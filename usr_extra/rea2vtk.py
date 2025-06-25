import numpy as np
import pyvista as pv
import sys
import struct
import matplotlib.pyplot as plt

assert len(sys.argv) == 3, '\n\nUsage: python3 ./%s in.rea out.vtk'%sys.argv[0]

fname = sys.argv[1]
fout = sys.argv[2]

def reader_rea(fname):

    def read_xyz1(f):
        line = f.readline()
        return np.array([float(s) for s in line.strip().split()], dtype=np.float64)

    def read_elem2d(f):
        line = f.readline() # group
        x = read_xyz1(f)
        y = read_xyz1(f)
        return x, y

    def read_elem3d(f):
        x = np.zeros((8,))
        y = np.zeros((8,))
        z = np.zeros((8,))

        line = f.readline() # group
        x[:4] = read_xyz1(f)
        y[:4] = read_xyz1(f)
        z[:4] = read_xyz1(f)

        x[4:] = read_xyz1(f)
        y[4:] = read_xyz1(f)
        z[4:] = read_xyz1(f)

        return x, y, z

    max_lines = 400
    with open(fname) as f:
        start = False
        nline = 0
        while not start and nline < max_lines:
            nline += 1
            line = f.readline()
            if "MESH DATA" in line:
                start = True
        assert start, 'ERROR: cannot find MESH in %s'%f
        
        line = f.readline()
        stmp = line.strip().split()

        nelt = int(stmp[0])
        dim = int(stmp[1])
        nelv = int(stmp[2])
        print('hdr:',nelt,dim,nelv)
        assert dim==2 or dim==3, 'invalid dim %d'%dim

        nv = pow(2,dim)
        xyz = np.zeros((nv, nelt, dim))
        if dim==2:
            for e in range(nelt):
                x,y = read_elem2d(f)
                xyz[:,e,0] = x
                xyz[:,e,1] = y
        else:
            for e in range(nelt):
                x,y,z = read_elem3d(f)
                xyz[:,e,0] = x
                xyz[:,e,1] = y
                xyz[:,e,2] = z

    return xyz

def reader_re2(fname):
    with open(fname, 'rb') as f:
        hdr = f.read(80).decode("utf-8").split()
        nelt, dim, nelv = int(hdr[1]), int(hdr[2]), int(hdr[3])
        print('hdr:',nelt,dim,nelv)

        wdsz = 8
        realtype = "d"

        # detect endianness
        etagb = f.read(4)
        etagL = struct.unpack("<f", etagb)[0]
        etagL = int(etagL * 1e5) / 1e5
        etagB = struct.unpack(">f", etagb)[0]
        etagB = int(etagB * 1e5) / 1e5

        if etagL == 6.54321:
            print('little-endian')
            emode = "<"
            endian = "little"
        elif etagB == 6.54321:
            print('big-endian')
            emode = ">"
            endian = "big"
        else:
            raise ValueError("Could not interpret endianness")

        nv = pow(2, dim)
        buf = f.read((dim * nv + 1) * wdsz * nelt)
        xyz = np.zeros((nv, nelt, dim))

        fi = np.frombuffer(
            buf,
            dtype = emode + realtype,
            count = (dim * nv + 1) * nelt,
            offset = 0
        )
        fi = np.reshape(fi, (nelt, dim * nv + 1 ))
        group = fi[:,0] # not used

        fi = np.reshape(fi[:,1:], (nelt, dim, nv)) # discard group id
        xyz = np.moveaxis(fi, [0,1,2], [1,2,0])
        return xyz

def read_mesh(fname):
    print('reading mesh section in %s ...'%fname)
    if fname.endswith('.rea'):
        return reader_rea(fname)
    elif fname.endswith('.re2'):
        return reader_re2(fname)
    else:
        print('Unrecognized extension of %s, expecting .rea/.re2'%fname)
        sys.exit(1)
    
def plot_quad(xyz):
    from matplotlib.patches import Polygon

    nv, E, dim = xyz.shape

    fig, ax = plt.subplots()
    for e in range(E):
        coords = xyz[:, e, :]
        polygon = Polygon(coords, closed=True, edgecolor='black', facecolor='lightgray')
        ax.add_patch(polygon)

    ax.set_aspect('equal')
    ax.autoscale()
    plt.xlabel("X")
    plt.ylabel("Y")

def plot_quad_v2(xyz):

    nv, E, dim = xyz.shape

    xyzl = np.concatenate([xyz, xyz[0:1]], axis=0)

    # Plot all lines at once
    plt.plot(xyzl[:,:,0], xyzl[:,:,1], 'k-')
    plt.gca().set_aspect('equal')
    plt.xlabel("X")
    plt.ylabel("Y")


def dump_vtk(xyz, fname):
    nv, E, dim = xyz.shape
    assert nv in (4, 8), '#vertex = 4 or 8'
    assert dim in (2, 3), 'bad dim'

    points = np.zeros((E * nv, 3))
    xyz_t = np.transpose(xyz, (1, 0, 2))
    points[:,:dim] = xyz_t.reshape(E * nv, dim)

    connectivity = np.arange(E * nv, dtype=np.int64).reshape(E, nv)
    cells = np.hstack([np.full((E, 1), nv, dtype=np.int64), connectivity])

    vtkCellType = 9 if dim==2 else 12
    cellType = np.full(E, vtkCellType, dtype=np.uint8)

    grid = pv.UnstructuredGrid(cells, cellType, points)
    grid.point_data["elementId"] = np.repeat(np.arange(1, E + 1), nv)

    grid.save(fname)

X = read_mesh(fname)
dump_vtk(X, fout)


#plot_quad(X)
#plot_quad_v2(X)
#plt.show()
