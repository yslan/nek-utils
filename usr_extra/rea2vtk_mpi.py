import numpy as np
import pyvista as pv
import sys
import os
import struct
import time
from mpi4py import MPI

assert len(sys.argv) == 3, '\n\nUsage: python3 ./%s in.rea out.pvtu'%sys.argv[0]

fname = sys.argv[1]
fout = sys.argv[2]
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

if rank==0:
    print('MPI %d ranks'%size, flush=True)

def tic(s,pad=2,end=''):
    comm.Barrier()
    s0=' '*pad
    if (rank==0):
        print(s0+'%-20s'%s, end=end, flush=True)
    return time.perf_counter()
def toc(t0):
    comm.Barrier()
    if (rank==0):
        print('  done! (%.4e sec)'% (time.perf_counter() - t0), flush=True)

def reader_rea(fname):
    assert size==1, 'mpi ver is not implemented'

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

def reader_re2_mpi(fname):
    t00 = time.perf_counter()

    t0 = tic('read header ...')
    if rank == 0:
        with open(fname, 'rb') as f:
            hdr = f.read(80).decode("utf-8").split()
            nelgt, dim, nelgv = int(hdr[1]), int(hdr[2]), int(hdr[3])
            print('  hdr:',nelgt,dim,nelgv)
    
            wdsz = 8
            realtype = "d"
    
            # detect endianness
            etagb = f.read(4)
            etagL = struct.unpack("<f", etagb)[0]
            etagL = int(etagL * 1e5) / 1e5
            etagB = struct.unpack(">f", etagb)[0]
            etagB = int(etagB * 1e5) / 1e5
    
            if etagL == 6.54321:
                print('  little-endian')
                emode = "<"
                endian = "little"
            elif etagB == 6.54321:
                print('  big-endian')
                emode = ">"
                endian = "big"
            else:
                raise ValueError("Could not interpret endianness")

            nv = 2 ** dim

            meta = {
                "nelgt": nelgt,
                "dim": dim,
                "nelgv": nelgv,
                "wdsz": wdsz,
                "realtype": realtype,
                "emode": emode,
                "endian": endian,
                "nv": nv,
            }

            assert nelgt > size, 'need more elements than #ranks'
    else:
        meta = None

    meta = comm.bcast(meta, root=0)

    nelgt = meta["nelgt"]
    dim = meta["dim"]
    nelgv = meta["nelgv"]
    wdsz = meta["wdsz"]
    realtype = meta["realtype"]
    emode = meta["emode"]
    endian = meta["endian"]
    nv = meta["nv"]
    toc(t0)

    # partition
    base = nelgt // size
    rem = nelgt % size

    nelt = base + (1 if rank < rem else 0)
    e0 = rank * base + min(rank, rem)
    e1 = e0 + nelt

    nrec = dim * nv + 1
    bytes_per_elem = nrec * wdsz
    data_offset = 80 + 4 # hdr + endian chk
    my_offset = data_offset + e0 * bytes_per_elem
    my_nbytes = nelt * bytes_per_elem

    if rank == 0:
        print(f"  partition: {size}", flush=True)
        print(f"    partition: base={base}, rem={rem}", flush=True)
        print(f"    bytes/elem: {bytes_per_elem}", flush=True)

    t0 = tic('mpiio read ...')
    fh = MPI.File.Open(comm, fname, MPI.MODE_RDONLY)
    buf = np.empty(my_nbytes, dtype=np.uint8)
    fh.Read_at_all(my_offset, buf)

    fh.Close()
    toc(t0)

    t0 = tic('read from buf...')
    fi = np.frombuffer(
        buf,
        dtype = emode + realtype,
        count = (dim * nv + 1) * nelt,
        offset = 0
    )
    toc(t0)
    
    t0=tic('reshape...')
    fi = fi.reshape((nelt, nrec))
    group = fi[:, 0].copy()
    fi = fi[:, 1:].reshape((nelt, dim, nv))
    xyz = np.moveaxis(fi, [0, 1, 2], [1, 2, 0]).copy()
    toc(t0)
   
    comm.Barrier()
    if (rank==0): 
        print('Done! (%.4e sec)'%(time.perf_counter()-t00))

    meta_out = {
        "dim", dim,
        "nelgt", nelgt,
        "nelgv", nelgv,
        "nelt", nelt,
        "e0", e0,
        "e1", e1,
    }
    return xyz, meta_out

def read_mesh_mpi(fname):
    if (rank ==0):
        print('Reading mesh section in %s ...'%fname)
    if fname.endswith('.rea'):
        return reader_rea(fname)
    elif fname.endswith('.re2'):
        return reader_re2_mpi(fname)
    else:
        print('Unrecognized extension of %s, expecting .rea/.re2'%fname)
        sys.exit(1)
    
def dump_pvtu(xyz, fname):
    def _write_pvtu_manifest(
        pvtu_path,
        piece_names,
        dim,
        has_point_element_id=True,
        has_cell_element_id=False,
    ):
        """
        Write a minimal .pvtu manifest for the rank-local .vtu pieces.
        """
    
        byte_order = "LittleEndian" if sys.byteorder == "little" else "BigEndian"
    
        lines = []
        lines.append('<?xml version="1.0"?>')
        lines.append(
            f'<VTKFile type="PUnstructuredGrid" version="1.0" byte_order="{byte_order}" header_type="UInt64">'
        )
        lines.append('  <PUnstructuredGrid GhostLevel="0">')
    
        if has_point_element_id:
            lines.append('    <PPointData>')
            lines.append('      <PDataArray type="Int64" Name="elementId"/>')
            lines.append('    </PPointData>')
        else:
            lines.append('    <PPointData/>')
    
        if has_cell_element_id:
            lines.append('    <PCellData>')
            lines.append('      <PDataArray type="Int64" Name="elementId"/>')
            lines.append('    </PCellData>')
        else:
            lines.append('    <PCellData/>')
    
        lines.append('    <PPoints>')
        lines.append('      <PDataArray type="Float64" NumberOfComponents="3"/>')
        lines.append('    </PPoints>')
    
        for piece in piece_names:
            lines.append(f'    <Piece Source="{os.path.basename(piece)}"/>')
    
        lines.append('  </PUnstructuredGrid>')
        lines.append('</VTKFile>')
    
        with open(pvtu_path, "w") as f:
            f.write("\n".join(lines))
            f.write("\n")

    assert fname.endswith('.pvtu'), 'require fname format: out.pvtu'
    if rank==0:
        print('Writing mesh to %s ...'%fname, flush=True)
    t00 = time.perf_counter()

    froot, fext = os.path.splitext(fname)
    outdir = os.path.dirname(froot)
    stem = os.path.basename(froot)

    if outdir == "":
        outdir = "."
    if rank == 0:
        os.makedirs(outdir, exist_ok=True)
    comm.Barrier()
    
    nv, E_local, dim = xyz.shape
    assert nv in (4, 8), '#vertex = 4 or 8'
    assert dim in (2, 3), 'bad dim'
    
    E_counts = comm.allgather(E_local)
    eoff = sum(E_counts[:rank])

    piece_name = f"{stem}_{rank:04d}.vtu"
    piece_path = os.path.join(outdir, piece_name)
    pvtu_path = os.path.join(outdir, f"{stem}.pvtu")

    t0 = tic('allocate ...')
    points = np.zeros((E_local * nv, 3))
    xyz_t = np.transpose(xyz, (1, 0, 2))
    points[:,:dim] = xyz_t.reshape(E_local * nv, dim)
    toc(t0)

    t0 = tic('connectivity ...')
    connectivity = np.arange(E_local * nv, dtype=np.int64).reshape(E_local, nv)
    cells = np.hstack([np.full((E_local, 1), nv, dtype=np.int64), connectivity])
    toc(t0)

    t0 = tic('cell type ...')
    vtkCellType = 9 if dim==2 else 12
    cellType = np.full(E_local, vtkCellType, dtype=np.uint8)
    toc(t0)

    t0 = tic('pv UnstructuredGrid ...')
    grid = pv.UnstructuredGrid(cells, cellType, points)
    grid.point_data["elementId"] = np.repeat(
        np.arange(eoff + 1, eoff + E_local + 1, dtype=np.int64), nv
    )
    #grid.cell_data["elementId"] = np.arange(eoff + 1, eoff + E_local + 1, dtype=np.int64)
    toc(t0)

    t0 = tic('write pvtu ...')
    grid.save(piece_path)
    toc(t0)

    t0 = tic('write manifest ...')
    if (rank==0):
        _write_pvtu_manifest(
            pvtu_path=pvtu_path,
            piece_names=[f"{stem}_{r:04d}.vtu" for r in range(size)],
            dim=dim,
            has_point_element_id=True,
            has_cell_element_id=False,
        )
    toc(t0)
    if (rank==0):
        print('Done! (%.4e sec)'%(time.perf_counter()-t00), flush=True)

X, meta = read_mesh_mpi(fname)
dump_pvtu(X, fout)

