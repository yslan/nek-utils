#!/bin/bash

if [ "$1" == "clean" ]; then
  set -x
  rm *.dra *.jou session.name fort.* 2> /dev/null
#  rm fluid.rea solid.rea combine.rea out.rea out.re2 2> /dev/null
  exit 0
fi

: ${NEK5000_PATH:="$HOME/Nek5000"}
if [ ! -d "$NEK5000_PATH" ]; then
  echo "Nek5000 not found at ${NEK5000_PATH}!"
  exit 1
fi

tools=(genbox pretex reatore2)
for f in ${tools[@]}; do
  str="$f=$NEK5000_PATH/bin/$f"
  eval $str
  if [ ! -f ${!f} ]; then  # ${!variable} is a variable variables.
    echo "Cannot find" ${!f}
    exit 1
  fi
done

# no gui mode on server
: ${GUI:=1}
if [ $GUI == 0 ]; then
  if ! command -v xvfb-run 2>&1 >/dev/null; then
    echo "A NO-GUI mode requires Xvfb (xvfb-run)!"
    exit 1
  fi
  pretex="xvfb-run "$pretex
fi

############################
# Main Meshing Starts Here #
############################

$genbox << EOF
fluid.box
EOF
mv box.rea fluid.rea

$genbox << EOF
solid.box
EOF
mv box.rea solid.rea

$pretex << EOF
combine
3
fluid
solid
2
EOF

$reatore2 << EOF
combine
out
EOF

