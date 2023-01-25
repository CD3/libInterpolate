#! /bin/bash

set -e
root=$(git rev-parse --show-toplevel)

cd $root

bindir="$PWD/build-and-test"
function cleanup()
{
  rm -r $bindir
}
trap cleanup EXIT

function copy_bindir()
{
  rm -rf $bindir.error
  mv $bindir $bindir.error
}
trap copy_bindir ERR

echo "Checking that project can be installed."
mkdir $bindir
cd $bindir
conan install .. -s build_type=Release --build missing
cmake .. -DCMAKE_TOOLCHAIN_FILE=conan_toolchain.cmake -DCMAKE_BUILD_TYPE=Release -G "Ninja"
cmake --build .
cmake --build . --target test

# test install
cmake --install . --prefix $bindir/install





echo "Checking that installed project can be detected and used."
mkdir app
cd app

cat << EOF > main.cpp
#include <iostream>
#include <libInterpolate/Interpolate.hpp>

int main()
{
  _1D::LinearInterpolator<double> i1;

  return 0;
}
EOF

cat << EOF > CMakeLists.txt
cmake_minimum_required(VERSION 3.1)
project(demo)
add_executable( main main.cpp )
find_package( libInterpolate REQUIRED )
target_link_libraries(main libInterpolate::Interpolate )
EOF

mkdir build1
cd build1
conan install ${root} -s build_type=Release --build missing
cmake .. -DCMAKE_TOOLCHAIN_FILE=conan_toolchain.cmake -DlibInterpolate_DIR=${bindir}/install/lib/cmake/ -DCMAKE_BUILD_TYPE=Release -G "Ninja"
cmake --build .
./main

cd ..

