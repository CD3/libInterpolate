  $ git clone $(dirname $(dirname $TESTDIR) )
  Cloning into 'libInterpolate'...
  done.
  $ cd libInterpolate
  $ mkdir build
  $ cd build
  $ conan install .. > /dev/null 2>&1
  $ bash -c "source activate.sh && cmake -DCMAKE_INSTALL_PREFIX=$PWD/INSTALL -DBUILD_TESTS=OFF .." > /dev/null
  $ cmake --build . --target install > /dev/null
  $ ls INSTALL
  cmake
  include
  $ ls testBin
  ls: cannot access 'testBin': No such file or directory
  [2]
  $ cd ..
  $ mkdir build2
  $ cd build2
  $ conan install .. > /dev/null 2>&1
  $ bash -c "source activate.sh && cmake -DCMAKE_INSTALL_PREFIX=$PWD/INSTALL .." > /dev/null
  $ cmake --build . --target install > /dev/null
  $ ls INSTALL
  cmake
  include
  $ ls testBin
  libInterpolate_CatchTests
