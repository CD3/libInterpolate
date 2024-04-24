install-deps:
  conan install . -s build_type=Debug
  conan install . -s build_type=Release

configure: install-deps
  cmake --preset conan-default

build: configure
  cmake --build build --config Debug
  cmake --build build --config Release

test: build
  cd build && testing/Debug/libInterpolate_CatchTests
  cd build && testing/Release/libInterpolate_CatchTests

format:
  fd . src --type f --exec clang-format -i
