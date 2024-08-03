list:
  just --list

install-deps:
  conan install . -s build_type=Debug --build missing
  conan install . -s build_type=Release --build missing

configure: install-deps
  cmake --preset conan-default

build: configure
  cmake --build build --config Debug
  cmake --build build --config Release

test: build
  cd build && testing/Debug/libInterpolate_CatchTests
  cd build && testing/Release/libInterpolate_CatchTests

install-deps-app:
  conan install applications -s build_type=Debug --build missing
  conan install applications -s build_type=Release --build missing

configure-app: install-deps-app
  cmake applications --preset conan-default -DBUILD_TESTS=OFF

build-app: configure-app
  cmake --build applications/build --config Release

format:
  fd . src --type f --exec clang-format -i

clean:
  rm -rf build CMakeUserPresets.json applications/build applications/CMakeUserPresets.json
