from conans import ConanFile, CMake, tools
from conan.tools.cmake  import CMakeToolchain,CMakeDeps,CMake
import os

class ConanBuild(ConanFile):
    settings = 'build_type'
    generators = "CMakeDeps", "CMakeToolchain"
    requires = 'boost/1.72.0', 'eigen/3.3.7'


    def build(self):
      cmake = CMake(self)
      cmake.configure()
      cmake.build()
