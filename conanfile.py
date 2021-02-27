from conans import ConanFile, CMake, tools
import os

class ConanBuild(ConanFile):
    generators = "cmake", "virtualenv"
    requires = 'boost/1.69.0@conan/stable', 'eigen/3.3.7@cd3/devel'


    def build(self):
      cmake = CMake(self)
      cmake.configure()
      cmake.build()
