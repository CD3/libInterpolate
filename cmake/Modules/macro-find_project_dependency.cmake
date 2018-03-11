#
# Searches for a dependency.
# 
#
macro( find_project_dependency DEPENDENCY )
  message(STATUS "Searching for dependency: ${DEPENDENCY}")

  # parse arguments
  # any extra arguments passed by the caller will be put in find_project_dependency_UNPARSED_ARGUMENTS
  # ALL unparsed arguments are passed on directly to find_package
  cmake_parse_arguments( find_project_dependency
                         "REQUIRED"
                         ""
                         "PATHS;SUBDIRECTORY_PATHS;OVERRIDE_PATHS"
                         ${ARGN} )
  # assume the dependency was not found
  set( ${DEPENDENCY}_FOUND 0 )

  # if the dependency wasn't found
  if( NOT ${${DEPENDENCY}_FOUND} )
    if( find_project_dependency_OVERRIDE_PATHS )
      message(STATUS "Checking for overrides for ${DEPENDENCY}.")
      find_and_add_subdirectory(${DEPENDENCY} "${find_project_dependency_OVERRIDE_PATHS}")
    endif()
  endif()

  # search for the dependency
  # look in the directories given by PATHS first, if they were given
  if( NOT ${${DEPENDENCY}_FOUND} )
  if( find_project_dependency_PATHS )
    find_package( ${DEPENDENCY} QUIET ${find_project_dependency_UNPARSED_ARGS} PATHS ${find_project_dependency_PATHS} )
  endif()
  endif()

  # look in the system directories if dependency wasn't found
  if( NOT ${${DEPENDENCY}_FOUND} )
    find_package( ${DEPENDENCY} QUIET ${find_project_dependency_UNPARSED_ARGS} )
  endif()

  # if the dependency wasn't found
  if( NOT ${${DEPENDENCY}_FOUND} )
    message(STATUS "Could not find install for ${DEPENDENCY}. Checking for directory to include.")
    # look for a subdirectory to include
    find_and_add_subdirectory(${DEPENDENCY} "${find_project_dependency_SUBDIRECTORY_PATHS}")
  endif()

  # if the dependency wasn't found
  if( NOT ${${DEPENDENCY}_FOUND} )
    # and REQUIRED was given
    if( find_project_dependency_REQUIRED )
      # throw an error
      message(FATAL_ERROR "Could not find required dependency: ${DEPENDENCY}")
    endif()
  endif()

endmacro(find_project_dependency)

macro(find_and_add_subdirectory DEPENDENCY SUBDIRECTORY_PATHS )
    foreach(SUBDIR ${SUBDIRECTORY_PATHS})
      if(EXISTS ${SUBDIR}/${DEPENDENCY}/CMakeLists.txt)
        message(STATUS "Found ${DEPENDENCY} in ${SUBDIR} directory to include.")
        add_subdirectory( ${SUBDIR}/${DEPENDENCY} )
        set( ${DEPENDENCY}_FOUND 1 )
        break()
      endif()
    endforeach()
endmacro(find_and_add_subdirectory)
