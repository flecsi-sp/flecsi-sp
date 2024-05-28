#------------------------------------------------------------------------------#
# Copyright (c) 2023, Triad National Security, LLC.
# All rights reserved.
#------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.19)

#------------------------------------------------------------------------------#
# CMake policies.
#------------------------------------------------------------------------------#

cmake_policy(SET CMP0144 NEW) # ``PACKAGE_NAME''_ROOT

#------------------------------------------------------------------------------#
# Add local module path.
#------------------------------------------------------------------------------#

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

#------------------------------------------------------------------------------#
# Project.
#------------------------------------------------------------------------------#

project(FleCSI-SP LANGUAGES C CXX)

#------------------------------------------------------------------------------#
# Find the FleCSI package. This must be done early so that the FleCSI cmake
# packages will be available to the project, e.g., 'include(FleCSI/format)'.
#
# CMake documentation for FleCSI client applications is included in the
# FleCSI userguide.
#------------------------------------------------------------------------------#

if(NOT FleCSI_FOUND)
  find_package(FleCSI 2 REQUIRED)
endif()

#------------------------------------------------------------------------------#
# Documentation
#------------------------------------------------------------------------------#

option(ENABLE_DOCUMENTATION "Enable documentation" OFF)
set(GITHUB_PAGES_REPO "git@github.com:flecsi-sp/flecsi-sp.git"
    CACHE STRING "GitHub Pages repository for gh-pages branch")
mark_as_advanced(ENABLE_DOCUMENTATION)
mark_as_advanced(GITHUB_PAGES_REPO)

if(ENABLE_DOCUMENTATION)
  include(FleCSI/documentation)
  include(FleCSI/doxygen)
  include(FleCSI/sphinx)
  include(subdirlist)

  # Avoid race condition by proactively creating output directory
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/.doxygen/xml)
  flecsi_add_doxygen_target(sp
    CONFIGS
      ${CMAKE_SOURCE_DIR}/doc/doxygen/sp.conf.in
      ${CMAKE_SOURCE_DIR}/doc/doxygen/conf.in
  )

  file(COPY ${CMAKE_SOURCE_DIR}/doc/sphinx
    DESTINATION ${CMAKE_BINARY_DIR}
  )
  file(REMOVE_RECURSE ${CMAKE_BINARY_DIR}/.sphinx)
  file(RENAME ${CMAKE_BINARY_DIR}/sphinx
    ${CMAKE_BINARY_DIR}/.sphinx
  )

  make_subdirlist(_SUBDIRECTORIES ${CMAKE_SOURCE_DIR}/sp False)
  foreach(_SUBDIR ${_SUBDIRECTORIES})
    file(COPY ${CMAKE_SOURCE_DIR}/sp/${_SUBDIR}/doc/sphinx/index.rst
      DESTINATION ${CMAKE_BINARY_DIR}/.sphinx/${_SUBDIR}
    )
    file(COPY ${CMAKE_SOURCE_DIR}/sp/${_SUBDIR}/doc/sphinx/${_SUBDIR}
      DESTINATION ${CMAKE_BINARY_DIR}/.sphinx/${_SUBDIR}
    )

    # Create output directory
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/.sphinx/${_SUBDIR}/${_SUBDIR}/src)
    execute_process(COMMAND
      ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_SOURCE_DIR}/sp/${_SUBDIR}/spec
        ${CMAKE_BINARY_DIR}/.sphinx/${_SUBDIR}/${_SUBDIR}/src/spec
    )
    execute_process(COMMAND
      ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_SOURCE_DIR}/sp/${_SUBDIR}/examples
        ${CMAKE_BINARY_DIR}/.sphinx/${_SUBDIR}/${_SUBDIR}/src/examples
    )
  endforeach()

  flecsi_add_sphinx_target(main
    CONFIG ${CMAKE_BINARY_DIR}/.sphinx
    OUTPUT ${CMAKE_BINARY_DIR}/doc
  )

  if(NOT GITHUB_PAGES_REPO STREQUAL "")
    flecsi_add_doc_deployment(deploy-docs
      GITHUB_PAGES_REPO ${GITHUB_PAGES_REPO})
  endif()
endif()

add_subdirectory(sp/fdm)
add_subdirectory(sp/fvm)
