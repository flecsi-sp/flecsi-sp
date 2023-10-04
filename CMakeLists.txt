#------------------------------------------------------------------------------#
# Copyright (c) 2023, Triad National Security, LLC.
# All rights reserved.
#------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.19)

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

set(PAGES_REPO "git@re-git.lanl.gov:flecsi-sp/flecsi-sp.git"
  CACHE STRING "Pages repository for pages branch")
option(ENABLE_DOCUMENTATION "Enable documentation" OFF)

if(ENABLE_DOCUMENTATION)
  include(FleCSI/documentation)
  include(deploy)

  option(ENABLE_DOXYGEN "Enable Doxygen documentation" ON)
  option(ENABLE_SPHINX  "Enable Sphinx documentation"  ON)
  mark_as_advanced(ENABLE_DOXYGEN)
  mark_as_advanced(ENABLE_SPHINX)

  if(ENABLE_DOXYGEN)
    include(FleCSI/doxygen)

    # Create Doxygen output directory eagerly to avoid race condition.
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/doc/api)

    flecsi_add_doxygen_target(core
      CONFIGS
        ${CMAKE_SOURCE_DIR}/doc/doxygen/core.conf.in
        ${CMAKE_SOURCE_DIR}/doc/doxygen/conf.in
    )
    flecsi_add_doxygen_target(user
      CONFIGS
        ${CMAKE_SOURCE_DIR}/doc/doxygen/user.conf.in
        ${CMAKE_SOURCE_DIR}/doc/doxygen/conf.in
    )
  endif()

  if(ENABLE_SPHINX)
    include(FleCSI/sphinx)

    flecsi_add_sphinx_target(main
      CONFIG ${CMAKE_SOURCE_DIR}/doc/sphinx
      OUTPUT ${CMAKE_BINARY_DIR}/doc
    )

    # This is for local install of the documentation. To deploy the
    # documentation to the flecsi website, use the "deploy-documentation"
    # target.
    install(
      DIRECTORY
        ${CMAKE_BINARY_DIR}/doc/
      DESTINATION
        ${CMAKE_INSTALL_PREFIX}/share/doc/Burton
      OPTIONAL
    )
  endif()

  if((NOT PAGES_REPO STREQUAL "") AND ENABLE_SPHINX AND ENABLE_DOXYGEN)
    add_doc_deployment(deploy-docs PAGES_REPO ${PAGES_REPO})
  endif()
endif()

