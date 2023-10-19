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
      ${CMAKE_INSTALL_PREFIX}/share/doc/flecsi-sp
    OPTIONAL
  )

  if(NOT PAGES_REPO STREQUAL "")
    add_doc_deployment(deploy-docs PAGES_REPO ${PAGES_REPO})
  endif()
endif()

