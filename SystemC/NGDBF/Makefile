#
# Makefile for small to medium SystemC applications
#
# Copyright (c) 2006-2010 OFFIS e.V. Institute for Information Technology
# All rights reserved.
#
# Author: Philipp A. Hartmann <philipp.hartmann@offis.de>
#
# -----------------------------------------------------------------------
#
# This Makefile uses, and thus requires the generic SystemC
# systemc.mk, which contains the actual build rules.
#
# Within this file, only configuration variables should be changed.
#
# -----------------------------------------------------------------------
# The main make targets, that are available by default, include:
#
# all
#    - build the application
# [MODULE]
#    - build and run the application
#      (name determined by MODULE variable below)
# clean
#    - cleanup generated files
#
# -----------------------------------------------------------------------
#
# Configuration section of this SystemC project
# (if a variable is omitted here, a default value is used)
#

# path (and name) of required "systemc.mk"
SYSTEMC_MAKE := ./systemc.mk

# general build options (yes/no)
DEBUG=yes
#PROFILE=yes
#VERBOSE=yes

# -----------------------------------------------------------------------
# SystemC installation (required)
#
# The variable SYSTEMC_HOME is required to determine the
# location of the SystemC installation. If set in the
# environment, this setting is used.
#SYSTEMC_HOME := /path/to/systemc/directory

# The variable SYSTEMC_LIB is required to locate the static
# library file libsystemc.a, matching the current commpiler
# and SystemC headers.
#SYSTEMC_LIB := $(SYSTEMC_HOME)/lib-$(TARGET_ARCH)
TARGET_ARCH := linux64
SYSTEMC_HOME := /usr/local/systemc-2.2.0
SYSTEMC_LIB := $(SYSTEMC_HOME)/lib-$(TARGET_ARCH)

#SYSTEMC_HOME := /Library/SystemC/Current
#SYSTEMC_LIB := $(SYSTEMC_HOME)/lib-macosx


# -----------------------------------------------------------------------
# Name of the SystemC application (required!)
MODULE=ldpcsim

# Default parameters to pass to test application
# (used, if "make run" is called)
#ARGS:=42

# -----------------------------------------------------------------------
# source file extension
#
# Set the extension of C++ source files.
# Defaults to "<name>.cpp"
#SRCEXT := cpp

# The actual source code files can be specified via a list of
# source directories (in the variable SRCDIRS) OR as a list of
# .cpp-files (in the variable SRCS). If set, the variable SRCS takes
# precedence, otherwise, all .cpp-files in SRCDIRS are taken as
# input files (see also: $(SRCEXT) above).
#
# Source-Directories
SRCDIRS := ./src
# Source files
#SRCS := example1.cpp example2.cpp

# -----------------------------------------------------------------------
#
# Other variables, that can be set here
#
# additional include options
EXTRA_INCLUDES := -I./inc/
#
# additional library directories and link directives
#EXTRA_LIBDIRS := -L/path/to/my/lib/dir
#EXTRA_LIBS    := -lmylib

# additional preprocessor symbols to define
# (as list of -Dmacro[=defn])
#EXTRA_DEFINES := -DSC_INCLUDE_DYNAMIC_PROCESSES -Dever=;;

# -----------------------------------------------------------------------
# compiler to use
#CXX   := =g++
#
# compiler flags (defaults shown below)
#
# always added (high warning level):
#CXXFLAGS := -arch i386
#CXXFLAGS := -Wall -Wextra -pedantic \
#            -Wno-variadic-macros -Wno-long-long
# only added to optimising builds
#CXXFLAGS_OPT   := -O3
# only added to DEBUG builds
#CXXFLAGS_DEBUG := -g
# always added (useful for command-line overrides)
#CXXFLAGS_EXTRA :=

# additional target to clean up current test application
#EXTRA_CLEAN=extra-clean
#
#extra-clean:
#	$(DEL) my-test-file.out

# -----------------------------------------------------------------------
# look for common build rules in generic places
#
SYSTEMC_MAKE += \
  ./systemc.mk   \
  ../systemc.mk  \
   $(SYSTEMC_HOME)/systemc.mk \
   $(SYSTEMC_HOME)/examples/systemc.mk
SYSTEMC_MAKE := $(word 1,$(wildcard $(SYSTEMC_MAKE)))
# include generic OSSS Makefile
ifneq (,$(strip $(SYSTEMC_MAKE)))
include $(SYSTEMC_MAKE)
else
$(error Common SystemC Makefile 'systemc.mk' not found! Bailing out)
endif

# Taf!
