#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([msp-boxes], [0.76])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/ordering/ordering-io.hpp])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
#AC_CHECK_HEADER_STDBOOL

# Checks for library functions.
AC_FUNC_STRTOD
AC_CHECK_FUNCS([sqrt strtol])

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/boxes/Makefile
                 src/ordering/Makefile
                 test/Makefile])
AC_OUTPUT
