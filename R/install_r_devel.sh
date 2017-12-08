## Output of ./configure alone:

R is now configured for x86_64-apple-darwin17.2.0

  Source directory:          .
  Installation directory:    /usr/local

  C compiler:                gcc  -g -O2
  Fortran 77 compiler:       gfortran  -g -O2

  Default C++ compiler:      g++   -g -O2
  C++98 compiler:            g++  -g -O2
  C++11 compiler:            g++ -std=gnu++11 -g -O2
  C++14 compiler:            g++ -std=gnu++14 -g -O2
  C++17 compiler:            g++ -std=gnu++1z -g -O2
  Fortran 90/95 compiler:    gfortran -g -O2
  Obj-C compiler:	     gcc -g -O2 -fobjc-exceptions

  Interfaces supported:      X11, aqua, tcltk
  External libraries:        readline, curl
  Additional capabilities:   NLS, ICU
  Options enabled:           shared BLAS, R profiling

  Capabilities skipped:      PNG, JPEG, TIFF, cairo
  Options not enabled:       memory profiling

  Recommended packages:      yes

## Relevant options:
# --with-libpng
# --with-jpeglib
# --with-libtiff
# --with-cairo
# --enable-R-profiling

## Relevant command:
# ./configure --with-libpng --with-jpeglib --with-libtiff --with-cairo --enable-R-profiling
