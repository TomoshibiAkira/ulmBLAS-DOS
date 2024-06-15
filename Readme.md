## Warning

You've been warned, it's very hacky. I don't plan to use LAPACK so I'm not compiling that for now.

ulmBLAS is compilable under DJGPP. Run ``make cblas`` to build.

An example is provided in ``dos_cblas_test`` folder: ``test_sgemm.c``. After ulmBLAS is built, run ``cd dos_cblas_test && make dos`` to built ``tsgemm.exe``.

Grab a copy of ``CWSDPMI`` and you're good to go.

## Benchmarks

The command is ``tsgemm 100 1000 1000``, which is around 200MFLOPs in total. All CPU tests are done in PCem.

* Intel 486DX4/120: 27.2s
* Cyrix 6x86L-PR133+: 13.5s
* Cyrix 6x86-PR200+: 9.9s
* Intel Pentium 100/66: 6.8s
* Intel Pentium 200: 3.4s
* AMD K6-2/266: 2.3s
* AMD K6-2/300: 2.0s
* Intel Pentium II/400: 0.7s
* Intel Celeron 533: 0.6s
* DOSBox 0.74-3 100% cycle: 0.7s

## Original README

This library is part of my lecture "Software Basics for High Performance
Computing" (MATH9367) at Ulm University:

   http://www.mathematik.uni-ulm.de/~lehn/ulmBLAS
   
   http://www.mathematik.uni-ulm.de/~lehn/sghpc

And yes, I am particularly proud of the section demonstrating how to achieve peak performance for the matrix-matrix product:

   http://www.mathematik.uni-ulm.de/~lehn/sghpc/gemm/index.html

Note: Further development will take place in ulmBLAS-core
