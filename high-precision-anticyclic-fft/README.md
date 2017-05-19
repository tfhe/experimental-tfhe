# High precision FFT on the anticyclic ring

This is an implementation of the Fast Fourier Transform for the
anticyclic ring R[X]/X^N+1 where N is a power of 2.

There are two versions: one 128-bit fixed-point version, using a custom
implementation based on int128 integers

One gmp version on arbitrary precision.

The aim of this project is to measure the actual running time of such
FFTS, which are the elementary operations in TFHE.
