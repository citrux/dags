#!/usr/bin/env python3

import sys
import numpy as np

if __name__ == '__main__':
    if len(sys.argv) > 1:
        fname = sys.argv[1]
        data = np.loadtxt(fname)
        signal = np.sqrt(np.sum(data[:,:3] ** 2, axis=1))
        spectrum = np.fft.fft(signal)
        signal_rev = np.fft.ifft(spectrum)
        np.savetxt("spectrum_py.dat", np.vstack((spectrum.real, spectrum.imag)).transpose())
        np.savetxt("signal_rev_py.dat", np.vstack((signal_rev.real, signal_rev.imag)).transpose())
