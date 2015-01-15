#!/usr/bin/env python3

import sys
import numpy as np
import struct

filename = sys.argv[1]

with open(filename, 'rb') as f:
    ndim, cell_size, ivar_min, ivar_max = struct.unpack('4i', f.read(16))
    
    cells = np.fromfile(f, dtype=np.float64)
    
cells_shape = [cell_size, cell_size, cell_size, ivar_max-ivar_min+1]
if (ndim==1):
    cells_shape[0] = 1
    cells_shape[1] = 1
if (ndim==2):
    cells_shape[0] = 1
cells = cells.reshape(cells_shape, order='F')

if filename.endswith('.dat'):
    filename = filename[:-4]

filename = filename + '.npy'

np.save(filename, cells)
