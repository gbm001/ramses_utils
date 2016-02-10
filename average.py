#!/usr/bin/env python3

import numpy as np
import sys
import os
import string
from io import StringIO, BytesIO

output_choices = sys.argv[1:]
N_outputs = len(output_choices)

pwd = os.getcwd()
print(pwd, output_choices)

strip_chars = string.digits + '_.'
if output_choices[0].endswith(('.txt', '.dat')):
    guess_base_name = output_choices[0][:-4].rstrip(strip_chars)
else:
    guess_base_name = output_choices[0].rstrip(strip_chars)
if guess_base_name:
    print('output to {}_average.dat'.format(guess_base_name))
else:
    print('output to average.dat')
    guess_base_name = 'average.dat'

filesets = []
nblocks = None

for i, output_choice in enumerate(output_choices):
    # Load blocks of text from file [i]
    line_blocks = []
    data_blocks = []
    with open(output_choice, 'rt') as f:
        string_data = f.readlines()
    current_line_block = []
    for line in string_data:
        if line.strip():
            current_line_block.append(line)
        else:
            if current_line_block:
                line_blocks.append(current_line_block)
            current_line_block = []
    if current_line_block:
        line_blocks.append(current_line_block)

    # Read block of text into numpy array
    for line_block in line_blocks:
        if line_block:
            strbuffer = StringIO(''.join(line_block))
            data = np.loadtxt(strbuffer)
            data_blocks.append(data)
    filesets.append(data_blocks)
    
    # Check number of blocks is consistent with previous files
    if nblocks is None:
        nblocks = len(data_blocks)
    else:
        if nblocks != len(data_blocks):
            raise ValueError('Inconsistent number of blocks!')

# Allocate space for sum of data in blocks
full_data = [None]*nblocks

for fileset in filesets:
    # Loop over each file
    for i, data in enumerate(fileset):
        if full_data[i] is None:
            # Allocate new storage for block
            N_x, N_y = data.shape
            full_data[i] = np.zeros((N_x, N_y))
        else:
            # Check block shape has not changed
            N_x, N_y = full_data[i].shape
            N_x_temp, N_y_temp = data.shape
            if (N_x != N_x_temp) or (N_y != N_y_temp):
                raise ValueError('Inconsistent shapes!')
        full_data[i] += data

output_buffer = BytesIO()

for i in range(nblocks):
    full_data[i] = full_data[i] / float(N_outputs)
    np.savetxt(output_buffer, full_data[i])
    if i != nblocks:
        output_buffer.write(b'\n')

with open('{}_average.dat'.format(guess_base_name), 'wb') as f:
    f.write(output_buffer.getvalue())
