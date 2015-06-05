#!/usr/bin/env python

from __future__ import print_function
from __future__ import division
range = xrange

import numpy as np
import sys
import os
import string

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

full_data = None

data_sets = []

for i, output_choice in enumerate(output_choices):
    data = np.loadtxt(output_choice)
    if full_data is None:
        N_x, N_y = data.shape
        full_data = np.zeros((N_x, N_y, N_outputs+1))
    else:
        N_x_temp, N_y_temp = data.shape
        if (N_x != N_x_temp) or (N_y != N_y_temp):
            raise ValueError('Inconsistent shapes!')
    full_data[..., i+1] = data
    full_data[..., 0] += data

full_data[..., 0] = full_data[..., 0] / float(N_outputs)

np.savetxt('{}_average.dat'.format(guess_base_name), full_data[..., 0])
