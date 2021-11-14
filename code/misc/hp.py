#!/usr/bin/env python
# Description: conversion of Counterr 3D numpy array data into a 2D dataframe

import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', default='hist_len_hp.npy')
parser.add_argument('-o', default='hp.csv')
args = parser.parse_args()

arr = np.load(args.i)
arr.transpose((0,1,2)).ravel()
arrReshaped = arr.reshape(arr.shape[0], -1)

np.savetxt(args.o, arrReshaped, fmt="%d", delimiter=",")