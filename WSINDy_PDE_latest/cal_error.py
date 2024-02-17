#!/usr/bin/env python
# coding: utf-8
import os; from pathlib import Path
import numpy as np
from scipy.io import loadmat, savemat
from argparse import ArgumentParser, Namespace

parser = ArgumentParser()
parser.add_argument('pfile', help="An output file from wsindy algo", type=str, default='', nargs='?') 
parser.add_argument('-pn', '--pde_name', help="PDE name", type=str, default="burgers")
args: Namespace = parser.parse_args()

fname = args.pfile
pde_name = args.pde_name

if pde_name == 'burgers':
    outputs = loadmat(fname)

    ground_coeffs = np.array([0.1, -0.5*2])
    est_coeffs = outputs["buffer_wsindy"][0, 0]
    est_coeffs = est_coeffs[np.abs(est_coeffs)>0]
    est_coeffs[1] = 2*est_coeffs[1]

    errs = np.abs(est_coeffs-ground_coeffs)/np.abs(ground_coeffs)
    errs = 100*errs
    print(errs.mean(), errs.std())
