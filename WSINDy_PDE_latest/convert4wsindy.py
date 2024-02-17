#!/usr/bin/env python
# coding: utf-8

import os; from pathlib import Path
import numpy as np
from scipy.io import loadmat, savemat
from argparse import ArgumentParser, Namespace

import sys; sys.path.append('../')
from misc import h5file

parser = ArgumentParser()

parser.add_argument('pfile', help="Estimated PDE solution", type=str, default='', nargs='?') 
parser.add_argument('-dp', '--data_path', help="Parent path storing datasets (template)", type=str, default="./datasets/")
parser.add_argument('-pn', '--pde_name', help="PDE name", type=str, default="burgers")
parser.add_argument('-lv', '--noise_lv', help="Noise level", type=int, default=0.0)
parser.add_argument('-s', '--save', help="decide weather to save or not?", action='store_true')

args: Namespace = parser.parse_args()

prediction_file = args.pfile; save = args.save; pde_name = args.pde_name; noise_lv = int(args.noise_lv); data_path = args.data_path
PDEFIND_data_path = "../Datasets/"
out_path = "./inputs"
if not Path(out_path).exists():
    os.mkdir(out_path)
    print(out_path, "created")

data = loadmat(os.path.join(PDEFIND_data_path, f"{pde_name}"))
template = loadmat(os.path.join(data_path, f"{pde_name}.mat"))
template["xs"][0, 0] = data["x"]
template["xs"][0, 1] = data["t"].reshape(1, -1)

if len(prediction_file) == 0:
    u_clean = data["usol"].real
    np.random.seed(0)
    if noise_lv > 0:
        un = u_clean + 0.01*np.abs(noise_lv)*(u_clean.std())*np.random.randn(u_clean.shape[0], u_clean.shape[1])
    out_fname = os.path.join(out_path, f"{pde_name}_noise{noise_lv}.mat") 
else:
    pf = Path(prediction_file)
    if pf.exists():
        suff = pf.suffix
        if '.mat' == suff:
            un = loadmat(prediction_file)["usol"]
        elif '.h5' == suff:
            un = h5file(prediction_file, mode='r', return_dict=True)["usol"]
        else:
            un = np.load(prediction_file)
    out_fname = os.path.join(out_path, Path(prediction_file).name.replace(suff, '.mat'))
    
template["U_exact"][0, 0] = un

keep_keys = set({"__header__", "__version__"", ""__globals__", 
                 "U_exact", "xs", "lhs", "tauhat", "tau"})

for k in template.copy().keys():
    if k not in keep_keys:
        template.pop(k)

if save:
    savemat(out_fname, template)
    print("save to", out_fname)
else:
    print("not save", out_fname)

# outputs = loadmat("./outputs/burgers_noise30_wsindy.mat")

# ground_coeffs = np.array([0.1, -0.5*2])
# est_coeffs = outputs["buffer_wsindy"][0, 0]
# est_coeffs = est_coeffs[np.abs(est_coeffs)>0]
# est_coeffs[1] = 2*est_coeffs[1]

# (1.7287405220410654, 1.2283569619815957)
# errs = np.abs(est_coeffs-ground_coeffs)/np.abs(ground_coeffs)
# errs = 100*errs
# errs.mean(), errs.std()
