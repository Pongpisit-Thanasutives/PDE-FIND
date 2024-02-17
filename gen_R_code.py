import argparse
import numpy as np
from misc import h5file

def gen_R(f_path, target_name, data_name):
    data = h5file(f_path, 'r', return_dict=True)
    best_subsets = data["best_subsets"]
    for i, bs in enumerate(best_subsets):
        efi = ['V' + str(e + 1) for e in np.where(bs > 0)[0].tolist()]
        inp_string = ' + '.join(efi)
        out_string = f"{target_name} ~ {inp_string}"
        print(f"form{i+1} <- {out_string} + 0")
        out_string = (
            f"fit{i + 1} = brm({out_string}, data={data_name}, family=gaussian(), save_pars=save_pars(all=TRUE))"
        )
        # print(out_string)
    return 

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate R code')
    parser.add_argument('f_path', type=str, help='Path to the h5 file')
    parser.add_argument('--target', type=str, default='u_t', help='Target name')
    parser.add_argument('--data', type=str, default='derivatives', help='Data name')
    args = parser.parse_args()
    gen_R(args.f_path, args.target, args.data)
