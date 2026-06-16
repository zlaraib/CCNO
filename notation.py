#!/bin/python
import sys
import argparse
import subprocess

dict_symbol2ascii = [
        ["Δ", "Delta_"],
        ["α", "alpha"],
        ["ν", "nu"],
        ["ₑ", "_e"],
        ["ₑ̄", "_anti_e"],
        ["ψ", "Psi"],
        ["ρ", "rho"],
        ["μ", "mu"],
        ["²", "_squared"],
        ["p̂", "phat"],
        ["ω", "omega"],
        ["π", "pi"],
        ["τ", "tau"],
        ["ξ", "xi"],
        ["ē", "antie"],
        
    ]

parser = argparse.ArgumentParser(
                    prog="notation.py",
                    description="Dictionary based notation swapper.",
                    epilog="WARNING: do NOT run recursively unless you are \
                    confident you will not overwrite important files! (hint: \
                    if the path will include this file itself, you are 100% \
                    going to break things!)")

# store_true defaults to storing false? okay.

parser.add_argument("-r", "--recurse", action="store_true",
                    help="recursively parse through given path")

# add this function later.
# known_dicts = ["symbol-ascii"]
# parser.add_argument("-d", "--dict", nargs=1, choices=known_dicts,
#                     help="specify unique dictionary to use")

parser.add_argument("-i", "--invert", action="store_true",
                    help="invert dictionary direction (A->B becomes B->A)")

parser.add_argument("filepath", help="file (dir) path")

args = parser.parse_args()


# This is adequate for now...
def sed_applydict(dictionary, path, recurse, invert):
    # for notational clarity, we convert from idx->fdx
    idx = invert # rename for clarity 0 if not inverted, 1 if inverted...
    fdx = 1 - invert # rename for clarity 0 if not inverted, 1 if inverted...
    for ent in dictionary:
        # Understanding this sed call takes... some effort.
        # Just, trust me it's only as dangerous as fireworks in a shed... :D
        sed_cmd = f"s/{ent[idx]}/{ent[fdx]}/g"
        if recurse:
            out = subprocess.run(["find", path, "-type", "f", "-name", "*",
                                  "-exec", "sed", "-i", sed_cmd, "{}", "+"],
                                 capture_output=True)
        else:
            out = subprocess.run(["sed", "-i", sed_cmd, path],
                                 capture_output=True)
            
        # did we fuck up?
        if out.returncode != 0:
            print("ERROR:")
            print(out)


if __name__ == "__main__":
    sed_applydict(dict_symbol2ascii, args.filepath, args.recurse, args.invert) # false is for "undo"
