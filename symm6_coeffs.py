# -*- coding: utf-8; mode: sage -*-

from sage.all import CuspForms, save, load
from kim_shahidi_lift_cong.symmetric_l_calculator import Symm6Calculator

import os

base_dir = os.path.join(os.getenv("HOME"), "Documents/misc")
dir_dict = {16: os.path.join(base_dir, "sym14_basis"),
            18: os.path.join(base_dir, "sym16_basis"),
            20: os.path.join(base_dir, "sym18_basis")}


def compute_coeffs(digits):
    for k in [16, 18, 20]:
        print k
        f = CuspForms(1, k).basis()[0]
        lcalc = Symm6Calculator(f, digits=digits)
        coeffs = lcalc.coeffs()
        save(coeffs, os.path.join(dir_dict[k], "coeffs" + str(k) + ".sobj"))


def compute_l(digits):
    for k in [16, 18, 20]:
        print k
        f = CuspForms(1, k).basis()[0]
        coeffs = load(os.path.join(dir_dict[k], "coeffs" + str(k) + ".sobj"))
        lcalc = Symm6Calculator(f, digits=digits, coeffs=coeffs)
        save(lcalc.algebraic_part_numeric(3*k - 2),
             os.path.join(dir_dict[k], "sym6_l_value" + str(k) + ".sobj"))
