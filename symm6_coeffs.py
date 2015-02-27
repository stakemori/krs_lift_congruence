# -*- coding: utf-8; mode: sage -*-

from sage.all import CuspForms, save, load
from kim_shahidi_lift_cong.symmetric_l_calculator import Symm6Calculator
import kim_shahidi_lift_cong.check_krs_cong
import os


def compute_coeffs(digits):
    for k in [16, 18, 20]:
        print k
        f = CuspForms(1, k).basis()[0]
        lcalc = Symm6Calculator(f, digits=digits)
        coeffs = lcalc.coeffs()
        save(coeffs, os.path.join(
            kim_shahidi_lift_cong.check_krs_cong.data_dir,
            "coeffs" + str(k) + ".sobj"))


def compute_l(digits):
    for k in [16, 18, 20]:
        print k
        f = CuspForms(1, k).basis()[0]
        coeffs = load(os.path.join(
            kim_shahidi_lift_cong.check_krs_cong.data_dir,
            "coeffs" + str(k) + ".sobj"))
        lcalc = Symm6Calculator(f, digits=digits, coeffs=coeffs)
        save(lcalc.algebraic_part_numeric(3*k - 2),
             os.path.join(kim_shahidi_lift_cong.check_krs_cong.data_dir,
                          "sym6_l_value" + str(k) + ".sobj"))
