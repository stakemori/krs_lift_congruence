# -*- coding: utf-8; mode: sage -*-

import os
from sage.all import (cached_method, QQ, mul, NumberField)

from degree2.const import ScalarModFormConst as SMFC
from degree2.const import ConstVectValued
from degree2.const import ConstVectValuedHeckeOp as CVH
from degree2.const import CalculatorVectValued

from degree2.basic_operation import number_of_procs
from degree2.vector_valued_smfs import VectorValuedSiegelModularForms
from degree2.scalar_valued_smfs import eisenstein_series_degree2
from degree2.utils import pmap


import krs_lift_congruence.utils

def cvv(cs, inc, tp=None):
    return ConstVectValued(14, cs, inc, tp)

sym14_21_consts = [cvv([SMFC([4]), SMFC([5]), SMFC([5, 6])], 1),
                   cvv([SMFC([4]), SMFC([6]), SMFC([4, 6])], 1),
                   cvv([SMFC([4]), SMFC([6]), SMFC([10])], 1),
                   cvv([SMFC([4]), SMFC([4, 4]), SMFC([4, 4])], 1),
                   cvv([SMFC([4]), SMFC([10]), SMFC([6])], 1),
                   cvv([SMFC([4]), SMFC([4, 6]), SMFC([6])], 1),
                   cvv([SMFC([4]), SMFC([5, 6]), SMFC([5])], 1),
                   cvv([SMFC([4]), SMFC([12]), SMFC([4])], 1),
                   cvv([SMFC([4]), SMFC([4, 4, 4]), SMFC([4])], 1),
                   cvv([SMFC([4]), SMFC([6, 6]), SMFC([4])], 1),
                   cvv([SMFC([5]), SMFC([6]), SMFC([4, 5])], 1),
                   cvv([SMFC([5]), SMFC([4, 5]), SMFC([6])], 1),
                   cvv([SMFC([5]), SMFC([4, 6]), SMFC([5])], 1),
                   cvv([SMFC([5]), SMFC([10]), SMFC([5])], 1),
                   cvv([SMFC([5]), SMFC([5, 6]), SMFC([4])], 1),
                   cvv([SMFC([6]), SMFC([4, 4]), SMFC([6])], 1),
                   cvv([SMFC([6]), SMFC([4, 5]), SMFC([5])], 1),
                   cvv([SMFC([6]), SMFC([4, 6]), SMFC([4])], 1),
                   cvv([SMFC([6]), SMFC([10]), SMFC([4])], 1),
                   cvv([SMFC([4]), SMFC([5]), SMFC([4, 5])], 3),
                   cvv([SMFC([4]), SMFC([6]), SMFC([4, 4])], 3),
                   cvv([SMFC([4]), SMFC([4, 4]), SMFC([6])], 3)]

def _sym14_21_consts1():
    l = sym14_21_consts[:2]
    return [CVH(c, 2) for c in l]

sym14_21_consts.extend(_sym14_21_consts1())

calculator = CalculatorVectValued(sym14_21_consts,
                                  krs_lift_congruence.utils.data_dir)
# calculator.calc_forms_and_save(6, verbose=True)

class VectorValuedSMFsSym14Wt17NonHol(VectorValuedSiegelModularForms):
    def __init__(self, prec):
        VectorValuedSiegelModularForms.__init__(self, 17, 14, prec)

    def dimension(self):
        return 24

    @cached_method
    def basis(self):
        d = calculator.forms_dict(self.prec)
        es4 = eisenstein_series_degree2(4, self.prec)
        f = es4**(-1)
        with number_of_procs(1):
            res = pmap(lambda c: d[c] * f, sym14_21_consts)
        return res

    def strum_bd_list(self):
        return 3

sym14_wt17_non_hol = VectorValuedSMFsSym14Wt17NonHol(6)



class VectorValuedSMFsSym14Wt17(VectorValuedSiegelModularForms):
    def __init__(self, prec):
        VectorValuedSiegelModularForms.__init__(self, 17, 14, prec)

    def dimension(self):
        return 13

    @cached_method
    def basis(self):
        M = VectorValuedSMFsSym14Wt17NonHol(self.prec)
        pl = mul([f for f, _ in M.hecke_charpoly(2).factor()
                  if f.degree() in [1, 12]])
        return M.basis_of_subsp_annihilated_by(pl, parallel=True)

def compute_lift_and_non_lift():
    krs_lift_congruence.utils.compute_lift_and_non_lift(
        sym14_wt17_non_hol, 6)


def mod_p(alpha, p):
    a = -4078080
    return krs_lift_congruence.utils.modulo_p(alpha, a, p)

def fname(f):
    return os.path.join(krs_lift_congruence.utils.data_dir, f)

def check_cong():
    p = 92467
    print "Checking when p = %s ..."%(p,)
    t2_eigenvalue = -4078080
    M = sym14_wt17_non_hol
    krs_lift_congruence.utils.check_cong(p, t2_eigenvalue, M)



# class VectorValuedSMFsSym14Wt17Lift(VectorValuedSiegelModularForms):
#     def __init__(self, prec):
#         VectorValuedSiegelModularForms.__init__(self, 17, 14, prec)

#     def dimension(self):
#         return 1

#     @cached_method
#     def basis(self):
#         M = VectorValuedSMFsSym14Wt17NonHol(self.prec)
#         pl = [f for f, _ in M.hecke_charpoly(2).factor() if f.degree() == 1][0]
#         return M.basis_of_subsp_annihilated_by(pl)

# def sym14wt17_lift(prec):
#     M = VectorValuedSMFsSym14Wt17Lift(prec)
#     return M.eigenform_with_eigenvalue_t2(QQ(-4078080))


# class VectorValuedSMFsSym14Wt17NonLift(VectorValuedSiegelModularForms):
#     def __init__(self, prec):
#         VectorValuedSiegelModularForms.__init__(self, 17, 14, prec)

#     def dimension(self):
#         return 12

#     @cached_method
#     def basis(self):
#         M = VectorValuedSMFsSym14Wt17NonHol(self.prec)
#         pl = [f for f, _ in M.hecke_charpoly(2).factor() if f.degree() == 12][0]
#         return M.basis_of_subsp_annihilated_by(pl, parallel=True)
