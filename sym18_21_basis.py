# -*- coding: utf-8; mode: sage -*-
from degree2.const import ScalarModFormConst as SMFC
from degree2.const import ConstVectValued
from degree2.const import ConstVectValuedHeckeOp as CVH
from degree2.const import CalculatorVectValued

from degree2.vector_valued_smfs import VectorValuedSiegelModularForms

import krs_lift_congruence.utils

def cvv(cs, inc, tp=None):
    return ConstVectValued(18, cs, inc, tp)

sym18_consts1 = [cvv([SMFC([4]), SMFC([6]), SMFC([4, 6])], 1),
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
                 cvv([SMFC([4]), SMFC([4, 4]), SMFC([6])], 3),
                 cvv([SMFC([4]), SMFC([4, 5]), SMFC([5])], 3),
                 cvv([SMFC([4]), SMFC([10]), SMFC([4])], 3),
                 cvv([SMFC([4]), SMFC([4, 6]), SMFC([4])], 3),
                 cvv([SMFC([5]), SMFC([4, 4]), SMFC([5])], 3),
                 cvv([SMFC([5]), SMFC([4, 5]), SMFC([4])], 3),
                 cvv([SMFC([6]), SMFC([4, 4]), SMFC([4])], 3),
                 cvv([SMFC([4]), SMFC([6]), SMFC([4]), SMFC([6])], 1, 'a'),
                 cvv([SMFC([5]), SMFC([6]), SMFC([4]), SMFC([5])], 1, 'a'),
                 cvv([SMFC([4]), SMFC([5]), SMFC([6]), SMFC([5])], 1, 'a'),
                 cvv([SMFC([4]), SMFC([5]), SMFC([6]), SMFC([5])], 1, 'b')]

sym18_consts1.extend([CVH(c, 2) for c in sym18_consts1[:8]])

calculator = CalculatorVectValued(sym18_consts1,
                                  krs_lift_congruence.utils.data_dir)

# calculator1 = CalculatorVectValued(sym18_consts1[10:31],
#                                    krs_lift_congruence.utils.data_dir)


class VectorValuedSMFsSym18Wt21(VectorValuedSiegelModularForms):
    def __init__(self, prec):
        VectorValuedSiegelModularForms.__init__(self, 21, 18, prec)

    def dimension(self):
        return 39

    def basis(self):
        d = calculator.forms_dict(self.prec)
        return [d[_c] for _c in sym18_consts1]

def compute_lift_and_non_lift():
    M = VectorValuedSMFsSym18Wt21(6)
    krs_lift_congruence.utils.compute_lift_and_non_lift(M, 6)

def check_cong():
    t2_eigenvalue = -383331840
    M = VectorValuedSMFsSym18Wt21(6)
    for p in [103, 5518029068479]:
        print "checking when p = %s ... "%(p)
        krs_lift_congruence.utils.check_cong(p, t2_eigenvalue, M)


# check_cong()                    # noerror!!

# with number_of_procs(4):
#     calculator.calc_forms_and_save(6, verbose=True)
