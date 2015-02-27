# -*- coding: utf-8; mode: sage -*-
from degree2.const import ScalarModFormConst as SMFC
from degree2.const import ConstVectValued
from degree2.const import ConstVectValuedHeckeOp as CVH
from degree2.const import CalculatorVectValued
from degree2.vector_valued_smfs import VectorValuedSiegelModularForms
import kim_shahidi_lift_cong.utils

def cvv(cs, inc, tp):
    return ConstVectValued(16, cs, inc, tp)


sym16_consts = [cvv([SMFC([4]), SMFC([4]), SMFC([10])], 1, None),
                cvv([SMFC([4]), SMFC([4]), SMFC([4, 6])], 1, None),
                cvv([SMFC([4]), SMFC([5]), SMFC([4, 5])], 1, None),
                cvv([SMFC([4]), SMFC([6]), SMFC([4, 4])], 1, None),
                cvv([SMFC([4]), SMFC([4, 4]), SMFC([6])], 1, None),
                cvv([SMFC([4]), SMFC([4, 5]), SMFC([5])], 1, None),
                cvv([SMFC([4]), SMFC([10]), SMFC([4])], 1, None),
                cvv([SMFC([4]), SMFC([4, 6]), SMFC([4])], 1, None),
                cvv([SMFC([5]), SMFC([5]), SMFC([4, 4])], 1, None),
                cvv([SMFC([5]), SMFC([4, 4]), SMFC([5])], 1, None),
                cvv([SMFC([5]), SMFC([4, 5]), SMFC([4])], 1, None),
                cvv([SMFC([6]), SMFC([6]), SMFC([6])], 1, None),
                cvv([SMFC([6]), SMFC([4, 4]), SMFC([4])], 1, None),
                cvv([SMFC([4]), SMFC([4]), SMFC([4, 4])], 3, None),
                cvv([SMFC([4]), SMFC([6]), SMFC([6])], 3, None),
                cvv([SMFC([4]), SMFC([4, 4]), SMFC([4])], 3, None),
                cvv([SMFC([5]), SMFC([5]), SMFC([6])], 3, None),
                cvv([SMFC([5]), SMFC([6]), SMFC([5])], 3, None),
                cvv([SMFC([6]), SMFC([6]), SMFC([4])], 3, None),
                cvv([SMFC([4]), SMFC([4]), SMFC([6]), SMFC([4])], 1, 'a'),
                cvv([SMFC([4]), SMFC([5]), SMFC([5]), SMFC([4])], 1, 'a')]

sym16_consts = sym16_consts + [CVH(sym16_consts[0], 2),
                               CVH(sym16_consts[1], 2)]

calculator = CalculatorVectValued(sym16_consts,
                                  kim_shahidi_lift_cong.utils.data_dir)
# prec is 6.
# calculator.calc_forms_and_save(6, verbose=True)


class VectorValuedSMFsSym16Wt19(VectorValuedSiegelModularForms):
    def __init__(self, prec):
        VectorValuedSiegelModularForms.__init__(self, 19, 16, prec)

    def dimension(self):
        return 23

    def basis(self):
        d = calculator.forms_dict(self.prec)
        return [d[c] for c in sym16_consts]


def check_cong():
    p = 37903031
    t2_eigenvalue = -8785920
    M = VectorValuedSMFsSym16Wt19(6)
    kim_shahidi_lift_cong.utils.check_cong(p, t2_eigenvalue, M)
