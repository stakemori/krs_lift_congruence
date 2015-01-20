# -*- coding: utf-8; mode: sage -*-
import os
from sage.all import flatten, PolynomialRing, FiniteField

from vector_valued_const.const import ScalarModFormConst as SMFC
from vector_valued_const.const import ConstVectValued
from vector_valued_const.const import ConstVectValuedHeckeOp as CVH
from vector_valued_const.const import CalculatorVectValued

from degree2.elements import SymWtModFmElt

from degree2.vector_valued_smfs import VectorValuedSiegelModularForms

import kim_shahidi_lift_cong.utils

sym18_data_dir = os.path.join(os.getenv("HOME"), "Documents/misc/sym18_basis")

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

calculator = CalculatorVectValued(sym18_consts1, sym18_data_dir)

calculator1 = CalculatorVectValued(sym18_consts1[10:31], sym18_data_dir)


class VectorValuedSMFsSym18Wt21(VectorValuedSiegelModularForms):
    def __init__(self, prec):
        VectorValuedSiegelModularForms.__init__(self, 21, 18, prec)

    def dimension(self):
        return 39

    def basis(self):
        d = calculator.forms_dict(self.prec)
        return [d[_c] for _c in sym18_consts1]

def mod_p(alpha, p):
    a = -383331840
    return sum([b * a**i for i, b in enumerate(alpha.list())])%p


def check_cong():
    def fname(f):
        return os.path.join(sym18_data_dir, f)

    lift = SymWtModFmElt.load_from(fname("lift_prec6.sobj"))
    non_lift = SymWtModFmElt.load_from(fname("non_lift_prec6.sobj"))

    M = VectorValuedSMFsSym18Wt21(6)
    g = lift - non_lift
    K = non_lift.base_ring
    v = M._to_vector(g)
    for p in [103, 5518029068479]:
        print "checking for p = %s ... "%(p)
        R = PolynomialRing(FiniteField(p), names="x")
        pl_modp = R(K.polynomial())
        # p is unramified.
        assert all((b == 1 for a, b in pl_modp().factor()))
        l = flatten([b.list() for b in v])
        l = [b.denominator() for b in l]
        assert all(b%p != 0 for b in l)
        assert all(mod_p(b, p)%p == 0 for b in v)

# check_cong()                    # noerror!!

# with number_of_procs(4):
#     calculator.calc_forms_and_save(6, verbose=True)
