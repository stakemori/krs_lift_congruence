# -*- coding: utf-8; mode: sage -*-
import os

from vector_valued_const.const import ScalarModFormConst as SMFC
from vector_valued_const.const import ConstVectValued
from vector_valued_const.const import ConstVectValuedHeckeOp as CVH
from vector_valued_const.const import CalculatorVectValued
from degree2.basic_operation import number_of_procs

from degree2.vector_valued_smfs import VectorValuedSiegelModularForms

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


# p = 5518029068479

# with number_of_procs(4):
#     calculator.calc_forms_and_save(6, verbose=True)

# sym18_consts = [cvv([SMFC([4]), SMFC([4]), SMFC([12])], 1),
#                 cvv([SMFC([4]), SMFC([4]), SMFC([4, 4, 4])], 1),
#                 cvv([SMFC([4]), SMFC([4]), SMFC([6, 6])], 1),
#                 cvv([SMFC([4]), SMFC([6]), SMFC([4, 6])], 1),
#                 cvv([SMFC([4]), SMFC([6]), SMFC([10])], 1),
#                 cvv([SMFC([4]), SMFC([4, 4]), SMFC([4, 4])], 1),
#                 cvv([SMFC([4]), SMFC([10]), SMFC([6])], 1),
#                 cvv([SMFC([4]), SMFC([4, 6]), SMFC([6])], 1),
#                 cvv([SMFC([4]), SMFC([5, 6]), SMFC([5])], 1),
#                 cvv([SMFC([4]), SMFC([12]), SMFC([4])], 1),
#                 cvv([SMFC([4]), SMFC([4, 4, 4]), SMFC([4])], 1),
#                 cvv([SMFC([4]), SMFC([6, 6]), SMFC([4])], 1),
#                 cvv([SMFC([5]), SMFC([5]), SMFC([10])], 1),
#                 cvv([SMFC([5]), SMFC([5]), SMFC([4, 6])], 1),
#                 cvv([SMFC([5]), SMFC([6]), SMFC([4, 5])], 1),
#                 cvv([SMFC([5]), SMFC([4, 5]), SMFC([6])], 1),
#                 cvv([SMFC([5]), SMFC([4, 6]), SMFC([5])], 1),
#                 cvv([SMFC([5]), SMFC([10]), SMFC([5])], 1),
#                 cvv([SMFC([5]), SMFC([5, 6]), SMFC([4])], 1),
#                 cvv([SMFC([6]), SMFC([6]), SMFC([4, 4])], 1),
#                 cvv([SMFC([6]), SMFC([4, 4]), SMFC([6])], 1),
#                 cvv([SMFC([6]), SMFC([4, 5]), SMFC([5])], 1),
#                 cvv([SMFC([6]), SMFC([4, 6]), SMFC([4])], 1),
#                 cvv([SMFC([6]), SMFC([10]), SMFC([4])], 1),
#                 cvv([SMFC([4, 4]), SMFC([4, 4]), SMFC([4])], 1),
#                 # inc 3
#                 cvv([SMFC([4]), SMFC([4]), SMFC([4, 6])], 3),
#                 cvv([SMFC([4]), SMFC([4]), SMFC([10])], 3),
#                 cvv([SMFC([4]), SMFC([5]), SMFC([4, 5])], 3),
#                 cvv([SMFC([4]), SMFC([6]), SMFC([4, 4])], 3),
#                 cvv([SMFC([4]), SMFC([4, 4]), SMFC([6])], 3),
#                 cvv([SMFC([4]), SMFC([4, 5]), SMFC([5])], 3),
#                 cvv([SMFC([4]), SMFC([10]), SMFC([4])], 3),
#                 cvv([SMFC([4]), SMFC([4, 6]), SMFC([4])], 3),
#                 cvv([SMFC([5]), SMFC([5]), SMFC([4, 4])], 3),
#                 cvv([SMFC([5]), SMFC([4, 4]), SMFC([5])], 3),
#                 cvv([SMFC([5]), SMFC([4, 5]), SMFC([4])], 3),
#                 cvv([SMFC([6]), SMFC([6]), SMFC([6])], 3),
#                 cvv([SMFC([6]), SMFC([4, 4]), SMFC([4])], 3),
#                 # quadruple es4, inc 1
#                 cvv([SMFC([4]), SMFC([4]), SMFC([4]), SMFC([4, 4])], 1, 'a'),
#                 cvv([SMFC([4]), SMFC([6]), SMFC([4]), SMFC([6])], 1, 'a'),
#                 cvv([SMFC([4]), SMFC([4, 4]), SMFC([4]), SMFC([4])], 1, 'a'),
#                 cvv([SMFC([5]), SMFC([5]), SMFC([4]), SMFC([6])], 1, 'a'),
#                 cvv([SMFC([5]), SMFC([6]), SMFC([4]), SMFC([5])], 1, 'a'),
#                 cvv([SMFC([6]), SMFC([6]), SMFC([4]), SMFC([4])], 1, 'a'),
#                 cvv([SMFC([4]), SMFC([4]), SMFC([4]), SMFC([4, 4])], 1, 'b'),
#                 cvv([SMFC([4]), SMFC([6]), SMFC([4]), SMFC([6])], 1, 'b'),
#                 cvv([SMFC([4]), SMFC([4, 4]), SMFC([4]), SMFC([4])], 1, 'b'),
#                 cvv([SMFC([5]), SMFC([5]), SMFC([4]), SMFC([6])], 1, 'b'),
#                 cvv([SMFC([5]), SMFC([6]), SMFC([4]), SMFC([5])], 1, 'b'),
#                 cvv([SMFC([6]), SMFC([6]), SMFC([4]), SMFC([4])], 1, 'b'),
#                 # quadruple es4, inc 3
#                 cvv([SMFC([4]), SMFC([5]), SMFC([4]), SMFC([5])], 1, 'a'),
#                 cvv([SMFC([4]), SMFC([6]), SMFC([4]), SMFC([4])], 1, 'a'),
#                 cvv([SMFC([4]), SMFC([5]), SMFC([4]), SMFC([5])], 1, 'b'),
#                 cvv([SMFC([4]), SMFC([6]), SMFC([4]), SMFC([4])], 1, 'b'),
#                 # quadruple es6, inc 1
#                 cvv([SMFC([4]), SMFC([5]), SMFC([6]), SMFC([5])], 1, 'a'),
#                 cvv([SMFC([4]), SMFC([6]), SMFC([6]), SMFC([4])], 1, 'a'),
#                 cvv([SMFC([4]), SMFC([5]), SMFC([6]), SMFC([5])], 1, 'b'),
#                 cvv([SMFC([4]), SMFC([6]), SMFC([6]), SMFC([4])], 1, 'b')]

# def const_to_str(c):
#     consts = c.consts
#     tp = c._type
#     inc = c.inc
#     wts_s = [c.wts for c in consts]
#     consts_str = ", ".join(["SMFC(%s)"%(wts) for wts in wts_s])
#     return "cvv([{consts_str}], {inc}{tp})".format(
#         consts_str=consts_str,
#         inc=str(inc),
#         tp=", '%s'"%(tp) if tp is not None else "")
