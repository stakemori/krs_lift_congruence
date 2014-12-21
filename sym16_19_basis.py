# -*- coding: utf-8; mode: sage -*-
import os

from vector_valued_const.const import ScalarModFormConst as SMFC
from vector_valued_const.const import ConstVectValued
from vector_valued_const.const import ConstVectValuedHeckeOp as CVH
from vector_valued_const.const import CalculatorVectValued

sym16_data_dir = os.path.join(os.getenv("HOME"), "Documents/misc/sym16_basis")

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

calcular = CalculatorVectValued(sym16_consts, sym16_data_dir)
# calcular.calc_forms_and_save(5, verbose=True)
