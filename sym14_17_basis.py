# -*- coding: utf-8; mode: sage -*-

import os

from vector_valued_const.const import ScalarModFormConst as SMFC
from vector_valued_const.const import ConstVectValued
from vector_valued_const.const import ConstVectValuedHeckeOp as CVH
from vector_valued_const.const import CalculatorVectValued

from degree2.vector_valued_smfs import VectorValuedSiegelModularForms

sym14_data_dir = os.path.join(os.getenv("HOME"), "Documents/misc/sym14_basis")

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
                   cvv([SMFC([4]), SMFC([4, 4]), SMFC([6])], 3),
                   cvv([SMFC([4]), SMFC([4, 5]), SMFC([5])], 3),
                   cvv([SMFC([4]), SMFC([10]), SMFC([4])], 3),
                   cvv([SMFC([4]), SMFC([4, 6]), SMFC([4])], 3),
                   cvv([SMFC([5]), SMFC([4, 4]), SMFC([5])], 3),
                   cvv([SMFC([5]), SMFC([4, 5]), SMFC([4])], 3),
                   cvv([SMFC([6]), SMFC([4, 4]), SMFC([4])], 3),
                   cvv([SMFC([4]), SMFC([6]), SMFC([4]), SMFC([6])], 1, 'a'),
                   cvv([SMFC([4]), SMFC([4, 4]), SMFC([4]), SMFC([4])], 1, 'a'),
                   cvv([SMFC([5]), SMFC([6]), SMFC([4]), SMFC([5])], 1, 'a'),
                   cvv([SMFC([4]), SMFC([5]), SMFC([6]), SMFC([5])], 1, 'a'),
                   cvv([SMFC([4]), SMFC([6]), SMFC([6]), SMFC([4])], 1, 'a'),
                   cvv([SMFC([4]), SMFC([6]), SMFC([4]), SMFC([6])], 1, 'b'),
                   cvv([SMFC([4]), SMFC([4, 4]), SMFC([4]), SMFC([4])], 1, 'b'),
                   cvv([SMFC([5]), SMFC([6]), SMFC([4]), SMFC([5])], 1, 'b'),
                   cvv([SMFC([4]), SMFC([5]), SMFC([6]), SMFC([5])], 1, 'b'),
                   cvv([SMFC([4]), SMFC([6]), SMFC([6]), SMFC([4])], 1, 'b')]

calculator = CalculatorVectValued(sym14_21_consts, sym14_data_dir)
# calculator.calc_forms_and_save(6, verbose=True)
