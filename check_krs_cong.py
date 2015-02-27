# -*- coding: utf-8; mode: sage -*-
import os
import kim_shahidi_lift_cong.sym14_17_basis
import kim_shahidi_lift_cong.sym16_19_basis
import kim_shahidi_lift_cong.sym18_21_basis

data_dir = os.path.join(os.getenv("HOME"), "ksr_lift_data_dir")

def compute_basis():
    '''
    Compute basis of M_{17, 14}, M_{19, 16}, M_{21, 18}
    with precision 6.
    '''
