# -*- coding: utf-8; mode: sage -*-
import krs_lift_congruence.sym14_17_basis
import krs_lift_congruence.sym16_19_basis
import krs_lift_congruence.sym18_21_basis

modules = [krs_lift_congruence.sym14_17_basis,
           krs_lift_congruence.sym16_19_basis,
           krs_lift_congruence.sym18_21_basis]

def compute_basis():
    '''
    Compute basis of M_{17, 14}, M_{19, 16}, M_{21, 18}
    with precision 6. Result will be saved to utils.data_dir.
    '''
    for m in modules:
        m.calculator.calc_forms_and_save(6, verbose=True)

def compute_lift_and_non_lift():
    """
    Compute KSR lift and non-lift for k = 16, 18 and 20
    with precision 6.
    """
    compute_basis()
    for m in modules:
        m.compute_lift_and_non_lift()

def check_congruence():
    """
    check congruence of Fourier coefficients of the lift
    and the non-lift for k = 16, 18 and 20.
    """
    for m in modules:
        m.check_cong()
    print "Checking done."
