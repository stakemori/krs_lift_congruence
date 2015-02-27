# -*- coding: utf-8; mode: sage -*-
from sage.all import QQ, PolynomialRing, FiniteField
from sage.rings.number_field.number_field import NumberField_absolute
import os

def modulo_p(alpha, a, p):
    return sum([b * a**i for i, b in enumerate(alpha.list())])%p

def is_p_integral(alpha, p):
    K = alpha.parent()
    if K == QQ:
        return alpha.denominator() % p != 0
    elif isinstance(K, NumberField_absolute):
        return all((a.denominator() % p != 0 for a in alpha.list()))

def check_cong(p, t2_eigenvalue, lift, non_lift, space):

    def mod_p(alpha):
        return modulo_p(alpha, t2_eigenvalue, p)

    v_lift = space._to_vector(lift)
    v_non_lift = space._to_vector(non_lift)
    # check p-integral
    assert all((is_p_integral(a, p) for a in v_lift))
    assert all((is_p_integral(a, p) for a in v_non_lift))
    K = non_lift.base_ring
    R = PolynomialRing(FiniteField(p), names="x")
    pl_modp = R(K.polynomial())
    # p is unramified.
    assert all((b == 1 for a, b in pl_modp().factor()))
    # check congruence
    v = v_lift - v_non_lift
    assert all(mod_p(b)%p == 0 for b in v)

data_dir = os.path.join(os.getenv("HOME"), "ksr_lift_data_dir")
