# -*- coding: utf-8; mode: sage -*-
from sage.all import QQ, PolynomialRing, FiniteField, NumberField
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

def _lift_chr_ply():
    x = PolynomialRing(QQ, names="x").gens()[0]
    d = {16: x + 4078080,
         18: x + 8785920,
         20: x + 383331840}
    return d


def fname(f):
    return os.path.join(data_dir, f)


def _lift_name(k):
    return fname("lift_%s_%s_prec6.sobj"%(str(k+1), str(k-2)))


def _nonlift_name(k):
    return fname("nonlift_%s_%s_prec6.sobj"%(str(k+1), str(k-2)))


def compute_lift_and_non_lift(space, prec):
    k = space.wt - 1
    lift_name = _lift_name(k)
    nonlift_name = _nonlift_name(k)
    lift_chr_ply = _lift_chr_ply()[k]
    lift_ev = -lift_chr_ply.constant_coefficient()
    R = PolynomialRing(QQ, names="x")
    K = NumberField(R((space.hecke_charpoly(2)/lift_chr_ply)), names="a")
    nonlift_ev = K.gens()[0]
    lift = space.eigenform_with_eigenvalue_t2(lift_ev)
    lift.save_as_binary(lift_name)
    nonlift = space.eigenform_with_eigenvalue_t2(nonlift_ev)
    nonlift.save_as_binary(nonlift_name)
