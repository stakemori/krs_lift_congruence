# -*- coding: utf-8; mode: sage -*-
from sage.all import QQ
from sage.rings.number_field.number_field import NumberField_absolute
def modulo_p(alpha, a, p):
    return sum([b * a**i for i, b in enumerate(alpha.list())])%p

def is_p_integral(alpha, p):
    K = alpha.parent()
    if K == QQ:
        return alpha.denominator() % p != 0
    elif isinstance(K, NumberField_absolute):
        return all((a.denominator() % p != 0 for a in alpha.list()))
