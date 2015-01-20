# -*- coding: utf-8; mode: sage -*-

def modulo_p(alpha, a, p):
    return sum([b * a**i for i, b in enumerate(alpha.list())])%p
