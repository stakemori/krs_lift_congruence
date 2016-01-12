# -*- coding: utf-8; mode: sage -*-

from sage.all import (eisenstein_series_qexp, bernoulli,
                      pi, gamma, prime_range, PowerSeriesRing,
                      QQ, mul, ceil, log, gp, cached_method, ComplexField,
                      Integer, cached_function, fork, floor, ZZ, PolynomialRing)

from sage.misc.all import verbose

import os


def to_bit_prec(prec):
    # a small overestimate of log(10,2)
    LOG_TEN_TWO_PLUS_EPSILON = 3.321928094887363
    return int((prec + 1) * LOG_TEN_TWO_PLUS_EPSILON) + 1


class DokchitserLCalc(object):

    '''
    The original implementation of Dokchitser L calculator
    spawns zombie python processes.
    '''

    def __init__(self,
                 digits=None,
                 conductor=None,
                 gammaV=None,
                 weight=None,
                 sgn=1,
                 Lresidues="automatic",
                 Lpoles=None):
        if Lpoles is None:
            Lpoles = []
        self.data_dct = {"conductor": conductor,
                         "gammaV": gammaV,
                         "weight": weight,
                         "sgn": sgn,
                         "Lresidues": Lresidues,
                         "Lpoles": Lpoles}
        self.digits = digits
        self._coeffs = None
        self._values = {}

    @property
    def _CC(self):
        return ComplexField(to_bit_prec(self.digits))

    def set_coeffs(self, coeffs):
        self._coeffs = coeffs

    def __to_CC(self, s):
        s = s.replace('.E', '.0E').replace(' ', '')
        # return self._CC(sage_eval(, locals={'I':self._CC.gen(0)}))
        return self._CC(s)

    def init_gp_code(self):
        fname = os.path.join(os.getenv("SAGE_ROOT"),
                             "src/ext/pari/dokchitser/computel.gp")
        d = self.data_dct
        return "; ".join(['read("{fname}")'.format(fname=fname),
                          "conductor=%s" % d["conductor"],
                          "gammaV=%s" % d["gammaV"],
                          "weight=%s" % d["weight"],
                          "default(realprecision,%s)" % self.digits,
                          "sgn=%s" % d["sgn"],
                          "Lpoles=%s" % d["Lpoles"],
                          "Lresidues=%s" % d["Lresidues"]])

    def init_gp(self):
        gp.eval(self.init_gp_code())

    @cached_method
    @fork
    def num_of_coeffs(self):
        self.init_gp()
        s = gp.eval("cflength()")
        try:
            return int(s)
        except:
            raise RuntimeError(s)

    def coeffs(self, m=None):
        if self._coeffs is not None:
            return self._coeffs

        if m is None:
            m = Integer(self.num_of_coeffs())
        primes = prime_range(m + 1)
        # Compute Dirichlet series from local Dirichlet sereis.
        res = [1] * m
        for p in primes:
            p_prec = ceil(log(m.n()) / log(p.n()))
            d = self.coeffs_p(p, p_prec)
            ord_ls = [0] * m
            for i in range(1, p_prec + 1):
                p_i = p ** i
                for a in range(1, m // p_i + 1):
                    ord_ls[a * p_i - 1] += 1
            for i, a in enumerate(ord_ls):
                if a != 0:
                    res[i] *= d[a]
        self.set_coeffs(res)
        return res

    def init_l_data_gp_code(self):
        coeffs = self.coeffs()
        return ";".join(['coeff = %s' % coeffs, 'initLdata("coeff[k]")'])

    def init_l_data_gp(self):
        gp.eval(self.init_l_data_gp_code())

    def coeffs_p(self, p, prec):
        '''Should return a list [a_0, a_1, .. , a_prec] that corresponds
        to the p-local Dirichlet series
        a_0 + a_1p^(-s) + a_2 p^(-2s) + ... + a_prec p^(-prec s).
        '''
        raise NotImplementedError

    def set_values(self, s, v):
        self._values[s] = v

    def l_value(self, s):
        s = self._CC(s)
        try:
            return self._values[s]
        except KeyError:
            pass

        @fork
        def computel():
            self.init_gp()
            self.init_l_data_gp()
            return gp.eval('L(%s)' % s)

        z = computel()
        if 'pole' in z:
            print z
            raise ArithmeticError
        elif '***' in z:
            print z
            raise RuntimeError(z)
        elif 'Warning' in z:
            i = z.rfind('\n')
            msg = z[:i].replace('digits', 'decimal digits')
            verbose(msg, level=-1)
            ans = self.__to_CC(z[i + 1:])
            self._values[s] = ans
            return ans
        ans = self.__to_CC(z)
        self.set_values(s, ans)
        return ans


class SymmLCalculator(DokchitserLCalc):

    def __init__(self, f, n,
                 digits=100,
                 gammaV=None,
                 eigenvalue_dct=None,
                 coeffs=None):
        self.n = n
        self.f = f
        k = self.f.weight()
        self.k = k
        weight = n * (k - 1) + 1
        if n % 2 == 0:
            r = n // 2
            gammaV = [- 2 * ((r * (k - 1)) // 2)]
        else:
            r = (n + 1) // 2
            gammaV = []

        for j in range(r):
            gammaV = [-j * (k - 1), -j * (k - 1) + 1] + gammaV

        DokchitserLCalc.__init__(self, digits=digits,
                                 conductor=1,
                                 gammaV=gammaV, weight=weight)

        self._eigenvalue_dct = eigenvalue_dct
        if coeffs is not None:
            self.set_coeffs(coeffs)

    def coeffs_p(self, p, prec):
        raise NotImplementedError

    def set_eigenvalue_dct(self, d):
        self._eigenvalue_dct = d

    def eigenvalue_dct(self):
        if self._eigenvalue_dct is not None:
            return self._eigenvalue_dct
        else:
            num = self.num_of_coeffs()
            f_exp = self.f.qexp(num + 1)
            d = {p: f_exp[p] for p in prime_range(num + 1)}
            self.set_eigenvalue_dct(d)
            return d

    def trancendential_part(self, m):
        return trance_part(self.f, m, n=self.n, digits=self.digits)

    def algebraic_part_numeric(self, m):
        return self.l_value(m) / self.trancendential_part(m)

    def critical_points(self):
        if self.n % 2 == 1:
            raise NotImplementedError
        r = self.n // 2
        k = self.k
        l1 = [
            s for s in range((k - 1) * (r - 1) + 1, (k - 1) * r + 1) if s % 2 == 1]
        l2 = [
            s for s in range((k - 1) * r + 1, (k - 1) * (r + 1) + 1) if s % 2 == 0]
        return l1 + l2


def gamma_inf_part(n, k, s):
    r, rsd = divmod(n, 2)
    if rsd == 0:
        return (pi ** (-ZZ(s) / ZZ(2)) *
                gamma(QQ(s) / QQ(2) - floor(QQ(r * (k - 1)) / QQ(2)))
                * gamma_inf_part(n - 1, k, s))
    else:
        r += 1
        return ((2 * pi) ** (-ZZ(r * s)) *
                mul([gamma(s - j * (k - 1)) for j in range(r)]))


def symml_euler_s_t_pol(j):
    R = PolynomialRing(QQ, names="a, b")
    U = PolynomialRing(QQ, names="s, t, x")
    T = PolynomialRing(R, names="x")
    x = T.gens()[0]
    a, b = R.gens()
    f = mul((1 - a ** i * b ** (j - i) * x) for i in range(j + 1))
    return sum(U(_to_s_t_pol(v)) * U("x") ** k for k, v in f.dict().items())


def _to_s_t_pol(f):
    '''f is a symmetric polynomial of a and b,
    s = a + b, t = a * b, return a polynomial of s and t,
    which is equal to f.
    '''
    R = PolynomialRing(QQ, names="a, b")
    a, b = R.gens()
    _s = a + b
    _t = a * b
    S = PolynomialRing(QQ, names="s, t")
    s, t = S.gens()
    res = 0
    while f != 0:
        i, j = list(sorted(f.dict().keys(), key=lambda x: x[0]))[0]
        i0, j0 = j - i, i
        v = f[(i, j)]
        res += v * s ** i0 * t ** j0
        f -= v * _s ** i0 * _t ** j0
    return res


class SymmNLCalculator(SymmLCalculator):

    def __init__(self, f, n, digits=100, eigenvalue_dct=None, coeffs=None):
        SymmLCalculator.__init__(
            self, f, n, digits=digits, eigenvalue_dct=eigenvalue_dct, coeffs=coeffs)
        self._s_t_pol = symml_euler_s_t_pol(self.n)

    def coeffs_p(self, p, prec):
        a_dct = self.eigenvalue_dct()
        R = PowerSeriesRing(QQ, names="q", default_prec=prec + 1)
        q = R.gens()[0]
        s = a_dct[p]
        t = p ** (self.k - 1)
        _s, _t, _x = self._s_t_pol.parent().gens()
        dnm = self._s_t_pol.subs({_s: s, _t: t, _x: q})
        return (1 / dnm).dict()


class Symm6Calculator(SymmLCalculator):

    def __init__(self, f, digits=100,
                 eigenvalue_dct=None, coeffs=None):
        SymmLCalculator.__init__(self, f, 6,
                                 digits=digits,
                                 eigenvalue_dct=eigenvalue_dct,
                                 coeffs=coeffs)

    def coeffs_p(self, p, prec):
        a_dct = self.eigenvalue_dct()
        R = PowerSeriesRing(QQ, names="q", default_prec=prec + 1)
        q = R.gens()[0]
        s = a_dct[p]
        t = p ** (self.k - 1)
        f1 = s ** 6 - 6 * t * s ** 4 + 9 * t ** 2 * s ** 2 - 2 * t ** 3
        f2 = t * (s ** 4 - 4 * t * s ** 2 + 2 * t ** 2)
        f3 = t ** 2 * (s ** 2 - 2 * t)
        dnm = (mul([1 - f * q + t ** 6 * q ** 2 for f in [f1, f2, f3]]) *
               (1 - t ** 3 * q))
        return (1 / dnm).dict()


class Symm2LCalculator(SymmLCalculator):

    def __init__(self, f, digits=100,
                 eigenvalue_dct=None, coeffs=None):
        SymmLCalculator.__init__(self, f, 2,
                                 digits=digits,
                                 eigenvalue_dct=eigenvalue_dct,
                                 coeffs=coeffs)

    def coeffs_p(self, p, prec):
        a_dct = self.eigenvalue_dct()
        k = self.k
        R = PowerSeriesRing(QQ, names="q", default_prec=prec + 1)
        q = R.gens()[0]
        s = a_dct[p]
        t = p ** (k - 1)
        f1 = s ** 2 - 2 * t
        dnm = (1 - f1 * q + t ** 2 * q ** 2) * (1 - t * q)
        return (1 / dnm).dict()


def coeff_of_f_in_es_prod(q, k):
    esq = eisenstein_series_qexp(q, normalization="constant")
    eskq = eisenstein_series_qexp(k - q, normalization="constant")
    g = esq * eskq - eisenstein_series_qexp(k, normalization="constant")
    a = (-1) ** (q // 2) * 2 ** (k - 3) * \
        bernoulli(q) * bernoulli(k - q) / (q * (k - q))
    return g[1] * a


@cached_function
def petersson_inner_prod(f, digits=100):
    '''
    Assumes that f is a normalized Hecke eigenform of level 1 and
    the Hecke field is eqaul to Q.
    '''
    k = f.weight()
    N = f.level()
    w = f.atkin_lehner_eigenvalue()
    e = (-1) ** (k // 2) * w
    l_calc = DokchitserLCalc(conductor=N, gammaV=[0, 1],
                             weight=k, sgn=e, digits=digits + 1,
                             Lpoles=[])

    num = l_calc.num_of_coeffs()
    f_qexp = f.qexp(prec=num + 2)
    l_calc.set_coeffs([f_qexp[a + 1] for a in range(num + 1)])

    q = [a for a in range(k // 2 + 2, k - 3) if a % 2 == 0][0]

    def l_func(s):
        return (2 * pi.n(digits=digits)) ** (-s) * gamma(s) * l_calc.l_value(s)

    return l_func(q) * l_func(k - 1) / coeff_of_f_in_es_prod(q, k)


def trance_part(f, m, n=6, digits=100):
    '''
    Assumes n is even.
    '''
    if n % 2 == 1:
        raise NotImplementedError
    pinprd = petersson_inner_prod(f, digits=digits)
    k = f.weight()
    r = n // 2
    a = pinprd ** (r * (r + 1) / 2)
    if m % 2 == 1:
        e = r * m - (r * (r - 1) / 2) * (k - 1)
    else:
        e = (r + 1) * m - r * (r + 1) / 2 * (k - 1)
    b = ((2 * pi) ** e).n(digits=digits)
    return a * b

# def f(x, k):
#     return gamma(x) * gamma(x - k + 1) * gamma(x - 2*k + 2)

# def g(x, k):
#     return f(x, k)/f((k - 1)*6 + 1 - x, k)


def complete_l_value(l_alg, f):
    k = f.weight()
    m = 3 * k - 2
    digits = 300
    trns = trance_part(f, m, digits=digits)
    cl_val = l_alg * trns * (gamma_inf_part(6, k, m)).n(digits=digits)
    return cl_val


def oposite_value_alg_part_using_feq(l_alg, f):
    k = f.weight()
    m = 3 * k - 2
    mop = 3 * k - 3
    digits = 300
    trns = trance_part(f, m, digits=digits)
    trns_op = trance_part(f, mop, digits=digits)
    cl_val = l_alg * trns * (gamma_inf_part(6, k, m)).n(digits=digits)
    l_val_op = cl_val / (gamma_inf_part(6, k, mop)).n(digits=digits)
    return l_val_op / trns_op
