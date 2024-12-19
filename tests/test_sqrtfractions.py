from sympy import factorint
from math import prod
from fractions import Fraction
from random import randint
from sqrtfractions import SqrtFraction


def test_sqrtfractions():
    from math import isclose as iscl
    
    def isclose(a, *b, rel_tol=1e-09, abs_tol=0.0):
        return all(iscl(a, bi, rel_tol=rel_tol, abs_tol=abs_tol) for bi in b)
    
    
    #creation & evaluation
    for n in range(1, 1000):
        f = factorint(n)
        s = prod(p**(e//2) for p, e in f.items())
        k = prod(p**(e%2) for p, e in f.items())
        assert s**2 * k == n
    
    for _ in range(100):
        n = randint(-100, +100)
        a = SqrtFraction(n)
        assert a==n and int(a)==n and isclose(float(a), n)
        
        n, d = randint(-100, +100), randint(-100, +100-1) or +100
        a = SqrtFraction(n, d)
        assert a.as_fraction()==Fraction(n, d) and isclose(float(a), n/d)
    
    
    #comparison
    for _ in range(100):
        a = SqrtFraction.random()
        assert isclose(float(abs(a)), abs(float(a)))
    
    for _ in range(100):
        a, b = SqrtFraction.random(), SqrtFraction.random()
        assert (a<b) == (float(a)<float(b))
    
    
    #add
    for _ in range(1000):
        a, b = SqrtFraction.random(), SqrtFraction.random()
        assert isclose(float(a)+float(b), float(a+b))
        
        a, b = SqrtFraction.random(), randint(-20, +20)
        assert isclose(float(a)+b, float(a+b), float(b+a))
    
    #sub
    for _ in range(1000):
        a, b = SqrtFraction.random(), SqrtFraction.random()
        assert isclose(float(a)-float(b), float(a-b))
        
        a, b = SqrtFraction.random(), randint(-20, +20)
        assert isclose(float(a)-b, float(a-b), -float(b-a))
    
    #mul
    for _ in range(1000):
        a, b = SqrtFraction.random(), SqrtFraction.random()
        assert isclose(float(a)*float(b), float(a*b))
        
        a, b = SqrtFraction.random(4), randint(-20, +20)
        assert isclose(float(a)*b, float(a*b), float(b*a))
    
    #invert
    for _ in range(10):
        a = SqrtFraction.random(5)
        assert a / a == 1
    
    #div
    for _ in range(10):
        a, b = SqrtFraction.random(5), SqrtFraction.random(5)
        assert isclose(float(a)/float(b), float(a/b), rel_tol=1e-5)
        
        a, b = SqrtFraction.random(5), randint(-20, +20-1) or +20
        assert isclose(float(a)/float(b), float(a/b))
        
        a, b = SqrtFraction.random(5), randint(-20, +20)
        assert isclose(float(b)/float(a), float(b/a), rel_tol=1e-5)
    
    #pow
    for _ in range(10):
        a = SqrtFraction.random(5)
        n = randint(-3, +3)
        assert isclose(float(a**n), float(a)**n)

