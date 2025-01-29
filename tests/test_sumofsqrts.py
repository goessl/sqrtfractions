from sympy import factorint
from math import prod
from fractions import Fraction
from random import randint
from sqrtfractions import SumOfSqrts




from math import isclose as iscl

def isclose(a, *b, rel_tol=1e-09, abs_tol=0.0):
    return all(iscl(a, bi, rel_tol=rel_tol, abs_tol=abs_tol) for bi in b)



def test_sumofsqrts():
    #creation & evaluation
    for _ in range(100):
        n = randint(-100, +100)
        a = SumOfSqrts(n)
        assert a==n and int(a)==n and isclose(float(a), n)
    
    
    #comparison
    for _ in range(100):
        a = SumOfSqrts.random()
        assert isclose(float(abs(a)), abs(float(a)))
    
    for _ in range(100):
        a, b = SumOfSqrts.random(), SumOfSqrts.random()
        assert (a<b) == (float(a)<float(b))



def test_sumofsqrts_neg():
    for _ in range(1000):
        a = SumOfSqrts.random()
        actual = -float(a)
        prediction = -a
        assert isclose(actual, float(prediction)) and prediction.validate()

def test_sumofsqrts_add():
    for _ in range(1000):
        a, b = SumOfSqrts.random(), SumOfSqrts.random()
        actual = float(a) + float(b)
        prediction = a + b
        assert isclose(actual, float(prediction)) and prediction.validate()
        
        a, b = SumOfSqrts.random(), randint(-20, +20)
        actual = float(a) + b
        prediction0 = a + b
        prediction1 = b + a
        assert isclose(actual, float(prediction0), float(prediction1)) \
                and prediction0.validate() and prediction1.validate()

def test_sumofsqrts_sub():
    for _ in range(1000):
        a, b = SumOfSqrts.random(), SumOfSqrts.random()
        actual = float(a) - float(b)
        prediction = a - b
        assert isclose(actual, float(prediction)) and prediction.validate()
        
        a, b = SumOfSqrts.random(), randint(-20, +20)
        actual = float(a) - b
        prediction0 = a - b
        prediction1 = -(b - a)
        assert isclose(actual, prediction0, prediction1) \
                and prediction0.validate() and prediction1.validate()

def test_sumofsqrts_mul():
    for _ in range(1000):
        a, b = SumOfSqrts.random(), SumOfSqrts.random()
        actual = float(a) * float(b)
        prediction = a * b
        assert isclose(actual, float(prediction)) and prediction.validate()
        
        a, b = SumOfSqrts.random(4), randint(-20, +20)
        actual = float(a) * b
        prediction0 = a * b
        prediction1 = b * a
        assert isclose(actual, prediction0, prediction1) \
                and prediction0.validate() and prediction1.validate()

def test_sumofsqrts_pow():
    for _ in range(10):
        a = SumOfSqrts.random(5)
        n = randint(0, +3)
        actual = float(a)**n
        prediction = a**n
        assert isclose(actual, float(prediction)) and prediction.validate()
