from math import sqrt, prod, gcd, lcm, sumprod
from collections import defaultdict
from fractions import Fraction
from operator import mul
import sympy as sp
from sympy.ntheory import factorint, primefactors
from itertools import product, islice, repeat, chain
from functools import total_ordering, reduce
from random import randint



@total_ordering
class SqrtFraction:
    r"""Class to handle numbers in the form $\sum_in_i\sqrt{r_i}/d$.
    
    The radicands are integers greater or equal 1,
    the factors in front of them are integers,
    the denominator is a integer greater or equal 1.
    This class is immutable.
    """
    
    #creation
    def __init__(self, n={}, d=1):
        """Create a new `SqrtFraction`.
        
        The numerator `n` is expected to be a dictionary with `ri:ni` as
        keys:values, where the radicands `ri` must be integers >=1
        and the factors `ni` must be integers.
        `n` can also be an integer for an integer fraction.
        The denominator `d` must be a non-zero integer.
        No arguments defaults to 0.
        Numerator and denominator will be simplified to shortest terms.
        """
        if isinstance(n, int):
            n = {1:n}
        elif isinstance(n, dict):
            if not (all(isinstance(k, int) for k in n.keys())
                    and all(isinstance(v, int) for v in n.values())):
                raise TypeError('Radicands and factors must be integers.')
            if not all(k>0 for k in n.keys()):
                raise ValueError('Radicands must be greater than zero.')
        else:
            raise TypeError('Numerator must be an integer or a dictionary.')
        
        if isinstance(d, int):
            if not d:
                raise ValueError('Denominator must be non-zero.')
        else:
            raise TypeError('Denominator must be an integer.')
        
        #simplify numerator
        _n = defaultdict(int)
        for k, v in n.items():
            f = factorint(k)
            s = prod(p**(e//2) for p, e in f.items())
            k = prod(p**(e%2) for p, e in f.items())
            _n[k] += s * v
        n = {k:v for k, v in _n.items() if v}
        #reduce fraction, make denominator positive
        cd = gcd(*n.values(), d) if d>0 else -gcd(*n.values(), d)
        n = {k:v//cd for k, v in n.items()}
        d //= cd
        
        self.n, self.d = dict(sorted(n.items())), d
    
    @staticmethod
    def random(N=10, precision=20):
        r"""Return a random `SqrtFraction.
        
        The factors from $\sqrt{1}$ up to $\sqrt{N}$ (incl.)
        will be initialised with random integers `ni`
        such that `-precision <= ni <= +precision`
        and with a random denominator `d` such that `1 <= d <= precision`.
        """
        n = {n:randint(-precision, +precision) for n in range(1, N+1)}
        d = randint(1, precision)
        return SqrtFraction(n, d)
    
    
    #evaluation
    def __float__(self):
        """Calculate the `float` approximation.
        
        The `float` value might overflow.
        """
        #dicts aren't hashable, so __float__ can't be made a cached_property
        if not hasattr(self, '_float'):
            self._float = sumprod(self.values(), map(sqrt, self.keys())) \
                    / self.d
            #empty sumprod=int(0) nevertheless becomes float due to division
        return self._float
    
    def is_integer(self):
        """Return if this is an integer.
        
        If `True`, `int()` will return an integer,
        if `False`, `int()` will raise a `ValueError`.
        """
        return set(self.keys())<={1} and self.d==1
    
    def is_fraction(self):
        """Return this is an integer fraction.
        
        If `True`, `as_fraction()` will return a fraction,
        if `False`, `as_fraction()` will raise a `ValueError`.
        """
        return set(self.keys()) <= {1}
    
    def __int__(self):
        """Return the integer value.
        
        If this doesn't represent an integer value a ValueError is raised.
        Can be checked beforehand with `is_integer()`.
        """
        if self.is_integer():
            return self.n.get(1, 0)
        else:
            raise ValueError('doesn\'t represent an integer')
    
    def as_fraction(self):
        """Return if this as an integer fraction.
        
        If this doesn't represent an integer fraction a ValueError is raised.
        Can be checked beforehand with `is_fraction()`.
        """
        if self.is_fraction():
            return Fraction(self.n.get(1, 0), self.d)
        else:
            raise ValueError('doesn\'t represent a fraction')
    
    def __bool__(self):
        """Return if this is non-zero."""
        #leave this to avoid `if obj` call len(obj) instead of int(obj)
        #https://docs.python.org/3/reference/datamodel.html#object.__bool__
        #this is a numeric class, not a container
        try:
            return bool(int(self))
        except:
            return True
    
    def sympify(self):
        """Return this as a `sympy`-expression."""
        return sumprod(self.values(), map(sp.sqrt, self.keys())) \
                / sp.Integer(self.d) #cast denominator in case of empty sum
    
    
    #collection
    def __len__(self):
        """Return the number of summands in the numerator."""
        return len(self.n)
    
    def keys(self):
        """Return the radicands `r_i`."""
        return self.n.keys()
    
    def values(self):
        """Return the factors `n_i` infront of the square roots."""
        return self.n.values()
    
    def items(self):
        """Return the radicands `r_i` with factors `n_i` as tuples."""
        return self.n.items()
    
    
    #ordering
    def __eq__(self, other):
        """Return if this equals another `SqrtFraction`, `int` or `Fraction`."""
        if isinstance(other, SqrtFraction):
            return self.n == other.n and self.d == other.d
        elif isinstance(other, int) or isinstance(other, Fraction):
            try: #Fraction handles int comparison
                return self.as_fraction() == other
            except:
                return False
        else:
            raise TypeError
    
    def __abs__(self):
        """Return the absolute value as a `SqrtFraction`."""
        return self if self>=0 else -self
    
    def __lt__(self, other):
        """Return if this is less than another `SqrtFraction` or `integer`."""
        #https://math.stackexchange.com/a/1076510
        if isinstance(other, SqrtFraction) or isinstance(other, int):
            l = self - other
            
            while not l.is_fraction():
                p = max(chain(*(primefactors(k) for k in l.keys())))
                
                #drop denominator
                r = SqrtFraction({k//p:-v for k, v in l.items() if k%p==0})
                l = SqrtFraction({k:v for k, v in l.items() if k%p!=0})
                #https://math.stackexchange.com/a/2347212
                l, r = l*abs(l), r*abs(r)*p
                
                l -= r
            
            return l.as_fraction() < 0
        else:
            raise TypeError
    
    
    #printing
    def __repr__(self):
        """Return a Unicode representation."""
        n = [f'{n:+d}{chr(0x221A)}{r}' for r, n in self.items()]
        if len(n) <= 1: #no parentheses needed
            return (n[0] if n else '0') + f'/{self.d}'
        else:
            return '('+''.join(n)+')' + f'/{self.d}'
    
    def _repr_latex_(self):
        """Return a Latex representation."""
        n = [f'{n:+d}\\sqrt{{{r}}}' for r, n in self.items()]
        return '$\\frac{' + (''.join(n) if n else '0') + f'}}{{{self.d}}}$'
    
    
    #arithmetic
    #Implement the main arithmetic operations (+, *) with type checks
    #Inverse arithmetic (-, /) is completely reused
    def __add__(self, other):
        """Return the sum with an other `SqrtFraction` or `int`."""
        if isinstance(other, SqrtFraction):
            d = lcm(self.d, other.d)
            fs, fo = d//self.d, d//other.d
            n = defaultdict(int)
            for k, v in self.items():
                n[k] += v * fs
            for k, v in other.items():
                n[k] += v * fo
            return SqrtFraction(n, d)
        elif isinstance(other, int):
            n = defaultdict(int, self.n)
            n[1] += other * self.d
            return SqrtFraction(n, self.d)
        else:
            raise TypeError('can only add SqrtFraction or int'
                    + f' (not "{type(other).__name__}") to SqrtFraction')
    __radd__ = __add__
    
    def __neg__(self):
        """Return the additive negation."""
        return SqrtFraction(self.n, -self.d)
    
    def __sub__(self, other):
        """Return the difference with an other `SqrtFraction` or `int`."""
        return self + (-other)
    
    def __rsub__(self, other):
        """Return the difference from an other `SqrtFraction` or `int`."""
        return (-self) + other
    
    def __mul__(self, other):
        """Return the product with an other `SqrtFraction` or `int`."""
        if isinstance(other, SqrtFraction):
            n = defaultdict(int)
            for ki, vi in self.items():
                for kj, vj in other.items():
                    n[ki*kj] += vi * vj
            return SqrtFraction(n, self.d*other.d)
        elif isinstance(other, int):
            return SqrtFraction({k:v*other for k, v in self.items()}, self.d)
        else:
            raise TypeError('can only multiply SqrtFraction or int'
                    + f' (not "{type(other).__name__}") with SqrtFraction')
    __rmul__ = __mul__
    
    def __invert__(self):
        """Return the multiplicative reciprocal."""
        #https://www.youtube.com/watch?v=SjP6Mer0aL8
        #https://en.wikipedia.org/wiki/Rationalisation_(mathematics)
        if self == 0: #repeat=-1 would also raise error, but this is clearer
            raise ZeroDivisionError
        if len(self) == 1:
            r, n = next(iter(self.items()))
            return SqrtFraction({r:self.d}, n*r)
        r = 1
        for p in islice(product((+1, -1), repeat=len(self)-1), 1, None):
            r *= SqrtFraction(
                    {k:s*v for (k, v), s in zip(self.items(), (+1,)+p)})
        return r * self.d / int(r * SqrtFraction(self.n))
    
    def __truediv__(self, other):
        """Return the quotient with an other `SqrtFraction` or `int`."""
        if other == 0:
            raise ZeroDivisionError
        if isinstance(other, SqrtFraction):
            return self * ~other
        elif isinstance(other, int):
            return SqrtFraction(self.n, self.d*other)
        else:
            raise TypeError('can only divide SqrtFraction by'
                    + f' SqrtFraction or int (not "{type(other).__name__}")')
    
    def __rtruediv__(self, other):
        """Return an other `SqrtFraction` or `int` divided by this."""
        return (~self) * other
    
    def __pow__(self, other):
        """Return this `SqrtFraction` raised to some integer power.
        
        The exponent may be negative.
        """
        #repeat does typecheck for int
        if other >= 0:
            return reduce(mul, repeat(self, other), SqrtFraction(1))
        else:
            return (~self)**(-other)



if __name__ == '__main__':
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
