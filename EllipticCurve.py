from functools import cache


class EllipticCurve:

    def __init__(self, a: int, b: int, field: int):
        self.a = a
        self.b = b
        self.field = field

    def weierstrass(self, x: int) -> int:
        """Returns the y for the weierstrass equation"""
        return (x ** 3 + self.a * x + self.b) % self.field

    def is_on_curve(self, point) -> bool:
        """Returns true if a point is on the curve"""
        x, y = get_xy(point)
        return y ** 2 % self.field == self.weierstrass(x)

    def __repr__(self):
        return f"EllipticCurve(field={self.field}, a={self.a}, b={self.b})"

    def __str__(self) -> str:
        return f"y^2 = x^3 + {self.a}x + {self.b} in F{self.field}"

    @property
    def latex_str(self) -> str:
        """ Latex repr, should be printed"""
        return f"y^2 = x^3 + {self.a}x + {self.b} in \\mathbb{{F}}_{self.field}"

    def multiplicative_inverse(self, divisor) -> int:
        """ Return multiplicative inverse with scalar and mod the field of this curve """
        return multiplicative_inverse(divisor, self.field)

    def __contains__(self, point) -> bool:
        """ Returns True if the point is ~~in~~ on the curve """
        return self.is_on_curve(point)

    def __eq__(self, curve) -> bool:
        """Returns True if a curve is the same"""
        return self.a == curve.a and self.b == curve.b and self.field == curve.field

    def __hash__(self):
        return hash(self.a * self.b * self.field)


class EllipticPoint:

    def __init__(self, x: int, y: int, curve: EllipticCurve):
        self.x = x
        self.y = y
        self.curve = curve

    @property
    def on_curve(self) -> bool:
        """Returns True if the point is on its curve"""
        return self.curve.is_on_curve(self)

    def slope(self, point):
        """ Calculate the slope between this point and another point. This is lambda or λ"""
        x, y = get_xy(point)
        field = self.curve.field
        if x != self.x:
            dividend = self.y - y
            divisor = self.x - x
        else:
            dividend = 3 * x ** 2 + self.curve.a
            divisor = (2 * y)
        inverse = multiplicative_inverse(divisor, field)
        return dividend * inverse % field

    def set_xy(self, x: int, y: int):
        """Set x and y attributes of the point in one go"""
        self.x = x
        self.y = y

    @property
    def inverse(self) -> 'EllipticPoint':
        """Return the inverse of the point. You can also use the ~ operator"""
        return self.__invert__()

    def __invert__(self) -> 'EllipticPoint':
        return EllipticPoint(x=self.x, y=(-self.y) % self.curve.field, curve=self.curve)

    def __add__(self, point) -> 'EllipticPoint':
        """Add points together"""
        x, y = get_xy(point)
        if isinstance(point, InfinityPoint) or x == float('inf') and y == float('inf'):  # Adding infinity return copy
            return EllipticPoint(x=self.x, y=self.y, curve=self.curve)
        if self == ~EllipticPoint(x, y, self.curve):  # Return point at infinity if adding inverse
            return InfinityPoint(curve=self.curve)

        slope = self.slope(point)
        x_result = (slope ** 2 - self.x - x) % self.curve.field
        y_result = (slope * (self.x - x_result) - self.y) % self.curve.field
        return EllipticPoint(x=x_result, y=y_result, curve=self.curve)

    def __iadd__(self, point):
        """Add points together and replace self"""
        result = self + point
        x, y = get_xy(result)
        self.set_xy(x, y)
        return result

    def __mul__(self, n: int):
        """Calculate scalar multiplication"""
        if n == 1:
            return self
        result = self + self
        while n - 2 > 0:
            result += self
            n -= 1
        if n == 1:
            result += self
        return result

    def __imul__(self, n):
        result = self.__mul__(n)
        x, y = get_xy(result)
        self.set_xy(x, y)
        return result

    def __repr__(self):
        return f"EllipticPoint(x={self.x}, y={self.y}, curve={self.curve.__repr__()})"

    def __str__(self):
        return f"{self.x, self.y}"

    @property
    def latex_str(self):
        return f"{self.x, self.y} \\in {self.curve.latex_str}"

    def __getitem__(self, index: int):
        if index == 0:
            return self.x
        elif index == 1:
            return self.y
        raise IndexError(f"{self.__repr__()} only supports index 0 and 1")

    def __iter__(self):
        return self.x, self.y

    def __eq__(self, point):
        return isinstance(point, EllipticPoint) \
               and self.x == point.x \
               and self.y == point.y \
               and self.curve == point.curve

    @property
    @cache
    def cycle(self) -> ['EllipticPoint']:
        """ Generate the cycle of the point. Keep adding the point to itself until the point at infinity is reached"""
        points = [self, self + self]
        if self not in self.curve:
            import warnings
            s = f"{self} is not on {self.curve} so .cycle might create an infinite loop!"
            warnings.warn(s)
        while (points[-1]) != InfinityPoint(self.curve):
            points.append(points[-1] + self)
        return points

    @property
    def order(self) -> int:
        """ The number of elements in the cycle of this point. """
        return len(self.cycle)

    def __hash__(self):
        return hash(self.x * self.y * hash(self.curve))


class InfinityPoint(EllipticPoint):
    def __init__(self, curve):
        super(InfinityPoint, self).__init__(x=float('inf'), y=float('inf'), curve=curve)

    def __add__(self, point):
        """ The infinite point is the identity element """
        x, y = get_xy(point)
        return EllipticPoint(x=x, y=y, curve=self.curve)


def get_xy(point):
    """ Return x and y from a point. This way we allow anything that satisfies below """
    return point[0], point[1]


def is_prime(p) -> bool:
    """ Basic prime algorithm returns True if p is prime """
    if p <= 3:
        return p > 1
    if p % 2 == 0 or p % 3 == 0:
        return False
    i = 5
    while i ** 2 <= p:
        if p % i == 0 or p % (i + 2) == 2:
            return False
    return True


def extended_gcd(a, b): 
    if a == 0:
        return b, 0, 1

    x, y, gcd = extended_gcd(b % a, a)

    return y - (b // a) * x, x, gcd


def multiplicative_inverse(n, mod):
    """ Give the multiplicative inverse """
    if is_prime(mod):  # If n is prime we can use Fermat.
        return n ** (mod - 2)

    x, y, gcd = extended_gcd(n, mod)
    if gcd == 1:
        return x
    else:
        raise ValueError(f"No inverse of {n}  mod {mod}!")


__author__ = "Quinten Cabo"


if __name__ == '__main__':
    # How to use:
    C13 = EllipticCurve(a=5, b=4, field=13)
    G = EllipticPoint(1, 7, C13)
    B = EllipticPoint(4, 7, C13)
    A = G * 3
    print(f"[3]{G} = {A} on {C13}")

    print(f'[3]{B} = {B * 3} on {C13}')
    print(f'[4]{B} = {B + B + B + B} on {C13}')

    print('Cycle of {G}:','\n'.join(map(str,G.cycle)))
    
    print("")

    C = EllipticCurve(11, 18, 23)
    D = EllipticCurve(11, 2, 22)
    P = EllipticPoint(15, 4, C)
    Q = EllipticPoint(14, 8, C)


    print(f"Curve: {C}")
    print(f"{P} in {C} = {P in C}")
    print(f"{Q} in {C} = {Q in C}")
    print(f"Slope: {P} and {Q} = {P.slope(Q)}")
    print(f"P + Q = {P} + {Q} = {P + Q}")
    print()
    print(f"P + P = {P + P}")  # (22,11)
    print(f"[3]P = P + P + P = {P * 3} = {P + P + P}")
    print()
    print(f"[1]P = [1]{P} = {P * 1}")
    print(f"[2]P = [2]{P} = {P * 2}")
    print("....")
    print(f"[30]P = [30]{P} = {P * 30}")
    print(f"[31]P = [31]{P} = {P * 31}")
    print(f"[32]P = [32]{P} = {P * 32}")
    print()
    print(f"Order of {P} is {P.order}")
