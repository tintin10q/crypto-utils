from functools import cache
from itertools import chain


class Step:
    """ A linked list of steps.
    The idea is that the text statements should be executable in python and that latex is the latex version of the step.
    Just print a step, and it will also print the steps after it.
    """

    def __init__(self, text: str = "", latex_text: str = "No latex version"):
        self.next_step = None
        self.text = text
        self.latex_text = latex_text

    def add_step(self, text: str, latex_text=None):
        if latex_text is None:
            latex_text = text

        # If my text is not set, set it
        if not self.text:
            self.text = text
            self.latex_text = latex_text
            return

        # If we are last in the chain set as our next step
        if self.next_step is None:
            self.next_step = Step(text, latex_text)
            return

        # Else set as the next step of the next step
        self.next_step.add_step(text, latex_text)

    def clear(self):
        """ Reset the steps in the linked list """
        if self.next_step:
            self.next_step.clear()
        del self.next_step
        self.next_step = None
        self.text = ""
        self.latex_text = ""

    def add_linebreak(self):
        self.add_step('\n', '\n')

    def __len__(self):
        if self.next_step is None:
            return 1
        return 1 + len(self.next_step)

    def __iter__(self):
        if self.next_step is None:
            raise StopIteration
        return chain([self.next_step], self.next_step.__iter__)

    def __str__(self):
        if self.next_step is None:
            return self.text
        return self.text + '\n' + self.next_step.__str__()

    @property
    def latex(self):
        """Print Full Latex output"""
        if self.next_step is None:
            return self.latex_text
        return self.latex_text + '\n' + self.next_step.latex

    @property
    def last(self):
        """Returns last step in the chain"""
        if self.next_step is None:
            return self
        return self.last

    def __iadd__(self, steps):
        """ Add a bunch of steps to the current step """
        if not isinstance(steps, Step):
            raise ValueError(f"You can only add steps to steps. You can't add {type(steps)} to Step.")
        self.last.next_step = steps
        return self

    def __add__(self, other):
        s = Step()
        s += self
        s += other
        return s


class EllipticCurve:

    def __init__(self, a: int, b: int, field: int):
        self.a = a
        self.b = b
        self.field = field

    def weierstrass(self, point: 'EllipticPoint') -> int:
        """Returns the y for the weierstrass equation"""
        x, y = get_xy(point)
        point.steps.add_step('Calculate weierstrass:', 'Calculate weierstrass:')
        point.steps.add_step('y**2 = (x ** 3 + a * x + b) % Fp', '$y^2 = x^3 + a \\times x + b \\mod \\mathbb{F}_p$')
        point.steps.add_step(f'y**2 = ({x} ** 3 + {self.a} * {x} + {self.b}) % {self.field}',
                             f'$y^2 = {x}^3 + {self.a} \\times {x} + {self.b} \\mod {self.field}$')
        point.steps.add_step(f'y**2 = ({x ** 3} + {self.a * x} + {self.b}) % {self.field}',
                             f'$y^2 = {x ** 3} + {self.a * x} + {self.b} \\mod {self.field}$')
        point.steps.add_step(f'y**2 = ({x ** 3 + self.a * x} + {self.b}) % {self.field}',
                             f'$y^2 = {x ** 3 + self.a * x} + {self.b} \\mod {self.field}$')
        point.steps.add_step(f'y**2 = {x ** 3 + self.a * x + self.b} % {self.field}',
                             f'$y^2 = {x ** 3 + self.a * x + self.b} \\mod {self.field}$')
        point.steps.add_step(f'y**2 = {(x ** 3 + self.a * x + self.b) % self.field}',
                             f'$y^2 = {(x ** 3 + self.a * x + self.b) % self.field}$')
        return (x ** 3 + self.a * x + self.b) % self.field

    def is_on_curve(self, point: 'EllipticPoint') -> bool:
        """Returns true if a point is on the curve"""
        x, y = get_xy(point)
        field = self.field
        point.steps.add_step(f'Is {point} on {self}?', f'Is ${point}$ on ${self.latex}$?')
        if point == InfinityPoint(point.curve):
            point.steps.add_step(f'Yes {point} on {self} as it is the point at infinity.',
                                 f'Yes, ${point}$ is on ${self.latex}$ as the point is $\mathcal{{O}}$?')
            return True
        weierstrass = self.weierstrass(point)
        point.steps.add_step(f'Now check if {y}^2 = {weierstrass}:', f'Now check if ${y}^2 = {weierstrass}$:')
        point.steps.add_step(f'{y}**2 % {field} = {weierstrass}', f'${y}^2 \\mod {field} = {weierstrass}$')
        point.steps.add_step(f'{y ** 2} % {field} = {weierstrass}', f'${y ** 2} \\mod {field} = {weierstrass}$')
        point.steps.add_step(f'{y ** 2 % field} = {weierstrass}', f'${y ** 2 % field} = {weierstrass}$')
        if y ** 2 % field == weierstrass:
            point.steps.add_step(f'{point} is on {self}!', f'Point is on the curve! ${point} \\in {self.latex}$')
        else:
            point.steps.add_step(f'{point} is not on {self}!', f'${point} \\notin {self.latex}$')

        return y ** 2 % field == weierstrass

    def __repr__(self):
        return f"EllipticCurve(field={self.field}, a={self.a}, b={self.b})"

    def __str__(self) -> str:
        return f"y^2 = x^3 + {self.a}x + {self.b} in F{self.field}"

    @property
    def latex(self) -> str:
        """ Latex repr, should be printed"""
        return f"y^2 = x^3 + {self.a}x + {self.b} \\in \\mathbb{{F}}_{{{self.field}}}"

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

    def __init__(self, x: int, y: int, curve: EllipticCurve, steps=None):
        self.x = x
        self.y = y
        self.curve = curve
        if steps is None:
            self.steps = Step()
        else:
            self.steps = steps

    @property
    def on_curve(self) -> bool:
        """Returns True if the point is on its curve"""
        return self.curve.is_on_curve(self)

    def slope(self, point, steps=None):
        """ Calculate the slope between this point and another point. This is lambda or λ"""
        x, y = get_xy(point)
        field = self.curve.field

        if steps is None:
            if x != self.x:
                dividend = self.y - y
                divisor = self.x - x
            else:
                dividend = (3 * x ** 2 + self.curve.a)
                divisor = (2 * y) % field
            inverse = multiplicative_inverse(divisor, field, steps=self.steps)
            the_slope = dividend * inverse % field
            return the_slope

        steps.add_step(f'We want to find to slope between P = {self} and Q = {point} on {self.curve}',
                       f'We want to find the slope ($\lambda$) between $P={self}$ and $Q={point}$ on ${self.curve.latex}$')
        if x != self.x:
            steps.add_step(f'x_y != x_q ({self.x} != {x}), slope = (y_p-y_q)\(x_p-q_x) % {field}',
                           f'$x_p \\neq x_q~~({self.x} \\neq {x}), \\lambda = \\frac{{y_p-y_q}}{{x_p-x_q}} \\mod {field}$')
            steps.add_step(f'slope = ({self.y}-{y})/({self.x}-{x}) % {field}',
                           f'$\\lambda = \\frac{{{self.y}-{y}}}{{{self.x}-{x}}} \\mod {field}$')
            steps.add_step(f'slope = ({self.y - y})/({self.x - x}) % {field}',
                           f'$\\lambda = \\frac{{{self.y - y}}}{{{self.x - x}}} \\mod {field}$')
            dividend = self.y - y
            divisor = self.x - x
        else:
            steps.add_step(f'x_p == x_q ({self.x} == {x}) so slope = (3 * x_p**2 + a)/(2 * y_p) % {field}',
                           f'$x_p == x_q~~({self.x} == {x})$ so: $\\lambda = \\frac{{3 * (x_p)^2 + a}}{{2 \\times y_p}} \\mod {field}$')
            steps.add_step(f'slope = (3 * {self.x}**2 + {self.curve.a})/(2 * {self.y}) % {field}',
                           f'$\\lambda = \\frac{{3 * {self.x}^2 + {self.curve.a}}}{{2 \\times {self.y}}} \\mod {field}$')
            steps.add_step(f'slope = (3 * {self.x ** 2} + {self.curve.a})/({2 * self.y}) % {field}',
                           f'$\\lambda = \\frac{{3 * {self.x ** 2} + {self.curve.a}}}{{{2 * self.y}}} \\mod {field}$')
            steps.add_step(f'slope = ({3 * self.x ** 2} + {self.curve.a})/({2 * self.y}) % {field}',
                           f'$\\lambda = \\frac{{{3 * self.x ** 2} + {self.curve.a}}}{{{2 * self.y}}} \\mod {field}$')
            steps.add_step(f'slope = ({3 * self.x ** 2 + self.curve.a})/({2 * self.y}) % {field}',
                           f'$\\lambda = \\frac{{{3 * self.x ** 2 + self.curve.a}}}{{{2 * self.y}}} \\mod {field}$')
            steps.add_step(
                f'slope = ({3 * self.x ** 2 + self.curve.a} % {field})/({2 * self.y} % {field}) % {field}',
                f'$\\lambda = \\frac{{{3 * self.x ** 2 + self.curve.a} \\mod {field}}}{{{2 * self.y} \\mod {field}}} \\mod {field}$')
            steps.add_step(
                f'slope = {(3 * self.x ** 2 + self.curve.a) % field}/{2 * self.y % field} % {field}',
                f'$\\lambda = \\frac{{{(3 * self.x ** 2 + self.curve.a) % field}}}{{{2 * self.y % field}}} \\mod {field}$')
            dividend = (3 * x ** 2 + self.curve.a)
            divisor = (2 * y) % field
        steps.add_linebreak()
        steps.add_step(f'To continue we need the multiplicative inverse of {divisor} % {field}',
                       f'To continue we need the multiplicative inverse of ${divisor} \\mod {field}$')
        inverse = multiplicative_inverse(n=divisor, mod=field, steps=steps)
        steps.add_linebreak()
        steps.add_step(f'Now we can continue with {divisor}**-1 % {field} = {inverse}',
                       f'Now we can continue with ${divisor}^{{-1}} \\mod {field} = {inverse}$')
        steps.add_step(f'slope = {dividend}/{divisor} % {field}',
                       f'$\\lambda = \\frac{{{dividend}}}{{{divisor}}} \\mod {field}$')
        steps.add_step(f'slope = {dividend}*{inverse} % {field}', f'$\\lambda = {dividend}*{inverse} \\mod {field}$')
        steps.add_step(f'slope = {dividend * inverse} % {field}', f'$\\lambda = {dividend * inverse} \\mod {field}$')
        steps.add_step(f'slope = {dividend * inverse % field}', f'$\\lambda = {dividend * inverse % field}$')
        the_slope = dividend * inverse % field
        # We don't have to return steps because you still have it when you gave it too the function
        return the_slope

    def set_xy(self, x: int, y: int):
        """Set x and y attributes of the point in one go"""
        self.x = x
        self.y = y

    @property
    def inverse(self) -> 'EllipticPoint':
        """Return the inverse of the point. You can also use the ~ operator"""
        return self.__invert__()

    def __invert__(self) -> 'EllipticPoint':
        point = EllipticPoint(x=self.x, y=(-self.y) % self.curve.field, curve=self.curve, steps=self.steps)
        point.steps.add_step(f'{self}**-1 = ({self.x},{-self.y})', f'{self}**-1 = ({self.x},{-self.y})')
        point.steps.add_step(f'{self}**-1 = ({self.x},{-self.y} % {self.curve.field})',
                             f'{self}**-1 = ({self.x},{-self.y} \\mod {self.curve.field})')
        point.steps.add_step(f'{self}**-1 = ({self.x},{(-self.y) % self.curve.field})',
                             f'{self}**-1 = ({self.x},{(-self.y) % self.curve.field})')
        return point

    def __add__(self, point) -> 'EllipticPoint':
        """Add points together"""
        x, y = get_xy(point)

        steps = point.steps
        steps.add_step(f'We want to do P + Q = R where P = {self} and Q=({x}, {y}) on {self.curve}',
                       f'We want to do $P + Q = R$ where $P = {self}$ and $Q = ({x},{y})$ on ${self.curve.latex}$')
        if isinstance(point, InfinityPoint) or x == float('inf') and y == float('inf'):  # Adding infinity return copy
            steps.add_step(f'Adding the point at infinity to a point gives the identity so P + Q = {P}',
                           f'Adding the point at infinity to a point gives the identity so $P + Q = {P}$')
            return EllipticPoint(x=self.x, y=self.y, curve=self.curve, steps=steps)
        if self == ~EllipticPoint(x, y, self.curve):  # Return point at infinity if adding inverse
            steps.add_step(
                f'Adding a point and its inverse gives the point at infinity so P + Q = {InfinityPoint(curve=self.curve)}',
                f'Adding the point at infinity to a point gives the identity so $P + Q = \\mathcal{{O}}$')
            return InfinityPoint(curve=self.curve, steps=steps)

        steps.add_step(f'We need to find the slope (lambda)', f'We need to find the slope ($\\lambda$)')
        steps.add_linebreak()
        slope = self.slope(point, steps)

        x_result = (slope ** 2 - self.x - x) % self.curve.field

        steps.add_linebreak()
        steps.add_step(f'Now that we have slope = {slope} we can calculate R',
                       f'Now that we have $\\lambda = {slope}$ we can calculate $R$')
        steps.add_step('First we calcualte x_r', 'First we calculate $x_r$')
        steps.add_step('x_r = slope**2 - x_p - x_q', '$x_r = \\lambda^2 - x_p - x_q$')
        steps.add_step(f'x_r = {slope}**2 - {self.x} - {x}', f'$x_r = {slope}^2 - {self.x} - {x}$')
        steps.add_step(f'x_r = {slope ** 2} + {-self.x - x}', f'$x_r = {slope ** 2} + {-self.x - x}$')
        steps.add_step(f'x_r = {slope ** 2 - self.x - x}', f'$x_r = {slope ** 2 - self.x - x}$')
        steps.add_step(f'x_r = {slope ** 2 - self.x - x} % {self.curve.field}',
                       f'$x_r = {slope ** 2 - self.x - x} \\mod {self.curve.field}$')
        steps.add_step(f'x_r = {(slope ** 2 - self.x - x) % self.curve.field}',
                       f'$x_r = {(slope ** 2 - self.x - x) % self.curve.field}$')
        steps.add_linebreak()
        steps.add_step(f'Now with x_r = {x_result} we can calculate y_r',
                       f'Now with $x_r={x_result}$, $y_r$ can be calculated')

        y_result = (slope * (self.x - x_result) - self.y) % self.curve.field

        steps.add_step('y_r = slope * (x_p - x_r) - y_p', '$y_r = \\lambda * (x_p - x_r) - y_p$')
        steps.add_step(f'y_r = {slope} * ({self.x} - {x_result}) - {self.y}',
                       f'$y_r = {slope} * ({self.x} - {x_result}) - {self.y}$')
        steps.add_step(f'y_r = {slope} * ({self.x - x_result}) - {self.y}',
                       f'$y_r = {slope} * ({self.x - x_result}) - {self.y}$')
        steps.add_step(f'y_r = {slope * (self.x - x_result)} - {self.y}',
                       f'$y_r = {slope * (self.x - x_result)} - {self.y}$')
        steps.add_step(f'y_r = {slope * (self.x - x_result) - self.y}',
                       f'$y_r = {slope * (self.x - x_result) - self.y}$')
        steps.add_step(f'y_r = {slope * (self.x - x_result) - self.y} % {self.curve.field}',
                       f'$y_r = {slope * (self.x - x_result) - self.y} \\mod {self.curve.field}$')
        steps.add_step(f'y_r = {(slope * (self.x - x_result) - self.y) % self.curve.field}',
                       f'$y_r = {(slope * (self.x - x_result) - self.y) % self.curve.field}$')
        steps.add_linebreak()
        steps.add_step(f'{self} + ({x},{y}) = ({x_result}, {y_result})',
                       f'${self} + ({x},{y}) = ({x_result}, {y_result})$')
        steps.add_linebreak()
        return EllipticPoint(x=x_result, y=y_result, curve=self.curve, steps=steps)

    def __iadd__(self, point):
        """Add points together and also add steps"""
        result = self + point
        x, y = get_xy(result)
        self.set_xy(x, y)
        self.steps.last += result.steps
        return result

    def __mul__(self, n: int):
        """Calculate scalar multiplication"""
        the_n = n
        steps = Step()
        steps.add_step(f'We want to do [{n}]{self}', f'We want to do $[{n}]{self}$'),
        if n == 1:
            steps.add_step(f'[1]{self} is the identity the the result is {self}', f'$[1]{self} = {self}$'),
            return self
        steps.add_step(f'Start with {self} + {self}', f'Start with {self} + {self}')
        result = self + self
        steps += result.steps
        steps.add_step(f'{self} + {self} = {result}', f'${self} + {self} = {result}$')
        while n - 2 > 0:
            steps.add_step(f'[{n}]: {result} += {self}', f'[{n}] ${result} += {self}$')
            result += self
            steps += result.steps
            # Line break here...
            n -= 1
        if n == 1:
            steps.add_step(f'[{the_n}] One final {result} += {self}', f'$[{the_n}] final ${result} += {self}$')
            result += self
        steps.add_step(f'[{the_n}]{self} = {result}', f'$[{the_n}]{self} = {result} += {self}$')
        return EllipticPoint(x=result.x, y=result.y, curve=self.curve, steps=steps)

    def __imul__(self, n):
        result = self.__mul__(n)
        x, y = get_xy(result)
        self.set_xy(x, y)
        return result

    def __repr__(self):
        return f"EllipticPoint(x={self.x}, y={self.y}, curve={self.curve.__repr__()})"

    def __str__(self):
        return f"{self.x, self.y}"

    def latex_with_curve(self):
        return f'{self.x, self.y} \\in {self.curve}'

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
    def __init__(self, curve, steps=None):
        super(InfinityPoint, self).__init__(x=float('inf'), y=float('inf'), curve=curve, steps=steps)

    def __add__(self, point):
        """ The infinite point is the identity element """
        x, y = get_xy(point)
        print("No step here!!! infinity point __add__ maybe not even needed...")
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


def extended_gcd(a, b=None, steps=None):
    """Compute extended_gcd (I think) steps is not implemented for this, but you can pass it."""
    if a == 0:
        return b, 0, 1

    x, y, gcd = extended_gcd(b % a, a, steps)

    return y - (b // a) * x, x, gcd


def multiplicative_inverse(n, mod, steps: Step):
    """ Give the multiplicative inverse. Pass steps as an argument because this doesn't take a point. """
    if is_prime(mod):  # If n is prime we can use Fermat.
        steps.add_step(f'{mod} is prime so we can use Fermat.')
        steps.add_step(f'{n}**-1 % {mod} = {n}**({mod} - 2)', f'${n}^{{-1}} \\mod {mod} = {n}^{{{mod} - 2}}$')
        steps.add_step(f'{n}**-1 % {mod} = {n}**{mod - 2} % {mod}', f'${n}^{{-1}} \\mod {mod} = {n}^{{{mod - 2}}}$')
        steps.add_step(f'{n}**-1 % {mod} = {n ** (mod - 2)} % {mod}',
                       f'${n}^{{-1}} \\mod {mod} = {n ** (mod - 2)} \\mod {mod}$')
        steps.add_step(f'{n}**-1 % {mod} = {n ** (mod - 2) % mod}',
                   f'${n}^{{-1}} \\mod {mod} = {n ** (mod - 2) % mod}$')
        return n ** (mod - 2) % mod


    steps.add_step(f"{mod} is not prime so Fermat can't be used. Extended gcd has to be used.")
    x, y, gcd = extended_gcd(n, mod, steps)
    if gcd == 1:
        steps.add_step(f"With extended gcd I found {n}**-1 = {x}", f"With extended gcd I found ${n}^{{-1}} = {x}$")
        return x
    else:
        steps.add_step(f"With extended gcd I found that there is no inverse!")
    raise ValueError(f"No inverse of {n}  mod {mod}!\n\nSteps to get here:\n{steps}")


__author__ = "Quinten Cabo"


if __name__ == '__main__':
    # How to use:

    C = EllipticCurve(a=5, b=4, field=13)
    G = EllipticPoint(1, 7, curve=C)
    B = EllipticPoint(4, 7, curve=C)
    
    S = B + B + B
    
    print(S.steps)

