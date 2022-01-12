#!/usr/bin/env python3

import sys

MIN_PYTHON = (3, 9)
if sys.version_info < MIN_PYTHON:
    sys.exit("Python %s.%s or later is required." % MIN_PYTHON)

from typing import Callable, Generator


def additive(value: int, adder: int, modulo: int) -> int:
    """Additive operation for group elements
        These operation functions have 3 inputs.
        :param value: The current value in the cycle
        :param adder: The value we always add
        :param modulo: The value to modulo with
    """
    return (value + adder) % modulo


additive.__name__ = '+'


__author__ = "Quinten Cabo"


def multiplicative(value, multiplier, modulo):
    """Additive operation for group elements
           These operation functions have 3 inputs.
           :param value: The current value in the cycle
           :param multiplier: The value we always add
           :param modulo: The value to modulo with
       """
    return (value * multiplier) % modulo


multiplicative.__name__ = 'x'


class Group:

    def __init__(self, n: int, operation: Callable, exclude=None):
        if exclude is None:
            exclude = set()
        self.n = n
        self.operation = operation
        # This is here because len needs exclude to only include real excludes it is better to do this only once here.
        exclude = set(filter(lambda x: x < n, exclude))  # don't exclude elements not in the set.
        self.exclude = exclude

    def __iter__(self) -> Generator:
        """ Returns group as python generator """
        for i in range(self.n):
            if i in self.exclude:
                continue
            yield i

    def __contains__(self, item: int) -> bool:
        if item in self.exclude:
            return False
        return self.n > item >= 0

    def __len__(self) -> int:
        return self.n - len(self.exclude)  # -1 because n is not in group?

    @property
    def order(self):
        """ Order of the group """
        return len(self)

    @property
    def neutral_element(self):
        """
        The neutral element of a group is an element that when the operation is applied to
        Another element of the group
        """
        for e in self:
            for item in self:
                if self.operation(item, e, self.n) != item:
                    break
            else:
                return e
        return None

    @property
    def generators(self) -> set:
        """Returns a set of all elements of this group that can generate the elements in the group"""
        return {i for i in self if self.is_generator(i)}

    @property
    def orders(self) -> set:
        """Returns a set of all the orders that an element in this group can have"""
        return {self.order_of(i) for i in self}

    def from_set(the_set: set, operation: Callable):
        """Returns group based on set"""
        n = max(the_set) + 1
        exclude = {i for i in range(n) if i not in the_set}
        return Group(n, operation, exclude)

    def order_of(self, item: int):
        """ Returns the order of an item in the group """
        if item not in self:
            raise ValueError(f"{item} is not in group!")
        past_values = set()
        value = 1
        order = 0
        while value not in past_values:
            past_values.add(value)
            value = self.operation(value, item, self.n)
            order += 1
        return order

    def cycle_of(self, g: int) -> list:
        """ Returns the cycle of an item in the group """
        if g not in self:
            raise ValueError(f"{g} is not in {self}")
        now = g
        past_values = []
        while now not in past_values:
            past_values.append(now)
            now = self.operation(past_values[-1], g, self.n)
        return past_values

    def inverse_of(self, item: int) -> int:
        """ Returns the inverse of an item in the group """
        if item not in self:
            raise ValueError(f"{item} not in {self}")
        identity = self.neutral_element
        for a in self:
            if self.operation(item, a, self.n) == identity:
                return a
        else:
            return None

    def is_generator(self, item: int) -> bool:
        """Returns true if item is a generator of every item in the group"""
        if self.order_of(item) == self.order:
            return True
        return False

    def is_subgroup(self, potential_subgroup: set) -> bool:
        """
        (B, ★) is a subgroup of (A, ★) if
        - e ∈ B neutral element is in B
        - B is a subset A
        - ∀a, b ∈ B: a ★ b ∈ B, is result of operation in with any value of in subgroup in group
        - ∀a ∈ B: the inverse of a is in B
        """
        # Order of B must be divisor of order of A (Lagrange)
        if len(self) % len(potential_subgroup) != 0:
            return False

        # is neutral element in subgroup
        if self.neutral_element not in potential_subgroup:
            return False

        # is potential_subgroup a subset A
        for b in potential_subgroup:
            if b not in self:
                return False

        # ∀a, b ∈ B: a ★ b ∈ B
        for b1 in potential_subgroup:
            for b2 in potential_subgroup:
                if self.operation(b1, b2, self.n) not in potential_subgroup:
                    return False

        # ∀a ∈ B, inverse of a is in B
        for b in potential_subgroup:
            if self.inverse_of(b) not in potential_subgroup:
                return False

        return True

    def __repr__(self):
        op = self.operation.__name__
        string = f"(Z/{self.n}Z, {op})"
        if self.exclude:
            string = "(" + string + f"\\{self.exclude})"
        return string

    @property
    def latex_repr(self):
        """ Latex representation of object. You have to print this. """
        op = r"\text{"f"{self.operation.__name__}""}"
        string = f"(\\mathbb{{Z}}/{self.n}\mathbb{{Z}},{op})"
        if self.exclude:
            exclude_str = f"\\{str(self.exclude)[:-1]}\\""}"
            string = "(" + string + f"\\backslash{exclude_str})"
        return string


class MultiplicativeGroup(Group):
    def __init__(self, n: int, excludes: set = None):
        super().__init__(n, multiplicative, excludes)


class AdditiveGroup(Group):
    def __init__(self, n: int, excludes: set = None):
        super().__init__(n, additive, excludes)


class MultiplicativeGroup_(MultiplicativeGroup):
    """ MultiplicativeGroup with {0} excluded """

    def __init__(self, n: int):
        super().__init__(n, excludes={0})
