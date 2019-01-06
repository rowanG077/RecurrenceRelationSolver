# -*- coding: utf-8 -*-

from .context import RecurrenceRelationParser

from .RecurrenceTestSuite import RecurrenceTestSuite

import unittest

class NonHomogeneousTestSuite(RecurrenceTestSuite):
    """Test cases for non homogeneous recurrence relations"""

    def test_simple_twice_theorem6(self):
        recurrence = """
            eqs :=
            [
            s(n) = s(n-1) + 2^n + 1,
            s(0) = 0
            ];
        """
        self.verify_range(self.parser.parse_recurrence(recurrence))

    def test_comass16(self):
        recurrence = """
            eqs :=
            [
            s(n) = 8*s(n-2)-16*s(n-4) +n^3,
            s(0) = 0,
            s(1) = 1,
            s(2) = 2,
            s(3) = 3
            ];
        """
        self.verify_range(self.parser.parse_recurrence(recurrence))

    def test_comass33(self):
        recurrence = """
            eqs :=
            [
            s(n) = (9/2)*s(n-2) +(3/2)*s(n-3)-5*s(n-4)-3*s(n-5) + (n-5)^2-3*(n-5)+7,
            s(0) = 2,
            s(1) = 4,
            s(2) = 8,
            s(3) = 1,
            s(4) = 3
            ];
        """
        self.verify_range(self.parser.parse_recurrence(recurrence))
    
    def test_comass36(self):
        recurrence = """
            eqs :=
            [
            s(n) = -2*s(n-1)+11*s(n-2)+12*s(n-3)-36*s(n-4) +41^(n-4)+3,
            s(0) = 1,
            s(1) = 1,
            s(2) = 1,
            s(3) = 1
            ];
        """
        self.verify_range(self.parser.parse_recurrence(recurrence))

    def test_week7_exercise4(self):
        recurrence = """
            eqs :=
            [
            s(n) = 2*s(n-1) + n + 5,
            s(0) = 4
            ];
        """
        self.verify_range(self.parser.parse_recurrence(recurrence))

    def test_week7_exercise5a(self):
        recurrence = """
            eqs :=
            [
            s(n) = 8*s(n-2) - 16*s(n-4) + (-2)^n,
            s(0) = 0,
            s(1) = 1,
            s(2) = 2,
            s(3) = 3,
            ];
        """
        self.verify_range(self.parser.parse_recurrence(recurrence))

    def test_week7_exercise5b(self):
        recurrence = """
            eqs :=
            [
            s(n) = 8*s(n-2) - 16*s(n-4) + n^2*4*n,
            s(0) = 0,
            s(1) = 1,
            s(2) = 2,
            s(3) = 3,
            ];
        """
        self.verify_range(self.parser.parse_recurrence(recurrence))

    def test_week7_exercise5c(self):
        recurrence = """
            eqs :=
            [
            s(n) = 8*s(n-2) - 16*s(n-4) + n^4 * 2^n,
            s(0) = 0,
            s(1) = 1,
            s(2) = 2,
            s(3) = 3,
            ];
        """
        self.verify_range(self.parser.parse_recurrence(recurrence))

    def test_week7_exercise6(self):
        recurrence = """
            eqs :=
            [
            s(n) = -5*s(n-1) - 6*s(n-2) + 42 * 4^n,
            s(0) = 56,
            s(1) = 278
            ];
        """
        self.verify_range(self.parser.parse_recurrence(recurrence))

    def test_week7_exercise7(self):
        recurrence = """
            eqs :=
            [
            s(n) = 4*s(n-1) - 3*s(n-2) + 2^n + n + 3,
            s(0) = 1,
            s(1) = 4
            ];
        """
        self.verify_range(self.parser.parse_recurrence(recurrence))

    def test_ex24_par8_2(self):
        recurrence = """
            eqs :=
            [
            s(n) = -5*s(n-1) - 6*s(n-2) + 42 * 4^n,
            s(1) = 56,
            s(2) = 278
            ];
        """
        self.verify_range(self.parser.parse_recurrence(recurrence))

    def test_example11(self):
        recurrence = """
            eqs :=
            [
            s(n) = 5*s(n-1) - 6*s(n-2) + 7^n,
            s(1) = 2,
            s(2) = 5
            ];
        """
        self.verify_range(self.parser.parse_recurrence(recurrence))

    def test_factorial(self):
        recurrence = """
            eqs :=
            [
            s(n) = s(n-1) + n!,
            s(1) = 1,
            s(2) = 2,
            s(3) = 2
            ];
        """
        self.verify_range(self.parser.parse_recurrence(recurrence))


if __name__ == '__main__':
    unittest.main()
