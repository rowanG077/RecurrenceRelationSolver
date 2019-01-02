# -*- coding: utf-8 -*-

from .context import RecurrenceRelation, RecurrenceRelationParser

from .RecurrenceTestSuite import RecurrenceTestSuite

import unittest


class HomogeneousTestSuite(RecurrenceTestSuite):
    """Test cases for homogeneous recurrence relations"""

    def test_comass03(self):
        recurrence = """
            eqs :=
            [
            s(n) = -4*s(n-2) + 4*s(n-1),
            s(0) = 6,
            s(1) = 8
            ];
        """

        self.verify_range(self.parser.parse_recurrence(recurrence))

    def test_comass07(self):
        recurrence = """
            eqs :=
            [
            s(n) = s(n-1)+s(n-2),
            s(0) = 1,
            s(1) = 1
            ];
        """
        self.verify_range(self.parser.parse_recurrence(recurrence))

    def test_week6_exercise4a(self):
        recurrence = """
            eqs :=
            [
            s(n) = s(n-1),
            s(0) = 2
            ];
        """
        self.verify_range(self.parser.parse_recurrence(recurrence))

    def test_week6_exercise4b(self):
        recurrence = """
            eqs :=
            [
            s(n) = -4*s(n-1) - 4*s(n-2),
            s(0) = 0,
            s(1) = 1
            ];
        """
        self.verify_range(self.parser.parse_recurrence(recurrence))

    def test_week6_exercise4c(self):
        recurrence = """
            eqs :=
            [
            s(n) = s(n-2)/4,
            s(0) = 1,
            s(1) = 0
            ];
        """
        self.verify_range(self.parser.parse_recurrence(recurrence))
    
    def test_week6_exercise6(self):
        recurrence = """
            eqs :=
            [
            s(n) = 10*s(n-1) - s(n-2),
            s(0) = 1,
            s(1) = 10
            ];
        """
        self.verify_range(self.parser.parse_recurrence(recurrence))


if __name__ == '__main__':
    unittest.main()