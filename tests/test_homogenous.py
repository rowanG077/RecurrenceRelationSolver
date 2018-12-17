# -*- coding: utf-8 -*-

from .context import RecurrenceRelation, RecurrenceRelationParser

import unittest


class HomogeneousTestSuite(unittest.TestCase):
    """Test cases for homogeneous recurrence relations"""

    def setUp(self):
        self.parser = RecurrenceRelationParser()
        self.verifyUpto = 50
        self.deviation = 1.0 / 100.0

    def test_comass03(self):
        recurrence = """
            eqs :=
            [
            s(n) = -4*s(n-2) + 4*s(n-1),
            s(0) = 6,
            s(1) = 8
            ];
        """
        relation = self.parser.parse_recurrence(recurrence)
        relation.verify_range(self.verifyUpto, self.deviation)

    def test_comass07(self):
        recurrence = """
            eqs :=
            [
            s(n) = s(n-1)+s(n-2),
            s(0) = 1,
            s(1) = 1
            ];
        """
        relation = self.parser.parse_recurrence(recurrence)
        relation.verify_range(self.verifyUpto, self.deviation)


if __name__ == '__main__':
    unittest.main()