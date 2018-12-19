# -*- coding: utf-8 -*-

from .context import RecurrenceRelationParser

import unittest

class NonHomogeneousTestSuite(unittest.TestCase):
    """Test cases for non homogeneous recurrence relations"""

    def setUp(self):
        self.parser = RecurrenceRelationParser()
        self.verifyUpto = 50
        self.deviation = 1.0 / 100.0

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
        relation = self.parser.parse_recurrence(recurrence)
        relation.verify_range(self.verifyUpto, self.deviation)

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
        relation = self.parser.parse_recurrence(recurrence)
        relation.verify_range(self.verifyUpto, self.deviation)
    
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
        relation = self.parser.parse_recurrence(recurrence)
        relation.verify_range(self.verifyUpto, self.deviation)

    def test_Wim(self):
            recurrence = """
                eqs :=
                [
                s(n) = s(n-2) + 0.5*n^2+0.5*n,
                s(0) = 0,
                s(1) = 1
                ];
            """
            relation = self.parser.parse_recurrence(recurrence)
            relation.verify_range(self.verifyUpto, self.deviation)




if __name__ == '__main__':
    unittest.main()
