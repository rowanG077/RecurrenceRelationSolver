# -*- coding: utf-8 -*-

from .context import RecurrenceRelation, RecurrenceRelationParser

import unittest
from nose.tools import assert_almost_equal

class RecurrenceTestSuite(unittest.TestCase):
    """Test cases for homogeneous recurrence relations"""

    def setUp(self):
        self.parser = RecurrenceRelationParser()

    def verify_range(self, relation):
        for i in range(0, 51):
            iterative_result = relation.calculateValueFromRecurrence(i)
            solved_result = relation.calculateValueFromSolved(i)
            msg = "Verification of solved recurrence failed at n = %d for relation: %s" % (i, str(relation._recurrence)) 
            assert_almost_equal(iterative_result, solved_result, 4, msg)
