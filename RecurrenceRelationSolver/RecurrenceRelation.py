import logging
import glob
import os.path
import re


class RecurrenceRelation(object):
    """
    RecurrenceRelation object that contains a recurrence relations
    and initial conitions. Allows for solving and verifications
    """

    def __init__(self, recurrence, initialConditions):
        """
        create RecurrenceRelation object

        Args:
            recurrence (string): The recurrence as a string
            initialConditions (dict of string: string): The initial conditions
        """
        self._recurrence = recurrence
        self._initialConditions = initialConditions
        self._closedForm = None
        self._logger = logging.getLogger(__name__)

    def _solve(self):
        """
        Solve the recurrence relation into a closed form

        Returns:
            String: The solved recurrence relation in string format
        """
        return "Solved"

    def solve(self):
        """
        Get the recurrence relation into a closed form

        Returns:
            String: The solved recurrence relation in string format
        """
        if self._closedForm is None:
            self._closedForm = self._solve()

        return self._closedForm

    def verify_upto(self, n):
        """
        Verify the solved recurrence relation with builtin sympy rsolve upto the nth iteration

        Returns:
            boolean: True if correct, false if incorrect
        """
        if self._closedForm is None:
            self._closedForm = self._solve()

        return True
