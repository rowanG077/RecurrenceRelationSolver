#!/usr/bin/env python3
# coding=utf-8
import logging
import re

from . import RecurrenceRelation


class RecurrenceRelationParser(object):
    """
    RecurrenceParser object that handles parsing of recurrence relations into RecurrenceRelationObject
    """

    def __init__(self):
        """
        create _RecurrenceParser object
        """
        self._pattern = re.compile("s\((\d+|n)\)\s*=\s*?(.+)")
        self._logger = logging.getLogger(__name__)

    def parse_recurrence(self, data):
        """
        parse and create a RecurrenceRelation object for it

        Args:
            filename (string): The file to read

        Returns:
            RecurrenceRelation: The parsed recurrence relation
        """
        parsed = {}
        for line in data.splitlines():
            line = line.strip()

            m = re.search(self._pattern, line)
            if not m:
                continue

            n = m.group(1)
            eq = m.group(2).rstrip(",").strip()
            if n in parsed:
                msg = 'Multiple equation found with condition %s. in data %s' % (
                    n, data)
                self._logger.warning(msg)

            # we use integers as keys in a dict for the initial conditions
            if n.isdigit():
                n = int(n)
            # if it is not an initial condition we append -s(n) so we can an
            # equal to 0
            else:
                eq += "-s(n)"

            # Remove the leading s defintion and the comma at the end
            parsed[n] = eq

        # remove the recurrence from the dictionary
        # and create a RecurrenceRelation object
        recurrence = parsed.pop("n")
        return RecurrenceRelation(recurrence, parsed)
