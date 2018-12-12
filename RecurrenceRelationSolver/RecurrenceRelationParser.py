import logging
import glob
import os.path
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
        self._pattern = re.compile("s\((\d+|n)\) = ")
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
            line = re.sub(' +', ' ', line).strip()

            m = re.search(self._pattern, line)
            if not m:
                continue

            n = m.group(1)
            if n in parsed:
                msg = 'Multiple equation found with condition %s. in data %s' % (
                    n, data)
                self._logger.warning(msg)

            if n.isdigit():
                n = int(n)

            # Remove the leading s defintion and the comma at the end
            parsed[n] = line[6:].strip().rstrip(",")
            # replace ^ with **
            parsed[n] = re.sub('^', '**', parsed[n])

        # remove the recurrence from the dictionary
        # and create a RecurrenceRelation object
        recurrence = parsed.pop("n")
        return RecurrenceRelation(recurrence, parsed)
