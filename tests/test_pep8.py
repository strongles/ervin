import logging
import sys
from flake8.api import legacy as flake8
from io import StringIO
from unittest import TestCase


class TestPep8(TestCase):

    def setUp(self) -> None:
        self.old_stdout = sys.stdout
        sys.stdout = self.new_stdout = StringIO()

        logging.getLogger("flake8").setLevel(logging.CRITICAL)

    def tearDown(self) -> None:
        sys.stdout = self.old_stdout
        logging.getLogger("flake8").setLevel(logging.INFO)

    def testPep8(self) -> None:

        style = flake8.get_style_guide(append_config="setup.cfg")
        report = style.check_files(["tests", "src"])
        errors = report.get_statistics("")
        error_string = f"\r\n{self.new_stdout.getvalue()}"

        self.assertEqual(len(errors), 0, error_string)
