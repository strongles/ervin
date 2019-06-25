from unittest import TestCase
from ervin.ervin_utils import read_and_sanitise_raw_data, read_from_fasta_file
from mock import patch, mock_open
import os


def dummy_raw_data():
    with open("tests/test_data/dummy_probe_blaster_source_data.txt") as file_in:
        print(os.getcwd())
        return file_in.read()


def dummy_expected_sanitised_data():
    return [
        ">something",
        "a",
        "multi",
        "line",
        "sequence",
        ">something else",
        "another",
        "shorter",
        "sequence"
    ]


class TestProbeBlaster(TestCase):

    @patch("builtins.open", new_callable=mock_open, read_data=dummy_raw_data())
    def test_raw_data_correctly_sanitised(self, *_):
        expected_data = dummy_expected_sanitised_data()
        actual_data = read_and_sanitise_raw_data("some_filename")
        self.assertListEqual(expected_data, actual_data)

    @patch("ervin.ervin_utils.read_and_sanitise_raw_data", return_value=dummy_expected_sanitised_data())
    def test_read_probe_records(self, mock_read_data):
        expected_records = [
            {
                "title": ">something",
                "seq": "amultilinesequence"
            },
            {
                "title": ">something else",
                "seq": "anothershortersequence"
            }
        ]
        actual_records = read_from_fasta_file("filename")
        mock_read_data.assert_called_once_with("filename")
        self.assertListEqual(expected_records, actual_records)
