import functools
import json
import os
import pathlib

from templates.config_template import CONFIG_TEMPLATE


DEFAULT_CONFIG_PATH = pathlib.Path.home() / ".config/ervin/config.json"


class InvalidConfigFile(Exception):
    def __init__(self, filename):
        message = f"Invalid config file provided: {filename}"
        super().__init__(self, message)
    pass


class Config:
    def __init__(self, config_dict):
        for attribute, value in config_dict.items():
            setattr(self, attribute, value)


def write_default_config():
    config_directory = DEFAULT_CONFIG_PATH.parent
    if not config_directory.exists():
        config_directory.mkdir(parents=True)
    with open(DEFAULT_CONFIG_PATH, "w") as config_file_out:
        json.dump(CONFIG_TEMPLATE, config_file_out)


@functools.lru_cache()
def get_config():
    config_filepath = os.environ.get("ERVIN_CONFIG_PATH", None)
    if not config_filepath:
        if not DEFAULT_CONFIG_PATH.exists():
            write_default_config()
        config_filepath = DEFAULT_CONFIG_PATH
    with open(config_filepath) as config_in:
        try:
            return Config(json.load(config_in))
        except json.JSONDecodeError:
            raise InvalidConfigFile(config_filepath)
