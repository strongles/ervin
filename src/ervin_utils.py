import datetime


def format_timestamp_for_filename():
    current_timestamp = datetime.datetime.now()
    return current_timestamp.strftime("%Y-%m-%d_%H-%M-%S")
