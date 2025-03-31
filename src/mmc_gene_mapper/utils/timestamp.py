import datetime


def get_timestamp():
    """
    Return a string with the current timestamp
    """
    now = datetime.datetime.now()
    result = f"{now.year:04d}-{now.month:02d}"
    result += f"-{now.day:02d}-{now.hour:02d}"
    result += f"-{now.minute:02d}-{now.second:02d}"
    return result
