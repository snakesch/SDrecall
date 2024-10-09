import uuid
from io import StringIO
import inspect
import logging
import time
import traceback

class ColoredFormatter(logging.Formatter):
    COLORS = {
        'DEBUG': '\033[96m',  # Cyan
        'INFO': '\033[92m',   # Green
        'WARNING': '\033[93m',  # Yellow
        'ERROR': '\033[91m',  # Red
        'CRITICAL': '\033[95m',  # Magenta
    }
    RESET = '\033[0m'  # Reset color

    def format(self, record):
        log_color = self.COLORS.get(record.levelname, self.RESET)
        # Get the current time as a formatted string
        time_str = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
        # Format the log message
        formatted_message = f"[{record.levelname}] - ({time_str}) - [{record.module} - {record.funcName}:{record.lineno}] -- {record.getMessage()}"
        return f"{log_color}{formatted_message}{self.RESET}"

def init_logger(handler = logging.StreamHandler(), tag = ""):
    logger = logging.getLogger("SDrecall")
    # logger = logging.getLogger(f"Process-{tag}")
    # handler.setFormatter(ColoredFormatter("%(levelname)s:%(asctime)s:%(module)s:%(funcName)s:%(lineno)s:%(message)s"))
    # logger.addHandler(handler)
    # logger.setLevel(logging.INFO)
    return logger

def log_command(func):
    def wrapper(*args, **kwargs):
        arg_names = list(inspect.signature(func).parameters.keys())
        arg_values = list(args)
        arg_dict = dict(zip(arg_names, arg_values))
        arg_dict.update(kwargs)

        tmp_tag = str(uuid.uuid4())
        log_stream = StringIO()
        ch = logging.StreamHandler(log_stream)
        logger = init_logger(handler = ch, tag = tmp_tag)

        # Get the default values of the function's keyword arguments
        sig = inspect.signature(func)
        defaults = {k: v.default for k, v in sig.parameters.items() if v.default is not inspect.Parameter.empty}

        # Fill in missing keyword arguments with their default values
        for arg_name, default_value in defaults.items():
            if arg_name not in arg_dict:
                arg_dict[arg_name] = default_value

        # Convert arguments and their values to a string
        args_str = ', '.join([f"{k}={arg_dict[k]}" for k in arg_names])
        defaults_str = ', '.join([f"{k}={arg_dict[k]}" for k in defaults])

        logger.debug(f"Executing: {func.__name__}({args_str}, {defaults_str})")
        start = time.time()
        try:
            result = func(*args, logger=logger, **kwargs)
            end = time.time()
            logger.debug(f"Finished: {func.__name__}({args_str}, {defaults_str}) in {end - start} seconds")
            log_contents = log_stream.getvalue()
            log_stream.close()
        except Exception as e:
            end = time.time()
            logger.debug(f"Failed: {func.__name__}({args_str}, {defaults_str}) in {end - start} seconds")
            tb_str = traceback.format_exc()
            log_contents = log_stream.getvalue()
            log_stream.close()
            return (False, (e, tb_str), log_contents)
        else:
            return (True, result, log_contents)

    return wrapper

def log_decorator(func):
    def wrapper(*args, **kwargs):
        tmp_tag = str(uuid.uuid4())
        log_stream = StringIO()
        ch = logging.StreamHandler(log_stream)
        logger = init_logger(handler = ch, tag = tmp_tag)
        result = func(*args, logger = logger, **kwargs)

        log_contents = log_stream.getvalue()
        log_stream.close()

        return result, log_contents
    return wrapper

def error_handling_decorator(func):
    def wrapper(*args, **kwargs):
        tmp_tag = str(uuid.uuid4())
        log_stream = StringIO()
        ch = logging.StreamHandler(log_stream)
        logger = init_logger(handler = ch, tag = tmp_tag)
        try:
            result = func(*args, logger=logger, **kwargs)
            log_contents = log_stream.getvalue()
            log_stream.close()
        except Exception as e:
            tb_str = traceback.format_exc()
            log_contents = log_stream.getvalue()
            log_stream.close()
            return (False, (e, tb_str), log_contents)
        else:
            return (True, result, log_contents)

    return wrapper