import logging

LOGGING_LEVEL = {
    "CRITICAL": 50, "ERROR": 40, "WARNING": 30,
    "INFO"    : 20, "DEBUG": 10, "NOTSET" : 0
}

def init_logger(level:str="INFO"):
    # Configurate logger
    logging_level = LOGGING_LEVEL[level]
    format = "[%(levelname)s] %(asctime)s, at %(filename)s:%(lineno)s - %(funcName)10s()\n%(message)s"
    logging.basicConfig(format=format, level=logging_level, datefmt='%Y/%m/%d %I:%M:%S %p')

    # Get logger
    logger = logging.getLogger(__name__)
    logger.info("Logger initialized.")


