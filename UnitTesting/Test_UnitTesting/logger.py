import unittest
import logging

class TestLogger(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        logging.basicConfig(level='INFO')
        logger = logging.getLogger(__name__)
        fh = logging.FileHandler('logging_output.log')
        fh.setLevel('DEBUG')
        logger.addHandler(fh)
        cls.logger = logger

    def test_log(self):

        self.logger.info(' Info message')
        self.logger.debug(' Debug message')

    def test_zchange_level(self):

        log_to_console_and_file(self, 'INFO', ' Hello')


def log_to_console_and_file(self, level, message):

    if level == 'INFO':
        level = 20

    self.logger.log(level, message)


if __name__ == '__main__':
    unittest.main()