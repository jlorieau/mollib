# -*- coding: utf-8 -*-

import unittest
try:
    import configparser
except ImportError:
    import ConfigParser as configparser

from mollib.core.settings import import_settings

class TestSettings(unittest.TestCase):
    "Tests the settings functions."

    def test_import_settings(self):
        "Tests the basic settings import functionality of import_settings."

        # Make a fake settings module using a class
        class Settings(object):
            def test_method(self): pass

        settings = Settings()
        settings.value1 = 'original'
        settings.value2 = 1.0
        settings.value3 = 1
        settings.value4 = False
        settings.value5 = 1

        # Setup the config
        config = configparser.ConfigParser()
        sample_config = """
        [settings.nonmatching]
        value1 = 'new'
        value2 = 2.0
        value3 = 1
        value4 = True
        value5 = 'different type'
        test_method = False
        """
        config.read_string(sample_config)

        # There are no matching sections, so this shouldn't change the values
        import_settings(config, 'settings', settings.__dict__)
        self.assertEqual(settings.value1, 'original')
        self.assertTrue(isinstance(settings.value1, str))

        self.assertEqual(settings.value2, 1.0)
        self.assertTrue(isinstance(settings.value2, float))

        self.assertEqual(settings.value3, 1)
        self.assertTrue(isinstance(settings.value3, int))

        self.assertEqual(settings.value4, False)
        self.assertTrue(isinstance(settings.value4, bool))

        self.assertTrue(callable(settings.test_method))

        # Update the config with matching values, which should replace the
        # values
        sample_config += """
        [settings]
        value1 = 'new'
        value2 = 2.0
        value3 = 2
        value4 = True
        value5 = 'different type'
        test_method = False
        """
        config.read_string(sample_config)

        # There is now a matching section and all specified values with a
        #  should be overwritten.
        import_settings(config, 'settings', settings.__dict__)
        self.assertEqual(settings.value1, 'new')
        self.assertTrue(isinstance(settings.value1, str))

        self.assertEqual(settings.value2, 2.0)
        self.assertTrue(isinstance(settings.value2, float))

        self.assertEqual(settings.value3, 2)
        self.assertTrue(isinstance(settings.value3, int))

        self.assertEqual(settings.value4, True)
        self.assertTrue(isinstance(settings.value4, bool))

        self.assertTrue(callable(settings.test_method))

        # Values without a matching type in the config should not be
        # overwritten
        self.assertEqual(settings.value5, 1)

        # Callables shouldn't be overwritten
        self.assertTrue(callable(settings.test_method))



