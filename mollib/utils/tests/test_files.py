"""
Unittests for file utility functions
"""
import unittest
import tempfile
import shutil
import os

from mollib.utils.files import write_file


class TestFileUtil(unittest.TestCase):
    def setUp(self):
        # Create a temporary directory
        self.tempdir = tempfile.mkdtemp()

    def tearDown(self):
        # Remove the directory after the test
        shutil.rmtree(self.tempdir)

    def test_write_ascii(self):
        """Test the write_file function for ASCII files."""
        txt = "Test text in temp file."
        test_path = os.path.join(self.tempdir, 'tmp.txt')
        returned_path = write_file(text=txt, filepath=test_path)

        # If successful, the write_file returns the path written to.
        self.assertEqual(test_path, returned_path)

        # Read the file and make sure it has the contents we wrote to it.
        with open(test_path, 'r') as f:
            self.assertEqual(txt, f.read())

        # Now see if we get a new file if we do not overwrite it.
        returned_path = write_file(text=txt, filepath=test_path,
                                   overwrite=False)

        # The returned path and the test path should be different because the
        # returned path now has a new version of the file
        self.assertNotEqual(test_path, returned_path)

        # The new file should still be readable though
        with open(returned_path, 'r') as f:
            self.assertEqual(txt, f.read())

        # Now try overwriting with a new message
        txt = "A new message"
        returned_path = write_file(text=txt, filepath=test_path,
                                   overwrite=True)

        # The returned path should be the same as the test_path because the
        # file was overwritten
        self.assertEqual(test_path, returned_path)

        # Read the file and make sure it has the contents we wrote to it.
        with open(test_path, 'r') as f:
            self.assertEqual(txt, f.read())

    def test_write_temp_ascii(self):
        """Test the writing of a temporary ASCII file."""
        txt = "Test text in another temp file."
        returned_path = write_file(text=txt, temporary=True)

        # The returned file should have the temporary path in it.
        temp_path = tempfile.gettempdir()
        self.assertIn(temp_path, returned_path)

        # Read the file and make sure it has the contents we wrote to it.
        with open(returned_path, 'r') as f:
            self.assertEqual(txt, f.read())
