import unittest

from mollib.utils import MDTable


class TestMarkdown(unittest.TestCase):

    maxDiff = None

    def test_MDTable(self):
        "Tests the Markdown table rendering."

        table = MDTable('one', 'two', 'three', 'four')

        # Check the header rendering
        self.assertEqual(table.content(),
                                  ("\n"
                                   "-----------------------\n"
                                   "one  two  three  four  \n"
                                   "---- ---- ------ ------\n"
                                   "-----------------------\n"
                                   "\n"))

        # Add rows, and check the rendering
        table.add_row('First', 'Second item', 'Third item', 'Fourth')
        table.add_row('Item 1', 'Item 2 and change', 'Item 3', 'Item Four')

        # Adding more or fewer columsn to a row produces an AssertionError
        with self.assertRaises(AssertionError):
            table.add_row('1', '2', '3')

        # This table does not render as a multiline
        self.assertEqual(table.content(),
                         ("\n"
                          "--------------------------------------------------\n"
                          "one     two                three       four       \n"
                          "------- ------------------ ----------- -----------\n"
                          "First   Second item        Third item  Fourth     \n"
                          "Item 1  Item 2 and change  Item 3      Item Four  \n"
                          "--------------------------------------------------\n"
                          "\n"))

        # This table will render as a multiline
        table.max_width = 40
        self.assertEqual(table.content(),
                         ("\n"
                          "----------------------------------------\n"
                          "one     two      three       four       \n"
                          "------- -------- ----------- -----------\n"
                          "First   Second   Third item  Fourth     \n"
                          "        item                            \n"
                          "\n"
                          "Item 1  Item 2   Item 3      Item Four  \n"
                          "        and                             \n"
                          "        change                          \n"
                          "\n"
                          "----------------------------------------\n"
                          "\n"))
