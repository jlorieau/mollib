"""
Utilities to format text for output.
"""

def word_list(words):
    """Group a list of words into a string with a list separated by commas
    and the word 'and'.
    
    Parameters
    ----------
    words: list of str
        A list of strings (words) to group together into a word list.
    
    Returns
    -------
    str
    
    Examples
    --------
    >>> word_list('February')
    'February'
    >>> word_list(['February', 'March'])
    'February and March'
    >>> word_list(['February', 'March', 'April', 'May'])
    'February, March, April and May'
    """
    # Convert words into a list, if it's a string
    if isinstance(words, str):
        words = [words,]

    # Generate the word list, depending on the number of items in the words
    # list
    num_words = len(words)

    if num_words == 0:
        return ''
    elif num_words == 1:
        return words[0]
    elif num_words == 2:
        return ' and '.join(words)
    else:
        return ', '.join(words[0:-1]) + ' and ' + words[-1]