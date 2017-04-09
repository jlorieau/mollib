Partial Alignment
=================
The ``pa`` command is used to fit residual dipolar couplings (RDCs) and residual
anisotropic chemical shifts (RACSs, sometimes known as RCSAs) from partially
aligned samples using NMR.

Usage
-----

    .. include:: output/cli_pa_help.html

Arguments
---------

    ``-a`` / ``--alignment``
        The file(s) with the RDC and RACS alignment data. These can be in
        either of the following formats:

        - The pa format. See :ref:`pa_format`


        - NMRPipe's DC format.

    ``-o`` / ``--out``
        (Optional) The filename for the output report. The output report is
        rendered in Markdown.

    ``--fix-sign`` / ``--nofix-sign``
        (Optional) Check to see if the sign of RDCs or RACSs of the same type
        need to be inverted to get a better fit. This operation is useful for
        automatically fixing the sign of couplings when the absolute value of
        the \|J+D\|- and \|J\|-couplings are measured. By default, the sign of
        couplings are checked and fixed.

    ``--fix-outliers`` / ``--nofix-outliers``
        (Optional) Check to see if there are outliers for each type of
        interaction. A warning outlier and a bad outlier are defined by those
        that give an alpha-critical cutoff of 95% and 99%, respectively,
        using a Grubbs test. If outliers are found, these will be removed from
        the fit and the reported statistics. By default, outliers are not
        removed.

.. _pa_format:

Partial Alignment Data File Format
----------------------------------

The file format has the following features:

1. The interaction labels for dipolar interactions refer to two atoms (ex:
   14N-H) and the interaction label for CSA interactions refer to one atom.

2. For dipolar interactions, redundant residue numbers and chain identifiers
   are not needed. For example, '14N-H' and '14N-14H' refer to the same dipole.

3. If the chain identifier is not specified, then the subunit 'A' is assumed.

4. Relative residue numbers are allowed. For example, '14N-C-1' is the same as
   the '14N-13C' dipole.

5. Errors are optional. If the error is not specified, a default value from
   the settings is used.

The partial alignment RDC and RACS data file has the following format:

::

    # Interaction   Value (Hz)   Error (optional)
    14N-H           -14.5        0.1
    15N-H             3.5
    A.16N-H          -8.5        0.2  # larger error

    A.16H-A.15C       0.5        0.1
    B.16H-B.15C       0.5        0.1

    # Residual anisotropic chemical shift data
    # Interaction   Value (ppb)   Error (optional)
    5C                112         1
    6C               -250

Examples
--------

