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

        - The pa format. ex:

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

        - NMRPipe's DC format.

    ``-o`` / ``--out`` (optional)
        The filename for the output report. The output report is rendered in
        Markdown.

    ``--fix-sign`` / ``--nofix-sign`` (optional)
        Check to see if the sign of RDCs or RACSs of the same type need to be
        inverted to get a better fit. This operation is useful for automatically
        fixing the sign of couplings when the absolute value of the \|J+D\|-
        and \|J\|-couplings are measured. By default, the sign of couplings are
        checked and fixed.

    ``--fix-outliers`` / ``--nofix-outliers`` (optional)
        Check to see if there are outliers for each type of interaction. A
        warning outlier and a bad outlier are defined by those that give an
        $alpha$-critical cutoff of 95% and 99%, respectively, using a Grubbs
        test. If outliers are found, these will be removed from the fit and the
        reported statistics. By default, outliers are not removed.

Examples
--------

