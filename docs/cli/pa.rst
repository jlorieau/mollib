.. _pa-command:

``pa`` command
==============
The ``pa`` command is used to fit residual dipolar couplings (RDCs) and residual
anisotropic chemical shifts (RACSs, sometimes known as RCSAs) from *partially
aligned* samples using NMR. The output table entries are colored for warning
outliers (yellow) and bad outliers (red).

Usage
-----

    .. include:: output/cli_pa_help.html

Arguments
---------

    ``-a`` / ``--alignment`` ``filename``
        The file(s) with the RDC and RACS alignment data. These can be in
        either of the following formats:

        - The pa format. See :ref:`pa_format`.

        - NMRPipe's DC format.

        - Magnetic resonance data files (``.mr``) submitted to the PDB. This
          function supports the automatic fetching and caching of magnetic
          resonance data files.

    ``-o`` / ``--out`` ``filename``
        (Optional) The filename for the output report. The output report is
        rendered in Markdown.

    ``-p`` / ``--pred`` ``filename``
        (Optional) The filename for the back-calculated RDCs and RACS from the
        SVD fit. The output report is rendered in Markdown.

    ``--summary``
        (Optional) Only display the fit summary.

    ``--exclude`` ``interactions-types``
        (Optional) Exclude one or more interactions types. ex: ``--exclude
        N-H CA-HA`` will exclude all N-H and CA-HA RDCs.

    ``--set`` ``id``
        (Optional) Use the given data set, if multiple data sets are available.
        This option is useful with ``.mr`` data from the PDB, which may contain
        mulitple alignment data sets from multiple alignment media. Sets can
        be selected from their alignment tensor value (ex: 500, 501, etc) or
        from their position within the data file, starting with 0. (ex: 0 for
        the first dataset, 1, for the second dataset and so on.)

    ``--project-methyls``
        (Optional) Use the C-C bond RDC values for the methyl ¹H-¹³C RDCs. This
        is the convention followed by X-plor NIH. By default, this is disabled.

    ``--methyl-scale`` ``number``
        (Optional) The scaling constant to use in fitting the methyl RDCs. This
        scaling may be needed if the contribution of the C3-rotational motion
        was not accounted for in the reported RDCs. By default, this value is
        1.0.

Fixer Arguments
^^^^^^^^^^^^^^^

    ``--fix-sign`` / ``--nofix-sign``
        (Optional) Check to see if the sign of RDCs or RACSs of the same type
        need to be inverted to get a better fit. This operation is useful for
        automatically fixing the sign of couplings when the absolute value of
        the \|J+D\|- and \|J\|-couplings are measured. By default, this fixer
        is **on**.

    ``--fix-outliers`` / ``--nofix-outliers``
        (Optional) Check to see if there are outliers for each type of
        interaction. A warning outlier and a bad outlier are defined by those
        that give an alpha-critical cutoff of 95% and 99%, respectively,
        using a Grubbs test. If outliers are found, these will be removed from
        the fit and the reported statistics. By default, this fixer is **off**.

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

.. code-block:: none

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

The following example fits the deposited RDCs for the hemagglutin fusion
peptide structure (``-a 2KXA``) to the deposited NMR structure
(``-i 2KXA``). The output table entries are colored for warning outliers
(yellow) and bad outliers (red).

.. include:: output/cli_pa_i_2KXA_a_2KXA.html

The following example fits the deposited RDCs for the first alignment
(``--set 0``) dataset of ubiquitin (``-a 2MJB``) to the deposited NMR structure
(``-i 2MJB``). The RDCs for methyl groups are projected onto the
corresponding C-C bonds (``--project-methyls``) and outliers are removed
from the fit (``--fix-outliers``).

.. include:: output/cli_pa_i_2MJB_a_2MJB_set_0_fix-outliers_project-methyls_summary.html

This example is the same as the last one, however 'CE-HE', 'CD-HD' and 'CE-SD'
RDCs are excluded (``--exclude``) from the fit.

.. include:: output/cli_pa_i_2MJB_a_2MJB_set_0_exclude_CE-HE_CD-HD_CE-SD_fix-outliers_project-methyls_summary.html

Likewise, the crystal structure of ubiquitin (``-i 1UBQ``) can be used in
the fit. In this case, the structure is missing hydrogen atoms, and these
must be added (``--hydrogenate``).

.. include:: output/cli_pa_i_1UBQ_a_2MJB_set_0_fix-outliers_project-methyls_hydrogenate_summary.html
