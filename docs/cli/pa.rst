.. _pa-command:

``pa`` command
==============
The ``pa`` command is used to fit residual dipolar couplings (RDCs) and residual
anisotropic chemical shifts (RACSs, sometimes known as RCSAs) from *partially
aligned* samples using NMR. The output table entries are colored for warning
outliers (yellow) and bad outliers (red).

.. include:: cmds/ml_pa_help.rst

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

.. note:: The models option (``-m``/``--models``) will load the models as
          multiple molecules to be fit together in the SVD rather than conduct
          a separate SVD for each.

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

``--fix-nh-scale`` / ``--nofix-nh-scale``
    (Optional) Check to see if RDCs and RACSs have been scaled to match the
    magnitude of N-H RDCs. If they have, scale them back down to their
    original values. By default, this fixer is **off**.

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

.. include:: cmds/ml_pa_2kxa_1.rst

The following example fits the deposited RDCs for the first alignment
(``--set 0``) dataset of ubiquitin (``-a 2MJB``) to the deposited NMR structure
(``-i 2MJB``). The RDCs for methyl groups are projected onto the
corresponding C-C bonds (``--project-methyls``) and outliers are removed
from the fit (``--fix-outliers``).

.. include:: cmds/ml_pa_2mjb_1.rst

This example is the same as the last one, however 'CE-HE', 'CD-HD' and 'CE-SD'
RDCs are excluded (``--exclude``) from the fit.

.. include:: cmds/ml_pa_2mjb_2.rst

Likewise, the crystal structure of ubiquitin (``-i 1UBQ``) can be used in
the fit. In this case, the structure is missing hydrogen atoms, and these
must be added (``--hydrogenate``).

.. include:: cmds/ml_pa_1ubq_1.rst

Tensor Conventions
------------------

In the absence of motion, dipolar tensors are axially symmetric (i.e.
:math:`\delta_{xx} = \delta_{yy}`) and the principal component
(:math:`\delta_{zz}`) is colinear with the internuclear vector in the principal
axis system (PAS).

Chemical shift tensors (CSA) may be axially asymmetric (i.e.
:math:`\delta_{xx} \neq \delta_{yy}`), and their geometries must be specified
in relation to internal atomic coordinates. We use the convention from
Cornilescu *et al.* [Cornilescu2000]_.

.. [Cornilescu2000] Cornilescu, G. & Bax, A. Measurement of Proton, Nitrogen,
    and Carbonyl Chemical Shielding Anisotropies in a Protein Dissolved in a
    Dilute Liquid Crystalline Phase. J. Am. Chem. Soc. 122, 10143–10154 (2000).

The literature reports both the chemical *shielding* tensor (\ :math:`\sigma`\ )
and the chemical *shift* tensor (\ :math:`\delta`\ ). The difference between
the two is an inversion of sign (i.e. :math:`\sigma = - \delta`\ ). As a result,
the ordering of components between different conventions will change. In the
Haeberlen convention, the chemical shift components are order by their
magnitudes.

.. math::

    | \delta_{zz} | \geq | \delta_{xx} | \geq | \delta_{yy} |

The isotropic component (\ :math:`\delta_{iso}`\ ) has already been subtracted
from the three components.

.. math::

    \delta_{iso} = \frac{1}{3} \left( \delta_{zz} + \delta_{xx} +
    \delta_{yy} \right)

.. figure:: img/backbone_tensor.*
    :align: right
    :width: 350
    :alt: Backbone CSA tensor conventions

    Backbone CSA tensor conventions

In the IUPAC convention, the components are normally ordered starting from the
largest component (with sign) as the '11' component. However, for chemical
*shielding* tensors, the '33' component is largest.

.. math::

    \sigma_{33} \geq \sigma_{22} \geq \sigma_{11}

Mollib uses the Haeberlen convention and chemical *shift* tensors. The backbone
H, C' and N CSA tensors are defined as follows:

The :sup:`13`\ C' tensor (blue) has the largest component
(\ :math:`\delta_{zz}`\ ) oriented orthogonal to the O-C-N plane, and it is
rotated about this component by the :math:`\alpha_z` angle.

The :sup:`15`\ N tensor (red) has the largest component
(\ :math:`\delta_{zz}`\ ) nearly colinear with the H-N bond, and it is rotated
away from the bond about the yy-component (orthogonal to the H-N-C' plane) with
an angle :math:`\beta_y`\ .

The :sup:`1`\ H tensor (green) has the largest component
(\ :math:`\delta_{zz}`\ ) nearly colinear with the H-N bond, and it is rotated
about the xx-component (orthogonal to the H-N-C' plane) by an angle
:math:`\gamma_x`\ .
