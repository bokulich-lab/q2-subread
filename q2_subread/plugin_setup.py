# ----------------------------------------------------------------------------
# Copyright (c) 2024, Michal Ziemski.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import Citations, Plugin
from q2_subread import __version__

citations = Citations.load("citations.bib", package="q2_subread")

plugin = Plugin(
    name="subread",
    version=__version__,
    website="https://github.com/mziemski/q2-subread",
    package="q2_subread",
    description="QIIME 2 plugin for wrapping the Subread software package.",
    short_description="Subread integration for QIIME 2.",
    # The plugin-level citation of 'Caporaso-Bolyen-2024' is provided as
    # an example. You can replace this with citations to other references
    # in citations.bib.
    citations=[citations['Caporaso-Bolyen-2024']]
)
