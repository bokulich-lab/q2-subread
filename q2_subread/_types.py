# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the GPL-3.0 license.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from rachis.core.type import SemanticType
from rachis.plugin import model


class SubreadIndexFileFmt(model.BinaryFileFormat):
    def validate(self, level="max"):
        pass


class SubreadIndexDirFmt(model.DirectoryFormat):
    array = model.File(r".+\.array$", format=SubreadIndexFileFmt)
    tab = model.File(r".+\.tab$", format=SubreadIndexFileFmt)
    lowinf = model.File(r".+\.lowinf$", format=SubreadIndexFileFmt)
    reads = model.File(r".+\.reads$", format=SubreadIndexFileFmt)

    # @array.set_path_maker
    # def array_path_maker(self):
    #     return "subread-index.00.b.array"
    #
    # @tab.set_path_maker
    # def tab_path_maker(self):
    #     return "subread-index.00.b.tab"
    #
    # @lowinf.set_path_maker
    # def lowinf_path_maker(self):
    #     return "subread-index.lowinf"
    #
    # @reads.set_path_maker
    # def reads_path_maker(self):
    #     return "subread-index.reads"


SubreadIndex = SemanticType("SubreadIndex")
