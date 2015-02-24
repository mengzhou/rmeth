# mapping.py: a module for mapping command calling of rmeth
# 
# Copyright (C) 2014 University of Southern California and
#                          Meng Zhou
# 
# Authors: Meng Zhou
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import sys, logging, subprocess

class Mapping:
  """The main class for read mapping.
  """
  def __init__(self, opt, files):
    # things passed from options
    self.info = opt.info
    self.warn = opt.warn
    self.error = opt.error
    self.mapper = opt.mapper
    self.directional = opt.directional
    self.isPairEnd = opt.isPairEnd
    # get all file names
    self.files = files
    self.file_list = self.files.get_file_list()
    self.index = opt.index
    self._indexCT = self.index + "/CTgenome"
    self._indexGA = self.index + "/GAgenome"

  def map_se(self):
    if self.directional:
      total = 4
    else:
      total = 2
    self.info("(1/%d) Mapping C to T converted read " + \
      "to C to T converted genome..."%total)
    self._call_mapper(self.file_list["read1CT"], self._indexCT)
    self.info("(2/%d) Mapping C to T converted read " + \
      "to G to A converted genome...")
    self._call_mapper(self.file_list["read1CT"], self._indexGA)
    if self.directional:
      self.info("(3/%d) Mapping G to A converted read " + \
        "to C to T converted genome..."%total)
      self._call_mapper(self.file_list["read1CT"], self._indexCT)
      self.info("(4/%d) Mapping G to A converted read " + \
        "to G to A converted genome...")
      self._call_mapper(self.file_list["read1CT"], self._indexGA)

  def _call_mapper(self, read, index, output):
    cmd = self._set_mapping_args(read, index)
    try:
      subprocess.check_call(cmd)
    except subprocess.CalledProcessError:
      self.error("An error occured for command:\n%s"%cmd)
      sys.exit(1)

  def _set_mapping_args(self, read, index):
    return []

class MappingWithTophat2(Mapping):
  def _set_mapping_args(self, read, index, output):
    args = {\
    }

class MappingWithSTAR(Mapping):
  def _set_mapping_args(self, read, index, output):
    args = {\
      "--runThreadN":"8", \
      "--readFilesIn":reads, \
      "--genomeDir":index, \
      "--outFileNamePrefix":output, \
      "--outFilterIntronMotifs":"RemoveNoncanonical", \
      #"--clip3pAdapterSeq":"", \
      "--outSAMstrandField":"intronMotif", \
      "--outSAMtype":"BAM Unsorted" \
    }
