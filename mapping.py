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

import sys, os, re
import logging, subprocess
from optparse import OptionParser
from random import choice

class Mapping:
  """The main class for read mapping.
  """
  def __init__(self, opt, files):
    self.files = files
    self.file_list = self.files.get_file_list()
    self.info = opt.info
    self.warn = opt.warn
    self.error = opt.error
    self.mapper = opt.mapper
    self.isPairEnd = opt.isPairEnd
    self.__get_file_paths(opt)
    self.__set_default_par()

  def map(self):
    """Start mapping.
    """
    #self.info("Mapping reads to C to T converted genome...")
    #self.info("Command: %s"%" ".join(self.args_map_to_CT))
    self.__map_to_CT()

    #self.info("Mapping reads to G to A converted genome...")
    #self.info("Command: %s"%" ".join(self.args_map_to_GA))
    self.__map_to_GA()

  def __map_to_CT(self):
    print "mapping to CT"

  def __map_to_GA(self):
    print "mapping to GA"

  def __get_file_paths(self, opt):
    self.mate1 = opt.mate1
    self.mate1_CT = opt.mate1_CT
    self.mate1_GA = opt.mate1_GA
    if opt.isPairEnd:
      self.mate2 = opt.mate2
      self.mate2_CT = opt.mate2_CT
      self.mate2_GA = opt.mate2_GA

  def __set_default_par(self):
    self.par = {}

class MappingWithTophat2(Mapping):
  def __set_default_par(self):
    self.par = {\
      "--readFilesIn":$READ, \
      "--outFileNamePrefix":$OUTPUT, \
      "--outFilterIntronMotifs":"RemoveNoncanonical", \
      "--outSAMstrandField":"intronMotif", \
      "--clip3pAdapterSeq":"", \
      "--outSAMtype":"BAM Unsorted" \
    }

class MappingWithSTAR(Mapping):
  def __set_default_par(self):
    if self.isPairEnd:
      reads = " ".join(self.files.get_mate1(), self.files.get_mate2())
    else:
      reads = self.files.get_mate1()
    self.par = {\
      "--readFilesIn":reads, \
      "--outFileNamePrefix":$OUTPUT, \
      "--outFilterIntronMotifs":"RemoveNoncanonical", \
      "--outSAMstrandField":"intronMotif", \
      "--clip3pAdapterSeq":"", \
      "--outSAMtype":"BAM Unsorted" \
    }
