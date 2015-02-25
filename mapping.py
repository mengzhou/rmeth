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

import sys, logging, subprocess, re

class Mapping:
  """The main class for read mapping.
  """
  def __init__(self, files, mapper, isDirectional, isPairEnd, \
    index):
    logging.basicConfig(level=20,
      format='[%(levelname)-5s][%(asctime)s] %(message)s ',
      datefmt='%H:%M:%S', stream=sys.stderr, filemode="w")
    self.info = logging.info
    self.warn = logging.warn
    self.error = logging.error
    self.mapper = mapper
    self.isDirectional = isDirectional
    self.isPairEnd = isPairEnd
    self.files = files
    self.file_list = self.files.get_file_list()
    self.index = index
    self._indexCT = self.index + "/CTgenome"
    self._indexGA = self.index + "/GAgenome"

  def map(self):
    total = 2
    if not self.isDirectional:
      total = 4

    if not self.isPairEnd:
      self.info("(1/%d) Mapping C to T converted read1 "%total + \
        "to C to T converted genome...")
      self._call_mapper_se(self.file_list["read1CT"], self._indexCT, \
        self.file_list["read1CT_to_CT"])
      self.info("(2/%d) Mapping C to T converted read1 "%total + \
        "to G to A converted genome...")
      self._call_mapper_se(self.file_list["read1CT"], self._indexGA, \
        self.file_list["read1CT_to_GA"])
      if not self.isDirectional:
        self.info("(3/%d) Mapping G to A converted read1 "%total + \
          "to C to T converted genome...")
        self._call_mapper_se(self.file_list["read1GA"], self._indexCT, \
          self.file_list["read1GA_to_CT"])
        self.info("(4/%d) Mapping G to A converted read1 "%total + \
          "to G to A converted genome...")
        self._call_mapper_se(self.file_list["read1GA"], self._indexGA, \
          self.file_list["read1GA_to_GA"])
    else:
      self.info("(1/%d) Mapping C to T converted read1 "%total + \
        "and G to A converted read 2 to C to T converted genome...")
      self._call_mapper_pe(self.file_list["read1CT"], \
        self.file_list["read2GA"], self._indexCT, \
        self.file_list["read1CTread2GA_to_CT"])
      self.info("(2/%d) Mapping C to T converted read1 "%total + \
        "and G to A converted read 2 to G to A converted genome...")
      self._call_mapper_pe(self.file_list["read1CT"], \
        self.file_list["read2GA"], self._indexGA, \
        self.file_list["read1CTread2GA_to_GA"])
      if not self.isDirectional:
        self.info("(3/%d) Mapping G to A converted read1 "%total + \
          "and C to T converted read 2 to C to T converted genome...")
        self._call_mapper_pe(self.file_list["read1GA"], \
          self.file_list["read2CT"], self._indexCT, \
          self.file_list["read1GAread2CT_to_CT"])
        self.info("(4/%d) Mapping G to A converted read1 "%total + \
          "and C to T converted read 2 to G to A converted genome...")
        self._call_mapper_pe(self.file_list["read1GA"], \
          self.file_list["read2CT"], self._indexGA, \
          self.file_list["read1GAread2CT_to_GA"])

  def _call_mapper_se(self, read, index, output):
    cmd = self._set_mapping_args_se(read, index, output)
    self.info(" ".join(cmd))
    #try:
    #  subprocess.check_call(cmd)
    #except subprocess.CalledProcessError:
    #  self.error("An error occured for command:\n%s"%cmd)
    #  sys.exit(1)

  def _call_mapper_pe(self, read1, read2, index, output):
    cmd = self._set_mapping_args_pe(read1, read2, index, output)
    self.info(" ".join(cmd))
    #try:
    #  subprocess.check_call(cmd)
    #except subprocess.CalledProcessError:
    #  self.error("An error occured for command:\n%s"%cmd)
    #  sys.exit(1)

  def _set_mapping_args_se(self, read, index, output):
    return ["-p", "value"]

  def _set_mapping_args_pe(self, read1, read2, index, output):
    return ["-p", "value"]

  def __str__(self):
    output = ""
    if not self.isPairEnd:
      output += " ".join(self._set_mapping_args_se(\
        self.file_list["read1CT"], self._indexCT, \
        self.file_list["wd_read1CT_to_CT"])) + "\n"
      output += " ".join(self._set_mapping_args_se(\
        self.file_list["read1CT"], self._indexGA, \
        self.file_list["wd_read1CT_to_GA"])) + "\n"
      if not self.isDirectional:
        output += " ".join(self._set_mapping_args_se(\
          self.file_list["read1GA"], self._indexCT, \
          self.file_list["wd_read1GA_to_CT"])) + "\n"
        output += " ".join(self._set_mapping_args_se(\
          self.file_list["read1GA"], self._indexGA, \
          self.file_list["wd_read1GA_to_GA"])) + "\n"
    else:
      output += " ".join(self._set_mapping_args_pe(\
        self.file_list["read1CT"], self.file_list["read2GA"], \
        self._indexCT, \
        self.file_list["wd_read1CTread2GA_to_CT"])) + "\n"
      output += " ".join(self._set_mapping_args_pe(\
        self.file_list["read1CT"], self.file_list["read2GA"], \
        self._indexGA, \
        self.file_list["wd_read1CTread2GA_to_GA"])) + "\n"
      if not self.isDirectional:
        output += " ".join(self._set_mapping_args_pe(\
          self.file_list["read1GA"], self.file_list["read2CT"], \
          self._indexCT, \
          self.file_list["wd_read1GAread2CT_to_CT"])) + "\n"
        output += " ".join(self._set_mapping_args_pe(\
          self.file_list["read1GA"], self.file_list["read2CT"], \
          self._indexGA, \
          self.file_list["wd_read1GAread2CT_to_GA"])) + "\n"

    return output
  
  def cleanup(self):
    """Delete temporary files.
    """
    return None

class MappingWithTophat2(Mapping):
  def _set_mapping_args_se(self, read, index, output):
    args = {\
    }

class MappingWithStar(Mapping):
  def _set_mapping_args_se(self, read, index, output):
    #prefix = re.match("(.+?\.)(Aligned.out.bam)?", output).group(1)
    args = {\
      "--runThreadN":"8", \
      "--readFilesIn":read, \
      "--genomeDir":index, \
      "--outFileNamePrefix":output + "/", \
      "--outFilterIntronMotifs":"RemoveNoncanonical", \
      #"--clip3pAdapterSeq":"", \
      "--outSAMstrandField":"intronMotif", \
      "--outSAMtype":"BAM Unsorted" \
    }
    cmd = [self.mapper]
    for i in args.keys():
      cmd.append(i)
      cmd.append(args[i])
    return cmd

  def _set_mapping_args_pe(self, read1, read2, index, output):
    #prefix = re.match("(.+?\.)(Aligned.out.bam)?", output).group(1)
    args = {\
      "--runThreadN":"8", \
      "--readFilesIn":" ".join((read1, read2)), \
      "--genomeDir":index, \
      "--outFileNamePrefix":output + "/", \
      "--outFilterIntronMotifs":"RemoveNoncanonical", \
      #"--clip3pAdapterSeq":"", \
      "--outSAMstrandField":"intronMotif", \
      "--outSAMtype":"BAM Unsorted" \
    }
    cmd = [self.mapper]
    for i in args.keys():
      cmd.append(i)
      cmd.append(args[i])
    return cmd
