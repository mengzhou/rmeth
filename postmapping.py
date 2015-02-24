# postmapping.py: a module for post mapping processing of rmeth
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

class CIGAR:
  def __init__(self, str):
    pattern = re.compile("([0-9]+)([DIHMNPSX=])")
    self.cigar = pattern.findall(str)

  def apply_with_qual(self, seq, qual):
    new_seq = ""
    new_qual = ""
    index = 0
    op_copy = ['M', '=', 'X']
    op_skip = ['S', 'I', 'H', 'P']
    op_pad = ['N', 'D']
    for i in self.cigar:
      length = int(i[0])
      operation = i[1]
      if operation in op_copy:
        new_seq += seq[index:index+length]
        new_qual += qual[index:index+length]
        index += length
      elif operation in op_skip:
        index += length
      elif operation in op_pad:
        new_seq += "N"*length
        new_qual += "B"*length

    return (new_seq, new_qual)

class PostMapping:
  def __init__(self, opt, files):
    self.files = files
    self.info = opt.info
    self.warn = opt.warn
    self.error = opt.error
    self.mapper = opt.mapper
    self.samtools = opt.samtools
    self.isPairEnd = opt.isPairEnd

  def build(self):
    """Build mapped read file from mapping results.
    """
    if self.isPairEnd:
      self.__build_pe()
    else:
      self.__build_se()

  def __build_pe():
    return

  def __build_se():
    """For single-end reads, only C->T conversion was done
    to the FASTQ file; and the converted reads were mapped
    to C->T and G->A genomes respectively.
    Reads that were mapped to C->T genome should only map
    to + strand, while mapping to G->A genome should only
    yield - strand. So proper filtering is needed.
    """
    extract_args = [self.samtools, "view", "-F", "256", "-q", "10"
      , "-o", self.file.CT_read_to_CT_primary_sam()
      , self.file.CT_read_to_CT()]
    self.info("Extracting mapped BAM file %s..."\
      %self.file.CT_read_to_CT())
    try:
      subprocess.check_call(extract_args)
    except subprocess.CalledProcessError:
      opt.error("An error occured in samtools merging.")
      sys.exit(1)


def read_fastq_all_reads(infh):
  """Read FASTQ file and return a dictionary with read name and 
  corresponding read sequence.
  """
  read_name_sequence = {}
  counter = 0
  for l in infh:
    ind = counter%4
    if ind == 0:
      name = l[1:].strip().replace(" ", "_")
    elif ind == 1:
      read_name_sequence[name] = l.strip()
    counter += 1

  return read_name_sequence

def rev_comp( str ):
  comp = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
  str = str.upper()
  return "".join([comp[str[::-1][i]] for i in xrange(len(str))])

def main():
  in_fastq_fh = open(sys.argv[1], 'r')
  in_sam_fh = open(sys.argv[2], 'r')

  read_name_seq = read_fastq_all_reads(in_fastq_fh)
  output = sys.stdout
  read_pool = []
  prev_name = ""

  for l in in_sam_fh:
    if l.startswith("@"):
      pass
    else:
      f = l.split()
      if f[0] == prev_name or prev_name == "":
        read_pool.append(f)
        prev_name = f[0]
      else:
        prev_name = f[0]
        f = choice(read_pool)
        flag = int(f[1])
        if flag & 0x10:
          strand = '-'
          seq = rev_comp(read_name_seq[f[0]])
          qual = f[10][::-1]
        else:
          strand = '+'
          seq = read_name_seq[f[0]]
          qual = f[10]
        chr = f[2]
        start = int(f[3]) - 1
        cigar_op = CIGAR(f[5])
        new_seq, new_qual = cigar_op.apply_with_qual(seq, qual)
        if flag & 0x10:
          new_seq = rev_comp(new_seq)
        outl = "\t".join([chr, str(start), str(start + len(new_seq)), \
          f[0], '0', strand, new_seq, new_qual])
        output.write(outl + '\n')
        read_pool = [l.split()]

  f = choice(read_pool)
  name = f[0]
  flag = int(f[1])
  if flag & 0x10:
    strand = '-'
    seq = rev_comp(read_name_seq[name])
    qual = f[10][::-1]
  else:
    strand = '+'
    seq = read_name_seq[name]
    qual = f[10]
  chr = f[2]
  start = int(f[3]) - 1
  cigar_op = CIGAR(f[5])
  new_seq, new_qual = cigar_op.apply_with_qual(seq, qual)
  if flag & 0x10:
    new_seq = rev_comp(new_seq)
  outl = "\t".join([chr, str(start), str(start + len(new_seq)), \
    name, '0', strand, new_seq, new_qual])
  output.write(outl + '\n')
  in_fastq_fh.close()
  in_sam_fh.close()

if __name__ == "__main__":
  main()
