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
    for i in cigar:
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

  def apply_with_qual_junction(self, seq, qual):
    new_seq = ""
    new_qual = ""
    index = 0
    junction = [[], []]
    op_copy = ['M', '=', 'X']
    op_skip = ['S', 'I', 'H', 'P']
    op_pad = ['D']
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
      else:
        # N for junction
        junction[0].append([new_seq, new_qual])
        junction[1].append(length)
        new_seq = ""
        new_qual = ""

    junction[0].append([new_seq, new_qual])

    return junction

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
  in_sam_fh = open(sys.argv[1], 'r')
  #in_fastq_fh = open(sys.argv[2], 'r')

  #read_name_seq = read_fastq_all_reads(in_fastq_fh)
  output = sys.stdout
  read_pool = []
  prev_name = ""

  for l in in_sam_fh:
    if l.startswith("@"):
      continue
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
          #seq = rev_comp(read_name_seq[f[0]])
          #seq = rev_comp(f[-1])
          seq = rev_comp(f[9])
          qual = f[10][::-1]
        else:
          strand = '+'
          #seq = read_name_seq[f[0]]
          #seq = f[-1]
          seq = f[9]
          qual = f[10]
        chr = f[2]
        start = int(f[3]) - 1
        cigar_op = CIGAR(f[5])
        junction = cigar_op.apply_with_qual_junction(seq, qual)
        if flag & 0x10:
          seq_order = range(len(junction[0])-1, -1, -1)
        else:
          seq_order = range(len(junction[0]))
        for i in seq_order:
          if len(junction[0]) > 1:
            name = f[0] + "_junc" + str(i)
          else:
            name = f[0]
          split_read = junction[0][i]
          if flag & 0x10:
            new_seq = rev_comp(split_read[0])
            new_qual = split_read[1][::-1]
            if i < len(seq_order)-1:
              new_start = start + junction[1][i]
            else:
              new_start = start
          else:
            new_seq = split_read[0]
            new_qual = split_read[1]
            if i > 0:
              new_start = start + junction[1][i-1]
            else:
              new_start = start
          outl = "\t".join([chr, str(new_start) \
            , str(new_start + len(new_seq)), \
            name, '0', strand, new_seq, new_qual])
          output.write(outl + '\n')

        read_pool = [l.split()]

  f = choice(read_pool)
  name = f[0]
  flag = int(f[1])
  if flag & 0x10:
    strand = '-'
    #seq = rev_comp(read_name_seq[name])
    #seq = rev_comp(f[-1])
    seq = rev_comp(f[9])
    qual = f[10][::-1]
  else:
    strand = '+'
    #seq = read_name_seq[name]
    #seq = f[-1]
    seq = f[9]
    qual = f[10]
  chr = f[2]
  start = int(f[3]) - 1
  cigar_op = CIGAR(f[5])
  junction = cigar_op.apply_with_qual_junction(seq, qual)
  if flag & 0x10:
    seq_order = range(len(junction[0])-1, -1, -1)
  else:
    seq_order = range(len(junction[0]))
  for i in seq_order:
    if len(junction[0]) > 1:
      name = f[0] + "_junc" + str(i)
    else:
      name = f[0]
    split_read = junction[0][i]
    if flag & 0x10:
      new_seq = rev_comp(split_read[0])
      new_qual = split_read[1][::-1]
      if i < len(seq_order)-1:
        new_start = start + junction[1][i]
      else:
        new_start = start
    else:
      new_seq = split_read[0]
      new_qual = split_read[1]
      if i > 0:
        new_start = start + junction[1][i-1]
      else:
        new_start = start
    outl = "\t".join([chr, str(new_start) \
      , str(new_start + len(new_seq)), \
      name, '0', strand, new_seq, new_qual])
    output.write(outl + '\n')

  read_pool = [l.split()]

  #in_fastq_fh.close()
  in_sam_fh.close()

if __name__ == "__main__":
  main()
