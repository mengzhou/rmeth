#!/usr/bin/env python
# rmeth.py: an RNA methylation analysis pipeline which uses tophat
# for split read mapping.
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

"""Rmeth is a pipeline for RNA methylation mapping and analysis.
"""

import sys, os, re
import logging, subprocess
from optparse import OptionParser
from random import choice

class Files:
  def __init__(opt):
    self._wd = opt.outdir
    self._wd_read1CT_to_CT = opt.mapped_dir_read1CT_to_CT
    self._wd_read1CT_to_GA = opt.mapped_dir_read1CT_to_GA
    name1 = re.match("(.+)(\.fastq|\.fq)", opt.read1).group(1)
    # original fastq files
    self._fn_read1 = opt.read1
    self._fn_read2 = opt.read2
    # converted fastq files
    self._fn_read1CT = self._wd + "/%s_CT.fastq"%name1
    # original mapped bam files
    self._fn_read1CT_to_CT = self._wd + opt.mapped_dir_read1CT_to_CT \
      + ".Aligned.out.bam"
    self._fn_read1CT_to_GA = self._wd + opt.mapped_dir_read1CT_to_GA \
      + ".Aligned.out.bam"
    # bam files with only primary mapping
    self._fn_read1CT_to_CT_primary = re.match("(.+)\.bam", \
      self._fn_read1CT_to_CT).group(1) + "primary.bam"
    self._fn_read1CT_to_GA_primary = re.match("(.+)\.bam", \
      self._fn_read1CT_to_GA).group(1) + "primary.bam"
    # bam files with compatible strand
    self._fn_read1CT_to_CT_compatible = re.match("(.+)\.bam", \
      self._fn_read1CT_to_CT).group(1) + "compatible.bam"
    self._fn_read1CT_to_GA_compatible = re.match("(.+)\.bam", \
      self._fn_read1CT_to_GA).group(1) + "compatible.bam"

    if not opt.directional:
      # converted fastq files
      self._fn_read1GA = self._wd + "/%s_GA.fastq"%name1
      # original mapped bam files
      self._fn_read1GA_to_CT = self._wd + opt.mapped_dir_read1GA_to_CT \
        + ".Aligned.out.bam"
      self._fn_read1GA_to_GA = self._wd + opt.mapped_dir_read1GA_to_GA \
        + ".Aligned.out.bam"
      # bam files with only primary mapping
      self._fn_read1GA_to_CT_primary = re.match("(.+)\.bam", \
        self._fn_read1GA_to_CT).group(1) + "primary.bam"
      self._fn_read1GA_to_GA_primary = re.match("(.+)\.bam", \
        self._fn_read1GA_to_GA).group(1) + "primary.bam"
      # bam files with compatible strand
      self._fn_read1GA_to_CT_compatible = re.match("(.+)\.bam", \
        self._fn_read1GA_to_CT).group(1) + "compatible.bam"
      self._fn_read1GA_to_GA_compatible = re.match("(.+)\.bam", \
        self._fn_read1GA_to_GA).group(1) + "compatible.bam"

    # merged bam file
    self._fn_read1_merged_bam = self._wd + "/%s.merged.bam"%(name1)

    if opt.isPairEnd:
      name2 = re.match("(.+)(\.fastq|\.fq)", opt.read2).group(1)
      # converted fastq files
      self._fn_read2GA = self._wd + "/%s_CT.fastq"%name2
      # original mapped bam files
      self._fn_read2GA_to_CT = self._wd + opt.mapped_dir_read2GA_to_CT \
        + ".Aligned.out.bam"
      self._fn_read2GA_to_GA = self._wd + opt.mapped_dir_read2GA_to_GA \
        + ".Aligned.out.bam"
      # bam files with only primary mapping
      self._fn_read2GA_to_CT_primary = re.match("(.+)\.bam", \
        self._fn_read2GA_to_CT).group(1) + "primary.bam"
      self._fn_read2GA_to_GA_primary = re.match("(.+)\.bam", \
        self._fn_read2GA_to_GA).group(1) + "primary.bam"
      # bam files with compatible strand
      self._fn_read2GA_to_CT_compatible = re.match("(.+)\.bam", \
        self._fn_read2GA_to_CT).group(1) + "compatible.bam"
      self._fn_read2GA_to_GA_compatible = re.match("(.+)\.bam", \
        self._fn_read2GA_to_GA).group(1) + "compatible.bam"

      # merged bam file
      self._fn_read2CT_merged_bam = self._wd + "/%s.merged.bam"%(name2)

      if not opt.directional:
        # converted fastq files
        self._fn_read2CT = self._wd + "/%s_GA.fastq"%name1
        # original mapped bam files
        self._fn_read2CT_to_CT = self._wd + opt.mapped_dir_read2CT_to_CT \
          + ".Aligned.out.bam"
        self._fn_read2CT_to_GA = self._wd + opt.mapped_dir_read2CT_to_GA \
          + ".Aligned.out.bam"
        # bam files with only primary mapping
        self._fn_read2CT_to_CT_primary = re.match("(.+)\.bam", \
          self._fn_read2CT_to_CT).group(1) + "primary.bam"
        self._fn_read2CT_to_GA_primary = re.match("(.+)\.bam", \
          self._fn_read2CT_to_GA).group(1) + "primary.bam"
        # bam files with compatible strand
        self._fn_read2CT_to_CT_compatible = re.match("(.+)\.bam", \
          self._fn_read2CT_to_CT).group(1) + "compatible.bam"
        self._fn_read2CT_to_GA_compatible = re.match("(.+)\.bam", \
          self._fn_read2CT_to_GA).group(1) + "compatible.bam"

      # merged bam file
      self._fn_read1_merged_bam = self._wd + "/%s.merged.bam"%(name1)

    # final output .mr
    self._fn_mr = self._wd + "/%s.mr"%name1

  def get_file_list(self):
    fn_list = {"read1":self._fn_read1, "read1CT":self._fn_read1CT, \
      "read1CT_to_CT":self._fn_read1CT_to_CT, \
      "read1CT_to_GA":self._fn_read1CT_to_GA, \
      "read1CT_to_CT_primary":self._fn_read1CT_to_CT_primary, \
      "read1CT_to_GA_primary":self._fn_read1CT_to_GA_primary, \
      "read1CT_to_CT_compatible":self._fn_read1CT_to_CT_compatible, \
      "read1CT_to_GA_compatible":self._fn_read1CT_to_GA_compatible, \
      "read1GA":self._fn_read1GA, \
      "read1GA_to_CT":self._fn_read1GA_to_CT, \
      "read1GA_to_GA":self._fn_read1GA_to_GA, \
      "read1GA_to_CT_primary":self._fn_read1GA_to_CT_primary, \
      "read1GA_to_GA_primary":self._fn_read1GA_to_GA_primary, \
      "read1GA_to_CT_compatible":self._fn_read1GA_to_CT_compatible, \
      "read1GA_to_GA_compatible":self._fn_read1GA_to_GA_compatible, \
      "read1_merged_bam":self._fn_read1_merged_bam,\
      "read2":self._fn_read2, "read2CT":self._fn_read2CT, \
      "read2GA":self._fn_read2GA, \
      "read2GA_to_CT":self._fn_read2GA_to_CT, \
      "read2GA_to_GA":self._fn_read2GA_to_GA, \
      "read2GA_to_CT_primary":self._fn_read2GA_to_CT_primary, \
      "read2GA_to_GA_primary":self._fn_read2GA_to_GA_primary, \
      "read2GA_to_CT_compatible":self._fn_read2GA_to_CT_compatible, \
      "read2GA_to_GA_compatible":self._fn_read2GA_to_GA_compatible, \
      "read2CT_to_CT":self._fn_read2CT_to_CT, \
      "read2CT_to_GA":self._fn_read2CT_to_GA, \
      "read2CT_to_CT_primary":self._fn_read2CT_to_CT_primary, \
      "read2CT_to_GA_primary":self._fn_read2CT_to_GA_primary, \
      "read2CT_to_CT_compatible":self._fn_read2CT_to_CT_compatible, \
      "read2CT_to_GA_compatible":self._fn_read2CT_to_GA_compatible, \
      "read2_merged_bam":self._fn_read2_merged_bam, \
      "mr":self._fn_mr
      }
    return fn_list

def read_fastq_one_read(infh):
  """Read FASTQ file and return read name and read sequence.
  """
  l = infh.readline()
  if not l:
    return (None, None)
  if not l.startswith("@"):
    opt.error("Corrupted FASTQ file!")
    sys.exit(1)
  else:
    name = l[1:].strip()

  l = infh.readline()
  seq = l.strip()
  # these two lines for quality score is not needed
  l = infh.readline()
  l = infh.readline()
  return (name, seq)

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

def SAM_isPairend_read1(flag):
  return flag & 0x1 and flag & 0x40

def SAM_isPairend_read2(flag):
  return flag & 0x1 and flag & 0x80

def filter_multi_mapping( list_of_candidates ):
  """Choose one from all ambiguously mapped reads.
  Now just randomly choose.
  """
  return choice(list_of_candidates)

def sam_to_mr( sam_fields ):
  """Convert SAM format to MappedRead. Input is list of SAM fields.
  """

def replace_sam_sequence_se(inf_sam, read_name_seq):
  """Replace read sequence of input SAM file by their original ones
  in FASTQ.
  """
  samfh = open(inf_sam, 'r')
  read_pool = []
  last_name = ""
  name = ""
  l = samfh.readline()
  while l:
    if l.startswith("@"):
      pass
    else:
      f = l.split("\t")
      name = f[0]
      if name != last_name:
        best_read = filter_multi_mapping(read_pool)
        best_read[9] = read_name_seq.pop(last_name)
        sam_to_mr(best_read)
        read_pool = f
      else:
        read_pool.append(f)
    last_name = name
    l = samfh.readline()

  sam_to_mr(read_pool, read_name_seq.pop(name))

def replace_sam_sequence_pe(inf_sam, read_name_seq, read_name_seq_read2):
  """Replace read sequence of input SAM file by their original ones
  in FASTQ.
  """
  samfh = open(inf_sam, 'r')
  fastqfh_read1 = open(inf_fastq_read1, 'r')
  fastqfh_read2 = open(inf_fastq_read2, 'r')
  (fastq_read1_name, fastq_read1_seq) = read_fastq_one_read(fastqfh_read1)
  fastq_name_read1 = fastq_read_name_process(fastq_name_read1, 1)
  (fastq_read2_name, fastq_read2_seq) = read_fastq_one_read(fastqfh_read2)
  fastq_name_read2 = fastq_read_name_process(fastq_name_read2, 2)

  for l in samfh:
    if l.startswith("@"):
      pass
    else:
      f = l.split("\t")
      name = f[0]
      seq = f[9]
      flag = int(f[1])
      while fastq_name and fastq_name != name:
        (fastq_name, fastq_seq) = read_fastq_one_read(fastqfh)
        fastq_name = fastq_read_name_process(fastq_name, mate)

def construct_mapped_reads(opt, files):
  """Filter mapped reads (BAM file) using preset parameters in samtools.
  Then replace the converted sequences with their original ones from FASTQ.
  """
  # set temporary file names. It seems that python does not provide a
  # good solution for real-time pipe handling, so I'm using temp file

  # 1. merge the mapped files
  files["mergedBAM"] = opt.outdir + "/merged.bam.tmp%d"%os.getpid()
  samtools_args = [opt.samtools, "merge", files["mergedBAM"],\
      files["mappedCT"], files["mappedGA"]]
  try:
    subprocess.check_call(samtools_args)
  except subprocess.CalledProcessError:
    opt.error("An error occured in samtools merging.")
    sys.exit(1)

  # 2. filter the merged file
  # -q INT: minimum mapping quality, set to 10 for < 0.1 error rate
  # -F INT: set to 256 to filter out secondary alignment specifically
  files["filteredSAM"] = opt.outdir + "/filtered.sam.tmp%d"%os.getpid()
  samtools_args = [opt.samtools, "view", "-q", "10", \
    "-F", "256", "-o", files["filteredSAM"], files["mergedBAM"]]
  try:
    subprocess.check_call(samtools_args)
  except subprocess.CalledProcessError:
    opt.error("An error occured in samtools filtering.")
    sys.exit(1)

  opt.info("Filtering mapped BAM file %s with samtools."%files["mappedCT"])
  try:
    subprocess.check_call(samtools_args_ct)
  except subprocess.CalledProcessError:
    opt.error("An error occured in samtools loading %s."%files["mappedCT"])
    sys.exit(1)
  os.remove(files["mergedBAM"])

  opt.info("Loading FASTQ file %s..."%files["read1"])
  read_name_seq = read_fastq_all_reads(files["read1"])
  if opt.isPairEnd:
    opt.info("Loading FASTQ file %s..."%files["read1"])
    read_name_seq_read2 = read_fastq_all_reads(files["read2"])

  opt.info("Processing mapped reads"%files["filteredSAM"])
  if opt.isPairEnd:
    replace_sam_sequence_pe(files["filteredSAM"], \
      read_name_seq, read_name_seq_read2)
  else:
    replace_sam_sequence_se(files["filteredSAM"], read_name_seq)

def get_tophat_args(opt, isCT):
  """Parse tophat parameters to a list that will be used for
  subprocess system call.
  --library-type will be added because the mapping will require
  bisulfite converted reads to be mapped to only one strand.
  """
  if not opt.pars:
    opt.pars = ""
  if opt.pars.find("--library-type") == -1:
    l = opt.pars.split()
    if isCT:
      l += ["--library-type", "fr-firststrand"]
    else:
      l += ["--library-type", "fr-secondstrand"]

  if opt.pars.find("-o") == -1 or opt.pars.find("--output-dir") == -1:
    if isCT:
      l += ["-o", opt.mapper_dir_CT]
    else:
      l += ["-o", opt.mapper_dir_GA]
  else:
    opt.error("Please do not specify -o for tophat parameter.")
    sys.exit(1)

  return l

def map_reads(opt, files):
  """Map reads using tophat.
  """
  args_ct = [opt.mapper]
  args_ct += get_tophat_args(opt, True)
  args_ct.append(opt.ctindex)
  args_ga = [opt.mapper]
  args_ga += get_tophat_args(opt, False)
  args_ga.append(opt.gaindex)
  if opt.isPairEnd:
    args_reads_ct = [files["read1CT"], files["read2GA"]]
    args_reads_ga = [files["read1GA"], files["read2CT"]]
  else:
    args_reads_ct = [files["read1CT"]]
    args_reads_ga = [files["read1GA"]]
  args_ct += args_reads_ct
  args_ga += args_reads_ga

  opt.info("Mapping reads to C to T converted genome...")
  opt.info("Command: %s"%" ".join(args_ct))
  try:
    subprocess.check_call(args_ct)
  except subprocess.CalledProcessError:
    opt.error("An error occured in tophat mapping. " + 
      "Please check your arguments:\n%s"%" ".join(args_ct))
    sys.exit(1)

  opt.info("Mapping reads to G to A converted genome...")
  opt.info("Command: %s"%" ".join(args_ga))
  try:
    subprocess.check_call(args_ga)
  except subprocess.CalledProcessError:
    opt.error("An error occured in tophat mapping. " +
      "Please check your arguments:\n%s"%" ".join(args_ga))
    sys.exit(1)

  # Try to locate mapped read files after tophat is done
  if os.path.isfile(opt.mapper_dir_CT + "/accepted_hits.bam"):
    files["mappedCT"] = opt.mapper_dir_CT + "/accepted_hits.bam"
  else:
    opt.error("No mapped file found for C to T converted genome" +
        "Some error might have occured during tophat run.")
    sys.exit(1)
  if os.path.isfile(opt.mapper_dir_GA + "/accepted_hits.bam"):
    files["mappedGA"] = opt.mapper_dir_GA + "/accepted_hits.bam"
  else:
    opt.error("No mapped file found for G to A converted genome" +
        "Some error might have occured during tophat run.")
    sys.exit(1)

def check_genome_index(opt):
  """Check if the bowtie2 genome indices for tophat are ready.
  Will check if CT and GA converted genomes exist.
  """
  opt.ctindex = opt.index.rstrip("/") + "/CTgenome/CTgenome"
  opt.gaindex = opt.index.rstrip("/") + "/GAgenome/GAgenome"

  if not os.path.isfile(opt.ctindex + ".1.bt2"):
    opt.error("C to T genome index is not found. Please check -i parameter!")
    opt.error("-i <dir> must be a directory including C_to_T and G_to_A directories.")
    sys.exit(1)
  if not os.path.isfile(opt.gaindex + ".1.bt2"):
    opt.error("G to A genome index is not found. Please check -i parameter!")
    opt.error("-i <dir> must be a directory including C_to_T and G_to_A directories.")
    sys.exit(1)

def fastq_read_name_process(name, mate):
  """Process read name in FASTQ files. Main purpose is to deal with spaces and
  add mate tag in pair-end data.
  """
  return name.replace(" ", "_")

def bs_conversion(infh, outfh, C_to_T, mate):
  """Convert FASTQ files using bisulfite rules.
  """
  infh.seek(0)
  counter = 0
  for l in infh:
    ind = counter%4
    if ind == 0:
      # It seems that read pairs without unique siffices also works. So why add them?
      #outfh.write(l.strip() + "/%d\n"%mate)
      # I think it is necessary to remove spaces in read name
      outfh.write(fastq_read_name_process(l, mate))
    elif ind == 1:
      if C_to_T:
        outfh.write(l.replace("C","T"))
      else:
        outfh.write(l.replace("G","A"))
    else:
      outfh.write(l)
    counter += 1

def init_files_paths(opt):
  """Initialize the list for files and file handlers.
  """
  if not os.path.exists(opt.outdir):
    os.makedirs(opt.outdir)
  opt.mapped_dir_read1CT_to_CT = opt.outdir + "/read1CT_to_CTgenome"
  opt.mapped_dir_read1CT_to_GA = opt.outdir + "/read1CT_to_GAgenome"
  os.makedirs(opt.mapped_dir_read1CT_to_CT)
  os.makedirs(opt.mapped_dir_read1CT_to_GA)
  if not opt.directional:
    opt.mapped_dir_read1GA_to_CT = opt.outdir + "/read1GA_to_CTgenome"
    opt.mapped_dir_read1GA_to_GA = opt.outdir + "/read1GA_to_GAgenome"
    os.makedirs(opt.mapped_dir_read1GA_to_CT)
    os.makedirs(opt.mapped_dir_read1GA_to_GA)
  if opt.isPairEnd:
    opt.mapped_dir_read2GA_to_CT = opt.outdir + "/read2GA_to_CTgenome"
    opt.mapped_dir_read2GA_to_GA = opt.outdir + "/read2GA_to_GAgenome"
    os.makedirs(opt.mapped_dir_read2GA_to_CT)
    os.makedirs(opt.mapped_dir_read2GA_to_GA)
    if not opt.directional:
      opt.mapped_dir_read2CT_to_CT = opt.outdir + "/read2CT_to_CTgenome"
      opt.mapped_dir_read2CT_to_GA = opt.outdir + "/read2CT_to_GAgenome"
      os.makedirs(opt.mapped_dir_read2CT_to_CT)
      os.makedirs(opt.mapped_dir_read2CT_to_GA)

  file_list = ["read1", "read2", "read1CT", "read1GA", "read2CT", "read2GA"]
  files = {}.fromkeys(file_list)
  fhs = {}.fromkeys(file_list)

  files["read1"] = opt.read1
  files["read2"] = opt.read2
  files["read1CT"] = re.sub("\.fastq|\.fq","",opt.read1) + "_CT.fastq"
  files["read1GA"] = re.sub("\.fastq|\.fq","",opt.read1) + "_GA.fastq"
  if opt.isPairEnd:
    files["read2CT"] = re.sub("\.fastq|\.fq","",opt.read2) + "_CT.fastq"
    files["read2GA"] = re.sub("\.fastq|\.fq","",opt.read2) + "_GA.fastq"
  try:
    fhs["read1"] = open(files["read1"], 'r')
  except:
    opt.error("Failed to open FASTQ file: %s."\
        %files["read1"])
    sys.exit(1)
  if opt.isPairEnd:
    try:
      fhs["read2"] = open(files["read2"], 'r')
    except:
      opt.error("Failed to open FASTQ file: %s."\
          %files["read2"])
      sys.exit(1)

  opt.info("Performing read conversion for mapping.")
  if os.path.isfile(files["read1CT"]):
    opt.warn("%s exists and is skipped."%files["read1CT"])
  else:
    fhs["read1CT"] = open(files["read1CT"], 'w')
    bs_conversion(fhs["read1"], fhs["read1CT"], True, 1)
  if os.path.isfile(files["read1GA"]):
    opt.warn("%s exists and is skipped."%files["read1GA"])
  else:
    fhs["read1GA"] = open(files["read1GA"], 'w')
    bs_conversion(fhs["read1"], fhs["read1GA"], False, 1)

  if opt.isPairEnd:
    if os.path.isfile(files["read2CT"]):
      opt.warn("%s exists and is skipped."%files["read2CT"])
    else:
      fhs["read2CT"] = open(files["read2CT"], 'w')
      bs_conversion(fhs["read2"], fhs["read2CT"], True, 2)
    if os.path.isfile(files["read2GA"]):
      opt.warn("%s exists and is skipped."%files["read2GA"])
    else:
      fhs["read2GA"] = open(files["read2GA"], 'w')
      bs_conversion(fhs["read2"], fhs["read2GA"], False, 2)

  for fh in fhs.values():
    if fh:
      fh.close()

  return files

def is_exe(file):
  return os.path.isfile(file) and os.access(file, os.X_OK)

def which(program):
  """Do the same thing as linux "which" command.
  """
  for path in os.environ["PATH"].split(os.pathsep):
    path = path.strip('"')
    exe_file = os.path.join(path, program)
    if is_exe(exe_file):
      return exe_file

  return None

def opt_validation(parser, opt):
  if not opt.read1 and not opt.mapper and not opt.samtools:
    parser.print_help()
    sys.exit(0)
  opt.info = opt.info
  opt.warn = opt.warn
  opt.error = opt.error
  if not opt.mapper:
    if which("tophat2"):
      opt.mapper = which("tophat2")
      opt.warn("--mapper is not specified but tophat2 is found. %s will be used."\
        %opt.mapper)
    else:
      if which("STAR"):
        opt.mapper = which("STAR")
        opt.warn("--mapper is not specified but STAR is found. %s will be used."\
          %opt.mapper)
      else:
        opt.error("Must set the path to mapper executable with --mapper!")
        sys.exit(1)
  if not is_exe(opt.mapper):
    opt.error(\
      "%s is not an executable file. Please check your path!"%opt.mapper)
    sys.exit(1)
  if not opt.samtools:
    if which("samtools"):
      opt.samtools = which("samtools")
      opt.warn("--samtools is not specified but found. %s will be used."\
        %opt.samtools)
    else:
      opt.error("Must set the path to samtools executable with --samtools!")
      sys.exit(1)
  if not is_exe(opt.samtools):
    opt.error(\
      "%s is not samtools executable file. Please check your path!"%opt.samtools)
    sys.exit(1)
  if not opt.read1:
    opt.error("Please specify at least one FASTQ file using -1.")
    sys.exit(1)
  if not opt.read2:
    opt.isPairEnd = False
  else:
    opt.isPairEnd = True

def main():
  usage = "Usage: %prog [-p \"tophat options\"] --tophat <path to tophat> " + \
      "--samtools <path to samtools> " + \
      "-i <bowtie2 index> -1 <read1> [-2 read2] [other options]"
  parser = OptionParser(usage=usage)
  parser.add_option("-o", "--output-dir", action="store", type="string", \
    default="%s/rmeth"%os.getcwd(), dest="outdir", \
    help="Output directory. Default: ./rmeth", metavar="<string>")
  parser.add_option("-p", "--parameters", action="store", type="string",
    dest="pars", help="Parameters for mapper." + \
      "Must be quoted because space is included.", \
    metavar="\"<parameter1 parameter2 ...>\"")
  parser.add_option("-1", "--read1", action="store", type="string",
    dest="read1", metavar="<FILE>",\
    help="FASTQ file for read. Mate 1 for pair-end or single-end file.")
  parser.add_option("-2", "--read2", action="store", type="string",
    dest="read2", metavar="<FILE>", \
    help="FASTQ file for read. Mate 2 for pair-end. Not neede for single-end.")
  parser.add_option("-i", "--index", action="store", type="string",
    dest="index", metavar="<DIR>",\
    help="Path to genome index directory (generated by rmeth) to be used in mapping.")
  parser.add_option("--mapper", action="store", type="string",
    dest="mapper", help="Path to tophat2/STAR executable program.", metavar="<FILE>")
  parser.add_option("--samtools", action="store", type="string",
    dest="samtools", help="Path to samtools.")
  parser.add_option("-d", "--directional", action="store_true",
    dest="directional", help="Library is directional. Default: false.") 
  (opt, args) = parser.parse_args(sys.argv)

  # Setup information format
  logging.basicConfig(level=20,
    format='[%(levelname)-5s][%(asctime)s] %(message)s ',
    datefmt='%H:%M:%S', stream=sys.stderr, filemode="w")

  # Option validation
  opt_validation(parser, opt)

  # Initialization of files names, file handlers, temporary directory,
  # option check, path check, etc. Also do BS convert for reads.
  check_genome_index(opt)
  opt.info("Output directory is set to %s"%opt.outdir)
  files = init_files_paths(opt)

  # Map reads using tophat
  map_reads(opt, files)

  # Things below need to be done with C++
  #######################################
  # Filter mapped reads with some quality threshold, and replace the
  # converted sequences to their original ones from FASTQ
  #construct_mapped_reads(opt, files)

  # Load original FASTQ file and the genome to determine methylation
  #GetMethylationLevel()
  #######################################

  # Output and cleanup

if __name__ == '__main__':
  main()
