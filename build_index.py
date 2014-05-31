#!/usr/bin/env python
# build_index.py: part of rmeth RNA methylation analysis pipeline 
# which is used for building BS converted bowtie2 genome index.
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

"""Script for building converted genome index using bowtie2.
Will generate two subdirectories CTgenome and GAgenome, including
CT and GA converted genome respectively.
"""

import sys, os
import logging, subprocess
from optparse import OptionParser

def build_index(opt, files):
  """Build bowtie2 index.
  """
  additional_pars = "-q"
  args_ct = [opt.bowtie2, additional_pars, \
    files["ctgenome"], opt.outdir_CT + "/CTgenome"]
  args_ga = [opt.bowtie2, additional_pars, \
    files["gagenome"], opt.outdir_GA + "/GAgenome"]

  logging.info("Building genome index of C to T converted genome...")
  logging.info("Command: %s"%" " .join(args_ct))
  try:
    subprocess.check_call(args_ct)
  except subprocess.CalledProcessError:
    logging.error("An error occured in building index. " + \
      "Please check your arguments:\n%s"%" ".join(args_ct))
    sys.exit(1)

  logging.info("Building genome index of G to A converted genome...")
  logging.info("Command: %s"%" " .join(args_ga))
  try:
    subprocess.check_call(args_ga)
  except subprocess.CalledProcessError:
    logging.error("An error occured in building index. " + \
      "Please check your arguments:\n%s"%" ".join(args_ga))
    sys.exit(1)

def bs_conversion(infh, outfh, C_to_T):
  """Convert FASTA file using bisulfite rules.
  """
  infh.seek(0)
  for l in infh:
    if l.startswith(">"):
      if C_to_T:
        outfh.write(l.strip() + "_CT\n")
      else:
        outfh.write(l.strip() + "_GA\n")
    else:
      l = l.upper()
      if C_to_T:
        outfh.write(l.replace("C","T"))
      else:
        outfh.write(l.replace("G","A"))

def init_files_paths(opt):
  """Initialize the list for files and file handlers.
  """
  if not os.path.exists(opt.outdir):
    os.makedirs(opt.outdir)
  opt.outdir_CT = opt.outdir + "/CTgenome"
  opt.outdir_GA = opt.outdir + "/GAgenome"
  if not os.path.exists(opt.outdir_CT):
    os.makedirs(opt.outdir_CT)
  if not os.path.exists(opt.outdir_GA):
    os.makedirs(opt.outdir_GA)

  file_list = ["genome", "ctgenome", "gagenome"]
  files = {}.fromkeys(file_list)
  fhs = {}.fromkeys(file_list)

  files["genome"] = opt.ref
  files["ctgenome"] = opt.outdir_CT + "/CTgenome.fa"
  files["gagenome"] = opt.outdir_GA + "/GAgenome.fa"

  logging.info("Performing bisulfite conversion for reference genome.")
  fhs["genome"] = open(files["genome"], 'r')
  if os.path.isfile(files["ctgenome"]):
    logging.warn("Converted genome %s exists and is skipped."%files["ctgenome"])
  else:
    fhs["ctgenome"] = open(files["ctgenome"], 'w')
    bs_conversion(fhs["genome"], fhs["ctgenome"], True)
  if os.path.isfile(files["gagenome"]):
    logging.warn("Converted genome %s exists and is skipped."%files["gagenome"])
  else:
    fhs["gagenome"] = open(files["gagenome"], 'w')
    bs_conversion(fhs["genome"], fhs["gagenome"], False)

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
  if not opt.bowtie2 and not opt.ref:
    parser.print_help()
    sys.exit(0)
  if not opt.bowtie2:
    if which("bowtie2-build"):
      opt.bowtie2 = which("bowtie2-build")
      logging.warn("--bowtie2-build is not specified but found. %s will be used."\
        %opt.bowtie2)
    else:
      logging.error("Must set the path to bowtie2-build executable with --bowtie2-build!")
      sys.exit(1)
  if not is_exe(opt.bowtie2):
    logging.error(\
      "%s is not bowtie2-build executable file. Please check your path!"%opt.bowtie2)
    sys.exit(1)
  if not opt.ref:
    logging.error("Please specify one FASTA file using -r.")
    sys.exit(1)
  if not os.path.isfile(opt.ref):
    logging.error("%s does not exist. Please check your path!"%opt.ref)
    sys.exit(1)
  opt.outdir = opt.outdir.rstrip("/")

def main():
  usage = "Usage: %prog --bowtie2-build <path_to_bowtie2-build> -r <reference_fasta> -n <index_name>"
  parser = OptionParser(usage=usage)
  parser.add_option("-b", "--bowtie2-build", action="store", type="string",
    dest="bowtie2", help="Path to bowtie2-build executable program.",
    metavar="<FILE>")
  parser.add_option("-r", "--reference", action="store", type="string",
    dest="ref", help="Path to reference genome FASTA file. Must be a single file.",\
    metavar="<FILE>")
  parser.add_option("-o","--output", action="store", type="string",
    dest="outdir", help="Path for the converted indices. Default: ./rmeth_genome",\
    default="./rmeth_genome")
  (opt, args) = parser.parse_args(sys.argv)

  # Setup information format
  logging.basicConfig(level=20,
    format='[%(levelname)-5s][%(asctime)s] %(message)s ',
    datefmt='%H:%M:%S', stream=sys.stderr, filemode="w")

  # Option validation
  opt_validation(parser, opt)

  # Initialize working directory and BS convert genome
  files = init_files_paths(opt)

  # Build bowtie2 indices for converted FASTA
  build_index(opt, files)

  logging.info("Finished.")

if __name__ == '__main__':
  main()
