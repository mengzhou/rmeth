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

def build_index_star(opt, files):
  """Build STAR index.
  """
  pars_ct = ["--runMode", "genomeGenerate",\
    "--runThreadN", "4",\
    "--genomeDir", opt.outdir_CT,\
    "--genomeFastaFiles", files["ctgenome"]\
    ]
  pars_ga = ["--runMode", "genomeGenerate",\
    "--runThreadN", "4",\
    "--genomeDir", opt.outdir_GA,\
    "--genomeFastaFiles", files["gagenome"]\
    ]
  args_ct = [opt.buildexe] + pars_ct
  args_ga = [opt.buildexe] + pars_ga

  opt.info("Building genome index of C to T converted genome...")
  opt.info("Command: %s"%" " .join(args_ct))
  try:
    os.chdir(opt.outdir_CT)
    subprocess.check_call(args_ct)
  except subprocess.CalledProcessError:
    opt.error("An error occured in building index. " + \
      "Please check your arguments:\n%s"%" ".join(args_ct))
    sys.exit(1)

  opt.info("Building genome index of G to A converted genome...")
  opt.info("Command: %s"%" " .join(args_ga))
  try:
    os.chdir(opt.outdir_GA)
    subprocess.check_call(args_ga)
  except subprocess.CalledProcessError:
    opt.error("An error occured in building index. " + \
      "Please check your arguments:\n%s"%" ".join(args_ga))
    sys.exit(1)

def build_index_bt2(opt, files):
  """Build bowtie2 index.
  """
  additional_pars = "-q"
  args_ct = [opt.buildexe, additional_pars, \
    files["ctgenome"], opt.outdir_CT + "/CTgenome"]
  args_ga = [opt.buildexe, additional_pars, \
    files["gagenome"], opt.outdir_GA + "/GAgenome"]

  opt.info("Building genome index of C to T converted genome...")
  opt.info("Command: %s"%" " .join(args_ct))
  try:
    subprocess.check_call(args_ct)
  except subprocess.CalledProcessError:
    opt.error("An error occured in building index. " + \
      "Please check your arguments:\n%s"%" ".join(args_ct))
    sys.exit(1)

  opt.info("Building genome index of G to A converted genome...")
  opt.info("Command: %s"%" " .join(args_ga))
  try:
    subprocess.check_call(args_ga)
  except subprocess.CalledProcessError:
    opt.error("An error occured in building index. " + \
      "Please check your arguments:\n%s"%" ".join(args_ga))
    sys.exit(1)

def bs_conversion(infh, outfh, C_to_T):
  """Convert FASTA file using bisulfite rules.
  """
  infh.seek(0)
  for l in infh:
    if l.startswith(">"):
      # cannot change chr name otherwise won't be able to use GTF annotation
      #if C_to_T:
      #  outfh.write(l.strip() + "_CT\n")
      #else:
      #  outfh.write(l.strip() + "_GA\n")
      outfh.write(l)
    else:
      l = l.upper()
      if C_to_T:
        outfh.write(l.replace("C","T"))
      else:
        outfh.write(l.replace("G","A"))

def init_files_paths(opt):
  """Initialize the list for files and file handlers.
  """
  opt.outdir = os.path.abspath(opt.outdir)
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

  opt.info("Performing bisulfite conversion for reference genome.")
  fhs["genome"] = open(files["genome"], 'r')
  if os.path.isfile(files["ctgenome"]):
    opt.warn("Converted genome %s exists and is skipped."%files["ctgenome"])
  else:
    fhs["ctgenome"] = open(files["ctgenome"], 'w')
    bs_conversion(fhs["genome"], fhs["ctgenome"], True)
  if os.path.isfile(files["gagenome"]):
    opt.warn("Converted genome %s exists and is skipped."%files["gagenome"])
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
  # Setup information format
  logging.basicConfig(level=20,
    format='[%(levelname)-5s][%(asctime)s] %(message)s ',
    datefmt='%H:%M:%S', stream=sys.stderr, filemode="w")
  opt.info = logging.info
  opt.error = logging.error
  opt.warn = logging.warn
  if not opt.buildexe and not opt.ref:
    parser.print_help()
    sys.exit(0)
  if not opt.buildexe:
    if which("bowtie2-build"):
      opt.buildexe = which("bowtie2-build")
      opt.warn("--build is not specified but bowtie2-build is found. "+\
        "%s will be used."%opt.buildexe)
      opt.isBowtie = True
    elif which("STAR"):
      opt.buildexe = which("STAR")
      opt.warn("--build is not specified but STAR is found. "+\
        "%s will be used."%opt.buildexe)
      opt.isBowtie = False
    else:
      opt.error("Must set the path to bowtie2-build executable "+\
        "or STAR executable with --build!")
      sys.exit(1)
  if not is_exe(opt.buildexe):
    opt.error(\
      "%s is not an executable file. Please check your path!"%opt.buildexe)
    sys.exit(1)
  else:
    if opt.buildexe.endswith("STAR"):
      opt.isBowtie = False
    elif opt.buildexe.endswith("bowtie2-build"):
      opt.isBowtie = True
    else:
      opt.error("Cannot determine which mapper %s belongs to. "%opt.buildexe)
      sys.exit(1)
  if not opt.ref:
    opt.error("Please specify one FASTA file using -r.")
    sys.exit(1)
  if not os.path.isfile(opt.ref):
    opt.error("%s does not exist or is not a file. "+\
      "Please check your path!"%opt.ref)
    sys.exit(1)
  opt.outdir = opt.outdir.rstrip("/")

  return opt

def main():
  usage = "Usage: %prog --build <path_to_build_program> -r <reference_fasta>"
  parser = OptionParser(usage=usage)
  parser.add_option("-b", "--build", action="store", type="string",
    dest="buildexe", help="Path to bowtie2-build executable program, or STAR executable.",
    metavar="<FILE>")
  parser.add_option("-r", "--reference", action="store", type="string",
    dest="ref", help="Path to reference genome FASTA file. Must be a single file.",\
    metavar="<FILE>")
  parser.add_option("-o","--output", action="store", type="string",
    dest="outdir", help="Path for the converted indices. Default: ./rmeth_genome",\
    default="./rmeth_genome")
  (opt, args) = parser.parse_args(sys.argv)

  # Option validation
  opt = opt_validation(parser, opt)

  # Initialize working directory and BS convert genome
  files = init_files_paths(opt)

  # Build bowtie2 indices for converted FASTA
  if opt.isBowtie:
    build_index_bt2(opt, files)
  else:
    build_index_star(opt, files)

  opt.info("Finished.")

if __name__ == '__main__':
  main()
