#!/usr/bin/python
import sys, os, getopt
import subprocess
import logging

'''
This program use vcflib to do the vcf file quallity control
and use annovar to do filter based variants annotation based on
differenct variant databases. It assume that all the vcf (uncompressed)
files are sitting under the input directory. The output files will
go to output directory
'''

#get the running script name as myapp
myapp, ext = os.path.splitext(__file__)

#configure the on screen message to the user during the program running
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s  %(message)s',
                    datefmt='%d %b %Y %H:%M:%S'
                 )
                    
def parseArg(argv):
   
   '''
   Configure the directory which holds all the input vcf files.
   It assmes all the vcf files are uncompressed already. 
   By default, it will look for a direcoty "vcfinfput" in your 
   home direcory. you can change it to the other place here 
   or configure the path through the command line option 
   -i <vcf file input directory>
   '''
   inputdir = '~/vcfinput'
   
   '''
   Configure the directory which holds all the output files.
   By default, it will look for a direcoty "output" in your
   home direcory. you can change it to the other place here
   or configure the path through the command line option
   -o <results output directory>
   '''
   outputdir = '~/output'

   '''
   Configure the path for finding vcffileter program. By default,
   it assumes the vcflib program in your home direcory. you can
   change it to the other place here or configure the path through
   the command line option -vf <vcffilter path>
   '''
   vcffilter = '~/vcflib/bin/vcffilter'
   
   '''
   Configure the path for finding annovar program directory. By default,
   it assumes the annovar program in your home direcory. you can
   change it to the other place here or configure the path through
   the command line options -a <annovar program directory>
   '''
   annovar = '~/annovar'
   
   '''
   Configure the path for finding annovar variation datbase directory.
   By default, it assumes databse will be sitting in your in your home
   direcory ~/humandb. you can change it to the other place here or 
   configure the path through the command line option -d <variation database direcotry>
   '''
   database_dir = '~/humandb'
   
   usage = """
   python varPrioritisation.py
   -i --idir    vcf files input direcotry. All the uncompressed vcf input
                files are sitting in this direcotry. By default, this
                has been configured to ~/vcfinput
                
   -o --odir    output directory. All the ouput files are writting into
                this directory. By default, it will write to ~/output
   
   -f          vcffilter program path. By default, it will look for the
                program in: ~/vcflib/bin/vcffilter
   -a --adir    annovar program direcotry. By default, it will look for
                the direcotry in: ~/annovar
   -d           annovar variation databse direcotry you have downloaded.
                By default, it will look for the database in: ~/humandb
   Example:
     python varPrioritisation.py -i <vcf file direcotry> -o <result output direcotory> -f <vcffilter program path> -a <annovar program direcotry> -d <downloaded annovar variation database direcotry>
   
   """

   try:
      opts, args = getopt.getopt(argv,"hi:o:f:a:d:",["idir=","odir=", "adir="])
   except getopt.GetoptError:
      print usage
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print usage
         sys.exit()
      elif opt in ("-i", "--idir"):
         inputdir = arg
      elif opt in ("-o", "--odir"):
         outputdir = arg
      elif opt == '-f':
         vcffilter = arg
      elif opt in ("-a", "--adir"):
         annovar = arg
      elif opt == '-d':
         database_dir = arg

         
   return inputdir, outputdir, vcffilter, annovar, database_dir

def is_exe(fpath):

   #return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
   return os.path.isfile(fpath)

def run_tool(cmd):
   """Perform the system call to the external tool"""

   logging.info('Start running: ' + cmd)
   p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
   output,error = p.communicate()
   if p.returncode:
      #print "Error: ",error.strip()
      raise Exception(errors)
   logging.info('Finish running: ' + cmd)

def qualityFilter(ifile, ofile, vcffilter):
   """Use vcflib to do the quality control"""
   
   #check whether vcffilter program exists or not
   if not is_exe(vcffilter):
      print vcffilter  + ' program does not exits'
      exit(1)
   
   #Filter by depth of coverage
   cmd = vcffilter + ' -f "DP > 10" ' + ifile + '> ' + ofile
   run_tool(cmd)


def vcf2AnnovarInput(ifile, ofile, annovar):
   """Covert vcf file to annovar input format"""

   annovar_convert = annovar + os.sep + 'convert2annovar.pl'
   #check whether convert2annovar.pl  program exists or not
   if not is_exe(annovar_convert):
      print annovar_convert  + ' program does not exits'
      exit(1)

   cmd = 'perl ' + annovar_convert + ' -format vcf4 ' + ifile + ' > ' + ofile
   run_tool(cmd)

def filter_by_dbtype(ifile, ofile, annovar, database_dir, dbtype, buildver, otherinfo):
   """Perform filter-based annotation. For each variants, filter it
   against a variation database, such as the 1000 Genomes Project
   database, to identify whether it has been reporte in the
   database. Exact match of nucleotide position and nucleotide
   composition are required."""
   
   annovar_annotate_variation = annovar + os.sep + 'annotate_variation.pl'
   #check whether annotate_variation.pl program exists or not
   if not is_exe(annovar_annotate_variation):
      print annovar_annotate_variation  + ' program does not exits'
      exit(1)

   cmd = 'perl ' + annovar_annotate_variation + ' -filter -dbtype ' + dbtype + ' -buildver ' + buildver + ' -comment -out ' + ofile + ' ' + ifile + ' ' +  database_dir + otherinfo
   run_tool(cmd)


if __name__ == "__main__":

   inputdir, outputdir, vcffilter, annovar, database_dir = parseArg(sys.argv[1:])
   print inputdir + ';' + outputdir + ';' + vcffilter + ';' + annovar + ';' + database_dir
   #check whether the input directory exists or not
   if not os.path.exists(inputdir):
      print inputdir  + ' directory does not exits'
      exit(1)
      
   '''
   check whether the output directory exists or not. Create one if
   it does not exist
   '''
   if not os.path.exists(outputdir):
      os.makedirs(outputdir)
      
   dirs = os.listdir(inputdir)
   if(len(dirs) == 0):
      print inputdir  + ' directory is empty'
      exit(1)
   
   logging.info(myapp + ' start running')
      
   #loop through the directory and process file one by one
   for file in dirs:
      if file.startswith('.'):
         continue
      prefix, ext = os.path.splitext(file)
      
      #quality filtering before variant annotation
      qualityFilter(inputdir + os.sep + file, outputdir + os.sep + prefix + '.qc.vcf', vcffilter)

      #convert vcf file to annovar input file format
      vcf2AnnovarInput(outputdir + os.sep + prefix + '.qc.vcf', outputdir + os.sep + prefix + '.annovar', annovar)

      #population frequency: annovar filter based annotation based on 1000g2015aug_eur db
      database_file = database_dir + os.sep + 'hg19_EUR.sites.2015_08.txt'
      if not os.path.isfile(database_file):
         print database_file + ' does not exist'
         exit(1)
      else:
         filter_by_dbtype(outputdir + os.sep + prefix + '.annovar', outputdir + os.sep + prefix, annovar, database_dir, '1000g2015aug_eur', 'hg19', ' -reverse -maf 0.05')
         
      
      database_file = database_dir + os.sep + 'hg19_exac03.txt'
      if not os.path.isfile(database_file):
         print database_file + ' does not exist'
         exit(1)
      else:
         #population frequency: annovar filter based annotation based on exac03 db
         filter_by_dbtype(outputdir + os.sep + prefix + '.annovar', outputdir + os.sep + prefix, annovar, database_dir, 'exac03', 'hg19', ' -reverse -otherinfo')
         
      
      database_file = database_dir + os.sep + 'hg19_snp138.txt'
      if not os.path.isfile(database_file):
         print database_file + ' does not exist'
         exit(1)
      else:
         #DBsnp annotated variants: annovar filter based annotation based on snp138 db
         filter_by_dbtype(outputdir + os.sep + prefix + '.annovar', outputdir + os.sep + prefix, annovar, database_dir, 'snp138', 'hg19','')
         
      database_file = database_dir + os.sep + 'hg19_ljb23_gerp++.txt'
      if not os.path.isfile(database_file):
         print database_file + ' does not exist'
         exit(1)
      else:
         #Use only coding variants to verify the variant effect: annovar filter based annotation based on ljb23_gerp++ db
         filter_by_dbtype(outputdir + os.sep + prefix + '.annovar', outputdir + os.sep + prefix, annovar, database_dir, 'ljb23_gerp++', 'hg19', ' -otherinfo')

      database_file = database_dir + os.sep + 'hg19_clinvar_20150629.txt'
      if not os.path.isfile(database_file):
         print database_file + ' does not exist'
         exit(1)
      else:
         #Clinvar Pathogenic:   annovar filter based annotation based on clinvar_20150629 db
         filter_by_dbtype(outputdir + os.sep + prefix + '.annovar', outputdir + os.sep + prefix, annovar, database_dir, 'clinvar_20150629', 'hg19', '')
