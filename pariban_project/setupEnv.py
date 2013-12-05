#!/bin/env python3

import subprocess, os, sys, urllib.request, tarfile, shlex, shutil
g_packageDir=os.path.abspath(os.path.dirname(__file__)) # get the absolute path of the script directory
sys.path.append(g_packageDir) # append to system path
g_binDir=os.path.join(g_packageDir, 'bin') # get the path of the bin directory
sys.path.append(g_binDir) # add to system path
os.environ["PATH"] += os.pathsep + g_packageDir # add to OS environment variable PATH
os.environ["PATH"] += os.pathsep + g_binDir # add to OS environment variable PATH

from utils import Assert # Assert class prints backtrace

class PackageRequirements:
  def __init__ (self, x_debug=False):
    self.debug= x_debug
    self.requiredPackages = {} # stores download path of each package in a map
    self.requiredPackages['scan_for_matches'] = "http://www.theseed.org/servers/downloads/scan_for_matches.tgz"
    self.requiredPackages['RNAfold'] = "http://www.tbi.univie.ac.at/~ronny/RNA/packages/source/ViennaRNA-2.1.5.tar.gz"
    self.requiredPackages['einverted'] = "ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz"
    self.checkPackages() # check if package is installed

  def checkPackages (self):
    for package in self.requiredPackages: # for all packages
      location = self.which(package) # check if the package is found in system path or bin directory
      if location: # found the package
        if self.debug:
          print("found", package, ":", location)
      else: # ask the user if she wants to install it
        choice = input("required package \""+package+"\" was not found. automatically install? (yes/no)?: ")
        if choice == "yes":
          self.installPackage(package) # install the package
        else: # the package was not found, so we cannot run successfully
          Assert(False, package+"  was not be found, please fix path and rerun or install manually if it is not installed")

        
  def installPackage (self, x_package):
    try:
      if x_package == "scan_for_matches" or x_package == "RNAfold" or x_package == "einverted": # we support these currently
        destFile=os.path.split(self.requiredPackages[x_package])[1] # name of the downloaded file
        if self.debug:
          print("downloading", self.requiredPackages[x_package])
        urllib.request.urlretrieve(self.requiredPackages[x_package], destFile) # download the file
        if self.debug:
          print("extracting", destFile)
        destDir = self.extract(destFile) # extract the file and get directory name
        if self.debug:
          print("extracted in", destDir, " ... now installing ...")
        global g_packageDir
        cmd = "cd %s && ./configure --prefix=%s && make install" % (destDir, g_packageDir) # configure, compile and install
        res=subprocess.check_output(["bash", "-c", cmd], stderr=subprocess.STDOUT) # run the command
        if self.debug:
          print(res.decode()) # output of the run was in unicode, so convert to string
        if self.debug:
          print("removing dirs and files", destFile, destDir)
        shutil.rmtree(destDir, True) # remove directory tree of the extracted package
        os.remove(destFile) # remove downloaded package
    except: # some error, the package was not installed
      Assert(False, x_package+"  could not be installed")
        
  def extract (self, tar_url, extract_path='.'):
    tar = tarfile.open(tar_url, 'r') # open tar file
    tar.extractall() # extract all files and directories
    dirName = tar.getnames()[0] # top level directory name
#    for item in tar:
#      tar.extract(item, extract_path)
#      if item.name.find(".tgz") != -1 or item.name.find(".tar") != -1:
#        extract(item.name, "./" + item.name[:item.name.rfind('/')])
    tar.close() # close tar file
    return dirName
  
  def which (self, x_program):
    def is_exe(fpath):# checks if a file exists and executable
      return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(x_program) # split into directory and file pair
    if fpath: # provided a path, so check if it exists
      if is_exe(x_program):
        return program
    else:
      for path in os.environ["PATH"].split(os.pathsep): # for all locations in PATH variable (includes bin directory)
        path = path.strip('"') # get rid of quotes
        exe_file = os.path.join(path, x_program) # create the name of the file as /path/to/dir/program
        if is_exe(exe_file): return exe_file # found the executable
    return None

if __name__ == "__main__": # stand alone tester
  PackageRequirements(True)
