# kyleaoman_utilities
Useful python functions for various tasks, especially as dependencies for other repositories/modules by kyleaoman.

**Installation:**
 - Download via web UI, or 'git clone https://github.com/kyleaoman/kyleaoman_utilities.git'
 - Install dependencies if necessary (see 'setup.py'), some may be found in other repositories by kyleaoman.
 - Global install (Linux): 
   - cd to directory with 'setup.py'
   - run 'sudo pip install -e .' (-e installs via symlink, so pulling repository will do a 'live' update of the installation)
 - User install (Linux):
   - cd to directory with 'setup.py'
   - ensure '~/lib/python3.6/site-packages' or similar is on your PYTHONPATH (e.g. 'echo $PYTHONPATH'), if not, add it (perhaps in .bash_profile or similar)
   - run 'pip install --prefix ~ -e .' (-e installs via symlink, so pulling repository will do a 'live' update of the installation)
 - cd to a directory outside the module and launch python; you should be able to do 'from kyleaoman_utilities import *'
