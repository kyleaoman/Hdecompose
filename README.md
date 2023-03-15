# Hdecompose
Decomposition of Hydrogen content of simulation particles into neutral/ionized and atomic/molecular. Implementations of Blitz & Rosolowsky (2006) and Rahmati et al (2013).

**Installation:**

From pypi:
 - 'pip install hdecompose'

From github:
 - Download via web UI, or 'git clone https://github.com/kyleaoman/Hdecompose.git'
 - Install dependencies if necessary (see 'setup.py'), some may be found in other repositories by kyleaoman.
 - Global install (Linux): 
   - cd to directory with 'setup.py'
   - run 'sudo pip install -e .' (-e installs via symlink, so pulling repository will do a 'live' update of the installation)
 - User install (Linux):
   - cd to directory with 'setup.py'
   - ensure '~/lib/python3.7/site-packages' or similar is on your PYTHONPATH (e.g. 'echo $PYTHONPATH'), if not, add it (perhaps in .bash_profile or similar)
   - run 'pip install --prefix ~ -e .' (-e installs via symlink, so pulling repository will do a 'live' update of the installation)
 - cd to a directory outside the module and launch python; you should be able to do 'from Hdecompose import *'
