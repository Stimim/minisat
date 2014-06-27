#Polybori
  * Download from [here](http://sourceforge.net/projects/polybori/files/)
  * install boost library and scons
  `sudo apt-get install libboost-all-dev scons`
  * get m4ri
  `cd polybori-0.8.3/ ; wget http://m4ri.sagemath.org/downloads/m4ri/m4ri-20121224.tar.gz`
  * build
  `scons devel-install DEVEL-PREFIX=../local PREFIX=../local`

#Minisat
  * edit Makefile (see git log)
