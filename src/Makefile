#################
### PLATO MAKEFILE
#################

all:
	make MDR

MDR:
	make libs
	cd sMDR; make
	cd pMDR; make

libs:
	cd mdrlib; make

clean: 
	cd lib; rm -f *.a
	cd mdrlib; make clean
	cd sMDR; make clean
	cd pMDR; make clean

tarball:
	make clean
	tar cvfz mdr_source_`date '+%Y%m%d'`.tar.gz \
	mdrlib sMDR pMDR
