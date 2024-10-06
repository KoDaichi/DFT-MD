include ./make.inc
all:
	cd tspace; make all
	cd src; make all
clean:
	cd tspace; make clean
	cd src; make clean
erase:
	cd tspace; make erase
	cd src; make erase
	make erase
package:
	git ls-files -z | xargs -0 tar cf band_symm.tar
