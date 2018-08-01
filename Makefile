all:

check:
	R CMD build .
	R CMD check .

install.deps:
	R -q -e 'source("http://bioconductor.org/biocLite.R") ; biocLite("faahKO")'

clean:
	$(RM) -r ..Rcheck

.PHONY: all clean check install.deps
