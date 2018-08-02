all:

check:
	R CMD build .
	R CMD check .

install.deps:
	R -q -e 'install.packages(c("Rcpp"))'
	R -q -e 'source("http://bioconductor.org/biocLite.R") ; biocLite(c("biocViews", "affy", "xcms", "faahKO"))'

clean:
	$(RM) -r ..Rcheck

.PHONY: all clean check install.deps
