all:

check:
	R CMD build .
	R CMD check .

install.deps:
	R -q -e 'install.packages(c("Rcpp"))'
	R -q -e 'source("http://bioconductor.org/biocLite.R") ; biocLite(c("biocViews", "affy", "xcms", "faahKO"))'

test:
	R -q -e "devtools::test('$(CURDIR)', reporter = c('progress', 'fail'))"

clean:
	$(RM) -r ..Rcheck

.PHONY: all clean check install.deps
