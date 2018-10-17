PKG_TAR=Risa_$(shell grep Version DESCRIPTION | sed 's/^Version: //').tar.gz

all:

build: $(PKG_TAR)

$(PKG_TAR): R/* man/* DESCRIPTION NAMESPACE NEWS README.md vignettes/*
	R CMD build .

check: $(PKG_TAR)
	R CMD check $(shell ls -t Risa*.tar.gz | head -n 1)

install.deps:
	R -q -e 'install.packages(c("Rcpp"))'
	R -q -e 'source("http://bioconductor.org/biocLite.R") ; biocLite(c("biocViews", "affy", "xcms", "faahKO"))'

test:
	R -q -e "devtools::test('$(CURDIR)', reporter = c('progress', 'fail'))"

clean:
	$(RM) -r ..Rcheck Risa_*.tar.gz

.PHONY: all clean check install.deps build
