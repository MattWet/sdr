


# runs checks with --run-donttest; can take a while
devcheck:
	Rscript -e "devtools::check()"

coverage: install
	Rscript -e 'covr::report(covr::package_coverage(), file = "stagewise_coverage.html")'

# R build and install
packageversion:=$(shell cat DESCRIPTION | egrep Version | sed 's/Version://g')

build:
	cd ../ && R CMD build stagewise

install: build
	cd ../ && R CMD INSTALL stagewise_$(shell printf "%s"${packageversion}).tar.gz

check:
	cd ../ && R CMD check --as-cran stagewise_$(shell printf "%s"${packageversion}).tar.gz

# Running tinytests
.PHONY: test
test: install
	Rscript -e "library('stagewise'); tinytest::test_all()"

.PHONY: all
all:
	make install
	make test
	make check
