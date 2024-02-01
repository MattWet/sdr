


# runs checks with --run-donttest; can take a while
devcheck:
	Rscript -e "devtools::check()"

coverage: install
	Rscript -e 'covr::report(covr::package_coverage(), file = "sdr_coverage.html")'

# R build and install
packageversion:=$(shell cat DESCRIPTION | egrep Version | sed 's/Version://g')

build:
	cd ../ && R CMD build sdr

install: build
	cd ../ && R CMD INSTALL sdr_$(shell printf "%s"${packageversion}).tar.gz

check:
	cd ../ && R CMD check --as-cran sdr_$(shell printf "%s"${packageversion}).tar.gz

# Running tinytests
.PHONY: test
test: install
	Rscript -e "library('sdr'); tinytest::test_all()"

.PHONY: all
all:
	make install
	make test
	make check
