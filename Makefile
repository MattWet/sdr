

#.PHONY: document
#document:
#	Rscript -e "devtools::document()"

#.PHONY: docs
#docs:
#	-rm -rf docs
#	Rscript -e "pkgdown::build_site()"

.PHONY: install
install:
	Rscript -e "devtools::install()"

test: 
	Rscript -e "devtools::load_all(); tinytest::test_all()"

check:
	Rscript -e "devtools::check()"

coverage: install
	Rscript -e 'covr::report(covr::package_coverage(), file = "sdr_coverage.html")'

.PHONY: all
all:
	make install
	make check
	make docs
