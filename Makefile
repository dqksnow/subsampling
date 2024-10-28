objects := $(wildcard R/*.R) DESCRIPTION
version := $(shell grep -E "^Version:" DESCRIPTION | awk '{print $$NF}')
pkg := $(shell  grep -E "^Package:" DESCRIPTION | awk '{print $$NF}')
tar := $(pkg)_$(version).tar.gz
checkLog := $(pkg).Rcheck/00check.log

.PHONY: check
check: $(tar)
	R CMD check --as-cran --no-manual $(tar)

.PHONY: build
build: $(tar)

.PHONY: install
install: $(tar)
	R CMD INSTALL $(tar)

# Target to build the package tarball
$(tar): $(objects)
	@$(RM) -rf src/RcppExports.cpp R/RcppExports.R
	@Rscript -e "library(methods);" \
	-e "Rcpp::compileAttributes()" \
	-e "devtools::document();";
	R CMD build .

# Target for generating vignettes with HTML output
.PHONY: preview
preview: $(vignettes)
	@Rscript -e "library(methods); pkgdown::build_site();"

# Check log generation
$(checkLog): $(tar)
	R CMD check --as-cran $(tar)


# Build each vignette to HTML format
vignettes/%.html: vignettes/%.Rmd
	Rscript -e "library(methods); rmarkdown::render('$?')"

# Render README file from RMarkdown
.PHONY: readme
readme: README.md
README.md: README.Rmd
	@Rscript -e "rmarkdown::render('$<')"

# Generate TAGS for easy code navigation
.PHONY: TAGS
TAGS:
	Rscript -e "utils::rtags(path = 'R', ofile = 'TAGS')"


.PHONY: clean
clean:
	@$(RM) -rf *~ */*~ *.Rhistroy *.tar.gz src/*.so src/*.o *.Rcheck/ .\#*