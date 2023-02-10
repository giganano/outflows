
TEXCOMPILER := pdflatex
BIBCOMPILER := bibtex

all: yields-onezone.pdf

%.pdf: %.tex %.bib
	@ $(TEXCOMPILER) $(basename $@)
	@ $(BIBCOMPILER) $(basename $@)
	@ $(TEXCOMPILER) $(basename $@)
	@ $(TEXCOMPILER) $(basename $@)

.PHONY: clean
clean:
	@ rm -f yields-onezone.aux
	@ rm -f yields-onezone.bbl
	@ rm -f yields-onezone.blg
	@ rm -f yields-onezone.log
	@ rm -f yields-onezone.out
	@ rm -f yields-onezone.pdf

