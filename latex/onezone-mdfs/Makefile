
TEXCOMPILER := pdflatex
BIBCOMPILER := bibtex

all: onezone-mdfs.pdf

onezone-mdfs.pdf: onezone-mdfs.tex onezone-mdfs.bib
	@ $(TEXCOMPILER) $(basename $@)
	@ $(BIBCOMPILER) $(basename $@)
	@ $(TEXCOMPILER) $(basename $@)
	@ $(TEXCOMPILER) $(basename $@)

.PHONY: clean
clean:
	@ rm -f *.aux
	@ rm -f *.bbl
	@ rm -f *.blg
	@ rm -f *.log
	@ rm -f *.out
	@ rm -f onezone-mdfs.pdf
