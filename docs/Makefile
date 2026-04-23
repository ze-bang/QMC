# Build the pedagogical LaTeX documents (SSE and DQMC).
#
#   make            -> both PDFs
#   make sse_qmc.pdf
#   make dqmc.pdf
#   make clean      -> remove auxiliary files
#   make distclean  -> also remove the PDFs

LATEX     ?= pdflatex
LATEXOPTS ?= -interaction=nonstopmode -halt-on-error

DOCS     = sse_qmc dqmc
PDFS     = $(addsuffix .pdf,$(DOCS))
AUX_EXT  = aux log out toc bbl blg fls fdb_latexmk synctex.gz

.PHONY: all clean distclean view

all: $(PDFS)

%.pdf: %.tex
	@$(LATEX) $(LATEXOPTS) $< >/dev/null
	@$(LATEX) $(LATEXOPTS) $< >/dev/null
	@printf "  PDF  %s (%s bytes, %s pages)\n" $@ "$$(wc -c < $@)" "$$(pdfinfo $@ 2>/dev/null | awk '/Pages/ {print $$2}')"

clean:
	@for d in $(DOCS); do rm -f $$(printf "$$d.%s " $(AUX_EXT)); done

distclean: clean
	@rm -f $(PDFS)

view: $(PDFS)
	@for p in $(PDFS); do xdg-open $$p >/dev/null 2>&1 & done
