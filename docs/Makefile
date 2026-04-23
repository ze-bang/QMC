# Build the pedagogical SSE QMC LaTeX document.
#
#   make           -> sse_qmc.pdf
#   make clean     -> remove auxiliary files
#   make distclean -> also remove the PDF

LATEX     ?= pdflatex
LATEXOPTS ?= -interaction=nonstopmode -halt-on-error

DOC      = sse_qmc
SRC      = $(DOC).tex
PDF      = $(DOC).pdf
AUX_EXT  = aux log out toc bbl blg fls fdb_latexmk synctex.gz

.PHONY: all clean distclean view

all: $(PDF)

$(PDF): $(SRC)
	@$(LATEX) $(LATEXOPTS) $(SRC) >/dev/null
	@$(LATEX) $(LATEXOPTS) $(SRC) >/dev/null
	@echo "  PDF  $(PDF) ($$(wc -c < $(PDF)) bytes, $$(pdfinfo $(PDF) 2>/dev/null | awk '/Pages/ {print $$2 " pages"}'))"

clean:
	@rm -f $(addprefix $(DOC)., $(AUX_EXT))

distclean: clean
	@rm -f $(PDF)

view: $(PDF)
	@xdg-open $(PDF) >/dev/null 2>&1 &
