OUTDIR=out
INDIR=.

all:
	mkdir -p $(OUTDIR)
	pdflatex -output-directory=$(OUTDIR) $(INDIR)/main.tex
	biber --output-directory=$(OUTDIR) $(INDIR)/main
	pdflatex -output-directory=$(OUTDIR) $(INDIR)/main.tex
	ln -sf $(OUTDIR)/main.pdf ./main.pdf

clean:
	rm -rf $(OUTDIR) main.pdf

.PHONY: all clean
