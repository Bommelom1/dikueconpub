
FUTHARK=futhark

all: dpsolve_demo.html

dpsolve_demo.md: dpsolve_demo.fut
	$(FUTHARK) literate dpsolve_demo.fut

dpsolve_demo.html: dpsolve_demo.md
	pandoc $< --css pandoc.css --citeproc --bibliography root.bib -s -t html5 '--mathjax=https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML' -H header --syntax-definition futhark.xml -o $@

dpsolve_demo.tex: dpsolve_demo.md
	pandoc $< --citeproc --bibliography root.bib -s -t latex --template=latex2.template --syntax-definition futhark.xml -o $@

dpsolve_demo.pdf: dpsolve_demo.tex
	pdflatex $<
	pdflatex $<
	pdflatex $<

.PHONY: clean
clean:
	rm -rf *~ dpsolve_demo-img dpsolve_demo.md dpsolve_demo.html dpsolve_demo doc *.c dpsolve_demo.pdf dpsolve_demo.tex *.toc *.aux *.log *.out

.PHONY: futversion
futversion:
	$(FUTHARK) --version

doc:
	$(FUTHARK) doc -v -o doc .
