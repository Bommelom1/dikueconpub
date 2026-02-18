
FUTHARK=futhark

#FUTHARK=$(HOME)/.local/bin/futhark

.PHONY: all
all:
	$(FUTHARK) c trmodel.fut

.PHONY: clean
clean:
	rm -rf *~ doc *.c trmodel trmodel_test trmodel_static trmodel_test.fut.*

.PHONY: futversion
futversion:
	$(FUTHARK) --version

.PHONY: test
test:
	$(FUTHARK) test trmodel_test.fut

doc:
	$(FUTHARK) doc -v -o doc .
