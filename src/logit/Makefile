
FUTBACK=c
INPUTDATA=123 10000 3 3 2

.PHONY: all
all: logit.out
	cat $<

logit.exe: logit.fut lib
	futhark $(FUTBACK) -o $@ $<

logit.in: logit.exe
	echo $(INPUTDATA) | ./$< -e generate_data -b > $@

logit.out: logit.exe logit.in
	./$< -e criterion < logit.in -t $@ -r 10

lib:
	futhark pkg sync

.PHONY: clean
clean:
	rm -rf *~ *.c *.in *.exe *.out
