
all:

logit.exe: logit.fut
	futhark opencl -o $@ $<

logit.in: logit.exe

clean:
	rm -rf *~ *.c *.in *.exe
