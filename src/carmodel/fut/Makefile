
FUTHARK=futhark

.PHONY: all
all: test
	$(MAKE) -C util dpsolve_demo.pdf

run_bellman_ad:
	$(FUTHARK) c run_bellman_ad.fut

run_bellman_j:
	$(FUTHARK) c run_bellman_j.fut

run_bellman:
	$(FUTHARK) c run_bellman.fut

.PHONY: clean
clean:
	rm -rf *~ *.c run_bellman_ad run_bellman_j run_bellman run_bellman.fut.* run_bellman_ad.fut.* run_bellman_j.fut.* *.exe *.tm *.tmavg *.pdf data1.dat data2.dat *.out data1j.dat
	$(MAKE) -C util clean
	$(MAKE) -C autotrade clean

.PHONY: test
test:
	$(FUTHARK) test -c run_bellman.fut run_bellman_ad.fut run_bellman_j.fut

.PHONY: futversion
futversion:
	$(FUTHARK) --version

%-c.exe: %.fut
	$(FUTHARK) c -o $@ $<

%-m.exe: %.fut
	$(FUTHARK) multicore -o $@ $<

%-o.exe: %.fut
	$(FUTHARK) opencl -o $@ $<

run_bellman-%-C.tm: run_bellman-c.exe
	echo "$*" | sed 's/-/ /g' | ./$< -e mainn -r 5 -t $@ > run_bellman-$*-c.out

run_bellman-%-M.tm: run_bellman-m.exe
	echo "$*" | sed 's/-/ /g' | ./$< -e mainn -r 5 -t $@ > run_bellman-$*-m.out

run_bellman-%-O.tm: run_bellman-o.exe
	echo "$*" | sed 's/-/ /g' | ./$< -d AMD -e mainn -r 5 -t $@ > run_bellman-$*-o.out

run_bellman_j-%-C.tm: run_bellman_j-c.exe
	echo "$*" | sed 's/-/ /g' | ./$< -e mainn -r 5 -t $@ > run_bellman_j-$*-c.out

run_bellman_j-%-M.tm: run_bellman_j-m.exe
	echo "$*" | sed 's/-/ /g' | ./$< -e mainn -r 5 -t $@ > run_bellman_j-$*-m.out

run_bellman_j-%-O.tm: run_bellman_j-o.exe
	echo "$*" | sed 's/-/ /g' | ./$< -d AMD -e mainn -r 5 -t $@ > run_bellman_j-$*-o.out

run_bellman_ad-%-C.tm: run_bellman_ad-c.exe
	echo "$*" | sed 's/-/ /g' | ./$< -e mainn -r 5 -t $@ > run_bellman_ad-$*-c.out

run_bellman_ad-%-M.tm: run_bellman_ad-m.exe
	echo "$*" | sed 's/-/ /g' | ./$< -e mainn -r 5 -t $@ > run_bellman_ad-$*-m.out

run_bellman_ad-%-O.tm: run_bellman_ad-o.exe
	echo "$*" | sed 's/-/ /g' | ./$< -d AMD -e mainn -r 5 -t $@ > run_bellman_ad-$*-o.out

%.tmavg: %.tm
	cat $< | python3 -c"import sys; nums = [int(n) for n in sys.stdin]; print(sum(nums)/len(nums)/1000000)" > $@

%.png: %.pdf
	convert -density 300 $< $@

plots.pdf: plot1.pdf plot1j.pdf plot1ad.pdf plot2.pdf
	pdfunite $^ $@

# N-C-Ax
DATA1=2-10-20 4-10-20 6-10-20 8-10-20 10-10-20 \
     4-5-20 4-15-20 4-20-20 4-25-20 4-30-20 4-35-20

DATA1_C_FILES=$(DATA1:%=run_bellman-%-C.tmavg)
DATA1_M_FILES=$(DATA1:%=run_bellman-%-M.tmavg)
DATA1_O_FILES=$(DATA1:%=run_bellman-%-O.tmavg)
#DATA1_FILES=$(DATA1_O_FILES) $(DATA1_M_FILES) $(DATA1_C_FILES)
DATA1_FILES=$(DATA1_M_FILES) $(DATA1_C_FILES)

data1.dat: $(DATA1_FILES)
	for n in $^ ; do \
             /bin/echo -n $$n | sed 's/[a-z_]*//g' | sed 's/[-.]/ /g' >> $@ ; \
             cat $$n >> $@ ; \
        done

plot1.pdf: plot1-noo.gp data1.dat matlab.dat
	gnuplot -e "load \"$<\""

DATA1_J=4-5-25 4-15-25 4-20-25 4-25-25 4-30-25 4-35-25

DATA1_J_C_FILES=$(DATA1_J:%=run_bellman_j-%-C.tmavg)
DATA1_J_M_FILES=$(DATA1_J:%=run_bellman_j-%-M.tmavg)
DATA1_J_O_FILES=$(DATA1_J:%=run_bellman_j-%-O.tmavg)
#DATA1_J_FILES=$(DATA1_J_O_FILES) $(DATA1_J_M_FILES) $(DATA1_J_C_FILES)
DATA1_J_FILES=$(DATA1_J_M_FILES) $(DATA1_J_C_FILES)

data1j.dat: $(DATA1_J_FILES)
	for n in $^ ; do \
             /bin/echo -n $$n | sed 's/[a-z_]*//g' | sed 's/[-.]/ /g' >> $@ ; \
             cat $$n >> $@ ; \
        done

plot1j.pdf: plot1j-noo.gp data1j.dat matlabj.dat
	gnuplot -e "load \"$<\""

DATA1_AD=4-5-25 4-15-25 4-20-25

DATA1_AD_C_FILES=$(DATA1_AD:%=run_bellman_ad-%-C.tmavg)
DATA1_AD_M_FILES=$(DATA1_AD:%=run_bellman_ad-%-M.tmavg)
DATA1_AD_O_FILES=$(DATA1_AD:%=run_bellman_ad-%-O.tmavg)
#DATA1_AD_FILES=$(DATA1_AD_O_FILES) $(DATA1_AD_M_FILES) $(DATA1_AD_C_FILES)
DATA1_AD_FILES=$(DATA1_AD_M_FILES) $(DATA1_AD_C_FILES)

data1ad.dat: $(DATA1_AD_FILES)
	for n in $^ ; do \
             /bin/echo -n $$n | sed 's/[a-z_]*//g' | sed 's/[-.]/ /g' >> $@ ; \
             cat $$n >> $@ ; \
        done

plot1ad.pdf: plot1ad-noo.gp data1ad.dat matlabj.dat
	gnuplot -e "load \"$<\""

DATA2=2-35-25 4-35-25 6-35-25 8-35-25 10-35-25
DATA2_C_FILES=$(DATA2:%=run_bellman-%-C.tmavg)
DATA2_M_FILES=$(DATA2:%=run_bellman-%-M.tmavg)
DATA2_O_FILES=$(DATA2:%=run_bellman-%-O.tmavg)
#DATA2_FILES=$(DATA2_O_FILES) $(DATA2_M_FILES) $(DATA2_C_FILES)
DATA2_FILES=$(DATA2_M_FILES) $(DATA2_C_FILES)
data2.dat: $(DATA2_FILES)
	for n in $^ ; do \
             /bin/echo -n $$n | sed 's/[a-z_]*//g' | sed 's/[-.]/ /g' >> $@ ; \
             cat $$n >> $@ ; \
        done

plot2.pdf: plot2-noo.gp data2.dat matlab2.dat
	gnuplot -e "load \"$<\""


.PHONY: ad-back-c
ad-back-c:
	$(FUTHARK) c run_bellman_ad.fut
	echo 10 25 200 | ./run_bellman_ad -r 5 -t /dev/stderr

.PHONY: back-opencl
back-opencl:
	$(FUTHARK) opencl run_bellman.fut
	echo 10 25 200 | ./run_bellman -d AMD -r 5 -t /dev/stderr

.PHONY: ad-back-opencl
ad-back-opencl:
	$(FUTHARK) opencl run_bellman_ad.fut
	echo 10 25 200 | ./run_bellman_ad -d AMD -r 5 -t /dev/stderr
