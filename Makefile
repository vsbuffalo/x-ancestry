NSIMS=5000
LNSIMS=300
DYLIB=_ancestrysim2.so
SEED=42

figs: images/x-graph-01.tex paper/images/x-tree.eps paper/images/x-arc.eps paper/images/prob-anc.eps paper/images/x-rm-tree.eps all-plots

clean: clean-x clean-auto clean-fig1 clean-xcuz

data: data-x data-auto data-fig1 data-xcuz inf

inf: data/x-halfcousin-blockcounts-rms.txt

#### Autosome Simulations
AUTODATA=data/auto-shared-ancestor-blockcounts-py.txt data/auto-shared-ancestor-blockcounts.txt data/auto-shared-ancestor-blocklens.txt
data-auto: $(AUTODATA)

clean-auto:
	rm -rf $(AUTODATA)

data/auto-shared-ancestor-blockcounts-py.txt: ancestrysim/ancestralsim.py
	# 5k sims too much for old routine
	python ancestrysim/ancestralsim.py auto --seed $(SEED) --ngens 1:11 --nsims 1000 > $@

data/auto-shared-ancestor-blockcounts.txt: ancestrysim/ancestrysim2.py $(DLIB)
	python ancestrysim/ancestrysim2.py ancestor --auto --ngens 1:10 --nsims $(NSIMS) --seed $(SEED) > $@

data/auto-shared-ancestor-blocklens.txt: ancestrysim/ancestrysim2.py $(DLIB)
	python ancestrysim/ancestrysim2.py ancestor --auto --ngens 1:10 --nsims $(LNSIMS) --seed $(SEED) --len > $@

#### X Chromosome Simulations
XDATA=data/x-shared-ancestor-blockcounts.txt data/x-shared-ancestor-blockcounts-py.txt data/x-shared-ancestor-blocklens.txt
data-x: $(XDATA)

clean-x:
	rm -rf $(XDATA)

data/x-shared-ancestor-blockcounts-py.txt: ancestrysim/ancestralsim.py
	# 5k sims too much for old routine
	python ancestrysim/ancestralsim.py x --ngens 1:11 --nsims 1000 --seed $(SEED) > $@

data/x-shared-ancestor-blockcounts.txt: ancestrysim/ancestrysim2.py $(DLIB)
	python ancestrysim/ancestrysim2.py ancestor --x --ngens 1:10 --nsims $(NSIMS) --seed $(SEED) > $@

data/x-shared-ancestor-blocklens.txt: ancestrysim/ancestrysim2.py
	python ancestrysim/ancestrysim2.py ancestor --x --ngens 1:10 --nsims $(NSIMS) --len --seed $(SEED) > $@


#### Figure 1 Simulations
FIG1DATA=data/genetic-auto-all-ancs-blockcounts.txt data/genetic-x-all-ancs-blockcounts.txt data/genetic-x-x-ancs-blockcounts.txt
data-fig1: $(FIG1DATA)

clean-fig1:
	rm -rf $(FIG1DATA)

data/genetic-auto-all-ancs-blockcounts.txt: ancestrysim/ancestrysim2.py
	python $< ancestor --auto --ngens 1:13 --nsims 100 --seed $(SEED) --gen-anc --all-ancestors --seed $(SEED) > $@

data/genetic-x-all-ancs-blockcounts.txt: ancestrysim/ancestrysim2.py
	python $< ancestor --x --ngens 1:13 --nsims 100 --seed $(SEED) --gen-anc --all-ancestors --seed $(SEED) > $@

data/genetic-x-x-ancs-blockcounts.txt: ancestrysim/ancestrysim2.py
	python $< ancestor --x --ngens 1:13 --nsims 100 --seed $(SEED) --gen-anc --seed $(SEED) > $@


#### X Chromosome Cousin Simulations
XCUZDATA=data/x-halfcousin-blockcounts.txt data/x-fullcousin-blockcounts.txt data/x-halfcousin-lens.txt
data-xcuz: $(XCUZDATA)

clean-xcuz:
	rm -rf $(XCUZDATA)

data/x-halfcousin-blockcounts.txt: ancestrysim/ancestrysim2.py $(DLIB)
	python ancestrysim/ancestrysim2.py cousin --x --ngens 1:10 --shared-random --nsims $(NSIMS) --seed $(SEED) > $@

data/x-fullcousin-blockcounts.txt: ancestrysim/ancestrysim2.py $(DLIB)
	python ancestrysim/ancestrysim2.py cousin --x --ngens 1:10 --shared-random --nsims $(NSIMS) --full-sibs --seed $(SEED) > $@

data/x-halfcousin-lens.txt: ancestrysim/ancestrysim2.py $(DLIB)
	python ancestrysim/ancestrysim2.py cousin --x --ngens 1:10 --nsims $(NSIMS) --half-sibs --len --seed $(SEED) > $@


## Misc Figures
figs: paper/images/x-tree.eps paper/images/x-arc.png paper/images/x-arc.eps

images/x-graph-01.tex: images/x-graph-01.dot
	dot2tex --figonly -tmath $< > $@

paper/images/x-rm-tree.eps:
	 python  ancestrysim/ancestrysim2.py dot --ngens 5 | dot -Tsvg -o paper/images/x-rm-tree.svg
	rsvg-convert -f eps -b white paper/images/x-rm-tree.svg > $@

paper/images/x-tree.eps: jsviz/x-tree.svg
	rsvg-convert -f eps -b white $^ > $@

paper/images/prob-anc.eps: jsviz/prob-anc.svg
	rsvg-convert -f eps -b white $^ > $@

paper/images/x-arc.png: jsviz/x-arc.svg
	# many PDF readers don't diplay the vector version well, so PNG is an alternate option.
	rm -f $@
	svg2png $^ --output $@

paper/images/x-arc.eps: jsviz/x-arc.svg
	rsvg-convert -f eps -b white $^ > $@

#### Paper Figures
all-plots: plots.R dist-functions.R data
	Rscript --vanilla --no-save --no-restore --no-init-file plots.R

## Inference
clean-inf:
	rm -rf data/x-halfcousin-blockcounts-rms.txt

data/x-halfcousin-blockcounts-rms.txt: ancestrysim/ancestrysim2.py $(DLIB)
		python ancestrysim/ancestrysim2.py cousin --x --ngens 1:10 --shared-random --nsims 30000 --rms --half-sibs --seed $(SEED) > $@

.PHONY: all clean data-fig1 data-x data-auto clean-fig1 clean-x clean-auto figs inf clean-inf
