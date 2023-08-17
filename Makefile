# PATH for libraries is guessed
RPACKAGE = smudgeplot
$(eval PATH := $(shell Rscript -e "noquote(.libPaths())" | tail -1 | cut -f 2 -d ' '))
CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

ifndef INSTAL_PREFIX
    INSTAL_PREFIX = /usr/local
endif

R_INSTALLATION = $(PATH)/$(RPACKAGE)
HET_KMERS_INST = $(INSTAL_PREFIX)/bin/smudgeplot.py $(INSTAL_PREFIX)/bin/PloidyPlot
SMUDGEPLOT_INST = $(INSTAL_PREFIX)/bin/smudgeplot_plot.R

.PHONY : install
install : $(R_INSTALLATION) $(HET_KMERS_INST) $(SMUDGEPLOT_INST) $(CUTOFF_INST)

$(INSTAL_PREFIX)/bin/% : exec/%
	install -C $< $(INSTAL_PREFIX)/bin

$(R_INSTALLATION) : R/*
	Rscript install.R

exec/PloidyPlot: src/PloidyPlot.c src/libfastk.c src/libfastk.h src/matrix.c src/matrix.h
	gcc $(CFLAGS) -o $@ src/PloidyPlot.c src/libfastk.c src/matrix.c -lpthread -lm


.PHONY : clean
clean :
	rm exec/PloidyPlot
	rm -r dist/
	rm -r smudgeplot_KamilSJaron.egg-info/
