# PATH for libraries is guessed
RPACKAGE = smudgeplot
$(eval RPATH := $(shell Rscript -e "noquote(.libPaths())" | tail -1 | cut -f 2 -d ' '))
CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

ifndef INSTALL_PREFIX
    INSTALL_PREFIX = /usr/local
endif

R_INSTALLATION = $(RPATH)/$(RPACKAGE)
HET_KMERS_INST = $(INSTALL_PREFIX)/bin/smudgeplot.py $(INSTALL_PREFIX)/bin/PloidyPlot
SMUDGEPLOT_INST = $(INSTALL_PREFIX)/bin/smudgeplot_plot.R

.PHONY : install
install : $(HET_KMERS_INST) $(SMUDGEPLOT_INST) $(CUTOFF_INST) $(R_INSTALLATION) 
# - removing the R package from the list

$(INSTALL_PREFIX)/bin/% : exec/%
	install -C $< $(INSTALL_PREFIX)/bin

$(R_INSTALLATION) : R/*
	Rscript install.R

exec/PloidyPlot: src_ploidyplot/PloidyPlot.c src_ploidyplot/libfastk.c src_ploidyplot/libfastk.h src_ploidyplot/matrix.c src_ploidyplot/matrix.h
	gcc $(CFLAGS) -o $@ src_ploidyplot/PloidyPlot.c src_ploidyplot/libfastk.c src_ploidyplot/matrix.c -lpthread -lm


.PHONY : clean
clean :
	rm -f exec/PloidyPlot
