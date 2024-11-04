# PATH for libraries is guessed
CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

ifndef INSTALL_PREFIX
    INSTALL_PREFIX = /usr/local
endif

HET_KMERS_INST = $(INSTALL_PREFIX)/bin/smudgeplot.py $(INSTALL_PREFIX)/bin/hetmers
SMUDGEPLOT_INST = $(INSTALL_PREFIX)/bin/smudgeplot_plot.R $(INSTALL_PREFIX)/bin/centrality_plot.R

.PHONY : install
install : $(HET_KMERS_INST) $(SMUDGEPLOT_INST) $(CUTOFF_INST) 

$(INSTALL_PREFIX)/bin/% : exec/%
	install -C $< $(INSTALL_PREFIX)/bin

exec/hetmers: src_ploidyplot/PloidyPlot.c src_ploidyplot/libfastk.c src_ploidyplot/libfastk.h src_ploidyplot/matrix.c src_ploidyplot/matrix.h
	gcc $(CFLAGS) -o $@ src_ploidyplot/PloidyPlot.c src_ploidyplot/libfastk.c src_ploidyplot/matrix.c -lpthread -lm


.PHONY : clean
clean :
	rm -f exec/hetmers
