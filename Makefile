# PATH for libraries is guessed
CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

ifndef INSTALL_PREFIX
    INSTALL_PREFIX = /usr/local
endif

HET_KMERS_INST = $(INSTALL_PREFIX)/bin/smudgeplot.py $(INSTALL_PREFIX)/bin/hetmers $(INSTALL_PREFIX)/bin/extract_kmer_pairs

.PHONY : install
install : $(HET_KMERS_INST) $(SMUDGEPLOT_INST)

$(INSTALL_PREFIX)/bin/% : exec/%
	install -C $< $(INSTALL_PREFIX)/bin

$(INSTALL_PREFIX)/bin/smudgeplot.py: smudgeplot.py
	install -C $< $(INSTALL_PREFIX)/bin

exec/hetmers: src_ploidyplot/PloidyPlot.c src_ploidyplot/libfastk.c src_ploidyplot/libfastk.h src_ploidyplot/matrix.c src_ploidyplot/matrix.h
	gcc $(CFLAGS) -o $@ src_ploidyplot/PloidyPlot.c src_ploidyplot/libfastk.c src_ploidyplot/matrix.c -lpthread -lm

exec/extract_kmer_pairs: src_ploidyplot/PloidyList.c src_ploidyplot/libfastk.c src_ploidyplot/libfastk.h src_ploidyplot/matrix.c src_ploidyplot/matrix.h
	gcc $(CFLAGS) -o $@ src_ploidyplot/PloidyList.c src_ploidyplot/libfastk.c src_ploidyplot/matrix.c -lpthread -lm

.PHONY : clean
clean :
	rm -f exec/hetmers
