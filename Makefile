# PATH for libraries is guessed
CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

ifndef INSTALL_PREFIX
    INSTALL_PREFIX = /usr/local
endif

HET_KMERS_INST = $(INSTALL_PREFIX)/bin/hetmers $(INSTALL_PREFIX)/bin/extract_kmer_pairs

.PHONY : default
default: exec/hetmers exec/extract_kmer_pairs

.PHONY : install
install : $(HET_KMERS_INST) 

$(INSTALL_PREFIX)/bin/% : exec/%
	install -C $< $(INSTALL_PREFIX)/bin

exec/hetmers: src/lib/PloidyPlot.c src/lib/libfastk.c src/lib/libfastk.h src/lib/matrix.c src/lib/matrix.h
	gcc $(CFLAGS) -o $@ src/lib/PloidyPlot.c src/lib/libfastk.c src/lib/matrix.c -lpthread -lm

exec/extract_kmer_pairs: src/lib/PloidyList.c src/lib/libfastk.c src/lib/libfastk.h src/lib/matrix.c src/lib/matrix.h
	gcc $(CFLAGS) -o $@ src/lib/PloidyList.c src/lib/libfastk.c src/lib/matrix.c -lpthread -lm

.PHONY : clean
clean :
	rm -f hetmers
	rm -f extract_kmer_pairs
