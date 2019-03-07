ifndef INSTAL_PREFIX
	INSTAL_PREFIX = /usr/local
endif

BIN = $(INSTAL_PREFIX)/bin/smudgeplot

.PHONY : install
install : $(BIN)

# alternativelly pytest
.PHONY : test
test :
	python3 -m unittest


$(BIN) : setup.py smudgeplot/* scripts/*
	python3 $< install

