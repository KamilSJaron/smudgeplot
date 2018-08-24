# PATH for libraries is guessed
PACKAGE = smudgeplot
$(eval PATH := $(shell Rscript -e "noquote(.libPaths())" | tail -1 | cut -f 2 -d ' '))

INSTALLATION = $(PATH)/$(PACKAGE)

.PHONY : install
install : $(INSTALLATION)

$(INSTALLATION) : R/*
	Rscript install.R