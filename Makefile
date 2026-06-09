SHELL := /bin/bash
XELATEX ?= xelatex

SUBDIRS := final_project/Qian week10 week9
SUMMARY := AI4Math_Project_Summary.pdf

.PHONY: all clean rebuild summary $(SUBDIRS)

all: final_project/Qian week10 week9

final_project/Qian:
	$(MAKE) -C final_project/Qian all

week9:
	$(MAKE) -C week9 all

week10:
	$(MAKE) -C week10 all

summary:
	$(XELATEX) -interaction=nonstopmode -jobname=AI4Math_Project_Summary project_summary.tex

clean:
	for dir in $(SUBDIRS); do $(MAKE) -C $$dir clean; done
	rm -f $(SUMMARY) AI4Math_Project_Summary.aux AI4Math_Project_Summary.log AI4Math_Project_Summary.out

rebuild: clean all
