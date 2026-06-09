PYTHON   ?= python
TEX      = xelatex
TEXFLAGS = -interaction=nonstopmode -output-directory=build
BIB      = biber
TARGET   = report.pdf
SRC      = report.tex
BIBFILE  = Qian.bib

.PHONY: all clean m1 m2 m3 m4 m5 m6 m7 m8

all: m8 $(TARGET)

m1:
	$(PYTHON) src/patched_conics.py --rp-au 0.2 --check

m2:
	$(PYTHON) src/nbody.py --method verlet --steps 2000 --check --convergence

m3:
	$(PYTHON) src/m3_ephemeris_check.py --dt 1800 --check

m4:
	$(PYTHON) src/flyby.py --v-inf 1.6 --periapsis-altitude 100 --check

m5:
	$(PYTHON) src/m5_single_day.py --check

m6: m5
	$(PYTHON) src/m6_launch_window.py --check

m7: m6
	$(PYTHON) src/m7_sensitivity.py --check

m8: m3 m7
	$(PYTHON) src/m8_visualize.py --check

$(TARGET): $(SRC) $(BIBFILE) | build
	$(TEX) $(TEXFLAGS) $<
	$(PYTHON) -c "from pathlib import Path; import shutil; shutil.copyfile('$(BIBFILE)', Path('build') / '$(BIBFILE)')"
	cd build && $(BIB) report
	$(TEX) $(TEXFLAGS) $<
	$(TEX) $(TEXFLAGS) $<
	$(PYTHON) -c "from pathlib import Path; import shutil; shutil.copyfile(Path('build') / '$@', '$@')"

build:
	$(PYTHON) -c "from pathlib import Path; Path('$@').mkdir(exist_ok=True)"

clean:
	$(PYTHON) -c "from pathlib import Path; import shutil; shutil.rmtree('build', ignore_errors=True); Path('$(TARGET)').unlink(missing_ok=True)"
