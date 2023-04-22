MAKEFLAGS += --warn-undefined-variables

SHELL := bash
.SHELLFLAGS := -eu -o pipefail -c

.DELETE_ON_ERROR:
.SUFFIXES:

package := nrpytutorial
pyseed ?= $(shell command -v python3 2> /dev/null)
venv := $(shell echo $${VIRTUAL_ENV-.venv})
VIRTUAL_ENV ?= ''
python := $(venv)/bin/python
pip := $(python) -m pip
FILE ?=
ARGS ?=
SHOULD_JUNIT_REPORT ?=

comma := ,
empty :=
space := $(empty) $(empty)
# now use to join elements of a list
# foo:= a b c
# bar:= $(subst $(space),$(comma),$(foo))

.DEFAULT_GOAL := help

.PHONY: help ### show this menu
help:
	@sed -nr '/#{3}/{s/\.PHONY:/--/; s/\ *#{3}/:/; p;}' ${MAKEFILE_LIST}

FORCE:

reveal-%: FORCE
	@echo $($*)

append-gitignore-%: FORCE
	@grep --quiet --line-regexp --fixed-strings $* .gitignore 2> /dev/null \
	|| echo $* >> .gitignore
	sort -o .gitignore{,}

remove-gitignore-%: FORCE
	sed --in-place '/$*/d' .gitignore
	sort -o .gitignore{,}

src-dir := src
test-dir := tests
bin-dir := bin
notebook-dir := $(bin-dir)/notebooks
build-dir := build
nbconvert-dir := $(build-dir)/nbconvert
stamp-dir := .stamps
doc-dir := docs
stamp-nrpy-dir := $(stamp-dir)/nrpy
diagram-dir := $(doc-dir)/diagrams
fig-dir := $(build-dir)/figs
dist-dir := dist
data-dir := data

maintainer-dirs := $(venv) $(build-dir) $(stamp-dir) $(dist-dir)

$(src-dir) $(test-dir) $(build-dir) $(stamp-dir) $(stamp-nrpy-dir) \
$(doc-dir) $(diagram-dir) $(fig-dir) $(data-dir) \
$(bin-dir) $(notebook-dir) $(nbconvert-dir):
	@mkdir --parents $@

.PHONY: stampdir
stampdir: $(stamp-dir) append-gitignore-$(stamp-dir)/

.PHONY: has-python-seed
has-python-seed:
	@if [[ -z $(pyseed) ]]; then \
		echo 'No python seed found, possible resolutions:'; \
		echo '-- pass inline, e.g. pyseed=/path/python make <goal>';\
		echo '-- overwrite <pyseed> variable in Makefile';\
		exit 2; \
	fi

.PHONY: is-venv-inactive
is-venv-inactive:
	@if [[ ! -z $(VIRTUAL_ENV) ]]; then \
		echo 'Python virtual environment ACTIVE'; \
		echo 'Run: deactivate'; \
		exit 2; \
	fi

venv-stamp := $(stamp-dir)/venv.stamp
$(venv-stamp): | $(stamp-dir)
	$(pyseed) -m venv $(venv)
	$(pip) install --upgrade pip > /dev/null
	$(pip) install --upgrade build > /dev/null
	touch $@

.PHONY: install-venv ### install virtual environment using python seed
install-venv: has-python-seed
install-venv: is-venv-inactive
install-venv: stampdir $(venv-stamp) append-gitignore-$(venv)/

.PHONY: uninstall-venv
uninstall-venv: is-venv-inactive remove-gitignore-$(venv)
	rm --force --recursive $(venv) $(venv-stamp)

.PHONY: is-venv-present
is-venv-present:
	@if [[ ! -e $(python) ]]; then \
		echo 'No virtual environment found'; \
		echo 'Run: make install-venv'; \
		exit 2; \
	fi

.PHONY: venv ### virtual environment control info
venv: is-venv-present
	@echo "Active shell: $$0"
	@echo "Command to activate virtual environment:"
	@echo "- bash/zsh: source $(venv)/bin/activate"
	@echo "- fish: source $(venv)/bin/activate.fish"
	@echo "- csh/tcsh: source $(venv)/bin/activate.csh"
	@echo "- PowerShell: $(venv)/bin/Activate.ps1"
	@echo "Exit: deactivate"

requirements-conf := requirements.txt
$(requirements-conf):
	@echo "pytest" >> $@
	@echo "pytest-mock" >> $@
	@echo "pytest-cov" >> $@
	@echo "pytest-datafiles" >> $@
	@echo "pylint" >> $@
	@echo "pylint-junit" >> $@
	@echo "autopep8" >> $@
	@echo "pynvim" >> $@

requirements-stamp := $(stamp-dir)/$(requirements-conf).stamp
$(requirements-stamp): $(requirements-conf) $(venv-stamp) | $(stamp-dir)
	$(pip) install --upgrade --requirement $<
	touch $@

.PHONY: install-requirements
install-requirements: is-venv-present $(requirements-stamp)

.PHONY: uninstall-requirements
uninstall-requirements: is-venv-present
	[[ ! -e $(requirements-conf) ]] || $(pip) uninstall --requirement $(requirements-conf) --yes
	rm --force $(requirements-stamp)

.PHONY: clean-requirements
clean-requirements: uninstall-requirements
	rm --force $(requirements-conf) $(requirements-stamp)

.PHONY: setup ### install venv and its requirements
setup: install-venv install-requirements

pyprojectrc := pyproject.toml README.md LICENSE
pyproject.toml:
	@echo '[build-system]' >> $@
	@echo 'requires = ["setuptools"]' >> $@
	@echo 'build-backend = "setuptools.build_meta"' >> $@
	@echo '' >> $@
	@echo '[project]' >> $@
	@echo 'name = "$(package)"' >> $@
	@echo 'version = "0.0.1"' >> $@
	@echo 'requires-python = ">=3.7"' >> $@
	@echo 'dependencies = []' >> $@

README.md:
	@echo '# $(package)' >> $@
	@echo 'Excellent package with much to offer' >> $@
	@echo '## Quickstart' >> $@
	@echo 'clone and `pyseed=/path/python make development`' >> $@

LICENSE:
	@echo 'MIT License' >> $@

pyproject-stamp := $(stamp-dir)/pyproject.stamp
$(pyproject-stamp): $(pyprojectrc) | $(stamp-dir)
	touch $@

.PHONY: install-pyproject
install-pyproject: $(pyproject-stamp)

.PHONY: uninstall-pyproject
uninstall-pyproject:
	rm --force $(pyprojectrc) $(pyproject-stamp)

package-stamp := $(stamp-dir)/package.stamp
$(package-stamp): $(pyproject-stamp) $(venv-stamp) | $(stamp-dir)
	$(pip) install --editable .
	touch $@

package-egg := $(package).egg-info
.PHONY: install-package
install-package: is-venv-present $(package-stamp) append-gitignore-$(package-egg)/

.PHONY: is-package-available
is-package-available: is-venv-present
	@$(python) -c 'import $(package)'

.PHONY: uninstall-package
uninstall-package: is-venv-present is-package-available remove-gitignore-$(package-egg)
	$(pip) uninstall $(package) --yes
	rm --force --recursive $(shell find . -type d -name "$(package-egg)") $(package-stamp)

.PHONY: development ### setup development package environment
development: setup install-package

.PHONY: install-dirs
install-dirs: $(src-dir) $(test-dir) $(stamp-dir)

.PHONY: uninstall-dirs
uninstall-dirs: remove-gitignore-$(stamp-dir)
	rm --recursive --force $(src-dir) $(test-dir) $(stamp-dir)

sample-package := $(src-dir)/$(package).py
$(sample-package): | $(src-dir)
	@echo "creating sample: $@"
	@echo "def sample(): return 0" >> $@

sample-tests := $(test-dir)/test_$(package).py
$(sample-tests): | $(test-dir)
	@echo "creating sample: $@"
	@echo "import pytest" >> $@
	@echo "from $(package) import sample" >> $@
	@echo "def test_scenario_1(): assert sample() == 0" >> $@
	@echo "def test_scenario_2(): assert not sample() != 0" >> $@

.PHONY: install-sample
install-sample: install-dirs $(sample-package) $(sample-tests)

.PHONY: uninstall-sample
uninstall-sample:
	rm --force $(sample-package) $(sample-tests)

.PHONY: sample ### create sample project and setup into development mode
sample: install-sample development

.PHONY: run ### run particular <FILE>, optionally pass some <ARGS>
run: is-package-available
	$(python) $(FILE) $(ARGS)

doctest-ignore := $(shell find -maxdepth 1 -mindepth 1 -type d -not -name $(src-dir))
.PHONY: doctest
doctest-module := pytest
doctest-module += --quiet
doctest-module += -rfE
doctest-module += --showlocals
doctest-module += $(addprefix --ignore=,$(doctest-ignore))
doctest-module += --doctest-modules

ifdef SHOULD_JUNIT_REPORT
doctest-module += --junit-xml=test-results/doctests/results.xml
endif

doctest: is-package-available
	$(python) -m $(doctest-module) $(FILE) || true

unittest-ignore := $(shell find -maxdepth 1 -mindepth 1 -type d -not -name $(test-dir))
.PHONY: unittest
unittest-module := pytest
unittest-module += --quiet
unittest-module += -rfE
unittest-module += --showlocals
unittest-module += $(addprefix --ignore=,$(unittest-ignore))

ifdef SHOULD_JUNIT_REPORT
unittest-module += --junit-xml=test-results/unittests/results.xml
endif

unittest: is-package-available
	$(python) -m $(unittest-module) $(FILE)

.PHONY: lint
lint-module := pylint
lint-module += --fail-under=7.5

ifdef FILE
lint-runner := $(python) -m $(lint-module) $(FILE)
else
lint-runner := $(python) -m $(lint-module) $(src-dir)/
endif

ifdef SHOULD_JUNIT_REPORT
lint-runner += --output-format=pylint_junit.JUnitReporter
lint-runner += > test-results/lint/results.xml
endif

lint: is-package-available
ifdef SHOULD_JUNIT_REPORT
	mkdir --parents test-results/lint/
endif
	$(lint-runner)

.PHONY: test
test: doctest unittest

.PHONY: tests-structure
tests-structure:
	find $(src-dir)/$(package) -type f -name '*.py' \
		| sed -rn "s/$(src-dir)\/$(package)/$(test-dir)/; s/.py//p" \
		| xargs mkdir --parents

.PHONY: coverage
coverage-module := pytest
coverage-module += --cov=$(src-dir)
coverage-module += --cov-branch
coverage-module += --cov-fail-under=75
coverage-module += --doctest-modules

ifdef SHOULD_JUNIT_REPORT
coverage-module += --cov-report=xml:test-results/coverage/report.xml
endif

ifdef SHOULD_HTML_REPORT
coverage-module += --cov-report=html
endif

coverage-runner := $(python) -m $(coverage-module)
coverage: is-package-available append-gitignore-.coverage
	$(coverage-runner)

.PHONY: check ### doctest, unittest, lint and coverage whole project
check: test coverage

.PHONY: format ### auto format into pep8 FILE=path/file or all
formatter-module := autopep8
formatter-module += --in-place
formatter-module += --aggressive --aggressive
formatter-module += --hang-closing --experimental

ifdef FILE
formatter-runner := $(python) -m $(formatter-module) $(FILE)
else
formatter-runner := $(python) -m $(formatter-module) --recursive $(src-dir)/ $(test-dir)/
endif

format: is-package-available
	$(formatter-runner)

.PHONY: dist ### create distribution package
build-module := build
build-runner := $(python) -m $(build-module) > /dev/null

dist: $(pyprojectrc) test append-gitignore-$(dist-dir)/
	$(build-runner)

.PHONY: distclean
distclean: remove-gitignore-$(dist-dir)
	rm --force --recursive $(dist-dir)

vpath %.puml $(diagram-dir)

diagram-srcs := $(wildcard $(diagram-dir)/*.puml)
diagram-goals := $(addprefix diagram-,$(notdir $(basename $(diagram-srcs))))
fig-dirname := $(notdir $(fig-dir))

$(fig-dir)/%.png: %.puml | $(fig-dir)
	@echo "-- creating $@"
	plantuml $< -o $(notdir $(fig-dir))
	touch $@

.PHONY: $(diagram-goals)
$(diagram-goals): diagram-%: $(fig-dir)/%.png | append-gitignore-$(fig-dirname)/

.PHONY: diagrams
diagrams: $(diagram-goals)

.PHONY: clean-$(fig-dirname)
clean-$(fig-dirname): remove-gitignore-$(fig-dirname)
	rm --force --recursive $(fig-dir)

ipython-stamp := $(stamp-dir)/ipython.stamp
$(ipython-stamp): | $(stamp-dir)
	$(pip) install ipython
	touch $@

.PHONY: run-ipython
run-ipython: is-venv-present $(ipython-stamp)
	$(python) -m IPython --colors Linux

notebook-stamp := $(stamp-dir)/notebook.stamp
$(notebook-stamp): | $(stamp-dir)
	$(pip) install notebook
	touch $@

jupyter-runner := $(venv)/bin/jupyter

vpath %.ipynb $(notebook-dir)

notebook-srcs := $(wildcard $(notebook-dir)/*.ipynb)
notebook-goals := $(addprefix run-notebook-,$(notdir $(basename $(notebook-srcs))))
nbconvert-dirname := $(notdir $(nbconvert-dir))

$(nbconvert-dir)/%.ipynb: %.ipynb | $(nbconvert-dir)
	@echo "-- creating $@"
	$(jupyter-runner) nbconvert --output-dir=$(nbconvert-dir) \
		--execute --to notebook --allow-errors $<
	touch $@

.PHONY: $(notebook-goals)
$(notebook-goals): run-notebook-%: is-venv-present $(notebook-stamp) $(nbconvert-dir)/%.ipynb | append-gitignore-$(nbconvert-dirname)/

.PHONY: convert-notebook
convert-notebook:
	@echo '-- converting into script $@'
	$(jupyter-runner) nbconvert --output-dir=$(nbconvert-dir) \
		--to script $(FILE)

.PHONY: clean-$(nbconvert-dirname)
clean-$(nbconvert-dirname): remove-gitignore-$(nbconvert-dirname)
	rm --force --recursive $(nbconvert-dir)

.PHONY: run-jupyter
run-jupyter: is-venv-present $(notebook-stamp)
	$(jupyter-runner) notebook

.PHONY: TAGS
maintainer-dirs += tags
TAGS: append-gitignore-tags
	ctags --languages=python --recurse

.PHONY: clean-TAGS
clean-TAGS:
	rm --force tags

nrpy-modules := outputC.py loop.py NRPy_param_funcs.py SIMD.py
nrpy-modules += expr_tree.py cse_helpers.py finite_difference.py
nrpy-modules += finite_difference_helpers.py grid.py indexedexp.py
nrpy-modules += functional.py reference_metric.py
nrpy-modules += cmdline_helper.py
nrpy-stamps := $(patsubst %,$(stamp-nrpy-dir)/%.stamp,$(nrpy-modules))

python-version := $(shell $(python) --version 2> /dev/null | grep -Po '\d+.\d+')
python-site-packages := $(venv)/lib/python$(python-version)/site-packages
$(stamp-nrpy-dir)/%.stamp: | $(stamp-nrpy-dir)
	@if [[ ! -e $(abspath ./$*) ]]; then \
		echo 'Not found nrpy module ./$*'; \
		exit 2; \
	fi
	ln -s $(abspath ./$*) $(python-site-packages)/.
	mkdir -p $(dir $@)
	touch $@

.PHONY: install-nrpy-modules ### install nrpy modules into project venv
install-nrpy-modules: is-venv-present $(nrpy-stamps)

.PHONY: uninstall-nrpy-modules
uninstall-nrpy-modules:
	rm --force $(addprefix $(python-site-packages)/,$(notdir $(nrpy-modules)))
	rm --force --recursive $(stamp-nrpy-dir)

nrpy-progress-sfcollapse-modules := Deprecated/deprecated_reference_metric.py
nrpy-progress-sfcollapse-modules += Deprecated/MoLtimestepping
nrpy-progress-sfcollapse-modules += Deprecated/BSSN
nrpy-progress-sfcollapse-modules += Deprecated/deprecated_NRPy_param_funcs.py
nrpy-progress-sfcollapse-modules += Deprecated/CurviBoundaryConditions
nrpy-progress-sfcollapse-stamps := $(patsubst %,$(stamp-nrpy-dir)/%.stamp,$(nrpy-progress-sfcollapse-modules))

nrpy-progress-sfcollapse-home := in_progress-NRPyCritCol
.PHONY: install-nrpy-progress-sfcollapse-simd
install-nrpy-progress-sfcollapse-simd:
	ln -sf $(abspath SIMD) $(nrpy-progress-sfcollapse-home)/.
	@grep --quiet --line-regexp --fixed-strings $(nrpy-progress-sfcollapse-home)/SIMD .gitignore 2> /dev/null \
				|| echo $(nrpy-progress-sfcollapse-home)/SIMD >> .gitignore
	sort -o .gitignore{,}

.PHONY: uninstall-nrpy-progress-sfcollapse-simd
uninstall-nrpy-progress-sfcollapse-simd:
	rm --force $(nrpy-progress-sfcollapse-home)/SIMD
	sed -ri "/$(nrpy-progress-sfcollapse-home)\/SIMD/d" .gitignore
	sort -o .gitignore{,}

nrpy-progress-sfcollapse-bssn-home := Deprecated/BSSN
nrpy-progress-sfcollapse-bssn-modules := BSSN_quantities.py BSSN_in_terms_of_ADM.py BSSN_RHSs.py
nrpy-progress-sfcollapse-bssn-modules += BSSN_gauge_RHSs.py ADM_in_terms_of_BSSN.py
nrpy-progress-sfcollapse-bssn-modules += ADMBSSN_tofrom_4metric.py BSSN_stress_energy_source_terms.py
nrpy-progress-sfcollapse-bssn-modules += Enforce_Detgammahat_Constraint.py
nrpy-progress-sfcollapse-bssn-modules += BSSN_constraints.py BSSN_T4UUmunu_vars.py
nrpy-progress-sfcollapse-bssn-modules := $(addprefix $(nrpy-progress-sfcollapse-bssn-home)/,$(nrpy-progress-sfcollapse-bssn-modules))

$(nrpy-progress-sfcollapse-bssn-modules):
	ln -sf $(abspath $(subst $(nrpy-progress-sfcollapse-bssn-home),BSSN,$@)) $(nrpy-progress-sfcollapse-bssn-home)/.
	echo $@ >> .gitignore
	sort -o .gitignore{,}

.PHONY: install-nrpy-progress-sfcollapse-bssn
install-nrpy-progress-sfcollapse-bssn: $(nrpy-progress-sfcollapse-bssn-modules)

nrpy-progress-sfcollapse-bssn-gitignore-clean := $(foreach module,$(nrpy-progress-sfcollapse-bssn-modules),/$(subst /,\/,$(module))/d;)
.PHONY: uninstall-nrpy-progress-sfcollapse-bssn
uninstall-nrpy-progress-sfcollapse-bssn:
	rm --force $(nrpy-progress-sfcollapse-bssn-modules)
	sed -ri '$(nrpy-progress-sfcollapse-bssn-gitignore-clean)' .gitignore
	sort -o .gitignore{,}

.PHONY: install-nrpy-progress-sfcollapse-modules
install-nrpy-progress-sfcollapse-modules: is-venv-present
install-nrpy-progress-sfcollapse-modules: install-nrpy-modules
install-nrpy-progress-sfcollapse-modules: install-nrpy-progress-sfcollapse-simd
install-nrpy-progress-sfcollapse-modules: install-nrpy-progress-sfcollapse-bssn
install-nrpy-progress-sfcollapse-modules: $(nrpy-progress-sfcollapse-stamps)

.PHONY: uninstall-nrpy-progress-sfcollapse-modules
uninstall-nrpy-progress-sfcollapse-modules: uninstall-nrpy-progress-sfcollapse-simd
uninstall-nrpy-progress-sfcollapse-modules: uninstall-nrpy-progress-sfcollapse-bssn
	rm --force $(addprefix $(python-site-packages)/,$(notdir $(nrpy-progress-sfcollapse-modules)))
	rm --force --recursive $(nrpy-progress-sfcollapse-stamps)

.PHONY: maintainer-clean ### clean Makefile can reconstruct
maintainer-clean:
	rm --force --recursive $(maintainer-dirs)

.PHONY: clean ### clean everything WITHOUT Makefile and .git
clean:
	find . -maxdepth 1 \
		-name ${MAKEFILE_LIST} -prune \
		-o -name .git -prune \
		-o -type d -name '.' \
		-o -print0 | xargs -0 rm --force --recursive
