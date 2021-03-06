.PHONY: inplace inplace-force clean docs test test-all test-cli build-data develop clean-data help docs profile
.DEFAULT_GOAL := help

PYTHON ?= python

PROFILE_FILES := $(wildcard analysis/profiling/*.pyopts)
PROFILE_TGT := $(PROFILE_FILES:.pyopts=.txt)

inplace:  ## Build extensions in place
	$(PYTHON) setup.py build_ext --inplace

inplace-force:  ## Build extensions in place (force0
	$(PYTHON) setup.py build_ext --inplace -f

test: inplace  ## Test the package with the current python version
	pip install 'pytest'
	pytest

test-cli:  ## Test the command-line interface
	@$(MAKE) --no-print-directory -C tests/cli test

test-all: clean test-cli  ## Test the package with multiple python environments using tox
	pip install 'tox>=2.7'
	tox

develop: inplace  ## Prepare the package for active development
	$(PYTHON) setup.py develop

clean:  ## Safely clean compiled package files, docs and test files
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name __pycache__ -type d -exec rm -rf {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '*.so' -exec rm -f {} +
	find analysis/ -name '*.txt' -exec rm -f {} +
	find analysis/ -name '*.svg' -exec rm -f {} +
	find analysis/ -name '*.png' -exec rm -f {} +
	rm -rf .tox
	$(MAKE) -C docs clean

build-data: inplace-force  ## Build the datasets
	pip install 'tqdm>=4.8'
	python setup.py build_data

clean-data:  ## Clean the datasets
	find mollib/data/*statistics/ -type f -exec rm -f {} +

build-analysis: build-data  ## Conduct the analysis of datasets
	pip install 'matplotlib>=2.0'
	cd analysis/datasets
	python ./ramachandran
	find analysis/ -name "*lowres.png" -exec cp '{}' docs/cli/img/ \;

docs: develop  ## Build the documentation
	pip install 'Sphinx>=1.5'
	pip install 'sphinxcontrib-napoleon2>=1.0'
	find analysis/ -name "*lowres.png" -exec cp '{}' docs/cli/img/ \;
	$(MAKE) -C docs cli html latexpdf
	@echo "\033[92m\n\nBuild successful!\033[0m"
	@echo "\033[92mView the html docs at docs/_build/html/index.html.\033[0m"
	@echo "\033[92mView the pdf docs docs/_build/latex/.\n\033[0m"

install:  ## Install mollib
	python setup.py install

#publish:  # Publish the sdist to pypi
# pip install 'twine>=1.5.0'

analysis/profiling/%.txt: analysis/profiling/%.pyopts
	python -m cProfile -s tottime `which mollib` `cat $<` > $@

profile: $(PROFILE_TGT)  ## Build the profile reports (analysis/profiling)
	@grep -A15 "function calls" analysis/profiling/*.txt | grep -v summary.txt > analysis/profiling/summary.txt
	@echo "\033[92mProfile reports built under analysis/profiling\033[0m"

help:  ## Print this help message
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-15s\033[0m %s\n", $$1, $$2}'
