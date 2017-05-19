.PHONY: inplace clean docs test test-all build-data develop clean-data help docs
.DEFAULT_GOAL := help

PYTHON ?= python

inplace: ## Build extensions in place
	$(PYTHON) setup.py build_ext --inplace -f

test: inplace  ## Test the package with the current python version
	pip install 'nose>=1.3'
	nosetests

test-all: clean  ## Test the package with multiple python environments using tox
	pip install 'tox>=2.7'
	tox

develop: inplace  ## Prepare the package for active development
	$(PYTHON) setup.py develop

clean:  ## Clean compiled package files, docs and test files
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '*.so' -exec rm -f {} +
	rm -rf .tox
	$(MAKE) -C docs clean

build-data: inplace  ## Build the datasets
	pip install 'tqdm>=4.8'
	python setup.py build_data

clean-data:  ## Clean the datasets
	find mollib/data/*statistics/ -type f -exec rm -f {} +

build-analysis: build-data  ## Conduct the analysis of datasets
	pip install 'matplotlib>=2.0'
	cd analysis/datasets
	python ./ramachandran

docs: develop build-analysis  ## Build the documentation
	pip install 'Sphinx>=1.5'
	pip install 'sphinxcontrib-napoleon2>=1.0'
	find analysis/ -name "*lowres.png" -exec cp '{}' docs/cli/img/ \;
	$(MAKE) -C docs cli html latexpdf
	@echo "\033[92m\n\nBuild successful!\033[0m"
	@echo "\033[92mView the html docs at docs/_build/html/index.html.\033[0m"
	@echo "\033[92mView the pdf docs docs/_build/latex/.\n\033[0m"

install:  ## Install mollib
	python setup.py install

#publish
# pip install 'twine>=1.5.0'

help:  ## Print this help message
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-15s\033[0m %s\n", $$1, $$2}'
