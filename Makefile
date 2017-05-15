.PHONY: inplace clean docs

PYTHON ?= python

inplace:
	$(PYTHON) setup.py build_ext --inplace -f

test: inplace
	pip install 'nose>=1.3'
	nosetests

testall: clean
	pip install 'tox>=2.7'
	tox

develop: inplace
	$(PYTHON) setup.py develop

clean:
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '*.so' -exec rm -f {} +
	rm -rf .tox
	$(MAKE) -C docs clean

builddata: inplace
	pip install 'tqdm>=4.8'
	python setup.py build_data

cleandata:
	find mollib/data/*statistics/ -type f -exec rm -f {} +

docs: develop
	pip install 'Sphinx>=1.5'
	pip install 'sphinxcontrib-napoleon2>=1.0'
	$(MAKE) -C docs cli html latexpdf
	@echo "\033[92m\n\nBuild successful!\033[0m"
	@echo "\033[92mView the html docs at docs/_build/html/index.html.\033[0m"
	@echo "\033[92mView the pdf docs docs/_build/latex/.\n\033[0m"

#publish
# pip install 'twine>=1.5.0'
