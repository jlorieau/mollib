.PHONY: inplace clean docs

PYTHON ?= python

inplace:
	$(PYTHON) setup.py build_ext --inplace

test: inplace
	nosetests

testfull: clean
	tox

develop: clean inplace
	$(PYTHON) setup.py develop

clean:
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '*.so' -exec rm -f {} +
	rm -rf .tox
	$(MAKE) -C docs clean

docs:
	$(MAKE) -C docs cli html latexpdf
	@echo "\033[92m\n\nBuild successful!\033[0m"
	@echo "\033[92mView the html docs at docs/_build/html/index.html.\033[0m"
	@echo "\033[92mView the pdf docs docs/_build/latex/.\n\033[0m"
