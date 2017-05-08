.PHONY: clean-pyc docs

clean-pyc:
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +

docs:
	$(MAKE) -C docs cli html latexpdf
	@echo "\033[95m\n\nBuild successful!\033[0m"
	@echo "\033[95mView the html docs at docs/_build/html/index.html.\033[0m"
	@echo "\033[95mView the pdf docs docs/_build/latex/.\n\033[0m"
