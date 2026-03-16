.PHONY: all lint test test-docker install dev clean distclean

PYTHON ?= python

all: ;

lint:
	q2lint
	flake8

test: all
	py.test

test-cov: all
	python -m pytest --cov=q2_subread -n 4 && coverage xml -o coverage.xml

test-docker: all
	qiime info
	qiime subread --help

install: all
	$(PYTHON) -m pip install -v .

dev: all
	pip install -e .

clean: distclean

distclean: ;
