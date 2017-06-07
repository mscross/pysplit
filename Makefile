all: clean develop test

test:
	nosetests -s -v pysplit

coverage: clean-cov
	nosetests pysplit --with-coverage --cover-package=pysplit

develop:
	python setup.py develop

clean-pyc:
	find pysplit -name "*.pyc" | xargs rm -f

clean-build:
	rm -rf ./build

clean-version:
	find pysplit -name "*version.py" | xargs rm -f

clean-cov:
	rm -rf ./coverage ./.coverage ./htmlcov

clean: clean-build clean-pyc clean-version clean-cov
