all: numeric

numeric:
	$(MAKE) -C nutils/numeric test_c

test:
	nosetests -v tests

clean:
	$(MAKE) -C nutils/numeric clean

.PHONY: all numeric test clean

# vim:noexpandtab
