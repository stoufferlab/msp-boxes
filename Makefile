
all:
	$(MAKE) -C src

clean:
	$(MAKE) clean -C src
	$(MAKE) clean -C test

check:
	$(MAKE) -C test
