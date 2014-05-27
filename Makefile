DEBUG = -D_$(DF)
export DEBUG

all: magrheol

.PHONY: magrheol

magrheol: 
	$(MAKE) -C src/

clean: 
	$(MAKE) clean -C src/ 
