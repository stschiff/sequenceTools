build/filterTrimFastq : filterTrimFastq.d gzip.d
	dmd -O filterTrimFastq.d gzip.d -odbuild -ofbuild/filterTrimFastq

clean :
	rm build/*
	
.PHONY: clean