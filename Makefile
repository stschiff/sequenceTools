all : build/filterTrimFastq build/vcf2eigenstrat build/genToVcf

build/filterTrimFastq : filterTrimFastq.d gzip.d
	dmd -O filterTrimFastq.d gzip.d -odbuild -ofbuild/filterTrimFastq

build/vcf2eigenstrat : vcf2eigenstrat.hs
	ghc -outputdir build -o build/vcf2eigenstrat vcf2eigenstrat.hs

build/genToVcf : genToVcf.hs
	ghc -outputdir build -o build/genToVcf genToVcf.hs

clean :
	rm build/*
	
.PHONY: clean