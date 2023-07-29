all:
	mamba env create -f SQANTI3.conda_env.yml -n SQ3
	mamba activate SQ3
	cd cDNA_Cupcake && $(MAKE)

