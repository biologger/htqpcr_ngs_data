
The full data analysis can be viewed as html export:

[htqpcr_ngs_comparison.html](https://github.com/biologger/htqpcr_ngs_data/blob/main/htqpcr_ngs_comparison.html)

and as html with embeded Supplementary Figures:

[htqpcr_ngs_comparison_R.html](https://github.com/biologger/htqpcr_ngs_data/blob/main/htqpcr_ngs_comparison_R.html)

Linux/Mac:

Clone repository

		$ git clone --recurse-submodules https://github.com/biologger/htqpcr_ngs_data

		# or alternatively
		$ git clone https://github.com/biologger/htqpcr_ngs_data.git
		$ git submodule init
		$ git submodule update


Install all dependencies with [Anaconda](https://www.anaconda.com/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html)

		$ conda env create -f conda_analysis_env.yml

Activate environment

		$ conda activate htqpcr_ngs_comparison

Start Jupyter Notebook

		$ jupyter-notebook

Open and Run the notebook

* htqpcr_ngs_comparison_R.ipynb

Windows:

We recommend a linux docker container or a linux OS in a virtual machine.
For Windows install R Studio and Rtools manually.
