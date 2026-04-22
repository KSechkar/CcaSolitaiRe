# CcaSolitaiRe
Code and experimental data accompanying Section 2 of Chapter 5 of 
Kirill Sechkar's D.Phil thesis 'Dynamical in vivo characterisation and design of resource-aware genetic circuits'. 
The repository includes experimental data for characterising different DNA cvonstructs on the way towards 
a single-plasmid implementation of CcaSR the optogenetic system.

## File organisation
The repository is organised as follows:
- Folder _Microplate experiments_:
  - _CcaSR_Ptac.ipynb_: Jupyter notebook for processing the data obtained by Kirill Sechkar
  whilst testing the CcaSR system's dependence on the transcriptional context
  by considering its implementations with and without an additional 
  inducible system for mScarlet-I expression controlled by a Ptac promoter. 
  The source data for this are found in the folder _CcaSR_Ptac_context_data_.
  - _014_testing.ipynb_: Jupyter notebook for processing the data obtained by Yuhang Xie
  whilst testing the pKS-02-014 constructs in non-DH5alpha strains in a microwell plate. 
  The source data for this are found in the folder _014_testing_.
  - _microplate_antools.py_: Python script containing functions for the analytical tools
  used in the Jupyter notebooks for processing microplate data.
- Folder _Screening experiments_:
  - _Screening.xlsx._: excel workbook sotring, processing and plotting data from
  construct screening experiments obtained by Kirill Sechkar,
  who incubated cells overnight in 2 mL cultures under different lighting conditions
  and measured the fluorescence of 200 uL samples 24 h after the start of the experiment.

## System requirements
Python code should be run using Python 3.12.3 with numpy 1.26.4, bokeh 3.4.1, jupyter 1.1.1,  and jupyterlab 4.2.5.
The xlsx file can be worked with using Microsoft Excel Version 2603.
