Some modification to FinaleToolkit are here. 

Please note that the code for delfi, delfi-gc-correct, delfi-merge-bins, and end-motifs still needs to be ported over.

Other than that, the intention for this branch is to optimize user experience, improve code efficiency, and make some minor changes. 

There is no documentation associated with this branch.

#### Installing
```
conda create -n r_finale python=3.11
conda activate r_finale
mkdir r_finale
cd r_finale
git clone -b ravi https://github.com/epifluidlab/FinaleToolkit.git
pip install -e ./FinaleToolkit/
finaletools
```
