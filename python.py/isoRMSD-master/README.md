# isoRMSD

- isoRMSD can calculate RMSD between 2 conformers with different atom names.
- The atom order of molecules also don't need to be same.
- Hydration or not is OK.
- This script is for small molecule only. It will be very slow for protein.

This script is base on RDKit cookbook -- [RMSD Calculation between N molecules](http://www.rdkit.org/docs/Cookbook.html).


### 1. Download and install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.continuum.io/)
Install command is simple!
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
# All enter 'yes' is OK.
```

### 2. Install [RDKit](http://rdkit.org/)
```
conda create -n rmsd -c conda-forge -c oddt rdkit oddt
conda activate rmsd

# Or conda install -c conda-forge -c oddt rdkit oddt
```

### 3. Usage
```console
(rmsd) [jcyang@x038:isoRMSD]$ python isoRMSD.py -r example/mol1.pdb -p example/mol2.pdb -o example/rmsd.csv
['Mol_Name', 'RMSD_Align', 'RMSD_NotAlign']

Assign bond orders from probe to reference.
mol0000 NaN 0.6661093347908585 1.7032691013398167

result save in example/rmsd.csv
(rmsd) [jcyang@x038:isoRMSD]$ cat rmsd.csv 
index,Mol_Name,RMSD_Align,RMSD_NotAlign
0,NaN,0.6661093347908585,1.7032691013398167
(rmsd) [jcyang@x038:isoRMSD]$ python isoRMSD.py -r example2/xtal-lig.pdb -p example2/test.mol2 -o example2/rmsd.txt
['Mol_Name', 'RMSD_Align', 'RMSD_NotAlign']

Assign bond orders from probe to reference.
[12:00:10] WARNING: More than one matching pattern found - picking one
mol0000 g1bl6-1_none 1.1324710650707734 1.5107986303156478

Assign bond orders from probe to reference.
[12:00:10] WARNING: More than one matching pattern found - picking one
mol0001 g1bl6-1_none 1.1327729438618293 1.6461493707768724

Assign bond orders from probe to reference.
[12:00:10] WARNING: More than one matching pattern found - picking one
mol0002 g1bl6-1_none 0.6510462210091607 1.4408565107412519

Assign bond orders from probe to reference.
[12:00:10] WARNING: More than one matching pattern found - picking one
mol0003 g1bl6-1_none 1.2676915455272233 8.289288499322483

result save in example2/rmsd.txt
(rmsd) [jcyang@x038:isoRMSD]$ cat example2/rmsd.txt 
index,Mol_Name,RMSD_Align,RMSD_NotAlign
0,g1bl6-1_none,1.1324710650707734,1.5107986303156478
1,g1bl6-1_none,1.1327729438618293,1.6461493707768724
2,g1bl6-1_none,0.6510462210091607,1.4408565107412519
3,g1bl6-1_none,1.2676915455272233,8.289288499322483
(rmsd) [jcyang@x038:isoRMSD]$ 
```
