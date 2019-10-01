# cif2lammps
## Authors

- Ryther Anderson

## Motivation
cif2lammps is a Python 3 program used to convert crystals (developed initially for metal-organic frameworks) to Large-scale Atomic Molecular Massively Parallel Simulator (LAMMPS) format. 

## Current Status
cif2lammps can be used to convert ToBaCCo (https://github.com/tobacco-mofs/tobacco_3.0) generated cifs LAMMPS data and input according to the UFF/UFF4MOF. More force-fields and the option to use custom force-fields will be added. Keep in mind this is the very first version of the code, and I (Ryther) wrote is quite quickly, so expect to need adapt it yourself to use in the current state. 

## Usage
Generally speaking just run:
```
python main_conversion.py directory_of_cifs
```
where "directory_of_cifs" is a directory with the cifs you want to convert. By default this will convert the cifs serially and add the data and in files to a new directory called unopt_lammps_data. The options currently are:
```
--parallel=[True|False]
```
If True this will run the conversion in parallel on all the available processors.
```
--force_field=UFF
```
Currently the only option for this is UFF, meaning UFF will be used with UFF4MOF parameters (when relevant).
```
--outdir=some_path
```
This will change to output location to the specified path.
```
--charges=[True|False]
```
Charges will be considered if set to True, this means the pair_style and such will be updated to include electrostatic interactions. The default is False. Only set to True if the the region _atom_site_charge is in your cif file(s).
```
--replication=[QxRxS | min_atoms:N ]
```
The CIF cell will be replicated to the shape QxRxS or to have atleast N atoms, in the latter case the most cubic possible shape is used.
## Requirements
Just download Anaconda. This should work for Py2 or Py3.

