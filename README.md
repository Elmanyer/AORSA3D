# AORSA 3D

The real 3D version of the AORSA code, forked from https://github.com/ORNL-Fusion/aorsa3d

This version is modified for its usage in MareNostrum 4 supercomputer.

Adapted initially by Dani Gallart and Pol Pastells.



## Compiling AORSA

To compile in MN4:

 - Load required modules: Intel compilers, NETcdf library.

```console
 foo@bar:~$ module load intel impi mkl netcdf
```

 - Compile

```console
 foo@bar:~$ make
 foo@bar:~$ ./bld_aorsaconv
```

## Running AORSA


### Execution

```console
 foo@bar:~$ sbatch /path/to/aorsa_utils/job.sh
```

### Inputs

- aorsa3d.in - contains Fortran namelists
- wout or wout.nc - contains VMEC output


### Outputs

- fpm, fpm.nc
- out15, out28, out38, transform.out, movie_wdot, out_edotj
- stdout and stderr


## Folder structure

- src: source code
- obj: compiled objects
- pgplot: pgplot related
- tables: needed tables to run the code, should be copied to execution folder (done by aorsa_utils/job.sh)
- inputs: commented aorsa3d.in file
