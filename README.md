# ACEMAKER

![https://github.com/IAEA-NDS/ACEMAKER/blob/master/LICENSE](https://img.shields.io/github/license/IAEA-NDS/ACEMAKER)

A code package to produce ACE-formatted files for MCNP calculations. This updated version allows preparing fast, dosimetry, thermal and photonuclear data for Monte Carlo simulations.

The ACEMAKER driver relies on specific versions of LINEAR, RECENT, LEGEND, SPECTRA, SIGMA1, FIXUP, GROUPIE, MERGER, and DICTIN from the [PREPRO] code suite.
The modules SIXLIN, DOACE, DODOS, and DOTSL prepare the required ACE-formatted files.
Please consult the [installation instructions](#installation) for guidance on how to set up this package including its dependencies.

The user should be aware that the comment lines within the source codes should always be considered as the most recent documentation and may supersede any earlier published report.

[PREPRO]: https://github.com/iaea-nds/prepro

## Authors

* **Daniel López Aldama** - Consultant, International Atomic Energy Agency, IAEA/NDS, Vienna, Austria
* **Andrej Trkov** - International Atomic Energy Agency, IAEA/NDS, Vienna, Austria

## PREPRO Package's Author

* **Dermott E. Cullen** - 1466 Hudson Way, Livermore, CA 94550, U.S.A.

## Installation

You can download the latest released (stable) version of ACEMAKER, either as `zip` or as `tar.gz`, from the releases page:
[https://github.com/IAEA-NDS/ACEMAKER/releases/latest](https://github.com/IAEA-NDS/ACEMAKER/releases/latest)

Alternatively, you can clone this repository with git, from the command line:
```
git clone https://github.com/IAEA-NDS/acemaker.git
```
If you want to use a version that corresponds to one of the releases, make sure to check out the corresponding tagged commit.

Binaries for Windows and Ubuntu are provided in this repository. Consult the [Usage](#usage) section for learning how to run ACEMAKER.
If you use another operating system, you need to compile the binaries yourself, which is described in the next section.

### Compiling the binaries from source
ACEMAKER depends on several codes of [PREPRO].
To compile these binaries with `gfortran`, first clone the PREPRO repository and change into the `src` directory:
```
git clone https://github.com/iaea-nds/prepro.git
cd prepro/source

```
At present, specific binaries need to be compiled from two version of [PREPRO] available as
the two branches `prepro_for_acemaker_2019` and `prepro_for_acemaker_2021`, respectively.
First, compile the binaries in the branch `prepro_for_acemaker_2021`:
```
git checkout prepro_for_acemaker_2021
make install graphics=no
```
This will compile the binaries from the source files and copy them over to the `bin` directory in the root of the repository.
Now check out `prepro_for_acemaker_2019`, compile the source files and copy the `fixup` and `spectra` binaries to the `bin` directory,
replacing the one from `prepro_for_acemaker_2021`:
```
git checkout prepro_for_acemaker_2019
make
cp fixup spectra ../bin
```

To compile the ACEMAKER binaries, make the ACEMAKER repository available:
```
git clone https://github.com/IAEA-NDS/acemaker.git
```
The `src` directory contains `.bat` files for Windows and two `.sh` files for Linux to perform the compilation.
Under Linux, run the commands:
```
cd acemaker/src
./gfortgroupie.sh
./gfortcomp.sh
```
This will create the binaries in the `build` subfolder.
Copy the binaries of the PREPRO suite created above and the binaries of the `build` folder into a common directory
to have all required binaries of ACEMAKER and PREPRO at one place.

For Windows call the corresponding `.bat` files for the compilation, instead.
In this case, the binaries will be created in the `src` directory.
Copy these binaries and the ones of PREPRO created above to a common directory.

## Usage

Copy all the binaries in `bin/windows` (or `bin/linux`) if you are using the precompiled versions
to a new folder where the calculation should be performed. If you have compiled the binaries yourself
according to the instructions in the previous section, copy those files instead to the calculation directory.


## References

* ACEMAKER: A code package to produce ACE-formatted files for MCNP calculations ([IAEA-NDS-0223](https://nds.iaea.org/publications/iaea-nds/iaea-nds-223.pdf))

## Disclaimer

Neither the author nor anybody else makes any warranty, express or implied, or assumes any legal liability or responsibility for the accuracy, completeness or usefulness of any information disclosed, or represents that its use would not infringe privately owned rights.
