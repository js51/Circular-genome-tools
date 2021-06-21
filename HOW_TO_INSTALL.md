# Installation Instructions
## Installing SageMath 9

The `cgt` package requires a version of `SageMath` that is based on Python 3 (as opposed to Python 2). If you're used to Python 2, the only thing that will need changing in your existing code will be print statements (they will need to be changed from `print x` to `print(x)`).

Follow the instructions located [here](https://doc.sagemath.org/html/en/installation/binary.html). Ignore the package manager solutions, the version of SageMath you will get from those is very outdated, at least on Ubuntu. 

The basic steps (for Linux) are:

- [Download binaries](https://www.sagemath.org/download-linux.html). You will most likely need the 64bit version. Choose the latest one.
- Make a directory on your machine called `Sage` or similar, where you will place the sage installation. Move the file you downloaded into this folder.
- Extract the files using `tar`. For example: `tar -xjf sage-9.3-Ubuntu_20.04-x86_64.tar.bz2.tar.bz2` After extracting, you can delete the `.bz2` file.
- `cd` into SageMath and run sage for the first time by executing the `sage` file. Wait for it to finish loading, to the point where you get the sage prompt and can type commands. Exit sage.
- Run the command `sudo ln -s /path/to/SageMath/sage /usr/local/bin/sage` You can now start SageMath from any directory by typing `sage`.

## Installing repsn

Some functions of Circular-genome-tools require the GAP package repsn to be installed within SageMath. To install repsn inside Sage,

- [Download repsn from this link](https://www.gap-system.org/Packages/repsn.html)
- Extract the files as explained above.
- move the directory `repsn-3.1.0` into the directory `SageMath/local/share/gap/pkg`

## Installing Circular Genome Tools

After completing the above steps, execute the command `sage -pip install cgt`. When you import the package, it will warn you of any missing python packages you might need to install, inclding `scipy`, `matplotlib` and `networkx`.

