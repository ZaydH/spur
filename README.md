# SPUR - Satisfying Perfectly Uniform Random Sampler

**SPUR** - Published at the [SAT 2018 Conference](http://sat2018.azurewebsites.net/)  
**Winner of Best Student Paper**

**Title**: Fast Sampling of Perfectly Uniform Satisfying Assignments  
**Authors**: [Dimitris Achlioptas](https://users.soe.ucsc.edu/~optas/), [Zayd Hammoudeh](https://users.soe.ucsc.edu/~zayd/), and Panos Theodoropoulos

SPUR is a fast uniform SAT assignment sampler.  It was presented at SAT-2018. 

# Building and Running the Binary

For user convenience, we have included in the root of this repository a build script, `build.sh`.  Running this script will build SPUR from source.

After building the binary files, SPUR can be run by calling from the command line:

`build\Release\spur -s <#Samples> -cnf <FormulaFile>` 

where `<#Samples>` is the number of samples to collect (>=1) and `<FormulaFile>` is a path to Boolean CNF formula in DIMACS file format.

The full set of command line options are:
* `-cnf [cnf_file]` Path to the CNF file
* `-s [s]` 	 Number of models "s" to sample uniformly at random
* `-tp`    	 Forces two-pass sampling. (Only applicable when s=1)
* `-out [out_file]`  Path to write the specified samples
* `-no-sample-write` Disable writing the final samples to a file.
* `-count-only`	 Perform only model counting.  Disable sampling.
* `-q`     	 Quiet mode
* `-v`     	 Verbose and trace mode
* `-d`     	 Debug mode
* `-t [s]` 	 Set time bound to s seconds
* `-cs [n]`	 Set max cache size to n MB

# Dependencies

* **Operating System**: This software is only tested on Ubuntu Linux and Mac OS High Sierra. Additional operating systems, but functionality is not guaranteed nor implied.
* **C++ Version**: C++11
* **Libraries**:
  - [GNU Multiple Precision Library](https://gmplib.org/)

# Output Format

As explained in our paper, the generated witnesses are represented as tuples.  The first item in the tuple is an integer between 1 and the number of samples (inclusive); this number represents the number of witnesses entailed by this tuple.  

The second item in the tuple is a variable assignment.  Variables are ordered from 1 to *N*, where *N* is the total number of variables.  Each variable is assigned to either "0", "1", or "\*".  Observe that the "\*" corresponds to an unconstrained variable that can be either "0" or "1".

**Example**: The tuple "3,10\*0\*" means that three witnesses are entailed by the assignment where variable one is 0b1, variables two and four are 0b0, while variables three and five can take either value.

# Python Libraries Available

Python libraries for running the SPUR binary are available upon request.  Email zayd.hammoudeh@gmail.com to request access.  Please include in the email request your name, email address, and institution.

# Special Thanks

The SPUR sampler is built on top of [sharpSAT](https://github.com/marcthurley/sharpSAT), which was developed by Marc Thurley.

