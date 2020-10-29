# BPTF2019 Routing Framework

## Building
### Requirements
Make sure to have the following libraries installed on your system:

- Tinyxml2
- GNU Scientific Library (GSL), including the BLAS
- OpenMP

### Structure of the project
Sourcecode is found in `src`. The subdirectory `core` contains the core router part, responsible for input/output/general datatypes. The other subdirectores are for each routing submodule ("strategy"). Used libraries are found in `lib`. When building, you can find the results in `build`.All header files should be in `include`. Stuff that does not belong here but was here is in `other`.

### Makefile Usage
You can build the project using `make` in the project root. By default, we use some debug and protection flags:

```
-g -fstack-protector-strong -fstack-clash-protection -fcf-protection -Wall -Wpedantic -Wextra -std=c++2a -D_GLIBCXX_ASSERTIONS
```

Furthermore, you can specify `SANITIZE=thread` or `SANITIZE=address` to either include the thread or address sanitizer in your build (default is address). When specifying `TYPE=DEBUG`, in addition to the flags above, the `-O0` flag is added. Without `TYPE=DEBUG`, we use `-O2 -D_FORTIFY_SOURCE=2` instead.

When specifying TYPE=RELEASE, all of these flags are omitted and `-Ofast` is added instead.

To select the used psychological model, you can use the PSYCHMOD variable. We default to user_equilibrium_2r.

To select which module you are building, you can use the STRATEGY variable. We default to SSOTD fulldisjoint. Possible strategies are:
- sstod. In this case, you can also specify SSOTD_VARIANT, which can either be onedisjoint,nodisjoint or fulldisjoint, which is our default.
- ssotdea. Leon's EA for SSOTD.
- ea. Maxi's EA.

You can also add more strategies just by creating more subfolders in src. Please remember to put your headers in include, as this is added to the include path. Also make sure, if you need more cpp files than one, to expand the Makefile accordingly. For this, you have to modify the ADDITIONALS variable. Just take a look at the example for the ea or ssotd and I'm sure you can figure it out.


## Usage
```
./router <input_file> <output_file>
```
