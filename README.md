# Column Parity Mixer trail search
This repository contains software to search truncated linear and differential trails for permutations/ciphers based on a column parity mixer, as is described in this [this academic research paper](https://tosc.iacr.org/index.php/ToSC/article/view/847/799). It is currently instantiated with [Mixifer](https://github.com/Ko-/mixifer).

## Compiling
Compiling requires a reasonably up-to-date C++ compiler that can deal with C++14. Just type 'make' to compile.

## What's in there
`main.cpp` contains the main function, which is currently empty. Uncomment one of the functions to test them. `distribution.cpp` and `bruteforce.cpp` are standalone tools that have been used at some point.

The software is loosely based on an earlier version of code that is now in the [KeccakTools repository](https://github.com/gvanas/KeccakTools). However, a lot has been modified to be able to deal with truncated trails and strong alignment.

Note that the code is not very polished/clean. It also still contains some old stuff that is not even touched/working anymore. Sorry about that.

The code contains several parts that are specific to the case of `m=4`, i.e. a state with 4 rows. So only changing `Permutation.cpp` is not going to be enough to try this on a different cipher. Although this was the goal originally, optimization were possible for this particular case that were deemed necessary to be able to cover a larger search space.
