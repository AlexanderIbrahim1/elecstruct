# TODO

Some things to improve the code:
- [DONE] speed up the creation of the Hermitian matrices in the HF method
- [DONE] perform a check to skip the calculation of two-electron integrals that will end up as 0
- [DONE] perform cleaner separation of source and header
- [DONE] add initial Fock matrix guess as option
- [DONE] add maximum number of iterations as an option
- [DONE] add density matrix difference tolerance as an option
- [DONE] add energy difference tolerance as an option
- add more molecule examples aside from just water
- create a special directory for the main executable
- update the README to be more friendly to incoming users
- replace the Yoshimine indexing with the indexing described in the Crawford group project
  - this will let me store the two-electron integrals in a 1D array
- implement MP2
- implement CCSD
- implement CCSD(T)
- add other basis sets aside from just STO-3G
