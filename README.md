# apt-vs-dift

Extension of the project PRISM-games [https://github.com/prismmodelchecker/prism-games] with the algorithms as described in CDC submission "Guaranteed Trade-Offs in Dynamic Information Flow Tracking Games".

## License:

For our modifications of the code, we give an MIT license. For the code that was originally part of PRISM-games, we refer to their license file in prism-games-3.0.beta-src/COPYING.txt and the "Licensing" section of their README in prism-games-3.0.beta-src/README.md.

## Organization:

- The folder final-plots contains the plots as they are used in Figure 4 of the paper for all case studies.
- The folders examples and props contain the case study files.
- The folder results contains all the results, including intermediate files we have to generate and final plots.
- The folder src contains our code, specifically two versions of PRISM-games: 
  1) The code for simulating models and estimating the probabilities in prism-fruit (extension of the code of this paper https://link.springer.com/chapter/10.1007%2F978-3-030-25540-4_29)
  2) The code for computing the Pareto frontier (which is the border of the set of achievable vectors) of a turn-based SG in prism-games.

To run the code and recreate the results, first install both prism-fruit as well as prism-games (the process is the same for both) and then run results.sh.

## Installation

### Dependencies
- Java 8 or more recent
- Python 3.7 or more recent
- C and C++ compiler
- PRISM-Games (https://prismmodelchecker.org/games/installation.php)
  - PPL

**We assume that you are using a Linux distribution.** For Windows and Mac OS users the guides we refer to also provide instructions, but you may encounter difficulties putting everything together.

### Setting up PRISM-games

Our implementation is an extension of PRISM-games.
Thus, we require PPL and PRISM to be installed.
Follow the guidelines given at:<br/>
https://prismmodelchecker.org/games/installation.php<br/>
Now the following command should execute without errors from the prism-games-3.0.beta-src/prism folder:
`./bin/prism`
At this point you should be able to run the scripts to generate the results.
