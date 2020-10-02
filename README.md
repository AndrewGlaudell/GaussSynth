# GaussSynth
Synthesis Package for two-qubit Gaussian Clifford + T matrices, i.e. the two-qubit Clifford + CS gate group

## Package Info

Package Version: 1.0


Written for Mathematica Version: 12.0


History:
1.0 - Initial version, completed 11/4/2019


Keywords: Quantum Compiling, Quantum Circuit Synthesis, Clifford Group, Controlled Phase Gate, Normal Forms, Exact Synthesis, Approximate Synthesis


## Description:
The GaussSynth.m package is a package for quantum compiling on two qubits using the Clifford group and the Controlled-Phase gate CS. The circuits which are exactly expressible over this gate set constitute every 4x4 unitary matrix which can be written as a matrix of Gaussian integers divided by some non-negative integer power of 2^(1/2) -- hence the package name. In this package, we supply a number of functions for performing quantum circuit synthesis on this gate set, both in the exact and approximate case. The algorithms in this package are based off of the work of Andrew Glaudell, Julien Ross, Matthew Amy, and Jake Taylor, and for details related to how these algorithms were developed, I suggest reading the articles [1-3] in the sources section below.


## Warnings:
I have used a fair amount of input checking so that functions only accept inputs of the appropriate form. This comes at the cost of some speed -- that being said, these checks cause constant overhead, and so their performance impact is worth it to prevent some erroneous calculation from being carried out. If you don't care about this input checking, one could relatively easily define their own functions from my own internal ones to slightly speed up their performance.


## Limitations:
This package is only intended for usage on two-qubit circuits. To perform circuit synthesis on larger circuits, I suggest loading this package and using these functions as subroutines.


## Algorithm Information:
Rather than describe these algorithms in detail here, I differ to the sources [1-3] listed above or the function descriptions.


## Requirements:
None


## Installation Instructions:
Install the file GaussSynth.m wherever you would like. To make things simplest, save to somewhere in your current Mathematica $Path or your working directory.


## Loading Instructions:
Make sure that GaussSynth.m is saved to somewhere on your current Mathematica $Path or in your current working directory. If it is not, either add that location to the $Path (see https://reference.wolfram.com/language/ref/$Path.html) or set your current directory to the one which includes GaussSynth.m with SetDirectory (see https://reference.wolfram.com/language/ref/SetDirectory.html).

Once this is done, use the following command to load the package:

Needs["GaussSynth`"];


## Usage Instructions:

We export a number of pre-defined variables and functions for the user. These can be checked directly in Mathematica after loading the package; for example, the command:

NormIt::usage

would list the use cases of the pre-defined NormIt function, as well as any options. Here is a list of the exported variables and functions:

- Id,X1,X2,Z1,Z2,W,H1,H2,S1,S2,CZ,CNOT12,CNOT21,EX,CS: Gates in U(4) representation

- R["P","Q"]: R[P,Q] gates in U(4) representation

- U4ToSO6: Map from U(4) to SO(6)

- IdSO6,X1SO6,X2SO6,Z1SO6,Z2SO6,\[CapitalIota]SO6,H1SO6,H2SO6,S1SO6,S2SO6,CZSO6,CNOT12SO6,CNOT21SO6,EXSO6,CSSO6: Gates in SO(6) representation

- RSO6["P","Q"]: R[P,Q] gates in SO(6) representation

- FromSequence: Read in string of operators as a U(4) matrix

- FromHexDec: ead in a Hexidecimal representation of a Clifford CS operator as a U(4) matrix

- CliffordQ: Check if an operator is a Clifford

- CliffordSynth: Find a Clifford circuit (in HexDec and String forms) for a Clifford operator

- RightCliffordSimilar, LeftCliffordSimilar: Check if two operators are right/left Clifford-similar

- GaussianQ: Check if an operator is a Clifford + CS operator

- LDE: Find the Least Denominator Exponent in either the U(4) or SO(6) representations

- OptimalCSCount: Find the optimal number of CS gates required to synthesize an operator in any representation

- SyllableList: List of canonical R[P,Q] Syllables for a normal form

- SyllableListAsymmetric: List of constructors for Clifford + CS for an alternative normal form using the usual generators

- NormIt: Find the normal form for a Clifford + CS operator

- FrobeniusDistance: Calculate the Frobenius distance between two matrices

- CandidateFinder: Find a two integer-component vector (a,b) and integer k such that 1/sqrt(2)^k (a,b) is close to some complex phase

- PauliRotation: Find a Clifford + CS operator which approximates a pauli-rotation on two qubits

- PauliRotationSequence: Find the corresponding normal-form sequence of a Clifford + CS operator which approximates a pauli-rotation

- PauliDecomposition: Find the Pauli rotation decomposition of a two-qubit unitary operator

- ApproximateOp: Find a Clifford + CS operator which approximates some operator on two qubits

- ApproximateSequence: Find the corresponding normal-form sequence of a Clifford + CS operator which approximates some operator on two qubits

- RandomCliffCS: Pseudorandomly sample from the uniform distribution of all Clifford + CS operators whose optimal CS-count is at-most some integer n.

## Copyright Info:
Copyright © 2019 Andrew Glaudell

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


## Sources:
[1] Amy, Matthew, Andrew N. Glaudell, and Neil J. Ross. "Number-theoretic characterizations of some restricted clifford+ t circuits." Quantum 4 (2020): 252.

[2] Glaudell, Andrew N., Neil J. Ross, and Jacob M. Taylor. "Optimal Two-Qubit Circuits for Universal Fault-Tolerant Quantum Computation." arXiv preprint arXiv:2001.05997 (2020). (under review)

[3] Glaudell, Andrew Noble. Quantum Compiling Methods for Fault-Tolerant Gate Sets of Dimension Greater than Two. Diss. 2019.
