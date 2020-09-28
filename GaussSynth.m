(* ::Package:: *)

(* ::Title:: *)
(*GaussSynth.m -- The Two-Qubit Clifford + CS Circuit Synthesis Package*)


(* ::Subtitle:: *)
(*Written and Maintained by Andrew Glaudell*)


(* ::Subtitle:: *)
(*The GaussSynth.m package is a package for quantum compiling on two qubits using the Clifford group and  the Controlled-Phase gate CS. The circuits which are exactly expressible over this gate set constitute every 4x4 unitary matrix which can be written as a matrix of Gaussian integers divided by some non-negative integer power of 2^(1/2) -- hence the package name. In this package, we supply a number of functions for performing quantum circuit synthesis on this gate set, both in the exact and approximate case. The algorithms in this package are based off of the work of Andrew Glaudell, Julien Ross, Matthew Amy, and Jake Taylor, and for details related to how these algorithms were developed, I suggest reading the articles [1-3] in the sources section below.*)


(* ::Section:: *)
(*Package Details*)


(* ::Text:: *)
(*Copyright \[Copyright] 2019 Andrew Glaudell*)
(**)
(*Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:*)
(**)
(*The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.*)
(**)
(*THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.*)
(**)
(*Package Version: 1.0*)
(**)
(*Written for Mathematica Version: 12.0*)
(**)
(*History:*)
(*1.0 - Initial version, completed 11/4/2019*)
(**)
(*Keywords: Quantum Compiling, Quantum Circuit Synthesis, Clifford Group, Controlled Phase Gate, Normal Forms, Exact Synthesis, Approximate Synthesis*)
(**)
(*Sources:*)
(*[1] Matt Amy, Andrew Glaudell, and Neil J. Ross. Number-theoretic characterizations of some restricted clifford+t circuits. Upcoming publication, preprint available from arXiv:1908.06076, 2019.*)
(*[2] Andrew Glaudell, Neil J. Ross, and Jacob M. Taylor. Optimal two-qubit circuits for universal fault-tolerant quantum computation. Upcoming publication, 2019.*)
(*[3] Andrew Glaudell. Inexact synthesis of pauli-rotations with the fault-tolerant clifford + controlled-phase gate set. Upcoming publication, in preparation, 2019.*)
(**)
(*Warnings:*)
(*I have used a fair amount of input checking so that functions only accept inputs of the appropriate form. This comes at the cost of some speed -- that being said, these checks cause constant overhead, and so their performance impact is worth it to prevent some erroneous calculation from being carried out. If you don't care about this input checking, one could relatively easily define their own functions from my own internal ones to slightly speed up their performance.*)
(**)
(*Limitations:*)
(*This package is only intended for usage on two-qubit circuits. To perform circuit synthesis on larger circuits, I suggest loading this package and using these functions as subroutines.*)
(**)
(*Discussion:*)
(*Rather than describe these algorithms in detail here, I differ to the sources [1-3] listed above or the function descriptions.*)
(**)
(*Requirements:*)
(*None*)


BeginPackage["GaussSynth`"];


(* ::Section:: *)
(*\[AliasDelimiter]\[AliasDelimiter]\[AliasDelimiter]Function Usage*)


GaussSynth::usage = "GaussSynth is a package for quantum compiling for two-qubit circuits using the Clifford + CS gate set";

Id::usage = "Id is the 4x4 Identity Matrix.";
X1::usage = "X1 is the unitary representation of the X\[CircleTimes]I gate.";
X2::usage = "X2 is the unitary representation of the I\[CircleTimes]X gate.";
Z1::usage = "Z1 is the unitary representation of the Z\[CircleTimes]I gate.";
Z2::usage = "Z2 is the unitary representation of the I\[CircleTimes]Z gate.";
W::usage = "W is the unitary representation of the primite 8th root of unity \[Omega].";
H1::usage = "H1 is the unitary representation of the H\[CircleTimes]I gate.";
H2::usage = "H2 is the unitary representation of the I\[CircleTimes]H gate.";
S1::usage = "S1 is the unitary representation of the S\[CircleTimes]I gate.";
S2::usage = "S2 is the unitary representation of the I\[CircleTimes]S gate.";
CZ::usage = "CZ is the unitary representation of the CZ gate.";
CNOT12::usage = "CNOT12 is the unitary representation of the CNOT gate with control qubit 1 and target qubit 2.";
CNOT21::usage = "CNOT21 is the unitary representation of the CNOT gate with control qubit 2 and target qubit 1.";
EX::usage = "EX is the unitary representation of the SWAP (Exchange) gate.";
CS::usage = "CS is the unitary representation of the CS gate.";

U4ToSO6::usage = "U4ToSO6[U] maps the 4x4 unitary U to the equivalent SO(6) representation for the operator U' = \!\(\*SuperscriptBox[\(e\), \(i\[Phi]\)]\) U \
	for some phase \[CurlyPhi] and U' an element of SU(4).";

IdSO6::usage = "Id is the 6x6 Identity Matrix.";
X1SO6::usage = "X1SO6 is the SO(6) representation of the X\[CircleTimes]I gate.";
X2SO6::usage = "X2SO6 is the SO(6) representation of the I\[CircleTimes]X gate.";
Z1SO6::usage = "Z1SO6 is the SO(6) representation of the Z\[CircleTimes]I gate.";
Z2SO6::usage = "Z2SO6 is the SO(6) representation of the I\[CircleTimes]Z gate.";
\[CapitalIota]SO6::usage = "\[CapitalIota]SO6 is the SO(6) representation of the complex phase I. As SU(4) is a double cover of SO(6), \
	we have \!\(\*SuperscriptBox[\(\[CapitalIota]SO6\), \(2\)]\) = IdSO6.";
H1SO6::usage = "H1SO6 is the SO(6) representation of the H\[CircleTimes]I gate.";
H2SO6::usage = "H2SO6 is the SO(6) representation of the I\[CircleTimes]H gate.";
S1SO6::usage = "S1SO6 is the SO(6) representation of the S\[CircleTimes]I gate.";
S2SO6::usage = "S2SO6 is the SO(6) representation of the I\[CircleTimes]S gate.";
CZSO6::usage = "CZSO6 is the SO(6) representation of the CZ gate. Note that this gate is not special unitary, and so we multiply \
	by a primitive 8th root of unity \!\(\*SuperscriptBox[\(\[Omega]\), \(\[Dagger]\)]\) before performing the transformation.";
CNOT12SO6::usage = "CNOT12SO6 is the SO(6) representation of the CNOT gate with control qubit 1 and target qubit 2. \
	Note that this gate is not special unitary, and so we multiply by a primitive 8th root of unity \!\(\*SuperscriptBox[\(\[Omega]\), \(\[Dagger]\)]\) \
	before performing the transformation.";
CNOT21SO6::usage = "CNOT21SO6 is the SO(6) representation of the CNOT gate with control qubit 2 and target qubit 1. \
	Note that this gate is not special unitary, and so we multiply by aa primitive 8th root of unity \!\(\*SuperscriptBox[\(\[Omega]\), \(\[Dagger]\)]\) \
	before performing the transformation.";
EXSO6::usage = "EXSO6 is the SO(6) representation of the SWAP (Exchange) gate. Note that this gate is not special unitary, \
	and so we multiply by a primitive 8th root of unity \!\(\*SuperscriptBox[\(\[Omega]\), \(\[Dagger]\)]\) before performing the transformation.";
CSSO6::usage = "CSSO6 is the SO(6) representation of the CS gate. Note that this gate is not special unitary, and so we multiply by \
	a primitive 16th root of unity \!\(\*SuperscriptBox[\(\[Omega]\), \(\(-1\)/2\)]\) before performing the transformation.";

FromSequence::usage = "FromSequnce[str] reads in the string str and interprets that string as a Clifford + CS operator in the U(4) representation.";

FromHexDec::usage = "FromHexDec[str] attempts to read in a string of a signed hexadecimal number as a Clifford + CS operator in the U(4) representation. \
	The sign indicates whether the operator corresponds to using symmetric or asymmetric syllables. The list of syllables takes hexadecimal values 1-f for \
	the fifteen unique syllables. After the syllables comes the marker 00000 which is followed by a hexidcimal integer which take the decimal values 1-92160.";

CliffordQ::usage = "CliffordQ[U] returns True if U is a Clifford and False otherwise.";
CliffordSynth::usage = "CliffordSynth[U] gives the index number (in Hexidecimal) and string of Clifford and Pauli operators \
	which constitute the Clifford which can be input as a string, an element of U(4), an element of SO(6), or a hexidcimal representation";
RightCliffordSimilar::usage = "CliffordSimilarRight[U,V] Returns True if there is a Clifford C such that U.C = V and False otherwise.";
LeftCliffordSimilar::usage = "CliffordSimilarLeft[U,V] Returns True if there is a Clifford C such that C.U = V and False otherwise.";

GaussianQ::usage = "GaussianQ[U] is a Boolean function which checks if U corresponds to a Clifford + CS circuit.";

SyllableList::usage = "A list of 15 syllables of CS-count one which are not right-Clifford similar. Each syllable is a \
	Clifford C which conjugates CS as C.CS.\!\(\*SuperscriptBox[\(C\), \(\[Dagger]\)]\). For each syllable, we supply the operator's \
	4x4 Unitary representation, the operator's 6x6 SO(6) representation, and its name according to the generators in string form.";
SyllableListAsymmetric::usage = "An alternative list of 15 syllables of CS-count one which are not right-Clifford similar. \
	Each syllable is equivalent to C.CS for C a Clifford. For each syllable, we supply the operator's 4x4 Unitary representation, \
	the operator's 6x6 SO(6) representation, and its name according to the generators in string form.";

NormIt::usage = "NormIt[U,options] takes as input a Gaussian Clifford + T operator U and outputs its equivalent normal form. \
	This normal form is output as a string of generators using the standard syllable list unless specified otherwise in the options. \
	The options for the \"OutputType\" are \"String\" or \"HexDec\", the options for the \"SyllableType\" are \"Normal\" or \"Asymmetric\", \
	and the options for \"UpToPhase\" are the booleans True and False. When either reading or outputting a string of generators, we use (only) the \
	following characters: \"W\", \"S1\", \"S2\", \"H1\", \"H2\", \"CZ\", \"EX\", \"X1\", \"X2\", \"Z1\", \"Z2\", \"CS\" where \"W\" is the phase \[Omega], \
	\"S1\" is the gate S1, and so on.";

FrobeniusDistance::usage = "FrobeniusDistance[U,V] computes the distance between U and V via the Frobenius Matrix Norm.";

PauliRotation::usage = "PauliRotation[\[CurlyPhi],\[Epsilon],Pauli] finds a unitary Gaussian Clifford + T operator which is within Frobenius Distance \[Epsilon] \
	of the Pauli rotation \!\(\*SuperscriptBox[\(e\), \(\(-i\[CurlyPhi]\)/2\\\ P\)]\) for the pauli matrix P which is specified by the string Pauli. \
	The Pauli can be one of the fifteen strings \"XI\", \"YI\", \"ZI\", \"IX\", \"IY\", \"IZ\", \"XX\", \"YX\", \"ZX\", \"XY\", \"YY\", \"ZY\", \
	\"XZ\", \"YZ\", or \"ZZ\".";
PauliRotationSequence::usage = "PauliRotationSequence[\[CurlyPhi],\[Epsilon],Pauli,options] finds a unitary Gaussian Clifford + T operator which is within Frobenius Distance \[Epsilon] \
	of the Pauli rotation \!\(\*SuperscriptBox[\(e\), \(\(-i\[CurlyPhi]\)/2\\\ P\)]\) for the pauli matrix P which is specified by the string Pauli. \
	The Pauli can be one of the fifteen strings \"XI\", \"YI\", \"ZI\", \"IX\", \"IY\", \"IZ\", \"XX\", \"YX\", \"ZX\", \"XY\", \"YY\", \"ZY\", \
	\"XZ\", \"YZ\", or \"ZZ\". It then outputs a normalized sequence of Clifford + CS operators via NormIt. The options for the \"OutputType\" are \
	\"String\" or \"HexDec\" and the options for the \"SyllableType\" are \"Normal\" or \"Asymmetric\".";
PauliDecomposition::usage = "PauliDecomp[U] finds a list of 15 angle parameters \!\(\*SubscriptBox[\(\[CurlyPhi]\), \(j\)]\) which constitutee a \
	decomposition of the form \!\(\*SubscriptBox[\(\[Product]\), \(1 \[LessEqual] j \[LessEqual] 15\)]\)\!\(\*SuperscriptBox[\(e\), \
	\(\(-\*SubscriptBox[\(i\[CurlyPhi]\), \(j\)]\)/2\\\ \*SubscriptBox[\(P\), \(j\)]\)]\) for the U(4) (up to a phase) or SO(6) operator U. The sequence of operators for \
	the Paulis \!\(\*SubscriptBox[\(P\), \(J\)]\) is (ZI,XI,ZI,IZ,IX,IZ,XX,YY,ZZ,ZI,XI,ZI,IZ,IX,IZ).";
ApproximateOp::usage = "Approximate[U,\[Epsilon]] finds an approximation within Frobenius Distance \[Epsilon] of the unitary or SO(6) operator U in the Clifford + CS gate set. \
	If U is an element of U(4), the result is a U(4) representation of a Clifford + CS circuit (up to a phase), and if U is an element of SO(6) the result is an \
	SO(6)  representation of a Clifford + CS circuit. Note that the Frobenius distance between U and its approximation is always calculated in the U(4) \
	representation.";
ApproximateSequence::usage = "ApproximateSequence[U,\[Epsilon],options] finds a normalized sequence of Clifford + CS operators which is within Frobenius Distance \[Epsilon] \
	(in the Unitary representation and up to an irrelevant phase) for the input U. U may be either an element of U(4) or SO(6). The output is given as a \
	string unless otherwise specified in the options. The options for the \"OutputType\" are \"String\" or \"HexDec\", the options for the \
	\"SyllableType\" are \"Normal\" or \"Asymmetric\", and the options for \"IfGaussianDoExact\" are the booleans True and False.";


(* ::Section::Closed:: *)
(*Possible Errors*)


General::invldopt = "Option `2` for function `1` received invalid value `3`";
U4ToSO6::notunitary = "The argument must be a 4x4 unitary matrix.";
FromSequence::notstring = "You have not entered a string.";
FromList::notagate = "The string `1` is not one of \"W\", \"S1\", \"S2\", \"H1\", \"H2\", \"CZ\", \"CS\", \"EX\", \"X1\", \"X2\", \"Z1\", or \"Z2\". 
	You may have forgotten a space between gate names or used a name for a gate which is not recognized.";
FromHexDec::invalidnumber = "The string `1` is not a valid hexadecimal representation. Make sure your string has the string \"00000\" 
	seperating the syllables from the Clifford. Otherwise, ensure the Clifford has an index from 1 to 92160 in decimal representation 
	before being written in hexadecimal representation, and that your syllables only take values of 1,2,3,4,5,6,7,8,9,a,b,c,d,e, or f.";
CliffordQ::notacircuit = "You have not entered a string of operators, a valid hexadecimal, an element of SO(6), or an element of U(4).";
CliffordSynth::notaclifford = "Your input is not a Clifford operator in string form, a valid hexadecimal, U(4) representation, or SO(6) representation.";
GaussianQ::notacircuit = "You have not entered a string of operators, a valid hexadecimal, an element of SO(6), or an element of U(4).";
NormIt::badopt = "The option `1` is not valid for `2`.";
NormIt::notacircuit = "You have not entered a string of operators, a valid hexadecimal, an element of SO(6), or an element of U(4).";
FrobeniusDistance::notequidimensionalmatrices = "Your inputs are not two matrices of equal size.";
CandidateFinder::notreals = "Your inputs are not two real numbers";
PauliRotation::notreals = "Your input does not include two real numbers";
PauliRotation::invldstring = "Your input does not include one string from the set of \"XI\", \"YI\", \"ZI\", \"IX\", \"IY\", \"IZ\", \"XX\", 
	\"YX\", \"ZX\", \"XY\", \"YY\", \"ZY\", \"XZ\", \"YZ\", or \"ZZ\".";
PauliDecomposition::notanoperator = "Your input is neither an element of U(4) or SO(6) and so it cannot be decomposed";
ApproximateOp::notreal = "Your error tolerance is not a real number.";
ApproximateOp::notanoperator = "Your input is neither an element of U(4) or SO(6) and so it cannot be approximated";


Begin["`Private`"];


(* ::Section:: *)
(*Function Definitions*)


(* ::Subsection:: *)
(*Functions for option checking*)


(* ::Text:: *)
(*These functions will be used to check options for functions which accept them. For each such function, we must define a test[f,op] function for a particular option type op of function f. Credit for this code snippet goes to Mr. Wizard in the Stack Overflow post https://mathematica.stackexchange.com/questions/116623.*)


optsMsg[f_][op_, val_] :=
  test[f, op][val] || Message[General::invldopt, f, op, val];

Attributes[optsCheck] = {HoldFirst};

optsCheck @ head_[___, opts : OptionsPattern[]] :=
  And @@ optsMsg[head] @@@ FilterRules[{opts}, Options @ head];


(* ::Subsection::Closed:: *)
(*Constants and Single-Qubit operators*)


(* ::Text:: *)
(*For internal use only.*)


\[Omega] = (1+I)/Sqrt[2];
\[Zeta] = Exp[I*Pi/8];
s = DiagonalMatrix[{1,I}];
h = 1/Sqrt[2]*{{1,1},{1,-1}};
x = PauliMatrix[1];
z = PauliMatrix[3];


(* ::Subsection::Closed:: *)
(*Unitary Representations of Two-qubit Clifford + CS operators*)


(* ::Text:: *)
(*These operators are exported to the user as 4x4 matrices in the standard Mathematica format.*)


Id = IdentityMatrix[4];
W = \[Omega]*Id;
S1 = KroneckerProduct[s,IdentityMatrix[2]];
S2 = KroneckerProduct[IdentityMatrix[2],s];
H1 = KroneckerProduct[h,IdentityMatrix[2]];
H2 = KroneckerProduct[IdentityMatrix[2],h];
CZ = DiagonalMatrix[{1,1,1,-1}];
CNOT12 = {
	{1,0,0,0},
	{0,1,0,0},
	{0,0,0,1},
	{0,0,1,0}
};
CNOT21 = {
	{1,0,0,0},
	{0,0,0,1},
	{0,0,1,0},
	{0,1,0,0}
};
EX = {
	{1,0,0,0},
	{0,0,1,0},
	{0,1,0,0},
	{0,0,0,1}
};
CS = DiagonalMatrix[{1,1,1,I}];
X1 = KroneckerProduct[x,IdentityMatrix[2]];
X2 = KroneckerProduct[IdentityMatrix[2],x];
Z1 = KroneckerProduct[z,IdentityMatrix[2]];
Z2 = KroneckerProduct[IdentityMatrix[2],z];


(* ::Subsection::Closed:: *)
(*Checks For U(4) and SO(6)*)


(* ::Text:: *)
(*These Boolean functions determine whether an operator is an element of U(4) or SO(6), respectively.*)


U4Q[U_]:= UnitaryMatrixQ[U] && (Dimensions[U] == {4,4});
SO6Q[O_] := OrthogonalMatrixQ[O] && (Dimensions[O] == {6,6}) && (Det[O] == 1);


(* ::Subsection::Closed:: *)
(*The SU(4)\[TildeFullEqual]SO(6) Isomorphism*)


(* ::Text:: *)
(*These definitions and functions allow one to compute the SO(6) representation of an element of U(4) (up to a phase).*)


(* ::Subsubsection:: *)
(*Rules for Inner Products of Wedge Products*)


(* ::Text:: *)
(*Defined using the unassigned Mathematica symbols of \[Wedge], \[LeftAngleBracket], and \[RightAngleBracket].*)


\[LeftAngleBracket]a_,O_,b_*c_\[RightAngleBracket]:=b*\[LeftAngleBracket]a,O,c\[RightAngleBracket];
\[LeftAngleBracket]a_*b_,O_,c_\[RightAngleBracket] := Conjugate[a]*\[LeftAngleBracket]b,O,c\[RightAngleBracket];
\[LeftAngleBracket]a_+b_,O_,c_\[RightAngleBracket]:=\[LeftAngleBracket]a,O,c\[RightAngleBracket]+\[LeftAngleBracket]b,O,c\[RightAngleBracket];
\[LeftAngleBracket]a_,O_,b_+c_\[RightAngleBracket]:=\[LeftAngleBracket]a,O,b\[RightAngleBracket]+\[LeftAngleBracket]a,O,c\[RightAngleBracket];
\[LeftAngleBracket]x_\[Wedge]y_,O_,u_\[Wedge]v_\[RightAngleBracket] := (Conjugate[x].O.u)*(Conjugate[y].O.v) - (Conjugate[x].O.v)*(Conjugate[y].O.u)


(* ::Subsubsection:: *)
(*Orthonormal Basis for Subscript[\[DoubleStruckCapitalC], 6]*)


(* ::Text:: *)
(*This basis is such that computing the above inner products for an element of U(4,\[DoubleStruckCapitalC]) will produce an element of SO(6,\[DoubleStruckCapitalR]). Moreover, the representations for Clifford + CS operators are easy to work with in this basis.*)


Basis6 ={
	(I/Sqrt[2])*(UnitVector[4,1]\[Wedge]UnitVector[4,2] - UnitVector[4,3]\[Wedge]UnitVector[4,4]),
	(1/Sqrt[2])*(UnitVector[4,1]\[Wedge]UnitVector[4,2] + UnitVector[4,3]\[Wedge]UnitVector[4,4]),
	(I/Sqrt[2])*(UnitVector[4,2]\[Wedge]UnitVector[4,3] - UnitVector[4,1]\[Wedge]UnitVector[4,4]),
	(1/Sqrt[2])*(UnitVector[4,2]\[Wedge]UnitVector[4,4] + UnitVector[4,3]\[Wedge]UnitVector[4,1]),
	(I/Sqrt[2])*(UnitVector[4,2]\[Wedge]UnitVector[4,4] - UnitVector[4,3]\[Wedge]UnitVector[4,1]),
	(1/Sqrt[2])*(UnitVector[4,2]\[Wedge]UnitVector[4,3] + UnitVector[4,1]\[Wedge]UnitVector[4,4])
};


(* ::Subsubsection:: *)
(*Calculating the SO(6) Representation for an element of SU(4) (and from U(4) up to a phase)*)


(* ::Text:: *)
(*Functions for mapping elements of SU(4) to SO(6) and U(4) to SO(6) (by converting that element of U(4) to an element of SU(4)).*)


SU4ToSO6[U_] := Table[Simplify[\[LeftAngleBracket]Basis6[[i]],U,Basis6[[j]]\[RightAngleBracket]],{i,1,6},{j,1,6}];

U4ToSO6[U_/;U4Q[U]] := SU4ToSO6[1/(Det[U])^(1/4)*U];
U4ToSO6[U_] := (
	Message[U4ToSO6::notunitary];
	$Failed
);


(* ::Subsection::Closed:: *)
(*SO(6) Representations of Two-qubit Clifford + CS operators*)


(* ::Text:: *)
(*Calculated using our transformations. These operators are exported to the user as 6x6 matrices in the standard Mathematica format. Note that we have to multiply by overall phases to ensure that the transformation uses an element of SU(4).*)


IdSO6 = SU4ToSO6[Id];
\[CapitalIota]SO6 = SU4ToSO6[W.W];
H1SO6 = SU4ToSO6[H1];
H2SO6 = SU4ToSO6[H2];
S1SO6 = SU4ToSO6[\[Omega]^(-1)*S1];
S2SO6 = SU4ToSO6[\[Omega]^(-1)*S2];
CZSO6 = SU4ToSO6[\[Omega]^(-1)*CZ];
CNOT12SO6 = SU4ToSO6[\[Omega]^(-1)*CNOT12];
CNOT21SO6 = SU4ToSO6[\[Omega]^(-1)*CNOT21];
EXSO6 = SU4ToSO6[\[Omega]^(-3)*EX];
CSSO6 = SU4ToSO6[\[Zeta]^(-1)*CS];
X1SO6 = SU4ToSO6[X1];
X2SO6 = SU4ToSO6[X2];
Z1SO6 = SU4ToSO6[Z1];
Z2SO6 = SU4ToSO6[Z2];


(* ::Subsection::Closed:: *)
(*Custom Representations of SO(6) Clifford + CS operators*)


(* ::Text:: *)
(*Our synthesis algorithms will use a custom data type for the SO(6) representation of a Clifford + CS operator. The basic data structure of this special representation is as follows:*)
(**)
(*{k,M} := 2^(-k/2) \[CenterDot] M*)
(**)
(*This allows easy tracking of the lde. They are packed in SparseArrays to help make things even a little faster, as every Clifford is just a permutation matrix. These representations are for internal use only.*)


(* ::Subsubsection:: *)
(*Switching between the standard representation for a 6x6 matrix and a representation specifically for Clifford + CS operators. *)


(* ::Text:: *)
(*Functions for converting to and from the special representation.*)


SO6ToSpecialRep[M_] := Module[{LDE,IntegerMat},
	LDE = FullSimplify[Max[Map[Simplify[Log[Sqrt[2],Denominator[#]]]&,Flatten[Simplify[M]]]]];
	IntegerMat = FullSimplify[(Sqrt[2]^LDE)*M];
	{LDE,SparseArray[IntegerMat]}
];
SpecialRepToSO6[{k_,M_}] := Simplify[1/Sqrt[2]^k*Normal[M]];


(* ::Subsubsection:: *)
(*Special Representations for SO(6) Clifford + CS operators*)


IdSp = SO6ToSpecialRep[IdSO6];
\[CapitalIota]Sp = SO6ToSpecialRep[\[CapitalIota]SO6];
H1Sp = SO6ToSpecialRep[H1SO6];
H2Sp = SO6ToSpecialRep[H2SO6];
S1Sp = SO6ToSpecialRep[S1SO6];
S2Sp = SO6ToSpecialRep[S2SO6];
CZSp = SO6ToSpecialRep[CZSO6];
CNOT12Sp = SO6ToSpecialRep[CNOT12SO6];
CNOT21Sp = SO6ToSpecialRep[CNOT21SO6];
EXSp = SO6ToSpecialRep[EXSO6];
CSSp = SO6ToSpecialRep[CSSO6];
X1Sp = SO6ToSpecialRep[X1SO6];
X2Sp = SO6ToSpecialRep[X2SO6];
Z1Sp = SO6ToSpecialRep[Z1SO6];
Z2Sp = SO6ToSpecialRep[Z2SO6];


(* ::Subsubsection:: *)
(*Basic Matrix Operations for the Special Representation*)


(* ::Text:: *)
(*These internal functions are used to reduce the denominator exponent to the lde, multiply operators in the special representation, and invert the special representation.*)


KReduceOnce[{k_,a_}]:=If[AllTrue[a,EvenQ,2] && k>1,{k - 2,a/2},{k,a}];
KReduce[o_]:=FixedPoint[KReduceOnce,o];

Dot2Sp[{k1_,a1_},{k2_,a2_}]:= KReduce[{k1 + k2,a1.a2}];
DotSp[x__]:=Fold[Dot2Sp,IdSp,{x}];

InvSp[{k_,a_}]:={k,Transpose[a]};


(* ::Subsubsection:: *)
(*Functions for Residues Modulo 2 and Finding Paired Matrix Rows*)


(* ::Text:: *)
(*These functions are used for finding paired rows in operators in the special representation.*)


PatternMats[{k_,a_}] := Mod[a,2];

MatrixRowPairs[list_] := GatherBy[
	Range@Length[Normal[list]],
	Normal[list][[#]]&
];

RowPairs[x_] := Sort@MatrixRowPairs@PatternMats@x


(* ::Subsection::Closed:: *)
(*String Reading*)


(* ::Text:: *)
(*Functions for reading in a string of operators. The string is always read in as an element of U(4).*)


FromSequence[str_String] := Module[{strlist},
	strlist = StringSplit[str];
	FromList[strlist]
];
FromSequence[x_] := (
	Message[FromSequence::notstring];
	$Failed
);

FromList[{}] = Id;
FromList[{"W"}] = W;
FromList[{"S1"}] = S1;
FromList[{"S2"}] = S2;
FromList[{"H1"}] = H1;
FromList[{"H2"}] = H2;
FromList[{"CZ"}] = CZ;
FromList[{"CS"}] = CS;
FromList[{"EX"}] = EX;
FromList[{"X1"}] = X1;
FromList[{"X2"}] = X2;
FromList[{"Z1"}] = Z1;
FromList[{"Z2"}] = Z2;
FromList[{str_}]:= (
	Message[FromList::notagate,str];
	$Failed
);
FromList[{h_,t__}] := FromList[{h}] . FromList[{t}];


(* ::Subsection::Closed:: *)
(*Hexadecimal Representations*)


(* ::Text:: *)
(*We shall use strings of signed hexadecimal integers to represent Clifford + CS operators. This form is much more compact than, for example, the string form of an operator. The string must be of the following form:*)
(**)
(*\"( | -)(1-9 | a-f)* 00000 C\"*)
(**)
(*where C is the hexadecimal representation of an integer from 1 to 92160 and the first option is simply an optional "-" character.*)


ValidHexDecQ[num_String] := Module[{seperator,typenum,syllabletype,syllablescliffordok,split,syllables,clifford,syllablescharacterlist,syllablesok,cliffordcharacterlist,cliffordpossible,cliffordok},
	seperator = StringCount[num,"00000"] == 1;
	typenum = StringCount[num,"-"];
	syllabletype = (typenum == 0) || ((typenum == 1) && StringStartsQ[num,"-"]);
	syllablescliffordok = If[
		seperator && syllabletype,
			split = StringSplit[num,{"00000","-"}];
			{syllables,clifford} = If[Length[split] > 1,split,{"",split[[1]]}];
			syllablescharacterlist = Union[Characters[syllables]];
			syllablesok = SubsetQ[{"1","2","3","4","5","6","7","8","9","a","b","c","d","e","f"},syllablescharacterlist];
			cliffordcharacterlist = Union[Characters[clifford]];
			cliffordpossible = (StringLength[clifford]) <= 5 && SubsetQ[{"0","1","2","3","4","5","6","7","8","9","a","b","c","d","e","f"},cliffordcharacterlist];
			cliffordok = If[
				cliffordpossible,
					1 <= FromDigits[clifford,16] <= 92160,
					False
			];
			syllablesok && cliffordok,
			False
	];
	syllablescliffordok	
];
ValidHexDecQ[U_] := False

PhaseSet = {Id,MatrixPower[W,4],W,MatrixPower[W,5]};
L1Set = {Id,H1,S1.H1.MatrixPower[W,7]};
L2Set = {Id,H2,S2.H2.MatrixPower[W,7]};
CZSet = {Id,CZ.MatrixPower[W,7]};
S1Set = {Id,S1.MatrixPower[W,7]};
P1Set = {Id,X1,X1.Z1,Z1};
S2Set = {Id,S2.MatrixPower[W,7]};
P2Set = {Id,X2,X2.Z2,Z2};
SWAPSet = {Id,EX.MatrixPower[W,5]};
\[CapitalIota]Set = {Id,MatrixPower[W,2]};

CliffordFromNumber[cliffnum_] := Module[{digits,cz,l1,l2,index},
	digits = PadLeft[IntegerDigits[cliffnum-1,MixedRadix[{4,10,3,3,2,4,2,4,2,2}]],10];
	cz = Unitize[digits[[2]]];
	l2 = Mod[digits[[2]]-cz,3];
	l1 = (digits[[2]]-l2-cz)/3;
	index = Join[{digits[[1]],l1,l2,cz},digits[[3;;-1]]] + 1;
	PhaseSet[[index[[1]]]] . L1Set[[index[[2]]]] . L2Set[[index[[3]]]] . CZSet[[index[[4]]]] . 
		L1Set[[index[[5]]]] . L2Set[[index[[6]]]] . S1Set[[index[[7]]]] . P1Set[[index[[8]]]] .
		S2Set[[index[[9]]]] . P2Set[[index[[10]]]] . SWAPSet[[index[[11]]]] . \[CapitalIota]Set[[index[[12]]]]
];

FromHexDec[str_String/;ValidHexDecQ[str]] := Module[{syllablelist,split,syllablestr,cliffstr,cliffnum,cliff,syllables},
	syllablelist = If[StringStartsQ[str,"-"],SyllableListAsymmetric[[All,1]],SyllableList[[All,1]]];
	split = StringSplit[str,{"00000","-"}];
	{syllablestr,cliffstr} = If[Length[split] > 1,split,{"",split[[1]]}];
	cliffnum = FromDigits[cliffstr,16];
	cliff = CliffordFromNumber[cliffnum];
	syllables = Map[syllablelist[[FromDigits[#,16]]]&,Characters[syllablestr]];
	Simplify[Dot @@ (Append[syllables,cliff])]
];
FromHexDec[str_] := (
	Message[FromHexDec::invalidnumber,str];
	$Failed
);


(* ::Subsection::Closed:: *)
(*The Two-Qubit Clifford Group*)


(* ::Text:: *)
(*Here we develop some basic functions for the two-qubit Clifford group.*)


(* ::Subsubsection::Closed:: *)
(*Identifying if an operator is a Clifford and if two operators are Clifford-similar*)


(* ::Text:: *)
(*Boolean functions for checking if an operator is a Clifford or if two operators are Clifford-similar.*)


CliffordQ[U_/;U4Q[U]] := (Det[U]^2 == 1) && AllTrue[Flatten[U4ToSO6[U]],IntegerQ];
CliffordQ[U_/;SO6Q[U]] := AllTrue[Flatten[U],IntegerQ];
CliffordQ[U_String/;ValidHexDecQ[U]] := CliffordQ[FromHexDec[U]];
CliffordQ[U_String] := CliffordQ[FromSequence[U]];
CliffordQ[U_] := (
	Message[CliffordQ::notacircuit];
	$Failed
);

RightCliffordSimilar[U_String/;ValidHexDecQ[U],V_String/;ValidHexDecQ[V]] := RightCliffordSimilar[FromHexDec[U],FromHexDec[V]];
RightCliffordSimilar[U_String,V_String] := RightCliffordSimilar[FromSequence[U],FromSequence[V]];
RightCliffordSimilar[U_,V_] := CliffordQ[ConjugateTranspose[U].V];

LeftCliffordSimilar[U_String/;ValidHexDecQ[U],V_String/;ValidHexDecQ[V]] := LeftCliffordSimilar[FromHexDec[U],FromHexDec[V]];
LeftCliffordSimilar[U_String,V_String] := LeftCliffordSimilar[FromSequence[U],FromSequence[V]];
LeftCliffordSimilar[U_,V_] := CliffordQ[V.ConjugateTranspose[U]];


(* ::Subsubsection:: *)
(*Synthesis of Clifford Circuits (with at most 1 CZ gate and 1 SWAP gate)*)


(* ::Text:: *)
(*Rather than carry around an explicit lookup table with 92160 elements in the U(4) case and 23040 elements in the SO(6) case, we instead will use a synthesis algorithm based on the special representation; this takes advantage of the sparsity of Cliffords in this representation. We implicitly define our regular Clifford synthesis algorithm, CliffordSynth, in terms of this other algorithm, to be defined below.*)


PhaseFixU4[{str_,power_,SUcliffnum_},U_] := Module[{unitary,phase,Wpower,diff,Wcosetnumber},
	unitary = FromSequence[str];
	phase = Simplify[(U. ConjugateTranspose[unitary])[[1,1]]];
	Wpower = Mod[Log[(1+I)/Sqrt[2],phase],8];
	diff = Mod[Wpower - power,8];
	Wcosetnumber = Which[
		(diff == 0) || (diff == 2),0,
		(diff == 4) || (diff == 6),1,
		(diff == 1) || (diff == 3),2,
		True,3
	];
	{IntegerString[23040*Wcosetnumber + FromDigits[SUcliffnum,16],16],str<>StringRepeat["W ",Wpower]}
];
PhaseFixSO6[{str_,power_,SUcliffnum_}] := {SUcliffnum,str<>StringRepeat[" W",power]};

CliffordSynth[U_String/;ValidHexDecQ[U]] := Module[{unitary},
	unitary = FromSequence[U];
	CliffordSynth[unitary]
];
CliffordSynth[U_String] := Module[{unitary},
	unitary = FromSequence[U];
	CliffordSynth[unitary]
];
CliffordSynth[U_/;(U4Q[U] && CliffordQ[U])] := Module[{orthogonal,synth},
	orthogonal = U4ToSO6[U];
	synth = CliffordSynthSp[{0,SparseArray[orthogonal]}];
	PhaseFixU4[synth,U]
];
CliffordSynth[U_/;CliffordQ[U]] := PhaseFixSO6[CliffordSynthSp[{0,SparseArray[U]}]];
CliffordSynth[U_] := (
	Message[CliffordSynth::notaclifford];
	$Failed
);


(* ::Subsubsection:: *)
(*Clifford Group in the Special representation*)


(* ::Text:: *)
(*Explicit CliffordSynthSp algorithm. We perform the synthesis by decomposing an element of the signed permutation group in terms of a generating set which is constructed from basic Clifford operators in the SO(6) representation.*)


L1Cosets = {
	{IdSp,"",0,0},
	{H1Sp,"H1 ",0,1},
	{DotSp[S1Sp,H1Sp],"S1 H1 ",7,2}
};
L2Cosets = {
	{IdSp,"",0,0},
	{H2Sp,"H2 ",0,1},
	{DotSp[S2Sp,H2Sp],"S2 H2 ",7,2}
};
CZSets = {
	{IdSp,"",0,0},
	{CZSp,"CZ ",7,1}
};
S1Sets = {
	{IdSp,"",0,0},
	{S1Sp,"S1 ",7,1},
};
Pauli1Sets = {
	{IdSp,"",0,0},
	{X1Sp,"X1 ",0,1},
	{DotSp[X1Sp,Z1Sp],"X1 Z1 ",0,2},
	{Z1Sp,"Z1 ",0,3}
};
S2Sets = {
	{IdSp,"",0,0},
	{S2Sp,"S2 ",7,1},
};
Pauli2Sets = {
	{IdSp,"",0,0},
	{X2Sp,"X2 ",0,1},
	{DotSp[X2Sp,Z2Sp],"X2 Z2 ",0,2},
	{Z2Sp,"Z2 ",0,3}
};
SwapSets = {
	{IdSp,"",0,0},
	{EXSp,"EX ",5,1}
};
\[CapitalIota]Sets = {
	{IdSp,"",0,0},
	{\[CapitalIota]Sp,"",2,1}
};

CliffordQSp[U_] := Module[{k,M},
	{k,M} = KReduce[U];
	(k == 0) && OrthogonalMatrixQ[M] && (Det[M] == 1)
];

CliffordSynthSp[{_,M_}] := Module[
	{
		CliffordList,SwapTest,MUnSwapped,UpperLeft,LowerRight,RemovedLeftCosets,CZRows,RemovedCZSet,
		UpperLeftCol,LowerRightCol,RemovedLeftCosetsNew,S1Test,S2Test,RemovedSSets,d1,d2,d3,
		RemovedPauli1Sets,d3new,d4,d5,d6,RemovedPauli2Sets,\[CapitalIota]Test,str,phase,list,firstdigit,cliffordnumber
	},
	
	CliffordList = ConstantArray[0,11];
	
	SwapTest = Total[Abs[M[[1;;3,1;;3]]],2];
	CliffordList[[10]] = If[SwapTest > 1,SwapSets[[1]],SwapSets[[2]]];
	MUnSwapped = M.Transpose[CliffordList[[10,1,2]]];
	
	UpperLeft = Total[Abs[MUnSwapped[[1;;3,1;;3]]],{2}];
	LowerRight = Total[Abs[MUnSwapped[[4;;6,4;;6]]],{2}];
	CliffordList[[1]] = Which[
		UpperLeft == {0,1,1},L1Cosets[[2]],
		UpperLeft == {1,0,1},L1Cosets[[3]],
		True,L1Cosets[[1]]
	];
	CliffordList[[2]] = Which[
		LowerRight == {0,1,1},L2Cosets[[2]],
		LowerRight == {1,0,1},L2Cosets[[3]],
		True,L2Cosets[[1]]
	];
	RemovedLeftCosets = Transpose[CliffordList[[1,1,2]].CliffordList[[2,1,2]]].MUnSwapped;
	
	CZRows = Total[Abs[RemovedLeftCosets[[3,1;;3]]],2];
	CliffordList[[3]] = If[CZRows == 1,CZSets[[1]],CZSets[[2]]];
	RemovedCZSet = Transpose[CliffordList[[3,1,2]]].RemovedLeftCosets;
	
	UpperLeftCol = Abs[RemovedCZSet[[1;;3,3]]];
	LowerRightCol = Abs[RemovedCZSet[[4;;6,6]]];
	CliffordList[[4]] = Which[
		UpperLeftCol == {0,0,1},L1Cosets[[1]],
		UpperLeftCol == {0,1,0},L1Cosets[[3]],
		True,L1Cosets[[2]]
	];
	CliffordList[[5]] = Which[
		LowerRightCol == {0,0,1},L2Cosets[[1]],
		LowerRightCol == {0,1,0},L2Cosets[[3]],
		True,L2Cosets[[2]]
	];
	RemovedLeftCosetsNew = Transpose[CliffordList[[4,1,2]].CliffordList[[5,1,2]]].RemovedCZSet;
	
	S1Test = Abs[RemovedLeftCosetsNew[[1,1]]];
	S2Test = Abs[RemovedLeftCosetsNew[[4,4]]];
	CliffordList[[6]] = If[S1Test == 1,S1Sets[[1]],S1Sets[[2]]];
	CliffordList[[8]] = If[S2Test == 1,S2Sets[[1]],S2Sets[[2]]];
	RemovedSSets = Transpose[CliffordList[[6,1,2]].CliffordList[[8,1,2]]].RemovedLeftCosetsNew;
	
	{d1,d2,d3} = Diagonal[RemovedSSets[[1;;3,1;;3]]];
	CliffordList[[7]] = Which[
		(d1 == d2) && (d2 == d3),Pauli1Sets[[1]],
		(d1 != d2) && (d2 == d3),Pauli1Sets[[2]],
		(d1 != d2) && (d2 != d3),Pauli1Sets[[3]],
		True,Pauli1Sets[[4]]
	];
	RemovedPauli1Sets = Transpose[CliffordList[[7,1,2]]].RemovedSSets;
	
	{d3new,d4,d5,d6} = Diagonal[RemovedPauli1Sets[[3;;6,3;;6]]];
	CliffordList[[9]] = Which[
		(d4 == d5) && (d5 == d6) && (d3new*d4 == 1),Pauli2Sets[[1]],
		(d4 != d5) && (d5 == d6),Pauli2Sets[[2]],
		(d4 != d5) && (d5 != d6),Pauli2Sets[[3]],
		True,Pauli2Sets[[4]]
	];
	RemovedPauli2Sets = Transpose[CliffordList[[9,1,2]]].RemovedPauli1Sets;
	
	\[CapitalIota]Test = RemovedPauli2Sets[[1,1]];
	CliffordList[[11]] = If[\[CapitalIota]Test == 1,\[CapitalIota]Sets[[1]],\[CapitalIota]Sets[[2]]];
	
	{str,phase,list} = Fold[{
			#1[[1]]<>#2[[2]],
			Mod[#1[[2]]+#2[[3]],8],
			Append[#1[[3]],#2[[4]]]
		}&,
		{"",0,{}},
		CliffordList
	];
	firstdigit = list[[3]] * (list[[1]]*3 + list[[2]] + 1);
	cliffordnumber = FromDigits[Prepend[list[[4;;-1]],firstdigit],MixedRadix[{10,3,3,2,4,2,4,2,2}]] + 1;
	
	{str,phase,IntegerString[cliffordnumber,16]}
];


(* ::Subsection::Closed:: *)
(*The Two-Qubit Clifford + CS Group*)


(* ::Text:: *)
(*Here we develop some basic constructors for Clifford + CS circuits. We also provide a function for checking if an operator corresponds to a Clifford + CS circuit.*)


(* ::Subsubsection:: *)
(*Check If an Operator Is a Gaussian Clifford + T Matrix*)


(* ::Text:: *)
(*We define a function for determining if an operator is a representation for a Gaussian Clifford + T matrix, i.e. a Clifford + CS operator.*)


DyadicQ[n_] := Module[{numerator,denominator},
	{numerator,denominator} = NumeratorDenominator[n];
	IntegerQ[numerator] && IntegerQ[Log2[denominator]]
];
GaussianDyadicQ[n_]:= AllTrue[ReIm[n],DyadicQ];

GaussianQ[U_/;U4Q[U]] := Module[{UGaussDyadicQ,UWGuassDyadicQ},
	UGaussDyadicQ = AllTrue[Flatten[U],GaussianDyadicQ];
	UWGuassDyadicQ = AllTrue[Flatten[W.U],GaussianDyadicQ];
	UGaussDyadicQ || UWGuassDyadicQ
];
GaussianQ[U_/;SO6Q[U]] := Module[{sp},
	sp = SO6ToSpecialRep[U];
	(Det[U] == 1) && IntegerQ[sp[[1]]] && AllTrue[Flatten[sp[[2]]],IntegerQ]
];
GaussianQ[U_/;ValidHexDecQ[U]] := True;
GaussianQ[U_String] := GaussianQ[FromSequence[U]];
GaussianQ[U_] := (Message[GaussianQ::notacircuit];
	$Failed
);


(* ::Subsubsection:: *)
(*Syllable Lists*)


(* ::Text:: *)
(*For the exported versions:*)


PrePostfixList = {
	{H1.H2,H1SO6.H2SO6,"H1 H2 ","H2 H1 "},
	{S1.H1.S2.H2,S1SO6.H1SO6.S2SO6.H2SO6,"S1 H1 S2 H2 ","H1 S1 S1 S1 H2 S2 S2 S2 "},
	{Id,IdSO6,"",""},
	{S1.H1,S1SO6.H1SO6,"S1 H1 ","H1 S1 S1 S1 "},
	{S2.H2,S2SO6.H2SO6,"S2 H2 ","H2 S2 S2 S2 "},
	{H2,H2SO6,"H2 ","H2 "},
	{H1,H1SO6,"H1 ","H1 "},
	{H1.S2.H2,H1SO6.S2SO6.H2SO6,"H1 S2 H2 ","H1 H2 S2 S2 S2 "},
	{S1.H1.H2,S1SO6.H1SO6.H2SO6,"S1 H1 H2 ","H1 S1 S1 S1 H2 "},
	{CNOT12.H1,CNOT12SO6.H1SO6,"H2 CZ H1 H2 ","H1 H2 CZ H2 "},
	{CZ.S1.H1.S2.H2,CZSO6.S1SO6.H1SO6.S2SO6.H2SO6,"CZ S1 H1 S2 H2 ","H1 S1 S1 S1 H2 S2 S2 S2 CZ "},
	{CZ.H1.H2,CZSO6.H1SO6.H2SO6,"CZ H1 H2 ","H1 H2 CZ "},
	{CNOT12.S1.H1,CNOT12SO6.S1SO6.H1SO6,"H2 CZ S1 H1 H2 ","H1 S1 S1 S1 H2 CZ H2 "},
	{CZ.S1.H1.H2,CZSO6.S1SO6.H1SO6.H2SO6,"CZ S1 H1 H2 ","H1 S1 S1 S1 H2 CZ "},
	{CZ.H1.S2.H2,CZSO6.H1SO6.S2SO6.H2SO6,"CZ H1 S2 H2 ","H1 H2 S2 S2 S2 CZ "}
};

SyllableList = Map[
	{#[[1]].CS.ConjugateTranspose[#[[1]]],#[[2]].CSSO6.Transpose[#[[2]]],#[[3]]<>"CS "<>#[[4]]}&,
	PrePostfixList
];
SyllableListAsymmetric = Map[
	{#[[1]].CS,#[[2]].CSSO6,#[[3]]<>"CS "}&,
	PrePostfixList
];


(* ::Text:: *)
(*And the internal data structure:*)


SyllableListSp = MapIndexed[
	{SO6ToSpecialRep[#1[[2]]],#1[[3]],IntegerString[#2,16]}&,
	SyllableList
];
SyllableListAsymmetricSp = MapIndexed[
	{SO6ToSpecialRep[#1[[2]]],#1[[3]],IntegerString[#2,16]}&,
	SyllableListAsymmetric
];


(* ::Subsection::Closed:: *)
(*Normalization*)


(* ::Text:: *)
(*In this section we develop a function for normalizing a circuit in any of our representations. This algorithm is described in detail in [2].*)


(* ::Subsubsection::Closed:: *)
(*Earliest Generator Ordering and associated row pairings*)


(* ::Text:: *)
(*We develop here an association which matches up each syllable to their row pairings under EGO.*)


TwoFourPaired = Map[
	{#,Complement[Range[1,6],#]}&,
	Subsets[Range[1,6],{4}]
];
TwoTwoTwoPaired = Flatten[Map[
	Table[{#[[1]],{#[[2,1]],#[[2,k]]},Complement[#[[2]],{#[[2,1]],#[[2,k]]}]},{k,2,4}]&,
	Table[{{1,j},Complement[Range[1,6],{1,j}]},{j,2,6}]
],1];
PairList = Map[Sort,Join[TwoTwoTwoPaired,TwoFourPaired]];

MatchingPairings[list_] := {
	list,
	{list[[3]],Union[list[[1]],list[[2]]]},
	{list[[2]],Union[list[[1]],list[[3]]]},
	{list[[1]],Union[list[[2]],list[[3]]]}
};
PossibleSyllableListPairings = Map[
	MatchingPairings[RowPairs[First[#]]]&,
	SyllableListSp
];
SyllableListPairings = Map[
	# -> First[FirstPosition[PossibleSyllableListPairings,#]]&,
	PairList
];
PairingKey = Association @@ SyllableListPairings;


(* ::Subsubsection::Closed:: *)
(*Finding a leftmost syllable and normalizing in the Special representation*)


(* ::Text:: *)
(*We develop separate synthesis algorithms here based on whether we want to synthesize a circuit using the symmetric or asymmetric syllables.*)


LeftmostSyllable[U_] := SyllableListSp[[PairingKey[RowPairs[U]]]];
RemoveLeftmost[{{k_,M_},str_,numberstr_}] := Module[{leftmost,clifford},
	Which[
		k > 0,
			leftmost = LeftmostSyllable[{k,M}];
			{Dot2Sp[InvSp[leftmost[[1]]],{k,M}],str<>leftmost[[2]],numberstr<>leftmost[[3]]},
		k == 0,
			clifford = CliffordSynthSp[{k,M}];
			{{-1,SparseArray[ConstantArray[0,{6,6}]]},str<>clifford[[1]]<>StringRepeat["W ",clifford[[2]]],numberstr<>"00000"<>clifford[[3]]},
		True,
			{{k,M},str,numberstr}
	]
];

LeftmostSyllableAsymmetric[U_] := SyllableListAsymmetricSp[[PairingKey[RowPairs[U]]]];
RemoveLeftmostAsymmetric[{{k_,M_},str_,numberstr_}] := Module[{leftmost,clifford},
	Which[
		k > 0,
			leftmost = LeftmostSyllableAsymmetric[{k,M}];
			{Dot2Sp[InvSp[leftmost[[1]]],{k,M}],str<>leftmost[[2]],numberstr<>leftmost[[3]]},
		k == 0,
			clifford = CliffordSynthSp[{k,M}];
			{{-1,SparseArray[ConstantArray[0,{6,6}]]},str<>clifford[[1]]<>StringRepeat["W ",clifford[[2]]],numberstr<>"00000"<>clifford[[3]]},
		True,
			{{k,M},str,numberstr}
	]
];

test[NormItSp,"SyllableType"] := MemberQ[{"Normal","Asymmetric"},#]&;

Options[NormItSp] = {"SyllableType" -> "Normal"};
NormItSp[U_,OptionsPattern[]]?optsCheck := Module[{type},
	type = If[OptionValue["SyllableType"] == "Normal",RemoveLeftmost,RemoveLeftmostAsymmetric];
	FixedPoint[type,{U,"",If[OptionValue["SyllableType"] == "Normal","","-"]}][[2;;3]]
];


(* ::Subsubsection::Closed:: *)
(*Normalizing from a string, a hexadecimal, U(4), or SO(6).*)


(* ::Text:: *)
(*We provide the function NormIt for normalizing operators in any of our forms.*)


Options[FixPhase] = {"UpToPhase" -> True};
FixPhase[{str_,numberstr_},U_,OptionsPattern[]] := Module[{syllablesplit,syllablenum,syllablepart,clifford,cliffnum,cliffstr,numberstrfixed,nophasestr,phasecount,strfixed},
	If[
		OptionValue["UpToPhase"] == True,
			{str,numberstr},
			syllablesplit = StringSplit[numberstr,"00000"];
			syllablenum = If[Length[syllablesplit] > 1, syllablesplit[[1]],""];
			syllablepart = FromHexDec[syllablenum<>"000001"];
			clifford = Simplify[ConjugateTranspose[syllablepart] . U];
			{cliffnum,cliffstr} = CliffordSynth[clifford];
			numberstrfixed = syllablenum<>"00000"<>cliffnum;
			nophasestr = StringDelete[str,"W "];
			phasecount = StringCount[cliffstr,"W "];
			strfixed = nophasestr<>StringRepeat["W ",phasecount];
			{strfixed,numberstrfixed}
	]
];

test[NormIt,"SyllableType"] := MemberQ[{"Normal","Asymmetric"},#]&;
test[NormIt,"OutputType"] := MemberQ[{"String","HexDec"},#]&;
test[NormIt,"UpToPhase"] := BooleanQ;

Options[NormIt] = {"OutputType" -> "String","SyllableType" -> "Normal","UpToPhase" -> True};
NormIt[U_/;(GaussianQ[U] && U4Q[U]),OptionsPattern[]]?optsCheck := Module[{index},
	index = If[OptionValue["OutputType"] == "String",1,2];
	FixPhase[NormItSp[SO6ToSpecialRep[U4ToSO6[U]],"SyllableType" -> OptionValue["SyllableType"]],U,"UpToPhase" -> OptionValue["UpToPhase"]][[index]]
];
NormIt[U_/;(GaussianQ[U] && SO6Q[U]),OptionsPattern[]]?optsCheck := Module[{index},
	index = If[OptionValue["OutputType"] == "String",1,2];
	NormItSp[SO6ToSpecialRep[U],"SyllableType" -> OptionValue["SyllableType"]][[index]]
];
NormIt[hexdec_String/;ValidHexDecQ[hexdec],opts:OptionsPattern[]] := NormIt[FromHexDec[hexdec],opts];
NormIt[str_String,opts:OptionsPattern[]] := NormIt[FromSequence[str],opts];
NormIt[U_] := (
	Message[NormIt::notacircuit,U];
	$Failed
);


(* ::Subsection:: *)
(*\[Epsilon]-Approximating Pauli Rotations*)


(* ::Text:: *)
(*This section describes algorithms used for finding approximations to Pauli Rotations*)


(* ::Subsubsection:: *)
(*Frobenius Distance of Matrices*)


(* ::Text:: *)
(*Computes the Frobenius distance between two matrices A and B (i.e. the Frobenius norm of the difference between A and B)*)


FrobeniusDistance[A_/;MatrixQ[A],B_/;MatrixQ[B]]/;(Dimensions[A] == Dimensions[B]) := Simplify[Sqrt[Re[Tr[(A-B).ConjugateTranspose[(A-B)]]]]];
FrobeniusDistance[_,_] := (
	Message[FrobeniusDistance::notequidimensionalmatrices];
	$Failed
);


(* ::Subsubsection:: *)
(*Continued fractions and affine transformation*)


(* ::Text:: *)
(*We describe a scheme for finding an appropriate affine transformation.*)


NextIterate[{aN_,rN_,pN_,qN_,pNm1_,qNm1_}] := Module[{aNp1},
	aNp1 = IntegerPart[1/rN];
	{aNp1,1/rN - aNp1,aNp1*pN + pNm1,aNp1*qN + qNm1,pN,qN}
];
ContinuedFractionToPrecision[\[Alpha]_,tol_] := Module[{a0,iterate},
	a0 = IntegerPart[\[Alpha]];
	iterate = NestWhile[NextIterate,{a0,N[\[Alpha] - a0,2*Log10[1/tol]],a0,1,1,0},#[[4]] <= 1/tol && #[[2]] != 0&];
	If[iterate[[2]] == 0 && iterate[[4]] <= 1/tol,
		iterate[[3;;4]],
		iterate[[5;;6]]
	]
];

AffineMatrix[\[CurlyPhi]_,\[Epsilon]_] := Module[{\[Alpha],q,r,s,t},
	\[Alpha] = Tan[\[CurlyPhi]/2];
	{r,q} = ContinuedFractionToPrecision[\[Alpha],Sqrt[\[Epsilon]/(8+\[Epsilon]*\[Alpha])]*Sec[\[CurlyPhi]/2]];
	{t,s} = ExtendedGCD[q,-r][[2]];
	{{q,r},{s,t}}
];


(* ::Subsubsection:: *)
(*Finding angle candidates for bounded and unbounded angles*)


(* ::Text:: *)
(*For calculating a candidate solution, we first map every angle into the interval [0,\[Pi]/2]. We then use our affine transformation to find solutions to the integer programming problem.*)


RealQ[x_] := Element[x,Reals];

CandidateFinder[\[CurlyPhi]_?RealQ,\[Epsilon]_?RealQ] := Module[{mod\[CurlyPhi]},
	mod\[CurlyPhi]=Mod[\[CurlyPhi],4*Pi];
	Which[
		0<=mod\[CurlyPhi]<=Pi/2,
			CandidateFinderBounded[mod\[CurlyPhi],\[Epsilon]],
		Pi/2<mod\[CurlyPhi]<=Pi,
			MapAt[Reverse,CandidateFinderBounded[Pi-mod\[CurlyPhi],\[Epsilon]],2],
		Pi<mod\[CurlyPhi]<=2*Pi,
			MapAt[{-#[[2]],#[[1]]}&,CandidateFinder[mod\[CurlyPhi]-Pi,\[Epsilon]],2],
		True,
			MapAt[{-#[[1]],-#[[2]]}&,CandidateFinder[mod\[CurlyPhi]-2*Pi,\[Epsilon]],2]
	]
];
CandidateFinder[_,_] := (
	Message[CandidateFinder::notreals];
	$Failed
);

CandidateFinderBounded[\[CurlyPhi]_,\[Epsilon]_] := Module[{root2k,c,s,\[Alpha],p,\[CapitalDelta],\[Delta],A,invA,scaledp\[Prime],scaled\[CapitalDelta]\[Prime],scaled\[Delta]\[Prime],\[Theta]\[CapitalDelta]\[Prime],m1,m2,x1,x2,x3,x4,y0,y1,transformed,valid},
	root2k = 1;
	c = Cos[\[CurlyPhi]/2];
	s = Sin[\[CurlyPhi]/2];
	\[Alpha] = 1-\[Epsilon]^2/8;
	p = {c*\[Alpha]+s*Sqrt[1-\[Alpha]^2],s*\[Alpha]-c*Sqrt[1-\[Alpha]^2]};
	\[CapitalDelta] = 2*Sqrt[1-\[Alpha]^2]*{-s,c};
	\[Delta] = \[Epsilon]^2/8*{c,s};
	A = AffineMatrix[\[CurlyPhi],\[Epsilon]];
	invA = Inverse[A];
	scaledp\[Prime] = A.p;
	scaled\[CapitalDelta]\[Prime] = A.\[CapitalDelta];
	scaled\[Delta]\[Prime] = A.\[Delta];
	\[Theta]\[CapitalDelta]\[Prime] = VectorAngle[scaled\[CapitalDelta]\[Prime],{1,0}];
	Which[
		\[Theta]\[CapitalDelta]\[Prime] > Pi/2,
			m1 = scaled\[CapitalDelta]\[Prime][[2]] / scaled\[CapitalDelta]\[Prime][[1]];
			m2 = scaled\[Delta]\[Prime][[2]] / scaled\[Delta]\[Prime][[1]];
			valid = Catch[Do[
				{x1,x2,x3,x4} = {
					Ceiling[scaledp\[Prime][[1]]+scaled\[CapitalDelta]\[Prime][[1]]],
					Floor[scaledp\[Prime][[1]]],
					Floor[scaledp\[Prime][[1]]+scaled\[CapitalDelta]\[Prime][[1]]+scaled\[Delta]\[Prime][[1]]],
					Floor[scaledp\[Prime][[1]]+scaled\[Delta]\[Prime][[1]]]
				};

				Do[
					y0 = If[x <= x2,
						Ceiling[m1*(x - scaledp\[Prime][[1]]) + scaledp\[Prime][[2]]],
						Ceiling[m2*(x - scaledp\[Prime][[1]]) + scaledp\[Prime][[2]]]
					];
					y1 = If[x <= x3,
						Floor[m2*(x - scaledp\[Prime][[1]]-scaled\[CapitalDelta]\[Prime][[1]]) + scaledp\[Prime][[2]] + scaled\[CapitalDelta]\[Prime][[2]]],
						Floor[m1*(x - scaledp\[Prime][[1]]-scaled\[Delta]\[Prime][[1]]) + scaledp\[Prime][[2]] + scaled\[Delta]\[Prime][[2]]]
					];
					Do[
						transformed = invA.{x,y};
						If[Norm[transformed] <= Sqrt[2]^k,Throw[{k,transformed}]],
						{y,y0,y1}
					];,
					{x,x1,x4}
				];
				{scaledp\[Prime],scaled\[CapitalDelta]\[Prime],scaled\[Delta]\[Prime]} *= Sqrt[2];,
				{k,0,Ceiling[4*Log2[1/\[Epsilon]]+6]}
			];];,
		\[Theta]\[CapitalDelta]\[Prime] == Pi/2,
			m1 = scaled\[Delta]\[Prime][[2]] / scaled\[Delta]\[Prime][[1]];
			valid = Catch[Do[
				{x1,x2} = {
					Ceiling[scaledp\[Prime][[1]]],
					Floor[scaledp\[Prime][[1]]+scaled\[Delta]\[Prime][[1]]]
				};
				Do[
					y0 = Ceiling[m1*(x - scaledp\[Prime][[1]]) + scaledp\[Prime][[2]]];
					y1 = Floor[m1*(x - scaledp\[Prime][[1]]-scaled\[CapitalDelta]\[Prime][[1]]) + scaledp\[Prime][[2]] + scaled\[CapitalDelta]\[Prime][[2]]];
					Do[
						transformed = invA.{x,y};
						If[Norm[transformed] <= Sqrt[2]^k,Throw[{k,transformed}]],
						{y,Floor[(y0+y1)/2],y1}
					];,
					{x,x1,x2}
				];
				{scaledp\[Prime],scaled\[CapitalDelta]\[Prime],scaled\[Delta]\[Prime]} *= Sqrt[2];,
				{k,0,Ceiling[4*Log2[1/\[Epsilon]]+6]}
			];];,
		True,
			m1 = scaled\[Delta]\[Prime][[2]] / scaled\[Delta]\[Prime][[1]];
			m2 = scaled\[CapitalDelta]\[Prime][[2]] / scaled\[CapitalDelta]\[Prime][[1]];
			valid = Catch[Do[
				{x1,x2,x3,x4} = {
					Ceiling[scaledp\[Prime][[1]]],
					Floor[scaledp\[Prime][[1]] + scaled\[Delta]\[Prime][[1]]],
					Floor[scaledp\[Prime][[1]] + scaled\[CapitalDelta]\[Prime][[1]]],
					Floor[scaledp\[Prime][[1]] + scaled\[Delta]\[Prime][[1]] + scaled\[CapitalDelta]\[Prime][[1]]]
				};
				Do[
					y0 = If[x <= x2,
						Ceiling[m1*(x - scaledp\[Prime][[1]]) + scaledp\[Prime][[2]]],
						Ceiling[m2*(x - scaledp\[Prime][[1]] - scaled\[Delta]\[Prime][[1]]) + scaledp\[Prime][[2]] + -scaled\[Delta]\[Prime][[2]]]
					];
					y1 = If[x <= x3,
						Floor[m2*(x - scaledp\[Prime][[1]]) + scaledp\[Prime][[2]]],
						Floor[m1*(x - scaledp\[Prime][[1]] - scaled\[CapitalDelta]\[Prime][[1]]) + scaledp\[Prime][[2]] + scaled\[CapitalDelta]\[Prime][[2]]]
					];
					Do[
						transformed = invA.{x,y};
						If[Norm[transformed] <= Sqrt[2]^k,Throw[{k,transformed}]],
						{y,y0,y1}
					];,
					{x,x1,x4}
				];
				{scaledp\[Prime],scaled\[CapitalDelta]\[Prime],scaled\[Delta]\[Prime]} *= Sqrt[2];,
				{k,0,Ceiling[4*Log2[1/\[Epsilon]]+6]}
			];];
	];
	valid
];

CandidateFinderBoundedNewApproach[\[CurlyPhi]_,\[Epsilon]_] := Module[{c,s,\[Alpha],p,\[CapitalDelta],\[Delta],A,v,invA,p\[Prime],\[CapitalDelta]\[Prime],\[Delta]\[Prime],v\[Prime],mInv,xL,xR,xB,xMin,xMax,yMin,yMax,kMax,k0,dx,dy,root2k,round,x0,x1,y0,y1,x0Int,x1Int,x2Int,y0Int,y1Int,y,vec},
	Catch[
		If[\[Epsilon] >= 2*Sqrt[2],Throw[{0,{1,0}}]];
		c = Cos[\[CurlyPhi]/2];
		s = Sin[\[CurlyPhi]/2];
		\[Alpha] = 1-\[Epsilon]^2/8;
		p = {c*\[Alpha]+s*Sqrt[1-\[Alpha]^2],s*\[Alpha]-c*Sqrt[1-\[Alpha]^2]};
		\[CapitalDelta] = 2*Sqrt[1-\[Alpha]^2]*{-s,c};
		\[Delta] = \[Epsilon]^2/8*{c,s};
		A = AffineMatrix[\[CurlyPhi],\[Epsilon]];
		v = Normalize[A[[1]]];
		invA = Inverse[A];
		p\[Prime] = A.p;
		\[CapitalDelta]\[Prime] = A.\[CapitalDelta];
		\[Delta]\[Prime] = A.\[Delta];
		v\[Prime] = A.v;
		mInv = \[CapitalDelta]\[Prime][[1]] / \[CapitalDelta]\[Prime][[2]];
		{xL,xR} = If[\[CapitalDelta]\[Prime][[1]] <= 0,
			{\[CapitalDelta]\[Prime] + p\[Prime],p\[Prime]},
			{p\[Prime],\[CapitalDelta]\[Prime] + p\[Prime]}
		];
		xB = If[(\[CurlyPhi]/2 - ArcCos[\[Alpha]]) < (ArcTan @@ v) < (\[CurlyPhi]/2 + ArcCos[\[Alpha]]),
			v\[Prime],
			xR
		];
		{xMin,xMax} = {xL[[1]],xR[[1]] + \[Delta]\[Prime][[1]]};
		{yMin,yMax} = If[\[CapitalDelta]\[Prime][[1]] <= 0,{xR[[2]],xL[[2]]},{xL[[2]],xR[[2]]}] + If[\[Delta]\[Prime][[2]]>0,{0,\[Delta]\[Prime][[2]]},{\[Delta]\[Prime][[2]],0}];
		kMax = Ceiling[4*Log2[1/\[Epsilon]]+6];
		k0 = Catch[Do[
			dx = Floor[xMax] - Ceiling[xMin] + 1;
			dy = Floor[yMax] - Ceiling[yMin] + 1;
			If[dx * dy > 0,Throw[k]];
			{xMin,xMax,yMin,yMax} *= Sqrt[2];,
			{k,0,kMax}
		]];
		root2k = Sqrt[2]^k0;
		round = If[\[CapitalDelta]\[Prime][[1]] <= 0,Ceiling,Floor];
		
		If[
			xB[[1]] > xR[[1]],
				If[mInv == 0,
					{x0,y0} = root2k * {xR[[1]],xB[[2]]};
					Do[
						{x0Int,y0Int} = {Ceiling[x0],Floor[y0]};
						Do[
							vec = invA.{x0Int,y};
							If[Norm[vec] <= root2k,Throw[{k,vec}]];,
							{y,y0Int,y0Int+1}
						];
						{root2k,x0,y0} *= Sqrt[2];,
						{k,k0,kMax}
					];,
					{x0,y0,x1,y1} = root2k * {xL[[1]],xL[[2]],xR[[1]],xB[[2]]};
					Do[
						{x0Int,x1Int} = {Ceiling[x0],Floor[x1]};
						Do[
							y = round[(x - x0)/mInv + y0];
							vec = invA.{x,y};
							If[Norm[vec] <= root2k,Throw[{k,vec}]];,
							{x,x0Int,x1Int}
						];
						x2Int = x1Int + 1;
						y1Int = Floor[y1];
						Do[
							vec = invA.{x2Int,y};
							If[Norm[vec] <= root2k,Throw[{k,vec}]];,
							{y,y1Int,y1Int+1}
						];
						{root2k,x0,y0,x1,y1} *= Sqrt[2];,
						{k,k0,kMax}
					];
				];,
				{x0,y0,x1} = root2k * {xL[[1]],xL[[2]],xR[[1]]};
				Do[
					{x0Int,x1Int} = {Ceiling[x0],Floor[x1]};
					Do[
						y = round[(x - x0)/mInv + y0];
						vec = invA.{x,y};
						If[Norm[vec] <= root2k,Throw[{k,vec}]];,
						{x,x0Int,x1Int}
					];
					{root2k,x0,y0,x1} *= Sqrt[2];,
					{k,k0,kMax}
				];
		];
	]
];


(* ::Subsubsection:: *)
(*Solving Lagrange Four-Squares*)


(* ::Text:: *)
(*Rather than implement our own solver based off of well known algorithms, we instead use a basic Mathematica function as the inputs for this problem never get too big.*)


Lagrange4[n_]:=Module[{x1,x2,x3,x4},({x1,x2,x3,x4}/.FindInstance[x1^2+x2^2+x3^2+x4^2==n,{x1,x2,x3,x4},Integers,1])[[1]]];


(* ::Subsubsection:: *)
(*SU(4) Approximations Using a Candidate Solution*)


(* ::Text:: *)
(*We provide an algorithm for finding a rotation by angle \[CurlyPhi] up to error \[Epsilon] for any of the 15 Pauli matrices.*)


SU4Z1Finder[\[CurlyPhi]_,\[Epsilon]_]:=Module[{k,x,y,a,b,c,d,Invroot2k,\[Alpha],\[Beta],\[Chi]},
	{k,{x,y}} = CandidateFinder[\[CurlyPhi],\[Epsilon]];
	{a,b,c,d} = Lagrange4[2^k-x^2-y^2];
	Invroot2k = 1/Sqrt[2]^k;
	\[Alpha] = Invroot2k*(x+I*y);
	\[Beta] = Invroot2k*(a+I*b);
	\[Chi] = Invroot2k*(c+I*d);
	{
		{\[Alpha],0,-Conjugate[\[Beta]],-Conjugate[\[Chi]]},
		{0,\[Alpha],\[Chi],-\[Beta]},
		{\[Beta],-Conjugate[\[Chi]],Conjugate[\[Alpha]],0},
		{\[Chi],Conjugate[\[Beta]],0,Conjugate[\[Alpha]]}
	}
];

PauliRotation[\[CurlyPhi]_/;Not[RealQ[\[CurlyPhi]]],\[Epsilon]_/;Not[RealQ[\[Epsilon]]],_] := (
	Message[PauliRotation::notreals];
	$Failed
);
PauliRotation[\[CurlyPhi]_,\[Epsilon]_,"ZI"] := SU4Z1Finder[\[CurlyPhi],\[Epsilon]];
PauliRotation[\[CurlyPhi]_,\[Epsilon]_,"XI"] := H1 . SU4Z1Finder[\[CurlyPhi],\[Epsilon]] . H1;
PauliRotation[\[CurlyPhi]_,\[Epsilon]_,"YI"] := S1 . H1 . SU4Z1Finder[\[CurlyPhi],\[Epsilon]] . H1 . ConjugateTranspose[S1];
PauliRotation[\[CurlyPhi]_,\[Epsilon]_,"IZ"] := EX . SU4Z1Finder[\[CurlyPhi],\[Epsilon]] . EX;
PauliRotation[\[CurlyPhi]_,\[Epsilon]_,"IX"] := EX . H1 . SU4Z1Finder[\[CurlyPhi],\[Epsilon]] . H1 . EX;
PauliRotation[\[CurlyPhi]_,\[Epsilon]_,"IY"] := EX . S1 . H1 . SU4Z1Finder[\[CurlyPhi],\[Epsilon]] . H1 . ConjugateTranspose[S1] . EX;
PauliRotation[\[CurlyPhi]_,\[Epsilon]_,"ZZ"] := CNOT21 . SU4Z1Finder[\[CurlyPhi],\[Epsilon]] . CNOT21;
PauliRotation[\[CurlyPhi]_,\[Epsilon]_,"XZ"] := H1 . CNOT21 . SU4Z1Finder[\[CurlyPhi],\[Epsilon]] . CNOT21 . H1;
PauliRotation[\[CurlyPhi]_,\[Epsilon]_,"YZ"] := S1 . H1 . CNOT21 . SU4Z1Finder[\[CurlyPhi],\[Epsilon]] . CNOT21 . H1 . ConjugateTranspose[S1];
PauliRotation[\[CurlyPhi]_,\[Epsilon]_,"ZX"] := H2 . CNOT21 . SU4Z1Finder[\[CurlyPhi],\[Epsilon]] . CNOT21 . H2;
PauliRotation[\[CurlyPhi]_,\[Epsilon]_,"XX"] := H1 . H2 . CNOT21 . SU4Z1Finder[\[CurlyPhi],\[Epsilon]] . CNOT21 . H2 . H1;
PauliRotation[\[CurlyPhi]_,\[Epsilon]_,"YX"] := S1 . H1 . H2 . CNOT21 . SU4Z1Finder[\[CurlyPhi],\[Epsilon]] . CNOT21 . H2 . H1 . ConjugateTranspose[S1];
PauliRotation[\[CurlyPhi]_,\[Epsilon]_,"ZY"] := S2 . H2 . CNOT21 . SU4Z1Finder[\[CurlyPhi],\[Epsilon]] . CNOT21 . H2 . ConjugateTranspose[S2];
PauliRotation[\[CurlyPhi]_,\[Epsilon]_,"XY"] := H1 . S2 . H2 . CNOT21 . SU4Z1Finder[\[CurlyPhi],\[Epsilon]] . CNOT21 . H2 . ConjugateTranspose[S2] . H1;
PauliRotation[\[CurlyPhi]_,\[Epsilon]_,"YY"] := S1 . H1 . S2 . H2 . CNOT21 . SU4Z1Finder[\[CurlyPhi],\[Epsilon]] . CNOT21 . H2 . ConjugateTranspose[S2] . H1 . ConjugateTranspose[S1];
PauliRotation[_,_,_] := (
	Message[PauliRotation::invldstring];
	$Failed
);

test[PauliRotationSequence,"SyllableType"] := MemberQ[{"Normal","Asymmetric"},#]&;
test[PauliRotationSequence,"OutputType"] := MemberQ[{"String","HexDec"},#]&;

Options[PauliRotationSequence] = {"OutputType" -> "String","SyllableType" -> "Normal"};

PauliRotationSequence[\[CurlyPhi]_,\[Epsilon]_,pauli_,opts:OptionsPattern[]]?optsCheck := NormIt[PauliRotation[\[CurlyPhi],\[Epsilon],pauli],opts,"UpToPhase" -> False];


(* ::Subsection:: *)
(*General Unitary Decomposition*)


(* ::Text:: *)
(*Here we develop an algorithm for the approximation of any U(4) matrix using the Clifford + CS gate set.*)


(* ::Subsubsection:: *)
(*Angles of rotation in the Pauli decomposition*)


(* ::Text:: *)
(*First, we develop a method for finding the Pauli decomposition of an operator in terms of its 15 Pauli-rotation angles.*)


SafeArcTan[a_,b_] := If[(a == 0) && (b == 0),0,ArcTan[a,b]];
AngleFind[{{a_,_},{b_,_}}] := SafeArcTan[a,b];
SmallPauliDecomposition[o_] := Module[{v1,v2,v3,\[CurlyPhi]1,\[CurlyPhi]2,\[CurlyPhi]3,M1,M2,o\[Prime],o\[Prime]\[Prime]},
	v1 = o[[1;;2,3]];
	\[CurlyPhi]1 = -SafeArcTan @@ Reverse[v1];
	M1 = {{Cos[\[CurlyPhi]1],Sin[\[CurlyPhi]1],0},{-Sin[\[CurlyPhi]1],Cos[\[CurlyPhi]1],0},{0,0,1}};
	o\[Prime] = Simplify[M1 . o];
	v2 = o\[Prime][[2;;3,3]];
	\[CurlyPhi]2 = -SafeArcTan @@ Reverse[v2];
	M2 = {{1,0,0},{0,Cos[\[CurlyPhi]2],Sin[\[CurlyPhi]2]},{0,-Sin[\[CurlyPhi]2],Cos[\[CurlyPhi]2]}};
	o\[Prime]\[Prime] = Simplify[M2 . o\[Prime]];
	v3 = o\[Prime]\[Prime][[1;;2,1]];
	\[CurlyPhi]3 = SafeArcTan @@ v3;
	{\[CurlyPhi]1,\[CurlyPhi]2,\[CurlyPhi]3}
];

PauliDecomposition[U_?U4Q] := PauliDecomposition @ U4ToSO6 @ U
PauliDecomposition[U_?SO6Q] := Module[{a,b,c,d,U1a\[Prime],Da,U2a\[Prime],U1d\[Prime],Dd,U2d\[Prime],U1a,U2a,U1d,U2d,Db,Dc,blocks,\[CurlyPhi]XX,\[CurlyPhi]YY,\[CurlyPhi]ZZ,\[CurlyPhi]ZI1,\[CurlyPhi]XI1,\[CurlyPhi]ZI2,\[CurlyPhi]ZI3,\[CurlyPhi]XI2,\[CurlyPhi]ZI4,\[CurlyPhi]IZ1,\[CurlyPhi]IX1,\[CurlyPhi]IZ2,\[CurlyPhi]IZ3,\[CurlyPhi]IX2,\[CurlyPhi]IZ4},
	a = U[[1;;3,1;;3]];
	b = U[[4;;6,1;;3]];
	c = U[[1;;3,4;;6]];
	d = U[[4;;6,4;;6]];
	{U1a\[Prime],Da,U2a\[Prime]} = SingularValueDecomposition[a];
	{U1d\[Prime],Dd,U2d\[Prime]} = SingularValueDecomposition[d];
	{U1a,U2a,U1d,U2d} = {
		FullSimplify[Det[U1a\[Prime]]] * U1a\[Prime],
		FullSimplify[Det[U2a\[Prime]]] * U2a\[Prime],
		FullSimplify[Det[U1d\[Prime]]] * U1d\[Prime],
		FullSimplify[Det[U2d\[Prime]]] * U2d\[Prime]
	};
	Db = Simplify[Transpose[U1d] . b . U2a];
	Dc = Simplify[Transpose[U1a] . c . U2d];
	blocks = Map[{{Da[[#,#]],Db[[#,#]]},{Dc[[#,#]],Dd[[#,#]]}}&,Range[1,3]];
	{\[CurlyPhi]XX,\[CurlyPhi]YY,\[CurlyPhi]ZZ} = -Map[AngleFind,blocks];
	{\[CurlyPhi]ZI1,\[CurlyPhi]XI1,\[CurlyPhi]ZI2} = -SmallPauliDecomposition[U1a];
	{\[CurlyPhi]ZI3,\[CurlyPhi]XI2,\[CurlyPhi]ZI4} = -SmallPauliDecomposition[Transpose[U2a]];
	{\[CurlyPhi]IZ1,\[CurlyPhi]IX1,\[CurlyPhi]IZ2} = -SmallPauliDecomposition[U1d];
	{\[CurlyPhi]IZ3,\[CurlyPhi]IX2,\[CurlyPhi]IZ4} = -SmallPauliDecomposition[Transpose[U2d]];
	{
		{\[CurlyPhi]ZI1,"ZI"},{\[CurlyPhi]XI1,"XI"},{\[CurlyPhi]ZI2,"ZI"},
		{\[CurlyPhi]IZ1,"IZ"},{\[CurlyPhi]IX1,"IX"},{\[CurlyPhi]IZ2,"IZ"},
		{\[CurlyPhi]XX,"XX"},{\[CurlyPhi]YY,"YY"},{\[CurlyPhi]ZZ,"ZZ"},
		{\[CurlyPhi]ZI3,"ZI"},{\[CurlyPhi]XI2,"XI"},{\[CurlyPhi]ZI4,"ZI"},
		{\[CurlyPhi]IZ3,"IZ"},{\[CurlyPhi]IX2,"IX"},{\[CurlyPhi]IZ4,"IZ"}
	}
];
PauliDecomposition[_] := (
	Message[PauliDecomposition::notanoperator];
	$Failed
);


(* ::Subsubsection:: *)
(*Reconstruction using Pauli-rotation angles*)


ApproximateOp[_,\[Epsilon]_/;Not[RealQ[\[Epsilon]]]] := (
	Message[Approximate::notreal];
	$Failed
);
ApproximateOp[U_?U4Q,\[Epsilon]_] := Module[{anglesaxes,pauliapproximations},
	anglesaxes = PauliDecomposition[U];
	pauliapproximations = Map[PauliRotation[#[[1]],\[Epsilon]/15,#[[2]]]&,anglesaxes];
	Simplify[Dot @@ pauliapproximations]
];
ApproximateOp[U_?SO6Q,\[Epsilon]_] := Module[{anglesaxes,pauliapproximations},
	anglesaxes = PauliDecomposition[U];
	pauliapproximations = Map[PauliRotation[#[[1]],\[Epsilon]/15,#[[2]]]&,anglesaxes];
	U4ToSO6 @ Simplify @ Dot @@ pauliapproximations
];
ApproximateOp[_,_] := (
	Message[Approximate::notanoperator];
	$Failed
);

test[ApproximateSequence,"IfGaussianDoExact"] := BooleanQ;
test[ApproximateSequence,"SyllableType"] := MemberQ[{"Normal","Asymmetric"},#]&;
test[ApproximateSequence,"OutputType"] := MemberQ[{"String","HexDec"},#]&;

Options[ApproximateSequence] = {"IfGaussianDoExact" -> True,"OutputType" -> "String","SyllableType" -> "Normal"};

ApproximateSequence[_,\[Epsilon]_/;Not[RealQ[\[Epsilon]]],OptionsPattern[]] := (
	Message[ApproximateSequence::notreal];
	$Failed
);
ApproximateSequence[U_,\[Epsilon]_,opts:OptionsPattern[]] := Module[{normitops},
	normitops = FilterRules[{opts},Options[NormIt]];
	If[
		OptionValue["IfGaussianDoExact"] && GaussianQ[U],
			NormIt[U,Append[normitops,"UpToPhase"->False]],
			NormIt[ApproximateOp[U,\[Epsilon]],normitops]
	]
];


(* ::Subsection:: *)
(*End of functions in the private context*)


End[];


(* ::Section:: *)
(*Exported Functions*)


Protect[Id,X1,X2,Z1,Z2,W,H1,H2,S1,S2,CZ,CNOT12,CNOT21,EX,CS,U4ToSO6,
	IdSO6,X1SO6,X2SO6,Z1SO6,Z2SO6,\[CapitalIota]SO6,H1SO6,H2SO6,S1SO6,S2SO6,CZSO6,CNOT12SO6,CNOT21SO6,EXSO6,CSSO6,
	FromSequence,FromHexDec,CliffordQ,CliffordSynth,RightCliffordSimilar,LeftCliffordSimilar,GaussianQ,SyllableList,SyllableListAsymmetric,NormIt,FrobeniusDistance,
	CandidateFinder,PauliRotation,PauliRotationSequence,PauliDecomposition,ApproximateOp,ApproximateSequence];


EndPackage[];
