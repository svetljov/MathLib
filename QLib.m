(* ::Package:: *)

<<"D:\\SvetPAPERS\\MathLib\\ShortCommands.m"


(*Applies the operator F on qubit k out of ntotal*)
applyF[F_,k_,ntotal_]:=KP[Sequence@@CArr[ID[2],k-1],F,Sequence@@CArr[ID[2],ntotal-k]]

{id, \[Sigma]x,\[Sigma]y,\[Sigma]z}=PauliMatrix[#]&/@{0,1,2,3};
Clear[x,y,z]
{\[Sigma][x],\[Sigma][y],\[Sigma][z]}={\[Sigma]x,\[Sigma]y,\[Sigma]z}

(*Applies the Pauli matrix on qubit k out of ntotal*)
\[Sigma][a_,k_,ntotal_]:=applyF[\[Sigma][a],k,ntotal]

(*Hadamard gate*)
H=1/Sqrt[2] {{1,1},{1,-1}};

(*Rotation of a single qubit*)
R[\[Theta]_,\[Phi]_]:=MatrixExp[-I \[Theta]/2 (\[Sigma]x Cos[\[Phi]]+\[Sigma]y Sin[\[Phi]])]//FS;

(*Toffoli gate*)
Toff={{1,0,0,0,0,0,0,0},{0,1,0,0,0,0,0,0},{0,0,1,0,0,0,0,0},{0,0,0,1,0,0,0,0},{0,0,0,0,1,0,0,0},{0,0,0,0,0,1,0,0},{0,0,0,0,0,0,0,1},{0,0,0,0,0,0,1,0}};

(*The discrete Fourier transform*)
DFT[Nions_]:=Normalize[#]&/@Table[Exp[2 I \[Pi] (j-1) (k-1)/2^Nions],{j,2^Nions},{k,2^Nions}]//N;

(*Swaps qubits i and j out of ntotal*)
swap[i_,j_,ntotal_]:=
Module[{F=ConstantArray[0,{2^ntotal,2^ntotal}],el,set,mm},
set=basis[ntotal];
Do[
el=set[[nn]];
el[[{i,j}]]=el[[{j,i}]];
{{mm}}=Position[set,el];
F[[nn,mm]]=1,
{nn,Length[set]}];
F]

(*Reverse the qubit order - implements a swap gate for mirror-positioned qubits - 12345 => 54321*)
swapAllQubits[set_]:=
Module[{F=ConstantArray[0,{Length[set],Length[set]}],mm},
Do[
{{mm}}=Position[set,Reverse[set[[nn]]]];
F[[nn,mm]]=1,{nn,Length[set]}];F
]


(*QUANTUM STATES*)


(*Here we generate some target states*)

basis[nn_]:=PadLeft[#,nn]&/@Table[IntegerDigits[k,2],{k,0,2^nn-1}];(*generates the basis of internal spin states for nn ions*)
set[nn_,n_]:=Permutations@PadRight[ConstantArray[1,{n}],nn];(*set of states forming the Dicke manifold of nn ions with n excitations*)
pos[nn_,n_]:=Flatten[Position[basis[nn],#]&/@set[nn,n]];


Dicke[nn_,n_]:=Normalize[Table[If[MemberQ[pos[nn,n],k],1,0],{k,2^nn}]];(*Dicke state of nn ions with n excitations*)


(*original definition of Cluster states*)
ClusterO[n_]:=Module[{st},
If[n==3,st=MatrixExp[I \[Pi]/4 (KP[\[Sigma]z,\[Sigma]z,id]+KP[id,\[Sigma]z,\[Sigma]z])].Normalize[Flatten[KP[{1,1},{1,1},{1,1}]]]];
If[n==4,st=MatrixExp[I \[Pi]/4 (KP[\[Sigma]z,\[Sigma]z,id,id]+KP[id,\[Sigma]z,\[Sigma]z,id]+KP[id,id,\[Sigma]z,\[Sigma]z])].Normalize[Flatten[KP[{1,1},{1,1},{1,1},{1,1}]]]];
If[n==5,st=MatrixExp[I \[Pi]/4 (KP[\[Sigma]z,\[Sigma]z,id,id,id]+KP[id,\[Sigma]z,\[Sigma]z,id,id]+KP[id,id,\[Sigma]z,\[Sigma]z,id]+KP[id,id,id,\[Sigma]z,\[Sigma]z])].Normalize[Flatten[KP[{1,1},{1,1},{1,1},{1,1},{1,1}]]]];
If[n==6,st=MatrixExp[I \[Pi]/4 (KP[\[Sigma]z,\[Sigma]z,id,id,id,id]+KP[id,\[Sigma]z,\[Sigma]z,id,id,id]+KP[id,id,\[Sigma]z,\[Sigma]z,id,id]+KP[id,id,id,\[Sigma]z,\[Sigma]z,id]+KP[id,id,id,id,\[Sigma]z,\[Sigma]z])].Normalize[Flatten[KP[{1,1},{1,1},{1,1},{1,1},{1,1},{1,1}]]]];
st]

(*transformed Cluster states*)
Cluster[n_]:=Module[{st},
If[n==3,st={1,0,0,0,0,0,0,1}/Sqrt[2]];
If[n==4,st={1,0,0,1,0,0,0,0,0,0,0,0,1,0,0,-1}/2];
If[n==5,st={1,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,-1,-1,0}/(2 Sqrt[2])];
If[n==6,st={1,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,-1,-1,0,0,0,0,0,0,0,0,0,1,0,0,1,0,-1,-1,0,1,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0}/4];
If[n==7,st={1,0,0,-1,0,-1,1,0,1,0,0,-1,0,1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,1,0,1,-1,0,1,0,0,-1,0,1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,-1,0,-1,1,0,-1,0,0,1,0,-1,1,0,-1,0,0,1,0,1,-1,0,-1,0,0,1,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}/(4 Sqrt[2])];
st]
