(* ::Package:: *)

GeneratePositions[Nions_]:=Module[{r,V,xset,eq,sn,UpLimit,rootset,min},
Do[
r=Table[If[i>Ceiling[NN/2],-x[NN-i+1],x[i]],{i,NN}];$Assumptions={Flatten[r]\[Element]Reals};
V=\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(NN\)]
\*SuperscriptBox[\(r[\([i]\)]\), \(2\)]\)+\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(NN\)]\(
\*UnderoverscriptBox[\(\[Sum]\), \(j = i + 1\), \(NN\)]
\*FractionBox[\(2\), \(Abs[r[\([i]\)] - r[\([j]\)]]\)]\)\)/.Abs[a_]->a;
xset=Array[x,Ceiling[NN/2]];eq=Table[D[V,x[i]]==0,{i,Ceiling[NN/2]}];
If[NN<2,
	sn=NSolve[eq,xset],
	UpLimit=Ceiling[(NN-1)/2];
	rootset=Table[{x[n],sn[[n]]},{n,UpLimit}];min=2.018/NN^0.559/2;
	If[OddQ[NN],rootset=Append[rootset,{x[UpLimit+1],$MachineEpsilon}],rootset[[UpLimit]]={x[UpLimit],min}];
	sn=xset/.FindRoot[eq,rootset]
],{NN,Nions}];
r/.x[n_]:>sn[[n]]]


GenerateHessian[Nions_]:=Module[{u},
u=GeneratePositions[Nions];
Table[If[n==m,1+2\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(p = 1\), \(m - 1\)]
\*FractionBox[\(1\), 
SuperscriptBox[\(Abs[u[\([m]\)] - u[\([p]\)]]\), \(3\)]]\)+2\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(p = m + 1\), \(Nions\)]
\*FractionBox[\(1\), 
SuperscriptBox[\(Abs[u[\([m]\)] - u[\([p]\)]]\), \(3\)]]\),-(2/Abs[u[[m]]-u[[n]]]^3)],{n,Nions},{m,Nions}]]


GenerateJ[Nions_]:=Module[{Ainv},
Ainv=Inverse[GenerateHessian[Nions]];
Ainv=Ainv-DiagonalMatrix[Diagonal[Ainv]];
Ainv(*/Ainv[[1,2]]*)]


(*This reverses the qubit order - implements a swap gate for pairs of qubits - 12345 => 54321*)
GenerateSwapMatrix[set_]:=Module[
{F=ConstantArray[0,{Length[set],Length[set]}]},
Do[
{{mm}}=Position[set,Reverse[set[[nn]]]];
F[[nn,mm]]=1,{nn,Length[set]}];F
]


DFT[Nions_]:=Normalize[#]&/@Table[Exp[2 I \[Pi] (j-1) (k-1)/2^Nions],{j,2^Nions},{k,2^Nions}]//N;


{\[Sigma]x,\[Sigma]y,\[Sigma]z,id}={PauliMatrix[1],PauliMatrix[2],PauliMatrix[3],PauliMatrix[4]};

(*rotation of a single spin*)
inv=U[\[Pi],0];
U[A_,\[Phi]_]=MatrixExp[-I A/2 (\[Sigma]x Cos[\[Phi]]+\[Sigma]y Sin[\[Phi]])]//FS;

(*spin-spin coupling interaction*)
\[CapitalPhi]2[T_]:=MatrixExp[I T/2 (J[[1,2]] KP[\[Sigma]z,\[Sigma]z])];
\[CapitalPhi]3[T_]:=MatrixExp[I T/2 (J[[1,2]] KP[\[Sigma]z,\[Sigma]z,id] + J[[1,3]] KP[\[Sigma]z,id,\[Sigma]z] + J[[2,3]] KP[id,\[Sigma]z,\[Sigma]z])];
\[CapitalPhi]4[T_]:=MatrixExp[I T/2 (J[[1,2]] KP[\[Sigma]z,\[Sigma]z,id,id] + J[[1,3]] KP[\[Sigma]z,id,\[Sigma]z,id] + J[[1,4]] KP[\[Sigma]z,id,id,\[Sigma]z] + J[[2,3]] KP[id,\[Sigma]z,\[Sigma]z,id] + J[[2,4]] KP[id,\[Sigma]z,id,\[Sigma]z] + J[[3,4]] KP[id,id,\[Sigma]z,\[Sigma]z])];
\[CapitalPhi]5[T_]:=MatrixExp[I T/2 (J[[1,2]] KP[\[Sigma]z,\[Sigma]z,id,id,id] + J[[1,3]] KP[\[Sigma]z,id,\[Sigma]z,id,id] + J[[1,4]] KP[\[Sigma]z,id,id,\[Sigma]z,id] + J[[1,5]] KP[\[Sigma]z,id,id,id,\[Sigma]z] + J[[2,3]] KP[id,\[Sigma]z,\[Sigma]z,id,id] + J[[2,4]] KP[id,\[Sigma]z,id,\[Sigma]z,id] + J[[2,5]] KP[id,\[Sigma]z,id,id,\[Sigma]z] + J[[3,4]] KP[id,id,\[Sigma]z,\[Sigma]z,id] + J[[3,5]] KP[id,id,\[Sigma]z,id,\[Sigma]z] + J[[4,5]] KP[id,id,id,\[Sigma]z,\[Sigma]z])];
\[CapitalPhi]6[T_]:=MatrixExp[I T/2 (J[[1,2]] KP[\[Sigma]z,\[Sigma]z,id,id,id,id] + J[[1,3]] KP[\[Sigma]z,id,\[Sigma]z,id,id,id] + J[[1,4]] KP[\[Sigma]z,id,id,\[Sigma]z,id,id] + J[[1,5]] KP[\[Sigma]z,id,id,id,\[Sigma]z,id] + J[[1,6]] KP[\[Sigma]z,id,id,id,id,\[Sigma]z] + J[[2,3]] KP[id,\[Sigma]z,\[Sigma]z,id,id,id] + J[[2,4]] KP[id,\[Sigma]z,id,\[Sigma]z,id,id] + J[[2,5]] KP[id,\[Sigma]z,id,id,\[Sigma]z,id] + J[[2,6]] KP[id,\[Sigma]z,id,id,id,\[Sigma]z] + J[[3,4]] KP[id,id,\[Sigma]z,\[Sigma]z,id,id] + J[[3,5]] KP[id,id,\[Sigma]z,id,\[Sigma]z,id] + J[[3,6]] KP[id,id,\[Sigma]z,id,id,\[Sigma]z] + J[[4,5]] KP[id,id,id,\[Sigma]z,\[Sigma]z,id] + J[[4,6]] KP[id,id,id,\[Sigma]z,id,\[Sigma]z] + J[[5,6]] KP[id,id,id,id,\[Sigma]z,\[Sigma]z])];


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