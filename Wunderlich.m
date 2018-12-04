(* ::Package:: *)

<<"D:\\SvetPAPERS\\MathLib\\ShortCommands.m"


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


(*spin-spin coupling interaction*)
\[CapitalPhi]2[T_]:=MatrixExp[I T/2 (J[[1,2]] KP[\[Sigma]z,\[Sigma]z])];
\[CapitalPhi]3[T_]:=MatrixExp[I T/2 (J[[1,2]] KP[\[Sigma]z,\[Sigma]z,id] + J[[1,3]] KP[\[Sigma]z,id,\[Sigma]z] + J[[2,3]] KP[id,\[Sigma]z,\[Sigma]z])];
\[CapitalPhi]4[T_]:=MatrixExp[I T/2 (J[[1,2]] KP[\[Sigma]z,\[Sigma]z,id,id] + J[[1,3]] KP[\[Sigma]z,id,\[Sigma]z,id] + J[[1,4]] KP[\[Sigma]z,id,id,\[Sigma]z] + J[[2,3]] KP[id,\[Sigma]z,\[Sigma]z,id] + J[[2,4]] KP[id,\[Sigma]z,id,\[Sigma]z] + J[[3,4]] KP[id,id,\[Sigma]z,\[Sigma]z])];
\[CapitalPhi]5[T_]:=MatrixExp[I T/2 (J[[1,2]] KP[\[Sigma]z,\[Sigma]z,id,id,id] + J[[1,3]] KP[\[Sigma]z,id,\[Sigma]z,id,id] + J[[1,4]] KP[\[Sigma]z,id,id,\[Sigma]z,id] + J[[1,5]] KP[\[Sigma]z,id,id,id,\[Sigma]z] + J[[2,3]] KP[id,\[Sigma]z,\[Sigma]z,id,id] + J[[2,4]] KP[id,\[Sigma]z,id,\[Sigma]z,id] + J[[2,5]] KP[id,\[Sigma]z,id,id,\[Sigma]z] + J[[3,4]] KP[id,id,\[Sigma]z,\[Sigma]z,id] + J[[3,5]] KP[id,id,\[Sigma]z,id,\[Sigma]z] + J[[4,5]] KP[id,id,id,\[Sigma]z,\[Sigma]z])];
\[CapitalPhi]6[T_]:=MatrixExp[I T/2 (J[[1,2]] KP[\[Sigma]z,\[Sigma]z,id,id,id,id] + J[[1,3]] KP[\[Sigma]z,id,\[Sigma]z,id,id,id] + J[[1,4]] KP[\[Sigma]z,id,id,\[Sigma]z,id,id] + J[[1,5]] KP[\[Sigma]z,id,id,id,\[Sigma]z,id] + J[[1,6]] KP[\[Sigma]z,id,id,id,id,\[Sigma]z] + J[[2,3]] KP[id,\[Sigma]z,\[Sigma]z,id,id,id] + J[[2,4]] KP[id,\[Sigma]z,id,\[Sigma]z,id,id] + J[[2,5]] KP[id,\[Sigma]z,id,id,\[Sigma]z,id] + J[[2,6]] KP[id,\[Sigma]z,id,id,id,\[Sigma]z] + J[[3,4]] KP[id,id,\[Sigma]z,\[Sigma]z,id,id] + J[[3,5]] KP[id,id,\[Sigma]z,id,\[Sigma]z,id] + J[[3,6]] KP[id,id,\[Sigma]z,id,id,\[Sigma]z] + J[[4,5]] KP[id,id,id,\[Sigma]z,\[Sigma]z,id] + J[[4,6]] KP[id,id,id,\[Sigma]z,id,\[Sigma]z] + J[[5,6]] KP[id,id,id,id,\[Sigma]z,\[Sigma]z])];
