#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Input::Initialization:: *)
vC1c=-( ((16 E^(-4 z))/(c^4 (E^((c E^z)/(4 \[Epsilon]))-S))-(24 E^(-4 z))/(c^4 (E^((c E^z)/(2 \[Epsilon]))-S))+(16 E^(-4 z))/(c^4 (E^((3 c E^z)/(4 \[Epsilon]))-S))-(4 E^(-4 z))/(c^4 (E^((c E^z)/\[Epsilon])-S))+(4 E^(-4 z))/(c^4 (-1+S))+E^(-3 z+(c E^z)/(4 \[Epsilon]))/(c^3 (E^((c E^z)/(4 \[Epsilon]))-S)^2 \[Epsilon])-(3 E^(-3 z+(c E^z)/(2 \[Epsilon])))/(c^3 (E^((c E^z)/(2 \[Epsilon]))-S)^2 \[Epsilon])+(3 E^(-3 z+(3 c E^z)/(4 \[Epsilon])))/(c^3 (E^((3 c E^z)/(4 \[Epsilon]))-S)^2 \[Epsilon])-E^(-3 z+(c E^z)/\[Epsilon])/(c^3 (E^((c E^z)/\[Epsilon])-S)^2 \[Epsilon]))/(E^(-4 z)/(c^4 (1-S))-(4 E^(-4 z))/(c^4 (E^((c E^z)/(4 \[Epsilon]))-S))+(6 E^(-4 z))/(c^4 (E^((c E^z)/(2 \[Epsilon]))-S))-(4 E^(-4 z))/(c^4 (E^((3 c E^z)/(4 \[Epsilon]))-S))+E^(-4 z)/(c^4 (E^((c E^z)/\[Epsilon])-S))));
vC2c=(-((24 E^(-4 z))/(c^4 (E^((c E^z)/(4 \[Epsilon]))-S)))+(36 E^(-4 z))/(c^4 (E^((c E^z)/(2 \[Epsilon]))-S))-(24 E^(-4 z))/(c^4 (E^((3 c E^z)/(4 \[Epsilon]))-S))+(6 E^(-4 z))/(c^4 (E^((c E^z)/\[Epsilon])-S))-(6 E^(-4 z))/(c^4 (-1+S))-E^(-2 z+(c E^z)/(2 \[Epsilon]))/(8 c^2 (E^((c E^z)/(4 \[Epsilon]))-S)^3 \[Epsilon]^2)+(3 E^(-2 z+(c E^z)/\[Epsilon]))/(4 c^2 (E^((c E^z)/(2 \[Epsilon]))-S)^3 \[Epsilon]^2)-(9 E^(-2 z+(3 c E^z)/(2 \[Epsilon])))/(8 c^2 (E^((3 c E^z)/(4 \[Epsilon]))-S)^3 \[Epsilon]^2)+E^(-2 z+(2 c E^z)/\[Epsilon])/(2 c^2 (E^((c E^z)/\[Epsilon])-S)^3 \[Epsilon]^2)-(E^(-2 z+(c E^z)/(4 \[Epsilon])) S)/(8 c^2 (E^((c E^z)/(4 \[Epsilon]))-S)^3 \[Epsilon]^2)+(3 E^(-2 z+(c E^z)/(2 \[Epsilon])) S)/(4 c^2 (E^((c E^z)/(2 \[Epsilon]))-S)^3 \[Epsilon]^2)-(9 E^(-2 z+(3 c E^z)/(4 \[Epsilon])) S)/(8 c^2 (E^((3 c E^z)/(4 \[Epsilon]))-S)^3 \[Epsilon]^2)+(E^(-2 z+(c E^z)/\[Epsilon]) S)/(2 c^2 (E^((c E^z)/\[Epsilon])-S)^3 \[Epsilon]^2)+E^(-3 z+(c E^z)/(2 \[Epsilon]))/(c^3 (E^((c E^z)/(4 \[Epsilon]))-S)^3 \[Epsilon])-(4 E^(-3 z+(c E^z)/(4 \[Epsilon])))/(c^3 (E^((c E^z)/(4 \[Epsilon]))-S)^2 \[Epsilon])-(3 E^(-3 z+(c E^z)/\[Epsilon]))/(c^3 (E^((c E^z)/(2 \[Epsilon]))-S)^3 \[Epsilon])+(12 E^(-3 z+(c E^z)/(2 \[Epsilon])))/(c^3 (E^((c E^z)/(2 \[Epsilon]))-S)^2 \[Epsilon])+(3 E^(-3 z+(3 c E^z)/(2 \[Epsilon])))/(c^3 (E^((3 c E^z)/(4 \[Epsilon]))-S)^3 \[Epsilon])-(12 E^(-3 z+(3 c E^z)/(4 \[Epsilon])))/(c^3 (E^((3 c E^z)/(4 \[Epsilon]))-S)^2 \[Epsilon])-E^(-3 z+(2 c E^z)/\[Epsilon])/(c^3 (E^((c E^z)/\[Epsilon])-S)^3 \[Epsilon])+(4 E^(-3 z+(c E^z)/\[Epsilon]))/(c^3 (E^((c E^z)/\[Epsilon])-S)^2 \[Epsilon])-(E^(-3 z+(c E^z)/(4 \[Epsilon])) S)/(c^3 (E^((c E^z)/(4 \[Epsilon]))-S)^3 \[Epsilon])+(3 E^(-3 z+(c E^z)/(2 \[Epsilon])) S)/(c^3 (E^((c E^z)/(2 \[Epsilon]))-S)^3 \[Epsilon])-(3 E^(-3 z+(3 c E^z)/(4 \[Epsilon])) S)/(c^3 (E^((3 c E^z)/(4 \[Epsilon]))-S)^3 \[Epsilon])+(E^(-3 z+(c E^z)/\[Epsilon]) S)/(c^3 (E^((c E^z)/\[Epsilon])-S)^3 \[Epsilon]))/(E^(-4 z)/(c^4 (1-S))-(4 E^(-4 z))/(c^4 (E^((c E^z)/(4 \[Epsilon]))-S))+(6 E^(-4 z))/(c^4 (E^((c E^z)/(2 \[Epsilon]))-S))-(4 E^(-4 z))/(c^4 (E^((3 c E^z)/(4 \[Epsilon]))-S))+E^(-4 z)/(c^4 (E^((c E^z)/\[Epsilon])-S))) ;
vC2c2=1/(32 \[Epsilon]^2) ((4 c^2 E^(2 z+(c E^z)/(2 \[Epsilon])) (E^((c E^z)/(2 \[Epsilon]))+S))/(E^((c E^z)/(2 \[Epsilon]))-S)^2+(9 c^2 E^(2 z+(3 c E^z)/(4 \[Epsilon])) (E^((3 c E^z)/(4 \[Epsilon]))+S))/(E^((3 c E^z)/(4 \[Epsilon]))-S)^2+(16 c^2 E^(2 z+(c E^z)/\[Epsilon]) (E^((c E^z)/\[Epsilon])+S))/(E^((c E^z)/\[Epsilon])-S)^2+(c^2 E^(2 z+(c E^z)/(4 \[Epsilon])) (36 E^((5 c E^z)/(4 \[Epsilon]))+80 E^((3 c E^z)/(4 \[Epsilon])) S+75 E^((c E^z)/\[Epsilon]) S+3 S^2+20 E^((c E^z)/(4 \[Epsilon])) S^2+27 E^((c E^z)/(2 \[Epsilon])) S (1+S)))/(E^((3 c E^z)/(2 \[Epsilon]))+5 E^((c E^z)/\[Epsilon]) S+3 E^((5 c E^z)/(4 \[Epsilon])) S+3 E^((c E^z)/(4 \[Epsilon])) S^2+5 E^((c E^z)/(2 \[Epsilon])) S^2+S^3+3 E^((3 c E^z)/(4 \[Epsilon])) S (1+S))-(6 c^2 E^(2 z+(c E^z)/\[Epsilon]) (E^((3 c E^z)/(4 \[Epsilon]))+E^((c E^z)/(2 \[Epsilon])) (3-2 S)-5 E^((c E^z)/(4 \[Epsilon])) S+S (-1+4 S)))/((-1+E^((c E^z)/(4 \[Epsilon]))) (E^((c E^z)/(4 \[Epsilon]))-S) (E^((c E^z)/(2 \[Epsilon]))-S) (E^((3 c E^z)/(4 \[Epsilon]))-S))+(c^2 E^(2 z+(c E^z)/(4 \[Epsilon])) (9 E^((3 c E^z)/(4 \[Epsilon]))+E^((c E^z)/(2 \[Epsilon])) (2-23 S)+S-4 S^2+E^((c E^z)/(4 \[Epsilon])) (1-2 S+16 S^2)))/((-1+E^((c E^z)/(4 \[Epsilon])))^2 (E^((c E^z)/(4 \[Epsilon]))-S)^2)-(32 c E^(z+(c E^z)/(4 \[Epsilon])) (1+3 E^((c E^z)/(4 \[Epsilon]))-4 S) \[Epsilon])/((-1+E^((c E^z)/(4 \[Epsilon]))) (E^((c E^z)/(4 \[Epsilon]))-S))+(96 c E^(z+(3 c E^z)/(4 \[Epsilon])) \[Epsilon])/(E^((3 c E^z)/(4 \[Epsilon]))-S)+320 \[Epsilon]^2+(4 c E^(z+(c E^z)/(2 \[Epsilon])) (-((c E^(z+(c E^z)/(4 \[Epsilon])) (1+3 E^((c E^z)/(4 \[Epsilon]))-4 S))/((-1+E^((c E^z)/(4 \[Epsilon]))) (E^((c E^z)/(4 \[Epsilon]))-S)))+16 \[Epsilon]))/(E^((c E^z)/(2 \[Epsilon]))-S)+(8 c E^(z+(c E^z)/\[Epsilon]) (c E^(z+(c E^z)/(4 \[Epsilon])) (2 E^((3 c E^z)/(2 \[Epsilon]))+5 E^((c E^z)/\[Epsilon]) S-E^((3 c E^z)/(4 \[Epsilon])) (-5+S) S-5 E^((c E^z)/(4 \[Epsilon])) S^2-E^((5 c E^z)/(4 \[Epsilon])) (6+S)+S^2 (-1+4 S)-E^((c E^z)/(2 \[Epsilon])) S (-3+5 S))+16 (E^((7 c E^z)/(4 \[Epsilon]))+E^((3 c E^z)/(4 \[Epsilon])) S+E^((c E^z)/\[Epsilon]) S^2+S^3-E^((3 c E^z)/(2 \[Epsilon])) (1+S)-E^((c E^z)/(4 \[Epsilon])) S^2 (1+S)) \[Epsilon]))/((-1+E^((c E^z)/(4 \[Epsilon]))) (E^((c E^z)/(4 \[Epsilon]))-S) (E^((c E^z)/(2 \[Epsilon]))-S) (E^((3 c E^z)/(4 \[Epsilon]))-S) (E^((c E^z)/\[Epsilon])-S))-(2 c E^(z+(c E^z)/(4 \[Epsilon])) (6 E^((5 c E^z)/(4 \[Epsilon]))+20 E^((3 c E^z)/(4 \[Epsilon])) S+15 E^((c E^z)/\[Epsilon]) S+3 S^2+10 E^((c E^z)/(4 \[Epsilon])) S^2+9 E^((c E^z)/(2 \[Epsilon])) S (1+S)) (4 c E^(z+(c E^z)/\[Epsilon]) (-1+E^((c E^z)/(4 \[Epsilon]))) (E^((c E^z)/(4 \[Epsilon]))-S) (E^((c E^z)/(2 \[Epsilon]))-S) (E^((3 c E^z)/(4 \[Epsilon]))-S)+(E^((c E^z)/\[Epsilon])-S) (c E^(z+(c E^z)/(4 \[Epsilon])) (2 E^((3 c E^z)/(2 \[Epsilon]))+5 E^((c E^z)/\[Epsilon]) S-E^((3 c E^z)/(4 \[Epsilon])) (-5+S) S-5 E^((c E^z)/(4 \[Epsilon])) S^2-E^((5 c E^z)/(4 \[Epsilon])) (6+S)+S^2 (-1+4 S)-E^((c E^z)/(2 \[Epsilon])) S (-3+5 S))+16 (E^((7 c E^z)/(4 \[Epsilon]))+E^((3 c E^z)/(4 \[Epsilon])) S+E^((c E^z)/\[Epsilon]) S^2+S^3-E^((3 c E^z)/(2 \[Epsilon])) (1+S)-E^((c E^z)/(4 \[Epsilon])) S^2 (1+S)) \[Epsilon])))/((-1+E^((c E^z)/(4 \[Epsilon]))) (E^((c E^z)/(4 \[Epsilon]))-S) (E^((c E^z)/(2 \[Epsilon]))-S) (E^((3 c E^z)/(4 \[Epsilon]))-S) (E^((c E^z)/\[Epsilon])-S) (E^((3 c E^z)/(2 \[Epsilon]))+5 E^((c E^z)/\[Epsilon]) S+3 E^((5 c E^z)/(4 \[Epsilon])) S+3 E^((c E^z)/(4 \[Epsilon])) S^2+5 E^((c E^z)/(2 \[Epsilon])) S^2+S^3+3 E^((3 c E^z)/(4 \[Epsilon])) S (1+S))));
Wc=1/c^4 (1/(1-S )-4/(E^(c/(4 \[Epsilon]))-S)+6/(E^(c/(2 \[Epsilon]))-S)-4/(E^((3 c)/(4 \[Epsilon]))-S)+1/(E^(c/\[Epsilon])-S))/.c->c w;
We=1/c^4 (1/(1-S )-4/(E^(c/(4 \[Epsilon]))-S)+6/(E^(c/(2 \[Epsilon]))-S)-4/(E^((3 c)/(4 \[Epsilon]))-S)+1/(E^(c/\[Epsilon])-S))/.\[Epsilon]->\[Epsilon] w;
Ws=w/c^4 (1/(1-S )-4/(E^(c/(4 \[Epsilon]))-S)+6/(E^(c/(2 \[Epsilon]))-S)-4/(E^((3 c)/(4 \[Epsilon]))-S)+1/(E^(c/\[Epsilon])-S))/.S->S w;
vT1e=(1/(1-S)-4/(E^((c E^-z)/(4 \[Epsilon]))-S)+6/(E^((c E^-z)/(2 \[Epsilon]))-S)-4/(E^((3 c E^-z)/(4 \[Epsilon]))-S)+1/(E^((c E^-z)/\[Epsilon])-S))/c^4;
vT2e=-((-((c E^(-z+(c E^-z)/(4 \[Epsilon])))/((E^((c E^-z)/(4 \[Epsilon]))-S)^2 \[Epsilon]))+(3 c E^(-z+(c E^-z)/(2 \[Epsilon])))/((E^((c E^-z)/(2 \[Epsilon]))-S)^2 \[Epsilon])-(3 c E^(-z+(3 c E^-z)/(4 \[Epsilon])))/((E^((3 c E^-z)/(4 \[Epsilon]))-S)^2 \[Epsilon])+(c E^(-z+(c E^-z)/\[Epsilon]))/((E^((c E^-z)/\[Epsilon])-S)^2 \[Epsilon]))/c^4);
vT3e=-(1/c^4)( (c  E^(-2 z+(c E^-z)/(4 \[Epsilon])) (-c E^((c E^-z)/(4 \[Epsilon]))-c S+4 E^(z+(c E^-z)/(4 \[Epsilon])) \[Epsilon]-4 E^z S \[Epsilon]))/(8 (E^((c E^-z)/(4 \[Epsilon]))-S)^3 \[Epsilon]^2)+(3 c E^(-2 z+(3 c E^-z)/(4 \[Epsilon])) (-3 c E^((3 c E^-z)/(4 \[Epsilon]))-3 c S+4 E^(z+(3 c E^-z)/(4 \[Epsilon])) \[Epsilon]-4 E^z S \[Epsilon]))/(8 (E^((3 c E^-z)/(4 \[Epsilon]))-S)^3 \[Epsilon]^2)-(3 c  E^(-2 z+(c E^-z)/(2 \[Epsilon])) (-c E^((c E^-z)/(2 \[Epsilon]))-c S+2 E^(z+(c E^-z)/(2 \[Epsilon])) \[Epsilon]-2 E^z S \[Epsilon]))/(4 (E^((c E^-z)/(2 \[Epsilon]))-S)^3 \[Epsilon]^2)-(c  E^(-2 z+(c E^-z)/\[Epsilon]) (-c E^((c E^-z)/\[Epsilon])-c S+E^(z+(c E^-z)/\[Epsilon]) \[Epsilon]-E^z S \[Epsilon]))/(2 (E^((c E^-z)/\[Epsilon])-S)^3 \[Epsilon]^2));
vT1s=(E^z (1/(1-E^z S)-4/(E^(c/(4 \[Epsilon]))-E^z S)+6/(E^(c/(2 \[Epsilon]))-E^z S)-4/(E^((3 c)/(4 \[Epsilon]))-E^z S)+1/(E^(c/\[Epsilon])-E^z S)))/c^4;
vT2s=(E^z (-(1/(1-E^z S))+(4 E^z S)/(E^(c/(4 \[Epsilon]))-E^z S)^2+4/(E^(c/(4 \[Epsilon]))-E^z S)-(6 E^z S)/(E^(c/(2 \[Epsilon]))-E^z S)^2-6/(E^(c/(2 \[Epsilon]))-E^z S)+(4 E^z S)/(E^((3 c)/(4 \[Epsilon]))-E^z S)^2+4/(E^((3 c)/(4 \[Epsilon]))-E^z S)-(E^z S)/(E^(c/\[Epsilon])-E^z S)^2-1/(E^(c/\[Epsilon])-E^z S)-(E^z S)/(-1+E^z S)^2))/c^4;
vT3s=1/c^4 E^z ((2 E^z S (E^(c/(4 \[Epsilon]))-3 E^z S))/(E^(c/(4 \[Epsilon]))-E^z S)^3-(4 E^z S)/(E^(c/(4 \[Epsilon]))-E^z S)^2-(3 E^z S (E^(c/(2 \[Epsilon]))-3 E^z S))/(E^(c/(2 \[Epsilon]))-E^z S)^3+(6 E^z S)/(E^(c/(2 \[Epsilon]))-E^z S)^2+(2 E^z S (E^((3 c)/(4 \[Epsilon]))-3 E^z S))/(E^((3 c)/(4 \[Epsilon]))-E^z S)^3-(4 E^z S)/(E^((3 c)/(4 \[Epsilon]))-E^z S)^2-(E^z S (E^(c/\[Epsilon])-3 E^z S))/(2 (E^(c/\[Epsilon])-E^z S)^3)+(E^z S)/(E^(c/\[Epsilon])-E^z S)^2+(E^z S)/(-1+E^z S)^2-(E^z S (-1+3 E^z S))/(2 (-1+E^z S)^3)+1/2 (-(1/(1-E^z S))+4/(E^(c/(4 \[Epsilon]))-E^z S)-6/(E^(c/(2 \[Epsilon]))-E^z S)+4/(E^((3 c)/(4 \[Epsilon]))-E^z S)-1/(E^(c/\[Epsilon])-E^z S)));

vC1e=vT2e/vT1e;
vC2e=-(vT3e/vT1e);
vC2e2=(c E^(-2 z+(c E^-z)/(4 \[Epsilon])) ((-c (E^((c E^-z)/(4 \[Epsilon]))+S)+8 E^z (E^((c E^-z)/(4 \[Epsilon]))-S) \[Epsilon])/(E^((c E^-z)/(4 \[Epsilon]))-S)^3+(6 E^((c E^-z)/(4 \[Epsilon])) (c (E^((c E^-z)/(2 \[Epsilon]))+S)-4 E^z (E^((c E^-z)/(2 \[Epsilon]))-S) \[Epsilon]))/(E^((c E^-z)/(2 \[Epsilon]))-S)^3-(3 E^((c E^-z)/(2 \[Epsilon])) (3 c (E^((3 c E^-z)/(4 \[Epsilon]))+S)-8 E^z (E^((3 c E^-z)/(4 \[Epsilon]))-S) \[Epsilon]))/(E^((3 c E^-z)/(4 \[Epsilon]))-S)^3+(4 E^((3 c E^-z)/(4 \[Epsilon])) (c (E^((c E^-z)/\[Epsilon])+S)-2 E^z (E^((c E^-z)/\[Epsilon])-S) \[Epsilon]))/(E^((c E^-z)/\[Epsilon])-S)^3))/( 8 (1/(1-S)-4/(E^((c E^-z)/(4 \[Epsilon]))-S)+6/(E^((c E^-z)/(2 \[Epsilon]))-S)-4/(E^((3 c E^-z)/(4 \[Epsilon]))-S)+1/(E^((c E^-z)/\[Epsilon])-S)) \[Epsilon]^2);

vC1s=-(vT2s/vT1s);
vC2s=vT3s/vT1s;
vC2s2=(E^z S (-(4/(E^(c/(4 \[Epsilon]))-E^z S)^2)+6/(E^(c/(2 \[Epsilon]))-E^z S)^2-4/(E^((3 c)/(4 \[Epsilon]))-E^z S)^2+1/(E^(c/\[Epsilon])-E^z S)^2+1/(-1+E^z S)^2+E^z S (1/(1-E^z S)^3-4/(E^(c/(4 \[Epsilon]))-E^z S)^3+6/(E^(c/(2 \[Epsilon]))-E^z S)^3-4/(E^((3 c)/(4 \[Epsilon]))-E^z S)^3+1/(E^(c/\[Epsilon])-E^z S)^3)))/(1/(1-E^z S)-4/(E^(c/(4 \[Epsilon]))-E^z S)+6/(E^(c/(2 \[Epsilon]))-E^z S)-4/(E^((3 c)/(4 \[Epsilon]))-E^z S)+1/(E^(c/\[Epsilon])-E^z S));

Growth=E^(-(c/\[Epsilon])(t+1))/c^4 (E^(c/(4 \[Epsilon]) (t+1))-1)^4;

(* Selection *)

pis = pi-Subscript[s, dir]h pi-Subscript[s, dir](1-h)dii-Subscript[s, dir](1-2h)dijj+Subscript[s, int](1-h)(1-2h)(dijij-dii djj);
pjs = pj-Subscript[s, dir]h pj-Subscript[s, dir](1-h)djj-Subscript[s, dir](1-2h)dijj+Subscript[s, int](1-h)(1-2h)(dijij-dii djj);

diis = dii-Subscript[s, dir](1-2h)(dijij-dii djj);
djjs = djj-Subscript[s, dir](1-2h)(dijij-dii djj);

dijjs = dijj-Subscript[s, dir](1-h)(dijij-dii djj);
dijijs= dijij;

(* Reproduction *)

pir=pis + u;

diir=F/(1+F) (diis+pis);
djjr=F/(1+F) (djjs+pjs);
dijjr=F/(1+F) dijjs;
dijijr=F/(1+F) ((1-R)(dijijs + pis pjs) + R(pis djjs + pjs diis));
\[CapitalDelta]dijij=Normal[Series[dijijr-dijij/.dijj-> dijj sm/.dijij-> dijij/.Subscript[s, dir]-> Subscript[s, dir]sm/.Subscript[s, int]->Subscript[s, int] sm^2/.dii-> F pi/.djj-> F pj,{sm,0,1}]]/.sm->1;
solDijij=Solve[\[CapitalDelta]dijij==0,dijij]//Simplify;
tempDijij=dijij/.solDijij[[1]];
temp2Dijij=Normal[Series[tempDijij/.pi->pi sm/.pj-> pj sm,{sm,0,2}]]/.sm->1//FullSimplify;
qdijij=Normal[Series[temp2Dijij,{Subscript[s, dir],0,1}]];
neutdijij=qdijij/.Subscript[s, dir]->0;

\[CapitalDelta]dii=Normal[Series[diir-dii/.Subscript[s, dir]->Subscript[s, dir]sm/.Subscript[s, int]->Subscript[s, int] sm^2/.dijj-> dijj sm/.dijij->neutdijij/.djj-> F pj,{sm,0,1}]]/.sm->1;
\[CapitalDelta]djj=Normal[Series[djjr-djj/.Subscript[s, dir]->Subscript[s, dir]sm/.Subscript[s, int]->Subscript[s, int]/.dijj-> dijj sm/.dijij->neutdijij/.dii-> F pi,{sm,0,1}]]/.sm->1;

sDii=Solve[\[CapitalDelta]dii==0,dii]//Simplify;
sDjj=Solve[\[CapitalDelta]djj==0,djj]//Simplify;

tempDii=dii/.sDii[[1]];
tempDjj=djj/.sDjj[[1]];

t1=Normal[Series[tempDii,{Subscript[s, dir],0,1}]];
t2=Normal[Series[t1,{pj,0,1}]];
qdii=t2;

t3=Normal[Series[tempDjj,{Subscript[s, dir],0,1}]];
t4=Normal[Series[t3,{pj,0,1}]];
qdjj=t4;

q1dii=qdii/.pj-> 0;
q1djj=qdjj/.pi-> 0;

\[CapitalDelta]dijj=Normal[Series[dijjr-dijj/.Subscript[s, dir]-> Subscript[s, dir]sm/.Subscript[s, int]-> Subscript[s, int] sm^2/.dijij->neutdijij/.dii-> F pi/.djj-> F pj/.dijj-> dijj sm,{sm,0,1}]]/.sm->1;
sDijj=Solve[\[CapitalDelta]dijj==0,dijj]//Simplify;
qdijj=Normal[Series[dijj/.sDijj[[1]],{Subscript[s, dir],0,1}]];

qle={dii->qdii,djj->qdjj,dijij->qdijij,dijj->qdijj};

dpi=Normal[Series[pir-pi/.qle/.Subscript[s, dir]-> Subscript[s, dir] sm/.Subscript[s, int]-> Subscript[s, int] sm^2,{sm,0,2}]]/.sm->1;
dpi2=Normal[Series[dpi/.pi-> pi sm/.pj-> pj sm,{sm,0,2}]]/.sm->1;

qself = {dii-> 1/2 (pi+dii),djj->1/2 (pj+djj),dijj-> 1/2 dijj,dijij-> 1/4 (pi pj + pj dii+pi djj+dijij)};
qout={dii-> 0,djj->0,dijj->0,dijij->0};

notID = {z-> 0,pj->pi};

qle={dii->qdii,djj->qdjj,dijij->qdijij,dijj->qdijj};
sqle={dii->F pi,djj->F pj,dijij->neutdijij,dijj->0,dij-> 0 ,divj-> 0};

\[CapitalDelta]piNoInt=pir-pi/.qle/.pj->0;
solNoInt=Solve[\[CapitalDelta]piNoInt==0,pi];
pNI=pi/.solNoInt[[1]];
solp=Solve[dpi2==0,pi];
peq=pi/.solp[[1]];
peq2=Normal[Series[peq,{pj,0,1}]];
peq3=Normal[Series[peq2/.z->0/.Subscript[s, dir]-> Subscript[s, dir] sm/.Subscript[s, int]-> Subscript[s, int] sm^2/.u-> u sm^2,{sm,0,2}]]/.sm->1;

(*Importing parameters table and removing names *)

pf=Import["math_par_lf.txt","Table"];
For[i=1,i<= Length[pf],i++,
pf[[i]]=Delete[pf[[i]],1];
];

cL= pf[[1]];
SL=pf[[2]];
eL=pf[[3]];
aL=Range[0,1,0.05];
uL=pf[[4]];
hL=pf[[5]];
selL=pf[[6]];
modelL={"c","e","s"};
imax=pf[[7]][[1]];
step=pf[[8]][[1]];
resume=pf[[9]][[1]];
lcorr=pf[[10]];

timetot=Length[cL]*Length[SL]*Length[eL]*Length[aL]*Length[uL]*Length[hL]*Length[selL]*Length[modelL]*Length[lcorr];

wS=E^(-2 s h pi-s(1-2h)F pi);
wG = E^(2 s h pi+s(1-2h)F pi);

If[resume == 0,

(*Create output files and initialize time counter *)

CreateFile["time_lf.txt"];
geni = 0;

CreateFile["data_lf.txt"];

(** Construct a parameter table **)

CreateFile["tab_lf.txt"];
nr=0;
For[imod=1,imod<= Length[modelL],imod++, 
For[icorr=1,icorr<= Length[lcorr],icorr++,
For[ic=1,ic<= Length[cL],ic++,
For[iS=1,iS<=Length[SL],iS++,
For[ie=1,ie<=Length[eL],ie++,
For[iu=1,iu<= Length[uL],iu++,
For[isel=1,isel<= Length[selL],isel++,
For[ih=1,ih<= Length[hL],ih++,
For[ia=1,ia<= Length[aL],ia++,

row= {modelL[[imod]],SL[[iS]],uL[[iu]],aL[[ia]],selL[[isel]],cL[[ic]],eL[[ie]],hL[[ih]],lcorr[[icorr]]};
riff=StringRiffle[row];
WriteLine["tab_lf.txt",riff];
startline= 1;
]]]]]]]]];

matpar=Import["tab_lf.txt","Table"];
nr=Length[matpar];,

(** Import existing parameter table **)

startline = resume;
geni=resume-1;
matpar=Import["tab_lf.txt","Table"];
nr=Length[matpar];

];


Close["time_cl.txt"];
Close["data_cl.txt"];

outputdata=OpenAppend["data_lf.txt"];
outputtime=OpenAppend["time_lf.txt"];

(** Let the calculations begin or resume ... **)
For[k=startline,k<= nr, k++,

(* Switching between models *)

model = matpar[[k]][[1]];

Switch[model,"c",
vsdir = s vC1c;
vsint = s^2 vC2c2;
cV=c wG;
\[Epsilon]V=\[Epsilon];
SV=S;
W=wG;
LF= 1/c^4 (1/(1-S )-4/(E^(c/(4 \[Epsilon]))-S)+6/(E^(c/(2 \[Epsilon]))-S)-4/(E^((3 c)/(4 \[Epsilon]))-S)+1/(E^(c/\[Epsilon])-S))/.c-> c w;,
"e",
vsdir = s vC1e;
vsint = s^2 vC2e2;
vsSurv = 0;
cV=c;
\[Epsilon]V=\[Epsilon] wG;
SV=S;
W=wG;
LF= 1/c^4 (1/(1-S )-4/(E^(c/(4 \[Epsilon]))-S)+6/(E^(c/(2 \[Epsilon]))-S)-4/(E^((3 c)/(4 \[Epsilon]))-S)+1/(E^(c/\[Epsilon])-S))/.\[Epsilon]-> \[Epsilon] w;,
"s",
vsdir = s vC1s;
vsint = s^2 vC2s2;
vsSurv = s;
cV=c;
\[Epsilon]V=\[Epsilon];
SV=S wS;
W=wS;
LF=w/c^4 (1/(1-S )-4/(E^(c/(4 \[Epsilon]))-S)+6/(E^(c/(2 \[Epsilon]))-S)-4/(E^((3 c)/(4 \[Epsilon]))-S)+1/(E^(c/\[Epsilon])-S))/.S-> S MeanW w;
];

(* Setting parameters *)

par = {c-> matpar[[k]][[6]], S->  matpar[[k]][[2]],  \[Epsilon]-> matpar[[k]][[7]],  u-> matpar[[k]][[3]], h->  matpar[[k]][[8]], \[Alpha]->  matpar[[k]][[4]], s->  matpar[[k]][[5]], r-> 0.5};
corr= matpar[[k]][[9]];

not1={pj-> pNI};
not1NI={pj->0};
not2={Subscript[s, dir]->vsdir,Subscript[s, int]-> vsint,F-> \[Alpha]/(2-\[Alpha]),R-> 2r(1-r),G-> (4\[Alpha](1-\[Alpha]))/((2-\[Alpha])^2 (4-\[Alpha]))};
not3={z-> 0, c-> cV,S-> SV, \[Epsilon]-> \[Epsilon]V};
not4={F-> \[Alpha]/(2-\[Alpha])};

Quiet[
tp=0;
For[j=0, j<= imax+2,j++,
For[i=1,i<= 9,i++,
ptp=tp;
ndpi=dpi2/.not1/.not2/.not3/.not4/.par/.pi-> ptp+i*10^(imax-j);
If[ndpi <0,
tp =ptp+(i-1)*10^(imax-j);
 Break[];,
l=dpi2/.not1/.not2/.not3/.not4/.par/.pi-> ptp+i*10^(imax-j)-Range[0,10^(imax-j)-10^(imax-j-1),10^(imax-j-1)] //Quiet;
If[Negative @ Min[DeleteCases[l, Indeterminate]],
 tp = ptp +(i-1)*10^(imax-j);
Break[];
];
];
If[And[ndpi > 0, i==9], 
tp =ptp+i*10^(imax-j);
Break[];
];
];
];
p=tp//N;
];

fit=W/.not1/.not2/.not3/.not4/.par/.pi-> p;
Ws=W/.F-> 1/2 (1+F)/.not1/.not2/.not3/.not4/.par/.pi-> p;
Wo=W/.F->0/.not1/.not2/.not3/.not4/.par/.pi-> p;

If[corr==1,
LFws=LF/.w->Ws/.MeanW-> fit;
LFwo=LF/.w-> Wo/.MeanW-> fit;,

LFws=LF/.w->Ws/.MeanW-> 1;
LFwo=LF/.w-> Wo/.MeanW-> 1;
];

\[Delta]=1-LFws/LFwo /.par //N;

Quiet[
tp=0;
For[j=0, j<= imax+2,j++,
For[i=1,i<= 9,i++,
ptp=tp;
ndpi=dpi2/.not1NI/.not2/.not3/.not4/.par/.pi-> ptp+i*10^(imax-j);
If[ndpi < 0,
tp =ptp+(i-1)*10^(imax-j);
 Break[];,
l=dpi2/.not1/.not2/.not3/.not4/.par/.pi-> ptp+i*10^(imax-j)-Range[0,10^(imax-j)-10^(imax-j-1),10^(imax-j-1)] //Quiet;
If[Negative @ Min[DeleteCases[l, Indeterminate]],
 tp = ptp +(i-1)*10^(imax-j);
Break[];
];
];
If[And[ndpi > 0, i==9], 
tp =ptp+i*10^(imax-j);
Break[];
];
];
];
pNoInteraction=tp;
];

row={model,S/.par,u/.par,\[Alpha]/.par,s/.par,c/.par,\[Epsilon]/.par,h/.par,p,pNoInteraction,corr, \[Delta],fit};
riff=StringRiffle[row];

WriteLine[outputdata,riff];

geni=geni+1;
rifftime=StringRiffle[{geni,timetot}," / "];
WriteLine[outputtime,rifftime];

];
