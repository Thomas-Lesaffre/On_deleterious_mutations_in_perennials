#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Input::Initialization:: *)
SetDirectory[NotebookDirectory[]];
vs\[Gamma]2s=(s^2 (E^(c/(4 \[Epsilon]))-E^z S) (E^(c/(2 \[Epsilon]))-E^z S) (E^((3 c)/(4 \[Epsilon]))-E^z S) (E^(c/\[Epsilon])-E^z S) ((E^(2 z) S^2 (3 E^(c/(4 \[Epsilon]))+5 E^(c/(2 \[Epsilon]))+3 E^((3 c)/(4 \[Epsilon]))+3 E^z S) Zi^2)/((E^(c/(4 \[Epsilon]))-E^z S) (E^(c/(2 \[Epsilon]))-E^z S) (E^((3 c)/(4 \[Epsilon]))-E^z S) (E^(c/\[Epsilon])-E^z S))-E^z S (3 E^((3 c)/(4 \[Epsilon]))+5 E^(c/\[Epsilon])+3 E^((5 c)/(4 \[Epsilon]))+6 E^(z+c/(4 \[Epsilon])) S+10 E^(z+c/(2 \[Epsilon])) S+6 E^(z+(3 c)/(4 \[Epsilon])) S+3 E^(2 z) S^2) Zi (-((E^z S Zi)/((E^(c/(4 \[Epsilon]))-E^z S) (E^(c/(2 \[Epsilon]))-E^z S) (E^((3 c)/(4 \[Epsilon]))-E^z S) (E^(c/\[Epsilon])-E^z S)^2))+(-((E^z S Zi)/((E^(c/(4 \[Epsilon]))-E^z S) (E^(c/(2 \[Epsilon]))-E^z S) (E^((3 c)/(4 \[Epsilon]))-E^z S)^2))+(-((E^z S Zi)/((E^(c/(4 \[Epsilon]))-E^z S) (E^(c/(2 \[Epsilon]))-E^z S)^2))-(E^z S Zi)/((E^(c/(4 \[Epsilon]))-E^z S)^2 (E^(c/(2 \[Epsilon]))-E^z S)))/(E^((3 c)/(4 \[Epsilon]))-E^z S))/(E^(c/\[Epsilon])-E^z S))+(E^((3 c)/(2 \[Epsilon]))+3 E^(z+(3 c)/(4 \[Epsilon])) S+5 E^(z+c/\[Epsilon]) S+3 E^(z+(5 c)/(4 \[Epsilon])) S+3 E^(2 z+c/(4 \[Epsilon])) S^2+5 E^(2 z+c/(2 \[Epsilon])) S^2+3 E^(2 z+(3 c)/(4 \[Epsilon])) S^2+E^(3 z) S^3) ((E^(2 z) S^2 Zi^2)/((E^(c/(4 \[Epsilon]))-E^z S) (E^(c/(2 \[Epsilon]))-E^z S) (E^((3 c)/(4 \[Epsilon]))-E^z S) (E^(c/\[Epsilon])-E^z S)^3)-(E^z S Zi (-((E^z S Zi)/((E^(c/(4 \[Epsilon]))-E^z S) (E^(c/(2 \[Epsilon]))-E^z S) (E^((3 c)/(4 \[Epsilon]))-E^z S)^2))+(-((E^z S Zi)/((E^(c/(4 \[Epsilon]))-E^z S) (E^(c/(2 \[Epsilon]))-E^z S)^2))-(E^z S Zi)/((E^(c/(4 \[Epsilon]))-E^z S)^2 (E^(c/(2 \[Epsilon]))-E^z S)))/(E^((3 c)/(4 \[Epsilon]))-E^z S)))/(E^(c/\[Epsilon])-E^z S)^2+1/(E^(c/\[Epsilon])-E^z S) ((E^(2 z) S^2 Zi^2)/((E^(c/(4 \[Epsilon]))-E^z S) (E^(c/(2 \[Epsilon]))-E^z S) (E^((3 c)/(4 \[Epsilon]))-E^z S)^3)-(E^z S Zi (-((E^z S Zi)/((E^(c/(4 \[Epsilon]))-E^z S) (E^(c/(2 \[Epsilon]))-E^z S)^2))-(E^z S Zi)/((E^(c/(4 \[Epsilon]))-E^z S)^2 (E^(c/(2 \[Epsilon]))-E^z S))))/(E^((3 c)/(4 \[Epsilon]))-E^z S)^2+((E^(2 z) S^2 Zi^2)/((E^(c/(4 \[Epsilon]))-E^z S) (E^(c/(2 \[Epsilon]))-E^z S)^3)+(E^(2 z) S^2 Zi^2)/((E^(c/(4 \[Epsilon]))-E^z S)^2 (E^(c/(2 \[Epsilon]))-E^z S)^2)+(E^(2 z) S^2 Zi^2)/((E^(c/(4 \[Epsilon]))-E^z S)^3 (E^(c/(2 \[Epsilon]))-E^z S)))/(E^((3 c)/(4 \[Epsilon]))-E^z S)))))/(E^((3 c)/(2 \[Epsilon]))+3 E^(z+(3 c)/(4 \[Epsilon])) S+5 E^(z+c/\[Epsilon]) S+3 E^(z+(5 c)/(4 \[Epsilon])) S+3 E^(2 z+c/(4 \[Epsilon])) S^2+5 E^(2 z+c/(2 \[Epsilon])) S^2+3 E^(2 z+(3 c)/(4 \[Epsilon])) S^2+E^(3 z) S^3)/.Zi-> 1;
vs\[Gamma]1s=s ((E^z) S ((4 (-1+E^(c/(4 \[Epsilon]))))/(E^(c/(4 \[Epsilon]))-E^z S)^2-(6 (-1+E^(c/(2 \[Epsilon]))))/(E^(c/(2 \[Epsilon]))-E^z S)^2+(4 (-1+E^((3 c)/(4 \[Epsilon]))))/(E^((3 c)/(4 \[Epsilon]))-E^z S)^2-(-1+E^(c/\[Epsilon]))/(E^(c/\[Epsilon])-E^z S)^2) )/(1+(6-6 E^z S)/(E^(c/(2 \[Epsilon]))-E^z S)+(1-E^z S)/(E^(c/\[Epsilon])-E^z S)+(4 (-1+E^z S))/(E^(c/(4 \[Epsilon]))-E^z S)+(4 (-1+E^z S))/(E^((3 c)/(4 \[Epsilon]))-E^z S));
vs\[Gamma]1e=(c dz E^-z (-E^(((4 c E^-z)/\[Epsilon]))-4 E^((13 c E^-z)/(4 \[Epsilon])) S-4 E^((7 c E^-z)/(2 \[Epsilon])) S+E^((4 c E^-z)/\[Epsilon]) S+3 E^((9 c E^-z)/(4 \[Epsilon])) S^2+5 E^((5 c E^-z)/(2 \[Epsilon])) S^2+3 E^((11 c E^-z)/(4 \[Epsilon])) S^2+5 E^((3 c E^-z)/\[Epsilon]) S^2+7 E^((13 c E^-z)/(4 \[Epsilon])) S^2+4 E^((7 c E^-z)/(2 \[Epsilon])) S^2-3 E^((9 c E^-z)/(4 \[Epsilon])) S^3-5 E^((5 c E^-z)/(2 \[Epsilon])) S^3-3 E^((11 c E^-z)/(4 \[Epsilon])) S^3-5 E^((3 c E^-z)/\[Epsilon]) S^3-3 E^((13 c E^-z)/(4 \[Epsilon])) S^3-3 E^((c E^-z)/\[Epsilon]) S^4-5 E^((5 c E^-z)/(4 \[Epsilon])) S^4-3 E^((3 c E^-z)/(2 \[Epsilon])) S^4-5 E^((7 c E^-z)/(4 \[Epsilon])) S^4-3 E^((2 c E^-z)/\[Epsilon]) S^4+4 E^((3 c E^-z)/(4 \[Epsilon])) S^5+7 E^((c E^-z)/\[Epsilon]) S^5+5 E^((5 c E^-z)/(4 \[Epsilon])) S^5+3 E^((3 c E^-z)/(2 \[Epsilon])) S^5+5 E^((7 c E^-z)/(4 \[Epsilon])) S^5+3 E^((2 c E^-z)/\[Epsilon]) S^5+E^((c E^-z)/(4 \[Epsilon])) S^6-4 E^((3 c E^-z)/(4 \[Epsilon])) S^6-4 E^((c E^-z)/\[Epsilon]) S^6-E^((c E^-z)/(4 \[Epsilon])) S^7))/((-1+E^((c E^-z)/(4 \[Epsilon]))) (E^((c E^-z)/(4 \[Epsilon]))-S) (E^((c E^-z)/(2 \[Epsilon]))-S) (E^((3 c E^-z)/(4 \[Epsilon]))-S) (E^((c E^-z)/\[Epsilon])-S) (E^((3 c E^-z)/(2 \[Epsilon]))+3 E^((3 c E^-z)/(4 \[Epsilon])) S+5 E^((c E^-z)/\[Epsilon]) S+3 E^((5 c E^-z)/(4 \[Epsilon])) S+3 E^((c E^-z)/(4 \[Epsilon])) S^2+5 E^((c E^-z)/(2 \[Epsilon])) S^2+3 E^((3 c E^-z)/(4 \[Epsilon])) S^2+S^3) \[Epsilon])/.dz-> -s;

vs\[Gamma]2e=(c dz^2 E^(-2 z+(c E^-z)/(4 \[Epsilon])) (-1+S) (c (E^((13 c E^-z)/(2 \[Epsilon]))+E^((25 c E^-z)/(4 \[Epsilon])) (-4+S)+17 E^((6 c E^-z)/\[Epsilon]) S-17 E^((c E^-z)/(2 \[Epsilon])) S^9-S^10+3 E^((3 c E^-z)/(4 \[Epsilon])) S^8 (7+S)+E^((c E^-z)/(4 \[Epsilon])) S^9 (-1+4 S)-3 E^((23 c E^-z)/(4 \[Epsilon])) S (1+7 S)+64 E^((13 c E^-z)/(4 \[Epsilon])) S^4 (-1+S^2)+4 E^((5 c E^-z)/(4 \[Epsilon])) S^7 (3-12 S+S^2)+3 E^((5 c E^-z)/\[Epsilon]) S^2 (11-15 S+2 S^2)-4 E^((21 c E^-z)/(4 \[Epsilon])) S (1-12 S+3 S^2)+E^((11 c E^-z)/(2 \[Epsilon])) S (-23+3 S+6 S^2)+5 E^((5 c E^-z)/(2 \[Epsilon])) S^5 (15-13 S+9 S^2)-3 E^((3 c E^-z)/(2 \[Epsilon])) S^6 (2-15 S+11 S^2)-3 E^((17 c E^-z)/(4 \[Epsilon])) S^3 (26-44 S+15 S^2)-5 E^((4 c E^-z)/\[Epsilon]) S^3 (9-13 S+15 S^2)+E^((19 c E^-z)/(4 \[Epsilon])) S^2 (33-93 S+20 S^2)+E^((c E^-z)/\[Epsilon]) S^7 (-6-3 S+23 S^2)+3 E^((9 c E^-z)/(4 \[Epsilon])) S^5 (15-44 S+26 S^2)-E^((7 c E^-z)/(4 \[Epsilon])) S^6 (20-93 S+33 S^2)-E^((2 c E^-z)/\[Epsilon]) S^6 (121-87 S+39 S^2)+E^((9 c E^-z)/(2 \[Epsilon])) S^2 (39-87 S+121 S^2)-E^((11 c E^-z)/(4 \[Epsilon])) S^4 (9-45 S+145 S^2+3 S^3)+E^((7 c E^-z)/(2 \[Epsilon])) S^3 (-9-33 S-195 S^2+5 S^3)+E^((15 c E^-z)/(4 \[Epsilon])) S^3 (3+145 S-45 S^2+9 S^3)+E^((3 c E^-z)/\[Epsilon]) S^4 (-5+195 S+33 S^2+9 S^3))-8 E^z (-1+E^((c E^-z)/(4 \[Epsilon]))) (E^((25 c E^-z)/(4 \[Epsilon]))+3 E^((23 c E^-z)/(4 \[Epsilon])) S-E^((6 c E^-z)/\[Epsilon]) S-3 E^((3 c E^-z)/(4 \[Epsilon])) (-2+S) S^8+E^((c E^-z)/(4 \[Epsilon])) S^9-3 E^((c E^-z)/(2 \[Epsilon])) S^9-S^10-3 E^((11 c E^-z)/(2 \[Epsilon])) S (-1+2 S)-6 E^((9 c E^-z)/(4 \[Epsilon])) S^6 (-4+3 S)+3 E^((5 c E^-z)/(4 \[Epsilon])) S^7 (-4+3 S)+3 E^((5 c E^-z)/\[Epsilon]) S^2 (-3+4 S)-6 E^((4 c E^-z)/\[Epsilon]) S^3 (-3+4 S)-3 E^((17 c E^-z)/(4 \[Epsilon])) S^3 (-9+5 S)-2 E^((5 c E^-z)/(2 \[Epsilon])) S^6 (-13+6 S)-3 E^((2 c E^-z)/\[Epsilon]) S^6 (-5+9 S)-2 E^((15 c E^-z)/(4 \[Epsilon])) S^3 (-6+13 S)-3 E^((19 c E^-z)/(4 \[Epsilon])) S^2 (4-6 S+S^2)+E^((c E^-z)/\[Epsilon]) S^7 (-3+12 S+S^2)+E^((21 c E^-z)/(4 \[Epsilon])) S (-1-12 S+3 S^2)+3 E^((3 c E^-z)/(2 \[Epsilon])) S^6 (1-6 S+4 S^2)+E^((7 c E^-z)/(4 \[Epsilon])) S^6 (8-27 S+6 S^2)+E^((13 c E^-z)/(4 \[Epsilon])) S^4 (-15+8 S^2)-E^((9 c E^-z)/(2 \[Epsilon])) S^2 (6-27 S+8 S^2)+E^((3 c E^-z)/\[Epsilon]) S^4 (-8+15 S^2)+3 E^((7 c E^-z)/(2 \[Epsilon])) S^3 (1-8 S+S^3)-3 E^((11 c E^-z)/(4 \[Epsilon])) S^4 (1-8 S^2+S^3)) \[Epsilon]))/(8 (-1+E^((c E^-z)/(4 \[Epsilon])))^2 (E^((c E^-z)/(4 \[Epsilon]))-S)^2 (E^((c E^-z)/(2 \[Epsilon]))-S)^2 (E^((3 c E^-z)/(4 \[Epsilon]))-S)^2 (E^((c E^-z)/\[Epsilon])-S)^2 (E^((3 c E^-z)/(2 \[Epsilon]))+5 E^((c E^-z)/\[Epsilon]) S+3 E^((5 c E^-z)/(4 \[Epsilon])) S+3 E^((c E^-z)/(4 \[Epsilon])) S^2+5 E^((c E^-z)/(2 \[Epsilon])) S^2+S^3+3 E^((3 c E^-z)/(4 \[Epsilon])) S (1+S)) \[Epsilon]^2)/.dz-> s;

vs\[Gamma]1c=-s (c E^(z+(c E^z)/(4 \[Epsilon])) (-E^(((15 c E^z)/(4 \[Epsilon]))) (-1+S)-4 E^((13 c E^z)/(4 \[Epsilon])) (-1+S) S+3 E^((2 c E^z)/\[Epsilon]) (-1+S) S^2+5 E^((9 c E^z)/(4 \[Epsilon])) (-1+S) S^2+3 E^((5 c E^z)/(2 \[Epsilon])) (-1+S) S^2+5 E^((11 c E^z)/(4 \[Epsilon])) (-1+S) S^2-5 E^((c E^z)/\[Epsilon]) (-1+S) S^4-3 E^((5 c E^z)/(4 \[Epsilon])) (-1+S) S^4-5 E^((3 c E^z)/(2 \[Epsilon])) (-1+S) S^4-3 E^((7 c E^z)/(4 \[Epsilon])) (-1+S) S^4+4 E^((c E^z)/(2 \[Epsilon])) (-1+S) S^5+(-1+S) S^6+E^((3 c E^z)/\[Epsilon]) S (4-7 S+3 S^2)+E^((3 c E^z)/(4 \[Epsilon])) S^4 (3-7 S+4 S^2))-4 (-1+E^((c E^z)/(4 \[Epsilon]))) (E^((4 c E^z)/\[Epsilon])+2 E^((15 c E^z)/(4 \[Epsilon])) S-6 E^((11 c E^z)/(4 \[Epsilon])) S^2-6 E^((5 c E^z)/(4 \[Epsilon])) S^5+2 E^((3 c E^z)/(4 \[Epsilon])) (-2+S) S^5+2 E^((c E^z)/(4 \[Epsilon])) S^6+S^7+5 E^((2 c E^z)/\[Epsilon]) S^3 (1+S)-E^((c E^z)/\[Epsilon]) S^5 (5+S)-2 E^((13 c E^z)/(4 \[Epsilon])) S (-1+2 S)-E^((7 c E^z)/(2 \[Epsilon])) S (-4+3 S)+E^((c E^z)/(2 \[Epsilon])) S^5 (-3+4 S)-E^((3 c E^z)/\[Epsilon]) S (1+5 S)-2 E^((7 c E^z)/(4 \[Epsilon])) S^3 (-2-2 S+S^2)+2 E^((9 c E^z)/(4 \[Epsilon])) S^2 (-1+2 S+2 S^2)+E^((5 c E^z)/(2 \[Epsilon])) S^2 (-7+3 S+3 S^2)-E^((3 c E^z)/(2 \[Epsilon])) S^3 (-3-3 S+7 S^2)) \[Epsilon])/((-1+E^((c E^z)/(4 \[Epsilon]))) (E^((c E^z)/(4 \[Epsilon]))-S) (E^((c E^z)/(2 \[Epsilon]))-S) (E^((3 c E^z)/(4 \[Epsilon]))-S) (E^((c E^z)/\[Epsilon])-S) (E^((3 c E^z)/(2 \[Epsilon]))+5 E^((c E^z)/\[Epsilon]) S+3 E^((5 c E^z)/(4 \[Epsilon])) S+3 E^((c E^z)/(4 \[Epsilon])) S^2+5 E^((c E^z)/(2 \[Epsilon])) S^2+S^3+3 E^((3 c E^z)/(4 \[Epsilon])) S (1+S)) \[Epsilon]);
vs\[Gamma]2c=s^2/(32 \[Epsilon]^2) ((4 c^2 E^(2 z+(c E^z)/(2 \[Epsilon])) (E^((c E^z)/(2 \[Epsilon]))+S))/(E^((c E^z)/(2 \[Epsilon]))-S)^2+(9 c^2 E^(2 z+(3 c E^z)/(4 \[Epsilon])) (E^((3 c E^z)/(4 \[Epsilon]))+S))/(E^((3 c E^z)/(4 \[Epsilon]))-S)^2+(16 c^2 E^(2 z+(c E^z)/\[Epsilon]) (E^((c E^z)/\[Epsilon])+S))/(E^((c E^z)/\[Epsilon])-S)^2+(c^2 E^(2 z+(c E^z)/(4 \[Epsilon])) (36 E^((5 c E^z)/(4 \[Epsilon]))+80 E^((3 c E^z)/(4 \[Epsilon])) S+75 E^((c E^z)/\[Epsilon]) S+3 S^2+20 E^((c E^z)/(4 \[Epsilon])) S^2+27 E^((c E^z)/(2 \[Epsilon])) S (1+S)))/(E^((3 c E^z)/(2 \[Epsilon]))+5 E^((c E^z)/\[Epsilon]) S+3 E^((5 c E^z)/(4 \[Epsilon])) S+3 E^((c E^z)/(4 \[Epsilon])) S^2+5 E^((c E^z)/(2 \[Epsilon])) S^2+S^3+3 E^((3 c E^z)/(4 \[Epsilon])) S (1+S))-(6 c^2 E^(2 z+(c E^z)/\[Epsilon]) (E^((3 c E^z)/(4 \[Epsilon]))+E^((c E^z)/(2 \[Epsilon])) (3-2 S)-5 E^((c E^z)/(4 \[Epsilon])) S+S (-1+4 S)))/((-1+E^((c E^z)/(4 \[Epsilon]))) (E^((c E^z)/(4 \[Epsilon]))-S) (E^((c E^z)/(2 \[Epsilon]))-S) (E^((3 c E^z)/(4 \[Epsilon]))-S))+(c^2 E^(2 z+(c E^z)/(4 \[Epsilon])) (9 E^((3 c E^z)/(4 \[Epsilon]))+E^((c E^z)/(2 \[Epsilon])) (2-23 S)+S-4 S^2+E^((c E^z)/(4 \[Epsilon])) (1-2 S+16 S^2)))/((-1+E^((c E^z)/(4 \[Epsilon])))^2 (E^((c E^z)/(4 \[Epsilon]))-S)^2)-(32 c E^(z+(c E^z)/(4 \[Epsilon])) (1+3 E^((c E^z)/(4 \[Epsilon]))-4 S) \[Epsilon])/((-1+E^((c E^z)/(4 \[Epsilon]))) (E^((c E^z)/(4 \[Epsilon]))-S))+(96 c E^(z+(3 c E^z)/(4 \[Epsilon])) \[Epsilon])/(E^((3 c E^z)/(4 \[Epsilon]))-S)+320 \[Epsilon]^2+(4 c E^(z+(c E^z)/(2 \[Epsilon])) (-((c E^(z+(c E^z)/(4 \[Epsilon])) (1+3 E^((c E^z)/(4 \[Epsilon]))-4 S))/((-1+E^((c E^z)/(4 \[Epsilon]))) (E^((c E^z)/(4 \[Epsilon]))-S)))+16 \[Epsilon]))/(E^((c E^z)/(2 \[Epsilon]))-S)+(8 c E^(z+(c E^z)/\[Epsilon]) (c E^(z+(c E^z)/(4 \[Epsilon])) (2 E^((3 c E^z)/(2 \[Epsilon]))+5 E^((c E^z)/\[Epsilon]) S-E^((3 c E^z)/(4 \[Epsilon])) (-5+S) S-5 E^((c E^z)/(4 \[Epsilon])) S^2-E^((5 c E^z)/(4 \[Epsilon])) (6+S)+S^2 (-1+4 S)-E^((c E^z)/(2 \[Epsilon])) S (-3+5 S))+16 (E^((7 c E^z)/(4 \[Epsilon]))+E^((3 c E^z)/(4 \[Epsilon])) S+E^((c E^z)/\[Epsilon]) S^2+S^3-E^((3 c E^z)/(2 \[Epsilon])) (1+S)-E^((c E^z)/(4 \[Epsilon])) S^2 (1+S)) \[Epsilon]))/((-1+E^((c E^z)/(4 \[Epsilon]))) (E^((c E^z)/(4 \[Epsilon]))-S) (E^((c E^z)/(2 \[Epsilon]))-S) (E^((3 c E^z)/(4 \[Epsilon]))-S) (E^((c E^z)/\[Epsilon])-S))-(2 c E^(z+(c E^z)/(4 \[Epsilon])) (6 E^((5 c E^z)/(4 \[Epsilon]))+20 E^((3 c E^z)/(4 \[Epsilon])) S+15 E^((c E^z)/\[Epsilon]) S+3 S^2+10 E^((c E^z)/(4 \[Epsilon])) S^2+9 E^((c E^z)/(2 \[Epsilon])) S (1+S)) (4 c E^(z+(c E^z)/\[Epsilon]) (-1+E^((c E^z)/(4 \[Epsilon]))) (E^((c E^z)/(4 \[Epsilon]))-S) (E^((c E^z)/(2 \[Epsilon]))-S) (E^((3 c E^z)/(4 \[Epsilon]))-S)+(E^((c E^z)/\[Epsilon])-S) (c E^(z+(c E^z)/(4 \[Epsilon])) (2 E^((3 c E^z)/(2 \[Epsilon]))+5 E^((c E^z)/\[Epsilon]) S-E^((3 c E^z)/(4 \[Epsilon])) (-5+S) S-5 E^((c E^z)/(4 \[Epsilon])) S^2-E^((5 c E^z)/(4 \[Epsilon])) (6+S)+S^2 (-1+4 S)-E^((c E^z)/(2 \[Epsilon])) S (-3+5 S))+16 (E^((7 c E^z)/(4 \[Epsilon]))+E^((3 c E^z)/(4 \[Epsilon])) S+E^((c E^z)/\[Epsilon]) S^2+S^3-E^((3 c E^z)/(2 \[Epsilon])) (1+S)-E^((c E^z)/(4 \[Epsilon])) S^2 (1+S)) \[Epsilon])))/((-1+E^((c E^z)/(4 \[Epsilon]))) (E^((c E^z)/(4 \[Epsilon]))-S) (E^((c E^z)/(2 \[Epsilon]))-S) (E^((3 c E^z)/(4 \[Epsilon]))-S) (E^((c E^z)/\[Epsilon])-S) (E^((3 c E^z)/(2 \[Epsilon]))+5 E^((c E^z)/\[Epsilon]) S+3 E^((5 c E^z)/(4 \[Epsilon])) S+3 E^((c E^z)/(4 \[Epsilon])) S^2+5 E^((c E^z)/(2 \[Epsilon])) S^2+S^3+3 E^((3 c E^z)/(4 \[Epsilon])) S (1+S))));

(* Survival *)
pis= pi-sSurv h pi-sSurv(1-h)dii-sSurv(1-2h)dijj+sSurv^2 (1-h)(1-2h)G pi pj;
pjs= pj-sSurv h pj-s(1-h)djj-sSurv(1-2h)dijj+sSurv^2 (1-h)(1-2h)G pi pj;

diis=dii-sSurv(1-2h)G pi pj;
djjs=djj-sSurv(1-2h)G pi pj;

dijjs=dijj-sSurv(1-h)G pi pj;
dijijs=dijij;

(* Gametes production *)
pig= pi-s\[Gamma]1 h pi-s\[Gamma]1(1-h)dii-s\[Gamma]1(1-2h)dijj+s\[Gamma]2(1-h)(1-2h)G pi pj;
pjg= pj-s\[Gamma]1 h pj-s\[Gamma]1(1-h)djj-s\[Gamma]1(1-2h)dijj+s\[Gamma]2(1-h)(1-2h)G pi pj;

diig=dii-s\[Gamma]1(1-2h)G pi pj;
djjg=djj-s\[Gamma]1(1-2h)G pi pj;

dijjg=dijj -s\[Gamma]1(1-h)G pi pj;
dijijg=dijij;


(* Syngamy *)

pir=pig + u;
pjr=pjg+u;

diir=F/(1+F) (diig+pig);
djjr=F/(1+F) (djjg+pjg);
dijjr=F/(1+F) (dijjg);
dijijr=F/(1+F) ((1-R)(dijijg + pig pjg) + R(pig djjg + pjg diig));

(* Juvenile survival *)

pij= pir-sSurv h pir-sSurv(1-h)diir-sSurv(1-2h)dijjr+sSurv^2 (1-h)(1-2h)G pi pj;
pjj= pjr-sSurv h pjr-sSurv(1-h)djjr-sSurv(1-2h)dijjr+sSurv^2 (1-h)(1-2h)G pi pj;

diij=diir-sSurv(1-2h)G pi pj;
djjj=djjr-sSurv(1-2h)G pi pj;

dijjj=dijjr -sSurv(1-h)G pi pj;
dijijj=dijijr;

(* Next Generation *)

pin=(1-S)pij+S pis;
pjn=(1-S)pjj+S pjs;

diin=(1-S)diij+S diis;
djjn=(1-S)djjj+S djjs;

dijjn=(1-S)dijjj+S dijjs;
dijijn=(1-S)dijijj+S dijijs;


\[CapitalDelta]dijij=Normal[Series[dijijn-dijij/.dijj-> dijj sm/.dijij-> dijij /.dij-> dij sm^2/.divj-> divj sm^2/.sSurv-> sSurv sm/.s\[Gamma]1->s\[Gamma]1 sm/.s\[Gamma]2-> s\[Gamma]2 sm^2/.dii-> F pi/.djj-> F pj,{sm,0,1}]]/.sm->1;
solDijij=Solve[\[CapitalDelta]dijij==0,dijij]//Simplify;
tempDijij=dijij/.solDijij[[1]];
temp2Dijij=Normal[Series[tempDijij/.pi->pi sm/.pj-> pj sm,{sm,0,2}]]/.sm->1//FullSimplify;
qdijij=Normal[Series[temp2Dijij,{Subscript[s, dir],0,1}]];
neutdijij=qdijij/.s->0 /.s\[Gamma]1-> 0/.s\[Gamma]2->0 //Simplify;

\[CapitalDelta]dii=Normal[Series[diin-dii/.sSurv-> sSurv sm/.s\[Gamma]1->s\[Gamma]1 sm/.s\[Gamma]2-> s\[Gamma]2 sm^2/.dij-> dij sm^2/.divj-> divj sm^2/.dijj-> dijj sm/.dijij->neutdijij/.djj-> F pj,{sm,0,1}]]/.sm->1;
\[CapitalDelta]djj=Normal[Series[djjn-djj/.sSurv-> sSurv sm/.s\[Gamma]1->s\[Gamma]1 sm/.s\[Gamma]2-> s\[Gamma]2 sm^2/.dij-> dij sm^2/.divj-> divj sm^2/.dijj-> dijj sm/.dijij->neutdijij/.dii-> F pi,{sm,0,1}]]/.sm->1;

sDii=Solve[\[CapitalDelta]dii==0,dii]//Simplify;
sDjj=Solve[\[CapitalDelta]djj==0,djj]//Simplify;

tempDii=dii/.sDii[[1]];
tempDjj=djj/.sDjj[[1]];

t1=Normal[Series[tempDii/.sSurv-> sSurv sm/.s\[Gamma]1->s\[Gamma]1 sm/.s\[Gamma]2-> s\[Gamma]2 sm^2,{sm,0,1}]]/.sm->1;
t2=Normal[Series[t1,{pj,0,1}]];
qdii=t2;

t3=Normal[Series[tempDjj/.sSurv-> sSurv sm/.s\[Gamma]1->s\[Gamma]1 sm/.s\[Gamma]2-> s\[Gamma]2 sm^2,{sm,0,1}]]/.sm->1;
t4=Normal[Series[t3,{pj,0,1}]];
qdjj=t4;

q1dii=qdii/.pj-> 0;
q1djj=qdjj/.pi-> 0;

\[CapitalDelta]dijj=Normal[Series[dijjn-dijj/.sSurv-> sSurv sm/.s\[Gamma]1->s\[Gamma]1 sm/.s\[Gamma]2-> s\[Gamma]2 sm^2/.dijij->neutdijij/.dii-> F pi/.djj-> F pj,{sm,0,1}]]/.sm->1;
sDijj=Solve[\[CapitalDelta]dijj==0,dijj]//Simplify;
qdijj=Normal[Series[dijj/.sDijj[[1]]/.sSurv-> sSurv sm/.s\[Gamma]1->s\[Gamma]1 sm/.s\[Gamma]2-> s\[Gamma]2 sm^2,{sm,0,1}]]/.sm->1;

qle={dii->qdii,djj->qdjj,dijij->qdijij,dijj->qdijj};

dpi=Normal[Series[pin-pi/.qle/.sSurv-> sSurv sm/.s\[Gamma]1-> s\[Gamma]1 sm/.s\[Gamma]2->s\[Gamma]2 sm^2,{sm,0,2}]]/.sm->1;
dpi2=Normal[Series[dpi/.pi-> pi sm/.pj-> pj sm,{sm,0,2}]]/.sm->1;

qself = {dii-> 1/2 (pi+dii),djj->1/2 (pj+djj),dijj-> 1/2 dijj,dijij-> 1/4 (pi pj + pj dii+pi djj+dijij)};
qout={dii-> 0,djj->0,dijj->0,dijij->0};

notID = {z-> 0,pj->pi};



qle={dii->qdii,djj->qdjj,dijij->qdijij,dijj->qdijj};
sqle={dii->F pi,djj->F pj,dijij->neutdijij,dijj->0,dij-> 0 ,divj-> 0};

\[CapitalDelta]piNoInt=pin-pi/.qle/.pj->0;
solNoInt=Solve[\[CapitalDelta]piNoInt==0,pi];
pNI=pi/.solNoInt[[1]];
solp=Solve[dpi2==0,pi];
peq=pi/.solp[[1]];
peq2=Normal[Series[peq,{pj,0,1}]];
peq3=Normal[Series[peq2/.z->0/.s-> s sm/.s\[Gamma]1-> s\[Gamma]1 sm/.s\[Gamma]2->s\[Gamma]2 sm^2/.u-> u sm^2,{sm,0,2}]]/.sm->1;

DF=List[{"mod","S","U","a","s","c","e","h","mc","pmc","del","mw"}];

(*Importing parameters table and removing names *)

pf=Import["math_par_cl.txt","Table"];
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

CreateFile["time_cl.txt"];
geni = 0;

CreateFile["data_cl.txt"];

(** Construct a parameter table **)

CreateFile["tab_cl.txt"];
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
WriteLine["tab_cl.txt",riff];
startline= 1;
]]]]]]]]];

matpar=Import["tab_cl.txt","Table"];
nr=Length[matpar];,

(** Import existing parameter table **)

startline = resume;
matpar=Import["tab_cl.txt","Table"];
nr=Length[matpar];
geni=resume - 1;

];

Close["time_cl.txt"];
Close["data_cl.txt"];

outputdata=OpenAppend["data_cl.txt"];
outputtime=OpenAppend["time_cl.txt"];

(** Let the calculations begin or resume ... **)
For[k=startline,k<= nr, k++,

(* Switching between models *)

model = matpar[[k]][[1]];

Switch[model,"c",
vs\[Gamma]1 = vs\[Gamma]1c;
vs\[Gamma]2 = vs\[Gamma]2c;
vsSurv = 0;
cV=c wG;
\[Epsilon]V=\[Epsilon];
SV=S;
W=wG;
LF= 1/c^4 (1/(1-S )-4/(E^(c/(4 \[Epsilon]))-S)+6/(E^(c/(2 \[Epsilon]))-S)-4/(E^((3 c)/(4 \[Epsilon]))-S)+1/(E^(c/\[Epsilon])-S))/.c-> c w;,
"e",
vs\[Gamma]1 = vs\[Gamma]1e;
vs\[Gamma]2 = vs\[Gamma]2e;
vsSurv = 0;
cV=c;
\[Epsilon]V=\[Epsilon] wG;
SV=S;
W=wG;
LF= 1/c^4 (1/(1-S )-4/(E^(c/(4 \[Epsilon]))-S)+6/(E^(c/(2 \[Epsilon]))-S)-4/(E^((3 c)/(4 \[Epsilon]))-S)+1/(E^(c/\[Epsilon])-S))/.\[Epsilon]-> \[Epsilon] w;,
"s",
vs\[Gamma]1 = vs\[Gamma]1s;
vs\[Gamma]2 = vs\[Gamma]2s;
vsSurv = s;
cV=c;
\[Epsilon]V=\[Epsilon];
SV=S wS;
W=wS;
LF=w/c^4 (1/(1-S )-4/(E^(c/(4 \[Epsilon]))-S)+6/(E^(c/(2 \[Epsilon]))-S)-4/(E^((3 c)/(4 \[Epsilon]))-S)+1/(E^(c/\[Epsilon])-S))/.S-> S w/.S-> S MeanW;
];

(* Setting parameters *)

par = {c-> matpar[[k]][[6]], S->  matpar[[k]][[2]],  \[Epsilon]-> matpar[[k]][[7]],  u-> matpar[[k]][[3]], h->  matpar[[k]][[8]], \[Alpha]->  matpar[[k]][[4]], s->  matpar[[k]][[5]], r-> 0.5};
corr= matpar[[k]][[9]];

not1={pj-> pNI};
not1NI={pj->0};
not2={s\[Gamma]1-> vs\[Gamma]1,s\[Gamma]2-> vs\[Gamma]2,sSurv-> vsSurv,F-> \[Alpha]/(2-\[Alpha]),R-> 2r(1-r),G-> (4\[Alpha](1-\[Alpha]))/((2-\[Alpha])^2 (4-\[Alpha]))};
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



