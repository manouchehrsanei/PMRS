
////////////////////////////////////////////////////////////////
// HydrostTest_Collapse
// Created 31/12/2018 by Manouchehr Sanei
// Labmec, State University of Campinas, Brazil
////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////
// Parameters
////////////////////////////////////////////////////////////////

lc =1.0e1;
r =0.0125;
h =0.05;
nh = 11;
nr = 2;


// 2D mesh


p1 = newp; Point(p1) = {-r,-h/2,0,lc};
p2 = newp; Point(p2) = {r,-h/2,0,lc};
p3 = newp; Point(p3) = {r,h/2,0,lc};
p4 = newp; Point(p4) = {-r,h/2,0,lc};


l1 = newl; Line(l1) = {p1,p2};
l2 = newl; Line(l2) = {p2,p3};
l3 = newl; Line(l3) = {p3,p4};
l4 = newl; Line(l4) = {p4,p1};

ll1 = newll; Line Loop(ll1) = {l1,l2,l3,l4};
s1 = news; Plane Surface(s1) = {ll1};


Transfinite Line {l2,l4} = nh;
Transfinite Line {l1,l3} = nr;
Transfinite Surface {s1};
Recombine Surface"*";

plug[] = {s1};
bottom[] = {l1};
rigth[] = {l2};
top[] = {l3};
left[] = {l4};



Physical Surface("plug") = {plug[]};
Physical Line("bottom") = {bottom[]};
Physical Line("top") = {top[]};
Physical Line("left") = {left[]};
Physical Line("rigth") = {rigth[]};

