
////////////////////////////////////////////////////////////////
// 2D wellbore with circular region
// Created 28/05/2018 by Manouchehr Sanei and Omar Duran
// Labmec, State University of Campinas, Brazil
////////////////////////////////////////////////////////////////

IsHexahedralQ = 1;
 
Mesh.ElementOrder = 1;
Mesh.SecondOrderLinear = 0;

wr = 0.1;
fr = 10.0;

nt = 8;
nr = 16;
h  = 5.0;
nh = 4;

radial_progression = 1.41875;

// center point
pc = newp; Point(pc) = {0,0,-h/2};

// internal circle
pi_1 = newp; Point(pi_1) = {wr,0,-h/2};
pi_2 = newp; Point(pi_2) = {0,wr,-h/2};
pi_3 = newp; Point(pi_3) = {-wr,0,-h/2};
pi_4 = newp; Point(pi_4) = {0,-wr,-h/2};

li_1 = newl; Circle(li_1) = {pi_1,pc,pi_2};
li_2 = newl; Circle(li_2) = {pi_2,pc,pi_3};
li_3 = newl; Circle(li_3) = {pi_3,pc,pi_4};
li_4 = newl; Circle(li_4) = {pi_4,pc,pi_1};

i_circle[] = {li_1,li_2,li_3,li_4};

// external circle
pe_1 = newp; Point(pe_1) = {fr,0,-h/2};
pe_2 = newp; Point(pe_2) = {0,fr,-h/2};
pe_3 = newp; Point(pe_3) = {-fr,0,-h/2};
pe_4 = newp; Point(pe_4) = {0,-fr,-h/2};

le_1 = newl; Circle(le_1) = {pe_1,pc,pe_2};
le_2 = newl; Circle(le_2) = {pe_2,pc,pe_3};
le_3 = newl; Circle(le_3) = {pe_3,pc,pe_4};
le_4 = newl; Circle(le_4) = {pe_4,pc,pe_1};

e_circle[] = {le_1,le_2,le_3,le_4};

// Auxiliary geometrical entities
lr1 = newl; Line(lr1) = {pi_1,pe_1};
lr2 = newl; Line(lr2) = {pi_2,pe_2};
lr3 = newl; Line(lr3) = {pi_3,pe_3};
lr4 = newl; Line(lr4) = {pi_4,pe_4};

radial_lines[] = {lr1,lr2,lr3,lr4};
azimuthal_lines[] = {i_circle[],e_circle[]};

Transfinite Line {azimuthal_lines[]} = nt;
Transfinite Line {radial_lines[]} = nr Using Progression radial_progression;

llq_1 = newll; Line Loop(llq_1) = {-li_1,lr1,le_1,-lr2};
llq_2 = newll; Line Loop(llq_2) = {-li_2,lr2,le_2,-lr3};
llq_3 = newll; Line Loop(llq_3) = {-li_3,lr3,le_3,-lr4};
llq_4 = newll; Line Loop(llq_4) = {-li_4,lr4,le_4,-lr1};

s1 = news; Plane Surface(s1) = {llq_1};
s2 = news; Plane Surface(s2) = {llq_2};
s3 = news; Plane Surface(s3) = {llq_3};
s4 = news; Plane Surface(s4) = {llq_4};

the_circle[] = {s1,s2,s3,s4};

Transfinite Surface"*";
If(IsHexahedralQ)
	Recombine Surface"*";
	out[] = Extrude {0,0,h} {
  		Surface{the_circle[]}; Layers{nh}; Recombine;
	};  
Else
	out[] = Extrude {0,0,h} {
  		Surface{the_circle[]}; Layers{nh};
	}; 	
EndIf

// external circle mid points
pem_1 = newp; Point(pem_1) = {fr,0,0.0};
pem_2 = newp; Point(pem_2) = {0,fr,0.0};
pem_3 = newp; Point(pem_3) = {-fr,0,0.0};
pem_4 = newp; Point(pem_4) = {0,-fr,0.0};


the_volume[] = {1,2,3,4};
NW_bc[] = {37};
NE_bc[] = {59};
SE_bc[] = {81};
SW_bc[] = {103};
Bottom_bc[] = {17,18,19,20};
Top_bc[] = {42,64,86,108};
Wellbore_bc[] = {29,51,73,95};

fixed_y_points[]={pem_1,pem_3};
fixed_x_points[]={pem_2,pem_4};

Point{fixed_y_points[],fixed_x_points[]} In Surface{37};
Point{fixed_y_points[],fixed_x_points[]} In Surface{59};
Point{fixed_y_points[],fixed_x_points[]} In Surface{81};
Point{fixed_y_points[],fixed_x_points[]} In Surface{103};

Physical Volume("Omega") = {the_volume[]};
Physical Surface("BCWest") = {NW_bc[]}; 
Physical Surface("BCSouth") = {NE_bc[]};
Physical Surface("BCTop") = {Top_bc[]};
Physical Surface("BCNorth") = {SE_bc[]};
Physical Surface("BCBottom") = {Bottom_bc[]};
Physical Surface("BCEast") = {SW_bc[]};
Physical Point("fixed_x") = {fixed_x_points[]};
Physical Point("fixed_y") = {fixed_y_points[]};


Coherence Mesh;


