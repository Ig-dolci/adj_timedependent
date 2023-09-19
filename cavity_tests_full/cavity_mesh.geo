lc = 0.5;

Point(1) = {0, 2, 0};
Point(2) = {0, 0, 0};
Point(3) = {2, 0,  0};
Point(4) = {2, 2,  0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1, 2, 3, 4};



Plane Surface(1) = {1};


Physical Line(0) = {1, 2, 3}; //wall
Physical Line(1) = {4}; //inflow

Physical Surface(10) = {1}; //domain


Transfinite Surface{1} = {1,2,3,4};


Transfinite Line{1} = 15;
Transfinite Line{2} = 15;
Transfinite Line{3} = 15;
Transfinite Line{4} = 15;


Recombine Surface{1};



Mesh.MshFileVersion = 2.2;


