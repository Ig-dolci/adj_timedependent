//+
Point(1) = {0.0, 0.0, 0, 1.0};
Point(2) = {1.0, 0.0, 0, 1.0};
Point(3) = {0.0, 1.0, 0, 1.0};
Point(4) = {1.0, 1.0, 0, 1.0};//+
Line(1) = {3, 4};
//+
Line(2) = {3, 1};
//+
Line(3) = {1, 2};
//+
Line(4) = {2, 4};

Point(5) = {0.2, 0.2, 0, 1.0};
Point(6) = {0.8, 0.2, 0, 1.0};
Point(7) = {0.2, 0.8, 0, 1.0};
Point(8) = {0.8, 0.8, 0, 1.0};//+//+
Line(5) = {7, 8};
//+
Line(6) = {7, 5};
//+
Line(7) = {5, 6};
//+
Line(8) = {8, 6};
//+
Curve Loop(1) = {5, 8, -7, -6};
//+
Plane Surface(1) = {1};
//+
Line(9) = {3, 7};
//+
Line(10) = {4, 8};
//+
Line(11) = {2, 6};
//+
Line(12) = {1, 5};
//+
Curve Loop(2) = {2, 12, -6, -9};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {7, -11, -3, 12};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {8, -11, 4, 10};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {1, 10, -5, -9};
//+
Plane Surface(5) = {5};
//+
Transfinite Curve {1, 2, 9, 10, 4, 11, 12, 3} = 20 Using Progression 1.1;

Transfinite Curve {1, 2, 3, 4} = 20 Using Progression 1;
//+
Transfinite Curve {6, 5, 8, 7} = 20 Using Progression 1;
//+
Transfinite Surface {2};
//+
Transfinite Surface {5};
//+
Transfinite Surface {4};
//+
Transfinite Surface {1};
//+
Transfinite Surface {3};
//+

//+
Physical Curve("wall", 13) = {2, 3, 4};
//+
Physical Curve("top", 14) = {1};
//+
Physical Surface("flow", 15) = {5, 2, 1, 4, 3};
//+
Recombine Surface {1, 5, 2, 4, 3};
