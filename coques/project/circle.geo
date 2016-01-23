Merge "circle.0-Windows\project\circle";
Point(1) = {0, 0, 0, 1.0};
Point(2) = {0, 0, 1.0, 1.0};
Point(3) = {0, 1.1, 0.5, 1.0};
Delete {
  Point{2};
}
Circle(1) = {3, 1, 3};
Physical Line(2) = {1};
Line Loop(3) = {1};
Ruled Surface(4) = {3};
Ruled Surface(5) = {3};
Delete {
  Surface{4};
}
Delete {
  Surface{5};
}
Plane Surface(5) = {3};
Plane Surface(6) = {3};
Delete {
  Surface{5};
}
Delete {
  Surface{6};
}
Ruled Surface(6) = {3};
Delete {
  Surface{6};
}
Plane Surface(6) = {3};
Physical Point(7) = {1};
Physical Point(8) = {3};
Delete {
  Surface{6};
}
Plane Surface(9) = {3};
Delete {
  Surface{9};
}
Plane Surface(9) = {3};
Physical Surface(10) = {9};
