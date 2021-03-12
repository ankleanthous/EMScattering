// Gmsh project created on Tue Nov 12 10:02:53 2019
SetFactory("OpenCASCADE");
Sphere(1) = {0, 0, 0, 0.5/2};
Dilate {{0, 0, 0}, {1, 1, 1/6}} { Volume{1}; }
Characteristic Length{:} = 0.02;
Mesh.Algorithm = 6;
