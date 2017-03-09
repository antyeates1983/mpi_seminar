% Read results from C++ or FORTRAN binary file
% and plot as a function of time.

close all; clear all

fid = fopen('u.dat','rb');
nx = fread(fid, 1, 'int32');
x = fread(fid, nx, 'float64');
u = fread(fid, nx, 'float64');
fclose(fid);

figure()
plot(x,u,'o-');

hold on;

fid = fopen('../fortran/u.dat','rb');
nx = fread(fid, 1, 'int32');
x = fread(fid, nx, 'float64');
u = fread(fid, nx, 'float64');
fclose(fid);
plot(x,u,'sr');

