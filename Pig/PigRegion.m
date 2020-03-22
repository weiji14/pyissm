%Getting the velocity in PIG vicinity for the ExpDraw

% Load Velocities
% http://nsidc.org/data/nsidc-0484.html
nsidc_vel='../Data/antarctic_ice_vel_phase_map_v01.nc';

% % Get necessary data to build up the velocity grid
% xmin    = ncreadatt(nsidc_vel,'/','xmin');
% ymax    = ncreadatt(nsidc_vel,'/','ymax');
% spacing = ncreadatt(nsidc_vel,'/','spacing');
%
% nx = double(ncreadatt(nsidc_vel,'/','nx'));
% ny = double(ncreadatt(nsidc_vel,'/','ny'));
% vx = double(ncread(nsidc_vel,'vx'));
% vy = double(ncread(nsidc_vel,'vy'));
%
% xmin = strtrim(xmin);  % this is a string, and we need to recover the double value
% xmin = str2num(xmin(1:end-2));  % get rid of the unit and convert to double
%
% ymax = strtrim(ymax);
% ymax = str2num(ymax(1:end-2));
%
% spacing = strtrim(spacing);
% spacing = str2num(spacing(1:end-2));
%
% % Build the coordinates
% x=xmin+(0:1:nx)'*spacing;
% y=(ymax)-(0:1:ny)'*spacing;

x = ncread(nsidc_vel,'x');
y = flipud(ncread(nsidc_vel,'y'));
vx = flipud(ncread(nsidc_vel, 'VX')');
vy = flipud(ncread(nsidc_vel, 'VY')');

%Limit the region to Pine Island
posx  = find(x>=-18.0e5 & x<=-12.0e5);
posx  = find(x>=-30.0e5 & x<=30.0e5);
x_pig = x(posx);
posy  = find(y>-4.0e5 & y<=1.0e5);
posy  = find(y>-30.0e5 & y<=20.0e5);
y_pig = y(posy);

vx_pig  = vx(:,:);
vy_pig  = vy(:,:);
vel_pig = sqrt(vx_pig.^2+vy_pig.^2);

imagesc(x_pig,y_pig,log(vel_pig+1));
axis xy equal tight;
