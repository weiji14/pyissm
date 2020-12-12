steps=[4];

if any(steps==1)   %Mesh Generation #1
  %Mesh parameters
  domain =['pig_domain.exp'];
  hinit=10000;   % element size for the initial mesh
  hmax=20000;    % maximum element size of the final mesh (originally 40000)
  hmin=250;     % minimum element size of the final mesh (originally 5000)
  gradation=1.7; % maximum size ratio between two neighboring elements
  err=8;         % maximum error between interpolated and control field

  % Generate an initial uniform mesh (resolution = hinit m)
  md=bamg(model,'domain',domain,'hmax',hinit);

  % Load Velocities
  nsidc_vel='Data/antarctic_ice_vel_phase_map_v01.nc';

  % Get necessary data to build up the velocity grid
  % xmin    = ncreadatt(nsidc_vel,'/','xmin');
  % ymax    = ncreadatt(nsidc_vel,'/','ymax');
  % spacing = ncreadatt(nsidc_vel,'/','spacing');
  % nx      = double(ncreadatt(nsidc_vel,'/','nx'));
  % ny      = double(ncreadatt(nsidc_vel,'/','ny'));
  % vx      = double(ncread(nsidc_vel,'vx'));
  % vy      = double(ncread(nsidc_vel,'vy'));
  %
  % % Read coordinates
  % xmin = strtrim(xmin);
  % xmin = str2num(xmin(1:end-2));
  % ymax = strtrim(ymax);
  % ymax = str2num(ymax(1:end-2));
  % spacing = strtrim(spacing);
  % spacing = str2num(spacing(1:end-2));
  %
  % % Build the coordinates
  % x=xmin+(0:1:nx)'*spacing;
  % y=(ymax-ny*spacing)+(0:1:ny)'*spacing;

  % % Interpolate velocities onto coarse mesh
  % vx_obs=InterpFromGridToMesh(x,y,flipud(vx'),md.mesh.x,md.mesh.y,0);
  % vy_obs=InterpFromGridToMesh(x,y,flipud(vy'),md.mesh.x,md.mesh.y,0);
  % vel_obs=sqrt(vx_obs.^2+vy_obs.^2);
  % clear vx vy x y;

  x = ncread(nsidc_vel,'x');
  y = flipud(ncread(nsidc_vel,'y'));
  vx = flipud(ncread(nsidc_vel, 'VX')');
  vy = flipud(ncread(nsidc_vel, 'VY')');

  % Interpolate velocities onto coarse mesh
  vx_obs=InterpFromGridToMesh(x,y,vx,md.mesh.x,md.mesh.y,0);
  vy_obs=InterpFromGridToMesh(x,y,vy,md.mesh.x,md.mesh.y,0);
  vel_obs=sqrt(vx_obs.^2+vy_obs.^2);
  clear vx vy x y;

  % Adapt the mesh to minimize error in velocity interpolation
  md=bamg(md,'hmax',hmax,'hmin',hmin,'gradation',gradation,'field',vel_obs,'err',err);

  %ploting
  plotmodel(md,'data','mesh')

  % Save model
  save ./Models/PIG_Mesh_generation md;
end

if any(steps==2)  %Masks #2

  md = loadmodel('Models/PIG_Mesh_generation');

  % Load SeaRISe dataset for Antarctica
  % http://websrv.cs.umt.edu/isis/index.php/Present_Day_Antarctica
  % searise='../Data/Antarctica_5km_withshelves_v0.75.nc';

  % Load ALBMAP dataset
  % https://doi.pangaea.de/10.1594/PANGAEA.734145
  albmap='Data/ALBMAPv1.nc';

  %read grounded/floating mask from ALBMAP
  x1=double(ncread(albmap,'x1'));
  y1=double(ncread(albmap,'y1'));
  mask=double(ncread(albmap,'mask'));

  %interpolate onto our mesh vertices
  groundedice=double(InterpFromGridToMesh(x1,y1,mask',md.mesh.x,md.mesh.y,0));
  groundedice(groundedice<=0)=-1;
  clear mask;

  %fill in the md.mask structure
  md.mask.ocean_levelset=groundedice; %ice is grounded for mask equal one
  md.mask.ice_levelset=-1*ones(md.mesh.numberofvertices,1);%ice is present when negative

  %ploting
  plotmodel(md,'data',md.mask.ocean_levelset,'title','grounded/floating','data',md.mask.ice_levelset,'title','ice/no-ice')

  % Save model
  save ./Models/PIG_SetMask md;
end

if any(steps==3)  %Parameterization #3
  for bedname = {'DeepBedMap','BedMachine'}
    fprintf('%s\n',bedname{1});
    md=loadmodel('./Models/PIG_SetMask');

    % DeepBedMap v1.1 https://doi.org/10.5281/zenodo.4054246
    if strcmp(bedname{1},'DeepBedMap')
      ncdata0='Data/deepbedmap_dem.nc';
      x     = single(ncread(ncdata0,'x'));
      y     = flipud(single(ncread(ncdata0,'y')));
      bed   = ncread(ncdata0,'z')';
    % BedMachine v2 https://doi.org/10.5067/E1QL9HFQ7A8M
    elseif strcmp(bedname{1},'BedMachine')
      ncdata1='Data/BedMachineAntarctica_2020-07-15_v02.nc';
      x     = single(ncread(ncdata1,'x'));
      y     = flipud(single(ncread(ncdata1,'y')));
      bed   = flipud(ncread(ncdata1,'bed')');
    end

    md.geometry.base=InterpFromGridToMesh(x,y,bed,md.mesh.x,md.mesh.y,0);
    md.friction = frictionschoof();  % Set to use Schoof (2005) type sliding law
    md=parameterize(md,'Pig/Pig.par');

    % Save model
    modelfilename=strcat('./Models/PIG_Parameterization_',bedname);
    save(modelfilename{1},'md');
  end
end

if any(steps==4)  %Control Method #4
  for bedname = {'BedMachine','DeepBedMap'}
    modelfilename=strcat('./Models/PIG_Parameterization_',bedname{1});
  	md=loadmodel(modelfilename);
    md.miscellaneous.name=strcat('pig_',bedname{1},'_control_drag_fullstokes');

    % Extrude mesh
    md=extrude(md,10,3); % 10 layers, with extrusion exponent of 3
    % md.friction=frictionschoof(md.friction); % Ensure still using Schoof-type

  	% Control general
  	md.inversion.iscontrol=1;
    % M1QN3 optimizer parameters
  	md.inversion.maxsteps=30; % maximum number of iterations (gradient computation)
  	md.inversion.maxiter=40; % maximum number of Function evaluation
  	md.inversion.dxmin=0.1; % convergence criterion: two points less than dxmin from each other (sup-norm) are considered identical
  	md.inversion.gttol=1.0e-4; % gradient relative convergence criterion 2
  	md.verbose=verbose('control',true);

  	% Cost functions
  	md.inversion.cost_functions=[101 103 501];
  	md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,3);
  	md.inversion.cost_functions_coefficients(:,1)=1;
  	md.inversion.cost_functions_coefficients(:,2)=100;
  	md.inversion.cost_functions_coefficients(:,3)=1e-7;

  	% Controls
  	md.inversion.control_parameters={'FrictionC'};
  	md.inversion.min_parameters=1*ones(md.mesh.numberofvertices,1);
  	md.inversion.max_parameters=200^2*ones(md.mesh.numberofvertices,1);

  	% Additional parameters
  	md.stressbalance.restol=0.01;
  	md.stressbalance.reltol=0.1;
  	md.stressbalance.abstol=NaN;

    % Set Flow Equation
    % https://issm.jpl.nasa.gov/documentation/approximations/
    % md.debug.valgrind=1;
    md = setflowequation(md,'FS','all');

  	% Solve
  	md.toolkits=toolkits;
  	md.cluster=generic('name',oshostname,'np',4,'executionpath',fullfile(pwd,'Models'));
  	md=solve(md,'Stressbalance');

  	% Update model friction fields accordingly
  	md.friction.coefficient=md.results.StressbalanceSolution.FrictionC;
    md.friction.effective_pressure=md.results.StressbalanceSolution.Pressure;

  	plotmodel(md,'data',md.friction.coefficient,'data',md.materials.rheology_B,'data',md.friction.effective_pressure)

  	% Save model
    modelfilename=strcat('./Models/',md.miscellaneous.name);
    save(modelfilename{1},'md');
  end
end

if any(steps==5) %Plot #5

  md = loadmodel('./Models/PIG_Control_drag');

  plotmodel(md,...
  'data',md.initialization.vel,'title','Observed velocity',...
  'data',md.results.StressbalanceSolution.Vel,'title','Modeled Velocity',...
  'data',md.geometry.base,'title','Bed elevation',...
  'data',md.results.StressbalanceSolution.FrictionCoefficient,'title','Friction Coefficient',...
  'colorbar#all','on','colorbartitle#1-2','(m/yr)',...
  'caxis#1-2',([1.5,4000]),...
  'colorbartitle#3','(m)', 'log#1-2',10);
end

if any(steps==6)  %Higher-Order #6

  % Load Model
  md = loadmodel('./Models/PIG_Control_drag');

  % Disable inversion
  md.inversion.iscontrol = 0

  % Extrude Mesh
  md = extrude(md,3,1)

  % Set Flowequation
  md = setflowequation(md,'FS','all')

  % Solve
  md=solve(md,'Stressbalance')

  % Save Model
  save ./Models/PIG_ModelHO md

end

if any(steps==7)  %Plot #7

  mdHO = loadmodel('./Models/PIG_ModelHO');
  mdSSA = loadmodel('./Models/PIG_Control_drag');

  basal=find(mdHO.mesh.vertexonbase);
  surf=find(mdHO.mesh.vertexonsurface);

  plotmodel(mdHO,'nlines',3,'ncols',2,'axis#all','equal',...
  'data',mdHO.initialization.vel,'title','Observed velocity',...
  'data',(mdHO.results.StressbalanceSolution.Vel(surf)-mdHO.initialization.vel(surf)),'title','(HO-observed) velocities',...
  'data',mdSSA.results.StressbalanceSolution.Vel,'title','Modeled SSA Velocity',...
  'data',(mdHO.results.StressbalanceSolution.Vel(surf)-mdSSA.results.StressbalanceSolution.Vel),'title','(HO-SSA) velocities',...
  'data',mdHO.results.StressbalanceSolution.Vel,'title','Modeled HO surface Velocities',...
  'data',(mdHO.results.StressbalanceSolution.Vel(surf)-mdHO.results.StressbalanceSolution.Vel(basal)),'title','(HOsurf-HO base) velocities',...
  'caxis#1',([1.5,4000]),'caxis#3',([1.5,4000]),'caxis#5',([1.5,4000]),...
  'colorbar#all','on','view#all',2,...
  'colorbartitle#all','(m/yr)',...
  'layer#5',1, 'log#1', 10,'log#3', 10,'log#5', 10);
end
