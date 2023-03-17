<!-- Title -->
<h1 align="center">
  SubZero
</h1>

<!-- description -->
<p align="center">
  <strong> Sea ice model to explicitly simulate individual floe life cycles using complex discrete elements with time-evolving shapes.</strong>
</p>
<!-- Information badges -->
<p align="center">
  <a href="https://www.repostatus.org/#active">
    <img alt="Repo status" src="https://www.repostatus.org/badges/latest/active.svg?style=flat-square" />
  <a href="https://doi.org/10.5281/zenodo.7199941">
    <img alt="Zenodo" src="https://zenodo.org/badge/doi/10.5281/zenodo.7199941.svg">
  </a>
</p>

This unique model uses parameterizations of floe-scale processes, such as collisions, fractures, ridging, and welding, to by pass resolving intra-floe bonded elements. We demonstrate the novel capabilities of the SubZero model in idealized experiments, including uniaxial compression, the summer-time sea ice flow through the Nares Strait and a winter-time equilibration of floe size and
ice thickness distributions.

For more information see paper: A Sea Ice Model with an Explicit Representation of the Floe Life Cycle by Georgy Manucharyan and Brandon Montemuro. 
This code is also saved on Zenodo

## Installation instructions

1. Clone SubZero contents to their own folder. Open folder within MATLAB.

2. Once within the MATLAB project, open the MATLAB command window.

3. Run command ```cd private``` to enter `private` folder within the command window.

4. Run command ```mex('-D__int64=__int64_t','clipper.cpp','mexclipper.cpp')``` in MATLAB command window. This will link the C++ code used for polygon clipping within the program. For more information on this linking see [here](https://www.mathworks.com/help/matlab/ref/mex.html).

5. If this does not work, you may have to install Xcode to compile the mex file on a Mac. Xcode is a large file and is available in App Store and online. Downloading the installer from the apple website may be quicker.

6. After installing Xcode, you may have to run these two lines in the terminal: 
    + ```sudo xcode-select -s /Applications/Xcode.app/Contents/Developer```
    + ```sudo xcodebuild -license accept```
   
 7. Shut down MATLAB, re-open, and re-do step 4.
 
 ## Running your first model

Model runs out of the box. Currently set up model creates ~10 floes that fill the entire domain with collisions and corner-griding on, and all other functionality turned off. Model will run for 7500 timesteps where each timestep is 10 seconds.

To change the model, the following steps can be taken:

1. Define the boundaries of the domain `initialize_boundaries.m` file. Three options are present, with two commented out. Choose one of the given options, or write your own.

2. Set up desired ocean forcings in `initialize_ocean.m` file. Specifically see the `psi_ocean`, `Uocn`, and `Vocn` variables. Note that your ocean must be bigger than the boundaries of the domain, in some cases significantly bigger to prevent errors with periodic boundaries. Additionally, the user is responsible for making the ocean doubly periodic outside of domain bounds. 

3. Define the atmospheric winds within the `SubZero.m` file. Specifically see the `winds.u` and `winds.v` variables.

4. Define the concentration and number of floes you wish to start with in the SubZero.m file. The concentration is the fraction of the domain that will initally be covered in ice. For this, see the variable `target_concentration` around line 57. This variable is repurposed below for concentration when packing (re-set to 1), so the initial value must be set for model initialization. The number of floes is the fourth argument in the following line:
  `[Floe, Nb] = initial_concentration(c2_boundary,target_concentration,height,10,min_floe_size)` - here 10 floes.

5. Within `initial_concentration`, topography can be added to the simulation. These will be represented as immovable, unbreakable floes. `Nb` is the number of these boundaries and their shape can be input here as well. See Nares Strait example within the `validation_cases` folder to see how to add these elements.

5. Set the flags for desired physical processes to `true` as well as the frequency with which you want them to execute in `SubZero.m` file. The flags are at the beginning of the file, while the frequencies are lower throughout the file. Look for lines of the form `if FRACTURES && mod(i_step,75)==0` with different flags. These lines give the frequency of each of these flag operations. In the example above, fractures happen every 75 timesteps.

7. Set the timestep `dt` directly under the flags.

6. Set the total number of snapshots of the system as well as how often that will be output in `SubZero.m` file. These are set under the `Set output parameters` comment.

7. Run the `SubZero.m` program to start a simulation. 

To see examples of different scenarios, see the examples included in `validation_cases` folder as tar files. These are the validation cases from the Subzero paper.  

## Validation Cases  
### 1. Uniaxial Compression:
Here we demonstrate the behavior of sea ice floes subject to uniaxial compression in a confined domain. The run is initialized with 200 floes in a fully-packed domain with the North/South boundaries moving towards the center of the domain, and stationary East/West boundaries. A relatively small time step, dt=5 s, is used to resolve the elastic waves in response to external boundary motion and changes in the floe configuration due to fractures. The atmospheric and oceanic stresses are set to zero for this simplified test. The floes are subject to Mohr-Coulomb fracture criteria, but there is no floe simplification, corner grinding, welding, ridging, rafting, or creation of new floes in this scenario. The boundaries move with a constant prescribed velocity, v_b = 0.1 m/s.

1. In `Subzero.m' set:
  1a. Flags to false except FRACTURES and COLLISION to true. 
  1b. Set dt = 5;
  1c. The variables `height.mean' is set to 1 and `height.delta' should be set to 0. 
  1d. Atospheric winds should be U0 = 0 and V0 = 0. 
  1e. Target_concentration = 1
  1f. The number of initial floes to 200
  1g. The modulus = 2.5e3*(mean(sqrt(cat(1,Floe).area))+min(sqrt(cat(1,Floe.area)))
  1h. Set all instances of doInt.flag = false so there are no forces from ocean/atmosphere and all movement is driven by boundaries.
  1i. The frequency for fractures is set to mod(i_step,200)
  1j. After updating time and i_step add the following lines to have the boundary move.
    if max(c2_boundary(2,:)) > 85000 && mod(i_step,30)==0
        xb = c2_boundary(1,:);
        yb = c2_boundary(2,:);
        yb = yb - 15*[-1 1 1 -1];
        c2_boundary = [xb; yb];
        Ly = max(c2_boundary(2,:));Lx = max(c2_boundary(1,:));
        c2_boundary_poly = polyshape(c2_boundary');
        c2_border = scale(c2_boundary_poly,2); c2_border = subtract(c2_border, c2_boundary_poly);
        floebound = initialize_floe_values(c2_border, height);
      end
  
2. In `initialize_boundaries.m' the box is set such that Lx=Ly=10 kilometers.
  
3. In `fracture.m' set Sig11 = 1.5e5 and make sure all the elliptical yield curve lines are commented out or deleted.  

4. In `floe_interactions.m' set mu = 0.3.
  
### 2. Nares Strait:
The Nares Strait simulation demonstrates the role of floe fractures in wind-driven sea ice transport through narrow straits. Nares Strait is a channel between Ellesmere Island (Canada) and Greenland. The simulation aims to reflect spring or summer-like conditions of Arctic sea ice export through Nares Strait after the breakup of its winter arches. Since the transport events are relatively short (order of days or less), the effects of thermodynamic sea ice melt could be considered secondary relative to mechanical floe processes such as collisions and fractures. We thus randomly initialize the model with relatively large floes of uniform thickness, covering only the area just north of the strait. The uniform 10 m/s southward winds generate stresses that push the floes through the strait, while the ocean is assumed to be stagnant. Coastal boundaries are prescribed using a series of static floes. All physical processes except collisions and fractures are turned off to model the spring/summer breakup of floes. To modify the program to run the Nares validation case on your own, the following changes must be made.
  
1. In `Subzero.m' set flags to false except FRACTURES and COLLISION to true. The variables `height.mean' is set to 1 and `height.delta' should be set to 0. Atospheric winds should be U0 = 0 and V0 = -10. Target_concentration = [1; 0], and the number of initial floes to 150. The frequency for fractures is set to mod(i_step,150).

2. In `initialize_ocean.m' the ocean grid is set such that dXo=20000 meters, Lx=Ly =2e6 meters, and ocean velocity field is set to zero.
  
3. In `initialize_boundaries.m' the box is set such that Lx=50 kilometers, the top y boundary is 500 kilometers and the bottom y boundary is -250 kilometers.
  
4. In `initialize_concentration.m' 
 4a. Add the following lines after `dy = y(2) -y(1)' to create the boundaries of Nares Strait: (Note that you will need to include the file `Nares_Strait_segments.mat' in a folder named Nares. This mat file can be found in the Nares validation zip file)
  
  Lx = abs(min(x));Ly = abs(min(y));
  Bx = [-Lx -Lx Lx    Lx];
  By = [-Ly  -1.5e5 -1.5e5  -Ly];
  B = polyshape(Bx',By');  

  %Create floes that act as boundaries and wont move
  load( './Nares/Nares_Strait_segments','R')

  %Remove islands
  if ~ISLANDS
    R(1:4) = [];
  end

  Floe = []; bound = c2_boundary_poly;
  for ii = 1:length(R)
    FloeNEW = initialize_floe_values(R(ii),height);
    bound = subtract(bound, FloeNEW.poly);
    Floe = [Floe FloeNEW];
  end
  N = length(Floe);
  for ii = 1:N
    Floe(ii).poly = polyshape(Floe(ii).c_alpha'+[Floe(ii).Xi Floe(ii).Yi]);
  end

  Bu = union([Floe.poly]);
  Bu = union(Bu,B);

 4b. Then, after `[~, b,~,~,~] = polybnd_voronoi([X Y],boundary.Vertices)' insert  the following lines to remove any overlap between the shapes generated by the Voronoi tesselation and Nares Strait topography.

  for kk = 1:length(b)
    poly(kk) = polyshape(b{kk});
  end
  poly = intersect([poly],polyB);
  polyout=subtract([poly],Bu);
  clear polyIce; polyIce = [];
  for kk = 1:length(polyout)
    R = regions(polyout(kk));
    polyIce = [polyIce; R];
  end
  
 4c. Replace the line `Nf = 1:length(b);' with the following:
  Nf = 1:length(polyIce);%randperm(length(b));
  
 4d. Replace the while loop with the following since we have already created the polyshapes;
  while Atot/area(polyB)<=c(jj,ii)
      floenew = initialize_floe_values(polyIce(Nf(count)),height);
      Floe = [Floe floenew];
      count = count+1;
      Atot = Atot+floenew.area;
      if count > length(Nf)
          Atot = area(polyB)+1;
      end
   end
  
 4d. Finally, after the lines `areas = cat(1,Floe.area)' add `areas(1:Nb) = min_floe_size;'

5. In `floe_interactions_all.m' update to include c2_boundary such that it reads 
  [tmp,Fx,Fy] =calc_trajectory(dt,ocean,winds,Floe(i),HFo,doInt,c2_boundary);
  
6. In `calc_trajectory.m':
  6a. The first line needs changed so the function can take in the addtional input and should read:
  function [floe,FxOA,FyOA] =calc_trajectory(dt,ocean,winds,floe,HFo, doInt,c2_boundary)
  
  6b. After calculating the new floe.c_alpha add the following if statement:
        if min(floe.c_alpha(2,:))+floe.Yi<min(c2_boundary(2,:))
            floe.alive = 0;
        end
                                                                
 7. In `fracture.m' set Pstar = 1e5 and make sure all the Mohr's cone lines are commented out or deleted.                                                    
 8. In `floe_interactions.m' set mu = 0.25. 
                                                                
### 3. Winter ITD and FSD equilibration:
Here we demonstrate an essential case of model equilibration in winter-like conditions, where all parameterizations are active. We subject sea ice to strong mechanical and thermodynamic forcing over a five week period to facilitate an accelerated model evolution away from the initialized floe shapes, sizes, and thicknesses towards typical winter-like distributions. Specifically, we prescribe idealized ice-ocean stresses in the form of four equal-strength counter-rotating gyres (arranged like mechanical gears) that create relative sea ice motion and facilitate floe fractures and ridging. Alternatively, one could prescribe atmosphere-ocean stresses to achieve the same goal, but in this run the winds are set to 0. To make this a winter-like simulation, we ensured continuous sea ice growth by specifying a fixed negative heat flux that increases the thickness of existing ice floes, the formation of new ice floes in open ocean regions, and welding between floes.
  
## Community Guidelines

If you experience issues using Subzero, please open a GitHub issue. If you're interested in contributing, please reach out to the paper authors. 
