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
  
### 2. Nares Strait:
The Nares Strait simulation demonstrates the role of floe fractures in wind-driven sea ice transport through narrow straits. Nares Strait is a channel between Ellesmere Island (Canada) and Greenland. The simulation aims to reflect spring or summer-like conditions of Arctic sea ice export through Nares Strait after the breakup of its winter arches. Since the transport events are relatively short (order of days or less), the effects of thermodynamic sea ice melt could be considered secondary relative to mechanical floe processes such as collisions and fractures. We thus randomly initialize the model with relatively large floes of uniform thickness, covering only the area just north of the strait. The uniform 10 m/s southward winds generate stresses that push the floes through the strait, while the ocean is assumed to be stagnant. Coastal boundaries are prescribed using a series of static floes. All physical processes except collisions and fractures are turned off to model the spring/summer breakup of floes.
  
### 3. Winter ITD and FSD equilibration:
Here we demonstrate an essential case of model equilibration in winter-like conditions, where all parameterizations are active. We subject sea ice to strong mechanical and thermodynamic forcing over a five week period to facilitate an accelerated model evolution away from the initialized floe shapes, sizes, and thicknesses towards typical winter-like distributions. Specifically, we prescribe idealized ice-ocean stresses in the form of four equal-strength counter-rotating gyres (arranged like mechanical gears) that create relative sea ice motion and facilitate floe fractures and ridging. Alternatively, one could prescribe atmosphere-ocean stresses to achieve the same goal, but in this run the winds are set to 0. To make this a winter-like simulation, we ensured continuous sea ice growth by specifying a fixed negative heat flux that increases the thickness of existing ice floes, the formation of new ice floes in open ocean regions, and welding between floes.
  
## Community Guidelines

If you experience issues using Subzero, please open a GitHub issue. If you're interested in contributing, please reach out to the paper authors. 
