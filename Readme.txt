Download SubZero contents to their own folder.

Within Matlab enter the private folder and run mex('-D__int64=__int64_t','clipper.cpp','mexclipper.cpp') in MATLAB command window.

*** How to set up and run SubZero ***

Define the boundaries of the domain initialize_boundaries.m file.
Set up desired ocean forcings in initialize_ocean.m file.
Define the atmospheric winds within the SubZero.m file.
Define the concentration and number of floes you wish to start the run with in SubZero.m file.
Set the flags for desired physical processes to true as well as the frequency with which you want them to execute in SubZero.m file.
Set the total number of snapshots of the system as well as how often that will be output in SubZero.m file.
Run the SubZero.m program to start a simulation. 

*** Note ***
You may have to install Xcode to compile the mex file on a Mac.

Xcode is a large file and is available in App Store and online. Downloading the installer from the apple website may be quicker .

After installing xcode, you may have to run these two lines in the terminal:
sudo xcode-select -s /Applications/Xcode.app/Contents/Developer
sudo xcodebuild -license accept

The validations from SubZero: A Sea Ice Model with an Explicit Representation of the Floe Life Cycle are included in validation_cases folder as tar files
