Download SubZero contents to their own folder.

Within Matlab enter the private folder and run mex('-D__int64=__int64_t','clipper.cpp','mexclipper.cpp') in MATLAB command window.

Run the SubZero.m program to start a simulation. 

To toggle any of the physical process set flags to true or false.

***Note ***
You may have to install Xcode to compile the mex file on a Mac.

Xcode is a large file and is available in App Store and online. Downloading the installer from the apple website may be quicker .

After installing xcode, you may have to run these two lines in the terminal:
sudo xcode-select -s /Applications/Xcode.app/Contents/Developer
sudo xcodebuild -license accept