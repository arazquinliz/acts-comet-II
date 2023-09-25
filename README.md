# acts-comet-II
Implementation of the [ACTS](https://github.com/acts-project/acts) tracking software into the [COMET](https://comet.kek.jp/) Phase-II experiment.

This folder contains the required modules to perform track reconstruction using [ACTS](https://acts.readthedocs.io/en/v28.2.0/index.html) in the context of the COMET Phase-II experiment. The files starting with uppercase are essential for the reconstruction process, files starting with "track" are programs that contain the essential structure to perform the process in the name. How to run them is explained inside each file. Files starting with "final" are the completed programs to perform seeding, fitting, finding and reconstruction with realistic COMET data. The files containing the actual data are not given for privacy purposes. Finally, files ending in ```.c``` are analysis tools that can give comprehensive plots from the output data of the aforementioned programs.

A [helloworld.cpp](https://github.com/arazquinliz/acts-comet-II/blob/main/helloworld.cpp) program is present to make sure that all the modules work fine. 

## Usage
Compile the code:
```g++ helloworld.cpp Geometry.cpp MagField.cpp Measurements.cpp Seeder.cpp 
        Logger.cpp KalmanFitting.cpp Writer.cpp Estimator.cpp -o helloworld.x 
        -L/path/to/acts_install/lib -lActsPluginTGeo -lActsPluginJson -lActsCore 
        -I/path/to/acts_install/include/ -I$(root-config --incdir) 
        -I$(root-config --evelibs) $(root-config --libs) 
        -lGeom -I/usr/include/eigen3
```

And run the code:
```./helloworld.x```

Each program requires input files in the form of ```.root``` files. Please do look into the required files and create them as needed.

## Requirements
ACTS 28.2.0 (and dependencies)

ROOT 6.28/04 or newer.

## Note
ACTS is a fast-changing toolkit, the code here presented might have some syntax differences from the current version. Some little code updates might be required to execute the programs.

## LICENSE 
Apache License 2.0

----------------------------------
Date: September 25 2023
