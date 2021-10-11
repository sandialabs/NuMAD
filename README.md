![](media/NuMAD.png) 
# Numerical Manufacturing And Design (NuMAD) Tool 
NuMAD (Numerical Manufacturing And Design) is an object-oriented, open-source software program written in Matlab which simplifies the process of creating a three-dimensional model of a wind turbine blade. The tool manages all blade information including aerodynamic, and material properties as well as, material placement. The blade information can be modified by a structural-optimizer in NuMAD or it can be used to run other tools in a stand-alone mode:

    The ANSYS® commercial finite element package
    The NWTC codes PreComp, BModes, and FAST
    PLOT3D file format for CFD mesh building

When used in conjunction with these other tools, NuMAD provides computation of blade cross section properties, various structural analyses, and estimation of blade aeroelastic instability. 

For any questions or to request more information, send an e-mail to: numad@sandia.gov


## Installation Instructions
1. Install [Matlab](https://www.mathworks.com/products/matlab.html)
2. Open `paths.m` in an editor and set the `NuMAD_path` to the `source` directory in the repo. Save and close `paths.m`.
3. Move the `paths.m` file to the MATLAB directory of your user account (e.g. C:\Users\username\Documents\MATLAB)

Every time a new MATLAB session is started type `paths` to be able to run NuMAD
### Models

Example YMAL files can be located in the [examples](https://github.com/sandialabs/NuMAD/examples) folder on GitHub. 

### Contributing, reporting bugs, and requesting help

To report issues or bugs please [create a new
issue](https://github.com/sandialabs/NuMAD/issues/new) on GitHub.

We welcome contributions from the community in form of bug fixes, feature
enhancements, documentation updates, etc. All contributions are processed
through pull-requests on GitHub. Please follow our [contributing
guidelines](https://github.com/sandialabs/NuMAD/blob/master/CONTRIBUTING.md)
when submitting pull-requests.
  
## License

NuMAD is licensed under BSD 3-clause license. Please see the
[LICENSE](https://github.com/sandialabs/NuMAD/LICENSE.md) included in
the source code repository for more details.

## Acknowledgements 

NuMAD is currently being developed with funding from Department of Energy's
(DOE) Energy Efficiency and Renewable Energy (EERE) Wind Energy Technology Office (WETO). 
