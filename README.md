# Artery Scaling Law Boundary Conditions

Scaling law boundary conditions for OpenFoam 2.1.1 based off the paper 
[Functional and anatomical measures for outflow boundary conditions in atherosclerotic coronary bifurcations](https://www.researchgate.net/publication/285045581_Functional_and_anatomical_measures_for_outflow_boundary_conditions_in_atherosclerotic_coronary_bifurcations).

A representative inlet flow rate for a case can be found by using a diameter flow
relation,

<img src="https://render.githubusercontent.com/render/math?math=q_{in}=1.43d^{2.55}">

For two diameters a scaling law for the set of outflow ratios can be used 

<img src="https://render.githubusercontent.com/render/math?math=\dfrac{q_2}{q_1}=\left(\dfrac{d_2}{d_1}\right)^{2.27}">

Using mass conservation we know that

<img src="https://render.githubusercontent.com/render/math?math=q_{1} %2B q_{2}=q_{in}">

Thus,

<img src="https://render.githubusercontent.com/render/math?math=q_{2}=q_{1}\left(\dfrac{d_2}{d_1}\right)^{2.27}">

and,

<img src="https://render.githubusercontent.com/render/math?math=q_{1}=\dfrac{q_{in}}{1 %2B \left(\dfrac{d_2}{d_1}\right)^{2.27}}">

Thus we can build our boundary conditions to satisfy  the above equations

## Building the Boundary Conditions

To install all boundary conditions source the OpenFOAM 2.1.1 .bashrc and run

``` shell
# Download the github repository
git clone https://github.com/Chr1sC0de/ArteryScalingLawsBC.git
# Provide permission to all the ./Allwmake and ./Allwclean files
chmod -R +x ArteryScalingLawsBC
cd  ArteryScalingLawsBC
# Install the boundary conditions
./Allwmake
```

## Running Examples

To run a single example on the gadi super computer

``` shell
cd examples/<case_folder>
qsub openfoam_job.sh
```

To change the number of processors used for all cases

``` shell
cd examples
# The following code will change the ncpus variable in the openfoam_job.sh
# As well as the numberOfSubdomains is the systems/decomposeParDict file
./change_ncpus.sh -n <number of cpus>
```

## Cleaning the Boundary Conditions

To uninstall all the boundary conditions 

``` shell
cd  ArteryScalingLawsBC
./Allwclean
```