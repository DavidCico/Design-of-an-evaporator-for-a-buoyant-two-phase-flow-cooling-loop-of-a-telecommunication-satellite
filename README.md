# Design-of-an-evaporator-for-a-buoyant-two-phases-flow-cooling-loop-of-a-telecommunication-satellite

<p align="justify">Satellites are constituted of electronic components that release heat flux. The heat released must be evacuated otherwise there are some risks to damage components in the satellite. In order to guarantee the correct working of the satellite, a thermal loop crossed by a refrigerant fluid must be added. The loop is constituted by a heater, a pressure regulator, a condenser, and a pump. The evaporator is the element of the loop in which the heat released by the electronic components is transferred to the fluid. The fluid, receiver of the components entailed energy, will then be heated and will evaporate itself little by little. This breeds the appearance of a two phase flow. The program, developed with the MATLAB code takes into account the different models of two phase flow that we can consider in our case. We are then able to follow the evolution of the physical variables which characterize the system like the components temperature, the drop loss, void fraction and the quality.</p>


## Getting Started

<p align="justify">These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.</p>

### Prerequisites

<p align="justify">The code being implemented in MATLAB, it requires the <a href="https://www.mathworks.com/products/matlab.html">MATLAB</a> software, which is licensed under the <a href="https://www.mathworks.com/">MathWorks</a> sofware company. MATLAB can be downloaded at the following link <a href="https://www.mathworks.com/downloads/">https://www.mathworks.com/downloads/</a>, and more information can be found about the license fee on the website.</p>


### File descriptions

<ul>
<li>"Conductance.m" (function that returns the equivalent conduction of our fin)</li> 
<li>"R123-R245FA.xls" (Excel file containing the coolant properties)</li> 
<li>"R245FA.m" (script loading the coolant properties)</li>
<li>"Rg_initial.m" (function that calculates Rg at an instant t)</li> 
<li>"calcul_Rg_comp.m" (script calling ode45 to solve the equations on Rg, x and P)</li> 
<li>"elbow_geometry.txt" (text file with geometry data, used by the code)</li>
<li>"elbow_geometry.xlsx" (Excel file containing the geometry data)</li> 
<li>"main.m" (main script)</li> 
<li>"pressure_Awad.m" (script on Awad approximation for pressure loss)</li>
<li>"pressure_Baroczy.m" (script on Baroczy approximation for pressure loss)</li>
<li>"pressure_L_M.m" (script on Lockhart & Martinelli approximation for pressure loss)</li>   
<li>"skin_friction.m" (function that returns the X parameter associated with skin friction)</li>
<li>"temp_G_W.m" (function that calculates temperature with Gunger & Winterton model)</li>
<li>"temp_K.m" (subroutine that calculates temperature with Kandlikar model)</li>
<li>"temp_S_G.m" (script that calculates temperature with Schrock & Grossman model)</li>
<li>"write_results.m" (script writing results of the model with gravity in a .txt)</li>
<li>"write_results_no_g.m" (script writing results of the model with no gravity in a .txt)</li>
</ul>

#### Geometry files

<p align="justify">There are normally two files describing the geometry of the case : a .txt and a .xlsx (excel 2007+). If modification is needed, then first the excel file should be modified and then be saved in .txt, in order to keep the MATLAB tabulations. Furthermore, the "." must be used for describing floating numbers. The geometry has been separated in different elements, each of these having thier own particularity (component, elbow, nothing). Here is the meaning of each column:</p>
 
	- initial heigt/position of the element (in mm)
	- final heigh of the element (in mm)
	- elbow (1 if the block is one, 0 else)
	- surface energy flux (if a component is present)
	- gravity value (0 in an elbow)
	- component power (if there is a component)
	- length of the element (in m)
	- component width (in m)
	- number of passages of the tube inside component (if a component is present)

<ul>
<li><div align="justify">"Bibliography" directory contains 4 research papers in pdf format for some background on the equations implemented inside the code.</div></li>
<li>"Example_plots" contains 2 figures that are shown below.</li>
<li><div align="justify">"Report.pdf" is the report of the project explaining the different equations used in the MATLAB code, and showing some of the results obtained under different conditions.</div></li>
</ul>


### Running the program

<p align="justify">To launch the program, it just requires to execute the script <b>main.m</b> on MATLAB. This is the main script of the code that regroups all the other subroutines, get the different data, plot the wanted curves, and write results in txt files. For good performance of the code, and if the geometry data and coolant properties need to be modified, then the same order of magnitude must be kept for the different parameters.</p>

<p align="justify">After starting the program, a choice will be asked to the user, to know which type of simulation he is willing to run:</p>
<ul>
	<li>2 with and without gravity</li>
	<li>1 with gravity only</li>
	<li>0 with no gravity</li>
</ul>

<p align="justify"The different results will be plotted and computed with MATLAB, and exported in .txt format: "results_with_gravity.txt" and "results_without_gravity.txt".</p>

<p align="justify">To compare the pressure loss results in the evaporator with different models, the scripts named "pressure_****.m" needs to be executed. However, this first requires running the main program once in order to get the data on the evolution of the title x (vapor mass fraction), and the geometry discretisation. Moreover, if a graphic comparison of the different models is required, the subroutine "pressure_L_M.m" must be executed before the others.</p>

<p align="justify">Below, are two of the figures that can be generated by the code. We can see here for instance, the evolution of both vapour mass fraction and velocity against the position inside the evaporator's tube. The reader is referred to the <a href="https://github.com/DavidCico/Design-of-an-evaporator-for-a-buoyant-two-phases-flow-cooling-loop-of-a-telecommunication-satellite/blob/master/Report.pdf">project report</a> for more information on the different results obtained through the numerical simulations.</p>

<p align="center">
<img src="https://github.com/DavidCico/Design-of-an-evaporator-for-a-buoyant-two-phases-flow-cooling-loop-of-a-telecommunication-satellite/blob/master/Example_plots/Vapour_fraction_tube.jpg" width="500" height="500"> <img src="https://github.com/DavidCico/Design-of-an-evaporator-for-a-buoyant-two-phases-flow-cooling-loop-of-a-telecommunication-satellite/blob/master/Example_plots/Velocity_comparison_tube.jpg" width="500" height="500" >
</p>


## Contributing

Please read [CONTRIBUTING.md](https://github.com/DavidCico/Design-of-an-evaporator-for-a-buoyant-two-phase-flow-cooling-loop-of-a-telecommunication-satellite/blob/master/CONTRIBUTING.md) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **David Cicoria** - *Initial work* - [DavidCico](https://github.com/DavidCico)

See also the list of [contributors](https://github.com/DavidCico/Design-of-an-evaporator-for-a-buoyant-two-phase-flow-cooling-loop-of-a-telecommunication-satellite/graphs/contributors) who participated in this project.
