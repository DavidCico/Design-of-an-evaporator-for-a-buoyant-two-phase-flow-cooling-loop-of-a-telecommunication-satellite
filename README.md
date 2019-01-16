# Design-of-an-evaporator-for-a-buoyant-two-phases-flow-cooling-loop-of-a-telecommunication-satellite

<p align="justify">Satellites are constituted of electronic components that release heat flux. The heat released must be evacuated otherwise there are some risks to damage components in the satellite. In order to guarantee the correct working of the satellite, a thermal loop crossed by a refrigerant fluid must be added. The loop is constituted by a heater, a pressure regulator, a condenser, and a pump. The evaporator is the element of the loop in which the heat released by the electronic components is transferred to the fluid. The fluid, receiver of the components entailed energy, will then be heated and will evaporate itself little by little. This breeds the appearance of a two phase flow. The program, developed with the MATLAB code takes into account the different models of two phase flow that we can consider in our case. We are then able to follow the evolution of the physical variables which characterize the system like the components temperature, the drop loss, void fraction and the quality.</p>


## Getting Started

<p align="justify">These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.</p>

### Prerequisites

<p align="justify">The code being implemented in MATLAB, it requires the <a href="https://www.mathworks.com/products/matlab.html">MATLAB</a> software, which is licensed under the <a href="https://www.mathworks.com/">MathWorks</a>sofware company. MATLAB can be downloaded at the following link <a href="https://www.mathworks.com/downloads/">https://www.mathworks.com/downloads/</a>, and more information can be found about the license fee on the website.</p>


## File descriptions

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

	* Geometry file
	There are normally two files describing the geometry of the case : a .txt and a .xlsx (excel 2007+). If modification
	is needed, then first the excel file should be modified and then be saved in .txt, in order to keep the MATLAB tabulations.
	Furthermore, the "." must be used for describing floating numbers.
	The geometry has been separated in different elements, each of these having thier own particularity (component, elbow, nothing)
 	Here is the meaning of each column :
		- initial heigt/position of the element (en mm)
		- final heigh of the element (en mm)
		- elbow (1 if the block is one, 0 else)
		- surface energy flux (if a component is present)
		- gravity value (0 in an elbow)
		- component power (if there is a component)
		- length of the element (en m)
		- component width (en m)
		- number of passages of the tube inside component (if a component is present)


	* Coolant data file
	Coolant data file is in .xls (Excel) format, and contains R245FA and R143 fluid properties
	for liquid and gaseous phase. The description of each column is inside the file itself

### Running the program

1. Input numerical values in the file physical_data.txt

        80            ! n  number of mesh cells in y    
        70            ! m  number of mesh cells in x    
        5             ! Length L1 of the domain    
        2             ! Height L2 (Left side of the domain)    
        2             ! Height L3 (Right side of the domain)     
        0.01          ! Diffusion coefficient    
        0             ! Velocity U (x direction)     
        0             ! Velocity V (y direction)          
        5             ! Final time    
        0.9           ! CFL number    
        0.4           ! Fourier number    

    Modifying and tuning these values in order to have an orthogonal mesh or not, diffusion and/or advection...

2. Use the **Makefile** to compile all the files and create the executable (run the command 'make' while being in the main directory of the program).

3. Launch the executable , which will create the mesh and run the discretised calculation on the latter.

4. Observe the concentration field on the domain using Paraview (open sol.pvd).

5. Remove the created files thanks to the commands 'make clean' and 'make solclean'.
