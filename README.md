# Design-of-an-evaporator-for-a-buoyant-two-phases-flow-cooling-loop-of-a-telecommunication-satellite

<p align="justify">Satellites are constituted of electronic components that release heat flux. The heat released must be evacuated otherwise there are some risks to damage components in the satellite. In order to guarantee the correct working of the satellite, a thermal loop crossed by a refrigerant fluid must be added. The loop is constituted by a heater, a pressure regulator, a condenser, and a pump. The evaporator is the element of the loop in which the heat released by the electronic components is transferred to the fluid. The fluid, receiver of the components entailed energy, will then be heated and will evaporate itself little by little. This breeds the appearance of a two phase flow. The program, developed with the MATLAB code takes into account the different models of two phase flow that we can consider in our case. We are then able to follow the evolution of the physical variables which characterize the system like the components temperature, the drop loss, void fraction and the quality.</p>


## Getting Started

<p align="justify">These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.</p>

### Prerequisites

<p align="justify">The code being implemented in MATLAB, it requires the <a href="https://www.mathworks.com/products/matlab.html">MATLAB</a> software, which is licensed under the <a href="https://www.mathworks.com/">MathWorks</a>sofware company. MATLAB can be downloaded at the following link <a href="https://www.mathworks.com/downloads/">https://www.mathworks.com/downloads/</a>, and more information can be found about the license fee on the website.</p>


## File descriptions

* '.f90' files in which the main code, as well as the different subroutines are programmed.
* 'physical_data.txt' which contains the different parameters to define the domain of computation, mesh size, and other parameters such as CFL or Fourier numbers.
* In the output directory 'ex_output_files', there are four files:     
-> 2 output *.vts files at t=0 et Tf/2    
-> 2 output *.txt files of velocity profile in x=1 for different data 

* 4 animations in the 'animations' directory.

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
