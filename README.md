# Panel Method: A Basic Application for Hand-Drawn Airfoils
<h2>About</h2>

This simple Python program has the objective to calculate common airfoil parameters by just setting a few inputs and drawing the geometry of the body in a dynamic way. There is certain freedom for the user to draw any airfoil (or arbitrary figure) by connecting dots, interpolating and later applying the Panel Method theory to do the calculations.

<h2>Setting the Airfoil</h2>
<h3>Input Parameters</h3>

The user must set four aerodynamic and geometric parameters for the calculation to be done. This is a very straightforward action, done by clicking the **SET** buttons as seen below. The default values are unit chord and unit velocity, zero angle of attack and air density at sea level.

![Program Image #1](https://i.ibb.co/gTv9GhB/1.png)

<h3>Drawing / Loading Airfoil</h3>

By clicking **DRAW**, the user is able to draw any arbitrary airfoil by connecting dots. However, some things should be noted in order for the calculation to return meaningful output:

![Program Image #2](https://i.ibb.co/YtFnfgy/2.png)

The program itself is going to use a spline method to interpolate the dots and make a much smoother airfoil.

Notes:
* Airfoil coordinates are in relation to the chord.
* The gray rectangle is used just for reference.

By clicking **LOAD**, the user is able to load several NACA airfoils. New files can be uploaded to the “SampleAirfoils” folders. The “x” and “y” coordinates in these files are also set in clockwise direction, with the last node being the same as the first (the trailing edge).

<h2>Output Parameters</h2>

Once the input parameters and airfoil are set, the user can calculate the output parameters by clicking the **Click to Calculate** button.
Some common parameters such as lift and moment coefficients are seen in the main window, while three plots are available for visualization (velocity stream plot, pressure coefficient over airfoil, and pressure coefficient contour plot).

For example:

![Program Image #3](https://i.ibb.co/8P5HNjQ/3.png)

<h2>References</h3>
This program uses the Source and Vortex Panel Method which can be found in John D. Anderson's Fundamentals of Aerodynamics.

The geometric parametrization of the airfoil was based on JoshTheEngineer's [Panel Method series](https://www.youtube.com/playlist?list=PLxT-itJ3HGuUDVMuWKBxyoY8Dm9O9qstP) on YouTube.
