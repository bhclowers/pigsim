pigsim_v2 Readme

Authors: 
Cullen Greer, Brian Clowers 
Dept. of Chemistry, Washington State University
Pullman, Washington

https://github.com/bhclowers/pigsim

Example Directories:
   To use the example directories "pigsim Tetramer example" or "pigsim Pickins example", simply open SIMION and choose "Run Lua Script" in the main view panel. In the dialog box that opens, navigate to where you've saved the example directory, and from that select "Extrude_SLIM_Setup.lua". This should open a small terminal window while SIMION generates the necessary potential arrays. For the tetramer example this should not exceed a few minutes of processing, though the Pickins example is substantially longer, and generation time will depend on computing power (I walk away from my machine).

pigsim.ipynb:
	To interact with the pigsim.ipynb (or iPython notebook), you need to install the jupyter Python plugin, as well as the dependencies of the notebook, though they are relatively few. The best way is through anaconda. Our preferred method is to install miniconda from their website: 
	<https://docs.conda.io/en/latest/miniconda.html>

	From there, create a conda environment in your terminal with:
		> conda create --name SLIM_Pickins python=3.10

	This environment also needs:
	numpy (1.23.5)     ## specific version probably not important, but this is what's in our implementation
	matplotlib (3.7.0)
	pillow (9.4.0)

	Activate your new environment with:
		> conda activate SLIM_Pickins

	and install these packages using:
		> conda install numpy=1.23.5
		etc.

	Finally, install the package for creating an ipynb kernel:
		> conda install -c anaconda ipykernel

	and create the kernel:
		> python -m ipykernel install --user --name=SLIM_Pickins
		# sometimes this fails, try running it with 'python3' instead of 'python'


	You should now be able to open a jupyter notebook in the directory of the .ipynb file:
		> jupyter notebook

	This should open a browser window where the pigsim.ipynb file should be interactive. Select it, and in the top menu, change the kernel to the one you've just created. You shuold now have a proper kernel to import the necessary packages and run the code blocks.

Be mindful of where you keep the scripts, running the ipynb generates a few small files. Running the setup lua generates several large files. The Pickins example requires approximately 6 GB of drive space to populate the full ion optics workbench. Additionally, the setup lua relies on the Required/ folder, so don't move that or edit anything in it. 