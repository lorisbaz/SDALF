SDALF
=====
In this package, you find a updated version of the MATLAB code for following paper:

**Person Re-identification by Symmetry-Driven Accumulation of Local Features**
M. Farenzena, L. Bazzani, A. Perina, V. Murino, and M. Cristani
*In IEEE Conference on Computer Vision and Pattern Recognition (CVPR), 2010*

Three demos have been provided
- L_main_tagging_SvsS.m : performs single-shot cross-validation on the given dataset, in order to emulate the 
		results reported on our paper. If you want to run the script from the beginning, 
		you can just change the experiment number ("expnum = 'XYZ'").
- L_main_tagging_MvsM.m : performs multi-shot cross-validation on the given dataset, in order to emulate the 
		results reported on our paper. 
- In the folder CMC_curves, we put the CMCs of our paper. So you can use them directly if you are 
		interested in making comparisons.
NOTE: when performing the tests on the ETHZ datasets, remember to set SUBfac = 0.5, because the original 
test of the PLS paper have been performed on half resolution.
		

##ADDITIONAL LIBS
We also provide with our package some additional libraries we used in our method. 
- Maximally Stable Colour Regions (addLibs/mscr) http://www.cs.ubc.ca/~perfo/software/
- CVlab (addLibs/CVlab_0.2) http://profs.sci.univr.it/~fusiello/teaching/visione/
- Local Binary Pattern (addLibs/LBP) http://www.scholarpedia.org/article/Local_Binary_Patterns
- GMM clustering (addLibs/SandG)


##COMPILE
Please note, some libraries contain mex-files that needs to be compiled for your machine. We provide 
a limited set of binary within the package. If you do not find the binary for your machine, please 
compile the source code (instructions in the external libraries folder "addLibs").
In particular, MSCR are partially in mex and also the functions detect_mscr_masked and whistcY in 
folder. whistcY does not need any dependence. For detect_mscr_masked please follow the same instruction 
for compiling the MSCR.


##BIBTEX AND ACKNOWLEDGMENTS
	@inproceedings{Farenzena2010CVPR,
	author = {Farenzena, Michela and Bazzani, Loris and Perina, Alessandro and Murino, Vittorio and Cristani, Marco},
	title = {Person Re-Identification by Symmetry-Driven Accumulation of Local Features},
	booktitle = {Proceedings of the 2010 IEEE Computer Society Conference on Computer Vision and Pattern Recognition},
	year = {2010},
	publisher = {IEEE Computer Society},
	address = {San Francisco, CA, USA},
	}

This research is founded by the EU-Project FP7 SAMURAI, grant FP7-SEC- 2007-01 No. 217899.
The code was written by Loris Bazzani and Michela Farenzena. For any problems contact 
the first author of the paper.
