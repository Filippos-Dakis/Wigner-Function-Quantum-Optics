Filippos Tzimkas-Dakis   Virginia Tech  November 2023

Any feedback and suggestions are much appreciated! 
    
    ----->  dakisfilippos@vt.edu  <-------
         
		 
		 
This file contains two main scripts:
	- CoherentBasis.m
	- FockBasis.m
which define two types of classes. The first one, CoherentBasis.m , refers to Coherent states basis 
of Qunatum Harmonic Oscillator (QHO), while the second, FockBasis.m, refers to the number/Fock states 
of QHO.

The file also includes two examples, one for each case. 

The key feature of these two classes (CoherentBasis, and FockBasis) is that they can compute the 
Wigner Distributio (Quasi-Probability distribution) for any given state. However, they also include 
other useful features such as state addition, normalization, Displacement operations etc.

This is the very first version and I intend to incorporate more features. Some of them are 
	- Q distribution 
	- time evolution 
	- plot number states in harmonic potential

I open to other suggestion. Please free to contact me. 


                ****  OPEN QUESTION TO COMMUNITY ****
The code slowes down when you are calculating the Wigner Function in the FockBasis. The main reason behind this
is the calculation of the Displacement operator in the fock basis. If you have any suggestions on that please
contact me! 