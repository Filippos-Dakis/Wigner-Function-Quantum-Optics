Filippos Tzimkas-Dakis   Virginia Tech  February 2024   

		 
All the scripts were developed on MATLAB 2023a (compatible with MATLAB 2022b or later).

This file contains two main scripts:
- CoherentBasis.m
- FockBasis.m

which define two types of classes. The first one, CoherentBasis.m , refers to Coherent states basis of Qunatum Harmonic Oscillator (QHO), while the second, FockBasis.m, refers to the number/Fock states 
of QHO.

The file also includes four examples
- CoherentBasis_Example_1.m
- FockBasis_Example_1.m
- Q_function_CoherentBasis_Example.m
- Q_function_FockBasis_Example.m
- Example_Squeezed_states.m

that help the user explore the features of the two main classes. 

The key feature of these two classes (CoherentBasis, and FockBasis) is the calculation of the **Wigner Quasi-Probability distribution** and the **Q-Husimi Distributions** for any given state! However, they also include 
other useful features such as quantum state addition, normalization, Displacement operations, annihilation and creation operators, etc.

For instance, **CoherentBasis_Example_1.m** produces the following Wigner Quasi-Probability distributions

![Wigner_Functions](https://github.com/Filippos-Dakis/Wigner-Function-Quantum-Optics/assets/114699564/686d66a3-1eba-4f42-acd4-ea8bdee7f206)

and **FockBasis_Example_1.m** produces these three Wigner Distributions

![Wigner_Functions_2](https://github.com/Filippos-Dakis/Wigner-Function-Quantum-Optics/assets/114699564/e9af7821-4a9c-4f5f-9c70-3f23885cf1e0)

Script **Q_function_CoherentBasis_Example.m** produces the following figure

![Q_function](https://github.com/Filippos-Dakis/Wigner-Function-Quantum-Optics/assets/114699564/a9cbdb41-c013-4341-b634-c2c1ac460a17)

A similar figure is produced by Q_function_FockBasis_Example.m

This is the very first version and I intend to incorporate more features. Some of them are  
- density matrix representation 
- qauntum Zeno gate for Cat-qubits
- switch between the two bases (connect the two classes)
 
Suggested textbooks:
- Exploring the Quantum Atoms, Cavities and Photons.  Serge Haroche Jean-Michel Raimond
- Introductory Quantum Optics. Christopher Gerry Peter Knight

I am open to other suggestion. Please feel free to contact me at dakisfilippos@vt.edu


If you found my code useful, please cite is as  https://github.com/Filippos-Dakis/Wigner-Function-Quantum-Optics


Any feedback and suggestions are much appreciated! 

More features are on the way! Stay tuned !!
