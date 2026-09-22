# FEMTool

This is the repository of the code used for the course on domain decomposition.
It provides basic infrastructure for creating a multi-dimensional P1-finite 
element code in C++. It assumes that the following are available on your machine:
- c++ compiler (up-to-date, covering c++23 features)
- make (compilation tool)
- [Vizir4](https://pyamg.saclay.inria.fr/vizir4.html) for vizualization
- [Gmsh](https://gmsh.info/) for mesh generation

An example file showcasing a typical usecase if located 
in the directory "example". To compile/execute this example:
- open a terminal and move to the directory "example"
- compile by typing make
- run the executable ./example

The example code should produce two output files named "output.mesh"
and "output.sol". You can vizualize them with vizir4.

# Organisation du Code
Les fichiers tpX-exoY sont chacun dans le répértoire tp/tpX/exoY 
Le code de l'exo 2 du TP3 est dans le fichier preconditioner.hpp et le code de l'exo3 TP1 dans le fichier iterativesolver.hpp
J'ai des doublons de fichiers .zone.identifier provenants de la  copie du projet ddm dans ma machine distante.
# Code visualisation
 Pour visulaliser les courbes j'ai utilisé un Jupiter Notebook  visu_projet2.ipynb python qui lit les fichiers .dat 

## Author
[Xavier Claeys](https://claeys.pages.math.cnrs.fr/)

## License
This software is under [LGPL](https://en.wikipedia.org/wiki/GNU_Lesser_General_Public_License)

