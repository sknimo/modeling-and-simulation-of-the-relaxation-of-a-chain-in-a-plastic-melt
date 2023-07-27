# Relaxation of a plastic melt during processing [Simulation]

![Image plastic](https://adamjeedurabuilt.com/wp-content/uploads/2021/08/010-800x313.jpg)

<sub><sup>https://adamjeedurabuilt.com/wp-content/uploads/2021/08/010-800x313.jpg</sup></sub>


## Introduction
---
>Modeling and simulation  are an integral part of engineering. It helps us test ideas and concepts before time and resources are investment in them. Eventhough commercial softwares are available to help with that, presently none of the software is geared towards polymer science. In the course of my study, I had to develop all the softwares that I required for my research. This current one is one of them. My hope is that, this initiative should end up being an open source project with contributions from others in the field.

This programs evaluates the manner in which a chain found in a plastic melt relaxes when the imposed stress is removed. During processing of plastics, products of desires size and properties are manufactured so as to meet certain engineering needs. When the plastic melt leaves the extruder, it has a tendency to swell, since the barrel from which it was extruded no longer imposes constrains on it. A clear understanding and evaluation of the relaxation the relaxation of the melt is therefore a must to produce materials of specific size and shapes. With this program, the relaxation modulus of a melt can be estimated. This program is an implementation of the theoretical works of [Pattamaprom et al(2000)](https://link.springer.com/article/10.1007/s003970000104). 

## Dataset
---
>The dataset was digitized from a paper published by [*kwakye-Nimo et al.(2022)*](https://pubs.acs.org/doi/abs/10.1021/acs.macromol.2c01102).

### Project structure

* the folder _reports_ contains, the executive summary
* all codes can be found in the folder _src_
* the file _notes_on_numerical_implementations.md_ gives details on all the numerical steps that were taken to develop this program

## Key insights
---
<!-- The subpopulations in a molecular weight distribution can be identified using Machine learning. With knowledge of the sub-population characteristics an optimization routine can be subsequently used to effectively perform peak deconvolution. This method eliminates the need for the user to know ahead of time the characteristics of the peaks he is seeking to identify. -->


<!-- ![Image result](reports/figures/deconvoluted.png) -->

## Tools used in this project
---
* Packages: Anaconda
* Libraries: Pandas, Numpy, Scipy, Matplotlib
* Programming languages: Python 3