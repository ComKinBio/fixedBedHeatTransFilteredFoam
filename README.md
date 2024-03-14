<div id="top"></div>
<!--
*** README template used
*** https://github.com/othneildrew/Best-README-Template
-->

<!-- PROJECT SHIELDS -->
<!--
*** Markdown "reference style" is used links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->


<!-- PROJECT -->
# fixedBedHeatTransFilteredFoam



<!-- PROJECT LOGO -->
The fixedBedHeatTransFilteredFoam is an extended solver based on the official "[DPMFoam](https://github.com/OpenFOAM/OpenFOAM-7/tree/master/applications/solvers/lagrangian/DPMFoam)" solver with thermal particles. 
<br />
<br />

[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![OpenFOAM v7](https://img.shields.io/badge/OpenFOAM-v7-brightgreen.svg)](https://openfoam.org/)
[![License: GPL v3][license-shield]][license-url]

<div align="center">
  <p align="center">
    <a href="https://github.com/ComKinBio/fixedBedHeatTransFilteredFoam/issues">Report Bug</a>
    ·
    <a href="https://github.com/ComKinBio/fixedBedHeatTransFilteredFoam/issues">Request Feature</a>
  </p>
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#about-the-project">Features</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#Contributing">Contributing</a></li>
    <li><a href="#Contact">Contact</a></li>
    <li><a href="#Publications">Publications</a></li>
  </ol>
</details>



<!-- Features -->
## Features

### Coupling method

IBM (interfaces-based model) and MBM (mesh-based model) are two typical thermally thick particle models. IBM was proposed by Thunman et al. [[3]](#3). In some literature, it also refers as sharp interface model or layer-based model. The IBM used in this repo is based on work of Ström et al. [[2]](#2), and certain modifications are also adopted [[1]](#1). The particle is divided into 4 layers (wet wood, dry wood, char, and ash) by 3 inifinite thin converting fronts (drying, devolatilization, and char burnout). Each layer is assumed to be uniform. The heat and mass transfer between the layers and fronts are calculated. 


### Colliding particle
The new lib has two new templates which are inheritaged from the official Lagrangian lib. The ReactingMultiphaseIBMCloud template is inheritaged from [ReactingMultiphaseCloud](https://github.com/OpenFOAM/OpenFOAM-7/tree/master/src/lagrangian/intermediate/clouds/Templates/ReactingMultiphaseCloud) template, and the ReactingMultiphaseIBMParcel is inheritaged from [ReactingMultiphaseParcel](https://github.com/OpenFOAM/OpenFOAM-7/tree/master/src/lagrangian/intermediate/parcels/Templates/ReactingMultiphaseParcel). 

Two submodels with RTS mechanism are added. The PyrolysisModel (which is modified from [DevolatilisationModel](https://github.com/OpenFOAM/OpenFOAM-7/tree/master/src/lagrangian/intermediate/submodels/ReactingMultiphase/DevolatilisationModel)) and the CharOxidizationModel (which is modified from [SurfaceReactionModel](https://github.com/OpenFOAM/OpenFOAM-7/tree/master/src/lagrangian/intermediate/submodels/ReactingMultiphase/SurfaceReactionModel)) are bundled with the reactingMultiphaseIBMCloud type. The reason to add these two submodels rather than using the original DevolatilisationModel and SurfaceReactionModel is that the pure vitrual functions in the original two has less access variables than wanted. 

### Solver

The case used in the paper is added as a test and tutorial case for this slover.

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- Contributing -->

## Contributing

This repo accepts updating. For example, correcting the coding style to the [OpenFOAM style](https://openfoam.org/dev/coding-style-guide/), adding particle shape submodel to [RTS mechanism](https://openfoamwiki.net/index.php/OpenFOAM_guide/runTimeSelection_mechanism), making submodels of thermally thick particle properties (currently hard coded)... ...

If you have any contribution to this repo, please fork the repo and create a pull request (to dev). You can also simply open an issue with the tag "improvement".

Besides coding, academic discussions through emails are most approciated.



<p align="right">(<a href="#top">back to top</a>)</p>



<!-- LICENSE -->
## License

Distributed under the [GPLv3 License](https://www.gnu.org/licenses/gpl-3.0.en.html). ([OpenFOAM license control](https://openfoam.org/licence/))

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- CONTACT -->
## Contact

Jingyuan Zhang - jingyuan.zhang@ntnu.no 

Tian Li - tian.li@ntnu.no / tian.li@risefr.no

Henrik Ström - henrik.strom@chalmers.se


Research group: [ComKin group at NTNU](https://www.ntnu.edu/comkin/)


<p align="right">(<a href="#top">back to top</a>)</p>

<!-- Publications -->
## Publications

If you want to use fixedBedHeatTransFilteredFoam in your research, you should cite the following papers:

* <a id="1">[1]</a> [Zhang J, Li T, Ström H, et al. A novel coupling method for unresolved CFD-DEM modeling[J]. International Journal of Heat and Mass Transfer, 2023, 203: 123817.](https://www.sciencedirect.com/science/article/pii/S0017931022012856)
 
<p align="right">(<a href="#top">back to top</a>)</p>



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/ComKinBio/fixedBedHeatTransFilteredFoam.svg?style=flat
[contributors-url]: https://github.com/ComKinBio/fixedBedHeatTransFilteredFoam/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/ComKinBio/fixedBedHeatTransFilteredFoam.svg?style=flat
[forks-url]: https://github.com/ComKinBio/fixedBedHeatTransFilteredFoam/network/members
[stars-shield]: https://img.shields.io/github/stars/ComKinBio/fixedBedHeatTransFilteredFoam.svg?style=flat
[stars-url]: https://github.com/ComKinBio/fixedBedHeatTransFilteredFoam/stargazers
[issues-shield]: https://img.shields.io/github/issues/ComKinBio/fixedBedHeatTransFilteredFoam.svg?style=flat
[issues-url]: https://github.com/ComKinBio/fixedBedHeatTransFilteredFoam/issues
[license-shield]: https://img.shields.io/badge/License-GPLv3-blue.svg
[license-url]: https://www.gnu.org/licenses/gpl-3.0

