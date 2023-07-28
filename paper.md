---
title: 'SubZero: a discrete element sea ice model that simulates floes as evolving concave polygons'
tags:
  - Sea Ice Modeling
  - Collisions, Fractures, Deformation
  - Discrete Element Methods
  - Deformable Polygonal Elements 
  - Floe Size Distribution
  - Ice Thickness Distribution
  
authors:
  - name: Brandon P. Montemuro
    orcid: 0000-0003-1946-4916
    equal-contrib: true
    affiliation: 1
  - name: Georgy E. Manucharyan
    orcid: 0000-0001-7959-2675
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 1
affiliations:
 - name: School of Oceanography, University of Washington, Seattle, WA, USA
   index: 1

date: 05 October 2022
bibliography: paper.bib

---

# Summary

SubZero is a conceptually new discrete element sea ice model geared to explicitly simulate the life cycles of individual floes by using polygonal elements with time-evolving boundaries. This unique model uses parameterizations of floe-scale processes, such as collisions, rafting, ridging, fracturing, and welding, to simulate the behavior on sea ice floes subject to mechanical and thermodynamic forcing in confined or periodic domains. SubZero enables the exploration of a wide range of floe interaction rules and fracture criteria to further our understanding of sea ice mechanics, including distributions of floe sizes, thicknesses, and shapes. \autoref{fig:Nares} shows snapshots from a validation study where SubZero simulates floes moving through a channel. A complete model description and example process studies demonstrating its capabilities can be found in @Manucharyan:2022. SubZero was developed using MATLAB [@MATLAB:2020], leveraging its built-in functions for polygonal operations. The model source code has been archived to Zenodo [@Montemuro:2022]. 

![The evolution of sea ice floes as they pass through the Nares Strait (between Canada and Greenland), including (a) initial floe state with the inset showing the location of Nares Strait, (b) floes shortly after sea ice breakup that occurred after about three days, and (c) floe state after ten days when many floes have passed through the Nares Strait. The initial distribution of floes was generated using a Voronoi tessellation, and the subsequent evolution of floe shapes is only subject to floe fractures. The green box in panel (a) shows the modeled part of Nares Strait. The blue arrows represent sea ice velocity after averaging floe momentum on an Eulerian grid.\label{fig:Nares}](Nares_Floes.png)



# Statement of need

Sea ice dynamics span a wide range of scales. At length scales O(10--100) km and smaller, sea ice exhibits granular behavior as individual floes and fracture networks become evident [@Rothrock:1984; @Zhang:2015; @Stern:2018]. Sea ice motion at relatively large scales, O(100 km), is commonly represented in climate models using continuous rheological models, like the viscous-plastic model [@Hibler:1979], which have not been formally derived from basic sea ice physics but have been postulated instead. The rheology defines a relationship between sea ice stress caused by floe-floe or floe-lead interactions, to the large scale deformation of ice cover, the material properties of sea ice, and the state of the ice cover. As such, the rheological models carry potentially large uncertainties in representing sea ice dynamics. They are also not designed to represent the scales of motion at which individual floes start to influence dynamics [@Coon:2007]. 

Developed initially in the context of granular assembles and rock dynamics [@Cundall:1979; @Potyondy:2004], Discrete Element Models (DEMs) are an alternative to continuous rheology models. DEMs can be computationally demanding because they represent media as a collection of many colliding bonded elements with specified shapes and contact laws. DEMs resort to setting the interaction laws between their elements and strive to calibrate them using micro-scale observations because the continuous equations of motion are often unknown. Existing floe-scale sea ice DEMs use bonded elements of simple preset shapes like disks [@Herman:2013; @Damsgaard:2018; @Chen:2021], polygons [@Kulchitsky:2017], or tetrahedra [@Liu:2018] to represent complex floe geometries. However, floe-scale modeling remains challenging due to difficulties reconciling discrete elements' idealized nature with complex floe-scale observations. Observations indicate that ice floes vary dramatically in size and shape and change over time due to various processes like fractures, rafting/ridging, lateral growth/melt, welding, etc. Therefore, using prescribed element shapes creates ambiguity about what elements and bonds between them represent. It is difficult to search for direct correspondence between the state variables of existing sea ice DEMs and real observations without a clear understanding of what an individual DEM element represents. 

SubZero was designed as an alternative to continuous rheology models and existing sea ice DEMs to improve the realism of sea ice simulations at floe scales. In contrast with existing sea ice models that use elements of pre-defined simple shapes, SubZero uses concave polygonal elements that are free to evolve in complexity in response to parameterized floe-scale processes. The model's capability of developing floe shapes naturally might bring us closer to direct model validation using floe-scale observations and advance our understanding of sea ice physics through idealized process studies. 




# Acknowledgements

B.P.M. and G.E.M. gratefully acknowledge support from the Office of Naval Research (ONR) grant N00014-19-1-2421. The authors highly appreciate the insightful discussions at the online workshop ``Modeling the Granular Nature of Sea Ice'' organized by the School of Oceanography, University of Washington, as part of the ONR grant N00014-19-1-2421. The authors thank Skylar Gering for reviewing the manuscript and the SubZero source code.

# References
