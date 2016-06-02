# oblimap-2.0
The second OBLIMAP release (2016)

Abstract. 
This paper accompanies the second OBLIMAP open source release. The package is developed to map climate fields between a general circulation model (GCM) and an ice sheet model (ISM) in both directions by using optimal aligned oblique projections, which minimize distortions. Both grids may be irregularly spaced and the ratio of the grids is allowed to differ largely. The stand-alone version of OBLIMAP is a powerful tool to map various differently gridded datasets on one uniform ISM grid with an optimal centered projection. This is demonstrated by relevant applicatons concerning the major ice caps. As this applies also for the reverse mapping direction, it can be used as an off-line coupler. Besides, OBLIMAP 2.0 is an embeddable GCM–ISM coupler, suited for high frequent on-line coupled experiments. A new fast scan method is presented as an alternative for the former time consuming grid search strategy, realising a performance gain of several orders of magnitude and enabling the mapping of high resolution datasets with a much larger number of grid nodes. Further a highly flexible masked mapping option is added. The limitations of the fast scan method with respect to unstructured and adaptive grids are discussed together with several proposed parallel implementations in order to achieve another performance gain.

See http://www.geosci-model-dev-discuss.net/gmd-2016-124/
