# oblimap-2.0
The second OBLIMAP release (2016)

This software is open source (See licensing details elsewhere). In addition to the formal licensing terms, we would greatly appreciate an acknowledgement. Preferably in the form of a citation and a link to the web-page.

Citation: Reerink et al. (2016,2010):

Reerink, T. J., van de Berg, W.J., van de Wal, R. S. W.  (2016), Fast gcm – ice sheet model coupling software oblimap 2.0, including on-line embeddable mapping routines, Geoscientific Model Development, 2016, 1–31, doi: 10.5194/gmd-2016-124, http://www.geosci-model-dev-discuss.net/gmd-2016-124/

Reerink, T. J., Kliphuis, M. A., van de Wal, R. S. W.: Mapping technique of climate fields between GCM’s and ice models, Geoscientific Model Development, 3, 13–41, doi:10.5194/gmd-3-13-2010, http://www.geosci-model-dev.net/3/13/2010/, 2010.


The OBLIMAP User Guide can be cited as follows:

Reerink, T. J.: OBLIMAP User Guide, version 1.0, accompanying OBLIMAP 2.0, Tech. rep., Institute for Marine and Atmospheric research Utrecht, Utrecht University, 3508 TA Utrecht, The Netherlands, https://github.com/oblimap/oblimap-2.0/blob/master/documentation/oblimap-user-guide.pdf, 2016.


Abstract. 
This paper accompanies the second OBLIMAP open source release. The package is developed to map climate fields between a general circulation model (GCM) and an ice sheet model (ISM) in both directions by using optimal aligned oblique projections, which minimize distortions. The curvature of the surfaces of the GCM and ISM grid differ, both grids may be irregularly spaced and the ratio of the grids is allowed to differ largely. OBLIMAP's stand-alone version is able to map data sets which differ in various aspects on the same ISM grid. Each grid may either coincide with the surface of a sphere, an ellipsoid or a flat plane, while the grid types might differ. Reprojection of e.g. ISM data sets is also facilitated. This is demonstrated by relevant applications concerning the major ice caps. As the stand-alone version applies also for the reverse mapping direction, it can be used as an off-line coupler. Besides, OBLIMAP 2.0 is an embeddable GCM - ISM coupler, suited for high-frequency on-line coupled experiments. A new fast scan method is presented for structured grids as an alternative for the former time consuming grid search strategy, realising a performance gain of several orders of magnitude and enabling the mapping of high resolution data sets with a much larger number of grid nodes. Further a highly flexible masked mapping option is added. The limitation of the fast scan method with respect to unstructured and adaptive grids is discussed together with a possible future parallel MPI implementation. 

See http://www.geosci-model-dev-discuss.net/gmd-2016-124/
