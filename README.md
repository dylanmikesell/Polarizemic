# 3 component seismic polarization analysis

Matt Haney wrote two polarization analysis codes. The first is based on the __covariance__ method. The second is based on the __coherency__ method from [this](http://www.bssaonline.org/content/76/5/1393.abstract) paper.

Run the poliriz.m file to see a comparison of the two methods using a syntethic input. 
 
---

Notes from Matt's original zip file with Dylan's modifications to run based on the Git repo structure.

polariz: a set of Matlab programs to perform polarization analysis on narrowband seismic data

This is a set of 4 Matlab programs - 3 functions and 1 script - that runs an example showing polarization analyses using a real-valued covariance matrix and a complex-valued coherency matrix. For real data, polarization is a function of frequency, so this should be used on data which has been band-pass filtered over a narrow frequency band.

To run the example, simply download the GIT repository, start Matlab, go to that directory in Matlab, and at the Matlab prompt type:

> polariz

This will generate 3 figures showing the results on a synthetic data set. 

The first figure plots the 3 components of motion (vertical, east, and north)  for the synthetic data. The synthetic data consists of 2 parts. The first part is rectilinearly polarized (all 3 components in phase) and the second part is elliptically polarized (vertical component 90 degrees out of phase relative to the horizontal components). In addition to the differences in polarization, the first part has relative amplitudes among the (Z,E,N) components set to (1,2,1). The second part has (Z,E,N) amplitudes set to (1,3,1). All polarization parameters estimated for this synthetic data use a time window of 2 cycles.  

The second figure shows the azimuth angle (in degrees clockwise from north), incidence angle (in degrees from vertical), and ellipticity (ratio of the intermediate axis of motion to the major axis of motion for the 3D ellipsoid) calculated using the complex-valued coherency matrix. For the first part, the rectilinear motion, these parameters should be:

	azimuth [deg] = arctan(2/1)*(180/pi) = 63.4 [deg]
	incidence [deg] = arctan(sqrt(1^2 + 2^2)/1)*(180/pi) = 65.9 [deg]
	ellipticity = 0

For the second part, the elliptical motion, these parameters should be:

	azimuth [deg] = arctan(3/1)*(180/pi) = 71.6 [deg]
	incidence [deg] = 90 [deg]
	ellipticity = 1/sqrt(3^2 + 1^2) = 0.316

Note that the incidence angle of elliptical motion is 90 degrees (horizontal incidence) and the ellipticity of rectilinear motion is 0. As can be seen in the second figure, the method based on the coherency matrix does a good job of estimating these parameters.

The third figure shows the same polarization parameters estimated using the covariance matrix instead of the coherency matrix. The covariance matrix method is seen to be adequate but inferior to the coherency matrix method in estimating the polarization of the synthetic data, in particular for the elliptically polarized part.

Matt Haney
Alaska Volcano Observatory
5/22/2009



