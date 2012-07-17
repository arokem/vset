The Visual System Engineering Toolbox (VSET) is a Matlab toolbox designed for calculating the properties of the front end of the visual system.  The code in VSET is a portion of the general toolbox, Image Systems Engineering Toolbox, that is sold commercially by Imageval Consulting, LLC.  That code is designed for industry to help design novel image sensor. The VSET portion of the ISET code is being freely distributed for use in modeling biological properties of image formation.

The VSET toolbox is written around several data structures.  Each of these is implemented using set/get/create syntax.  Major computations are implemented by computing with these data structures.

The scene data structure describes the scene radiance (photons).  It is set up to permit depth encoding, though nearly all of the current examples are based on a radiance field originating from a single plane.

The optical image transforms the scene radiance into the irradiance distribution at the sensor.  There are several optics models that implement the transformation either as diffraction limited, shift invariant, or ray trace (shift variant) depending on the level of information you have about the optics.  There is an implementation of the optics of the human eye, based on work from Marimont and Wandell.  This toolbox also works seamlessly with Wavefront Toolbox (also in github) that uses data from adaptive optics to simulate the defocus based on measurements of many different human eyes.

The sensor transforms the irradiance into a spatial array of cone absorptions. The sensor pixels and spectral quantum efficiency can be set to model the human cones at various eccentricities and with various types of inert pigments (macular pigment, lens density, optical density).  The sensor simulates the photon absorptions in the cone (and rod) mosaics.

For a general introduction to vision, please see:

  http://foundationsofvision.stanford.edu/

For examples of the code and some tutorials on human vision, please see

  https://www.stanford.edu/group/vista/cgi-bin/FOV/computational-examples/
