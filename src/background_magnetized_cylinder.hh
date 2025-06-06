/*!
\file background_magnetized_cylinder.hh
\brief Declares a background of an infinite cylinder magnetized perpendicular to its axis
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_MAGNETIZED_CYLINDER_HH
#define SPECTRUM_BACKGROUND_MAGNETIZED_CYLINDER_HH

#include "background_cylindrical_obstacle.hh"

namespace Spectrum {

//! Method for computing derivatives of B (0: analytical, 1: numerical)
#define MAGNETIZED_CYLINDER_DERIVATIVE_METHOD 0

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundMagnetizedCylinder class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the class
const std::string bg_name_magnetized_cylinder = "BackgroundMagnetizedCylinder";

/*!
\brief Magnetic field of a uniformly magnmetized cylinder, inside and outside
\author Vladimir Florinski
\author Juan G Alonso Guzman

Parameters: (BackgroundCylindricalObstacle)
*/
class BackgroundMagnetizedCylinder : public BackgroundCylindricalObstacle {

protected:

//! Compute the internal u, B, and E fields
   void EvaluateBackground(void) override;

//! Compute the internal u, B, and E derivatives
   void EvaluateBackgroundDerivatives(void) override;

public:

//! Default constructor
   BackgroundMagnetizedCylinder(void);

//! Copy constructor
   BackgroundMagnetizedCylinder(const BackgroundMagnetizedCylinder& other);

//! Destructor
   ~BackgroundMagnetizedCylinder() override = default;

//! Clone function
   CloneFunctionBackground(BackgroundMagnetizedCylinder);
};

};

#endif
