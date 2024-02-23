/*!
\file background_solarwind.hh
\brief Implements a plasma background class for the constant speed supersonic wind of a rotating star
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_solarwind_termshock.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundSolarWind methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 02/22/2023
*/
BackgroundSolarWindTermShock::BackgroundSolarWindTermShock(void)
                            : BackgroundSolarWind(bg_name_solarwind_termshock, 0, MODEL_STATIC)
{
};

/*!
\author Juan G Alonso Guzman
\date 02/22/2023
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBackground()" with the argument of "true".
*/
BackgroundSolarWindTermShock::BackgroundSolarWindTermShock(const BackgroundSolarWindTermShock& other)
                            : BackgroundSolarWind(other)
{
   RAISE_BITS(_status, MODEL_STATIC);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBackground(true);
};

/*!
\author Juan G Alonso Guzman
\date 02/22/2023
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void BackgroundSolarWindTermShock::SetupBackground(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) BackgroundSolarWind::SetupBackground(false);
   container.Read(&r_TS);
   container.Read(&w_TS);
   container.Read(&s_TS);

   s_TS_inv = 1.0 / s_TS;
   dmax_TS = dmax_fraction * w_TS;
};

/*!
\author Juan G Alonso Guzman
\date 02/22/2023
*/
void BackgroundSolarWindTermShock::EvaluateBackground(void)
{
// First compute regular Parker Spiral fields
   BackgroundSolarWind::EvaluateBackground();

// Scale U (radial) and phi component of B if position is past the termination shock. The radial component of B and E should be the same.
   posprime.ChangeFromBasis(eprime);
   double r = r_ref * posprime.Norm();
   GeoVector Br, Bp, rhat = UnitVec(posprime);
   if(r > r_TS) {
      Br = (_spdata.Bvec * rhat) * rhat;
      Bp = (_spdata.Bvec - Br);
      if(r > r_TS + w_TS) {
         _spdata.Uvec *= s_TS_inv * Sqr((r_TS + w_TS) / r);
         _spdata.Bvec = Br + s_TS * Bp * Sqr(r / (r_TS + w_TS));
      }
      else {
         _spdata.Uvec *= 1.0 + (s_TS_inv - 1.0) * (r - r_TS) / w_TS;
         _spdata.Bvec = Br + Bp * (1.0 + (s_TS - 1.0) * (r - r_TS) / w_TS);
      };
   };
};

/*!
\author Juan G Alonso Guzman
\date 02/22/2023
*/
void BackgroundSolarWindTermShock::EvaluateBackgroundDerivatives(void)
{
#if SOLARWIND_DERIVATIVE_METHOD == 0

   if(BITS_RAISED(_spdata._mask, BACKGROUND_gradU)) {
//TODO: complete
   };
   if(BITS_RAISED(_spdata._mask, BACKGROUND_gradB)) {
//TODO: complete
   };
   if(BITS_RAISED(_spdata._mask, BACKGROUND_gradE)) {
      _spdata.gradEvec = -((_spdata.gradUvec ^ _spdata.Bvec) + (_spdata.Uvec ^ _spdata.gradBvec)) / c_code;
   };
   if(BITS_RAISED(_spdata._mask, BACKGROUND_dUdt)) _spdata.dUvecdt = gv_zeros;
   if(BITS_RAISED(_spdata._mask, BACKGROUND_dBdt)) _spdata.dBvecdt = gv_zeros;
   if(BITS_RAISED(_spdata._mask, BACKGROUND_dEdt)) _spdata.dEvecdt = gv_zeros;

#else
   NumericalDerivatives();
#endif
};

/*!
\author Juan G Alonso Guzman
\date 02/23/2024
*/
void BackgroundSolarWindTermShock::EvaluateDmax(void)
{
   BackgroundSolarWind::EvaluateDmax();
   double r = posprime.Norm();
   if(r_TS - dmax0 < r && r < r_TS + w_TS + dmax0) {
      if(r < r_TS) _spdata.dmax += (dmax_TS - dmax0) * (r - r_TS + dmax0 / dmax0);
      else if(r > r_TS + w_TS) _spdata.dmax -= (dmax_TS - dmax0) * (r - r_TS - w_TS - dmax0 / dmax0);
      else _spdata.dmax = dmax_TS;
   };
};

};
