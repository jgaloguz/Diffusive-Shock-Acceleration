/*!
\file distribution_base.cc
\brief Implements a base class to record distributions
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "distribution_base.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistributionBase methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/14/2021
*/
DistributionBase::DistributionBase(void)
                : Params("", 0, STATE_NONE)
{
};

/*!
\author Vladimir Florinski
\date 06/14/2021
\param[in] name_in   Readable name of the class
\param[in] specie_in Particle's specie
\param[in] status_in Initial status
*/
DistributionBase::DistributionBase(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
                : Params(name_in, specie_in, status_in)
{
};

/*!
\author Vladimir Florinski
\date 06/14/2021
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDistribution()" with the argument of "true".
*/
DistributionBase::DistributionBase(const DistributionBase& other)
                : Params(other)
{
// Params' constructor resets all flags
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistribution(true);
};

/*!
\author Vladimir Florinski
\date 09/13/2022
\param [in] construct Whether called from a copy constructor or separately
*/
void DistributionBase::SetupDistribution(bool construct)
{
};

/*!
\author Vladimir Florinski
\date 06/18/2021
\param[in] cont_in Container with parameters

This is the default method to set up an object. It should only be defined in the base class (XXXXBase). Derived classes should _not_ modify it! This version always calls the correct virtual "SetupDistribution()" method.
*/
void DistributionBase::SetupObject(const DataContainer& cont_in)
{
   Params::SetContainer(cont_in);
   SetupDistribution(false);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 09/13/2022
*/
void DistributionBase::AddEvent(void)
{
};

/*!
\author Juan G Alonso Guzman
\date 09/13/2022
*/
void DistributionBase::AddRecord(void)
{
};

/*!
\author Juan G Alonso Guzman
\date 09/13/2022
\param[in] n_records_in Total number of events
*/
void DistributionBase::SetNRecords(int n_records_in)
{
};

/*!
\author Vladimir Florinski
\date 09/13/2022
*/
void DistributionBase::ResetDistribution(void)
{
};

/*!
\author Juan G Alonso Guzman
\date 09/13/2022
*/
void DistributionBase::ResetRecords(void)
{
};

/*!
\author Vladimir Florinski
\date 09/13/2022
\param[in] other Second distribution
\return Reference to this object
*/
DistributionBase& DistributionBase::operator +=(const DistributionBase& other)
{
   return *this;
};

/*!
\author Juan G Alonso Guzman
\date 09/13/2022
\param[in] other Second distribution
*/
void DistributionBase::CopyRecords(const DistributionBase& other)
{
};

/*!
\author Vladimir Florinski
\date 09/13/2022
\param[out] size Size of an element of the "distro" array
\return Pointer to the distribution storage
*/
void* DistributionBase::GetDistroAddress(size_t& size)
{
   size = 0;
   return nullptr;
};

/*!
\author Juan G Alonso Guzman
\date 09/13/2022
\param[out] size Size of an element of the "weights_record" array
\return Pointer to the weights storage
*/
void* DistributionBase::GetWeightsRecordAddress(size_t& size)
{
   size = 0;
   return nullptr;
};

/*!
\author Vladimir Florinski
\date 05/05/2022
\param[in] t1        First time value
\param[in] pos1      First position vector
\param[in] mom1      First momentum vector
\param[in] spdata1   First spatial data
\param[in] t2        Second time value
\param[in] pos2      Second position vector
\param[in] mom2      Second momentum vector
\param[in] spdata2   Second spatial data
\param[in] action_in Action to calculate the weight
*/
void DistributionBase::ProcessTrajectory(double t1, const GeoVector& pos1, const GeoVector& mom1, const SpatialData& spdata1,
                                         double t2, const GeoVector& pos2, const GeoVector& mom2, const SpatialData& spdata2,
                                         double amp2, double logwgt2, int action_in)
{
   SetState(t1, pos1, mom1);
   _spdata._mask = spdata1._mask;
   _spdata = spdata1;
   _t2 = t2;
   _pos2 = pos2;
   _mom2 = mom2;
   _spdata2._mask = spdata2._mask;
   _spdata2 = spdata2;
   _spdata2.Bmag_min = spdata2.Bmag_min;
   _spdata2.Bmag_max = spdata2.Bmag_max;
   _amp = amp2;
   _logwgt = logwgt2;
   EvaluateValue();
   EvaluateWeight(action_in);
   AddEvent();
   if (keep_records) AddRecord();
};

/*!
\author Vladimir Florinski
\date 09/13/2022
*/
void DistributionBase::EvaluateValue(void)
{
};

/*!
\author Vladimir Florinski
\date 09/13/2022
\param[in] action_in Action index from "ActionTable"
*/
void DistributionBase::EvaluateWeight(int action_in)
{
};

/*!
\author Vladimir Florinski
\date 09/13/2022
\param[in] file_name Distribution file name
*/
void DistributionBase::Dump(const std::string& file_name) const
{
};

/*!
\author Vladimir Florinski
\date 09/13/2022
\param[in] file_name Distribution file name
*/
void DistributionBase::Restore(const std::string& file_name)
{
};
   
/*!
\author Vladimir Florinski
\date 09/13/2019
\param[in] ijk        Which component to print (the remaining dimensions are collapsed)
\param[in] dist_name  Distribution file name
\param[in] phys_units Use physical units for output
*/
void DistributionBase::Print1D(int ijk, const std::string& dist_name, bool phys_units) const
{
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 06/29/2022
\param[in] ijk1       Component 1 to print
\param[in] ijk2       Component 2 to print
\param[in] dist_name  Distribution file name
\param[in] phys_units Use physical units for output
*/
void DistributionBase::Print2D(int ijk1, int ijk2, const std::string& dist_name, bool phys_units) const
{
};

/*!
\author Juan G Alonso Guzman
\date 09/13/2022
\param[in] dist_name  Distribution file name
\param[in] phys_units Use physical units for output
*/
void DistributionBase::PrintRecords(const std::string& dist_name, bool phys_units) const
{
};

};
