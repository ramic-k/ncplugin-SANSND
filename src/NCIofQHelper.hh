#ifndef NCPlugin_IofQHelper_hh
#define NCPlugin_IofQHelper_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2025 NCrystal developers                                   //
//                                                                            //
//  Licensed under the Apache License, Version 2.0 (the "License");           //
//  you may not use this file except in compliance with the License.          //
//  You may obtain a copy of the License at                                   //
//                                                                            //
//      http://www.apache.org/licenses/LICENSE-2.0                            //
//                                                                            //
//  Unless required by applicable law or agreed to in writing, software       //
//  distributed under the License is distributed on an "AS IS" BASIS,         //
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  //
//  See the License for the specific language governing permissions and       //
//  limitations under the License.                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "NCrystal/NCPluginBoilerplate.hh"
#include "NCrystal/core/NCTypes.hh"
#include "NCPointwiseDist.hh"

namespace NC = NCrystal;

namespace NCPluginNamespace {

  //We implement the actual physics model in this completely custom C++ helper
  //class. That decouples it from NCrystal interfaces (which is nice in case the
  //NCrystal API changes at some point), and it makes it easy to directly
  //instantiate and test the modelling implementation from standalone C++ code.
  //
  //We mark the class as MoveOnly, to make sure it doesn't get copied around by
  //accident (since it could easily end up having large data members).
  class IofQHelper {
  public:

    //The class is similar to NC::IofQHelper but we redifine it here to implement
    //theta min feature
    //The class is constructed from Q and I(Q) values:

    IofQHelper( const NC::VectD& Q, const NC::VectD& IofQ, double thetaMin = 0);
    //Calculate the integral of Q*I(Q) from Q=0 to Qmax=2k, where k is the
    //wavenumber of the neutron of the provided energy. Note that to convert it
    //to a cross section one must still multiply it with a factor of c/E where c
    //is an appropriate constant and E is the neutron energy:
    double calcQIofQIntegral( NC::NeutronEnergy ) const;
    //Calculate the integral of Q*I(Q) from Q=2k*sin(theta_min) to Qmax=2k
    double calcQIofQIntegralMin( NC::NeutronEnergy ) const;

    //Sample a Q value according to Q*I(Q) over the interval from Q=0 to
    //Qmax=2k, where k is the wavenumber of the neutron of the provided energy:
    double sampleQValue( NC::RNG&, NC::NeutronEnergy ) const;
    double sampleQValueTrunc( NC::RNG&, NC::NeutronEnergy ) const;

    //QMax value:
    double getQMax() const;

  private:
    PointwiseDist m_pwdist;
    NC::NeutronEnergy m_ekinMax;
    double m_normFact;
    double m_thetaMin;
    struct internal_t;
    IofQHelper( internal_t );
  };
}

////////////////////////////
// Inline implementations //
////////////////////////////


inline double NCPluginNamespace::IofQHelper::calcQIofQIntegral( NC::NeutronEnergy ekin ) const
{
  if ( ekin >= m_ekinMax )
    return m_normFact;
  constexpr double kkk = 4.0 * NC::ekin2ksq(1.0);
  const double twok = std::sqrt( kkk * ekin.dbl() );
  return m_pwdist.commulIntegral( twok ) * m_normFact;
}

inline double NCPluginNamespace::IofQHelper::calcQIofQIntegralMin( NC::NeutronEnergy ekin ) const
{
  if ( ekin >= m_ekinMax )
    return m_normFact;
  
  constexpr double kkk = 4.0 * NC::ekin2ksq(1.0);
  const double twok = std::sqrt( kkk * ekin.dbl() );
  double fullInt = m_pwdist.commulIntegral( twok ) * m_normFact;
  
  double lowerInt = m_pwdist.commulIntegral( twok*std::sin(m_thetaMin) ) * m_normFact;
  return fullInt - lowerInt;
}

inline double NCPluginNamespace::IofQHelper::sampleQValue( NC::RNG& rng, NC::NeutronEnergy ekin ) const
{
  constexpr double kkk = 4.0 * NC::ekin2ksq(1.0);
  const double twok = std::sqrt( kkk * std::min<double>(m_ekinMax.dbl(),ekin.dbl()) );
  return m_pwdist.sampleBelow( rng, twok );
}
inline double NCPluginNamespace::IofQHelper::sampleQValueTrunc( NC::RNG& rng, NC::NeutronEnergy ekin ) const
{
  constexpr double kkk = 4.0 * NC::ekin2ksq(1.0);
  const double twok = std::sqrt( kkk * std::min<double>(m_ekinMax.dbl(),ekin.dbl()) );
  return m_pwdist.sampleBelowTrunc( rng, twok, twok*std::sin(m_thetaMin) );
}

#endif
