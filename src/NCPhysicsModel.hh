#ifndef NCPlugin_PhysicsModel_hh
#define NCPlugin_PhysicsModel_hh


#include "NCrystal/NCPluginBoilerplate.hh"//Common stuff (includes NCrystal
                                          //public API headers, sets up
                                          //namespaces and aliases)
#include "NCIofQHelper.hh"

namespace NCPluginNamespace {

  class PhysicsModel final : public NC::MoveOnly {
  public:
    //after debugging
    //enum class model_def : unsigned { FILE=0, PPF=1, GPF=2, HSFBA=3 };
    //A few static helper functions which can extract relevant data from NCInfo
    //objects (the createFromInfo function will raise BadInput exceptions in
    //case of syntax errors in the @CUSTOM_ section data):

    static bool isApplicable( const NC::Info& );
    static PhysicsModel createFromInfo( const NC::Info& );//will raise BadInput in case of syntax errors
    enum class Model : unsigned { FILE=0, PPF=1, GPF=2, HSFBA=3 };
    //Constructor gets the models string and vector of parameters:
    PhysicsModel( Model model, NC::VectD param, double thetaMin = 0.0 );
    //Constructor gets the models string and the dist R file or input I(q) file:
    PhysicsModel( Model model, std::string filename, double thetaMin = 0.0 );

    //Provide cross sections for a given neutron:
    double calcCrossSection( double neutron_ekin ) const;
    //Sample scattering vector from inverse CDF (rng is random number stream).
    double sampleScatteringVector( NC::RNG& rng, double neutron_ekin ) const;

    //Sample scattering event. Results are given
    //as the final ekin of the neutron and scat_mu which is cos(scattering_angle).
    struct ScatEvent { double ekin_final, mu; };
    ScatEvent sampleScatteringEvent( NC::RNG& rng, double neutron_ekin ) const;

  private:
    //Data members:
    Model m_model;
    NC::Optional<NC::VectD> m_param;
    NC::Optional<NCP::IofQHelper> m_helper;
    double m_thetaMin{0.0};
  };

}
#endif
