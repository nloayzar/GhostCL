#ifndef GHOST_C2G2_H //Usual macro guard to prevent multiple inclusion
#define GHOST_C2G2_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Daniel G. Figueroa, Adrien Florio, Francisco Torrenti,  Year: 2020

#include "CosmoInterface/cosmointerface.h"

//Include cosmointerface to have access to all of the library.

namespace TempLat
{
    /////////
    // Model name and number of fields
    /////////

    // In the following class, we define the defining parameters of your model:
    // number of fields of each species and the type of tinteractions.

    struct ModelPars : public TempLat::DefaultModelPars {
    	  static constexpr size_t NScalars = 1;
        static constexpr size_t NGhostScalars = 1;
        // In our phi4 example, we only want 2 scalar fields.
        static constexpr size_t NPotTerms = 1;
        static constexpr size_t NGhostMassTerms = 1;
        static constexpr size_t NGhostPotTerms = 1;
        // Our potential naturaly splits into two terms: the inflaton potential
        // and the interaction with the daughter field.

        // All the numbers of fields are 0 by default, so we need only
        // to specify that we want two scalar fields.
        // See the model with gauge fields to have an example of how to turn
        // them on and specify interactions.
    };

  #define MODELNAME ghost_c2g2
  // Here we define the name of the model. This should match the name of your file.

  template<class R>
  using Model = MakeModel(R, ModelPars);
  // In this line, we define an appropriate generic model, with the correct
  // number of fields, ready to be customized.
  // If you are curious about what this is doing, the macro is defined in
  // the "CosmoInterface/abstractmodel.h" file.

  class MODELNAME : public Model<MODELNAME>
  // Declaration of our model. It inherits from the generic model defined above.
  {
 //...
private:

  double q, mass_phi, mass_g, planckmass;
// Here are the declaration of the model specific parameters. They are 'private'
// to force you using them only within your model and not outside.

// Some parameters which are declared in the class "Model" and which are useful (they are all 'public'):

// fldS0, piS0 : arrays which should contain the initial homogeneous values of
//               the scalar fields
//
// alpha, fStar, omegaStar : time and field rescaling to go to program units.
//
// fldS : The actual object which contains the scalar fields.

  public:

    MODELNAME(ParameterParser& parser, RunParameters<double>& runPar, std::shared_ptr<MemoryToolBox> toolBox): //Constructor of our model.
    Model<MODELNAME>(parser,runPar.getLatParams(), toolBox, runPar.dt, STRINGIFY(MODELLABEL)) //MODELLABEL is defined in the cmake.
    {

      /////////
      // Independent parameters of the model (read from parameters file)
      /////////

      mass_phi = parser.get<double>("mass_phi");
      mass_g = parser.get<double>("mass_g");
      planckmass = parser.get<double>("planckmass",{2.435e18});
      //  We start by initializing our model paramteters. We read them from the
      // input file/command line.  Effectively, by calling 'par.get<double>("lambda")'
      // we declare a new parameter which needs to be in the input data.  Its name is
      // "lambda" and we specify it is a 'double'.

      q = parser.get<double>("q");
      // In the same way, we declare an input parameter 'q'.
      grav_ON = parser.get<bool>("grav_ON",{false});

      //g = sqrt(q*lambda);
      //For convenience, we also define g as a function of lambda and q.


        /////////
        // Initial homogeneous components of the fields
        // (read from parameters file, or specified here if not)
        /////////

        fldS0 = parser.get<double, 1>("initial_amplitudes");
        piS0 = parser.get<double, 1>("initial_momenta", {0});
        
        // Then, we need to specify the initial homogeneous
        // value of our fields. We read them again from the input file. The int '2' means
        // that we actually expect two values and that we will get an array of
        // double of size two.
        // Contrary to the "initial_amplitudes" parameter and the others above,
        //, the "initial_momenta" is an optional parameter. It can still be specified through
        //  command line or input file as initial_momenta=value1 value2 ... valueNs,
        // but it can also be omitted, as we specified a default value of '{0, 0}'.


        /////////
        // Rescaling for program variables
        /////////

        alpha = 0;
        fStar = mass_phi;
        omegaStar = mass_phi;
        PlanckMass = planckmass;
        // We now need to specify the rescaling from physical units to program units.
        // This consists of the  time rescaling exponent alpha, the field rescaling fStar
        // and the velocity rescaling omegaStar.
        // See the paper for more information on how to fix them.

        setInitialPotentialAndMassesFromPotential();
        // Here we call this function to compute the value of the potential on the homogeneous
        // initial condition  (useful to set the initial Hubble rate). We also compute
        // in this function the masses from the second derivative of the potential
        // evaluated on the homogeneous initial conditions. If you want to do something else,
        // uncomment the next section and do whatever suits your needs.

        /*
          masses2S = {..., ...};
          setInitialPotentialFromPotential();
         */
    }

   /////////
   // Program potential (add as many functions as terms are in the potential)
   /////////

    auto potentialTerms(Tag<0>) // Scalar mass term 
    {
        return 0.5 * pow<2>(fldS(0_c));
    }

    auto ghostpotentialTerms(Tag<0>) // Scalar Ghost interaction term 
    {
        return 0.5 * q * pow<2>(fldS(0_c)) * pow<2>(fldGS(0_c));
    }

    auto ghostMassTerms(Tag<0>) // Ghost mass term
    {
        return 0.5 * pow<2>(mass_g/omegaStar) * pow<2>(fldGS(0_c));
    }

   /////////
   // Derivatives of the program potential with respect fields
   // (add one function for each field).
   /////////

    auto potDeriv(Tag<0>)  // Derivative with respect of scalar of scalar mass term and interaction term
    {
      return   fldS(0_c) + q * fldS(0_c) * pow<2>(fldGS(0_c)) ;
    }

    auto potDerivGS(Tag<0>)  // Derivative with respect of ghost of interaction term
    {
      return  q * fldGS(0_c) * pow<2>(fldS(0_c));
    }

    auto potMassDerivGS(Tag<0>)  // Derivative respect of ghost of ghost mass term
    {
      return  pow<2>(mass_g/omegaStar) * (fldGS(0_c));
    }
	
   /////////
   //  Second derivatives of the program potential with respect fields
   // (add one function for each field)
   /////////

    auto potDeriv2(Tag<0>) // Second derivative with respect inflaton
    // Finally, for the purpose of initializing the masses, the user needs to define
    // in the same fashion the second derivatives of the potential
    // (put 'return 0' if you are not using this feature).
    {
      return  1.0 +  q * pow<2>(fldGS(0_c)) ;
    }

    auto potDerivG2(Tag<0>) // Second derivative with respect daughter field
    {
      return  pow<2>(mass_g/omegaStar) + q * pow<2>(fldS(0_c)) ;
    }
		
	
    };
}

#endif 
