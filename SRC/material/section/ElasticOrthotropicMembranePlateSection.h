/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.0 $
// $Date: 2024-01-29$
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/ElasticOrthotropicMembranePlateSection.h,v $
//  Elastic Orthotropic Plate Section with membrane
//

// Implemented by Jose Gallardo 

#ifndef ElasticOrthotropicMembranePlateSection_h
#define ElasticOrthotropicMembranePlateSection_h

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <Vector.h>
#include <Matrix.h>
#include <ID.h>

#include <SectionForceDeformation.h>


class ElasticOrthotropicMembranePlateSection : public SectionForceDeformation{

//-------------------Declarations-------------------------------

  public : 

    //null constructor
    ElasticOrthotropicMembranePlateSection( ) ;

    //full constructor
    ElasticOrthotropicMembranePlateSection(int tag,
        double E1,
        double E2,
        double nu12,
        double nu21,
        double G12,
        double G13,
        double G23,
        double h = 1.0,
        double rho = 0.0);


    //destructor
    ~ElasticOrthotropicMembranePlateSection( ) ;

    //make a clone of this material
    SectionForceDeformation *getCopy( ) ;

    const char *getClassType(void) const {return "ElasticOrthotropicMembranePlate";};

    //send back order of strain in vector form
    int getOrder( ) const ;

    //send back order of strain in vector form
    const ID& getType( ) ;

    //swap history variables
    int commitState( ) ; 

    //revert to last saved state
    int revertToLastCommit( ) ;

    //revert to start
    int revertToStart( ) ;

    //get the strain and integrate plasticity equations
    int setTrialSectionDeformation( const Vector &strain_from_element ) ;

    //send back the strain
    const Vector& getSectionDeformation( ) ;

    //send back the stress 
    const Vector& getStressResultant( ) ;

    //send back the tangent 
    const Matrix& getSectionTangent( ) ;

    //send back the initial tangent 
    const Matrix& getInitialTangent( ) ;

    //print out data
    void Print( OPS_Stream &s, int flag ) ;

    //density per unit area
    double getRho() ;

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    // parameters
    int setParameter(const char** argv, int argc, Parameter& param);
    int updateParameter(int parameterID, Information& info);

  private :

    double E1;      // elastic modulus in direction 1
    double E2;      // elastic modulus in direction 2
    double nu12;    // poisson ratio in dir 12
    double nu21;    // poisson ratio in dir 21
    double G12;     // Shear modulus in dir 12
    double G13;     // Shear modulus in dir 13
    double G23;     // Shear modulus in dir 23
    double h;       // MembranePlate thickness
    double rho;     // density
    double rhoH;    //mass per unit 2D area

    static const double Ks ; // =5/6 = shear correction factor

    Vector strain ;
    static Vector stress ;
    static Matrix tangent ;
    static ID array ;  

} ; //end of ElasticOrthotropicMembranePlateSection declarations


#endif
