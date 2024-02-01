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
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/ElasticOrthotropicMembranePlateSection.cpp,v $
//  Elastic Orthotropic Plate Section with membrane
//

// Implemented by Jose Gallardo


#include <ElasticOrthotropicMembranePlateSection.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <Parameter.h>

//parameters
const double ElasticOrthotropicMembranePlateSection::Ks = 5.0/6.0 ; //shear correction

//static vector and matrices
Vector  ElasticOrthotropicMembranePlateSection::stress(8) ;
Matrix  ElasticOrthotropicMembranePlateSection::tangent(8,8) ;
ID      ElasticOrthotropicMembranePlateSection::array(8) ;

void* OPS_ElasticOrthotropicMembranePlateSection()
{
    if (OPS_GetNumRemainingInputArgs() < 8) {
        opserr << "WARNING insufficient arguments\n";
        opserr << "Want: section ElasticOrthotropicMembranePlateSection tag? E1? E2? nu12? nu21? G12? G13? G23? h? <rho?>\n";
        return 0;
    }

    int tag;
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
        opserr << "WARNING invalid tag\n";
        return 0;
    }

    double data[9]  = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    numdata = OPS_GetNumRemainingInputArgs();

    if (numdata > 8) numdata = 8;

    if (numdata < 7) {
        opserr << "WARNING invalid double values\n";
        return 0;
    }

    if (OPS_GetDoubleInput(&numdata, data) < 0) {
        opserr << "WARNING invalid double values\n";
        return 0;
    }

    return new ElasticOrthotropicMembranePlateSection(tag, data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7]);
}

//null constructor
ElasticOrthotropicMembranePlateSection::ElasticOrthotropicMembranePlateSection( ) :
SectionForceDeformation( 0, SEC_TAG_ElasticOrthotropicMembranePlateSection )
{ 

}

//full constructor
ElasticOrthotropicMembranePlateSection::ElasticOrthotropicMembranePlateSection(int tag,
    double _E1, double _E2, double _nu12, double _G12, double _G13, double _G23,
    double _h, double _rho)
    : SectionForceDeformation(tag, SEC_TAG_ElasticOrthotropicMembranePlateSection),
    E1(_E1), E2(_E2), nu12(_nu12), G12(_G12), G13(_G13), G23(_G23), h(_h), rho(_rho), strain(8)
{
    rhoH  = rho * h;
	nu21 = nu12*E2/E1;

}



//destructor
ElasticOrthotropicMembranePlateSection::~ElasticOrthotropicMembranePlateSection( )
{ 

} 



//make a clone of this material
SectionForceDeformation*  ElasticOrthotropicMembranePlateSection::getCopy( )
{
  ElasticOrthotropicMembranePlateSection *clone ;

  clone = new ElasticOrthotropicMembranePlateSection(this->getTag(), E1, E2, nu12, G12, G13, G23, h, rho) ; //new instance of this class

  //    *clone = *this ; //assignment to make copy
  clone->rhoH = this->rhoH ;
  clone->strain = this->strain;

  return clone ;
}

//density per unit area
double ElasticOrthotropicMembranePlateSection::getRho( )
{
  return rhoH ;
}


//send back order of strain in vector form
int ElasticOrthotropicMembranePlateSection::getOrder( ) const
{
  return 8 ;
}


//send back order of strain in vector form
const ID& ElasticOrthotropicMembranePlateSection::getType( )
{
    static bool initialized = false;
    if (!initialized) {
        array(0) = SECTION_RESPONSE_FXX;
        array(1) = SECTION_RESPONSE_FYY;
        array(2) = SECTION_RESPONSE_FXY;
        array(3) = SECTION_RESPONSE_MXX;
        array(4) = SECTION_RESPONSE_MYY;
        array(5) = SECTION_RESPONSE_MXY;
        array(6) = SECTION_RESPONSE_VXZ;
        array(7) = SECTION_RESPONSE_VYZ;
        initialized = true;
    }
    return array;
}



//swap history variables
int ElasticOrthotropicMembranePlateSection::commitState( )
{
  return 0 ;
}

//revert to last saved state
int ElasticOrthotropicMembranePlateSection::revertToLastCommit( )
{
  return 0 ;
}

//revert to start
int ElasticOrthotropicMembranePlateSection::revertToStart( )
{
  return 0 ;
}

//get the strain 
int ElasticOrthotropicMembranePlateSection::setTrialSectionDeformation( const Vector &strain_from_element)
{
  this->strain = strain_from_element ;

  return 0 ;
}


//send back the strain
const Vector& ElasticOrthotropicMembranePlateSection::getSectionDeformation( )
{
  return this->strain ;
}


//send back the stress 
const Vector&  ElasticOrthotropicMembranePlateSection::getStressResultant( )
{

    //membrane behavior
    double A11 = h * (E1 / (1.0 - nu12 * nu21));
    double A12 = h * (E1 * nu12 / (1 - nu12 * nu21));
    double A22 = h * (E2 / (1.0 - nu12 * nu21));
    double A66 = h * G12;

    stress(0) = A11 * strain(0) + A12 * strain(1);
    stress(1) = A12 * strain(0) + A22 * strain(1);
    stress(2) = A66 * strain(2);

    // bending 
    double D11 = (h * h * h / 12.0) * (E1 / (1.0 - nu12 * nu21));
    double D12 = (h * h * h / 12.0) * (nu12 * E2 / (1.0 - nu12 * nu21));
    double D22 = (h * h * h / 12.0) * (E2 / (1.0 - nu12 * nu21));
    double D66 = (h * h * h / 12.0) * G12;

    stress(3) = -(D11 * strain(3) + D12 * strain(4)); 
    stress(4) = -(D12 * strain(3) + D22 * strain(4));
    stress(5) = -D66 * strain(5);

    // Shear
    stress(6) = G13 * strain(6);
    stress(7) = G23 * strain(7);

 
  return this->stress ;
}


//send back the tangent 
const Matrix&  ElasticOrthotropicMembranePlateSection::getSectionTangent( )
{

    double A11 = h * (E1 / (1.0 - nu12 * nu21));
    double A12 = h * (E1 * nu12 / (1 - nu12 * nu21));
    double A22 = h * (E2 / (1.0 - nu12 * nu21));
    double A66 = h * G12;

    tangent.Zero();

  // Membrane tangent terms
  tangent(0, 0) = A11;
  tangent(1, 1) = A22;
  tangent(0, 1) = A12;
  tangent(1, 0) = tangent(0, 1);
  tangent(2, 2) = A66;

  // Bending tangent terms

  double D11 = (h * h * h / 12.0) * (E1 / (1.0 - nu12 * nu21));
  double D12 = (h * h * h / 12.0) * (nu12 * E2 / (1.0 - nu12 * nu21));
  double D22 = (h * h * h / 12.0) * (E2 / (1.0 - nu12 * nu21));
  double D66 = (h * h * h / 12.0) * G12;

  

  tangent(3, 3) = -D11;
  tangent(4, 4) = -D22;

  tangent(3, 4) = -D12;
  tangent(4, 3) = tangent(3, 4);

  tangent(5, 5) = -D66;

  // Shear tangent terms
  tangent(6, 6) = G13;
  tangent(7, 7) = G23;

  return this->tangent ;
}


//send back the initial tangent 
const Matrix&  ElasticOrthotropicMembranePlateSection::getInitialTangent( )
{

    double A11 = h * (E1 / (1.0 - nu12 * nu21));
    double A12 = h * (E1 * nu12 / (1 - nu12 * nu21));
    double A22 = h * (E2 / (1.0 - nu12 * nu21));
    double A66 = h * G12;

  tangent.Zero() ;

  // Membrane tangent terms

  tangent(0, 0) = A11;
  tangent(1, 1) = A22;

  tangent(0, 1) = A12;
  tangent(1, 0) = tangent(0, 1);

  tangent(2, 2) = A66;


  double D11 = (h * h * h / 12.0) * (E1 / (1.0 - nu12 * nu21));
  double D12 = (h * h * h / 12.0) * (nu12 * E2 / (1.0 - nu12 * nu21));
  double D22 = (h * h * h / 12.0) * (E2 / (1.0 - nu12 * nu21));
  double D66 = (h * h * h / 12.0) * G12;
  // Bending tangent terms

  tangent(3, 3) = -D11;
  tangent(4, 4) = -D22;
  tangent(3, 4) = -D12;
  tangent(4, 3) = tangent(3, 4);
  tangent(5, 5) = -D66;

  // Shear tangent terms
  tangent(6, 6) = G13;
  tangent(7, 7) = G23;

  return this->tangent ;
}


//print out data
void  ElasticOrthotropicMembranePlateSection::Print( OPS_Stream &s, int flag )
{
    if (flag == OPS_PRINT_PRINTMODEL_SECTION) {
        s << "ElasticOrthotropicMembranePlateSection: \n ";
        s << "  Modulus of elasticity: dir 1, E1 = " << E1 << endln;
        s << "  Modulus of elasticity: dir 2, E2 = " << E2 << endln;
        s << "  Poisson's Ratio: dir 12, nu12 = " << nu12 << endln;
        s << "  Poisson's Ratio: dir 21, nu21 = " << nu21 << endln;
        s << "  Shear modulus: dir 12, G12 = " << G12 << endln;
        s << "  Shear modulus: dir 13, G13 = " << G13 << endln;
        s << "  Shear modulus: dir 23, G23 = " << G23 << endln;
        s << "  Thickness, h = " << h << endln;
        s << "  Density, rho = " << rho << endln;
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"ElasticOrthotropicMembranePlateSection\", ";
        s << "\"E1\": " << E1 << ", ";
        s << "\"E1\": " << E2 << ", ";
        s << "\"nu12\": " << nu12 << ", ";
        s << "\"nu21\": " << nu21 << ", ";
        s << "\"G12\": " << G12 << ", ";
        s << "\"G13\": " << G13 << ", ";
        s << "\"G23\": " << G23 << ", ";
        s << "\"thickness\": " << h << ", ";
        s << "\"density\": " << rho << "}";
    }
}


int ElasticOrthotropicMembranePlateSection::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(11);
  data(0) = this->getTag();
  data(1) = E1;
  data(2) = E2;
  data(3) = nu12;
  data(4) = nu21;
  data(5) = G12;
  data(6) = G13;
  data(7) = G23;
  data(8) = h;
  data(9) = rho;
  data(10) = rhoH;

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "ElasticOrthotropicMembranePlateSection::sendSelf() - failed to send data\n";

  return res;
}


int ElasticOrthotropicMembranePlateSection::recvSelf(int cTag, Channel &theChannel,
				      FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(6);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "ElasticOrthotropicMembranePlateSection::recvSelf() - failed to recv data\n";
  else {
    this->setTag(data(0));
    E1 = data(1);
    E2 = data(2);
    nu12 = data(3);
    nu21 = data(4);
    G12 = data(5);
    G13 = data(6);
    G23 = data(7);
    h = data(8);
    rho = data(9);
    rhoH = data(10);
  }

  return res;
}

int ElasticOrthotropicMembranePlateSection::setParameter(const char** argv, int argc, Parameter& param)
{
    if (argc < 1)
        return -1;
    if (strcmp(argv[0], "E1") == 0) {
        param.setValue(E1);
        return param.addObject(1, this);
    }
    if (strcmp(argv[0], "E2") == 0) {
        param.setValue(E2);
        return param.addObject(2, this);
    }
    if (strcmp(argv[0], "nu12") == 0) {
        param.setValue(nu12);
        return param.addObject(3, this);
    }
    if (strcmp(argv[0], "nu21") == 0) {
        param.setValue(nu21);
        return param.addObject(4, this);
    }
    if (strcmp(argv[0], "G12") == 0) {
        param.setValue(G12);
        return param.addObject(5, this);
    }
    if (strcmp(argv[0], "G13") == 0) {
        param.setValue(G13);
        return param.addObject(6, this);
    }
    if (strcmp(argv[0], "G23") == 0) {
        param.setValue(G23);
        return param.addObject(7, this);
    }
    if (strcmp(argv[0], "h") == 0) {
        param.setValue(h);
        return param.addObject(8, this);
    }
    if (strcmp(argv[0], "rho") == 0) {
        param.setValue(rho);
        return param.addObject(9, this);
    }
    return -1;
}

int ElasticOrthotropicMembranePlateSection::updateParameter(int parameterID, Information& info)
{
    if (parameterID == 1) {
        E1 = info.theDouble;
    }
    else if (parameterID == 2) {
        E2 = info.theDouble;
    }
    else if (parameterID == 3) {
        nu12 = info.theDouble;
    }
    else if (parameterID == 4) {
        nu21 = info.theDouble;
    }
    else if (parameterID == 5) {
        G12 = info.theDouble;
    }
    else if (parameterID == 6) {
        G13 = info.theDouble;
    }
    else if (parameterID == 7) {
        G23 = info.theDouble;
    }
    else if (parameterID == 8) {
        h = info.theDouble;
    }
    else if (parameterID == 9) {
        rho = info.theDouble;
    }
    return 0;
}
