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
                                                                        

// $Source: /usr/local/cvs/OpenSees/SRC/element/TSM_2D.h,v $
                                                                        
#ifndef TSM_2D_h
#define TSM_2D_h



#include <Element.h>
#include <Matrix.h>
#include <Vector.h>

class Channel;
class UniaxialMaterial;
class Response;

//#define ELE_TAG_TSM_2D 100002

 class TSM_2D : public Element
{
  public:

      // Constructors
      TSM_2D(int tag, int Nd1, int Nd2,
          double Kvi, double fas, double Kb, double Kr, double Do, double Di,
          double Tr, double fc, double a, double k, double pm,
          UniaxialMaterial** theMaterials, const Vector y = 0,
          double mass = 0.0);

      TSM_2D();
      ~TSM_2D();

      // method to get class type
	  const char *getClassType(void) const {return "TSM_2D";};

      // public methods to obtain information about dof & connectivity    
      int getNumExternalNodes(void) const;
      const ID &getExternalNodes(void);
      Node **getNodePtrs(void);
      int getNumDOF(void);	
      void setDomain(Domain *theDomain);

      // public methods to set the state of the element    
      int commitState(void);
      int revertToLastCommit(void);        
      int revertToStart(void);        
      int update(void);
    
      // public methods to obtain stiffness, mass, damping and 
      // residual information    
      const Matrix &getTangentStiff();
      const Matrix &getInitialStiff();  
      const Matrix &getDamp();
      const Matrix &getMass();

      void zeroLoad();	
      int addLoad(const Vector &addP);
      int addInertiaLoadToUnbalance(const Vector &accel);
      const Vector &getResistingForce(void);
      const Vector &getResistingForceIncInertia(void);            

      // public methods for element output
      int sendSelf(int commitTag, Channel &theChannel);
      int recvSelf(int commitTag, Channel &theChannel,
          FEM_ObjectBroker &theBroker);
      int displaySelf(Renderer &, int mode, float fact, const char **displayModes=0, int numModes=0);
      void Print(OPS_Stream &s, int flag =0);    

      Response *setResponse(const char **argv, int argc, OPS_Stream &);
      int getResponse(int responseID, Information &eleInformation);

      int setParameter(const char **argv, int argc, Parameter &param);
	  int updateParameter (int parameterID, Information &info);
	  
  protected:
  
  private:
      void setUp();

      ID  connectedExternalNodes;   // contains the tags of the end nodes
      Node* theNodes[2];                  // array of nodes
      UniaxialMaterial* theMaterials[1];  // array of uniaxial materials

      // parameters

      double abs(double);
      double max(double, double);
      double expo(double);

      double Kvi;             // initial vertical stiffness
      double fas;             // factor of vertical stiffness softening due lateral deformation
      double Kb;              // bending stiffness
      double Kr;              // Stiffness of the rotational spring
      double Tr;              // rubber height
      double fc;              // cavitation force
      double a;               // rate of cavitation damage
      double k;               // cavitation parameter
      double pm;              // maximum damage

      double mass;            // mass of the element
      Vector y, x;            // local y, and x direction
      double L;               // element length

      double Do, Di;           // outer and inner diameter


      // internal variables
      double s;
      double theta;
      double um;
      double fpc;

      //trial state

      double t_s;
      double t_theta;
      double t_um;
      double t_fpc;

      // state variables
    
      Vector ub;          // displacements in basic system
      Vector qb;          // forces in basic system
      Matrix kt;          // stiffness matrix in basic system
      Vector ul;          // displacements in local system
      Matrix Tgl;         // transformation matrix from global to local system
      Matrix Tlb;         // transformation matrix from local to basic system

      // committed history variables
      Vector ub_c;         // displacement in basic system

      // initial stiffness matrix in basic system
      Matrix kinit;

      static Matrix theMatrix;             // matrix to return stiff, damp & mass
      Vector theVector;             // vector to return the residual
      Vector theLoad;
};

#endif

