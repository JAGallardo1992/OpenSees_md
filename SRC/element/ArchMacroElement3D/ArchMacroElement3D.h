#ifndef ArchMacroElement3D_h // Proveri da li treba da bude istog naziva kao i file? A: Ne treba, tu definises macro. 
#define ArchMacroElement3D_h

#define PY_SSIZE_T_CLEAN
#include <Python.h>


#ifndef _bool_h
#include "bool.h"
#endif

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
//#include "C:/Users/ivana/BackUp/VaultMacroElement_ENAC-IT/ESSD-OpenSees/SRC/element/forceBeamColumn/BeamIntegration.h" TODO: Promenila si 
// da imas standardne vrednosti intLength, i sa dva sectiona fiksno
#include <vector>
#include <string>

class Node;
class SectionForceDeformation;
class CrdTransf; // Da li ovo treba CrdTransf, ili moja matrica za geom transformaciju? A: Posto imam normalne dof u cvorovima, moze ova default CrdTransf
class Response;
class Channel;

class ArchMacroElement3D : public Element
{
public:
    ArchMacroElement3D(int tag,
        int nd1, int nd2,
        SectionForceDeformation* sI, SectionForceDeformation* sJ,
        Vector intLength, Vector intLengthMasses,
        CrdTransf& coordTransf,
        double dt, double theta, const char* type);

    ArchMacroElement3D();
    ~ArchMacroElement3D(); //ovo  ~ je destructor

    const char* getClassType(void) const { return "ArchMacroElement3D"; }; // Is it needed? A:  Seems that is for debugging

    // public methods to obtain information about dof & connectivity    
    int getNumExternalNodes(void) const;
    const ID& getExternalNodes(void);
    Node** getNodePtrs(void);

    int getNumDOF(void);
    void setDomain(Domain* theDomain); // initialization

    //  ---- methods dealing with committed state and update - u ME nisu jednake nuli! DA LI TREBA = 0??
    int commitState(void);  // called when a converged solution has been obtained for a time step
    int revertToLastCommit(void);

    //-----
    int revertToStart(void); //FIXME: Ova je dodata, da bude isto kao FV, iako mi nije potrebna
    //-----
    // public methods to obtain stiffness, mass, damping and residual information
    int update(void);  // called when a new trial step has been set at the nodes



    const Matrix& getTangentStiff(void);

    //const Matrix &getBasicStiff(Matrix &kb, int initial);
    const Matrix& getInitialStiff(void);
    const Matrix& getMass(void);
    const Matrix& getInitialBasicStiff(void);


    void zeroLoad();

    //-----
    int addLoad(ElementalLoad* theLoad, double loadFactor);  //FIXME: Ova je dodata, da bude isto kao FV, iako mi nije potrebna
    //-----


    int addInertiaLoadToUnbalance(const Vector& accel);

    const Vector& getResistingForce(void);
    const Vector& getResistingForceIncInertia(void);

    // public methods for element output

    int sendSelf(int commitTag, Channel& theChannel);
    int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);
    void Print(OPS_Stream& s, int flag = 0);

    Response* setResponse(const char** argv, int argc, OPS_Stream& s);
    int getResponse(int responseID, Information& eleInfo);

    // damping methods - FIXME: ovde si koment Damping 
    int setRayleighDampingFactors(double alphaM, double betaK, double betaK0, double betaKc);
    const Matrix& getDamp(void);

    //-----
    //const Matrix &getSecantStiff(void); //FIXME: Ova je dodata, da bude isto kao FV, iako mi nije potrebna
    //-----




protected:
    const Vector& getRayleighDampingForces(void); // set as zero. Here or in cpp? Probably cpp


    // In your .h file
private:

    //const Matrix &getInitialBasicStiff(void);

    // explicit matrix operations to transform from local to global accounting for node offset(s)
      //int trasformMatrixToGlobal(Matrix& A);  

    //-----
    const int numSections;                      // number of sections = 2
    SectionForceDeformation** theSections;      // pointers to sectional models
    //BeamIntegration* beamInt; // Promenila si da imas intLength
    CrdTransf* crdTransf;           // pointer to coordinate tranformation object
    const char* type; // Add this line to hold the type name


    ID  connectedExternalNodes;                 // tags of element nodes (i, j)
    Node* theNodes[2];                          // pointers to element nodes (i,j) 

    static Matrix K;                            // stores element stiffness, damping, or mass Matrix
    static Vector P;                            // Element resisting force vector - Treba za getResistingForce - u getResponse

    Vector Q;                                   // Applied nodal loads
    Vector q;                                   //  Basic force
    double q0[6];                                // Fixed end forces in basic system
    double p0[6];                                // Reactions in basic system of distributed loads 



    // DO I need this? I think yes, for connectedExternalNode
    int parameterID;   //ovo je isto nesto za sensitivity analysis
    static double workArea[];  // not used here

    Vector intLength, intLengthMasses;           // integration lengths for sectional responses and lumped masses
    double dt;
    double theta;




    /* int dampingModel;*/


     //--------This from VaultMacroELement: --------------
     // static const std::string acc_file_name;
     // static const std::string result_file_name;
     // static const std::string python_script_name;

     // std::vector<std::vector<double>> result;

     // Load output in ::result
     //void loadResult();
   //---------------------------------------------------    

   // // Link with Python
    int initializePython();
    PyObject* pModule;
    PyObject* pUpdateFunc;
    PyObject* pGetResistingForceFunc;
    PyObject* pGetInitialStiffFunc;
    PyObject* pCommitFunc;
    PyObject* pIsStateCommitted;
    PyObject* pGetMassMomInertiaFunc;
    PyObject* pSetDomainFunc;
    PyObject* pPassDtThetaFunc;
    

};



#endif