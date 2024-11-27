#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "ArchMacroElement3D.h"
#include <Channel.h>
#include <CrdTransf.h>
#include <CompositeResponse.h>
#include <Domain.h>
#include <ID.h>
#include <Information.h>
#include <Node.h>
#include <Parameter.h>
#include <Renderer.h>
#include <TransientIntegrator.h>
#include <Vector.h>
#include <elementAPI.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <FEM_ObjectBroker.h>
#include <Matrix.h>
#include <math.h>
#include <map>
#include <ElementIter.h>
#include <SectionForceDeformation.h>
#include <MovableObject.h>
#include <TaggedObject.h>
#include <ElasticSection3d.h>
//#include <ElasticSection2d.h>   TODO: Proveri da li ti ovo treba
/* -------------------------------------------------------------------------- */
#include <algorithm>
#include <limits>
#include <string.h>
#include <cmath>
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
#ifndef DBL_EPSILON
#define DBL_EPSILON (std::numeric_limits<double>::epsilon())
#endif
#include <typeinfo>
#include <string>

// -----This from ME: -------------------
#ifdef _USRDLL
#include <windows.h>
#define OPS_Export extern "C" _declspec(dllexport)
#elif _MACOSX
#define OPS_Export extern "C" __attribute__((visibility("default")))
#else
#define OPS_Export extern "C"
#endif
//------------------------------------------

Matrix ArchMacroElement3D::K(12, 12); //Iz ME je 18 x 18, u 2d bi trebalo ba bude 6x6,  u 3d 12x12
Vector ArchMacroElement3D::P(12);    // Iz ME je 18 - ti ovde imas 2 node = 12
double ArchMacroElement3D::workArea[200];
// initialize static variables

static int numMyMacroelement = 0;



OPS_Export void*
OPS_ArchMacroElement3D()
{
    opserr << "OPS f START";

    Vector intLength(2);
    intLength(0) = 1.0 / 2;
    intLength(1) = 1.0 / 2;

    Vector intLengthMasses(2);
    intLengthMasses(0) = 0.50;
    intLengthMasses(1) = 0.50;
    char filePath[1024] = { 0 };  // Ensure it's clear

    if (numMyMacroelement == 0) {
        opserr << "ArchMacroElement3D loaded from external library - Written by IB, EPFL, 2023" << endln;
        numMyMacroelement++;
    }

    if (OPS_GetNumRemainingInputArgs() < 7) {
        opserr << "insufficient arguments:eleTag,iNode,jNode,transfTag, dt, theta , type \n";
        return 0;
    }

    // inputs: 
    int iData[4];
    int numData = 4;
    if (OPS_GetIntInput(&numData, &iData[0]) < 0) {
        opserr << "WARNING: invalid integer inputs" << endln;
        return 0;
    }

    double dData[2]; //dt, theta
    numData = 2;
    if (OPS_GetDoubleInput(&numData, &dData[0]) < 0) {
        opserr << "WARNING: invalid double inputs" << endln;
        return 0;
    }

    double dt = dData[0];
    double theta = dData[1];
   
    const char* type = OPS_GetString();
    if (type == nullptr) {
        opserr << "Failed to retrieve string or no string to retrieve.\n";
        return nullptr;
    }

    // options
   /* int dampingModel = 0;*/
    ///////////////////////////////////////////////////////////////////
    // CODE FOR PASSING THETA AND DT, UNFORTUNATELY IT CANT WOEK WELL IF PUT HERE
    /*opserr << "OPS dt je: " << dt << endln;
    opserr << "OPS theta je: " << theta << endln;*/

    // Pas dt to Python code

    //initializePython();
    //Vector paramTHAnalysis(2); // vector for dt and theta

    //paramTHAnalysis(0) = dt;
    //paramTHAnalysis(1) = theta;
    ////

    ////// Convert Vectors to Python lists
    //PyObject* pparamTHAnalysisVec;
    //convertVectorToPyList(paramTHAnalysis, &pparamTHAnalysisVec);

    ////// Create a Python tuple from the lists
    //PyObject* pArgss = PyTuple_Pack(1, pparamTHAnalysisVec);

    ////if (!pArgss) {
    ////    PyErr_Print();
    ////    err = 1;
    ////}
    ////// Call python function
    //PyObject_CallObject(ppassParamTHAnalysis, pArgss);

    ////// Clean up
    //Py_DECREF(pArgss);
    //Py_DECREF(pparamTHAnalysisVec);

    //// Print the Vector contents using opserr for debugging or verification
    //opserr << "Vector paramTHAnalysis Contents: [";
    //opserr << paramTHAnalysis(0) << ", " << paramTHAnalysis(1);
    //opserr << "]" << endln;
    ///////////////////////////////////////////////////////////////////
    // 
    // Used for reading a single double value

  //while (OPS_GetNumRemainingInputArgs() > 0) {
  //    const char* type = OPS_GetString();
  //    if (strcmp(type, "-density") == 0) {

  //        numData = 1;
  //        // Read the next argument as the rho (density) value
  //        if (OPS_GetDoubleInput(&numData, &rho) < 0) {
  //            opserr << "WARNING: invalid density value" << endln;
  //            return 0;  // Exit if rho (density) value is not provided correctly
  //        }
  //    }
  //   
  //}
  //opserr << "Parsed rho value2: " << rho << endln; //DEBUG:Ovde vrati dobru vrednost rho. 

// check transf
    CrdTransf* theTransf = OPS_getCrdTransf(iData[3]);
    if (theTransf == 0) {
        opserr << "WARNING: invalid density value" << endln;
        opserr << "coord transfomration not found\n";
        return 0;
    }

    //FIXME: Da li treba nullpts u zagradi ili ne? ME3d ima to u ()
    // SectionForceDeformation* theSectionI{nullptr};
      // SectionForceDeformation* theSectionJ{nullptr};
    SectionForceDeformation* theSectionI{ nullptr };
    SectionForceDeformation* theSectionJ{ nullptr };

    // oov mi ne treba
    int E_ = 5;
    int A = 10;
    int Iz = 15;
    int Iy = 15;
    int G = 2;
    int J = 20;

    theSectionI = new ElasticSection3d(0, E_, A, Iz, Iy, G, J);
    theSectionJ = new ElasticSection3d(0, E_, A, Iz, Iy, G, J);


    Element* theEle = new ArchMacroElement3D(iData[0], iData[1], iData[2], theSectionI, theSectionJ, intLength, intLengthMasses, *theTransf, dData[0], dData[1], type);

    //FIXME: Da li treba i ovo za deleteSection ? 

       //  if (theSectionI!=NULL)  delete theSectionI;
       //  if (theSectionJ!=NULL)  delete theSectionJ;
   //ZA DEBUG:
// Assuming 'theEle' is a pointer to an existing VaultMacroElement2D object
// and 'crdTransf' is public or has a public getter method

    // Assuming you have a valid Domain object named 'theDomain'
        // Set the domain of the element

    // Now, you want to check the initial length
    // If 'crdTransf' has been initialized properly in 'setDomain', you can access the initial length
    //VaultMacroElement2D* vaultElement = dynamic_cast<VaultMacroElement2D*>(theEle);

    //if (vaultElement != nullptr) {
    //    // Call the Print function
    //    OPS_Stream& s = opserr;
    //    vaultElement->Print(s, 0);
    //}
    //else {
    //    // Handle the error if the cast fails
    //    opserr << "ERROR: Casting to VaultMacroElement2D failed." << endln;
    //}

    // Now, set the domain of the element
    //Domain* theDomain = OPS_GetDomain();
    //if (!theDomain) {
    //    opserr << "WARNING: theDomain is null" << endln;
    //    return 0;
    //}

    //// Set the domain of the element
    //vaultElement->setDomain(theDomain);


    const Matrix& K = static_cast<ArchMacroElement3D*>(theEle)->getMass();

    //// Print the mass matrix to the screen
    //opserr << "Mass Matrix:" << endln;
    //for (int i = 0; i < K.noRows(); ++i) {
    //    for (int j = 0; j < K.noCols(); ++j) {
    //        opserr << K(i, j) << " ";
    //    }
    //    opserr << endln; // Use 'endln' instead of 'std::endl' since you're using the OpenSees opserr stream
    //}
    // Check if the dynamic cast was successful
    // //DEBUG STIFFNESS MATRICES
        // Call the getInitialBasicStiff function
    //const Matrix& kb = vaultElement->getInitialBasicStiff();

    //opserr << "Matrix kb:" << endln;
    //for (int i = 0; i < kb.noRows(); ++i) {
    //    for (int j = 0; j < kb.noCols(); ++j) {
    //        opserr << kb(i, j) << " ";
    //    }
    //    opserr << endln;
    //}
    //const Matrix& K = vaultElement->getInitialStiff();

    //// Print the stiffness matrix to the screen
    //opserr << "Stiffness Matrix:" << endln;
    //for (int i = 0; i < K.noRows(); ++i) {
    //    for (int j = 0; j < K.noCols(); ++j) {
    //        opserr << K(i, j) << " ";
    //    }
    //    opserr << endln;
    //}
    //double initialLength = theTransf->getInitialLength();
    //opserr << "Initial Length: " << initialLength << endln;
    // Print the retrieved string to check it
    opserr << "Retrieved string: " << type << "\n";

   // opserr << "Dovde" << endln;
    //// Now you want to call getMass and print the mass matrix

    opserr << "OPS f END";
    delete theSectionI;
    delete theSectionJ;
    // If the element is successfully created, return it
    return theEle;

    // Catch-all return to ensure there's always a return statement
    //return nullptr; // Or return 0; if you prefer, but nullptr is clearer for pointers

}

ArchMacroElement3D::ArchMacroElement3D(int tag, int nd1, int nd2,
    SectionForceDeformation* sI, SectionForceDeformation* sJ,
    Vector _intLength, Vector _intLengthMasses,
    CrdTransf& coordTransf, double dt, double theta, const char* type)
    : Element(tag, 0), numSections(2), theSections(0),
    crdTransf(0), intLength(_intLength), intLengthMasses(_intLengthMasses), connectedExternalNodes(2),
    Q(12), q(6), dt(dt), theta(theta), type(type), parameterID(0)

{

    // Allocate arrays of pointers to SectionForceDeformations
    theSections = new SectionForceDeformation * [numSections];

    if (theSections == 0) {
        opserr << "ArchMacroElement3D::ArchMacroElement3D - failed to allocate section model pointer\n";
        exit(-1);
    }

    theSections[0] = sI->getCopy();       if (theSections[0] == 0) { opserr << "ArchMacroElement3D::ArchMacroElement3D -- failed to get a copy of section model\n";  exit(-1); }
    theSections[1] = sJ->getCopy();       if (theSections[1] == 0) { opserr << "ArchMacroElement3D::ArchMacroElement3D -- failed to get a copy of section model\n";  exit(-1); }

    //FIXME: OVde, u ME imas comittedDeformations, da li ti treba ovde? 

    crdTransf = coordTransf.getCopy3d();

    if (crdTransf == 0) {
        opserr << "ArchMacroElement3D::ArchMacroElement3D - failed to copy coordinate transformation\n";
        exit(-1);
    }

    // Set connected external node IDs
    connectedExternalNodes(0) = nd1;
    connectedExternalNodes(1) = nd2;

    theNodes[0] = 0;
    theNodes[1] = 0;

    q0[0] = 0.0; // Zasto ovo nisam imala pre? 
    q0[1] = 0.0;
    q0[2] = 0.0;
    q0[3] = 0.0;
    q0[4] = 0.0;

    p0[0] = 0.0;
    p0[1] = 0.0;
    p0[2] = 0.0;
    p0[3] = 0.0;
    p0[4] = 0.0;

    initializePython();
}

// TODO: Da li damping u zagradi 0? - Ovde nije isti red kao pre
ArchMacroElement3D::ArchMacroElement3D()
    : Element(0, 0), numSections(0), theSections(0), crdTransf(0), intLength(2), intLengthMasses(2),
    connectedExternalNodes(2), Q(6), q(3), dt(0.0), theta(0.0), type(type), parameterID(0)
{
    q0[0] = 0.0;
    q0[1] = 0.0;
    q0[2] = 0.0;
    q0[3] = 0.0;
    q0[4] = 0.0;

    p0[0] = 0.0;
    p0[1] = 0.0;
    p0[2] = 0.0;
    p0[3] = 0.0;
    p0[4] = 0.0;

    theNodes[0] = 0;
    theNodes[1] = 0;


    initializePython();
}


ArchMacroElement3D::~ArchMacroElement3D()
{
    for (int i = 0; i < numSections; i++) {
        if (theSections[i])
            delete theSections[i];
    }


    // Delete the array of pointers to SectionForceDeformation pointer arrays
    if (theSections)
        delete[] theSections;

    if (crdTransf)
        delete crdTransf;
    //crdTransf = nullptr;



    Py_Finalize(); //TODO: DA li mi rteba ovo u deconstructoru? 
}


int
ArchMacroElement3D::getNumExternalNodes() const {
    return 2;
}

const ID&
ArchMacroElement3D::getExternalNodes() {
    return connectedExternalNodes;
}

Node**
ArchMacroElement3D::getNodePtrs() {
    return theNodes;
}

int
ArchMacroElement3D::getNumDOF() {
    return 12;
}

void convertStringToPyString(const char* str, PyObject** pStr) {
    *pStr = PyUnicode_FromString(str);
    if (!*pStr) {
        PyErr_Print();
        opserr << "Error: Failed to convert C-string to PyObject" << endln;
    }
}


int getPythonFunction(PyObject* pModule, const char* functionName, PyObject** pFunction) {
    *pFunction = PyObject_GetAttrString(pModule, functionName);

    // Check the returned object directly.
    if (!*pFunction || !PyCallable_Check(*pFunction)) {
        PyErr_Print();
        Py_Finalize();
        return 1;
    }

    return 0;
}


void convertVectorToPyList(const Vector& v, PyObject** pList) {
    *pList = PyList_New(v.Size());
    if (!*pList) {
        return; // Or handle the error
    }

    for (int i = 0; i < v.Size(); i++) {
        PyObject* item = PyFloat_FromDouble(v[i]);
        if (!item) {
            // Handle the error, maybe release the list and return.
        }
        PyList_SetItem(*pList, i, item);
    }
}


int ArchMacroElement3D::initializePython() {
    // Setup Python process
    Py_SetPythonHome(L"C:\\Users\\ivana\\AppData\\Local\\Programs\\Python\\Python38");

    Py_Initialize();

    // Print the Python version
    PyRun_SimpleString("import sys; print('Python version:', sys.version)");

    // Add current working directory to path
    //PyRun_SimpleString("import sys");

    PyRun_SimpleString("sys.path.append(\".\")");
    // Modify this path as needed
    // Import Python module
    std::cout << "Importing module" << std::endl;
    pModule = PyImport_ImportModule("ArchMacroElement3D");
    if (!pModule) {
        PyErr_Print();
        Py_Finalize();
        return 1;
    }

    // Get references to Python functions
    if (getPythonFunction(pModule, "update", &pUpdateFunc) != 0) return 1;
    if (getPythonFunction(pModule, "get_resisting_force", &pGetResistingForceFunc) != 0) return 1;
    if (getPythonFunction(pModule, "get_initial_stiff", &pGetInitialStiffFunc) != 0) return 1;
    //FIXME: Ja sam dodala liniju ispod, ne ENACIT
    if (getPythonFunction(pModule, "commit_state_in_python", &pCommitFunc) != 0) return 1;
    if (getPythonFunction(pModule, "get_mass_mom_inertia", &pGetMassMomInertiaFunc) != 0) return 1;
    if (getPythonFunction(pModule, "set_domain", &pSetDomainFunc) != 0) return 1;

    if (getPythonFunction(pModule, "pass_dt_theta", &pPassDtThetaFunc) != 0) return 1;
    return 0; // Success: This ensures that the function always has a return value at the end. 
}


void
ArchMacroElement3D::setDomain(Domain* theDomain) {
    // Check Domain is not null - invoked when object removed from a domain
    opserr << "setDomain f START";
    if (theDomain == 0) {
        theNodes[0] = 0;
        theNodes[1] = 0;
        return;
    }

    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);


    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);

    if (theNodes[0] == 0 || theNodes[1] == 0) {
        opserr << "WARNING ArchMacroElement3D (tag: %d), node not found in domain" << this->getTag() << endln;;
        return;
    }

    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();

    if (dofNd1 != 6 || dofNd2 != 6) {
        opserr << "FATAL ERROR ArchMacroElement3D (tag: %d), has differing number of DOFs at its nodes",
            this->getTag();
        return;
    }

    if (crdTransf->initialize(theNodes[0], theNodes[1])) {
        opserr << "Error in initializing crdTransf" << endln;
    }

    double L = crdTransf->getInitialLength();

    if (L == 0.0) {
        opserr << "Init.length is zero. Check crdTransf" << endln;
    }


    opserr << "Set domain print intLength frm OPS: " << intLength << "\n";
    opserr << "Set domain print type frm OPS: " << type << "\n";

    ////////// PASS DT I THETA
              // Create arguments for python function
    opserr << "Set domain print dt frm OPS: " << dt << "\n";
    opserr << "Set domain print theta frm OPS: " << theta << "\n";

    Vector dt_vector(2);  // Create a vector of size 2
    dt_vector[0] = dt;    // Assign dt to the first element
    dt_vector[1] = 0;     // Set the second element to 0

    Vector theta_vector(2);  // Similarly for theta
    theta_vector[0] = theta;
    theta_vector[1] = 0;


    PyObject* pDt;
    PyObject* pTheta;
    convertVectorToPyList(dt_vector, &pDt);  // Convert the dt vector
    convertVectorToPyList(theta_vector, &pTheta);  // Convert the theta vector

    PyObject* pArgss = PyTuple_Pack(2, pDt, pTheta);
    if (!pArgss) {
        opserr << "Error: Failed to pack dt and theta into Python tuple. Check dt and theta values." << endln;
        PyErr_Print();  // Print the detailed Python error message
        // Handle error, possibly return or clean up
    }

    // Call python function
    PyObject* result = PyObject_CallObject(pPassDtThetaFunc, pArgss);
    if (!result) {
        PyErr_Print();  // Handle failure
        // Clean up
    }

    // Clean up
    Py_DECREF(pArgss);
    Py_DECREF(pDt);
    Py_DECREF(pTheta);





    ////////////// PASS PATH WHERE GEOMETRY IS
              // Create arguments for python function
    PyObject* pType;
    
    convertStringToPyString(type, &pType);
    
    if (!pType) {
        opserr << "Failed to convert type to Python string." << endln;
         // or handle error appropriately
    }

    // Assume pUpdateFunc is already defined and is a PyObject pointing to a Python function
    PyObject* pArgs = PyTuple_Pack(1, pType); // Packing just one argument here, expand as needed

   
    if (!pArgs) {
        opserr << "SetDomain error. Check geometry path." << endln;
        
        PyErr_Print();
        //err = 1;
    }

    // Call the Python function
    PyObject_CallObject(pSetDomainFunc, pArgs);
   

    // Clean up Python objects to prevent memory leaks
    Py_DECREF(pArgs);
    Py_DECREF(pType);
    
   

    ////////////////// PASS THETA AND DT




    this->DomainComponent::setDomain(theDomain);
    opserr << "setDomain f END -before this->update";
    this->update();
    opserr << "setDomain f END";
}






int
ArchMacroElement3D::commitState() {
    // I changed this to be as in ElasticBeam2d, not as before when it looped over int points to get stresses, bcs we do not
    // go to the level of stresses.
    opserr << "commitState f START";
    int retVal = 0;

    // call element commitState to do any base class stuff
    // Original commitState functionality
    if ((retVal = this->Element::commitState()) != 0) {
        opserr << "ArchMacroElement3D::commitState() - failed in base class";
        //return retVal; // Handle the error as needed //FIXME: In original code first it goes crdTransf and then return.
    }
    retVal += crdTransf->commitState();

    // Define success and failure Vectors
    Vector successVec(2); // Vector of size 2
    Vector failureVec(2); // Vector of size 2

    // Populate the Vectors with appropriate values
    if (retVal == 0) {
        successVec(0) = 1;
        successVec(1) = 1;
        opserr << "ArchMacroElement3D::commitState() - success Commited";
    }
    else {
        successVec(0) = 0;
        successVec(1) = 0;
    }

    // Convert Vectors to Python lists
    PyObject* pSuccessVec;
    convertVectorToPyList(successVec, &pSuccessVec);
    //convertVectorToPyList(failureVec, &pFailureVec);

    // Create a Python tuple from the lists
    PyObject* pArgs = PyTuple_Pack(1, pSuccessVec);

    // Check for errors in creating Python arguments
    if (!pArgs) {
        PyErr_Print();
        retVal = 1; // Handle error as needed
    }
    else {
        // Call Python function with the arguments
        PyObject_CallObject(pCommitFunc, pArgs);

        // Clean up Python objects
        Py_DECREF(pArgs);
        Py_DECREF(pSuccessVec);
        //Py_DECREF(pFailureVec);
    }
    opserr << "commitState f END";
    return retVal;
}



int
ArchMacroElement3D::revertToLastCommit()
{
    opserr << "revertToLastCommit f START";
    int retVal = 0;

    // Loop over the integration points and revert to last committed state
    for (int i = 0; i < numSections; i++)
        retVal += theSections[i]->revertToLastCommit();

    retVal += crdTransf->revertToLastCommit();

    opserr << "revertToLastCommit f END";
    return retVal;
}

// FIXME: Dodata da bude isto kao FV, iako mi ne treba u sustini
int
ArchMacroElement3D::revertToStart()
{
    opserr << "revertToStart f START";
    int retVal = 0;

    // Loop over the integration points and revert states to start
    for (int i = 0; i < numSections; i++)
        retVal += theSections[i]->revertToStart();

    retVal += crdTransf->revertToStart();
    opserr << "revertToStart f END";
    return retVal;
}



void convertPyListToVector(PyObject* pList, Vector& v) {
    v.resize(PyList_Size(pList));
    for (int i = 0; i < v.Size(); i++) {
        PyObject* item = PyList_GetItem(pList, i);
        if (!item) {
            opserr << "Error: Failed to get list item" << endln;
            return;
        }
        v[i] = PyFloat_AsDouble(item);
        if (PyErr_Occurred()) {
            opserr << "Error: Failed to convert list item to float" << endln;
            return;
        }
    }
}









int
ArchMacroElement3D::update(void)
{
    int err = 0;
    opserr << "Update f START";
    //// First update the coordinate transformation
    int result = crdTransf->update();
    if (result != 0) {
        opserr << "Error updating coordinate transformation\n";
        return result;  // Return the error code if update was not successful
    }

    // Update the transformation OROGINAL
   // crdTransf->update();

    /// Get basic deformations
    //const Vector& disp = crdTransf->getBasicTrialDisp();
    //const Vector& accel = crdTransf->getBasicTrialAccel();

    //// Print the values of disp and accel
    //opserr << "Basic Trial Displacements: " << disp;
    //opserr << "Basic Trial Accelerations: " << accel;
    // 
    const Vector& dispI = theNodes[0]->getTrialDisp();
    const Vector& dispJ = theNodes[1]->getTrialDisp();

    //// Print the values of disp and accel
    //opserr << "Basic Trial Displacements: " << disp;
    //opserr << "Basic Trial Accelerations: " << accel;
    const Vector& accelI = theNodes[0]->getTrialAccel();
    const Vector& accelJ = theNodes[1]->getTrialAccel();

    opserr << "Trial Accelerations - accelI: " << accelI;
    opserr << "Trial Accelerations - accelJ: " << accelJ;

    //opserr << "Basic Trial Accelerations: " << accel;
    // Calculate the norm of the trial accelerations
    double accelNorm = accelI.Norm();

    // Check if the norm is greater than DBL_EPSILON
    if (accelNorm > DBL_EPSILON) {
        opserr << "The norm of trial accelerations is bigger than DBL_EPSILON.\n";
    }
    else {
        opserr << "The norm of trial accelerations is too small (less than or equal to DBL_EPSILON).\n";
    }

    double L = crdTransf->getInitialLength();
    double oneOverL = 1.0 / L;


    //FIXME: IT4R: Pass to my python getBasicTrialAcc/Displ 

          // Create arguments for python function
    PyObject* pDisp;
    PyObject* pAccel;
    convertVectorToPyList(dispI, &pDisp);
    convertVectorToPyList(accelI, &pAccel);
    PyObject* pArgs = PyTuple_Pack(2, pDisp, pAccel);
    if (!pArgs) {
        PyErr_Print();
        err = 1;
    }

    // Call python function
    PyObject_CallObject(pUpdateFunc, pArgs);

    // Clean up
    Py_DECREF(pArgs);
    Py_DECREF(pDisp);
    Py_DECREF(pAccel);

    //opserr << "Update invoked  " << disp;
    opserr << "Update f END";
    return 0;
}


const Matrix&
ArchMacroElement3D::getInitialBasicStiff()
{

    opserr << "getInitialBasicStiff f START";
    static Matrix kb(6, 6); //TODO: PROMENILA SI NA 6,6

    kb.Zero();

    // Call Python function - Ja sam dodala isto kao sto su ENACIT u INitial stiffness fji. Proveri da li radi - rteba da bude 6x6
    PyObject* pResult = PyObject_CallObject(pGetInitialStiffFunc, NULL);
    if (!pResult) {
        PyErr_Print();
        return kb;
    }

    // Convert Python object to flat Vector
    Vector kbFlat;
    convertPyListToVector(pResult, kbFlat); //FIXME: IT4R passes vector, and crdTransf needs matrix; so let's make a matrix kbMatr

    // Check if the size matches the expected size of the kb matrix
    if (kbFlat.Size() != 36) {
        opserr << "Error: The size of the flattened array does not match the expected size of the kb matrix." << endln;
        return kb;
    }

    // Populate the kb matrix with values from kbFlat //TODO: PROVERI
    int counter = 0;
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            kb(i, j) = kbFlat(counter++);
        }
    }

    // Print the kb matrix
    opserr << "kb matrix:" << endln;
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            opserr << kb(i, j) << " ";
        }
        opserr << endln; // Move to the next line after printing each row
    }

    opserr << "getInitialBasicStiff f END";
    return kb; //TODO: Ovde treba drugu matricu vratiti - zasto???
}



const Matrix&
ArchMacroElement3D::getInitialStiff()
{

    opserr << "getInitialStiff f START";
    const Matrix& kb = this->getInitialBasicStiff();
    // Call computeElemtLengthAndOrient before getting the initial length
    crdTransf->initialize(theNodes[0], theNodes[1]);

    // Now you can get the initial length
    double L = crdTransf->getInitialLength();
    //opserr << "length is " << L << endln;
    // Transform to global stiffness
    K = crdTransf->getInitialGlobalStiffMatrix(kb);

    // Print the K matrix
    opserr << "K matrix:" << endln;
    for (int i = 0; i < K.noRows(); ++i) { // Assuming K.noRows() gives the number of rows
        for (int j = 0; j < K.noCols(); ++j) { // Assuming K.noCols() gives the number of columns
            opserr << K(i, j) << " ";
        }
        opserr << endln; // New line at the end of each row
    }


    opserr << "getInitialStiff f END";
    return K;
}

//FIXME: Za sada, TangStiff je ista kao InitStiff, u globalnom koord sist ovde. Promeniti? 
//Ova f=ja  mora postojati implementirana, from Element
const Matrix&
ArchMacroElement3D::getTangentStiff()
{
    opserr << "getTangentStiff f START";
    const Matrix& kb = this->getInitialBasicStiff();
    crdTransf->initialize(theNodes[0], theNodes[1]);

    // Now you can get the initial length
    double L = crdTransf->getInitialLength();
    //opserr << "length is " << L << endln;

    // Transform to global stiffness
    K = crdTransf->getInitialGlobalStiffMatrix(kb);

    //Print K tangent matrix
    opserr << "K tangent matrix:" << endln;
    for (int i = 0; i < K.noRows(); ++i) { // Assuming K.noRows() gives the number of rows
        for (int j = 0; j < K.noCols(); ++j) { // Assuming K.noCols() gives the number of columns
            opserr << K(i, j) << " ";
        }
        opserr << endln; // New line at the end of each row
    }
    opserr << "getTangentStiff f END";
    return K;
}


const Matrix&
ArchMacroElement3D::getMass()
{
    //opserr << "Python sys.path before calling function:\n";

    // Import the sys module and print the Python path to check
    // PyRun_SimpleString("import sys");
    // PyRun_SimpleString("print(sys.path)");

    //NEW VERSION:
    opserr << "getMass f START";
    K.Zero();

    // Call Python function to get mass and moment of inertia
    PyObject* pResult = PyObject_CallObject(pGetMassMomInertiaFunc, NULL);
    if (!pResult) {
        PyErr_SetString(PyExc_RuntimeError, "getMass f - Error calling Python function.");
        PyErr_Print();
        return K;
    }

    // Convert Python object to Vector
    Vector massMomInertia;
    convertPyListToVector(pResult, massMomInertia);

    // Print the vector's elements
    opserr << "massMomInertia contents: ";
    for (int i = 0; i < massMomInertia.Size(); ++i) {
        opserr << massMomInertia[i] << " ";
    }
    opserr << endln;

    // Check that the vector has the correct size (6 elements)
    if (massMomInertia.Size() != 12) {
        PyErr_SetString(PyExc_RuntimeError, "Python function did not return 12 elements.");
        PyErr_Print();
        return K;
    }

    // Populate the K matrix with the values returned from Python
    /*for (int i = 0; i < 12; ++i) {
        K(i, i) = massMomInertia(i);
    }*/
    K(0, 0) = massMomInertia(0);
    K(1, 1) = massMomInertia(1);
    K(2, 2) = massMomInertia(2);

    K(6, 6) = massMomInertia(6);
    K(7, 7) = massMomInertia(7);
    K(8, 8) = massMomInertia(8);

    // Print the K matrix
    opserr << "K matrix - actually massMatrix:\n";
    for (int i = 0; i < K.noRows(); ++i) {
        for (int j = 0; j < K.noCols(); ++j) {
            opserr << K(i, j) << " ";
        }
        opserr << endln; // New line at the end of each row
    }
    opserr << "getMass f END";
    return K;

}
// Ovde imas pet sila externih, jer nemas za torziju
void
ArchMacroElement3D::zeroLoad(void)
{
    Q.Zero();

    q0[0] = 0.0;
    q0[1] = 0.0;
    q0[2] = 0.0;
    q0[3] = 0.0;
    q0[4] = 0.0;

    p0[0] = 0.0;
    p0[1] = 0.0;
    p0[2] = 0.0;
    p0[3] = 0.0;
    p0[4] = 0.0;

    return;
}


int
ArchMacroElement3D::addLoad(ElementalLoad* theLoad, double loadFactor) {

    opserr << "addLoad f START";
    int type;
    const Vector& data = theLoad->getData(type, loadFactor);
    double L = crdTransf->getInitialLength();
    opserr << "Load factor bego of addLoad f: " << loadFactor << endln;

    if (type == LOAD_TAG_Beam3dUniformLoad) {
        double wy = data(0) * loadFactor;  // Transverse
        double wz = data(1) * loadFactor;  // Transverse
        double wx = data(2) * loadFactor;  // Axial (+ve from node I to J)

        double Vy = 0.5 * wy * L;
        double Mz = Vy * L / 6.0; // wy*L*L/12
        double Vz = 0.5 * wz * L;
        double My = Vz * L / 6.0; // wz*L*L/12
        double P = wx * L;

        // Reactions in basic system
        p0[0] -= P;
        p0[1] -= Vy;
        p0[2] -= Vy;
        p0[3] -= Vz;
        p0[4] -= Vz;

        // Fixed end forces in basic system
        q0[0] -= 0.5 * P;
        q0[1] -= Mz;
        q0[2] += Mz;
        q0[3] += My;
        q0[4] -= My;
    }

    else if (type == LOAD_TAG_Beam3dPointLoad) {
        double Py = data(0) * loadFactor;
        double Pz = data(1) * loadFactor;
        double N = data(2) * loadFactor;
        double aOverL = data(3);
        opserr << "Load factor inside point load: " << loadFactor << endln;

        if (aOverL < 0.0 || aOverL > 1.0)
            return 0;

        double a = aOverL * L;
        double b = L - a;

        // Reactions in basic system
        p0[0] -= N; //Subtract the value of V1 from the current value of p0[3]
        double V1, V2;
        V1 = Py * (1.0 - aOverL);
        V2 = Py * aOverL;
        p0[1] -= V1;
        p0[2] -= V2;
        V1 = Pz * (1.0 - aOverL);
        V2 = Pz * aOverL;
        p0[3] -= V1;
        p0[4] -= V2;

        double L2 = 1.0 / (L * L);
        double a2 = a * a;
        double b2 = b * b;

        // Fixed end forces in basic system
        q0[0] -= N * aOverL;
        double M1, M2;
        M1 = -a * b2 * Py * L2;
        M2 = a2 * b * Py * L2;
        q0[1] += M1;
        q0[2] += M2;
        M1 = -a * b2 * Pz * L2;
        M2 = a2 * b * Pz * L2;
        q0[3] -= M1;
        q0[4] -= M2;
    }
    else {
        opserr << "ArchMacroElement3D::addLoad()  -- load type unknown for element with tag: " << this->getTag() << endln;
        return -1;
    }
}



int
ArchMacroElement3D::addInertiaLoadToUnbalance(const Vector& accel) {
    opserr << "addInertiaLoadToUnbalance f START";
    const Matrix& MM = this->getMass();
    double m = MM(0, 0);

    // Print the value of m to check it
    opserr << "Value of m  in addInertLoadtoUn function: " << m << endln;

    // Get R * accel from the nodes
    const Vector& Raccel1 = theNodes[0]->getRV(accel);
    const Vector& Raccel2 = theNodes[1]->getRV(accel);

    // Print the values of disp and accel
    opserr << "addInLoadToUnb getRV: " << Raccel1;
    opserr << "addInLoadToUnb getRV: " << Raccel2;


    if (6 != Raccel1.Size() || 6 != Raccel2.Size()) {
        opserr << "ArchMacroElement3D::addInertiaLoadToUnbalance matrix and vector sizes are incompatible\n";
        return -1;
    }

    //TODO: Proveri jer si ovo uzela iz DispBasedBeam2d; nemas u rotacijama nista ovako. Da li treba isto i u 3d? 
    Q(0) -= m * Raccel1(0);
    Q(1) -= m * Raccel1(1);
    Q(2) -= m * Raccel1(2);

    Q(6) -= m * Raccel2(0);
    Q(7) -= m * Raccel2(1);
    Q(8) -= m * Raccel2(2);
    // Print the Q vector
    opserr << "Q vector: ";
    for (int i = 0; i < Q.Size(); ++i) {
        opserr << Q(i) << " ";
    }
    opserr << endln; // New line at the end

    opserr << "addInertiaLoadToUnbalance f END";

    return 0;
}


const Vector&
ArchMacroElement3D::getResistingForceIncInertia()
{
    opserr << "getResistingForceIncInertia f START";

    P = this->getResistingForce();

    // Subtract other external nodal loads ... P_res = P_int - P_ext
    P.addVector(1.0, Q, -1.0);

    // add the damping forces if rayleigh damping
    if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
        P.addVector(1.0, this->getRayleighDampingForces(), 1.0);



    const Vector& accel1 = theNodes[0]->getTrialAccel();
    const Vector& accel2 = theNodes[1]->getTrialAccel();

    // Print the values of disp and accel
    opserr << "addInLoadToUnb getRV: " << accel1;
    opserr << "addInLoadToUnb getRV: " << accel2;


    const Matrix& MM = this->getMass();
    double m = MM(0, 0);

    // Print the value of m to check it
    opserr << "Value of m inside getResForceIncInertia: " << m << endln;
    //double m = 477.0;//Stavi da procita iz mass MAtrix
    P(0) += m * accel1(0);
    P(1) += m * accel1(1);
    P(2) += m * accel1(2);

    P(6) += m * accel2(0);
    P(7) += m * accel2(1);
    P(8) += m * accel2(2);


    opserr << "getResistingForceIncInertia f END";
    return P;
}



const Vector&
ArchMacroElement3D::getResistingForce()
{
    opserr << "getResistingForce f START";
    PyObject* pResult = PyObject_CallMethod(pModule, "get_resisting_force", NULL);
    if (!pResult) {
        PyErr_Print();
    }

    // Convert Python object to Vector
    Vector resistingForce;
    convertPyListToVector(pResult, resistingForce);
    opserr << "resisting Force is - from Py " << resistingForce << endln;
    //FIXME: OVO JE BILO PRE< KAO GLOBAL FORCES IZ PY, RADILO JE
    // transform resisting forces  from local to global coordinates

    //Matrix plMatrix(12, 1); // A Matrix with 12 rows and 1 column.
    //// Assuming indices for resistingForce are correctly mapped to plMatrix elements:
    //plMatrix(0, 0) = resistingForce[6];
    //plMatrix(1, 0) = resistingForce[7];
    //plMatrix(2, 0) = resistingForce[8];
    //plMatrix(3, 0) = resistingForce[9];
    //plMatrix(4, 0) = resistingForce[10];
    //plMatrix(5, 0) = resistingForce[11];

    //plMatrix(6, 0) = resistingForce[0];
    //plMatrix(7, 0) = resistingForce[1];
    //plMatrix(8, 0) = resistingForce[2];
    //plMatrix(9, 0) = resistingForce[3];
    //plMatrix(10, 0) = resistingForce[4];
    //plMatrix(11, 0) = resistingForce[5];


    //// Print each element of P_matrix
    //opserr << "plMatrix Contents:" << endln;
    //for (int i = 0; i < 12; ++i) {
    //    opserr << "plMatrix(" << i << ", 0): " << plMatrix(i, 0) << endln;
    //}


    //Matrix P_matrix = crdTransf->getGlobalMatrixFromLocal(plMatrix);

    //// Print each element of P_matrix
    //opserr << "P_matrix Contents:" << endln;
    //for (int i = 0; i < 12; ++i) {
    //    opserr << "P_matrix(" << i << ", 0): " << P_matrix(i, 0) << endln;
    //}


    //// Initialize a Vector to hold the result, assuming 'P' has 12 rows and 1 column
    //Vector resultVector(12);

    //// Copy each element from 'P' back into the vector
    //for (int i = 0; i < 12; ++i) {
    //    resultVector(i) = P_matrix(i, 0); // Assuming this is the correct way to access elements in your Vector and Matrix classes
    //}
    //Vector P = resultVector;

    Vector pg(12);


    pg(0) = -resistingForce[6];
    pg(1) = -resistingForce[7];
    pg(2) = -resistingForce[8];
    pg(3) = -resistingForce[9];
    pg(4) = -resistingForce[10];
    pg(5) = -resistingForce[11];

    pg(6) = -resistingForce[0];
    pg(7) = -resistingForce[1];
    pg(8) = -resistingForce[2];
    pg(9) = -resistingForce[3];
    pg(10) = -resistingForce[4];
    pg(11) = -resistingForce[5];

    P = pg;

    opserr << "P  is " << P << endln;

    // Check for NaN values in P
    for (int i = 0; i < P.Size(); ++i) {
        if (std::isnan(P(i))) {
            opserr << "Error: Resisting force contains NaN values at index " << i << ". Stopping analysis." << endln;
            return -1;  // Return an error code when NaN values are found
        }
    }

    opserr << "getResistingForce f END";


    return P;
    //FIXME: LOCAL FORCES IZ PY:
   //Vector pl(12);

   //pl(0) = resistingForce[3];
   //pl(1) = resistingForce[4];
   //pl(2) = resistingForce[5];
   //pl(3) = resistingForce[0];
   //pl(4) = resistingForce[1];
   //pl(5) = resistingForce[2];
   //// transform resisting forces  from local to global coordinates
   //Vector pg(6);
   //double Angle = 0.0;

   //pg(0) = cos(Angle) * pl[0] - sin(Angle) * pl[1];
   //pg(1) = sin(Angle) * pl[0] + cos(Angle) * pl[1];

   //pg(3) = cos(Angle) * pl[3] - sin(Angle) * pl[4];
   //pg(4) = sin(Angle) * pl[3] + cos(Angle) * pl[4];

   //pg(2) = pl[2];
   //pg(5) = pl[5];

   //P = pg;
   //opserr << "P from getResistingForces is " << P << endln;
   //opserr << "getResistingForce f END";
   //return P;

}

// sendSelf, recvSelf, There are two methods provided which are required is the user uses to use the database or parallel procesing features of the OpenSees applications.If neither are to be used, the developer need simply return a negative value in both methods


int
ArchMacroElement3D::sendSelf(int /*commitTag*/, Channel&/*theChannel*/)
{
    opserr << "VaultMacroElement2D::sendSelf: method not implemented\n";
    return -1;
}

int
ArchMacroElement3D::recvSelf(int /*commitTag*/, Channel&/*theChannel*/,
    FEM_ObjectBroker&/*theBroker*/)
{
    opserr << "VaultMacroElement2D::recvSelf: method not implemented\n";
    return -1;
}


void
ArchMacroElement3D::Print(OPS_Stream& s, int /*flag*/)
{
    s << "\nVaultMacroElement2D, element id:  " << this->getTag() << endln;
    s << "\tConnected external nodes:  " << connectedExternalNodes;

}


Response*
ArchMacroElement3D::setResponse(const char** argv, int argc, OPS_Stream& output)
{

    Response* theResponse = 0;

    output.tag("ElementOutput");
    output.attr("eleType", "DispBeamColumn3d");
    output.attr("eleTag", this->getTag());
    output.attr("node1", connectedExternalNodes[0]);
    output.attr("node2", connectedExternalNodes[1]);

    //
    // we compare argv[0] for known response types 
    //

    // global force - 
    if (strcmp(argv[0], "forces") == 0 || strcmp(argv[0], "force") == 0
        || strcmp(argv[0], "globalForce") == 0 || strcmp(argv[0], "globalForces") == 0) {

        output.tag("ResponseType", "Px_1");
        output.tag("ResponseType", "Py_1");
        output.tag("ResponseType", "Pz_1");
        output.tag("ResponseType", "Mx_1");
        output.tag("ResponseType", "My_1");
        output.tag("ResponseType", "Mz_1");
        output.tag("ResponseType", "Px_2");
        output.tag("ResponseType", "Py_2");
        output.tag("ResponseType", "Pz_2");
        output.tag("ResponseType", "Mx_2");
        output.tag("ResponseType", "My_2");
        output.tag("ResponseType", "Mz_2");


        theResponse = new ElementResponse(this, 1, P);

        // local force -
    }
    else if (strcmp(argv[0], "localForce") == 0 || strcmp(argv[0], "localForces") == 0) {

        output.tag("ResponseType", "N_1");
        output.tag("ResponseType", "Vy_1");
        output.tag("ResponseType", "Vz_1");
        output.tag("ResponseType", "T_1");
        output.tag("ResponseType", "My_1");
        output.tag("ResponseType", "Mz_1");
        output.tag("ResponseType", "N_2");
        output.tag("ResponseType", "Vy_2");
        output.tag("ResponseType", "Vz_2");
        output.tag("ResponseType", "T_2");
        output.tag("ResponseType", "My_2");
        output.tag("ResponseType", "Mz_2");

        theResponse = new ElementResponse(this, 2, P);
        output.endTag();
        if (theResponse == 0)
            return Element::setResponse(argv, argc, output);
        else
            return theResponse;
    }
}


int
ArchMacroElement3D::getResponse(int responseID, Information& eleInfo)
{
    double N, V, M1, M2, T;
    double L = crdTransf->getInitialLength();
    double oneOverL = 1.0 / L;

    if (responseID == 1)
        return eleInfo.setVector(this->getResistingForce());

    else if (responseID == 2) {
        // Axial
        N = q(0);
        P(6) = N;
        P(0) = -N + p0[0];

        // Torsion
        T = q(5);
        P(9) = T;
        P(3) = -T;

        // Moments about z and shears along y
        M1 = q(1);
        M2 = q(2);
        P(5) = M1;
        P(11) = M2;
        V = (M1 + M2) * oneOverL;
        P(1) = V + p0[1];
        P(7) = -V + p0[2];

        // Moments about y and shears along z
        M1 = q(3);
        M2 = q(4);
        P(4) = M1;
        P(10) = M2;
        V = (M1 + M2) * oneOverL;
        P(2) = -V + p0[3];
        P(8) = V + p0[4];

        return eleInfo.setVector(P);
    }
    else
        return -1;
} // Return -1 to indicate failure or that no action was performed

int
ArchMacroElement3D::setRayleighDampingFactors(double alpham, double betak, double betak0, double betakc)
{
    opserr << "setRayleighDampingFactors f START";
    alpham = 0.0;
    betak = 0.0;
    betak0 = 0.0;
    betakc = 0.0;

    // alphaM = 0.0;
       //betaK = 0.0;
       //betaK0 = 0.0;
    // betaKc = 0.0;
    alphaM = alpham;
    betaK = betak;
    betaK0 = betak0;
    betaKc = betakc;

    // check that memory has been allocated to store compute/return
    // damping matrix & residual force calculations
    if (index == -1) {
        int numDOF = this->getNumDOF();

        for (int i = 0; i < numMatrices; i++) {
            Matrix* aMatrix = theMatrices[i];
            if (aMatrix->noRows() == numDOF) {
                index = i;
                i = numMatrices;
            }
        }
        if (index == -1) {
            Matrix** nextMatrices = new Matrix * [numMatrices + 1];
            if (nextMatrices == 0) {
                opserr << "Element::getTheMatrix - out of memory\n";
            }
            int j;
            for (j = 0; j < numMatrices; j++)
                nextMatrices[j] = theMatrices[j];
            Matrix* theMatrix = new Matrix(numDOF, numDOF);
            if (theMatrix == 0) {
                opserr << "Element::getTheMatrix - out of memory\n";
                exit(-1);
            }
            nextMatrices[numMatrices] = theMatrix;

            Vector** nextVectors1 = new Vector * [numMatrices + 1];
            Vector** nextVectors2 = new Vector * [numMatrices + 1];
            if (nextVectors1 == 0 || nextVectors2 == 0) {
                opserr << "Element::getTheVector - out of memory\n";
                exit(-1);
            }

            for (j = 0; j < numMatrices; j++) {
                nextVectors1[j] = theVectors1[j];
                nextVectors2[j] = theVectors2[j];
            }

            Vector* theVector1 = new Vector(numDOF);
            Vector* theVector2 = new Vector(numDOF);
            if (theVector1 == 0 || theVector2 == 0) {
                opserr << "Element::getTheVector - out of memory\n";
                exit(-1);
            }

            nextVectors1[numMatrices] = theVector1;
            nextVectors2[numMatrices] = theVector2;

            if (numMatrices != 0) {
                delete[] theMatrices;
                delete[] theVectors1;
                delete[] theVectors2;
            }
            index = numMatrices;
            numMatrices++;
            theMatrices = nextMatrices;
            theVectors1 = nextVectors1;
            theVectors2 = nextVectors2;
        }
    }

    // if need storage for Kc go get it
    if (betaKc != 0.0) {
        if (Kc == 0)
            Kc = new Matrix(this->getTangentStiff());
        if (Kc == 0) {
            opserr << "WARNING - ELEMENT::setRayleighDampingFactors - out of memory\n";
            betaKc = 0.0;
        }

        // if don't need storage for Kc & have allocated some for it, free the memory
    }
    else if (Kc != 0) {
        delete Kc;
        Kc = 0;
    }
    opserr << "setRayleighDampingFactors f END";
    return 0;
}




const Matrix&
ArchMacroElement3D::getDamp(void)
{
    //// Retrieve the current matrix
    //Matrix *theMatrix = theMatrices[index];
    //
    //// Clear the matrix
    //theMatrix->Zero();

    opserr << "getDamp f START";
    if (index == -1) {
        this->setRayleighDampingFactors(alphaM, betaK, betaK0, betaKc);
    }

    // now compute the damping matrix
    Matrix* theMatrix = theMatrices[index];
    theMatrix->Zero();
    if (alphaM != 0.0)
        theMatrix->addMatrix(0.0, this->getMass(), alphaM);
    if (betaK != 0.0)
        theMatrix->addMatrix(1.0, this->getTangentStiff(), betaK);
    if (betaK0 != 0.0)
        theMatrix->addMatrix(1.0, this->getInitialStiff(), betaK0);
    if (betaKc != 0.0)
        theMatrix->addMatrix(1.0, *Kc, betaKc);

    //// Print the damping matrix
    //opserr << "Damping Matrix: \n";
    //for (int i = 0; i < theMatrix->noRows(); i++) {
    //    for (int j = 0; j < theMatrix->noCols(); j++) {
    //        opserr << (*theMatrix)(i, j) << " ";
    //    }
    //    opserr << "\n"; // New line for each row
    //}

    opserr << "getDamp f END";
    // return the computed matrix
    return *theMatrix;
}





const Vector&
ArchMacroElement3D::getRayleighDampingForces(void)
{
    opserr << "getRayleighDampingForces f START";
    if (index == -1) {
        this->setRayleighDampingFactors(alphaM, betaK, betaK0, betaKc);
    }

    Matrix* theMatrix = theMatrices[index];
    Vector* theVector = theVectors2[index];
    Vector* theVector2 = theVectors1[index]; //velocity data from nodes stored in theVector2 


    // determine the vel vector from ele nodes
    Node** theNodes = this->getNodePtrs();
    int numNodes = this->getNumExternalNodes();
    int loc = 0;
    for (int i = 0; i < numNodes; i++) {
        const Vector& vel = theNodes[i]->getTrialVel();
        for (int i = 0; i < vel.Size(); i++) {
            (*theVector2)(loc++) = vel[i];
        }
    }

    // now compute the damping matrix
    theMatrix->Zero();
    if (alphaM != 0.0)
        theMatrix->addMatrix(0.0, this->getMass(), alphaM);
    if (betaK != 0.0)
        theMatrix->addMatrix(1.0, this->getTangentStiff(), betaK);
    if (betaK0 != 0.0)
        theMatrix->addMatrix(1.0, this->getInitialStiff(), betaK0);
    if (betaKc != 0.0)
        theMatrix->addMatrix(1.0, *Kc, betaKc);

    // finally the D * v
    // computes the product of the damping matrix and the velocity vector (theVector2) and stores the result in theVector. 
    // This operation simulates the effect of damping forces on the system.

    theMatrix->Zero();

    theVector->addMatrixVector(0.0, *theMatrix, *theVector2, 1.0);

    //// Print theVector contents (damping forces)
    //opserr << "Damping Forces - from getRaylDampForces: ";
    //for (int i = 0; i < theVector->Size(); i++) {
    //    opserr << (*theVector)(i) << " ";
    //}
    //opserr << "\n";


    //// Return the cleared vector
    //return *theVector;
    opserr << "getRayleighDampingForces f END";
    return *theVector;
}