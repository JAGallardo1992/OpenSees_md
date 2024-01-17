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
                                                                        
// $Revision: 1.8 $
// $Date: 2007-02-02 01:30:47 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/NewElement.cpp,v $
                                                                        
// Written: fmk 
// Created: 08/01
//
// Description: This file contains the implementation for the NewElement class.
//
// What: "@(#) NewElement.cpp, revA"


#include "TSM_2D.h"

#include <Information.h>
#include <ElementResponse.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <OPS_Globals.h>
#include <OPS_Stream.h>

#include <elementAPI.h>
#include <UniaxialMaterial.h>

#define PI 3.14159l

static int numMyBearing = 0;

Matrix TSM_2D::theMatrix(6, 6);
// Vector TSM_2D::theVector(6);

void* OPS_TSM_2D()
{
	if (numMyBearing == 0) {
		opserr << "TSM element - Written by Jose Gallardo, PUC, 2022\n";
		//opserr << "Based on the paper of Gallardo et al. 2022\n";
		numMyBearing++;
	}

	int ndf = OPS_GetNDF();
	if (ndf != 3) {
		opserr << "WARNING invalid ndf: " << ndf;
		opserr << ", for plane problem need 3 dfs\n";
		return 0;
	}


	if (OPS_GetNumRemainingInputArgs() < 14) {
		opserr << "WARNING insufficient arguments\n";
		opserr << "Want: TSM_2D eleTag iNode jNode kvi fas Kb Kr Do Di Tr fc a kappa phi_m -shear matTag <-orient y1 y2 y3>  <-mass m> \n";
		return 0;
	}

	// tags
	int idata[3];
	int num = 3;
	if (OPS_GetIntInput(&num, idata) < 0) {
		opserr << "WARNING: invalid integer inputs\n";
		return 0;
	}

	// data
	double data[11];
	num = 11;
	if (OPS_GetDoubleInput(&num, data) < 0) {
		opserr << "WARNING: invalid double inputs\n";
		return 0;
	}

	// materials
	UniaxialMaterial* mats[1] = { 0 };
	const char* type = OPS_GetString();
	if (strcmp(type, "-shear") != 0) {
		opserr << "WARNING: want -shear\n";
		return 0;
	}
	int matTag;
	num = 1;
	if (OPS_GetIntInput(&num, &matTag) < 0) {
		opserr << "WARNING: invalid matTag for shear\n";
		return 0;
	}
	mats[0] = OPS_getUniaxialMaterial(matTag);
	if (mats[0] == 0) {
		opserr << "WARNING: material not found\n";
		return 0;
	}

	// options
	Vector y(3); y(0) = 1.0; y(1) = 0.0; y(2) = 0.0;
	double mass = 0.0;

	if (OPS_GetNumRemainingInputArgs() < 1) {
		return new TSM_2D(idata[0], idata[1], idata[2], data[0], data[1], data[2],
			data[3], data[4], data[5], data[6], data[7], data[8], data[9], data[10],
			mats, y, mass);
	}
	
	while (OPS_GetNumRemainingInputArgs() > 0) {
		type = OPS_GetString();
		if (strcmp(type, "-orient") == 0) {
			if (OPS_GetNumRemainingInputArgs() < 2) {
				opserr << "WARNING: insufficient arguments after -orient\n";
				return 0;
			}
			num = 2;
			if (OPS_GetDoubleInput(&num, &y(0)) < 0) {
				opserr << "WARNING: invalid orient value\n";
				return 0;
			}
		}
		
		else if (strcmp(type, "-mass") == 0) {
			if (OPS_GetNumRemainingInputArgs() < 1) {
				opserr << "WARNING: insufficient args\n";
				return 0;
			}
			num = 1;
			if (OPS_GetDoubleInput(&num, &mass) < 0) {
				opserr << "WARNING: invalid mass\n";
				return 0;
			}
		}
		
	}
	
	return new TSM_2D(idata[0], idata[1], idata[2], data[0], data[1], data[2],
		data[3], data[4], data[5], data[6], data[7], data[8], data[9], data[10],
		mats, y, mass);
}

// constructors:
TSM_2D::TSM_2D(int tag, int Nd1, int Nd2,
    double _Kvi, double _fas, double _Kb, double _Kr, double _Do, double _Di,
	double _Tr, double _fc, double _a, double _k, double _pm,
    UniaxialMaterial** Material, const Vector _y, double _mass)
 : Element(tag,ELE_TAG_TSM_2D),
  connectedExternalNodes(2), Kvi(_Kvi), fas(_fas), y(_y), Do(_Do), Di(_Di),
	Tr(_Tr), fc(_fc), a(_a), k(_k), pm(_pm),
    Kb(_Kb), Kr(_Kr), L(0.0), ub(3), qb(3), kt(3, 3), mass(_mass),
    ul(6), Tgl(6, 6), Tlb(3, 6), ub_c(3), kinit(3,3), theLoad(6)
{
	
	// ensure the connectedExternalNode ID is of correct size & set values
	if (connectedExternalNodes.Size() != 2) {
		opserr << "TSM_2D::TSM_2D() - element: "
			<< this->getTag() << " - failed to create an ID of size 2.\n";
		exit(-1);
	}

	connectedExternalNodes(0) = Nd1;
	connectedExternalNodes(1) = Nd2;

	// set node pointers to NULL
	for (int i = 0; i < 2; i++)
		theNodes[i] = 0;

		// check material input
	if (Material == 0) {
		opserr << "TSM_2D::TSM_2D() - null material for shear behavior \n ";
		exit(-1);
	}

	// get copies of the uniaxial materials

	if (Material[0] == 0) {
		opserr << "TSM_2D::TSM_2D() - "
			"null uniaxial material pointer passed.\n";
		exit(-1);
	}


	theMaterials[0] = Material[0]->getCopy();
	// initialize initial stiffness matrix
	kinit.Zero();
	kinit.resize(3, 3);
	kinit(0, 0) = Kvi;
	kinit(1, 1) = theMaterials[0]->getInitialTangent();
	kinit(2, 2) = Kb;

	ub.resize(3);
	qb.resize(3);
	ub_c.resize(3);
	
	// initialize other variables
	this->revertToStart();
	
}

TSM_2D::TSM_2D()
 :Element(0,ELE_TAG_TSM_2D),
  connectedExternalNodes(2), Kvi(0.0), fas(0.0), y(0), Do(0.0), Di(0.0),
	Tr(0.0), fc(0.0), a(0.0), k(0.0), pm(0.0),
	Kb(0.0), Kr(0.0), L(0.0), ub(3), qb(3), kt(3, 3), mass(0.0),
	ul(6), Tgl(6, 6), Tlb(3, 6), ub_c(3), kinit(3, 3), theLoad(6),
	s(0.0), t_s(0.0), theta(0.0), t_theta(0.0), 
	fpc(0.0), t_fpc(0.0), um(0.0), t_um(0.0)
{
	// ensure the connectedExternalNode ID is of correct size
	if (connectedExternalNodes.Size() != 2) {
		opserr << "TSM_2D::TSM_2D() - element: "
			<< this->getTag() << " - failed to create an ID of size 2.\n";
		exit(-1);
	}

	// set node pointers to NULL
	for (int i = 0; i < 2; i++)
		theNodes[i] = 0;

	// set material pointers to NULL
		theMaterials[0] = 0;

}

//  destructor:
TSM_2D::~TSM_2D()
{

}

int TSM_2D::getNumExternalNodes(void) const
{
	//opserr << "getNumExternalNodes \n";
    return 2;
}

const ID & TSM_2D::getExternalNodes(void)
{
	//opserr << "getExternalNodes \n";
    return connectedExternalNodes;
}

Node** TSM_2D::getNodePtrs(void)
{
	//opserr << "getNodePtrs \n";
  return theNodes;
}

int TSM_2D::getNumDOF(void)
{
	//opserr << "getNumDoF \n";
    return 6;
}


void TSM_2D::setDomain(Domain *theDomain)
{
	//opserr << "setDomain \n";
	// check Domain is not null - invoked when object removed from a domain
	if (!theDomain) {
		theNodes[0] = 0;
		theNodes[1] = 0;
		return;
	}

	// first set the node pointers
	theNodes[0] = theDomain->getNode(connectedExternalNodes(0));
	theNodes[1] = theDomain->getNode(connectedExternalNodes(1));

	// if can't find both - send a warning message
	if (!theNodes[0] || !theNodes[1]) {
		if (!theNodes[0]) {
			opserr << "WARNING TSM_2D::setDomain() - Nd1: "
				<< connectedExternalNodes(0)
				<< " does not exist in the model for";
		}
		else {
			opserr << "WARNING TSM_2D::setDomain() - Nd2: "
				<< connectedExternalNodes(1)
				<< " does not exist in the model for";
		}
		opserr << " element: " << this->getTag() << ".\n";

		return;
	}

	// now determine the number of dof and the dimension
	int dofNd1 = theNodes[0]->getNumberDOF();
	int dofNd2 = theNodes[1]->getNumberDOF();

	// if differing dof at the ends - print a warning message
	if (dofNd1 != 3) {
		opserr << "TSM_2D::setDomain() - node 1: "
			<< connectedExternalNodes(0)
			<< " has incorrect number of DOF (not 3).\n";
		return;
	}
	if (dofNd2 != 3) {
		opserr << "TSM_2D::setDomain() - node 2: "
			<< connectedExternalNodes(1)
			<< " has incorrect number of DOF (not 3).\n";
		return;
	}

	// call the base class method
    this->DomainComponent::setDomain(theDomain);
	// set up the transformation matrix for orientation

	this->setUp();
}   	 

int TSM_2D::commitState()
{
	//opserr << "CommitState \n";
	int retVal = 0;

	s = t_s;
	theta = t_theta;
	ub_c = ub;
	um = t_um;
	fpc = t_fpc;

	retVal += theMaterials[0]->commitState();

	// call the base class method
	retVal += this->Element::commitState();
	if (retVal < 0) {
		opserr << "TSM_2D::commitState() - failed in base class\n";
		return retVal;
	}
	return retVal;
}

int TSM_2D::revertToLastCommit()
{
	//opserr << "revertToLastCommit \n";
	int retVal = 0;
  
	retVal += theMaterials[0]->revertToLastCommit();
	
	return retVal;
}

int TSM_2D::revertToStart()
{
	//opserr << "revertToStart \n";
	int errCode = 0;

	ub.Zero();
	qb.Zero();
	ub_c.Zero();

	s = 0.0;
	theta = 0.0;
	t_s = 0.0;
	t_theta = 0.0;
	um = fc / Kvi;
	t_um = um;
	fpc = fc;
	t_fpc = fpc;

	// reset stiffness matrix in basic system
	kt = kinit;

	errCode += theMaterials[0]->revertToStart();

	return errCode;
}

int TSM_2D::update(void)
{
	//opserr << "update \n";
	int errCode = 0;
	// get global trial displacements and velocities
	const Vector& dsp1 = theNodes[0]->getTrialDisp();
	const Vector& dsp2 = theNodes[1]->getTrialDisp();
	const Vector& vel1 = theNodes[0]->getTrialVel();
	const Vector& vel2 = theNodes[1]->getTrialVel();

	static Vector ug(6), ugdot(6), uldot(6), ubdot(3);
	for (int i = 0; i < 3; i++) {
		ug(i) = dsp1(i);  ugdot(i) = vel1(i);
		ug(i + 3) = dsp2(i);  ugdot(i + 3) = vel2(i);
	}

	// transform response from the global to the local system
	ul.addMatrixVector(0.0, Tgl, ug, 1.0);
	uldot.addMatrixVector(0.0, Tgl, ugdot, 1.0);

	// transform response from the local to the basic system
	ub.addMatrixVector(0.0, Tlb, ul, 1.0);
	ubdot.addMatrixVector(0.0, Tlb, uldot, 1.0);

	// 1) get axial force and stiffness in basic x-direction

	double uh = ub(1);
	double ac = Di / Do;
	double R = sqrt(1 + pow(ac, 2)) * Do / 4;

	double Kv;
	Kv = Kvi / (1 + (fas / pow(PI, 2)) * pow(uh / R, 2));
	kt.resize(3, 3);
	qb.resize(3);
	qb.Zero();
	
	
	if (ub(0) > 0.0) {
		t_um = max(ub(0), um);
		double uc = fc / Kv;
		double phi = pm * (1 - expo(-ac * ((t_um - uc) / uc)));
		double fcn = fc * (1 - phi);
		double ucn = fcn / Kv;
		double kd = (fpc - fcn) / (t_um - ucn);

		if (ub(0) <= ucn) {
			qb(0) = Kv * ub(0);
			kt(0, 0) = Kv;
		}
		else if (ub(0) < um) {
			qb(0) = fcn + kd * (ub(0) - ucn);
			kt(0, 0) = kd;
		}
		else {
			qb(0) = fc * (1 + (1 / (k * Tr)) * (1 - expo(-k * (ub(0) - uc))));
			kt(0, 0) = (fc/k)* expo(-k * (ub(0) - uc));
			t_fpc = max(fpc, qb(0));
		}
	}
	else {
		qb(0) = ub(0) * Kv;
		kt(0, 0) = Kv;
	}

	// 2) evaluate the shear force and stiffness

	double s_p = s;
	double theta_p = theta;

	t_s = s;
	t_theta = theta;

	int iter = 0;
	double tol = 1.0e-15;
	double error = 1.0;
	int max_iter = 50;
	double fs_p = 0;

	double kshear;
	double F;

	do {
		t_s = uh - L * t_theta;
		//opserr << "s: " << t_s << "\n";
		//opserr << "v: " << ubdot(1) << "\n";
		theMaterials[0]->setTrialStrain(t_s, ubdot(1));

		double fs = theMaterials[0]->getStress();
		kshear = theMaterials[0]->getTangent();

		//opserr << "ks: " << kshear << "\n";
		F = fs + qb(0) * t_theta;

		t_theta = (F * L - qb(0) * t_s) / (Kr + qb(0) * L);

		error = pow(abs(fs_p - fs) + abs(s_p - t_s) + abs(theta_p - t_theta), 2);

		s_p = t_s;
		theta_p = t_theta;
		fs_p = fs;
		iter++;
	} while ((error > tol) && (iter < max_iter));

	qb(1) = F;

	kt(1, 1) = kshear -qb(0) / (Kr+ qb(0)*L);


	// 3) bending shear and stiffness

	qb(2) = ub(2) * Kb;
	kt(2, 2) = Kb;
	
  return 0;
}

const Matrix& TSM_2D::getTangentStiff(void)
{
	//opserr << "getTangentStiff \n";
	// zero the matrix
	theMatrix.Zero();

	// transform from basic to local system
	static Matrix kl(6, 6);
	kl.addMatrixTripleProduct(0.0, Tlb, kt, 1.0);

	// // add geometric stiffness to local stiffness
	// double kGeo1 = 0.5 * qb(0);
	// kl(2, 1) -= kGeo1;
	// kl(2, 4) += kGeo1;
	// kl(5, 1) -= kGeo1;
	// kl(5, 4) += kGeo1;
	// double kGeo2 = kGeo1 * 0.5 * L;
	// kl(2, 2) += kGeo2;
	// kl(5, 2) -= kGeo2;
	// double kGeo3 = kGeo1 * (1.0 - 0.5) * L;
	// kl(2, 5) -= kGeo3;
	// kl(5, 5) += kGeo3;

	// transform from local to global system
	static Matrix kg(6, 6);
	kg.addMatrixTripleProduct(0.0, Tgl, kl, 1.0);
	
	return kg;
}

const Matrix& TSM_2D::getInitialStiff(void)
{
	//opserr << "getInitialStiff \n";
	// zero the matrix
	theMatrix.Zero();

	// transform from basic to local system
	static Matrix klInit(6, 6);
	klInit.addMatrixTripleProduct(0.0, Tlb, kinit, 1.0);

	// transform from local to global system
	theMatrix.addMatrixTripleProduct(0.0, Tgl, klInit, 1.0);

	return theMatrix;
}
    
void TSM_2D::zeroLoad(void)
{
	//opserr << "zeroLoad \n";
	theLoad.Zero();
}

int TSM_2D::addLoad(const Vector &addP)
{
	// opserr << "addLoad \n";
	opserr << "TSM_2D::addLoad() - "
		<< "load type unknown for element: "
		<< this->getTag() << ".\n";

	return -1;
}

int TSM_2D::addInertiaLoadToUnbalance(const Vector &accel)
{
	//opserr << "addInertiaLoadToUnbalance \n";
	// check for quick return
	if (mass == 0.0)
		return 0;

	// get R * accel from the nodes
	const Vector& Raccel1 = theNodes[0]->getRV(accel);
	const Vector& Raccel2 = theNodes[1]->getRV(accel);

	if (3 != Raccel1.Size() || 3 != Raccel2.Size()) {
		opserr << "TSM_2D::addInertiaLoadToUnbalance() - "
			<< "matrix and vector sizes are incompatible.\n";
		return -1;
	}

	// want to add ( - fact * M R * accel ) to unbalance
	// take advantage of lumped mass matrix
	double m = 0.5 * mass;
	for (int i = 0; i < 2; i++) {
		theLoad(i) -= m * Raccel1(i);
		theLoad(i + 3) -= m * Raccel2(i);
	}

	return 0;
}

const Vector& TSM_2D::getResistingForce()
{	
	//opserr << "getResistingForce\n";
	// zero the residual
	theVector.Zero();
	theVector.resize(6);
	// determine resisting forces in local system
	static Vector ql(6);
	ql.addMatrixTransposeVector(0.0, Tlb, qb, 1.0);

	// add P-Delta moments to local forces
	// double kGeo1 = 0.5 * qb(0);
	// double MpDelta1 = kGeo1 * (ul(4) - ul(1));
	// ql(2) += MpDelta1;
	// ql(5) += MpDelta1;
	// double MpDelta2 = kGeo1 * 0 * L * ul(2);
	// ql(2) += MpDelta2;
	// ql(5) -= MpDelta2;
	// double MpDelta3 = kGeo1 * (1.0 - 0) * L * ul(5);
	// ql(2) -= MpDelta3;
	// ql(5) += MpDelta3;

	// determine resisting forces in global system
	theVector.addMatrixTransposeVector(0.0, Tgl, ql, 1.0);

  return theVector;
}

const Vector& TSM_2D::getResistingForceIncInertia()
{
	//opserr << "getResistingForceIncInertia \n";
	// this already includes damping forces from materials
	theVector = this->getResistingForce();

	// subtract external load
	theVector.addVector(1.0, theLoad, -1.0);

	// add inertia forces from element mass
	if (mass != 0.0) {
		const Vector& accel1 = theNodes[0]->getTrialAccel();
		const Vector& accel2 = theNodes[1]->getTrialAccel();

		double m = 0.5 * mass;
		for (int i = 0; i < 2; i++) {
			theVector(i) += m * accel1(i);
			theVector(i + 3) += m * accel2(i);
		}
	}

	return theVector;
}

int TSM_2D::sendSelf(int commitTag, Channel &theChannel)
{
	//opserr << " sendSelf \n";
	static Vector data(15);
	data(0) = this->getTag();
	data(1) = Kvi;
	data(2) = fas;
	data(3) = Kb;
	data(4) = Kr;
	data(5) = Do;
	data(6) = Di;
	data(7) = Tr;
	data(8) = fc;
	data(9) = a;
	data(10) = k;
	data(11) = pm;
	data(12) = mass;
	data(13) = x.Size();
	data(14) = y.Size();
	theChannel.sendVector(0, commitTag, data);

	// send the two end nodes
	theChannel.sendID(0, commitTag, connectedExternalNodes);

	// send the material class tags
	ID matClassTags(1);
	matClassTags(0) = theMaterials[0]->getClassTag();
	theChannel.sendID(0, commitTag, matClassTags);


	theMaterials[0]->sendSelf(commitTag, theChannel);

	// send remaining data
	if (x.Size() == 3)
		theChannel.sendVector(0, commitTag, x);
	if (y.Size() == 3)
		theChannel.sendVector(0, commitTag, y);

  return 0;
}

int TSM_2D::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	//opserr << "recvSelf \n";
	// delete material memory
	if (theMaterials[0] != 0)
		delete theMaterials[0];

	// receive element parameters
	static Vector data(10);
	theChannel.recvVector(0, commitTag, data);
	this->setTag((int)data(0));
	Kvi = data(1);
	fas = data(2);
	Kb = data(3);
	Kr = data(4);
	Do = data(5);
	Di = data(6);
	Tr = data(7);
	fc = data(8);
	a = data(9);
	k = data(10);
	pm = data(11);
	mass=data(12) = mass;

	// receive the two end nodes
	theChannel.recvID(0, commitTag, connectedExternalNodes);

	// receive the material class tags
	ID matClassTags(1);
	theChannel.recvID(0, commitTag, matClassTags);


	theMaterials[0] = theBroker.getNewUniaxialMaterial(matClassTags(0));
	if (theMaterials[0] == 0) {
		opserr << "TSM_2D::recvSelf() - "
			<< "failed to get blank uniaxial material.\n";
		return -2;
	}
	theMaterials[0]->recvSelf(commitTag, theChannel, theBroker);


	// receive remaining data
	if ((int)data(13) == 3) {
		x.resize(3);
		theChannel.recvVector(0, commitTag, x);
	}
	if ((int)data(14) == 3) {
		y.resize(3);
		theChannel.recvVector(0, commitTag, y);
	}
	

	// initialize initial stiffness matrix
	kinit.Zero();
	kinit.resize(3, 3);
	kinit(0, 0) = Kvi;
	kinit(1, 1) = theMaterials[0]->getInitialTangent();
	kinit(2, 2) = Kb;

	// initialize other variables
	this->revertToStart();

	return 0;
}


int TSM_2D::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
	//opserr << "displaySelf \n";
	// first determine the end points of the element based on
	// the display factor (a measure of the distorted image)
	const Vector& end1Crd = theNodes[0]->getCrds();
	const Vector& end2Crd = theNodes[1]->getCrds();

	static Vector v1(3);
	static Vector v2(3);

	if (displayMode >= 0) {
		const Vector& end1Disp = theNodes[0]->getDisp();
		const Vector& end2Disp = theNodes[1]->getDisp();

		for (int i = 0; i < 2; i++) {
			v1(i) = end1Crd(i) + end1Disp(i) * fact;
			v2(i) = end2Crd(i) + end2Disp(i) * fact;
		}
	}
	else {
		int mode = displayMode * -1;
		const Matrix& eigen1 = theNodes[0]->getEigenvectors();
		const Matrix& eigen2 = theNodes[1]->getEigenvectors();

		if (eigen1.noCols() >= mode) {
			for (int i = 0; i < 2; i++) {
				v1(i) = end1Crd(i) + eigen1(i, mode - 1) * fact;
				v2(i) = end2Crd(i) + eigen2(i, mode - 1) * fact;
			}
		}
		else {
			for (int i = 0; i < 2; i++) {
				v1(i) = end1Crd(i);
				v2(i) = end2Crd(i);
			}
		}
	}

	return theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag(), 0);
}

void TSM_2D::Print(OPS_Stream &s, int flag)
{
	//opserr << "Print \n";
	if (flag == OPS_PRINT_CURRENTSTATE) {
		// print everything
		s << "Element: " << this->getTag() << endln;
		s << "  type: TSM_2D\n";
		s << "  iNode: " << connectedExternalNodes(0);
		s << "  jNode: " << connectedExternalNodes(1) << endln;
		s << "  Shear Material: " << theMaterials[0]->getTag() << endln;
		s << "  mass: " << mass << endln;
		// determine resisting forces in global system
		s << "  resisting force: " << this->getResistingForce() << endln;
	}

	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": " << this->getTag() << ", ";
		s << "\"type\": \"TSM_2D\", ";
		s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
		s << "\"Shear material\": [\"";
		s << theMaterials[0]->getTag() << "\"], ";
		s << "\"mass\": " << mass << "}";
	}
}

Response* TSM_2D::setResponse(const char **argv, int argc, OPS_Stream &output)
{
	//opserr << "setRespone\n";
	Response* theResponse = 0;

	output.tag("ElementOutput");
	output.attr("eleType", "ElastomericBearingBoucWen2d");
	output.attr("eleTag", this->getTag());
	output.attr("node1", connectedExternalNodes[0]);
	output.attr("node2", connectedExternalNodes[1]);

	// global forces
	if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0 ||
		strcmp(argv[0], "globalForce") == 0 || strcmp(argv[0], "globalForces") == 0)
	{
		output.tag("ResponseType", "Px_1");
		output.tag("ResponseType", "Py_1");
		output.tag("ResponseType", "Mz_1");
		output.tag("ResponseType", "Px_2");
		output.tag("ResponseType", "Py_2");
		output.tag("ResponseType", "Mz_2");

		theResponse = new ElementResponse(this, 1, theVector);
	}
	// local forces
	else if (strcmp(argv[0], "localForce") == 0 || strcmp(argv[0], "localForces") == 0)
	{
		output.tag("ResponseType", "N_1");
		output.tag("ResponseType", "V_1");
		output.tag("ResponseType", "M_1");
		output.tag("ResponseType", "N_2");
		output.tag("ResponseType", "V_2");
		output.tag("ResponseType", "M_2");

		theResponse = new ElementResponse(this, 2, theVector);
	}
	// basic forces
	else if (strcmp(argv[0], "basicForce") == 0 || strcmp(argv[0], "basicForces") == 0)
	{
		output.tag("ResponseType", "qb1");
		output.tag("ResponseType", "qb2");
		output.tag("ResponseType", "qb3");

		theResponse = new ElementResponse(this, 3, Vector(3));
	}
	// local displacements
	else if (strcmp(argv[0], "localDisplacement") == 0 ||
		strcmp(argv[0], "localDisplacements") == 0)
	{
		output.tag("ResponseType", "ux_1");
		output.tag("ResponseType", "uy_1");
		output.tag("ResponseType", "rz_1");
		output.tag("ResponseType", "ux_2");
		output.tag("ResponseType", "uy_2");
		output.tag("ResponseType", "rz_2");

		theResponse = new ElementResponse(this, 4, theVector);
	}
	// basic displacements
	else if (strcmp(argv[0], "deformation") == 0 || strcmp(argv[0], "deformations") == 0 ||
		strcmp(argv[0], "basicDeformation") == 0 || strcmp(argv[0], "basicDeformations") == 0 ||
		strcmp(argv[0], "basicDisplacement") == 0 || strcmp(argv[0], "basicDisplacements") == 0)
	{
		output.tag("ResponseType", "ub1");
		output.tag("ResponseType", "ub2");
		output.tag("ResponseType", "ub3");

		theResponse = new ElementResponse(this, 5, Vector(3));
	}
	
	// material output
	else if (strcmp(argv[0], "material") == 0) {
		if (argc > 2) {
			int matNum = atoi(argv[1]);
			if (matNum == 1)
				theResponse = theMaterials[0]->setResponse(&argv[2], argc - 2, output);
		}
	}

	output.endTag(); // ElementOutput

	return theResponse;
}

int TSM_2D::getResponse(int responseID, Information &eleInfo)
{
	//opserr << "getRespone\n";
	// double kGeo1, MpDelta1, MpDelta2, MpDelta3;

	switch (responseID) {
	case 1:  // global forces
		return eleInfo.setVector(this->getResistingForce());

	case 2:  // local forces
		theVector.Zero();
		// determine resisting forces in local system
		theVector.addMatrixTransposeVector(0.0, Tlb, qb, 1.0);
		// add P-Delta moments
		// kGeo1 = 0.5 * qb(0);
		// MpDelta1 = kGeo1 * (ul(4) - ul(1));
		// theVector(2) += MpDelta1;
		// theVector(5) += MpDelta1;
		// MpDelta2 = kGeo1 * 0.5 * L * ul(2);
		// theVector(2) += MpDelta2;
		// theVector(5) -= MpDelta2;
		// MpDelta3 = kGeo1 * (1.0 - 0.5) * L * ul(5);
		// theVector(2) -= MpDelta3;
		// theVector(5) += MpDelta3;

		return eleInfo.setVector(theVector);

	case 3:  // basic forces
		return eleInfo.setVector(qb);

	case 4:  // local displacements
		return eleInfo.setVector(ul);

	case 5:  // basic displacements
		return eleInfo.setVector(ub);

	default:
		return -1;
	}
}

int TSM_2D::setParameter(const char **argv, int argc, Parameter &param)
{
	//opserr << "setParameter\n";
  return 1;
}   

int TSM_2D::updateParameter(int parameterID, Information &info)
{
	// opserr << "updateParameter\n";
  return -1;
}

void TSM_2D::setUp()
{
	//opserr << "setUp\n";
	const Vector& end1Crd = theNodes[0]->getCrds();
	const Vector& end2Crd = theNodes[1]->getCrds();
	Vector xp = end2Crd - end1Crd;

	L = xp.Norm();

	if (L > DBL_EPSILON) {
		if (x.Size() == 0) {
			x.resize(3);
			x(0) = xp(0);  x(1) = xp(1);  x(2) = 0.0;
			y.resize(3);
			y(0) = -x(1);  y(1) = x(0);  y(2) = 0.0;
		}
		else  {
			opserr << "Error TSM_2D::setUp() - "
				<< "element: " << this->getTag()
				<< " - incorrct length of the element.\n";
			exit(-1);
		}
	}
	// check that vectors for orientation are of correct size
	if (x.Size() != 3 || y.Size() != 3) {
		opserr << "TSM_2D::setUp() - "
			<< "element: " << this->getTag()
			<< " - incorrect dimension of orientation vectors.\n";
		exit(-1);
	}

	// establish orientation of element for the transformation matrix
	// z = x cross y
	static Vector z(3);
	z(0) = x(1) * y(2) - x(2) * y(1);
	z(1) = x(2) * y(0) - x(0) * y(2);
	z(2) = x(0) * y(1) - x(1) * y(0);

	// y = z cross x
	y(0) = z(1) * x(2) - z(2) * x(1);
	y(1) = z(2) * x(0) - z(0) * x(2);
	y(2) = z(0) * x(1) - z(1) * x(0);

	// compute length(norm) of vectors
	double xn = x.Norm();
	double yn = y.Norm();
	double zn = z.Norm();

	// check valid x and y vectors, i.e. not parallel and of zero length
	if (xn == 0 || yn == 0 || zn == 0) {
		opserr << "TSM_2D::setUp() - "
			<< "element: " << this->getTag()
			<< " - invalid orientation vectors.\n";
		exit(-1);
	}

	// create transformation matrix from global to local system
	Tgl.Zero();
	Tgl(0, 0) = Tgl(3, 3) = x(0) / xn;
	Tgl(0, 1) = Tgl(3, 4) = x(1) / xn;
	Tgl(1, 0) = Tgl(4, 3) = y(0) / yn;
	Tgl(1, 1) = Tgl(4, 4) = y(1) / yn;
	Tgl(2, 2) = Tgl(5, 5) = z(2) / zn;

	// create transformation matrix from local to basic system (linear)
	Tlb.Zero();
	Tlb(0, 0) = Tlb(1, 1) = Tlb(2, 2) = -1.0;
	Tlb(0, 3) = Tlb(1, 4) = Tlb(2, 5) = 1.0;
	Tlb(1, 2) = - L;
	Tlb(1, 5) = 0.0;

}

double TSM_2D::abs(double x)
{
	if (x < 0) return -x;
	else return x;
}

double TSM_2D::max(double x,double y) {
	if (x < y) return y;
	else return x;
}

double TSM_2D::expo(double x) {
	double e = 2.718281828459046;
	return pow(e, x);
}

const Matrix& TSM_2D::getDamp()
{
	// opserr << "getDamp \n";
	// zero the matrix
	theMatrix.Zero();

	// call base class to setup Rayleigh damping
	double factThis = 0.0;


	// now add damping tangent from materials
	static Matrix cb(3, 3);
	cb.Zero();
	cb(0, 0) = 0.0;
	cb(1, 1) = 0.0;
	cb(2, 2) = 0.0;

	// transform from basic to local system
	static Matrix cl(6, 6);
	cl.addMatrixTripleProduct(0.0, Tlb, cb, 1.0);

	// transform from local to global system and add to cg
	theMatrix.addMatrixTripleProduct(factThis, Tgl, cl, 1.0);

	return theMatrix;
}

const Matrix& TSM_2D::getMass()
{
	// zero the matrix
	theMatrix.Zero();

	// check for quick return
	if (mass == 0.0) {
		return theMatrix;
	}

	double m = 0.5 * mass;
	for (int i = 0; i < 2; i++) {
		theMatrix(i, i) = m;
		theMatrix(i + 3, i + 3) = m;
	}

	return theMatrix;
}