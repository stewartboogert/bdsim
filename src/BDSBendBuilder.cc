#include "globals.hh" // geant4 globals / types
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"

#include "BDSDebug.hh"
#include "BDSAcceleratorComponent.hh"
#include "BDSBendBuilder.hh"
#include "BDSComponentFactory.hh"
#include "BDSMagnet.hh"
#include "BDSSimpleComponent.hh"
#include "BDSElement.hh"
#include "BDSBeamPipeInfo.hh"
#include "BDSMagnetOuterInfo.hh"
#include "BDSTiltOffset.hh"
#include "BDSTransform3D.hh"
#include "BDSUtilities.hh"
#include "parser/element.h"
#include "parser/elementtype.h"


#include <algorithm>
#include <iterator>
#include <ostream>
#include <utility>  // for std::pair
#include <vector>

using namespace GMAD;

BDSBendBuilder* BDSBendBuilder::_instance = nullptr;

BDSBendBuilder* BDSBendBuilder::Instance()
{
    if(_instance == nullptr)
    {_instance = new BDSBendBuilder();}
    return _instance;
}


BDSBendBuilder::BDSBendBuilder()
{
  outerDiameter = BDSGlobalConstants::Instance()->OuterDiameter();
}

BDSBendBuilder::~BDSBendBuilder()
{
    _instance = nullptr;
}

BDSLine* BDSBendBuilder::SBendLine(Element*  element,
                                         G4double angleIn,
                                         G4double angleOut,
                                         BDSMagnetStrength* st)
{

    // Calculate number of sbends to split parent into
    G4int nSBends = CalculateNSBendSegments(element);

    //Zero angle bend only needs one element.
    std::string thename = element->name + "_1_of_1";

    G4double length = element->l * CLHEP::m;
    //copy of angle
    G4double angle = (*st)["angle"];

    BDSLine* sbendline  = new BDSLine(element->name);

    // Single element if no poleface and zero bend angle or dontSplitSBends=1, therefore nSBends = 1
    if (!BDS::IsFinite(angle) || (nSBends == 1))
    {
        BDSFieldInfo* vacuumField = new BDSFieldInfo(BDSFieldType::dipole,
                                                     brho,
                                                     BDSIntegratorType::dipole,
                                                     st);
        // prepare one sbend segment
        BDSMagnet* oneBend = new BDSMagnet(BDSMagnetType::sectorbend,
                                           thename,
                                           length,
                                           PrepareBeamPipeInfo(element, -angleIn, -angleOut),
                                           PrepareMagnetOuterInfo(element, -angleIn, -angleOut),
                                           vacuumField);

        oneBend->SetBiasVacuumList(element->biasVacuumList);
        oneBend->SetBiasMaterialList(element->biasMaterialList);

        sbendline->AddComponent(oneBend);

        return sbendline;
    }

    // prepare one name for all that makes sense
    thename = element->name + "_1_of_" + std::to_string(nSBends);
    //calculate their angles and length
    G4double semiangle  = -element->angle / (G4double) nSBends;
    G4double semilength = length / (G4double) nSBends;
    G4double rho        = element->l*CLHEP::m/element->angle;

    angle      = element->angle;
    angleIn    = element->e1*CLHEP::rad;
    angleOut   = element->e2*CLHEP::rad;

    BDSMagnetType magType = BDSMagnetType::sectorbend;
    // check magnet outer info
    BDSMagnetOuterInfo* magnetOuterInfoCheck = PrepareMagnetOuterInfo(element,angleIn,angleOut);
    CheckBendLengthAngleWidthCombo(semilength, semiangle, magnetOuterInfoCheck->outerDiameter, thename);
    // clean up
    delete magnetOuterInfoCheck;

    // angle increment for sbend elements with poleface rotation(s) specified
    G4double deltastart = -element->e1/(0.5*(nSBends-1));
    G4double deltaend   = -element->e2/(0.5*(nSBends-1));

    // first element should be fringe if poleface specified
    if (BDS::IsFinite(angleIn) && includeFringe)
    {
        BDSMagnetStrength* fringeStIn  = new BDSMagnetStrength();
        (*fringeStIn)["field"]         = (*st)["field"];
        (*fringeStIn)["length"]        = thinElementLength;
        (*fringeStIn)["angle"]         = -thinElementLength/rho;
        (*fringeStIn)["polefaceangle"] = element->e1;
        thename                        = element->name + "_e1_fringe";
        angle                          = element->e1 + 0.5*((*fringeStIn)["angle"]);
        BDSMagnet* startfringe = DipoleFringe(element, -angle, angle, thename, magType, fringeStIn);
        sbendline->AddComponent(startfringe);
    }

    for (int i = 0; i < nSBends; ++i)
    {
        thename = element->name + "_"+std::to_string(i+1)+"_of_" + std::to_string(nSBends);

        // subtract thinElementLength from first and last elements if fringe & poleface specified
        length = semilength;
        if ((BDS::IsFinite(element->e1)) && (i == 0) && includeFringe)
        {length -= thinElementLength;}
        if ((BDS::IsFinite(element->e2)) && (i == nSBends-1) && includeFringe)
        {length -= thinElementLength;}
        semiangle = -length/rho;

        // Default angles for all segments
        angleIn = -semiangle*0.5;
        angleOut = -semiangle*0.5;

        // Input and output angles added to or subtracted from the default as appropriate
        // Note: case of i == 0.5*(nSBends-1) is just the default central wedge.
        // More detailed methodology/reasons in developer manual
        if ((BDS::IsFinite(element->e1))||(BDS::IsFinite(element->e2)))
        {
            if (i < 0.5*(nSBends-1))
            {
                angleIn -= (element->e1 + (i*deltastart));
                angleOut -= ((0.5*(nSBends-3)-i)*deltastart);
            }
            else if (i > 0.5*(nSBends-1))
            {
                angleIn  +=  (0.5*(nSBends+1)-i)*deltaend;
                angleOut += -(0.5*(nSBends-1)-i)*deltaend;
            }
        }
        if ((BDS::IsFinite(element->e1)) && (i == 0) && includeFringe)
        {angleIn += thinElementLength/rho;}
        if ((BDS::IsFinite(element->e2)) && (i == nSBends-1) && includeFringe)
        {angleOut += thinElementLength/rho;}

        // Check for intersection of angled faces.
        G4double intersectionX = BDS::CalculateFacesOverlapRadius(angleIn,angleOut,semilength);
        BDSMagnetOuterInfo* magnetOuterInfo = PrepareMagnetOuterInfo(element,angleIn,angleOut);
        G4double magnetRadius= 0.625*magnetOuterInfo->outerDiameter*CLHEP::mm;
        // Every geometry type has a completely arbitrary factor of 1.25 except cylindrical
        if (magnetOuterInfo->geometryType == BDSMagnetGeometryType::cylindrical)
        {magnetRadius= 0.5*magnetOuterInfo->outerDiameter*CLHEP::mm;}

        //Check if intersection is within radius
        if ((BDS::IsFinite(intersectionX)) && (std::abs(intersectionX) < magnetRadius))
        {
            G4cerr << __METHOD_NAME__ << "Angled faces of element "<< thename
            << " intersect within the magnet radius." << G4endl;
            exit(1);
        }

        BDSMagnetStrength* stSemi = new BDSMagnetStrength(*st); // copy field strength - ie B
        (*stSemi)["length"] = length;
        (*stSemi)["angle"]  = semiangle;  // override copied length and angle

        BDSFieldInfo* vacuumField = new BDSFieldInfo(BDSFieldType::dipole,
                                                     brho,
                                                     BDSIntegratorType::dipole,
                                                     stSemi);

        BDSMagnet* oneBend = new BDSMagnet(magType,
                                           thename,
                                           length,
                                           PrepareBeamPipeInfo(element, angleIn, angleOut),
                                           magnetOuterInfo,
                                           vacuumField,
                                           semiangle);

        oneBend->SetBiasVacuumList(element->biasVacuumList);
        oneBend->SetBiasMaterialList(element->biasMaterialList);
        sbendline->AddComponent(oneBend);

#ifdef BDSDEBUG
        G4cout << "---->creating sbend line,"
     << " element= " << thename
	 << " angleIn= " << angleIn
	 << " angleOut= " << angleOut << "m"
	 << G4endl;
#endif
    }
    //Last element should be fringe if poleface specified
    if (BDS::IsFinite(element->e2) && includeFringe)
    {
        BDSMagnetStrength* fringeStOut  = new BDSMagnetStrength();
        (*fringeStOut)["angle"]         = -thinElementLength/rho;
        (*fringeStOut)["field"]         = (*st)["field"];
        (*fringeStOut)["polefaceangle"] = element->e2;
        (*fringeStOut)["length"]        = thinElementLength;
        angle                           = element->e2+ 0.5*((*fringeStOut)["angle"]);
        thename                         = element->name + "_e2_fringe";

        BDSMagnet* endfringe = DipoleFringe(element, angle, -angle, thename, magType, fringeStOut);
        sbendline->AddComponent(endfringe);
    }
    return sbendline;
}



BDSLine* BDSBendBuilder::RBendLine(Element* element,
                                   Element* prevElement,
                                   Element* nextElement,
                                   G4double angleIn,
                                   G4double angleOut,
                                   G4double brho,
                                   BDSMagnetStrength* st)
{
    BDSLine* rbendline  = new BDSLine(element->name);

    PoleFaceRotationsNotTooLarge(element);

    G4double angle       = element->angle;
    G4double length      = element->l*CLHEP::m;
    G4String thename     = element->name;
    G4double rho         = element->l*CLHEP::m/element->angle;
    G4bool prevModifies  = false;
    G4bool nextModifies  = false;

    BDSMagnetType magType = BDSMagnetType::rectangularbend;

    // booleans for modification by previous/next element
    if ((prevElement) && (prevElement->type == ElementType ::_RBEND))
    {prevModifies = true;}
    if ((nextElement) && (nextElement->type == ElementType ::_RBEND))
    {nextModifies = true;}

    // poleface angles
    G4double polefaceAngleIn = element->e1 + 0.5*(length-thinElementLength)/rho;
    G4double polefaceAngleOut = element->e2 + 0.5*(length-thinElementLength)/rho;

    // poleface angles and main element angles are modified if next/previous is an rbend
    if ((prevElement) && (prevElement->type == ElementType::_RBEND)){
        polefaceAngleIn -= 0.5*element->angle;
        angleIn += 0.5*(thinElementLength)/rho;}
    if ((nextElement) && (nextElement->type == ElementType::_RBEND)){
        polefaceAngleOut -= 0.5*element->angle;
        angleOut += 0.5*(thinElementLength)/rho;}

    // first element should be fringe if poleface specified
    if (BDS::IsFinite(element->e1) && includeFringe &&(!prevModifies))
    {
        BDSMagnetStrength* fringeStIn  = new BDSMagnetStrength();
        (*fringeStIn)["field"]         = (*st)["field"];
        (*fringeStIn)["polefaceangle"] = element->e1;
        (*fringeStIn)["length"]        = thinElementLength;
        (*fringeStIn)["angle"]         = -thinElementLength/rho;
        thename                        = element->name + "_e1_fringe";
        angle                          = polefaceAngleIn;

        BDSMagnet* startfringe = DipoleFringe(element, -angle, angle, thename, magType, fringeStIn);
        rbendline->AddComponent(startfringe);
    }

    // subtract thinElementLength from main rbend element if fringe & poleface(s) specified
    if (BDS::IsFinite(element->e1) && includeFringe && (!prevModifies))
    {length   -= thinElementLength;
     angleIn  += 0.5*(thinElementLength)/rho;
     angleOut -= 0.5*(thinElementLength)/rho;}
    if (BDS::IsFinite(element->e2) && includeFringe && (!nextModifies))
    {length   -= thinElementLength;
     angleOut += 0.5*(thinElementLength)/rho;
     angleIn  -= 0.5*(thinElementLength)/rho;}
    angle = -length/rho;

    if (nextModifies){
        angleOut  -= 0.5*(thinElementLength)/rho;
    }
    if (prevModifies){
        angleIn  -= 0.5*(thinElementLength)/rho;
    }

    // override copied length and angle
    (*st)["length"] = length;
    (*st)["angle"]  = angle;

    BDSFieldInfo* vacuumField = new BDSFieldInfo(BDSFieldType::dipole,
                                                 brho,
                                                 BDSIntegratorType::dipole,
                                                 st);

    BDSMagnet* oneBend = new BDSMagnet(magType,
                                       element->name,
                                       length,
                                       PrepareBeamPipeInfo(element, angleIn, angleOut),
                                       PrepareMagnetOuterInfo(element, angleIn, angleOut),
                                       vacuumField,
                                       angle,
                                       nullptr);

    rbendline->AddComponent(oneBend);

    //Last element should be fringe if poleface specified
    if (BDS::IsFinite(element->e2) && includeFringe && (!nextModifies)){
        BDSMagnetStrength *fringeStOut = new BDSMagnetStrength();
        (*fringeStOut)["field"] = (*st)["field"];
        (*fringeStOut)["polefaceangle"] = element->e2;
        (*fringeStOut)["length"] = thinElementLength;
        (*fringeStOut)["angle"] = -thinElementLength / rho;
        thename = element->name + "_e2_fringe";
        angle = polefaceAngleOut;

        BDSMagnet *endfringe = DipoleFringe(element, angle, -angle, thename, magType, fringeStOut);
        rbendline->AddComponent(endfringe);
    }

    return rbendline;
}



BDSMagnet* BDSBendBuilder::DipoleFringe(GMAD::Element* element,
                                                G4double angleIn,
                                                G4double angleOut,
                                                G4String name,
                                                BDSMagnetType magType,
                                                BDSMagnetStrength* st)
{
    BDSBeamPipeInfo* beamPipeInfo = PrepareBeamPipeInfo(element, angleIn, angleOut);
    BDSMagnetOuterInfo* magnetOuterInfo = PrepareMagnetOuterInfo(element, angleIn, angleOut);
    magnetOuterInfo->geometryType = BDSMagnetGeometryType::none;

    BDSFieldInfo* vacuumField = new BDSFieldInfo(BDSFieldType::dipole,
                                                 brho,
                                                 BDSIntegratorType::fringe,
                                                 st);

    return new BDSMagnet(magType,
                         name,
                         (*st)["length"],
                         beamPipeInfo,
                         magnetOuterInfo,
                         vacuumField,
                         (*st)["angle"],
                         nullptr);
}


G4int BDSBendBuilder::CalculateNSBendSegments(GMAD::Element const* element,
                              const G4double aperturePrecision)
{
    //if maximum distance between arc path and straight path larger than 1mm, split sbend into N chunks,
    //this also works when maximum distance is less than 1mm as there will just be 1 chunk!

    G4double length = element->l*CLHEP::m;
    // from formula: L/2 / N tan (angle/N) < precision. (L=physical length)
    // add poleface rotations onto angle as absolute number (just to be safe)
    G4double totalAngle = std::abs(element->angle) + std::abs(element->e1) + std::abs(element->e2);
    G4int nSBends = (G4int) ceil(std::sqrt(length*totalAngle/2/aperturePrecision));
    if (nSBends==0)
    {nSBends = 1;} // can happen in case angle = 0
    if (BDSGlobalConstants::Instance()->DontSplitSBends())
    {nSBends = 1;}  // use for debugging
    if (nSBends % 2 == 0)
    {nSBends += 1;} // always have odd number of poles for poleface rotations
#ifdef BDSDEBUG
    G4cout << __METHOD_NAME__ << " splitting sbend into " << nSBends << " sbends" << G4endl;
#endif
    return nSBends;
}