/* 
Beam Delivery Simulation (BDSIM) Copyright (C) BDSIM Collaboration, 2001 - 2026.

This file is part of BDSIM.

BDSIM is free software: you can redistribute it and/or modify 
it under the terms of the GNU General Public License as published 
by the Free Software Foundation version 3 of the License.

BDSIM is distributed in the hope that it will be useful, but 
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BDSIM.  If not, see <http://www.gnu.org/licenses/>.
*/
//
// Created by Stewart Boogert on 07/06/2025.
//
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"

#include "BDSLinkTrackerInterface.hh"
#include "BDSIMLink.hh"
#include "BDSLinkBunch.hh"
#include "BDSSamplerCustom.hh"
#include "BDSParser.hh"
#include "BDSException.hh"
#include "BDSIonDefinition.hh"
#include "BDSParticleDefinition.hh"
#include "BDSAcceleratorComponent.hh"
#include "BDSAcceleratorComponentRegistry.hh"
#include "BDSLinkComponent.hh"
#include "BDSLinkOpaqueBox.hh"

#include "G4Electron.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4Transform3D.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"

#include "G4FastSimulationHelper.hh"
#include "G4FastSimulationManagerProcess.hh"
#include "G4Proton.hh"
#include "G4coutDestination.hh"
#include "G4ios.hh"
#include "G4GenericIon.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "CLHEP/Units/SystemOfUnits.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>


template <typename T>
T* make_ptr(py::array_t<T> &arr) {
  auto buf = arr.request();
  T* ptr = static_cast<T*>(buf.ptr);
  return ptr;
}

template <typename T>
inline void set_element(py::array_t<T> &arr, int index, T value) {
  T* ptr = make_ptr<T>(arr);
  ptr[index] = value;
}

void TrackXSuite(BDSLinkTrackerInterface *tracker_interface,
                 int iElement,
                 std::string elementName,
                 py::object particles,
                 int verbosity = 0);

void TrackRFTrack(BDSLinkTrackerInterface *tracker_interface,
                  py::object particles);


namespace {
// Geant4 produces excited nuclei (isomers), with a non-zero level I in their PDG
// ID (10LZZZAAAI), but G4IonTable::GetIon() does not find most of these when they
// are sent back into the link in a later element ("Ion ... not found"). As a
// workaround, ions are handled as their ground state (I = 0): when they are
// returned to xtrack, and when an excited ion that cannot be found is sent in.
// The mass and charge returned by BDSIM are not changed (the excitation energy is
// at most a few MeV).
// TODO: fix this in BDSIM, such that excited ions (isomers) can be sent into the
// link as well, and remove this workaround.
inline int64_t GroundStatePDGID(int64_t pdgID)
{
  if (pdgID > 1000000000)
    {return pdgID - pdgID % 10;}
  return pdgID;
}

// Particle definition for a particle coming from xsuite, given its momentum.
// BDSLinkTrackerInterface::PrepareBDSParticleDefinition is not used: it does not
// include the electron mass for partially stripped ions, has no check that the
// particle exists, and its declaration and definition have the momentum and
// kinetic energy arguments in a different order.
BDSParticleDefinition* PrepareBDSParticleDefinition(long long int pdgIDIn, double momentumIn,
                                                    double kineticEnergyIn, double ionChargeIn)
{
  G4int pdgID = (G4int) pdgIDIn;

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particleDefGeant;

  // Wrap in our class that calculates momentum and kinetic energy.
  // Requires that one of E, Ek, P be non-zero (only one).
  BDSParticleDefinition* particleDefinition = nullptr;
  BDSIonDefinition* ionDef = nullptr;

  // PDG for ions = 10LZZZAAAI
  if (pdgID > 1000000000) // is an ion
  {

    G4IonTable* ionTable = particleTable->GetIonTable();
    particleDefGeant = ionTable->GetIon(pdgID);
    if (!particleDefGeant && pdgID % 10 != 0)
    {
      // Excited ion (isomer) that BDSIM cannot find: use the ground state instead
      // (see GroundStatePDGID; TODO: fix in BDSIM)
      pdgID = (G4int) GroundStatePDGID(pdgID);
      particleDefGeant = ionTable->GetIon(pdgID);
    }
    if (!particleDefGeant)
    {throw BDSException("BDSXtrackInterface> Ion \"" + std::to_string(pdgID) + "\" not found");}

    G4int ionCharge = ionChargeIn==0 ? particleDefGeant->GetAtomicNumber() : (G4int) ionChargeIn;
    ionDef = new BDSIonDefinition(particleDefGeant->GetAtomicMass(),
                                  particleDefGeant->GetAtomicNumber(),
                                  ionCharge);

    // correct the particle mass in the case of partially-stripped ions
    G4double mass   = ionTable->GetIonMass(ionDef->Z(), ionDef->A());
    G4double charge = ionDef->Charge(); // correct even if overridden
    mass += ionDef->NElectrons()*G4Electron::Definition()->GetPDGMass();
    // The constructor with a custom name requires the bdsim name of the ion
    std::string bdsimPartName = "ion " + std::to_string(static_cast<int>(ionDef->A()))
                                + " " + std::to_string(static_cast<int>(ionDef->Z()))
                                + " " + std::to_string(static_cast<int>(charge));

    particleDefinition = new BDSParticleDefinition(bdsimPartName, mass, charge, 0,
                                                   kineticEnergyIn, momentumIn, 1, ionDef, pdgID);
  }
  else
  {
    particleDefGeant = particleTable->FindParticle(pdgID);
    if (!particleDefGeant)
    {throw BDSException("BDSXtrackInterface> Particle \"" + std::to_string(pdgID) + "\" not found");}
    particleDefinition = new BDSParticleDefinition(particleDefGeant, 0, kineticEnergyIn, momentumIn, 1, nullptr);
  }

  return particleDefinition;
}
} // namespace

PYBIND11_MODULE(bdslinktrackerinterface, m) {
  py::class_<BDSLinkTrackerInterface>(m,"BDSLinkTrackerInterface")
      .def_static("GetInstance", [](std::string bdsimConfigFile,
                                    int referenceParticlePDG,
                                    double referenceKineticEnergy,
                                    double relativeEnergyCut,
                                    int seed,
                                    int referenceIonCharge,
                                    bool batchMode,
                                    bool no_neutral_particles) {
        auto* obj = BDSLinkTrackerInterface::GetInstance(bdsimConfigFile,
                                                         referenceParticlePDG,
                                                         referenceKineticEnergy,
                                                         relativeEnergyCut,
                                                         seed,
                                                         referenceIonCharge,
                                                         batchMode);
        // Reference particle with the electron mass for partially stripped ions
        obj->SetReferenceParticleDefinition(PrepareBDSParticleDefinition(referenceParticlePDG,
                                                                         0,
                                                                         referenceKineticEnergy,
                                                                         referenceIonCharge));
        obj->SetNoNeutralParticles(no_neutral_particles);
        return obj;
    },            py::arg("bdsimConfigFileIn") = "trackerInterface.gmad",
                  py::arg("referenceParticlePDG") =11,
                  py::arg("referenceKineticEnergy") = 100,
                  py::arg("relativeEnergyCut") = 0.01,
                  py::arg("seed") = 1234,
                  py::arg("referenceIonCharge") = 1,
                  py::arg("batchMode") = true,
                  py::arg("no_neutral_particles") = true,
                  py::return_value_policy::reference)

    .def_static("GetInstance", [](BDSParser* parser,
                                  int referenceParticlePDG,
                                  double referenceKineticEnergy,
                                  double relativeEnergyCut,
                                  int seed,
                                  int referenceIonCharge,
                                  bool batchMode,
                                  bool no_neutral_particles) {
      auto* obj = BDSLinkTrackerInterface::GetInstance(parser,
                                                       referenceParticlePDG,
                                                       referenceKineticEnergy,
                                                       relativeEnergyCut,
                                                       seed,
                                                       referenceIonCharge,
                                                       batchMode);
      obj->SetNoNeutralParticles(no_neutral_particles);
      return obj;
    }, py::arg("parser"),
       py::arg("referenceParticlePDG") =11,
       py::arg("referenceKineticEnergy") = 100,
       py::arg("relativeEnergyCut") = 0.01,
       py::arg("seed") = 1234,
       py::arg("referenceIonCharge") = 1,
       py::arg("batchMode") = true,
       py::arg("no_neutral_particles") = true,
       py::return_value_policy::reference)

    .def_static("GetInstance", []() {return BDSLinkTrackerInterface::GetInstance();},
                py::return_value_policy::reference)
    .def("Reset", &BDSLinkTrackerInterface::Reset)
    .def("PrepareBDSParticleDefinition", &BDSLinkTrackerInterface::PrepareBDSParticleDefinition)
    .def("GetReferenceParticleDefinition", [](BDSLinkTrackerInterface &ti) {
      return ti.GetReferenceParticleDefinition();
    }, py::return_value_policy::reference)
    .def("SetReferenceParticleDefinition",&BDSLinkTrackerInterface::SetReferenceParticleDefinition, py::arg("particleDefinition"))
    .def("GetBDSIMConfigFile", &BDSLinkTrackerInterface::GetBDSIMConfigFile)
    .def("GetReferenceParticlePDG", &BDSLinkTrackerInterface::GetReferenceParticlePDG)
    .def("GetReferenceParticleKineticEnergy", &BDSLinkTrackerInterface::GetReferenceParticleKineticEnergy)
    .def("GetRelativeEnergyCut", &BDSLinkTrackerInterface::GetRelativeEnergyCut)
    .def("GetSeed", &BDSLinkTrackerInterface::GetSeed)
    .def("GetReferenceIonCharge", &BDSLinkTrackerInterface::GetReferenceIonCharge)
    .def("GetBatchMode", &BDSLinkTrackerInterface::GetBatchMode)
    .def("GetMinimumKineticEnergy", &BDSLinkTrackerInterface::GetMinimumKineticEnergy)
    .def("SetNoNeutralParticles", &BDSLinkTrackerInterface::SetNoNeutralParticles)
    .def("GetNoNeutralParticles", &BDSLinkTrackerInterface::GetNoNeutralParticles)
    .def("GetBunchLink",&BDSLinkTrackerInterface::GetBunchLink,py::return_value_policy::reference)
    .def("GetBDSIMLink",&BDSLinkTrackerInterface::GetBDSIMLink,py::return_value_policy::reference)
    .def("AddParticle",[](BDSLinkTrackerInterface &ti, double x, double y, double px, double py,
                          double ct, double deltap, double chi, double chargeRatio,
                          double s, int trackid, int pdgid) {
      ti.AddParticle(x,y, px, py, ct, deltap, chi, chargeRatio, s, trackid, pdgid);
    }, py::arg("x"), py::arg("y"), py::arg("xp"), py::arg("yp"),
    py::arg("ct"), py::arg("deltap"), py::arg("chi"), py::arg("chargeRatio"),
    py::arg("s"), py::arg("trackid"), py::arg("pdgid"))
    .def("AddParticle",[](BDSLinkTrackerInterface &ti, double x, double y, double px, double py,
                          double pz, double t, double s,
                          int trackid, int pdgid) {
      ti.AddParticle(x, y, px, py, pz, t, s, trackid, pdgid);
    }, py::arg("x"), py::arg("y"), py::arg("px"), py::arg("py"), py::arg("pz"),
    py::arg("t"), py::arg("s"), py::arg("trackid"), py::arg("pdgid"))
    .def("AddParticles",[](BDSLinkTrackerInterface &ti,
                           std::vector<double> &x,
                           std::vector<double> &y,
                           std::vector<double> &xp,
                           std::vector<double> &yp,
                           std::vector<double> &ct,
                           std::vector<double> &deltap,
                           std::vector<double> &chi,
                           std::vector<double> &chargeRatio,
                           std::vector<double> &s,
                           std::vector<int> &trackid,
                           std::vector<int> &pdgid) {
      ti.AddParticles(x,y,xp,yp, ct, deltap, chi, chargeRatio, s, trackid, pdgid);
    }, py::arg("x"), py::arg("y"), py::arg("xp"), py::arg("yp"),
    py::arg("ct"), py::arg("deltap"), py::arg("chi"), py::arg("chargeRatio"),
    py::arg("s"), py::arg("trackid"), py::arg("pdgid"))
    .def("AddParticles",[](BDSLinkTrackerInterface &ti,
                           std::vector<double> x,
                           std::vector<double> y,
                           std::vector<double> px,
                           std::vector<double> py,
                           std::vector<double> pz,
                           std::vector<double> t,
                           std::vector<double> s,
                           std::vector<int> trackid,
                           std::vector<int> pdgid) {
      ti.AddParticles(x,y,px,py,pz,t,s,trackid,pdgid);
    }, py::arg("x"), py::arg("y"), py::arg("px"), py::arg("py"), py::arg("pz"),
    py::arg("t"), py::arg("s"), py::arg("trackid"), py::arg("pdgid"))
    .def("ClearData",&BDSLinkTrackerInterface::ClearData)
    .def("GetParticleDefinition", &BDSLinkTrackerInterface::GetParticleDefinition,py::keep_alive<1, 2>())
    .def("GetParticlePDGMass",&BDSLinkTrackerInterface::GetParticlePDGMass)
    .def("GetParticlePDGCharge", &BDSLinkTrackerInterface::GetParticlePDGCharge)
    .def("GetChargeRatio", &BDSLinkTrackerInterface::GetChargeRatio)
    .def("GetMassRatio",&BDSLinkTrackerInterface::GetMassRatio)
    .def("GetChi", &BDSLinkTrackerInterface::GetChi)
    .def("TrackXSuite",[](BDSLinkTrackerInterface *tracker_interface,
        int iElement,
        std::string elementName,
        py::object particles,
        int verbosity) {
      // Start from a clean link state, and never leave a half-filled bunch, sampler
      // hits or active-state flags behind when anything throws: they would be
      // picked up by the next call and corrupt the particles object.
      tracker_interface->ClearData();
      try {
        TrackXSuite(tracker_interface, iElement, elementName, particles, verbosity);
      } catch (...) {
        tracker_interface->ClearData();
        throw;
      }
    }, py::arg("iElement"), py::arg("elementName"), py::arg("particles"),
       py::arg("verbosity") = 0)
    .def("TrackRFTrack",[](BDSLinkTrackerInterface *tracker_interface, py::object bunch6d) {
      TrackRFTrack(tracker_interface, bunch6d);
    });
}

// Print all volumes with a dense material (i.e. not vacuum or air) inside lv,
// with their bounding box in the frame of the link box (beam axis at x = y = 0).
namespace {
void PrintSolidVolumes(const G4LogicalVolume* lv, const G4Transform3D& toBox,
                       int depth, std::ostream& out)
{
  if (depth > 6) {return;}
  for (size_t i = 0; i < lv->GetNoDaughters(); ++i)
    {
      const G4VPhysicalVolume* pv = lv->GetDaughter((G4int)i);
      const G4LogicalVolume* dlv  = pv->GetLogicalVolume();
      // daughter frame -> mother frame -> box frame
      G4Transform3D local(pv->GetObjectRotationValue(), pv->GetObjectTranslation());
      G4Transform3D tr = toBox * local;
      const G4Material* mat = dlv->GetMaterial();
      if (mat && mat->GetDensity() > 0.01*CLHEP::g/CLHEP::cm3)
        {
          G4ThreeVector pMin, pMax;
          dlv->GetSolid()->BoundingLimits(pMin, pMax);
          double xmin = std::numeric_limits<double>::max(), xmax = -xmin;
          double ymin = xmin, ymax = -xmin;
          for (int c = 0; c < 8; ++c)
            {
              HepGeom::Point3D<double> corner((c & 1) ? pMax.x() : pMin.x(),
                                              (c & 2) ? pMax.y() : pMin.y(),
                                              (c & 4) ? pMax.z() : pMin.z());
              corner = tr * corner;
              xmin = std::min(xmin, corner.x()); xmax = std::max(xmax, corner.x());
              ymin = std::min(ymin, corner.y()); ymax = std::max(ymax, corner.y());
            }
          out << "             volume " << pv->GetName() << ":  material " << mat->GetName()
              << ",  x [" << xmin/CLHEP::mm << ", " << xmax/CLHEP::mm << "] mm"
              << ",  y [" << ymin/CLHEP::mm << ", " << ymax/CLHEP::mm << "] mm\n";
        }
      PrintSolidVolumes(dlv, tr, depth + 1, out);
    }
}

// Print what BDSIM itself has for the selected link element: this is independent
// of the bookkeeping on the tracker side, to verify that particles are sent to the
// intended element.
void PrintLinkElementInfo(BDSIMLink* link, int iElement, const std::string& elementName,
                          int nActive, double arcLength)
{
  std::ostringstream out;
  out << std::setprecision(5);
  const int linkIndex = link->GetLinkIndex(elementName);
  out << "BDSIM link: tracking " << std::setw(7) << nActive << " particles through "
      << elementName << " (link ID " << linkIndex;
  if (linkIndex != iElement)
    {out << ", WARNING: tracker requested link ID " << iElement;}
  out << "):  length " << arcLength/CLHEP::m << " m";

  // Type and material from the component registry (BDSIM component names are the
  // link element names). An element that BDSIM failed to build shows up as drift.
  const BDSAcceleratorComponent* component = nullptr;
  for (const auto& kv : *BDSAcceleratorComponentRegistry::Instance())
    {
      if (kv.second && kv.second->GetName() == elementName)
        {component = kv.second; break;}
    }
  if (component)
    {out << ",  type " << component->GetType() << ",  material " << component->Material();}
  else
    {out << ",  type unknown (not found in the component registry)";}
  out << "\n";

  // The volumes actually placed in Geant4, with their position: the jaws show up here
  const BDSLinkComponent* linkComp = nullptr;
  for (int i = 0; ; ++i)
    {
      // the beamline holds link components only; nullptr past its end
      const BDSLinkComponent* lc = link->GetLinkComponent(i);
      if (!lc) {break;}
      if (lc->LinkName() == elementName) {linkComp = lc; break;}
    }
  if (linkComp && linkComp->Component() && linkComp->Component()->GetContainerLogicalVolume())
    {PrintSolidVolumes(linkComp->Component()->GetContainerLogicalVolume(), G4Transform3D(), 0, out);}
  else
    {out << "             (geometry of " << elementName << " not found in the link beamline)\n";}
  std::cout << out.str() << std::flush;
}

// Crystal channelling (G4ChannelingFastSimModel) is done with Geant4 fast
// simulation, which BDSIM's channelling physics lists only activate for e-, e+ and
// protons (BDS::ChannellingPhysicsComplete). The model itself works for any charged
// particle, but ions (which all use the process manager of GenericIon) never reach
// it, so they do not channel. When fast simulation is active for protons, it is
// activated here for ions as well, in the same way as G4FastSimulationPhysics does
// (once, after the physics is built).
// TODO: fix this in BDSIM (activate fast simulation for GenericIon as well), and
// remove this workaround.
bool HasFastSimulation(G4ProcessManager* pm)
{
  if (!pm)
    {return false;}
  G4ProcessVector* pv = pm->GetProcessList();
  for (G4int i = 0; i < (G4int)pv->size(); ++i)
    {
      if (dynamic_cast<G4FastSimulationManagerProcess*>((*pv)[i]))
        {return true;}
    }
  return false;
}

void ActivateIonChannelling()
{
  static bool done = false;
  if (done)
    {return;}
  done = true;
  G4ProcessManager* pmProton = G4Proton::Proton()->GetProcessManager();
  G4ProcessManager* pmIon    = G4GenericIon::GenericIon()->GetProcessManager();
  if (pmIon && HasFastSimulation(pmProton) && !HasFastSimulation(pmIon))
    {G4FastSimulationHelper::ActivateFastSimulation(pmIon);}
}

// G4Fragment::CalculateExcitationEnergy() prints a warning (five lines) whenever
// the excitation energy of a fragment is below -10 eV. For heavy ions at LHC
// energies this happens all the time through rounding alone: the excitation energy
// is sqrt(E^2 - p^2) - M, and at E ~ 5e8 MeV the rounding error on E^2 - p^2 gives
// errors of ~1e-4 MeV. Geant4 then sets the excitation energy to 0, which is the
// correct result, so these warnings are dropped here (only these: all other output
// is passed on unchanged).
void InstallG4coutFilter()
{
  static bool done = false;
  if (done)
    {return;}
  done = true;
  static G4coutDestination destination;  // default destination: prints to std::cout
  static int skip = 0;                   // lines of a warning still to be dropped
  destination.AddCoutTransformer([](G4String& msg) -> G4bool {
      if (msg.find("G4Fragment::CalculateExcitationEnergy(): WARNING") != std::string::npos)
        {
          skip = 4;  // the fragment printout (3 lines) and an empty line
          return false;
        }
      if (skip > 0)
        {
          --skip;
          const bool blank = msg.find_first_not_of(" \t\r\n") == std::string::npos;
          if (msg.rfind("Fragment: A =", 0) == 0 || msg.rfind("          P = (", 0) == 0
              || msg.rfind("    #spin=", 0) == 0 || blank)
            {return false;}
          skip = 0;
        }
      return true;
    });
  G4iosSetDestination(&destination);
}
} // namespace

void TrackXSuite(BDSLinkTrackerInterface *tracker_interface,
                            int iElement,
                            std::string elementName,
                            py::object particles,
                            int verbosity) {

    py::array_t<double> x = py::cast<py::array_t<double>>(particles.attr("x"));
    py::array_t<double> y = py::cast<py::array_t<double>>(particles.attr("y"));
    py::array_t<double> px = py::cast<py::array_t<double>>(particles.attr("px"));
    py::array_t<double> py = py::cast<py::array_t<double>>(particles.attr("py"));
    py::array_t<double> ct = py::cast<py::array_t<double>>(particles.attr("zeta"));
    py::array_t<double> deltap = py::cast<py::array_t<double>>(particles.attr("delta"));
    py::array_t<double> chi = py::cast<py::array_t<double>>(particles.attr("chi"));
    py::array_t<double> charge_ratio = py::cast<py::array_t<double>>(particles.attr("charge_ratio"));
    py::array_t<double> massratio = py::cast<py::array_t<double>>(particles.attr("mass_ratio"));
    py::array_t<double> s = py::cast<py::array_t<double>>(particles.attr("s"));
    py::array_t<int64_t> pdgid = py::cast<py::array_t<int64_t>>(particles.attr("pdg_id"));
    py::array_t<int64_t> trackid = py::cast<py::array_t<int64_t>>(particles.attr("particle_id"));
    py::array_t<int64_t> parentid = py::cast<py::array_t<int64_t>>(particles.attr("parent_particle_id"));
    py::array_t<int64_t> state = py::cast<py::array_t<int64_t>>(particles.attr("state"));
    py::array_t<int64_t> at_element = py::cast<py::array_t<int64_t>>(particles.attr("at_element"));
    py::array_t<int64_t> at_turn = py::cast<py::array_t<int64_t>>(particles.attr("at_turn"));

    // fill BDSLinkBunch
    auto bunch = tracker_interface->GetBunchLink();
    auto ref = tracker_interface->GetReferenceParticleDefinition();
    // bunch->AddParticle();

    // Fast pointers to numpy arrays
    auto get = [&](const char* name) {
        return py::cast<py::array_t<double>>(particles.attr(name)).mutable_data();
    };
    auto geti = [&](const char* name) {
        return py::cast<py::array_t<int64_t>>(particles.attr(name)).mutable_data();
    };

    double* x_ptr = get("x");
    double* y_ptr = get("y");
    double* px_ptr = get("px");
    double* py_ptr = get("py");
    double* zeta_ptr = get("zeta");
    double* delta_ptr = get("delta");
    double* chi_ptr = get("chi");
    double* charge_ratio_ptr = get("charge_ratio");
    double* s_ptr = get("s");
    int64_t* pdgid_ptr = geti("pdg_id");
    int64_t* trackid_ptr = geti("particle_id");
    int64_t* state_ptr = geti("state");
    int64_t* at_element_ptr = geti("at_element");
    int64_t* at_turn_ptr = geti("at_turn");

    // Add particles to the bunch to be tracked
    ssize_t n = py::cast<py::array_t<double>>(particles.attr("x")).size();
    auto& active_state = tracker_interface->GetParticleActiveState();
    int n_active = 0;
    size_t n_unphysical = 0;
    for (ssize_t i = 0; i < n; ++i) {
        if (state_ptr[i] < 1) {  // there can be different types of alive states in xsuite
            active_state.push_back(false);
            continue;
        }
        else {
            G4double q = charge_ratio_ptr[i] * ref->Charge();
            G4double mass_ratio = charge_ratio_ptr[i] / chi_ptr[i];
            G4double p = ref->Momentum() * (delta_ptr[i] + 1) * mass_ratio;

            // xsuite normalizes px, py and delta to the particle's own mass
            // (px = Px / (P0 m/m0)), so the slope is px / (1 + delta), without mass ratio
            G4double oneplusdelta = 1 + delta_ptr[i];
            G4double xp = px_ptr[i] / oneplusdelta;
            G4double yp = py_ptr[i] / oneplusdelta;
            // A particle whose transverse momentum is not smaller than its total momentum
            // cannot be built (BDSBunch::CalculateZp would throw and abort the whole
            // call). This can only result from unphysical upstream tracking: mark it as
            // lost in this element (same marker as particles lost in BDSIM) and skip it.
            if (!std::isfinite(xp) || !std::isfinite(yp) || xp*xp + yp*yp >= 1) {
                set_element<int64_t>(state, i, -333);
                active_state.push_back(false);
                ++n_unphysical;
                continue;
            }
            G4double zp = BDSBunch::CalculateZp(xp, yp, 1);

            auto partDef = PrepareBDSParticleDefinition(pdgid_ptr[i], p, 0, q);
            G4double t = - zeta_ptr[i] * CLHEP::m / (ref->Beta() * CLHEP::c_light);

            BDSParticleCoordsFull coords(x_ptr[i] * CLHEP::m,
                                         y_ptr[i] * CLHEP::m,
                                         0,
                                         xp, yp, zp,
                                         t,
                                         0,
                                         partDef->TotalEnergy(),
                                         1);
            active_state.push_back(true);
            ++n_active;
            bunch->AddParticle(partDef, coords, trackid_ptr[i], trackid_ptr[i]);
        }
    }

    // One active-state flag per slot: the write-back below relies on this
    if (active_state.size() != (size_t)n)
      {throw std::runtime_error("TrackXSuite> Inconsistent link state (stale data from a previous call?)");}

    if (n_unphysical > 0) {
        std::cerr << "TrackXSuite> WARNING: " << n_unphysical << " particles with unphysical "
                  << "angles (xp^2 + yp^2 >= 1) entering element " << elementName
                  << " were marked as lost." << std::endl;
    }

    // run n particles
    auto link = tracker_interface->GetBDSIMLink();
    try {
        link->SelectLinkElement(elementName);  // throws if the element doesn't exist
    } catch (...) {
        tracker_interface->ClearData(); // don't leave a filled bunch behind for the next call
        throw;
    }
    // Length of the selected element. Look it up by name: iElement is a link ID
    // (output sampler ID), which is not the index in the link beamline.
    const double collLength = link->GetArcLengthOfLinkElement(elementName);
    if (collLength < 0)
      {throw std::runtime_error("TrackXSuite> Could not determine the length of element " + elementName);}
    if (verbosity > 0)
      {PrintLinkElementInfo(link, iElement, elementName, n_active, collLength);}
    ActivateIonChannelling();
    InstallG4coutFilter();
    link->BeamOn((G4int)bunch->Size());
    // get sampler data
    const BDSHitsCollectionSamplerLink* hits = link->SamplerHits();
    size_t hitsCount = 0;
    if (hits)
    {
        hitsCount = hits->GetSize();
    }
    // Count the number of secondary particles that will be written back
    // (neutral secondaries are skipped below if requested, so don't count them)
    size_t secondaryCount = 0;
    for (size_t i = 0; i < hitsCount; i++)
    {
        auto hit = (*hits)[i];
        if (hit->externalParticleID == hit->externalParentID) { continue; }
        if (hit->charge == 0 && tracker_interface->GetNoNeutralParticles()) { continue; }
        secondaryCount++;
    }

    // The output arrays have slots for all particles, regardless of lost or not, and for secondary particles
    size_t output_size = x.size();

    // Secondaries can only be written into unallocated slots of the particles object
    // (xsuite marks these with state == LAST_INVALID_STATE). Lost particles from upstream
    // elements also sit after the active ones, so we must not simply write from n_active.
    // Check the capacity *before* writing anything, so that a failure leaves the particles
    // object untouched.
    const int64_t LAST_INVALID_STATE = -999999999;
    std::vector<size_t> free_slots;
    free_slots.reserve(secondaryCount);
    for (size_t j = 0; j < output_size && free_slots.size() < secondaryCount; j++)
    {
        if (state_ptr[j] == LAST_INVALID_STATE) { free_slots.push_back(j); }
    }
    if (free_slots.size() < secondaryCount)
    {
        tracker_interface->ClearData();
        throw std::runtime_error("TrackXSuite> Not enough free slots in the particles object for the "
                                 + std::to_string(secondaryCount) + " secondaries returned by BDSIM (only "
                                 + std::to_string(free_slots.size()) + " free). Increase _capacity.");
    }

    // Loop through the particles in the *original* bunch - the primaries
    size_t hits_index = 0;
    bool prim_survied = false;
    // double sum_deltaplusone_sec = 0.0;
    double sum_secondary_energy = 0.0;

    size_t free_slot_index = 0;
    size_t n_multiple_primary_hits = 0;
    for (size_t i=0; i < active_state.size(); i++){
        if (!active_state.at(i)){
            continue; // This was an inactive particle that hasn't been processed, do not change it
        }

        // ID of the primary particle. This is the external ID it was added to the bunch
        // with (see AddParticle above), so read it directly instead of advancing the
        // bunch: BeamOn already consumed it, so GetNextParticle would wrap around
        // (printing "looping to start of bunch" on every call).
        const int prim_part_id = static_cast<int>(trackid_ptr[i]);

        // Now start looping over the hits - the particles to be returned to the tracker
        // These can be primary or secondary particles. Each primary can produce 0, 1, or 2+ products
        // The products need to be sorted to keep the array order - surviving primary particles are all
        // filled in first. If a primary didn't survive, keep the original coordinates and make it inactive.
        // The hits are ordered by primary event, so just need one loop.
        while (hits_index < hitsCount)
        {
            BDSHitSamplerLink* hit = (*hits)[hits_index];
            if (hit->externalParentID != prim_part_id) { // The hits corresponding to the current primary are exhausted
                break;
            }

            const BDSParticleCoordsFull &coords = hit->coords;

            double qratio = hit->charge / ref->Charge();
            double mratio = hit->mass / ref->Mass();
            double dp = (hit->momentum / mratio - ref->Momentum()) / ref->Momentum();

            /// Need to compensate for the geometry construction in BDSIM
            /// There is a safety margin that is added to the collimator legnth
            double collMargin = 2.5 * BDSSamplerCustom::ChordLength();
            double zt = ref->Beta() * CLHEP::c_light * ((collLength + collMargin) / (CLHEP::c_light * ref->Beta()) - coords.T);
            double oneplusdelta = 1 + dp;  // xsuite px = xp (1 + delta), see above

            auto track_id = hit->externalParticleID;
            auto parent_id = hit->externalParentID;
            const int64_t pdg_id = GroundStatePDGID(hit->pdgID);

            if (track_id == parent_id){
                // This is a primary particle as its parent is itself
                if (prim_survied) {
                    // The primary crossed more than one output sampler: it was not stopped
                    // at the end of the element (e.g. missing user limits in the physics
                    // list). Keep the first hit; later ones belong to other elements.
                    ++n_multiple_primary_hits;
                    ++hits_index;
                    continue;
                }
                prim_survied = true;

                set_element<double>(s, i, s_ptr[i] + collLength / CLHEP::m);
                set_element<double>(x, i, coords.x / CLHEP::m);
                set_element<double>(px, i, coords.xp * oneplusdelta);
                set_element<double>(y, i, coords.y / CLHEP::m);
                set_element<double>(py, i, coords.yp * oneplusdelta);
                set_element<double>(ct, i, zt / CLHEP::m);
                set_element<double>(deltap, i, dp);
                set_element<int64_t>(pdgid, i, pdg_id);
                set_element<double>(charge_ratio, i, qratio);
                set_element<double>(chi, i, qratio / mratio);
                // parent_particle_id of a surviving primary is left untouched: if it was itself
                // a secondary from an upstream element, it must keep pointing to its own parent.
                set_element<int64_t>(state, i, 1);

            }
            else
            {
                if (qratio == 0 && tracker_interface->GetNoNeutralParticles()) {
                    ++hits_index;
                    continue;
                }
                // Secondary particles are written into free (unallocated) slots. Capacity was checked
                // above. Only the fields below are set here; the remaining per-particle fields
                // (p0c, beta0, gamma0, weight, particle_id, rng, ...) are filled in on the Python side.
                size_t prod_write_index = free_slots[free_slot_index++];
                set_element<double>(s, prod_write_index, s_ptr[i] + collLength / CLHEP::m);
                set_element<double>(x, prod_write_index, coords.x / CLHEP::m);
                set_element<double>(px, prod_write_index, coords.xp * oneplusdelta);
                set_element<double>(y, prod_write_index, coords.y / CLHEP::m);
                set_element<double>(py, prod_write_index, coords.yp * oneplusdelta);
                set_element<double>(ct, prod_write_index, zt / CLHEP::m);
                set_element<double>(deltap, prod_write_index, dp);
                set_element<int64_t>(pdgid, prod_write_index, pdg_id);
                set_element<int64_t>(at_element, prod_write_index, at_element_ptr[i]); // same element as the parent
                set_element<int64_t>(at_turn, prod_write_index, at_turn_ptr[i]);
                set_element<double>(charge_ratio, prod_write_index, qratio);
                set_element<int64_t>(parentid, prod_write_index, parent_id);
                set_element<int64_t>(state, prod_write_index, 1);
                double chi_val = qratio / mratio;
                set_element<double>(chi, prod_write_index, chi_val);

                sum_secondary_energy += std::sqrt(std::pow(hit->momentum,2) + std::pow(hit->mass,2));
            }

            hits_index++;
        }

        if (!prim_survied) // Primary didn't survive - set inactive
        {
            // i < output_size is guaranteed, as active_state has exactly one entry per slot.
            // -333 is a generic "lost in BDSIM" marker; the calling code (e.g. xcoll) maps it
            // to its own loss states.
            set_element<int64_t>(state, i, -333); // inactive
            // Correct the energy of the lost primary particle to account for the production of secondaries
            // The effective delta is such that the lost particle has the effective delta
            // which corresponds to the energy in - energy out for this primary

            // reconstruct the incoming primary particle energy
            G4double mass_ratio_prim = charge_ratio_ptr[i] / chi_ptr[i];
            G4double p_prim = ref->Momentum() * (delta_ptr[i] + 1) * mass_ratio_prim;
            G4double mass_prim = mass_ratio_prim * ref->Mass();
            G4double energy_prim = std::sqrt(std::pow(p_prim, 2) + std::pow(mass_prim, 2));

            // compute the effective delta
            G4double energy_diff = energy_prim - sum_secondary_energy;
            G4double p_eff;
            G4double squared_energy_mass_diff = std::pow(energy_diff, 2) - std::pow(mass_prim, 2);
            if (squared_energy_mass_diff < 0){
                // This means that the total energy escaping includes part of the rest mass of
                // the primary. Tolerate the error for now, as otherwise need to adjust also the
                // mass and PDG id of the lost primary particle
                p_eff = 0;
            }
            else
            {
                p_eff = std::sqrt(squared_energy_mass_diff);
            }
            G4double delta_eff = (p_eff / mass_ratio_prim - ref->Momentum()) / ref->Momentum();
            set_element<double>(deltap, i, delta_eff);
        }
        prim_survied = false; // reset for next particle
        sum_secondary_energy = 0.0;
    }
    if (n_multiple_primary_hits > 0) {
        std::cerr << "TrackXSuite> WARNING: " << n_multiple_primary_hits << " extra sampler hits of "
                  << "primaries in element " << elementName << " were ignored: particles are not "
                  << "stopped at the end of the element. Check that the physics list includes "
                  << "the BDSIM cuts and limits process." << std::endl;
    }
    // clean BDSLinkBunch of particles
    tracker_interface->ClearData();
}

void TrackRFTrack(BDSLinkTrackerInterface *tracker_interface, py::object bunch6d) {
  py::print("TrackRFTrack> Bunch6d::", bunch6d);
  auto size_method = bunch6d.attr("size");
  py::print("TrackRFTrack> Bunch6d::size",size_method());

  auto bdsim_link = tracker_interface->GetBDSIMLink();

  auto refPDG = tracker_interface->GetReferenceParticleDefinition()->PDGID();
  // clear sampler hits (do this first and not at end as sampler data will
  // no longer available in python)
  bdsim_link->ClearSamplerHits();

  int nparticle = bunch6d.attr("size")().cast<int>();

  for(int i = 0; i < nparticle ;i++) {
    // py::print("TrackRFTrack> ",i);

    auto p = bunch6d.attr("get_particle")(i);
    auto x = py::cast<double>(p.attr("x"));
    auto y = py::cast<double>(p.attr("y"));
    auto p4v = (py::cast<py::array_t<double>>(p.attr("get_four_momentum")())).unchecked<2>();
    // TODO time needs to calculated wrt to reference particle
    // auto t = py::cast<double>(p.attr("t"));
    auto pdgID = py::cast<int>(p.attr("pdg_id")); // initiall this will be zero and add particle will set to reference particle

    tracker_interface->AddParticle(x*CLHEP::mm, y*CLHEP::mm, // x, y
                                   p4v(1,0), p4v(2,0), p4v(3,0), // px, px, pz
                                   0, // ct
                                   0, // s
                                   i, pdgID); // trackID, parent;
    if (pdgID == 0) {
      p.attr("pdg_id") = py::cast(refPDG);
    }
  }
  bdsim_link->BeamOn(nparticle);

  // loop over sampler hits and update bunch
  auto sh = bdsim_link->SamplerHits();

  // set all bunch particles as if they didn't make it (s)
  //auto endpoint = (*sh)[0]->coords.s;
  auto endpoint = 0.1;
  for(int i=0; i<nparticle; i++) {
    auto p = bunch6d.attr("get_particle")(i);
    p.attr("S_lost") = py::cast(endpoint);
  }


  auto append_method = bunch6d.attr("append");
  std::vector<double> new_particle_data = {0, 0, 0, 0, 0, 0, 0, 0, 0};

  for(std::size_t i=0; i<sh->entries(); i++) {
    auto h = (*sh)[i];
    // std::cout << i << " " << h->externalParentID << std::endl;
    auto externalID = h->externalParentID;
    auto trackID = h->trackID;
    auto p = bunch6d.attr("get_particle")(externalID);
    if(trackID == 1) { // existing partucke
      p.attr("x") = py::cast(h->coords.x);
      p.attr("y") = py::cast(h->coords.y);
      p.attr("xp") = py::cast(h->coords.xp);
      p.attr("yp") = py::cast(h->coords.yp);
      p.attr("S_lost") = py::cast(std::nan(""));
      // TODO colulate absolute T for particle
    }
    else { // new particle

      new_particle_data[0] = h->coords.x;
      new_particle_data[1] = h->coords.xp;
      new_particle_data[2] = h->coords.y;
      new_particle_data[3] = h->coords.yp;

      // look up particle mass and charge
      G4ParticleDefinition* particle = G4ParticleTable::GetParticleTable()->FindParticle(h->pdgID);
      new_particle_data[6] = particle->GetPDGMass();
      new_particle_data[7] = particle->GetPDGCharge();
      new_particle_data[8] = py::cast<double>(p.attr("N"))*h->coords.weight;

      // append bdsim generated particle
      int idx_insert = bunch6d.attr("size")().cast<int>();
      py::array_t<double> new_particle_arr(new_particle_data.size(), new_particle_data.data());
      append_method(new_particle_arr);

      // need to set particle as not lost
      p = bunch6d.attr("get_particle")(idx_insert);
      p.attr("S_lost") = py::cast(std::nan(""));

      // set particle momentum
      p.attr("Pc") = py::cast(h->momentum);

      // set particle pdg
      p.attr("pdg_id") = py::cast(h->pdgID);

      // TODO colulate absolute T for new particles particle
      // set particle time

      // set created particle lifetime for RF Track
      auto lifetime = particle->GetPDGLifeTime();
      auto mass = particle->GetPDGMass();
      if (lifetime < 0) { // set infinite lifetime
        p.attr("lifetime") = std::numeric_limits<double>::infinity();
      }
      else if(mass > 0) { // need to have non zero mass
        //vauto boosted_lifetime = h->coords.totalEnergy/mass*lifetime;
        // auto sampled_lifetime = -boosted_lifetime * std::log(G4UniformRand());
        // std::cout << lifetime << " " << h->coords.totalEnergy << " " << mass << std::endl;
        // p.attr("lifetime") = py::cast(boosted_lifetime);

        p.attr("lifetime") = std::numeric_limits<double>::infinity();
      }
    }
  }
}