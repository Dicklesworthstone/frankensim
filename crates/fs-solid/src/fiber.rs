//! Fiber-section beam-columns (bead tfz.14): cross-sections
//! discretized into fibers, each carrying a uniaxial fs-material law —
//! Mander confined/unconfined concrete plus Menegotto–Pinto steel for
//! the reinforced-concrete frame use-case. The section state update
//! `ε(y) = ε₀ − y·κ → (N, M, 2×2 tangent)` is the seismic flagship's
//! inner hot loop; [`update_sections_batched`] runs MANY sections per
//! call with the 2×2 tangent solves through fs-la's batched Cholesky
//! (the §15.2 pairing), and the battery ledgers the measured
//! throughput.
//!
//! Cyclic correctness is the point: the hysteresis battery drives an
//! RC section through growing curvature cycles and checks loop
//! closure, positive per-cycle dissipation, and the peak moment
//! against a hand capacity estimate.

use crate::SolidError;
use fs_la::batched::{BatchMat, BatchVec, batch_cholesky, batch_cholesky_solve};
use fs_material::creep::plastic::uniaxial::{UniaxialNortonJ2, UniaxialNortonJ2State};
use fs_material::creep::recovery::plastic::uniaxial::{
    UniaxialRecoveringNortonJ2, UniaxialRecoveringNortonJ2State,
};
use fs_material::creep::recovery::uniaxial::{
    UniaxialRecoveringNorton, UniaxialRecoveringNortonState,
};
use fs_material::creep::thermal::uniaxial::{
    UniaxialThermalNorton, UniaxialThermalNortonJ2, UniaxialThermalNortonJ2State,
    UniaxialThermalNortonState,
};
use fs_material::creep::thermal::uniaxial::recovery::{
    UniaxialThermalRecoveringNorton, UniaxialThermalRecoveringNortonJ2,
    UniaxialThermalRecoveringNortonJ2State, UniaxialThermalRecoveringNortonState,
};
use fs_material::creep::uniaxial::{UniaxialNorton, UniaxialNortonState};
use fs_material::fiber::{ManderConcrete, ManderState, MenegottoPintoSteel, MpState, Uniaxial};
use fs_material::ductile::uniaxial::UniaxialDuctileJ2State;
use fs_material::plastic::chaboche::axial::{UniaxialChaboche, UniaxialChabocheState};
use fs_material::plastic::tabulated::rate::axial::UniaxialTabulatedOverstress;
use fs_material::plastic::thermal::tabulated::axial::{
    UniaxialThermalTabulated, UniaxialThermalTabulatedState,
};
use fs_material::plastic::reduction::{
    UniaxialJ2, UniaxialJ2State, UniaxialTabulated, UniaxialVoce,
};

/// Bounded runtime dispatch of exact native multi-backstress axial laws.
pub mod chaboche;
/// Native regularized ductile axial response with an explicit physical interval.
pub mod ductile;
/// Material-bound fiber laws, RC section geometry and coherent member mass.
pub mod material;
use chaboche::{ChabocheFiberLaw, ChabocheFiberState};
use chaboche::thermal::{ThermalChabocheFiberLaw, ThermalChabocheFiberState};
use ductile::DuctileFiberLaw;

/// A fiber's constitutive card.
#[derive(Debug, Clone)]
pub enum FiberLaw {
    /// Free-sided axial Norton creep for one explicit physical interval.
    /// Original material receipts remain attached to source-bound laws.
    Norton(UniaxialNorton),
    /// Permanent primary creep with independent resistance recovery during dwell.
    RecoveringNorton(UniaxialRecoveringNorton),
    /// Simultaneous axial Norton creep and isotropic J2 with separate memories.
    NortonJ2(UniaxialNortonJ2),
    /// Simultaneous permanent creep/plasticity with separate recovering resistance.
    RecoveringNortonJ2(UniaxialRecoveringNortonJ2),
    /// Sourced thermal expansion plus Norton creep and reciprocal heat.
    ThermalNorton(UniaxialThermalNorton),
    /// Sourced thermal expansion plus simultaneous Norton/J2 response.
    ThermalNortonJ2(UniaxialThermalNortonJ2),
    /// Sourced thermal expansion, primary-resistance recovery and reciprocal heat.
    ThermalRecoveringNorton(UniaxialThermalRecoveringNorton),
    /// Sourced thermal expansion with joint recovering creep/J2 and native heat.
    ThermalRecoveringNortonJ2(UniaxialThermalRecoveringNortonJ2),
    /// Free-sided axial restriction of the shared three-dimensional J2 law.
    /// Isotropic and Prager hardening retain separate reversal histories.
    J2(UniaxialJ2),
    /// Native axial J2 initiation and regularized damage with retained energy.
    /// This prepared law consumes one explicitly owned physical interval.
    DuctileJ2(DuctileFiberLaw),
    /// Native Voce isotropic saturation with retained linear Prager translation.
    Voce(UniaxialVoce),
    /// Native bounded measured-curve hardening, with retained Prager translation.
    Tabulated(UniaxialTabulated),
    /// Original bounded hardening with associative Perzyna flow over one physical interval.
    TabulatedOverstress(UniaxialTabulatedOverstress),
    /// Original temperature-dependent hardening, expansion and optional Perzyna
    /// flow over one prepared physical interval.
    ThermalTabulated(UniaxialThermalTabulated),
    /// Source-bound thermal cyclic flow with one through eight independent branches.
    ThermalChaboche(ThermalChabocheFiberLaw),
    /// Independent nonlinear backstress recovery with native Voce isotropic hardening.
    Chaboche(UniaxialChaboche),
    /// Two through eight independently evolving native recovery branches.
    ChabocheBranches(ChabocheFiberLaw),
    /// Menegotto–Pinto steel.
    Steel(MenegottoPintoSteel),
    /// Mander concrete (confined or unconfined by parameters).
    Concrete(ManderConcrete),
}

/// A fiber's committed state.
#[derive(Debug, Clone)]
pub enum FiberState {
    /// Signed permanent axial creep and accumulated equivalent creep.
    Norton(UniaxialNortonState),
    /// Permanent signed/path histories and a separate recovering resistance.
    RecoveringNorton(UniaxialRecoveringNortonState),
    /// Signed creep/plastic strains and their separate accumulated path lengths.
    NortonJ2(UniaxialNortonJ2State),
    /// Both permanent histories and independently recovering creep resistance.
    RecoveringNortonJ2(UniaxialRecoveringNortonJ2State),
    /// Signed creep history plus the accepted sourced temperature.
    ThermalNorton(UniaxialThermalNortonState),
    /// Signed creep/plastic histories plus the accepted sourced temperature.
    ThermalNortonJ2(UniaxialThermalNortonJ2State),
    /// Permanent creep, recovering resistance and accepted source temperature.
    ThermalRecoveringNorton(UniaxialThermalRecoveringNortonState),
    /// Both permanent histories, recovering resistance and accepted temperature.
    ThermalRecoveringNortonJ2(UniaxialThermalRecoveringNortonJ2State),
    /// Signed axial plastic strain and accumulated equivalent plastic strain.
    J2(UniaxialJ2State),
    /// Native plasticity, irreversible initiation/damage and separate work accounts.
    DuctileJ2(UniaxialDuctileJ2State),
    /// The same exact axial plastic subspace, with a saturating yield radius.
    Voce(UniaxialJ2State),
    /// Signed axial plastic strain and accumulated strain on the supplied table.
    Tabulated(UniaxialJ2State),
    /// Permanent axial/accumulated plastic flow; accepted viscous overstress may persist.
    TabulatedOverstress(UniaxialJ2State),
    /// Permanent signed/path plastic history and its accepted source temperature.
    ThermalTabulated(UniaxialThermalTabulatedState),
    /// Accepted temperature, plastic flow and every independent backstress.
    ThermalChaboche(ThermalChabocheFiberState),
    /// Signed axial plastic history and independently retained backstress.
    Chaboche(UniaxialChabocheState),
    /// Ordered native histories, with the same branch count as the fiber law.
    ChabocheBranches(ChabocheFiberState),
    /// Steel state.
    Steel(MpState),
    /// Concrete state.
    Concrete(ManderState),
}

/// One fiber: position, area, law.
#[derive(Debug, Clone)]
pub struct Fiber {
    /// Distance from the section reference axis (bending arm).
    pub y: f64,
    /// Tributary area.
    pub area: f64,
    /// The constitutive card.
    pub law: FiberLaw,
}

impl From<UniaxialNorton> for FiberLaw {
    fn from(law: UniaxialNorton) -> Self { Self::Norton(law) }
}
impl From<UniaxialRecoveringNorton> for FiberLaw {
    fn from(law: UniaxialRecoveringNorton) -> Self { Self::RecoveringNorton(law) }
}
impl From<UniaxialNortonJ2> for FiberLaw {
    fn from(law: UniaxialNortonJ2) -> Self { Self::NortonJ2(law) }
}
impl From<UniaxialRecoveringNortonJ2> for FiberLaw {
    fn from(law: UniaxialRecoveringNortonJ2) -> Self { Self::RecoveringNortonJ2(law) }
}
impl From<UniaxialThermalNorton> for FiberLaw {
    fn from(law: UniaxialThermalNorton) -> Self { Self::ThermalNorton(law) }
}
impl From<UniaxialThermalNortonJ2> for FiberLaw {
    fn from(law: UniaxialThermalNortonJ2) -> Self { Self::ThermalNortonJ2(law) }
}
impl From<UniaxialThermalRecoveringNorton> for FiberLaw {
    fn from(law: UniaxialThermalRecoveringNorton) -> Self { Self::ThermalRecoveringNorton(law) }
}
impl From<UniaxialThermalRecoveringNortonJ2> for FiberLaw {
    fn from(law: UniaxialThermalRecoveringNortonJ2) -> Self { Self::ThermalRecoveringNortonJ2(law) }
}
impl From<UniaxialTabulated> for FiberLaw {
    fn from(law: UniaxialTabulated) -> Self { Self::Tabulated(law) }
}
impl From<DuctileFiberLaw> for FiberLaw {
    fn from(law: DuctileFiberLaw) -> Self { Self::DuctileJ2(law) }
}
impl From<UniaxialTabulatedOverstress> for FiberLaw {
    fn from(law: UniaxialTabulatedOverstress) -> Self { Self::TabulatedOverstress(law) }
}
impl From<UniaxialThermalTabulated> for FiberLaw {
    fn from(law: UniaxialThermalTabulated) -> Self { Self::ThermalTabulated(law) }
}

impl FiberLaw {
    /// These native returns consume their selected physical interval exactly once.
    pub(crate) fn requires_physical_duration(&self) -> bool {
        matches!(self,
            Self::Norton(_) | Self::RecoveringNorton(_) | Self::NortonJ2(_)
                | Self::RecoveringNortonJ2(_) | Self::ThermalNorton(_)
                | Self::ThermalNortonJ2(_) | Self::ThermalRecoveringNorton(_)
                | Self::ThermalRecoveringNortonJ2(_) | Self::TabulatedOverstress(_)
                | Self::ThermalTabulated(_) | Self::ThermalChaboche(_)
                | Self::DuctileJ2(_)
        )
    }
}

/// A section: fibers plus committed states.
#[derive(Debug, Clone)]
pub struct Section {
    /// The fibers.
    pub fibers: Vec<Fiber>,
    /// Committed states, one per fiber.
    pub states: Vec<FiberState>,
}

/// A section's response at a strain state.
#[derive(Debug, Clone, Copy, Default)]
pub struct SectionState {
    /// Axial force.
    pub n: f64,
    /// Bending moment.
    pub m: f64,
    /// Tangent [[dN/dε₀, dN/dκ], [dM/dε₀, dM/dκ]].
    pub tangent: [[f64; 2]; 2],
}

impl Section {
    /// Build with virgin states.
    #[must_use]
    pub fn new(fibers: Vec<Fiber>) -> Section {
        let states = fibers
            .iter()
            .map(|f| match &f.law {
                FiberLaw::Norton(l) => FiberState::Norton(l.initial_state()),
                FiberLaw::RecoveringNorton(l) => FiberState::RecoveringNorton(l.initial_state()),
                FiberLaw::NortonJ2(l) => FiberState::NortonJ2(l.initial_state()),
                FiberLaw::RecoveringNortonJ2(l) => FiberState::RecoveringNortonJ2(l.initial_state()),
                FiberLaw::ThermalNorton(l) => FiberState::ThermalNorton(l.initial_state()),
                FiberLaw::ThermalNortonJ2(l) => FiberState::ThermalNortonJ2(l.initial_state()),
                FiberLaw::ThermalRecoveringNorton(l) => FiberState::ThermalRecoveringNorton(l.initial_state()),
                FiberLaw::ThermalRecoveringNortonJ2(l) => FiberState::ThermalRecoveringNortonJ2(l.initial_state()),
                FiberLaw::J2(l) => FiberState::J2(l.initial_state()),
                FiberLaw::DuctileJ2(l) => FiberState::DuctileJ2(l.initial_state()),
                FiberLaw::Voce(l) => FiberState::Voce(l.initial_state()),
                FiberLaw::Tabulated(l) => FiberState::Tabulated(l.initial_state()),
                FiberLaw::TabulatedOverstress(l) => FiberState::TabulatedOverstress(l.initial_state()),
                FiberLaw::ThermalTabulated(l) => FiberState::ThermalTabulated(l.initial_state()),
                FiberLaw::ThermalChaboche(l) => FiberState::ThermalChaboche(l.initial_state()),
                FiberLaw::Chaboche(l) => FiberState::Chaboche(l.initial_state()),
                FiberLaw::ChabocheBranches(l) => FiberState::ChabocheBranches(l.initial_state()),
                FiberLaw::Steel(l) => FiberState::Steel(l.initial_state()),
                FiberLaw::Concrete(l) => FiberState::Concrete(l.initial_state()),
            })
            .collect();
        Section { fibers, states }
    }

    /// Response at (ε₀, κ) from the COMMITTED states (pure — commit
    /// separately).
    #[must_use]
    pub fn respond(&self, eps0: f64, kappa: f64) -> SectionState {
        self.response(eps0, kappa, false)
            .expect("law/state pairing fixed at construction")
    }

    /// Checked response for nonlinear structural solvers. J2 admission uses
    /// the full compatible tensor strain, including lateral contraction. Concrete
    /// retains the material's inclusive crushing-strain bound and history checks.
    /// Failure leaves every fiber's accepted history unchanged.
    pub fn try_respond(&self, eps0: f64, kappa: f64) -> Result<SectionState, SolidError> {
        self.response(eps0, kappa, true)
    }

    fn response(&self, eps0: f64, kappa: f64, checked: bool) -> Result<SectionState, SolidError> {
        let refused = |what: &str| SolidError::MaterialRefused { what: what.into() };
        if checked
            && (!eps0.is_finite() || !kappa.is_finite() || self.fibers.len() != self.states.len())
        {
            return Err(refused(
                "fiber strains and law/state cardinality must be valid",
            ));
        }
        let mut out = SectionState::default();
        for (f, st) in self.fibers.iter().zip(&self.states) {
            let eps = eps0 - f.y * kappa;
            if checked
                && (!eps.is_finite() || !f.y.is_finite() || !(f.area.is_finite() && f.area > 0.0))
            {
                return Err(refused(
                    "finite fiber strain and positive finite area required",
                ));
            }
            // Section convention is tension-positive; the Mander card
            // is compression-positive (tension returns zero) — flip in
            // and out for concrete.
            let (sig, tan) = match (&f.law, st) {
                (FiberLaw::Norton(l), FiberState::Norton(s)) => {
                    let response = l
                        .response(eps, s)
                        .map_err(|error| refused(&error.to_string()))?;
                    (response.stress_pa, response.tangent_pa)
                }
                (FiberLaw::RecoveringNorton(l), FiberState::RecoveringNorton(s)) => {
                    let response = l.response(eps, s)
                        .map_err(|error| refused(&error.to_string()))?;
                    (response.stress_pa, response.tangent_pa)
                }
                (FiberLaw::NortonJ2(l), FiberState::NortonJ2(s)) => {
                    let response = l.response(eps, s)
                        .map_err(|error| refused(&error.to_string()))?;
                    (response.stress_pa, response.tangent_pa)
                }
                (FiberLaw::RecoveringNortonJ2(l), FiberState::RecoveringNortonJ2(s)) => {
                    let response = l.response(eps, s)
                        .map_err(|error| refused(&error.to_string()))?;
                    (response.stress_pa, response.tangent_pa)
                }
                (FiberLaw::ThermalNorton(l), FiberState::ThermalNorton(s)) => {
                    let response = l.response(eps, s)
                        .map_err(|error| refused(&error.to_string()))?;
                    (response.stress_pa, response.tangent_pa)
                }
                (FiberLaw::ThermalNortonJ2(l), FiberState::ThermalNortonJ2(s)) => {
                    let response = l.response(eps, s)
                        .map_err(|error| refused(&error.to_string()))?;
                    (response.stress_pa, response.tangent_pa)
                }
                (FiberLaw::ThermalRecoveringNorton(l), FiberState::ThermalRecoveringNorton(s)) => {
                    let response = l.response(eps, s)
                        .map_err(|error| refused(&error.to_string()))?;
                    (response.stress_pa, response.tangent_pa)
                }
                (FiberLaw::ThermalRecoveringNortonJ2(l), FiberState::ThermalRecoveringNortonJ2(s)) => {
                    let response = l.response(eps, s)
                        .map_err(|error| refused(&error.to_string()))?;
                    (response.stress_pa, response.tangent_pa)
                }
                (FiberLaw::J2(l), FiberState::J2(s)) => {
                    let response = l.response(eps, s);
                    if checked {
                        let norm_squared =
                            fs_material::tensor::contract(&response.strain, &response.strain);
                        let bound = l.backbone().elastic.strain_limit.min(0.05);
                        if !norm_squared.is_finite() || norm_squared > bound * bound {
                            return Err(refused(
                                "J2 fiber exceeds its full-tensor small-strain domain",
                            ));
                        }
                    }
                    (response.stress_pa, response.tangent_pa)
                }
                (FiberLaw::DuctileJ2(l), FiberState::DuctileJ2(s)) => {
                    let response = l.response(eps, s)
                        .map_err(|error| refused(&error.to_string()))?;
                    (response.stress_pa, response.tangent_pa)
                }
                (FiberLaw::Steel(l), FiberState::Steel(s)) => (l.stress(eps, s), l.tangent(eps, s)),
                (FiberLaw::Voce(l), FiberState::Voce(s)) => {
                    let response = l.response(eps, s);
                    if checked {
                        let norm_squared =
                            fs_material::tensor::contract(&response.strain, &response.strain);
                        let bound = l.backbone().backbone().elastic.strain_limit.min(0.05);
                        if !norm_squared.is_finite() || norm_squared > bound * bound {
                            return Err(refused(
                                "Voce fiber exceeds its full-tensor small-strain domain",
                            ));
                        }
                    }
                    (response.stress_pa, response.tangent_pa)
                }
                (FiberLaw::TabulatedOverstress(l), FiberState::TabulatedOverstress(s)) => {
                    let response = l.response(eps, s)
                        .map_err(|error| refused(&error.to_string()))?;
                    (response.stress_pa, response.tangent_pa)
                }
                (FiberLaw::ThermalTabulated(l), FiberState::ThermalTabulated(s)) => {
                    let response = l.response(eps, s)
                        .map_err(|error| refused(&error.to_string()))?;
                    (response.stress_pa, response.tangent_pa)
                }
                (FiberLaw::ThermalChaboche(l), FiberState::ThermalChaboche(s)) => {
                    l.response(eps, s).map_err(|error| refused(&error.to_string()))?
                }
                (FiberLaw::Tabulated(l), FiberState::Tabulated(s)) => match l.response(eps, s) {
                    // The native checked reduction admits the full compatible
                    // tensor strain and the table's accumulated-strain support.
                    Ok(response) => (response.stress_pa, response.tangent_pa),
                    Err(error) if checked => return Err(refused(&error.to_string())),
                    Err(_) => (f64::NAN, f64::NAN),
                },
                (FiberLaw::Chaboche(l), FiberState::Chaboche(s)) => match l.response(eps, s) {
                    Ok(response) => (response.stress_pa, response.tangent_pa),
                    Err(error) if checked => {
                        return Err(SolidError::MaterialRefused {
                            what: error.to_string(),
                        });
                    }
                    Err(_) => (f64::NAN, f64::NAN),
                },
                (FiberLaw::ChabocheBranches(l), FiberState::ChabocheBranches(s)) => {
                    match l.response(eps, s) {
                        Ok(response) => response,
                        Err(error) if checked => return Err(refused(&error.to_string())),
                        Err(_) => (f64::NAN, f64::NAN),
                    }
                }
                (FiberLaw::Chaboche(_), FiberState::ChabocheBranches(_))
                | (FiberLaw::ChabocheBranches(_), FiberState::Chaboche(_)) => {
                    if checked {
                        return Err(refused(
                            "fiber law and backstress-history branch counts differ",
                        ));
                    }
                    (f64::NAN, f64::NAN)
                }
                (FiberLaw::Concrete(l), FiberState::Concrete(s)) => match l.response(-eps, s) {
                    Ok(response) => (-response.stress_pa, response.tangent_pa),
                    Err(error) if checked => return Err(refused(&error.to_string())),
                    Err(_) => (f64::NAN, f64::NAN),
                },
                _ => return Err(refused("fiber law and committed state must match")),
            };
            out.n += sig * f.area;
            out.m -= sig * f.y * f.area;
            let ea = tan * f.area;
            out.tangent[0][0] += ea;
            out.tangent[0][1] -= ea * f.y;
            out.tangent[1][0] -= ea * f.y;
            out.tangent[1][1] += ea * f.y * f.y;
        }
        if checked
            && (![out.n, out.m].into_iter().all(f64::is_finite)
                || out.tangent.iter().flatten().any(|v| !v.is_finite()))
        {
            return Err(refused(
                "fiber section force and tangent must be representable",
            ));
        }
        Ok(out)
    }

    /// Commit an admitted strain pair. Prefer [`Self::try_commit`] when the input
    /// has not already passed [`Self::try_respond`].
    ///
    /// # Panics
    /// Panics on a refused strain pair or invalid public fiber/history data,
    /// before changing any accepted history.
    pub fn commit(&mut self, eps0: f64, kappa: f64) {
        self.try_commit(eps0, kappa)
            .expect("fiber section commit requires an admitted response");
    }

    /// Admit every fiber and atomically commit the proposed section history.
    /// A failure in any fiber leaves all accepted states unchanged, including
    /// fibers earlier in the section that would otherwise have loaded plastically.
    pub fn try_commit(&mut self, eps0: f64, kappa: f64) -> Result<(), SolidError> {
        self.try_respond(eps0, kappa)?;
        let mut proposed = self.states.clone();
        for (f, st) in self.fibers.iter().zip(proposed.iter_mut()) {
            let eps = eps0 - f.y * kappa;
            match (&f.law, st) {
                (FiberLaw::Norton(l), FiberState::Norton(s)) => *s = l.update_state(eps, s),
                (FiberLaw::RecoveringNorton(l), FiberState::RecoveringNorton(s)) => {
                    *s = l.update_state(eps, s)
                }
                (FiberLaw::NortonJ2(l), FiberState::NortonJ2(s)) => *s = l.update_state(eps, s),
                (FiberLaw::RecoveringNortonJ2(l), FiberState::RecoveringNortonJ2(s)) => {
                    *s = l.update_state(eps, s)
                }
                (FiberLaw::ThermalNorton(l), FiberState::ThermalNorton(s)) => {
                    *s = l.update_state(eps, s)
                }
                (FiberLaw::ThermalNortonJ2(l), FiberState::ThermalNortonJ2(s)) => {
                    *s = l.update_state(eps, s)
                }
                (FiberLaw::ThermalRecoveringNorton(l), FiberState::ThermalRecoveringNorton(s)) => {
                    *s = l.update_state(eps, s)
                }
                (FiberLaw::ThermalRecoveringNortonJ2(l), FiberState::ThermalRecoveringNortonJ2(s)) => {
                    *s = l.update_state(eps, s)
                }
                (FiberLaw::J2(l), FiberState::J2(s)) => *s = l.update_state(eps, s),
                (FiberLaw::DuctileJ2(l), FiberState::DuctileJ2(s)) => {
                    *s = l.response(eps, s)
                        .map_err(|error| SolidError::MaterialRefused {
                            what: error.to_string(),
                        })?.state;
                }
                (FiberLaw::Voce(l), FiberState::Voce(s)) => *s = l.update_state(eps, s),
                (FiberLaw::Tabulated(l), FiberState::Tabulated(s)) => {
                    *s = l.response(eps, s)
                        .map_err(|error| SolidError::MaterialRefused {
                            what: error.to_string(),
                        })?.state;
                }
                (FiberLaw::TabulatedOverstress(l), FiberState::TabulatedOverstress(s)) => {
                    *s = l.response(eps, s)
                        .map_err(|error| SolidError::MaterialRefused {
                            what: error.to_string(),
                        })?.state;
                }
                (FiberLaw::ThermalTabulated(l), FiberState::ThermalTabulated(s)) => {
                    *s = l.response(eps, s)
                        .map_err(|error| SolidError::MaterialRefused {
                            what: error.to_string(),
                        })?.state;
                }
                (FiberLaw::ThermalChaboche(l), FiberState::ThermalChaboche(s)) => {
                    l.commit(eps, s).map_err(|error| SolidError::MaterialRefused {
                        what: error.to_string(),
                    })?;
                }
                (FiberLaw::Chaboche(l), FiberState::Chaboche(s)) => *s = l.update_state(eps, s),
                (FiberLaw::ChabocheBranches(l), FiberState::ChabocheBranches(s)) => {
                    l.commit(eps, s)
                }
                (FiberLaw::Steel(l), FiberState::Steel(s)) => *s = l.update_state(eps, s),
                (FiberLaw::Concrete(l), FiberState::Concrete(s)) => {
                    *s = l
                        .response(-eps, s)
                        .map_err(|error| SolidError::MaterialRefused {
                            what: error.to_string(),
                        })?
                        .state;
                }
                _ => unreachable!("law/state pairing fixed at construction"),
            }
        }
        self.states = proposed;
        Ok(())
    }
}

/// Batched section update: respond to per-section strain pairs and
/// solve each 2×2 tangent system `K·d = rhs` through fs-la's batched
/// Cholesky — the flagship's inner loop shape (state determination
/// needs exactly these solves). Returns the responses and the batch
/// solutions.
#[must_use]
pub fn update_sections_batched(
    sections: &[Section],
    strains: &[(f64, f64)],
    rhs: &[(f64, f64)],
) -> (Vec<SectionState>, Vec<(f64, f64)>) {
    assert_eq!(sections.len(), strains.len(), "one strain pair per section");
    assert_eq!(sections.len(), rhs.len(), "one rhs per section");
    let m = sections.len();
    let responses: Vec<SectionState> = sections
        .iter()
        .zip(strains)
        .map(|(s, &(e, k))| s.respond(e, k))
        .collect();
    let mut a = BatchMat::zeros(2, m);
    let mut b = BatchVec::zeros(2, m);
    for (i, r) in responses.iter().enumerate() {
        // SPD guard: fiber tangents can momentarily lose definiteness
        // on crushing branches; a tiny diagonal floor keeps the batch
        // factorization honest without masking the response itself.
        let floor = 1e-9 * (r.tangent[0][0].abs() + r.tangent[1][1].abs()).max(1e-9);
        a.set(i, 0, 0, r.tangent[0][0] + floor);
        a.set(i, 0, 1, r.tangent[0][1]);
        a.set(i, 1, 0, r.tangent[1][0]);
        a.set(i, 1, 1, r.tangent[1][1] + floor);
        b.set(i, 0, rhs[i].0);
        b.set(i, 1, rhs[i].1);
    }
    let (chol, failures) = batch_cholesky(&a);
    // Indefinite outliers fall back to their own 2×2 direct solve.
    let mut out = Vec::with_capacity(m);
    if failures.is_empty() {
        batch_cholesky_solve(&chol, &mut b);
        for i in 0..m {
            out.push((b.get(i, 0), b.get(i, 1)));
        }
    } else {
        for (i, r) in responses.iter().enumerate() {
            let det = r.tangent[0][0] * r.tangent[1][1] - r.tangent[0][1] * r.tangent[1][0];
            let d = if det.abs() > 1e-30 { det } else { 1e-30 };
            out.push((
                (r.tangent[1][1] * rhs[i].0 - r.tangent[0][1] * rhs[i].1) / d,
                (-r.tangent[1][0] * rhs[i].0 + r.tangent[0][0] * rhs[i].1) / d,
            ));
        }
    }
    (responses, out)
}

/// A reinforced-concrete rectangular section fixture: confined core +
/// unconfined cover Mander fibers plus top/bottom Menegotto–Pinto
/// steel layers. Depth `d`, width `b`, `layers` concrete fibers.
///
/// # Panics
/// On invalid material parameters (fixture constants are valid).
#[must_use]
#[allow(clippy::cast_precision_loss)]
pub fn rc_section(d: f64, b: f64, layers: usize, steel_area: f64) -> Section {
    let core = ManderConcrete::new(42e6, 0.004, 32e9, 0.015).expect("core card");
    let cover = ManderConcrete::new(30e6, 0.002, 30e9, 0.005).expect("cover card");
    let steel = MenegottoPintoSteel::new(200e9, 450e6, 0.01).expect("steel card");
    let cover_t = 0.1 * d;
    let mut fibers = Vec::new();
    let t = d / layers as f64;
    for i in 0..layers {
        let y = -0.5 * d + (i as f64 + 0.5) * t;
        let is_cover = y.abs() > 0.5 * d - cover_t;
        fibers.push(Fiber {
            y,
            area: b * t,
            law: FiberLaw::Concrete(if is_cover {
                cover.clone()
            } else {
                core.clone()
            }),
        });
    }
    for y in [-(0.5 * d - cover_t), 0.5 * d - cover_t] {
        fibers.push(Fiber {
            y,
            area: steel_area,
            law: FiberLaw::Steel(steel.clone()),
        });
    }
    Section::new(fibers)
}
