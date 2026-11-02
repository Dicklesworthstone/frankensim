//! Simultaneous adhesive sliding and spin on shared physical bodies.
//!
//! Every junction keeps its actual solid-contact radius and uniform static and
//! kinetic shear strengths. The unchanged full radius-scaled compliance couples
//! all three coordinates of every patch. The common circular-support engine
//! admits positive-semidefinite responses, redistributes dependent static loads,
//! and checks every final patch against its own uniform-disk support function.
//! Its torsional bound is `2*tau*A*a/3`; no pressure law or normal load is used.
//!
//! Block sweeps are numerical iterations, not physical timesteps. Work and loss
//! use the final common motion against the complete original matrix. A moving
//! junction retains kinetic capacity while a stationary junction retains its
//! static set. The mechanical owner supplies phases, physical ordering and
//! stopping events. An unresolved static feasibility decision refuses instead
//! of admitting a lower-strength kinetic alternative. Heat allocation, normal
//! evolution, partial slip and interface peeling remain with their owners.

use super::{AdhesiveWrenchContact, AdhesiveWrenchPhase};
use crate::sliding_patch::contact::multiple::{CircularContactInput, PreparedCircularContacts};

/// Shared full-response, sweep, evaluation and dimensional closure budgets.
pub use crate::sliding_patch::contact::multiple::MultipleHertzConfig as MultipleAdhesiveConfig;
/// Complete-candidate cancellation, unresolved feasibility and physical refusals.
pub use crate::sliding_patch::contact::multiple::MultipleHertzError as MultipleAdhesiveError;
/// All original patch wrenches and their common full-matrix physical work.
pub use crate::sliding_patch::contact::multiple::MultipleHertzStep as MultipleAdhesiveStep;

/// One original adhesive junction and its released physical motion.
#[derive(Clone, Copy, Debug)]
pub struct AdhesiveContactInput<'a> {
    /// Original real radius and independently supplied static/kinetic strengths.
    pub law: &'a AdhesiveWrenchContact,
    /// Released tangential displacement in this contact's original frame [m].
    pub free_displacement_m: [f64; 2],
    /// Released relative rotation about this contact's original normal [rad].
    pub free_spin_angle_rad: f64,
}

/// Prepared physical response; source laws are borrowed for each trial.
/// Zero-area detached junctions have no wrench and must be omitted rather than
/// assigned an invented radius to enter the scaled mechanical system.
#[derive(Clone, Debug, PartialEq)]
pub struct PreparedAdhesiveContacts {
    response: PreparedCircularContacts,
}

impl PreparedAdhesiveContacts {
    /// Admit the complete row-major `(3*count)^2` compliance [m/N] in coordinates
    /// `[u_x,u_y,a*spin]` and dual efforts `[F_x,F_y,M/a]` for each real radius.
    /// Dependent rows are allowed; each local block must be positive definite
    /// or exactly zero. The original physical matrix is never regularized.
    pub fn new(
        count: usize,
        scaled_compliance_m_per_n: Vec<f64>,
        config: MultipleAdhesiveConfig,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MultipleAdhesiveError> {
        Ok(Self {
            response: PreparedCircularContacts::new(
                count,
                scaled_compliance_m_per_n,
                config,
                cancelled,
            )?,
        })
    }

    /// Original number of attached junctions, in the mechanical owner's order.
    #[must_use]
    pub const fn contact_count(&self) -> usize {
        self.response.contact_count()
    }

    /// Original complete physical response; no factor or diagonal substitute.
    #[must_use]
    pub fn scaled_compliance_m_per_n(&self) -> &[f64] {
        self.response.scaled_compliance_m_per_n()
    }

    /// Original physical acceptance and bounded numerical work controls.
    #[must_use]
    pub const fn config(&self) -> MultipleAdhesiveConfig {
        self.response.config()
    }

    /// Resolve stationary-admissible junctions from one common released state.
    /// Static redistribution precedes kinetic selection on dependent responses.
    /// No individual patch is accepted before the final full-matrix audit.
    pub fn advance(
        &self,
        sites: &[AdhesiveContactInput<'_>],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MultipleAdhesiveStep, MultipleAdhesiveError> {
        let sites = self.inputs(sites, &mut cancelled)?;
        self.response.advance(&sites, cancelled)
    }

    /// Resolve the caller's original moving/static phases together. A larger
    /// static strength cannot prematurely capture a junction identified as
    /// moving. Zero mean motion still requires native endpoint/event checks.
    pub fn advance_phases(
        &self,
        sites: &[AdhesiveContactInput<'_>],
        phases: &[AdhesiveWrenchPhase],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MultipleAdhesiveStep, MultipleAdhesiveError> {
        let sites = self.inputs(sites, &mut cancelled)?;
        self.response.advance_phases(&sites, phases, cancelled)
    }

    fn inputs<'a>(
        &self,
        sites: &[AdhesiveContactInput<'a>],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<CircularContactInput<'a, AdhesiveWrenchContact>>, MultipleAdhesiveError> {
        if cancelled() {
            return Err(MultipleAdhesiveError::Cancelled);
        }
        if sites.len() != self.contact_count() {
            return Err(MultipleAdhesiveError::Invalid("contact input count"));
        }
        sites
            .iter()
            .map(|site| {
                if cancelled() {
                    return Err(MultipleAdhesiveError::Cancelled);
                }
                Ok(CircularContactInput {
                    law: site.law,
                    free_displacement_m: site.free_displacement_m,
                    free_spin_angle_rad: site.free_spin_angle_rad,
                })
            })
            .collect()
    }
}
