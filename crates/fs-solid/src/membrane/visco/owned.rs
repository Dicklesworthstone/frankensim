//! Owned native motion with immutable material bindings and shared histories.
//!
//! An accepted native motion can outlive its caller's patch and law storage.
//! Borrowed views share its original accepted buffers; native advancement still
//! owns all geometry, inertia, local material solves and work gates. Only an
//! accepted native result can supply new histories. No public snapshot restore
//! can attach arbitrary histories to a different source or constitutive law.
//!
//! A kinetics-only event changes finite positive relaxation times and nothing
//! else. Since neither the energy nor natural reference changes, every accepted
//! stress, pressure, energy ledger, metric, displacement and velocity survives
//! exactly. No transverse solve, work, heat, impulse or elapsed time is created.

/// Owned prescribed-reference proposals with the complete native work accounts.
pub mod reference_change;
pub use reference_change::OwnedFiniteViscoMembraneMotionReferenceChange;

/// Per-triangle source temperatures with retained natural-strain history.
pub mod temperature;

use std::sync::Arc;

use fs_material::sheet::SheetError;
use fs_material::visco_hyper::ViscoHyperelastic;
use fs_material::visco_hyper::membrane::ViscoMembraneState;

use super::dynamics::{FiniteViscoMembraneMotion, FiniteViscoMotionControls, InertialEnergy};
use super::{
    FiniteViscoMembrane, FiniteViscoMembraneEnergy, FiniteViscoMembraneError, MembranePatch,
    MembraneSolution, MembraneSolveSettings, ViscoelasticLawField, invalid, material_error, poll,
};
use crate::membrane::inertia::MembraneMassImage;

#[derive(Debug)]
enum OwnedLaws {
    Uniform(ViscoHyperelastic),
    PerElement(Vec<ViscoHyperelastic>),
}

impl OwnedLaws {
    fn from_native(laws: ViscoelasticLawField<'_>) -> Self {
        match laws {
            ViscoelasticLawField::Uniform(law) => Self::Uniform(law.clone()),
            ViscoelasticLawField::PerElement(laws) => Self::PerElement(laws.to_vec()),
        }
    }

    fn as_native(&self) -> ViscoelasticLawField<'_> {
        match self {
            Self::Uniform(law) => ViscoelasticLawField::Uniform(law),
            Self::PerElement(laws) => ViscoelasticLawField::PerElement(laws),
        }
    }
}

// Private accepted buffers contain no patch/law borrows. Capturing an already
// admitted native body does no computation and never constructs virgin history.
#[derive(Clone, Debug)]
struct AcceptedMembrane {
    displacement_m: Arc<[[f64; 3]]>,
    states: Arc<[ViscoMembraneState]>,
    first_piola_pa: Arc<[[f64; 9]]>,
    pressure_pa: Arc<[Option<f64>]>,
    resultants_n_m: Arc<[[f64; 3]]>,
    stored_energy_j: f64,
    time_s: f64,
}

impl AcceptedMembrane {
    fn from_native(membrane: &FiniteViscoMembrane<'_>) -> Self {
        Self {
            displacement_m: Arc::clone(&membrane.displacement_m),
            states: Arc::clone(&membrane.states),
            first_piola_pa: Arc::clone(&membrane.first_piola_pa),
            pressure_pa: Arc::clone(&membrane.pressure_pa),
            resultants_n_m: Arc::clone(&membrane.resultants_n_m),
            stored_energy_j: membrane.stored_energy_j,
            time_s: membrane.time_s,
        }
    }

    fn bind<'a>(
        &self,
        patch: &'a MembranePatch,
        material: ViscoelasticLawField<'a>,
    ) -> FiniteViscoMembrane<'a> {
        FiniteViscoMembrane {
            patch,
            material,
            displacement_m: Arc::clone(&self.displacement_m),
            states: Arc::clone(&self.states),
            first_piola_pa: Arc::clone(&self.first_piola_pa),
            pressure_pa: Arc::clone(&self.pressure_pa),
            resultants_n_m: Arc::clone(&self.resultants_n_m),
            stored_energy_j: self.stored_energy_j,
            time_s: self.time_s,
        }
    }
}

/// Accepted native motion independent of the caller's patch and law lifetimes.
/// Clones and borrowed views share every immutable history/velocity buffer.
#[derive(Clone, Debug)]
pub struct OwnedFiniteViscoMembraneMotion {
    patch: Arc<MembranePatch>,
    laws: Arc<OwnedLaws>,
    accepted: AcceptedMembrane,
    velocity_mps: Arc<[[f64; 3]]>,
    mass_image: MembraneMassImage,
}

/// The complete native physical interval with its owned accepted successor.
/// Physical material heat and numerical kinetic loss remain separate accounts.
#[derive(Debug)]
pub struct OwnedFiniteViscoMembraneMotionStep {
    pub motion: OwnedFiniteViscoMembraneMotion,
    pub equilibrium: MembraneSolution,
    pub material_energy: FiniteViscoMembraneEnergy,
    pub element_dissipated_j: Vec<f64>,
    pub inertia: InertialEnergy,
    pub balance_error_j: f64,
}

impl OwnedFiniteViscoMembraneMotion {
    /// Capture a validated native motion. The original immutable patch and law
    /// field are copied once; accepted history, responses and velocity are shared.
    /// Uniform laws remain uniform rather than being duplicated per triangle.
    pub fn from_native(motion: &FiniteViscoMembraneMotion<'_>) -> Self {
        Self {
            patch: Arc::new(motion.membrane.patch.clone()),
            laws: Arc::new(OwnedLaws::from_native(motion.membrane.material)),
            accepted: AcceptedMembrane::from_native(&motion.membrane),
            velocity_mps: Arc::clone(&motion.velocity_mps),
            mass_image: motion.mass_image,
        }
    }

    /// O(1) native view. Its bindings borrow this owner while accepted buffers
    /// are shared; constructing the view performs no mass or material solve.
    pub fn as_native(&self) -> FiniteViscoMembraneMotion<'_> {
        FiniteViscoMembraneMotion {
            membrane: self.accepted.bind(&self.patch, self.native_laws()),
            velocity_mps: Arc::clone(&self.velocity_mps),
            mass_image: self.mass_image,
        }
    }

    /// The same original reference geometry, source bundles and graded inventory.
    pub fn patch(&self) -> &MembranePatch {
        &self.patch
    }

    /// Immutable current kinetics, without source resolution or point cloning.
    pub fn native_laws(&self) -> ViscoelasticLawField<'_> {
        self.laws.as_native()
    }

    /// Replace only branch relaxation times in original triangle order.
    ///
    /// Supply exactly one law per triangle, including for a previously uniform
    /// field. Equilibrium energy/domain, ordered branch energies, natural
    /// reference and local solver controls must equal the accepted law. Every
    /// new relaxation time must be finite and positive. No energy, stress or
    /// pressure is recomputed: a kinetics-only event cannot change them.
    ///
    /// The entire candidate is private until all laws and the final cancellation
    /// check pass. Subsequent native advancement uses these new times with the
    /// retained metrics. Temperature adapters must shift their original
    /// calibrated spectrum, never repeatedly shift this current endpoint.
    pub fn with_kinetics(
        &self,
        target_laws: &[ViscoHyperelastic],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, FiniteViscoMembraneError> {
        poll(&mut cancelled)?;
        if target_laws.len() != self.accepted.states.len() {
            return Err(invalid(
                "one kinetic target law per original membrane triangle required",
            ));
        }
        let current = self.native_laws();
        let mut proposed = Vec::with_capacity(target_laws.len());
        for (element, target) in target_laws.iter().enumerate() {
            poll(&mut cancelled)?;
            let previous = current.at(element);
            if target.branches.len() != previous.branches.len() {
                return Err(material_error(
                    element,
                    SheetError::Invalid("kinetic events retain the original ordered branches"),
                ));
            }
            // Clone only the already admitted, bounded law. An incompatible
            // target cannot force copying an arbitrary new energy payload.
            let mut admitted = previous.clone();
            for (branch, target) in admitted.branches.iter_mut().zip(&target.branches) {
                if !target.tau_s.is_finite() || target.tau_s <= 0.0 {
                    return Err(material_error(
                        element,
                        SheetError::Invalid("finite positive target relaxation times required"),
                    ));
                }
                branch.tau_s = target.tau_s;
            }
            // Whole-law equality includes private reference/solver fields and
            // any future native fields; the admitted clone changed only times.
            if admitted != *target {
                return Err(material_error(
                    element,
                    SheetError::Invalid(
                        "kinetic events cannot change energy, domain, natural reference or solver controls",
                    ),
                ));
            }
            proposed.push(admitted);
        }
        poll(&mut cancelled)?;
        Ok(Self {
            laws: Arc::new(OwnedLaws::PerElement(proposed)),
            ..self.clone()
        })
    }

    /// Advance through the original native momentum/material solver. This owner
    /// adds no constitutive solve, time integrator or separate force recovery.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_pressure(
        &self,
        dt_s: f64,
        dead_load_n: &[[f64; 3]],
        pressure_pa: f64,
        prescribed: &[(usize, f64)],
        settings: MembraneSolveSettings,
        controls: FiniteViscoMotionControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<OwnedFiniteViscoMembraneMotionStep, FiniteViscoMembraneError> {
        let native = self.as_native();
        let step = native.advance_with_pressure(
            dt_s,
            dead_load_n,
            pressure_pa,
            prescribed,
            settings,
            controls,
            &mut cancelled,
        )?;
        let motion = Self {
            patch: Arc::clone(&self.patch),
            laws: Arc::clone(&self.laws),
            accepted: AcceptedMembrane::from_native(&step.motion.membrane),
            velocity_mps: Arc::clone(&step.motion.velocity_mps),
            mass_image: step.motion.mass_image,
        };
        poll(&mut cancelled)?;
        Ok(OwnedFiniteViscoMembraneMotionStep {
            motion,
            equilibrium: step.equilibrium,
            material_energy: step.material_energy,
            element_dissipated_j: step.element_dissipated_j,
            inertia: step.inertia,
            balance_error_j: step.balance_error_j,
        })
    }
}
