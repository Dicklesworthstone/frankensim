//! Live collisional and radiative ionization with diffuse spectral feedback.
//!
//! Each fixed-volume slab carries one original adjacent charge-stage channel.
//! Impact ionization, three-body capture, optional electron/heavy exchange and
//! held pool powers alternate symmetrically with the shared spectral diffuse
//! transfer kernel. Every radiation stage refreshes opacity and the Milne
//! spectrum at its actual populations and electron temperature. Prescribed
//! boundary beams share the same reaction/heat update. Only photons crossing
//! the two boundaries enter or leave the column; self-recycled photons
//! cancel before rounding in net species and energy updates.
//!
//! A single full/two-half controller estimates temporal error in populations,
//! both temperatures, reaction and energy integrals, and each directional face
//! exposure. The diffuse update is explicit midpoint with positivity/source
//! admission and step rejection, not an unconditionally positive or L-stable
//! method. Thermal and collisional substeps reuse the positive exponential
//! midpoint kernel. Spectral/angular/tail controls are separate from time error.
//! No state/bank cloning or species-vector allocation occurs during trial steps.
//!
//! Geometry is homogeneous plane-parallel slabs, transversely infinite, with
//! quasi-static photon propagation and immediate electron thermalization. There
//! is no finite photon inventory, radiation-pressure work,
//! scattering, stimulated emission, excited-level cascade, spatial gas flow or
//! multichannel opacity within a slab. The local/source and spectral estimates
//! do not certify the evolving solution, spatial resolution or source accuracy.

use super::{DiffuseRecombinationBeam, DiffuseRecombinationError, DiffuseRecombinationSolve};
use crate::chemistry::collisional_radiative::{
    CollisionalRadiativeBalance, CollisionalRadiativeError, CollisionalRadiativeSolve,
    CollisionalRadiativeThermal,
};
use fs_qty::{Amount, Area, Time};
use fs_thermochem::equilibrium::composition::FrozenIonizedGasStateV1;
use fs_thermochem::kinetics::atomic_ionization::radiative_recombination::AtomicRadiativeRecombinationV1;

mod solve;

/// One original atomic channel and optional thermal exchange in a mixed slab.
#[derive(Clone, Copy, Debug)]
pub struct DiffuseRecombinationEvolutionCell<'a> {
    /// Actual initial populations, separate temperatures and fixed volume.
    pub initial: &'a FrozenIonizedGasStateV1,
    /// Original reciprocal emission/absorption and atomic impact source.
    pub source: &'a AtomicRadiativeRecombinationV1,
    /// Original neutral/Coulomb exchange laws, held powers and heavy bounds.
    pub thermal: Option<&'a CollisionalRadiativeThermal>,
}
impl<'a> DiffuseRecombinationEvolutionCell<'a> {
    /// Keep the heavy temperature fixed and apply no external pool power.
    pub const fn new(
        initial: &'a FrozenIonizedGasStateV1,
        source: &'a AtomicRadiativeRecombinationV1,
    ) -> Self {
        Self {
            initial,
            source,
            thermal: None,
        }
    }
    /// Evolve both thermal pools through the existing original-source owner.
    pub const fn with_thermal(mut self, thermal: &'a CollisionalRadiativeThermal) -> Self {
        self.thermal = Some(thermal);
        self
    }
}

/// Combined temporal controls and independently bounded radiation queries.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DiffuseRecombinationEvolutionSolve {
    /// Population/thermal source bounds and temporal error/step controls.
    /// The joule tolerance applies per cell and once to the whole column.
    pub chemistry: CollisionalRadiativeSolve,
    /// Spectral, angular and omitted-tail controls for each refreshed field.
    pub radiation: DiffuseRecombinationSolve,
    /// Total refreshed fields, including rejected trials and predictors.
    pub maximum_radiation_evaluations: usize,
    /// Total cell/source/angular work across every refreshed field.
    pub maximum_radiation_work: usize,
}
impl DiffuseRecombinationEvolutionSolve {
    /// Explicit chemistry controls with the existing radiation defaults.
    pub fn new(chemistry: CollisionalRadiativeSolve) -> Self {
        let mut radiation = DiffuseRecombinationSolve::default();
        radiation.relative_tolerance = radiation
            .relative_tolerance
            .min(0.05 * chemistry.relative_tolerance);
        radiation.source_quadrature.relative_tolerance = radiation
            .source_quadrature
            .relative_tolerance
            .min(0.05 * radiation.relative_tolerance);
        Self {
            chemistry,
            radiation,
            maximum_radiation_evaluations: 20_000,
            maximum_radiation_work: 1_000_000_000,
        }
    }
}

/// Time-integrated directional traffic through one common-area interface.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DiffuseRecombinationExposure {
    /// Number of photons moving toward larger slab indices.
    pub right_photons: f64,
    /// Their integrated energy [J].
    pub right_energy_j: f64,
    /// Number of photons moving toward smaller slab indices.
    pub left_photons: f64,
    /// Their integrated energy [J].
    pub left_energy_j: f64,
}

/// Accepted local endpoint with original immutable inputs and radiation accounts.
#[derive(Clone, Debug, PartialEq)]
pub struct DiffuseRecombinationEvolutionCellResult {
    initial: FrozenIonizedGasStateV1,
    state: FrozenIonizedGasStateV1,
    source: AtomicRadiativeRecombinationV1,
    thermal: Option<CollisionalRadiativeThermal>,
    extent: Amount,
    balance: CollisionalRadiativeBalance,
    radiation: DiffuseRecombinationIntegratedRates,
}
impl DiffuseRecombinationEvolutionCellResult {
    /// Original bank with actual accepted populations and temperatures.
    pub const fn state(&self) -> &FrozenIonizedGasStateV1 {
        &self.state
    }
    /// Unmodified original state.
    pub const fn initial_state(&self) -> &FrozenIonizedGasStateV1 {
        &self.initial
    }
    /// Original reciprocal atomic channel, without a fitted replacement.
    pub const fn source(&self) -> &AtomicRadiativeRecombinationV1 {
        &self.source
    }
    /// Retained original optional exchange sources and held pool powers.
    pub const fn thermal_inputs(&self) -> Option<&CollisionalRadiativeThermal> {
        self.thermal.as_ref()
    }
    /// Net pair creation by collisions and radiation [mol].
    pub const fn reaction_extent(&self) -> Amount {
        self.extent
    }
    /// Original local reaction/energy accounts. Here emitted photons enter the
    /// transported field; only the column boundary exposure counts as escape.
    pub const fn balance(&self) -> CollisionalRadiativeBalance {
        self.balance
    }
    /// Stable net radiation and resolved self/nonlocal recycling integrals.
    pub const fn radiation(&self) -> DiffuseRecombinationIntegratedRates {
        self.radiation
    }
}

/// Integrated transfer quantities, with net accounts evaluated independently.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DiffuseRecombinationIntegratedRates {
    /// Photons absorbed in their emitting slab.
    pub self_absorbed_photons: f64,
    /// Photons absorbed after emission by another slab.
    pub nonlocal_absorbed_photons: f64,
    /// Photons absorbed from the prescribed external boundary beams.
    pub incident_absorbed_photons: f64,
    /// Their absorbed energy [J].
    pub incident_absorbed_energy_j: f64,
    /// Radiation-only net pair creation [mol]; self-recycling cancels first.
    pub net_ionization_extent_mol: f64,
    /// Radiation-only net electron sensible-energy input [J].
    pub net_electron_energy_j: f64,
    /// Radiation-only net gas-energy input [J].
    pub net_gas_energy_j: f64,
}

/// Actual whole-column conservation and bounded integration/query work.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DiffuseRecombinationEvolutionBalance {
    /// Integrated prescribed photons entering both boundaries.
    pub incident_photons: f64,
    /// Their integrated input energy [J].
    pub incident_energy_j: f64,
    /// Integrated outgoing photons, including transmitted boundary beams.
    pub escaped_photons: f64,
    /// Their integrated energy [J].
    pub escaped_energy_j: f64,
    /// Sum of independently integrated net radiation pair creation [mol].
    pub net_radiative_extent_mol: f64,
    /// Sum of direct net gas radiation inputs [J].
    pub net_radiation_energy_j: f64,
    /// Actual endpoint chemical + electron + heavy energy change [J].
    pub gas_energy_change_j: f64,
    /// Held external electron and heavy work [J].
    pub external_work_j: f64,
    /// Net radiative pair creation times N_A + outgoing - incoming photons.
    pub photon_residual: f64,
    /// Actual gas energy change + outgoing - incoming radiation - external work [J].
    /// The boundary comparison includes `boundary_energy_roundoff_j`.
    pub energy_residual_j: f64,
    /// Actual gas energy change - direct net radiation input - external work [J].
    /// This independent storage/source check keeps the absolute chemistry tolerance.
    pub gas_source_residual_j: f64,
    /// Additional allowance [J] only for comparisons of rounded interface energies.
    /// This is 128 epsilon times total incoming plus outgoing energy; it does not
    /// increase the chemistry tolerance or certify a physical transport error.
    pub boundary_energy_roundoff_j: f64,
    /// Complete accepted steps; two-half-step endpoints are retained.
    pub accepted_steps: usize,
    /// Rejected temporal or source/positivity trials.
    pub rejected_steps: usize,
    /// Largest admitted dimensionless temporal error estimate.
    pub maximum_normalized_local_error: f64,
    /// Refreshed fields including predictor and rejected-trial work.
    pub radiation_evaluations: usize,
    /// Total cell/source/angular work, including initial preparation.
    pub radiation_work: usize,
    /// Largest spectral estimate in any successfully evaluated trial field.
    pub maximum_normalized_spectral_estimate: f64,
    /// Largest angular estimate in any successfully evaluated trial field.
    pub maximum_normalized_angular_estimate: f64,
    /// Largest conditional omitted-tail estimate in any evaluated trial field.
    pub maximum_normalized_tail_bound: f64,
}

/// Complete accepted evolution; failed solves publish no partial endpoint.
#[derive(Clone, Debug, PartialEq)]
pub struct DiffuseRecombinationEvolutionResult {
    cells: Vec<DiffuseRecombinationEvolutionCellResult>,
    interfaces: Vec<DiffuseRecombinationExposure>,
    balance: DiffuseRecombinationEvolutionBalance,
    area: Area,
    elapsed: Time,
    controls: DiffuseRecombinationEvolutionSolve,
    incident_beams: Vec<DiffuseRecombinationBeam>,
}
impl DiffuseRecombinationEvolutionResult {
    /// Ordered source-retained gas endpoints.
    pub fn cells(&self) -> &[DiffuseRecombinationEvolutionCellResult] {
        &self.cells
    }
    /// N+1 directional integrated exposures, including illuminated boundaries.
    pub fn interfaces(&self) -> &[DiffuseRecombinationExposure] {
        &self.interfaces
    }
    /// Whole-column physical accounts and explicit work.
    pub const fn balance(&self) -> DiffuseRecombinationEvolutionBalance {
        self.balance
    }
    /// Original common slab face area [m²].
    pub const fn area(&self) -> Area {
        self.area
    }
    /// Complete requested duration [s].
    pub const fn elapsed_time(&self) -> Time {
        self.elapsed
    }
    /// Applied source, temporal and radiation controls.
    pub const fn controls(&self) -> DiffuseRecombinationEvolutionSolve {
        self.controls
    }
    /// Original prescribed boundary beams, held throughout this exposure.
    pub fn incident_beams(&self) -> &[DiffuseRecombinationBeam] {
        &self.incident_beams
    }
}

/// Source, numerical or conservation refusal without caller mutation.
#[derive(Clone, Debug, PartialEq)]
pub enum DiffuseRecombinationEvolutionError {
    /// Invalid geometry, duration or work controls.
    Invalid(&'static str),
    /// Indexed original local source or state admission refused.
    Cell {
        index: usize,
        error: CollisionalRadiativeError,
    },
    /// Shared spectral/angular radiation field refused.
    Radiation(DiffuseRecombinationError),
    /// Accepted plus rejected steps or minimum representable step exhausted.
    StepBudget,
    /// Total field evaluations or source/characteristic work exhausted.
    RadiationBudget,
    /// Required finite arithmetic or bounded allocation refused.
    Arithmetic(&'static str),
    /// Integrated photon conservation failed.
    PhotonBalance { residual_photons: f64 },
    /// Actual gas, boundary radiation and external work did not close.
    EnergyBalance { residual_j: f64 },
    /// Cancellation before final publication.
    Cancelled,
}
impl core::fmt::Display for DiffuseRecombinationEvolutionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "evolving diffuse recombination radiation: {self:?}")
    }
}
impl std::error::Error for DiffuseRecombinationEvolutionError {}

/// Evolve original atomic populations and thermal pools with diffuse recycling.
///
/// Neutral cells may acquire their first physical carriers from a luminous
/// neighbor. The emitted spectrum is recomputed from the changing emitter, and
/// absorption changes the receiver's next opacity and electron temperature.
/// Empty columns and completely inactive dark cells return identity endpoints.
pub fn advance_diffuse_recombination(
    cells: &[DiffuseRecombinationEvolutionCell<'_>],
    area: Area,
    duration: Time,
    controls: DiffuseRecombinationEvolutionSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<DiffuseRecombinationEvolutionResult, DiffuseRecombinationEvolutionError> {
    solve::advance(cells, area, &[], duration, controls, cancelled)
}

/// Evolve ionization, shielding and thermal pools with external illumination.
///
/// Every radiation stage transports the prescribed discrete boundary beams
/// through its current populations together with spontaneous diffuse emission.
/// Absorption may create the first electrons without an artificial seed, and
/// all incident energy is charged to the incoming photon boundary accounts.
/// The beams are held over the interval; no drive time profile or photon flight
/// delay is inferred. Empty columns retain their complete transmitted exposure.
pub fn advance_diffuse_recombination_with_illumination(
    cells: &[DiffuseRecombinationEvolutionCell<'_>],
    area: Area,
    incident_beams: &[DiffuseRecombinationBeam],
    duration: Time,
    controls: DiffuseRecombinationEvolutionSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<DiffuseRecombinationEvolutionResult, DiffuseRecombinationEvolutionError> {
    solve::advance(cells, area, incident_beams, duration, controls, cancelled)
}
