//! Spectral, bidirectional reabsorption of spontaneous recombination radiation.
//!
//! Ordered, individually mixed plane-parallel slabs have a common face area
//! and vacuum boundaries. Their ORIGINAL Milne spectra emit isotropically;
//! original bound-free cross sections absorb at every photon energy. Exact
//! uniform-slab characteristics include self-absorption and illumination from
//! both directions. Positive angular quadrature and adaptive spectral
//! quadrature retain frequency-dependent opacity, never a mean emitted photon.
//!
//! This is a frozen-state, quasi-static radiation query. Its per-cell reaction
//! and electron-energy rates are useful inputs to coupled evolution; this
//! owner does not itself change populations or temperatures. No incoming beam,
//! scattering, excited-level cascade, stimulated emission, finite propagation
//! delay, line radiation or radiation-pressure work is inferred. Geometry is
//! transversely infinite; the area scales inventories and face throughput.

use fs_qty::Area;
use fs_thermochem::electron_collision::MaxwellianRateControlsV1;
use fs_thermochem::equilibrium::composition::FrozenIonizedGasStateV1;
use fs_thermochem::kinetics::atomic_ionization::photoionization::PhotoionizationErrorV1;
use fs_thermochem::kinetics::atomic_ionization::radiative_recombination::{
    AtomicRadiativeRecombinationV1, RadiativeRecombinationErrorV1,
};

mod solve;

/// One homogeneous slab with a retained atomic channel and actual populations.
#[derive(Clone, Copy, Debug)]
pub struct DiffuseRecombinationCell<'a> {
    /// Original thermodynamic bank, amounts, volume and temperatures.
    pub state: &'a FrozenIonizedGasStateV1,
    /// Original reciprocal emission and absorption source.
    pub source: &'a AtomicRadiativeRecombinationV1,
}

/// Numerical controls for all transported number and energy moments.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DiffuseRecombinationSolve {
    /// Relative spectral/angular estimates and conditional omitted-tail bound.
    pub relative_tolerance: f64,
    /// Absolute allowance for each photon-rate observable [s⁻¹].
    pub absolute_photon_rate_tolerance_per_s: f64,
    /// Absolute allowance for each energy-rate observable [W].
    pub absolute_power_tolerance_w: f64,
    /// Original frozen-state quasineutral admission policy.
    pub relative_charge_tolerance: f64,
    /// Maximum admitted slab count.
    pub maximum_cells: usize,
    /// Total attempted energy panels, including all refinements/restarts.
    pub maximum_spectral_panels: usize,
    /// Maximum composite angular subdivisions under mu=t⁴ on 0<t<1.
    pub maximum_angular_subdivisions: usize,
    /// Total cell/source or angular cell-characteristic evaluations.
    pub maximum_work: usize,
    /// Existing original-source total-moment quadrature controls.
    pub source_quadrature: MaxwellianRateControlsV1,
}
impl Default for DiffuseRecombinationSolve {
    fn default() -> Self {
        Self {
            relative_tolerance: 1e-6,
            absolute_photon_rate_tolerance_per_s: 1.0,
            absolute_power_tolerance_w: 1e-15,
            relative_charge_tolerance: 1e-12,
            maximum_cells: 256,
            maximum_spectral_panels: 32_768,
            maximum_angular_subdivisions: 64,
            maximum_work: 100_000_000,
            source_quadrature: MaxwellianRateControlsV1 {
                relative_tolerance: 1e-10,
                maximum_subdivisions: 4096,
            },
        }
    }
}

/// Directional photon throughput across an interface, integrated over angles.
/// These are total rates through the common area, not angle-integrated intensity.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DiffuseRecombinationFlux {
    /// Photons traveling toward larger cell indices [s⁻¹].
    pub right_photons_per_s: f64,
    /// Corresponding power [W].
    pub right_power_w: f64,
    /// Photons traveling toward smaller cell indices [s⁻¹].
    pub left_photons_per_s: f64,
    /// Corresponding power [W].
    pub left_power_w: f64,
}

/// Independent source, absorption and gas-energy rates in one frozen cell.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DiffuseRecombinationCellBalance {
    /// Gross spontaneous photons from the represented spectrum [s⁻¹].
    pub emitted_photons_per_s: f64,
    /// Gross emitted chemical plus captured-electron energy [W].
    pub emitted_power_w: f64,
    /// Captured electron kinetic-energy loss [W], not 3kT/2 per event.
    pub captured_electron_power_w: f64,
    /// Absorption of this same cell's emitted photons [s⁻¹].
    pub self_absorbed_photons_per_s: f64,
    /// Absorption of photons emitted in other cells [s⁻¹].
    pub nonlocal_absorbed_photons_per_s: f64,
    /// Total absorbed photon energy [W].
    pub absorbed_power_w: f64,
    /// Photoelectron kinetic energy deposited after local threshold cost [W].
    pub photoelectron_power_w: f64,
    /// (Absorptions - captures)/N_A [mol/s]; no state update is implied.
    pub net_ionization_rate_mol_per_s: f64,
    /// Photoelectron heating minus captured-electron cooling [W].
    pub net_electron_power_w: f64,
    /// Original chemical gap times net ionization [W].
    pub chemical_storage_rate_w: f64,
    /// Gas-energy source: absorption minus emission [W].
    pub net_gas_power_w: f64,
    /// Source plus both incoming interfaces minus absorption and both exits.
    pub photon_residual_per_s: f64,
    /// Independently transported energy version of that balance [W].
    pub energy_residual_w: f64,
    /// Chemical storage plus electron energy minus actual gas power [W].
    pub gas_energy_residual_w: f64,
}

/// A source-retained frozen cell and its complete radiative rates.
#[derive(Clone, Debug, PartialEq)]
pub struct DiffuseRecombinationCellResult {
    state: FrozenIonizedGasStateV1,
    source: AtomicRadiativeRecombinationV1,
    balance: DiffuseRecombinationCellBalance,
    emission_cutoff_j: f64,
}
impl DiffuseRecombinationCellResult {
    /// Original unmodified state.
    pub const fn state(&self) -> &FrozenIonizedGasStateV1 {
        &self.state
    }
    /// Original reciprocal source.
    pub const fn source(&self) -> &AtomicRadiativeRecombinationV1 {
        &self.source
    }
    /// Gas and radiation rates.
    pub const fn balance(&self) -> DiffuseRecombinationCellBalance {
        self.balance
    }
    /// Upper energy represented for THIS emitter; absorption uses common support.
    pub const fn emission_cutoff_j(&self) -> f64 {
        self.emission_cutoff_j
    }
}

/// Whole-column radiation conservation and explicit numerical limitations.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DiffuseRecombinationBalance {
    /// Gross emission and reabsorption are separate, including self-absorption.
    pub emitted_photons_per_s: f64,
    /// Gross emitted photon power [W].
    pub emitted_power_w: f64,
    /// Total photons absorbed across all cells [s⁻¹].
    pub absorbed_photons_per_s: f64,
    /// Total absorbed photon power [W].
    pub absorbed_power_w: f64,
    /// Escapes through both vacuum boundaries [s⁻¹].
    pub escaped_photons_per_s: f64,
    /// Power escaping through both boundaries [W].
    pub escaped_power_w: f64,
    /// Emission minus absorption and escape [s⁻¹].
    pub photon_residual_per_s: f64,
    /// Emission minus absorption and escape [W].
    pub energy_residual_w: f64,
    /// Conditional upper bound on all omitted emitted photon rates [s⁻¹].
    pub omitted_photon_rate_bound_per_s: f64,
    /// Conditional upper bound on all omitted emitted photon power [W].
    pub omitted_power_bound_w: f64,
    /// Largest spectral Gauss 8/16 estimate divided by its requested allowance.
    pub maximum_normalized_spectral_estimate: f64,
    /// Largest angular Gauss 8/16 estimate divided by its requested allowance.
    pub maximum_normalized_angular_estimate: f64,
    /// Largest conditional tail bound divided by a transport observable's allowance.
    pub maximum_normalized_tail_bound: f64,
    /// Total attempted energy panels.
    pub spectral_panels: usize,
    /// Accepted angular subdivisions under the grazing-angle transform.
    pub angular_subdivisions: usize,
    /// Actual bounded work counter.
    pub work: usize,
}

/// Steady diffuse radiation field with unmodified source cells.
#[derive(Clone, Debug, PartialEq)]
pub struct DiffuseRecombinationResult {
    cells: Vec<DiffuseRecombinationCellResult>,
    interfaces: Vec<DiffuseRecombinationFlux>,
    balance: DiffuseRecombinationBalance,
    area: Area,
    controls: DiffuseRecombinationSolve,
}
impl DiffuseRecombinationResult {
    /// Original cells and their radiative reaction/energy rates.
    pub fn cells(&self) -> &[DiffuseRecombinationCellResult] {
        &self.cells
    }
    /// N+1 interfaces, including the two vacuum boundaries.
    pub fn interfaces(&self) -> &[DiffuseRecombinationFlux] {
        &self.interfaces
    }
    /// Whole-column accounts and numerical estimates.
    pub const fn balance(&self) -> DiffuseRecombinationBalance {
        self.balance
    }
    /// Original plane-parallel face area.
    pub const fn area(&self) -> Area {
        self.area
    }
    /// Retained numerical controls.
    pub const fn controls(&self) -> DiffuseRecombinationSolve {
        self.controls
    }
}

/// Refusal before publishing a complete radiation field or any partial state.
#[derive(Clone, Debug, PartialEq)]
pub enum DiffuseRecombinationError {
    /// Invalid geometry, controls, state/source identity or arithmetic input.
    Invalid(&'static str),
    /// Original emission source refused in the indexed cell.
    Source {
        index: usize,
        error: RadiativeRecombinationErrorV1,
    },
    /// Original absorption source refused at a transported photon energy.
    Absorption {
        index: usize,
        error: PhotoionizationErrorV1,
    },
    /// Cell or quadrature allocation exceeded its admitted bounds.
    Allocation,
    /// Total cell/source/characteristic work exhausted.
    WorkBudget,
    /// Total attempted spectral panels exhausted.
    SpectralBudget,
    /// Angular refinement exhausted before its transfer estimates converged.
    AngularBudget,
    /// Original shared source support cannot resolve the omitted spectrum.
    UnresolvedTail { maximum_normalized_bound: f64 },
    /// Photon or independently evaluated energy account did not close.
    Balance {
        photon_residual_per_s: f64,
        energy_residual_w: f64,
    },
    /// A required finite quantity cannot be represented.
    Arithmetic(&'static str),
    /// Cancellation before publication.
    Cancelled,
}
impl core::fmt::Display for DiffuseRecombinationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "diffuse recombination radiation: {self:?}")
    }
}
impl std::error::Error for DiffuseRecombinationError {}

/// Evaluate isotropic spontaneous emission, self/nonlocal absorption and escape.
///
/// Each cell uses its own original spectrum and threshold. Empty columns and
/// absent emitters return exact darkness. The numerical field is a frozen,
/// quasi-static input to evolution, not a finite photon inventory or a claim
/// that the gas stays fixed under the returned nonzero reaction/energy rates.
pub fn evaluate_diffuse_recombination(
    cells: &[DiffuseRecombinationCell<'_>],
    area: Area,
    controls: DiffuseRecombinationSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<DiffuseRecombinationResult, DiffuseRecombinationError> {
    solve::evaluate(cells, area, controls, cancelled)
}
