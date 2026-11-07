//! Finite photon transport coupled to local collisional and radiative chemistry.
//!
//! Ordered, individually mixed cells share a fixed beam area. Each original
//! reciprocal atomic source supplies photoionization, electron impact, its
//! three-body reverse and spontaneous radiative capture. Their changing
//! populations determine opacity, and their energy transfers determine the
//! actual electron temperature used by subsequent rate evaluations.
//!
//! Symmetric composition alternates the allocation-free local reaction kernel
//! with exact monochromatic whole-column subflows. The spectrum itself uses
//! symmetric band composition. Full/two-half comparison controls the combined
//! second-order local error in populations, temperature, reaction integrals
//! and every band photon account. Cross sections prepare once; trial steps do
//! not clone thermodynamic states, allocate species vectors or invoke nested
//! adaptive solves. No LTE projection or artificial electron seed is used.
//!
//! Incident bands have constant number flux over the exposure. Propagation is
//! quasi-static; geometry and volume remain fixed. Optional original collision
//! sources exchange heat between electron/heavy pools, and explicit held powers
//! can drive either pool through the same local reaction kernel. Without these
//! inputs the heavy temperature is fixed. Each
//! radiative capture emits one photon whose chemical plus captured-electron
//! energy LEAVES the modeled system (unit escape probability). This explicit
//! loss closure does not transport diffuse recombination radiation or represent
//! radiation trapping. Scattering, stimulated emission, excited levels,
//! motion and light-travel delay are absent.
//! Source support is required at numerical split stages; unsupported trial
//! temperatures trigger step reduction. Local temporal estimates do not certify
//! spatial discretization, spectral quadrature or empirical model error.

use super::collisional_radiative::{
    CollisionalRadiativeBalance, CollisionalRadiativeError, CollisionalRadiativeSolve,
    CollisionalRadiativeThermal,
};
use super::photoionization::PhotoionizationError;
use fs_qty::{Amount, Area, Time};
use fs_thermochem::equilibrium::composition::FrozenIonizedGasStateV1;
use fs_thermochem::kinetics::atomic_ionization::photoionization::{
    MAX_PHOTON_GROUPS_V1, PhotonNumberFluxV1,
};
use fs_thermochem::kinetics::atomic_ionization::radiative_recombination::AtomicRadiativeRecombinationV1;

mod solve;

/// Immutable gas cell and its original, reciprocal bound-free atomic source.
#[derive(Clone, Copy, Debug)]
pub struct RadiationTransportCell<'a> {
    /// Original gas state, including charge-balancing spectator species.
    pub initial: &'a FrozenIonizedGasStateV1,
    /// One shared atomic channel for collisions, photoionization and capture.
    pub source: &'a AtomicRadiativeRecombinationV1,
    /// Optional original electron/heavy collision sources and held pool powers.
    pub thermal: Option<&'a CollisionalRadiativeThermal>,
}
impl<'a> RadiationTransportCell<'a> {
    /// Fixed-heavy cell with no external pool power or elastic heat exchange.
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
    /// Enable the existing coupled thermal owner without changing photon transport.
    pub const fn with_thermal(mut self, thermal: &'a CollisionalRadiativeThermal) -> Self {
        self.thermal = Some(thermal);
        self
    }
}

/// Explicit local source support, combined integration and spatial work budgets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RadiationTransportSolve {
    /// Applied to the combined evolution. The joule allowance applies to every
    /// cell and once to the whole column; it is not multiplied by cell count.
    pub chemistry: CollisionalRadiativeSolve,
    /// Numerical per-band photon closure relative to that band's incoming count.
    pub relative_photon_tolerance: f64,
    /// Maximum cells admitted before allocation or source evaluation.
    pub maximum_cells: usize,
    /// Maximum discrete bands, also bounded by the L1 source policy.
    pub maximum_groups: usize,
    /// Maximum cell/group product; bounds all band work and account buffers.
    pub maximum_cell_groups: usize,
}
impl RadiationTransportSolve {
    /// Keep the caller's source, quadrature and integration controls explicitly.
    pub fn new(chemistry: CollisionalRadiativeSolve) -> Self {
        Self {
            chemistry,
            relative_photon_tolerance: 1e-12,
            maximum_cells: 4096,
            maximum_groups: MAX_PHOTON_GROUPS_V1,
            maximum_cell_groups: 1_048_576,
        }
    }
}

/// Finite incoming/absorbed/escaped photons for one incident spectral band.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct RadiationTransportBandBalance {
    /// Photons entering the cell or the complete column during the exposure.
    pub incident_photons: f64,
    /// Photons consumed by photoionization, including repeated turnover.
    pub absorbed_photons: f64,
    /// Photons leaving toward the next cell or the external boundary.
    pub transmitted_photons: f64,
    /// Actual absorbed photons times this original band's energy [J].
    pub absorbed_energy_j: f64,
    /// Compensated absorbed plus transmitted minus incoming count.
    pub photon_residual: f64,
}

/// Source-retained cell endpoint and its distinct photon and reaction accounts.
#[derive(Clone, Debug, PartialEq)]
pub struct RadiationTransportCellResult {
    initial: FrozenIonizedGasStateV1,
    state: FrozenIonizedGasStateV1,
    source: AtomicRadiativeRecombinationV1,
    thermal: Option<CollisionalRadiativeThermal>,
    extent: Amount,
    group_balances: Vec<RadiationTransportBandBalance>,
    balance: CollisionalRadiativeBalance,
}
impl RadiationTransportCellResult {
    /// Actual endpoint for material, electrical, collision and optical owners.
    pub const fn state(&self) -> &FrozenIonizedGasStateV1 {
        &self.state
    }
    /// Original gas state and thermodynamic bank.
    pub const fn initial_state(&self) -> &FrozenIonizedGasStateV1 {
        &self.initial
    }
    /// Original reciprocal atomic source, preserved without rebuilding a fit.
    pub const fn source(&self) -> &AtomicRadiativeRecombinationV1 {
        &self.source
    }
    /// Original optional collision/held-power inputs for this cell.
    pub const fn thermal_inputs(&self) -> Option<&CollisionalRadiativeThermal> {
        self.thermal.as_ref()
    }
    /// Signed net ion/electron-pair creation [mol]; gross turnover is separate.
    pub const fn reaction_extent(&self) -> Amount {
        self.extent
    }
    /// Band accounts in boundary-group order. Downstream flux varies with time;
    /// these are integrated photon counts, not a fixed local illumination.
    pub fn group_balances(&self) -> &[RadiationTransportBandBalance] {
        &self.group_balances
    }
    /// Actual gas energy, gross radiation events and signed collisional progress.
    pub const fn balance(&self) -> CollisionalRadiativeBalance {
        self.balance
    }
}

/// Complete spatial solution; no partial endpoints are published on refusal.
#[derive(Clone, Debug, PartialEq)]
pub struct RadiationTransportResult {
    cells: Vec<RadiationTransportCellResult>,
    groups: Vec<PhotonNumberFluxV1>,
    group_balances: Vec<RadiationTransportBandBalance>,
    balance: CollisionalRadiativeBalance,
    area: Area,
    elapsed: Time,
    controls: RadiationTransportSolve,
}
impl RadiationTransportResult {
    /// Ordered endpoints from the illuminated boundary toward the exit.
    pub fn cells(&self) -> &[RadiationTransportCellResult] {
        &self.cells
    }
    /// Original constant boundary photon-energy/number-flux groups.
    pub fn groups(&self) -> &[PhotonNumberFluxV1] {
        &self.groups
    }
    /// Whole-column input, summed absorption and exit counts for each band.
    pub fn group_balances(&self) -> &[RadiationTransportBandBalance] {
        &self.group_balances
    }
    /// Compensated sums of physical accounts and maxima of numerical estimates.
    /// Step counts describe the common column integration, not a sum over cells.
    pub const fn balance(&self) -> CollisionalRadiativeBalance {
        self.balance
    }
    /// Fixed common area normal to every incident band [m²].
    pub const fn illuminated_area(&self) -> Area {
        self.area
    }
    /// Completed physical exposure [s].
    pub const fn elapsed_time(&self) -> Time {
        self.elapsed
    }
    /// Applied source and numerical work controls.
    pub const fn controls(&self) -> RadiationTransportSolve {
        self.controls
    }
}

/// Input, source, budget or conservation refusal, leaving caller cells intact.
#[derive(Clone, Debug, PartialEq)]
pub enum RadiationTransportError {
    /// Invalid common geometry or numerical allowance.
    Invalid(&'static str),
    /// Original local source or numerical integration controls refused.
    Chemistry(CollisionalRadiativeError),
    /// Spatial work budget exceeded before allocation.
    CellLimit { requested: usize, maximum: usize },
    /// Discrete spectrum budget exceeded before allocation.
    GroupLimit { requested: usize, maximum: usize },
    /// Cell/group product exceeds its bound or cannot be represented.
    CellGroupLimit {
        cells: usize,
        groups: usize,
        maximum: usize,
    },
    /// Indexed local source, temperature or state refusal.
    Cell {
        index: usize,
        error: CollisionalRadiativeError,
    },
    /// Indexed photon source/exposure refusal; None identifies boundary input.
    Band {
        cell: Option<usize>,
        group: usize,
        error: PhotoionizationError,
    },
    /// Step count, minimum step or representable time increment exhausted.
    StepBudget,
    /// A physical account, subflow or bounded allocation is unrepresentable.
    Arithmetic(&'static str),
    /// A band fails photon closure at a cell (Some) or the whole column (None).
    PhotonBalance {
        cell: Option<usize>,
        group: usize,
        residual_photons: f64,
    },
    /// Actual gas plus escaped radiation minus absorbed beam and external work [J].
    EnergyBalance { residual_j: f64 },
    /// Cancellation before the complete result is published.
    Cancelled,
}
impl core::fmt::Display for RadiationTransportError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "reacting radiation transport: {self:?}")
    }
}
impl std::error::Error for RadiationTransportError {}

/// Transport finite photon bands while collisions and radiative capture compete.
///
/// Incoming flux is specified at the first cell only. Recombination can restore
/// neutral opacity during illumination; repeated photoionization consumes new
/// incident photons each time. Electron sensible energy and the original
/// chemical gap close against absorbed and escaped radiation. Empty columns
/// transmit all photons; darkness still evolves nonzero charge populations.
/// Source and caloric admission remain mandatory even for dark downstream cells.
pub fn advance_radiation_transport(
    cells: &[RadiationTransportCell<'_>],
    groups: &[PhotonNumberFluxV1],
    area: Area,
    duration: Time,
    controls: RadiationTransportSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<RadiationTransportResult, RadiationTransportError> {
    solve::advance(cells, groups, area, duration, controls, cancelled)
}
