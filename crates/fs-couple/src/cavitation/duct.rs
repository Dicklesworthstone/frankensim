//! Event-resolved liquid-column separation on an original material DG duct.
//!
//! This is a discrete, isothermal, pure-vapor endpoint cavity, not a spherical
//! bubble or a two-phase bulk equation of state. The retained liquid uses the
//! native small-strain acoustic operator, including any original wall memory.
//! The caller supplies the absolute background pressure and a small allowed
//! cavity volume. Trapped gas, nucleation delay, interface capillarity, moving
//! meshes and latent-heat feedback are outside this equilibrium cavity model.
//!
//! With native energy coordinates `x`, generator `L` and inflow column `b`,
//! `p_abs = p_ref + b^T x`. A filled endpoint has `Q=0`. While the cavity is
//! present, `Q = -b^T L x / (b^T b)` keeps the INSTANTANEOUS pressure at vapor
//! pressure, and `V' = Q`: positive flow into the liquid grows the cavity.
//! Rupture and collapse are localized on the native RK45 dense trajectory.
//! Neither the acoustic state nor the integrated volume is reset or clipped.
//! Tiny signed root residuals remain visible under explicit numerical bounds.
//!
//! The original duct adopts the actual integrated `p Q` work, pressure impulse
//! and original quadratic dissipation. Together with cavity potential
//! `(p_ref-p_sat) V`, its energy closes without inventing heat from void volume.

use fs_duct::DuctError;
use fs_duct::dg::{DgBoundary, DgDuct, DgFluidMaterial, DgState, DgStepOutput, PortId};
use fs_material::MaterialError;
use fs_material::cavitation::{CavitationProperties, CavitationState};
use fs_material::fluid_acoustics::FluidAcoustics;
use fs_time::adaptive::events::EventError;
use fs_time::adaptive::events::multiple::EventSetError;

/// Simultaneous discrete vapor pockets at original endpoint/interior ports.
pub mod pockets;

/// Complete original acoustic/wall state retained by this bounded dense lane.
pub const MAX_CAVITATING_DUCT_DOFS: usize = 128;

/// Physical admission of the frozen liquid-column approximation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DuctCavityLimits {
    /// Absolute reference for every native DG pressure perturbation [Pa].
    pub reference_pressure_pa: f64,
    /// Total cavity-volume ceiling [m³]. Each station additionally enforces a
    /// small volume relative to its neighboring liquid-column spans.
    pub maximum_cavity_volume_m3: f64,
    /// Bound on |nodal pressure perturbation|/(rho c²), in (0, 0.1].
    pub maximum_liquid_strain: f64,
    /// Bound on each actual liquid-interface velocity / c, in (0, 0.1].
    /// Internal pockets check both native through-flow traces separately.
    pub maximum_interface_mach: f64,
    /// Known saturated rho_v/rho_l must not exceed this bound, in (0, 0.1].
    /// Missing optional vapor density is an explicit dilute-vapor assumption.
    pub maximum_vapor_density_ratio: f64,
}

/// Tolerances and finite work for one complete physical interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DuctCavityNumerics {
    pub relative_tolerance: f64,
    pub absolute_energy_tolerance_j: f64,
    pub absolute_pressure_tolerance_pa: f64,
    pub absolute_volume_tolerance_m3: f64,
    /// Bound the accepted step and resolve the shortest wave/contact time.
    pub maximum_step_s: f64,
    pub event_time_tolerance_s: f64,
    /// Ordered samples per accepted dense step, in 1..=4096. Finite scans do
    /// not certify tangencies or multiple crossings between adjacent samples.
    pub event_scan_substeps: usize,
    pub maximum_attempts: usize,
    pub maximum_events: usize,
}

/// Which complementarity branch currently owns the endpoint boundary.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DuctCavityMode {
    LiquidFilled,
    VaporCavity,
}

/// Accepted native state, cavity coordinate and one physical clock.
#[derive(Clone, Debug, PartialEq)]
pub struct DuctCavityState {
    duct: DgState,
    time_s: f64,
    volume_m3: f64,
    mode: DuctCavityMode,
    pockets: Vec<pockets::DuctPocketState>,
}
impl DuctCavityState {
    pub const fn duct(&self) -> &DgState {
        &self.duct
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    /// Integrated cavity volume, WITHOUT a positivity clip. A localized
    /// collapse can leave a signed residual bounded by the declared absolute
    /// volume tolerance; a filled branch retains that residual verbatim.
    pub const fn volume_m3(&self) -> f64 {
        self.volume_m3
    }
    pub const fn mode(&self) -> DuctCavityMode {
        self.mode
    }
    /// Every pocket in the original native port order. The legacy volume/mode
    /// getters above observe the first pocket; this slice is the multi-pocket view.
    pub fn pockets(&self) -> &[pockets::DuctPocketState] {
        &self.pockets
    }
}

/// A source/domain/numerical refusal never replaces the accepted history.
#[derive(Clone, Debug, PartialEq)]
pub enum DuctCavityError {
    Material(MaterialError),
    Duct(DuctError),
    Integration(EventError),
    EventSetIntegration(EventSetError),
    ActiveConstraints(&'static str),
    InvalidInput(&'static str),
    OutsideDomain(&'static str),
    RootResolution(&'static str),
    AttemptLimit,
    EventLimit,
    Cancelled,
    EnergyBalance,
}
impl core::fmt::Display for DuctCavityError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "cavitating liquid duct: {self:?}")
    }
}
impl std::error::Error for DuctCavityError {}
impl From<MaterialError> for DuctCavityError {
    fn from(value: MaterialError) -> Self {
        Self::Material(value)
    }
}
impl From<DuctError> for DuctCavityError {
    fn from(value: DuctError) -> Self {
        Self::Duct(value)
    }
}
impl From<EventError> for DuctCavityError {
    fn from(value: EventError) -> Self {
        Self::Integration(value)
    }
}

impl From<EventSetError> for DuctCavityError {
    fn from(value: EventSetError) -> Self {
        Self::EventSetIntegration(value)
    }
}

/// Native port and complete coupled-storage accounting for an accepted step.
#[derive(Clone, Debug)]
pub struct DuctCavityStep {
    pub duct: DgStepOutput,
    pub mode: DuctCavityMode,
    pub end_absolute_pressure_pa: f64,
    pub end_inflow_m3_s: f64,
    pub cavity_volume_m3: f64,
    /// Integral Q dt = V_end-V_start, positive when liquid retreats.
    pub displaced_liquid_volume_m3: f64,
    pub cavity_potential_j: f64,
    pub total_stored_energy_j: f64,
    /// Change in duct+cavity storage plus independently integrated native loss.
    pub energy_balance_residual_j: f64,
    pub rupture_events: usize,
    pub collapse_events: usize,
    pub last_rupture_time_s: Option<f64>,
    pub last_collapse_time_s: Option<f64>,
    pub accepted_substeps: usize,
    pub rejected_substeps: usize,
}

/// One original liquid duct with one cavitating ENDPOINT and no extra ports.
/// The other endpoint may retain its native closed, matched or pressure-release
/// condition. The source, every DG coordinate and wall memory remain original.
pub struct CavitatingDuct<F: CavitationProperties + FluidAcoustics + Clone> {
    inner: pockets::CavitatingDuctPockets<F>,
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn fraction(value: f64) -> bool {
    positive(value) && value <= 0.1
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn multiply(matrix: &[f64], x: &[f64], out: &mut [f64]) {
    for (row, value) in matrix.chunks_exact(x.len()).zip(out) {
        *value = dot(row, x);
    }
}
fn quadratic_loss(generator: &[f64], x: &[f64]) -> f64 {
    let n = x.len();
    let mut loss = 0.0;
    for i in 0..n {
        loss -= generator[i * n + i] * x[i] * x[i];
        for j in 0..i {
            loss -= (generator[i * n + j] + generator[j * n + i]) * x[i] * x[j];
        }
    }
    loss
}

impl DuctCavityNumerics {
    fn admit(&self) -> Result<(), DuctCavityError> {
        if !positive(self.relative_tolerance)
            || self.relative_tolerance >= 0.1
            || ![
                self.absolute_energy_tolerance_j,
                self.absolute_pressure_tolerance_pa,
                self.absolute_volume_tolerance_m3,
                self.maximum_step_s,
                self.event_time_tolerance_s,
            ]
            .into_iter()
            .all(positive)
            || self.event_time_tolerance_s > self.maximum_step_s
            || !(1..=4096).contains(&self.event_scan_substeps)
            || self.maximum_attempts == 0
            || self.maximum_events == 0
        {
            return Err(DuctCavityError::InvalidInput(
                "positive finite tolerances and bounded nonzero work required",
            ));
        }
        Ok(())
    }
}

impl<F: CavitationProperties + FluidAcoustics + Clone> CavitatingDuct<F> {
    /// Bind one original endpoint cavity. Single and multiple pockets use the
    /// SAME coupled constraints, dense events, trajectory and native energy adoption.
    pub fn new(
        duct: DgDuct<DgFluidMaterial<F>>,
        initial: DgState,
        time_s: f64,
        initial_cavity_volume_m3: f64,
        limits: DuctCavityLimits,
    ) -> Result<Self, DuctCavityError> {
        let view = duct.continuous_operator();
        if view.ports().len() != 1 {
            return Err(DuctCavityError::InvalidInput(
                "one original endpoint port required by the single-cavity facade",
            ));
        }
        let port = view.ports()[0].0;
        if duct.inlet() != DgBoundary::Port(port) && duct.outlet() != DgBoundary::Port(port) {
            return Err(DuctCavityError::InvalidInput(
                "use CavitatingDuctPockets for an interior vapor pocket",
            ));
        }
        let stations = [pockets::DuctPocket {
            port,
            initial_volume_m3: initial_cavity_volume_m3,
            maximum_volume_m3: limits.maximum_cavity_volume_m3,
        }];
        Ok(Self {
            inner: pockets::CavitatingDuctPockets::new(duct, initial, time_s, &stations, limits)?,
        })
    }

    pub const fn duct(&self) -> &DgDuct<DgFluidMaterial<F>> {
        self.inner.duct()
    }
    pub fn source(&self) -> &F {
        self.inner.source()
    }
    pub const fn material(&self) -> CavitationState {
        self.inner.material()
    }
    pub const fn state(&self) -> &DuctCavityState {
        self.inner.state()
    }
    pub const fn limits(&self) -> DuctCavityLimits {
        self.inner.limits()
    }

    /// Atomic advancement through the shared multi-pocket trajectory owner.
    pub fn advance<Cancel: FnMut() -> bool>(
        &mut self,
        end_time_s: f64,
        numerics: &DuctCavityNumerics,
        cancelled: &mut Cancel,
    ) -> Result<DuctCavityStep, DuctCavityError> {
        let report = self.inner.advance(end_time_s, numerics, cancelled)?;
        let pocket = report.pockets[0];
        Ok(DuctCavityStep {
            duct: report.duct,
            mode: pocket.mode,
            end_absolute_pressure_pa: pocket.end_absolute_pressure_pa,
            end_inflow_m3_s: pocket.end_inflow_m3_s,
            cavity_volume_m3: pocket.volume_m3,
            displaced_liquid_volume_m3: pocket.displaced_liquid_volume_m3,
            cavity_potential_j: report.cavity_potential_j,
            total_stored_energy_j: report.total_stored_energy_j,
            energy_balance_residual_j: report.energy_balance_residual_j,
            rupture_events: pocket.rupture_events,
            collapse_events: pocket.collapse_events,
            last_rupture_time_s: pocket.last_rupture_time_s,
            last_collapse_time_s: pocket.last_collapse_time_s,
            accepted_substeps: report.accepted_substeps,
            rejected_substeps: report.rejected_substeps,
        })
    }
}

#[cfg(test)]
mod tests;
