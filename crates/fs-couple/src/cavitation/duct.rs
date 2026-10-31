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
use fs_duct::dg::continuous::state::DgTrajectoryPort;
use fs_duct::dg::{DgBoundary, DgDuct, DgFluidMaterial, DgState, DgStepOutput, PortId};
use fs_material::MaterialError;
use fs_material::cavitation::{CavitationProperties, CavitationState};
use fs_material::fluid_acoustics::FluidAcoustics;
use fs_time::adaptive::events::{
    EventAdvance, EventDirection, EventError, EventOptions, InitialEvent, rk45_until_event,
};
use fs_time::adaptive::{AdaptiveState, AdaptiveStatus, PiController};
use std::cell::RefCell;

/// Complete original acoustic/wall state retained by this bounded dense lane.
pub const MAX_CAVITATING_DUCT_DOFS: usize = 128;

/// Physical admission of the frozen liquid-column approximation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DuctCavityLimits {
    /// Absolute reference for every native DG pressure perturbation [Pa].
    pub reference_pressure_pa: f64,
    /// Positive cavity-volume ceiling [m³], at most 10% of endpoint area ×
    /// column length. A smaller ceiling is needed when geometry varies rapidly.
    pub maximum_cavity_volume_m3: f64,
    /// Bound on |nodal pressure perturbation|/(rho c²), in (0, 0.1].
    pub maximum_liquid_strain: f64,
    /// Bound on |cavity inflow|/(endpoint area × liquid c), in (0, 0.1].
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
}

/// A source/domain/numerical refusal never replaces the accepted history.
#[derive(Clone, Debug, PartialEq)]
pub enum DuctCavityError {
    Material(MaterialError),
    Duct(DuctError),
    Integration(EventError),
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
    duct: DgDuct<DgFluidMaterial<F>>,
    material: CavitationState,
    limits: DuctCavityLimits,
    port: PortId,
    port_area_m2: f64,
    cavity_node: usize,
    state: DuctCavityState,
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
    /// Bind an existing native material operator and its OWN initial state.
    /// A nonempty initial cavity requires p_abs=p_sat to floating-point roundoff;
    /// a filled endpoint requires p_abs>=p_sat. No pressure is projected.
    pub fn new(
        duct: DgDuct<DgFluidMaterial<F>>,
        initial: DgState,
        time_s: f64,
        initial_cavity_volume_m3: f64,
        limits: DuctCavityLimits,
    ) -> Result<Self, DuctCavityError> {
        let material = duct.material().cavitation_state()?;
        if material.acoustic_state()? != duct.acoustic_material() {
            return Err(DuctCavityError::InvalidInput(
                "cavitation and native duct acoustic coefficients must share one material state",
            ));
        }
        let view = duct.continuous_operator();
        let coordinates = view.energy_coordinates(&initial)?;
        if view.dimension() > MAX_CAVITATING_DUCT_DOFS || view.ports().len() != 1 {
            return Err(DuctCavityError::InvalidInput(
                "one original endpoint port and at most 128 complete native coordinates required",
            ));
        }
        let (port, column) = &view.ports()[0];
        let area = if duct.inlet() == DgBoundary::Port(*port) {
            duct.inlet_area_m2()
        } else if duct.outlet() == DgBoundary::Port(*port) {
            duct.outlet_area_m2()
        } else {
            return Err(DuctCavityError::InvalidInput(
                "a side-port is not an endpoint liquid-column cavity",
            ));
        };
        if !time_s.is_finite()
            || time_s < 0.0
            || !initial_cavity_volume_m3.is_finite()
            || initial_cavity_volume_m3 < 0.0
            || !limits.reference_pressure_pa.is_finite()
            || limits.reference_pressure_pa < material.saturation_pressure_pa()
            || !positive(limits.maximum_cavity_volume_m3)
            || limits.maximum_cavity_volume_m3 > 0.1 * area * duct.total_length()
            || initial_cavity_volume_m3 > limits.maximum_cavity_volume_m3
            || !fraction(limits.maximum_liquid_strain)
            || !fraction(limits.maximum_interface_mach)
            || !fraction(limits.maximum_vapor_density_ratio)
            || coordinates.iter().any(|x| !x.is_finite())
            || !positive(dot(column, column))
        {
            return Err(DuctCavityError::InvalidInput(
                "invalid liquid cavity clock, reference, volume or physical limits",
            ));
        }
        if material.saturated_vapor().is_some_and(|v| {
            v.density_kg_m3() / material.liquid_density_kg_m3() > limits.maximum_vapor_density_ratio
        }) {
            return Err(DuctCavityError::OutsideDomain(
                "vapor mass is not negligible relative to liquid",
            ));
        }
        let pressure = limits.reference_pressure_pa + dot(column, coordinates);
        let roundoff = 128.0 * f64::EPSILON * limits.reference_pressure_pa.max(1.0);
        let gap = pressure - material.saturation_pressure_pa();
        if gap < -roundoff || (initial_cavity_volume_m3 > 0.0 && gap.abs() > roundoff) {
            return Err(DuctCavityError::InvalidInput(
                "initial pressure and cavity volume violate complementarity",
            ));
        }
        let mut derivative = vec![0.0; coordinates.len()];
        multiply(view.generator_per_s(), coordinates, &mut derivative);
        let mode = if initial_cavity_volume_m3 > 0.0
            || (gap.abs() <= roundoff && dot(column, &derivative) < 0.0)
        {
            DuctCavityMode::VaporCavity
        } else {
            DuctCavityMode::LiquidFilled
        };
        let port = *port;
        let cavity_node = if duct.inlet() == DgBoundary::Port(port) {
            0
        } else {
            duct.field_dofs() - 1
        };
        if mode == DuctCavityMode::VaporCavity {
            let inflow = -dot(column, &derivative) / dot(column, column);
            if !inflow.is_finite()
                || (inflow / (area * material.liquid_sound_speed_m_s())).abs()
                    > limits.maximum_interface_mach
            {
                return Err(DuctCavityError::OutsideDomain(
                    "initial interface speed exceeds declared liquid Mach limit",
                ));
            }
        }
        let model = Self {
            duct,
            material,
            limits,
            port,
            port_area_m2: area,
            cavity_node,
            state: DuctCavityState {
                duct: initial,
                time_s,
                volume_m3: initial_cavity_volume_m3,
                mode,
            },
        };
        model.admit_field(
            model
                .duct
                .continuous_operator()
                .energy_coordinates(&model.state.duct)?,
            1.0,
            roundoff,
            false,
        )?;
        Ok(model)
    }

    pub const fn duct(&self) -> &DgDuct<DgFluidMaterial<F>> {
        &self.duct
    }
    pub fn source(&self) -> &F {
        self.duct.material()
    }
    pub const fn material(&self) -> CavitationState {
        self.material
    }
    pub const fn state(&self) -> &DuctCavityState {
        &self.state
    }
    pub const fn limits(&self) -> DuctCavityLimits {
        self.limits
    }

    fn admit_field(
        &self,
        normalized: &[f64],
        scale: f64,
        pressure_tolerance: f64,
        allow_cavity_overshoot: bool,
    ) -> Result<(), DuctCavityError> {
        // Original physical mass maps, including wall-memory ordering. This is
        // not a reconstructed acoustic state or a second liquid mass operator.
        let (pressures, _) = self.duct.from_energy_coords(normalized);
        let bulk = self.duct.acoustic_material().adiabatic_bulk_modulus_pa();
        if pressures
            .iter()
            .any(|p| !p.is_finite() || (p * scale / bulk).abs() > self.limits.maximum_liquid_strain)
        {
            return Err(DuctCavityError::OutsideDomain(
                "native liquid pressure exceeds the declared small-strain limit",
            ));
        }
        if pressures.iter().enumerate().any(|(i, p)| {
            !(allow_cavity_overshoot && i == self.cavity_node)
                && self.limits.reference_pressure_pa + p * scale
                    < self.material.saturation_pressure_pa() - pressure_tolerance
        }) {
            return Err(DuctCavityError::OutsideDomain(
                "another native liquid node reaches vapor pressure; distributed cavitation is outside the endpoint model",
            ));
        }
        Ok(())
    }

    /// Prepare and atomically adopt one complete coupled interval. Cancellation,
    /// unresolved roots, work limits, domain failure or energy failure preserve
    /// the accepted cavity AND the complete native duct history.
    pub fn advance<Cancel: FnMut() -> bool>(
        &mut self,
        end_time_s: f64,
        numerics: &DuctCavityNumerics,
        cancelled: &mut Cancel,
    ) -> Result<DuctCavityStep, DuctCavityError> {
        self.advance_kernel(end_time_s, numerics, cancelled)
    }

    fn advance_kernel(
        &mut self,
        end_time_s: f64,
        numerics: &DuctCavityNumerics,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<DuctCavityStep, DuctCavityError> {
        numerics.admit()?;
        if cancelled() {
            return Err(DuctCavityError::Cancelled);
        }
        let dt = end_time_s - self.state.time_s;
        if !end_time_s.is_finite() || !positive(dt) {
            return Err(DuctCavityError::InvalidInput(
                "the physical endpoint must be finite and strictly later",
            ));
        }
        let view = self.duct.continuous_operator();
        let original = view.energy_coordinates(&self.state.duct)?;
        let n = original.len();
        let generator = view.generator_per_s();
        let column = &view.ports()[0].1;
        let column_norm = dot(column, column);
        let vapor = self.material.saturation_pressure_pa();
        let pressure_scale = self.limits.reference_pressure_pa.max(1.0);
        let volume_scale = self.limits.maximum_cavity_volume_m3;
        let energy_scale = (2.0 * self.state.duct.stored_energy_j())
            .max(pressure_scale * volume_scale)
            .max(numerics.absolute_energy_tolerance_j);
        let coordinate_scale = energy_scale.sqrt();
        let impulse_scale = pressure_scale * dt;
        let atol = (numerics.absolute_energy_tolerance_j / energy_scale)
            .min(numerics.absolute_volume_tolerance_m3 / volume_scale)
            .min(numerics.absolute_pressure_tolerance_pa / pressure_scale);
        if !positive(energy_scale)
            || !positive(impulse_scale)
            || !positive(atol)
            || numerics.absolute_volume_tolerance_m3 >= volume_scale
        {
            return Err(DuctCavityError::InvalidInput(
                "numerical scales overflow or cannot resolve the admitted cavity",
            ));
        }
        let mut y: Vec<_> = original.iter().map(|x| x / coordinate_scale).collect();
        y.extend([self.state.volume_m3 / volume_scale, 0.0, 0.0, 0.0]);
        let mut integration =
            AdaptiveState::new(self.state.time_s, &y, numerics.maximum_step_s.min(dt));
        let mut mode = self.state.mode;
        let mut rupture_events = 0;
        let mut collapse_events = 0;
        let mut last_rupture_time_s = None;
        let mut last_collapse_time_s = None;
        let options = EventOptions {
            direction: EventDirection::Falling,
            initial: InitialEvent::Ignore,
            max_step: numerics.maximum_step_s,
            scan_substeps: numerics.event_scan_substeps,
            time_tolerance: numerics.event_time_tolerance_s,
            max_iterations: 128,
        };
        let pi = PiController::default();
        while integration.t < end_time_s {
            if cancelled() {
                return Err(DuctCavityError::Cancelled);
            }
            let attempts = integration
                .accepted
                .checked_add(integration.rejected)
                .ok_or(DuctCavityError::AttemptLimit)?;
            if attempts >= numerics.maximum_attempts {
                return Err(DuctCavityError::AttemptLimit);
            }
            let failure = RefCell::new(None);
            let rhs = |_: f64, u: &[f64], out: &mut [f64]| {
                multiply(generator, &u[..n], &mut out[..n]);
                let inflow = if mode == DuctCavityMode::VaporCavity {
                    -coordinate_scale * dot(column, &out[..n]) / column_norm
                } else {
                    0.0
                };
                let pressure = coordinate_scale * dot(column, &u[..n]);
                let domain = self
                    .admit_field(
                        &u[..n],
                        coordinate_scale,
                        numerics.absolute_pressure_tolerance_pa,
                        mode == DuctCavityMode::LiquidFilled,
                    )
                    .and_then(|()| {
                        if (inflow / (self.port_area_m2 * self.material.liquid_sound_speed_m_s()))
                            .abs()
                            > self.limits.maximum_interface_mach
                        {
                            return Err(DuctCavityError::OutsideDomain(
                                "interface speed exceeds declared liquid Mach limit",
                            ));
                        }
                        if u[n] * volume_scale > self.limits.maximum_cavity_volume_m3 {
                            return Err(DuctCavityError::OutsideDomain(
                                "cavity exceeds the fixed-column volume limit",
                            ));
                        }
                        if mode == DuctCavityMode::VaporCavity
                            && (self.limits.reference_pressure_pa + pressure - vapor).abs()
                                > numerics.absolute_pressure_tolerance_pa
                        {
                            return Err(DuctCavityError::RootResolution(
                                "active pressure is outside the vapor-pressure root tolerance",
                            ));
                        }
                        Ok(())
                    });
                if let Err(error) = domain {
                    if failure.borrow().is_none() {
                        *failure.borrow_mut() = Some(error);
                    }
                    out.fill(f64::NAN);
                    return;
                }
                for (value, b) in out[..n].iter_mut().zip(column) {
                    *value += b * (inflow / coordinate_scale);
                }
                out[n] = inflow / volume_scale;
                out[n + 1] = pressure / impulse_scale;
                out[n + 2] = pressure * inflow / energy_scale;
                out[n + 3] = quadratic_loss(generator, &u[..n]);
            };
            let guard = |_: f64, u: &[f64]| {
                if mode == DuctCavityMode::VaporCavity {
                    u[n]
                } else {
                    (self.limits.reference_pressure_pa + coordinate_scale * dot(column, &u[..n])
                        - vapor)
                        / pressure_scale
                }
            };
            let start = integration.t;
            let result = rk45_until_event(
                &mut integration,
                &rhs,
                &guard,
                end_time_s,
                numerics.relative_tolerance,
                atol,
                &pi,
                numerics.maximum_attempts - attempts,
                &options,
                &mut || cancelled(),
            );
            if let Some(error) = failure.into_inner() {
                return Err(error);
            }
            match result? {
                EventAdvance::Stopped(report) => match report.status {
                    AdaptiveStatus::ReachedEnd => break,
                    AdaptiveStatus::StepLimit => return Err(DuctCavityError::AttemptLimit),
                    AdaptiveStatus::Cancelled => return Err(DuctCavityError::Cancelled),
                },
                EventAdvance::Event(_) => {
                    if integration.t <= start {
                        return Err(DuctCavityError::RootResolution(
                            "a distinct contact time could not be represented",
                        ));
                    }
                    if rupture_events + collapse_events >= numerics.maximum_events {
                        return Err(DuctCavityError::EventLimit);
                    }
                    let mut derivative = vec![0.0; n];
                    multiply(generator, &integration.u[..n], &mut derivative);
                    let constrained_flow =
                        -coordinate_scale * dot(column, &derivative) / column_norm;
                    if mode == DuctCavityMode::VaporCavity {
                        if (integration.u[n] * volume_scale).abs()
                            > numerics.absolute_volume_tolerance_m3
                            || constrained_flow >= 0.0
                        {
                            return Err(DuctCavityError::RootResolution(
                                "collapse root or direction is unresolved",
                            ));
                        }
                        collapse_events += 1;
                        last_collapse_time_s = Some(integration.t);
                        mode = DuctCavityMode::LiquidFilled;
                    } else {
                        let gap = self.limits.reference_pressure_pa
                            + coordinate_scale * dot(column, &integration.u[..n])
                            - vapor;
                        if gap.abs() > numerics.absolute_pressure_tolerance_pa
                            || constrained_flow <= 0.0
                        {
                            return Err(DuctCavityError::RootResolution(
                                "rupture root or direction is unresolved",
                            ));
                        }
                        rupture_events += 1;
                        last_rupture_time_s = Some(integration.t);
                        mode = DuctCavityMode::VaporCavity;
                    }
                }
            }
        }
        let volume = integration.u[n] * volume_scale;
        let endpoint: Vec<_> = integration.u[..n]
            .iter()
            .map(|x| x * coordinate_scale)
            .collect();
        let pressure = self.limits.reference_pressure_pa + dot(column, &endpoint);
        if volume < -numerics.absolute_volume_tolerance_m3
            || volume > self.limits.maximum_cavity_volume_m3
            || (mode == DuctCavityMode::LiquidFilled
                && (volume.abs() > numerics.absolute_volume_tolerance_m3
                    || pressure < vapor - numerics.absolute_pressure_tolerance_pa))
            || (mode == DuctCavityMode::VaporCavity
                && (pressure - vapor).abs() > numerics.absolute_pressure_tolerance_pa)
        {
            return Err(DuctCavityError::RootResolution(
                "endpoint complementarity is outside numerical root bounds",
            ));
        }
        let loss = integration.u[n + 3] * energy_scale;
        let receipt = DgTrajectoryPort {
            id: self.port,
            pressure_impulse_pa_s: integration.u[n + 1] * impulse_scale,
            work_into_duct_j: integration.u[n + 2] * energy_scale,
        };
        let (candidate, duct_report) = view.candidate_from_trajectory(
            &self.state.duct,
            endpoint.clone(),
            dt,
            &[receipt],
            loss,
            numerics.absolute_energy_tolerance_j,
            numerics.relative_tolerance,
        )?;
        let potential = (self.limits.reference_pressure_pa - vapor) * volume;
        let before = self.state.duct.stored_energy_j()
            + (self.limits.reference_pressure_pa - vapor) * self.state.volume_m3;
        let total = candidate.stored_energy_j() + potential;
        let residual = total - before + loss;
        let energy_tolerance = numerics.absolute_energy_tolerance_j
            + numerics.relative_tolerance * (total.abs() + before.abs() + loss.abs());
        if !residual.is_finite() || residual.abs() > energy_tolerance {
            return Err(DuctCavityError::EnergyBalance);
        }
        let mut derivative = vec![0.0; n];
        multiply(generator, &endpoint, &mut derivative);
        let inflow = if mode == DuctCavityMode::VaporCavity {
            -dot(column, &derivative) / column_norm
        } else {
            0.0
        };
        let report = DuctCavityStep {
            duct: duct_report,
            mode,
            end_absolute_pressure_pa: pressure,
            end_inflow_m3_s: inflow,
            cavity_volume_m3: volume,
            displaced_liquid_volume_m3: volume - self.state.volume_m3,
            cavity_potential_j: potential,
            total_stored_energy_j: total,
            energy_balance_residual_j: residual,
            rupture_events,
            collapse_events,
            last_rupture_time_s,
            last_collapse_time_s,
            accepted_substeps: integration.accepted,
            rejected_substeps: integration.rejected,
        };
        if cancelled() {
            return Err(DuctCavityError::Cancelled);
        }
        self.state = DuctCavityState {
            duct: candidate,
            time_s: end_time_s,
            volume_m3: volume,
            mode,
        };
        Ok(report)
    }
}

#[cfg(test)]
mod tests;
