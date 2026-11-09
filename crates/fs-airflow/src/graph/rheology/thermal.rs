//! Original shear-dependent hydraulic resistance coupled to sensible heat.
//!
//! This is an explicitly uniform, recirculating thermal inventory: the caller's
//! total C(T) includes the pipes AND all returning supply fluid. Every hydraulic
//! pressure/injection port exchanges fluid at that same temperature. Pump work
//! enters through the existing pressure ports; motor losses, external donor
//! enthalpy, mass, density, heat capacity and thermal equilibration are not
//! inferred from a viscosity law. An open circuit needs an enthalpy owner.
//!
//! For held pressure/flow ports, imposed power Q and bath conductance H, the
//! implicit midpoint balance is
//! `E(T1)-E(T0) = dt * [D((T0+T1)/2) + Q + H*(Tb-(T0+T1)/2)]`.
//! Every D is the actual original-law Carreau graph dissipation. Each thermal
//! trial re-prepares the same source laws, geometry and physical rate ceilings;
//! the native nonlinear pressure solve retains its original work bounds.
//! One cloned trial topology is reused. The existing hydraulic solver and its
//! returned profiles still allocate O(nodes+pipes) storage for each trial.
//!
//! Source/caloric refusal is terminal for the proposed interval. Initial,
//! midpoint and final hydraulics are admitted before an immutable next owner
//! is published. A bounded bracket/refinement solve and endpoint power-sign
//! guard reject unresolved intervals and detected equilibrium overshoots.
//! These checks do not isolate all equilibria, certify uniqueness or source
//! support along an arbitrary intermediate trajectory, or replace temporal
//! refinement. Hydraulics are quasisteady: radial inertia, hydraulic storage,
//! spatial temperature fields, open-stream enthalpy, phase change and thermal
//! expansion remain outside this reduction.

/// Resolved radial momentum, hydraulic storage and one shared caloric inventory.
pub mod transient;

use fs_material::thermal::transport::{ThermalTransportError, storage::TabulatedHeatStorage};

use super::{
    CarreauGraphControls, CarreauGraphError, CarreauGraphSolution, CarreauPipeGraph, Cx,
    FixedPressure, FlowInjection, GraphError, LaminarCarreauPipe, Sum, checkpoint, finite,
};

/// Includes every attempted full hydraulic solve, including final admission.
pub const MAX_THERMAL_HYDRAULIC_SOLVES: usize = 256;

/// Independent thermal accuracy/work bounds and the original hydraulic controls.
#[derive(Clone, Copy, Debug)]
pub struct ThermalCarreauControls {
    pub hydraulics: CarreauGraphControls,
    pub maximum_hydraulic_solves: usize,
    pub maximum_temperature_change_k: f64,
    pub temperature_tolerance_k: f64,
    pub absolute_energy_tolerance_j: f64,
    pub relative_energy_tolerance: f64,
}

impl Default for ThermalCarreauControls {
    fn default() -> Self {
        Self {
            hydraulics: CarreauGraphControls::default(),
            maximum_hydraulic_solves: 96,
            maximum_temperature_change_k: 5.0,
            temperature_tolerance_k: 1e-8,
            absolute_energy_tolerance_j: 1e-10,
            relative_energy_tolerance: 1e-9,
        }
    }
}

impl ThermalCarreauControls {
    fn validate(self) -> Result<(), ThermalCarreauError> {
        self.hydraulics.validate()?;
        if !(2..=MAX_THERMAL_HYDRAULIC_SOLVES).contains(&self.maximum_hydraulic_solves)
            || !positive(self.maximum_temperature_change_k)
            || !positive(self.temperature_tolerance_k)
            || !nonnegative(self.absolute_energy_tolerance_j)
            || !nonnegative(self.relative_energy_tolerance)
            || self.relative_energy_tolerance >= 1.0
            || self.absolute_energy_tolerance_j == 0.0 && self.relative_energy_tolerance == 0.0
        {
            return Err(GraphError::InvalidInput("invalid Carreau thermal controls").into());
        }
        Ok(())
    }
}

/// Held heat input and whole-inventory exchange with one prescribed bath.
/// Native hydraulic heat is already deposited; do not include it in Q again.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalCarreauBoundary {
    pub prescribed_power_w: f64,
    pub bath_conductance_w_k: f64,
    pub bath_temperature_k: f64,
}

impl ThermalCarreauBoundary {
    fn validate(self) -> Result<(), ThermalCarreauError> {
        if !self.prescribed_power_w.is_finite()
            || !nonnegative(self.bath_conductance_w_k)
            || !positive(self.bath_temperature_k)
        {
            return Err(GraphError::InvalidInput("invalid Carreau thermal boundary").into());
        }
        Ok(())
    }

    fn bath_power(self, temperature_k: f64) -> Result<f64, GraphError> {
        finite(
            self.bath_conductance_w_k * (self.bath_temperature_k - temperature_k),
            "Carreau bath power",
        )
    }
}

#[derive(Debug)]
pub enum ThermalCarreauError {
    Hydraulic(CarreauGraphError),
    TransientHydraulic(super::transient::TransientCarreauError),
    Caloric(ThermalTransportError),
    TemperatureDomain { temperature_k: f64 },
    HydraulicSolveBudget,
    UnbracketedTemperature,
    TemperatureResolution,
    EquilibriumCrossing,
    EnergyBalance,
}

impl From<GraphError> for ThermalCarreauError {
    fn from(value: GraphError) -> Self {
        Self::Hydraulic(value.into())
    }
}
impl From<CarreauGraphError> for ThermalCarreauError {
    fn from(value: CarreauGraphError) -> Self {
        Self::Hydraulic(value)
    }
}
impl From<super::transient::TransientCarreauError> for ThermalCarreauError {
    fn from(value: super::transient::TransientCarreauError) -> Self {
        Self::TransientHydraulic(value)
    }
}
impl From<ThermalTransportError> for ThermalCarreauError {
    fn from(value: ThermalTransportError) -> Self {
        Self::Caloric(value)
    }
}
impl core::fmt::Display for ThermalCarreauError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "thermal Carreau graph: {self:?}")
    }
}
impl std::error::Error for ThermalCarreauError {}

/// Accepted original pipe graph and the declared recirculating caloric inventory.
/// The authored storage is borrowed unchanged; it supplies no inferred mass or
/// evidence-qualified heat capacity. All original pipes share one temperature.
#[derive(Clone, Debug)]
pub struct ThermalCarreauGraph<'a> {
    graph: CarreauPipeGraph,
    storage: &'a TabulatedHeatStorage,
    temperature_span_k: [f64; 2],
    time_s: f64,
}

/// Physical integrated exchanges and separate native numerical work defects.
#[derive(Clone, Debug)]
pub struct ThermalCarreauReport {
    pub duration_s: f64,
    pub initial_temperature_k: f64,
    pub midpoint_temperature_k: f64,
    pub final_temperature_k: f64,
    pub stored_energy_change_j: f64,
    /// The only hydraulic contribution deposited in caloric storage.
    pub viscous_heat_j: f64,
    pub prescribed_heat_j: f64,
    /// Signed bath heat into the inventory.
    pub bath_heat_j: f64,
    pub reservoir_work_j: f64,
    pub injection_work_j: f64,
    pub constitutive_residual_work_j: f64,
    pub continuity_residual_work_j: f64,
    /// Native reservoir/injection work minus original viscous heat.
    pub hydraulic_work_residual_j: f64,
    /// Storage minus viscous, imposed and bath heat.
    pub thermal_residual_j: f64,
    /// Storage minus reservoir/injection work, imposed and bath heat.
    pub combined_energy_residual_j: f64,
    pub thermal_residual_allowance_j: f64,
    pub temperature_bracket_width_k: f64,
    pub hydraulic_solves: usize,
    /// Actual accepted midpoint profiles, ports and native work accounting.
    pub midpoint_hydraulics: CarreauGraphSolution,
}

/// A complete proposal. The borrowed initial owner is never mutated on refusal.
#[derive(Clone, Debug)]
pub struct ThermalCarreauStep<'a> {
    pub next: ThermalCarreauGraph<'a>,
    pub report: ThermalCarreauReport,
}

#[derive(Clone, Copy)]
struct Point {
    temperature: f64,
    residual: f64,
    allowance: f64,
}

struct Trial<H = CarreauGraphSolution> {
    point: Point,
    stored: f64,
    heat: f64,
    imposed: f64,
    bath: f64,
    hydraulics: H,
}

impl<'a> ThermalCarreauGraph<'a> {
    pub fn new(
        graph: CarreauPipeGraph,
        storage: &'a TabulatedHeatStorage,
        time_s: f64,
    ) -> Result<Self, ThermalCarreauError> {
        if !nonnegative(time_s) {
            return Err(GraphError::InvalidInput("invalid Carreau thermal clock").into());
        }
        let span = shared_temperature_span(&graph, storage)?;
        Ok(Self {
            graph,
            storage,
            temperature_span_k: span,
            time_s,
        })
    }

    pub const fn graph(&self) -> &CarreauPipeGraph {
        &self.graph
    }
    pub const fn storage(&self) -> &TabulatedHeatStorage {
        self.storage
    }
    pub fn temperature_k(&self) -> f64 {
        self.graph.branches[0].pipe.temperature_k()
    }
    pub const fn temperature_span_k(&self) -> [f64; 2] {
        self.temperature_span_k
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }

    /// Advance held hydraulic/thermal loads with caloric implicit midpoint.
    /// Reservoir and injection streams must return within the explicitly shared
    /// inventory. This is not an open-stream or radial-momentum time integrator.
    /// Original source errors are terminal; only a complete next owner publishes.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &self,
        duration_s: f64,
        boundaries: &[FixedPressure],
        injections: &[FlowInjection],
        boundary: ThermalCarreauBoundary,
        controls: ThermalCarreauControls,
        cx: &Cx<'_>,
    ) -> Result<ThermalCarreauStep<'a>, ThermalCarreauError> {
        checkpoint(cx)?;
        controls.validate()?;
        boundary.validate()?;
        if !positive(duration_s) {
            return Err(
                GraphError::InvalidInput("Carreau thermal duration must be positive").into(),
            );
        }
        let next_time = finite(self.time_s + duration_s, "Carreau thermal clock")?;
        if next_time <= self.time_s {
            return Err(GraphError::InvalidInput("Carreau thermal clock lost resolution").into());
        }
        let initial_temperature = self.temperature_k();
        let mut candidate = self.graph.clone();
        let mut evaluations = 0;
        let mut solve = |temperature: f64| -> Result<Trial, ThermalCarreauError> {
            checkpoint(cx)?;
            // Reserve one complete pressure solve for physical endpoint admission.
            if evaluations >= controls.maximum_hydraulic_solves - 1 {
                return Err(ThermalCarreauError::HydraulicSolveBudget);
            }
            if !temperature.is_finite()
                || temperature < self.temperature_span_k[0]
                || temperature > self.temperature_span_k[1]
            {
                return Err(ThermalCarreauError::TemperatureDomain {
                    temperature_k: temperature,
                });
            }
            let midpoint = 0.5 * initial_temperature + 0.5 * temperature;
            let stored = self
                .storage
                .energy_change_j(initial_temperature, temperature)?;
            prepare_temperature(&self.graph, &mut candidate, midpoint, cx)?;
            evaluations += 1;
            let hydraulics =
                candidate.solve_with_injections(boundaries, injections, controls.hydraulics, cx)?;
            let heat = finite(
                duration_s * hydraulics.dissipated_power_w,
                "Carreau viscous heat",
            )?;
            if heat < 0.0 {
                return Err(ThermalCarreauError::EnergyBalance);
            }
            let imposed = finite(
                duration_s * boundary.prescribed_power_w,
                "Carreau imposed heat",
            )?;
            let bath = finite(
                duration_s * boundary.bath_power(midpoint)?,
                "Carreau bath heat",
            )?;
            let residual = sum(&[stored, -heat, -imposed, -bath])?;
            let scale = stored.abs().max(heat).max(imposed.abs()).max(bath.abs());
            let allowance = finite(
                controls.absolute_energy_tolerance_j + controls.relative_energy_tolerance * scale,
                "Carreau thermal energy tolerance",
            )?;
            checkpoint(cx)?;
            Ok(Trial {
                point: Point {
                    temperature,
                    residual,
                    allowance,
                },
                stored,
                heat,
                imposed,
                bath,
                hydraulics,
            })
        };

        let initial_slope = finite(
            self.storage.heat_capacity_j_k(initial_temperature)?
                + (0.5 * duration_s) * boundary.bath_conductance_w_k,
            "Carreau thermal initial slope",
        )?;
        let (accepted, width, initial_point) = solve_temperature(
            initial_temperature,
            self.temperature_span_k,
            controls.maximum_temperature_change_k,
            controls.temperature_tolerance_k,
            initial_slope,
            &mut solve,
        )?;
        drop(solve);
        if accepted.point.residual.abs() > accepted.point.allowance {
            return Err(ThermalCarreauError::EnergyBalance);
        }
        prepare_temperature(&self.graph, &mut candidate, accepted.point.temperature, cx)?;
        evaluations += 1;
        let endpoint =
            candidate.solve_with_injections(boundaries, injections, controls.hydraulics, cx)?;
        let endpoint_net = sum(&[
            endpoint.dissipated_power_w,
            boundary.prescribed_power_w,
            boundary.bath_power(accepted.point.temperature)?,
        ])?;
        let endpoint_heat = finite(duration_s * endpoint_net, "Carreau endpoint thermal power")?;
        if initial_point.residual < 0.0 && endpoint_heat < -accepted.point.allowance
            || initial_point.residual > 0.0 && endpoint_heat > accepted.point.allowance
        {
            return Err(ThermalCarreauError::EquilibriumCrossing);
        }
        let reservoir = finite(
            duration_s * accepted.hydraulics.boundary_power_w,
            "Carreau reservoir work",
        )?;
        let injection = finite(
            duration_s * accepted.hydraulics.injection_power_w,
            "Carreau injection work",
        )?;
        let constitutive = finite(
            duration_s * accepted.hydraulics.constitutive_residual_power_w,
            "Carreau constitutive residual work",
        )?;
        let continuity = finite(
            duration_s * accepted.hydraulics.node_residual_power_w,
            "Carreau continuity residual work",
        )?;
        let hydraulic_residual = sum(&[reservoir, injection, -accepted.heat])?;
        let combined = sum(&[
            accepted.stored,
            -reservoir,
            -injection,
            -accepted.imposed,
            -accepted.bath,
        ])?;
        checkpoint(cx)?;
        Ok(ThermalCarreauStep {
            next: Self {
                graph: candidate,
                storage: self.storage,
                temperature_span_k: self.temperature_span_k,
                time_s: next_time,
            },
            report: ThermalCarreauReport {
                duration_s,
                initial_temperature_k: initial_temperature,
                midpoint_temperature_k: 0.5 * initial_temperature
                    + 0.5 * accepted.point.temperature,
                final_temperature_k: accepted.point.temperature,
                stored_energy_change_j: accepted.stored,
                viscous_heat_j: accepted.heat,
                prescribed_heat_j: accepted.imposed,
                bath_heat_j: accepted.bath,
                reservoir_work_j: reservoir,
                injection_work_j: injection,
                constitutive_residual_work_j: constitutive,
                continuity_residual_work_j: continuity,
                hydraulic_work_residual_j: hydraulic_residual,
                thermal_residual_j: accepted.point.residual,
                combined_energy_residual_j: combined,
                thermal_residual_allowance_j: accepted.point.allowance,
                temperature_bracket_width_k: width,
                hydraulic_solves: evaluations,
                midpoint_hydraulics: accepted.hydraulics,
            },
        })
    }
}

fn solve_temperature<H>(
    initial_temperature: f64,
    temperature_span_k: [f64; 2],
    maximum_temperature_change_k: f64,
    temperature_tolerance_k: f64,
    initial_slope: f64,
    solve: &mut impl FnMut(f64) -> Result<Trial<H>, ThermalCarreauError>,
) -> Result<(Trial<H>, f64, Point), ThermalCarreauError> {
    let initial = solve(initial_temperature)?;
    let initial_point = initial.point;
    let (accepted, width) = if initial.point.residual == 0.0 {
        (initial, 0.0)
    } else {
        // Bracket only in the physical initial heating/cooling direction.
        // Keep scalar endpoints; native profile vectors need not be retained.
        drop(initial);
        let direction = -initial_point.residual.signum();
        let bound = if direction > 0.0 {
            (initial_temperature + maximum_temperature_change_k).min(temperature_span_k[1])
        } else {
            (initial_temperature - maximum_temperature_change_k).max(temperature_span_k[0])
        };
        let limit = (bound - initial_temperature).abs();
        if !positive(limit) {
            return Err(ThermalCarreauError::UnbracketedTemperature);
        }
        let scale = initial_slope;
        let mut distance = (initial_point.residual.abs() / scale).min(limit);
        if !positive(distance) {
            return Err(ThermalCarreauError::TemperatureResolution);
        }
        let outer = loop {
            let temperature = if distance == limit {
                bound
            } else {
                initial_temperature + direction * distance
            };
            if temperature == initial_temperature {
                return Err(ThermalCarreauError::TemperatureResolution);
            }
            let trial = solve(temperature)?;
            if trial.point.residual == 0.0
                || trial.point.residual.signum() != initial_point.residual.signum()
            {
                break trial;
            }
            if distance == limit {
                return Err(ThermalCarreauError::UnbracketedTemperature);
            }
            distance = (2.0 * distance).min(limit);
        };
        if outer.point.residual == 0.0 {
            (outer, 0.0)
        } else {
            let (mut low, mut high) = if direction > 0.0 {
                (initial_point, outer.point)
            } else {
                (outer.point, initial_point)
            };
            drop(outer);
            loop {
                let width = high.temperature - low.temperature;
                let a = low.residual.abs();
                let b = high.residual.abs();
                let fraction = if a >= b {
                    1.0 / (1.0 + b / a)
                } else {
                    let ratio = a / b;
                    ratio / (1.0 + ratio)
                };
                let fraction = if (0.1..=0.9).contains(&fraction) {
                    fraction
                } else {
                    0.5
                };
                let temperature = low.temperature + fraction * width;
                if temperature <= low.temperature || temperature >= high.temperature {
                    return Err(ThermalCarreauError::TemperatureResolution);
                }
                let trial = solve(temperature)?;
                if trial.point.residual == 0.0
                    || width <= temperature_tolerance_k
                        && trial.point.residual.abs() <= trial.point.allowance
                {
                    let accepted_width = if trial.point.residual == 0.0 {
                        0.0
                    } else {
                        width
                    };
                    break (trial, accepted_width);
                }
                if trial.point.residual.signum() == low.residual.signum() {
                    low = trial.point;
                } else {
                    high = trial.point;
                }
            }
        }
    };
    Ok((accepted, width, initial_point))
}

fn shared_temperature_span(
    graph: &CarreauPipeGraph,
    storage: &TabulatedHeatStorage,
) -> Result<[f64; 2], ThermalCarreauError> {
    // The hydraulic owner admits one original law and a nonempty topology,
    // but may carry distinct prescribed branch temperatures. This caloric
    // owner has one shared inventory and therefore requires a uniform state.
    let reference = graph.branches[0].pipe;
    if graph
        .branches
        .iter()
        .any(|branch| branch.pipe.temperature_k() != reference.temperature_k())
    {
        return Err(GraphError::InvalidInput(
            "one-inventory Carreau thermal graph requires uniform temperature",
        )
        .into());
    }
    let source_span = reference.law().temperature_shift().validity_k();
    let storage_span = storage.temperature_span_k();
    let span = [
        source_span[0].max(storage_span[0]),
        source_span[1].min(storage_span[1]),
    ];
    let temperature = reference.temperature_k();
    if span[0] >= span[1] || temperature < span[0] || temperature > span[1] {
        return Err(ThermalCarreauError::TemperatureDomain {
            temperature_k: temperature,
        });
    }
    storage.heat_capacity_j_k(temperature)?;
    Ok(span)
}

fn prepare_temperature(
    original: &CarreauPipeGraph,
    candidate: &mut CarreauPipeGraph,
    temperature_k: f64,
    cx: &Cx<'_>,
) -> Result<(), ThermalCarreauError> {
    for (index, (source, target)) in original
        .branches
        .iter()
        .zip(&mut candidate.branches)
        .enumerate()
    {
        checkpoint(cx)?;
        let pipe = source.pipe;
        target.pipe = LaminarCarreauPipe::new(
            *pipe.law(),
            temperature_k,
            pipe.radius_m(),
            pipe.length_m(),
            pipe.maximum_laminar_wall_rate_per_s(),
        )
        .map_err(|source| CarreauGraphError::Pipe {
            branch: index,
            source,
        })?;
    }
    checkpoint(cx)?;
    Ok(())
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
fn sum(values: &[f64]) -> Result<f64, GraphError> {
    let mut sum = Sum::default();
    for &value in values {
        sum.add(value)?;
    }
    sum.value()
}
