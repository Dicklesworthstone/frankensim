//! Pump head curves coupled to the actual material-dependent cooling network.
//!
//! A fixed-speed head/volume-flow curve supplies dp = rho_in * g * H(Q_in).
//! Q_in uses the resolved COMMON INLET density, never a hot branch density.
//! Every pressure trial re-solves all branch material/enthalpy operating points.
//! The same matcher is used for fixed walls and common-solid feedback.
//!
//! This is a nominal, low-compressibility operating point, not a uniqueness,
//! stability, cavitation, efficiency or pump-viscosity-correction certificate.
//! The caller must supply a curve applicable to the declared inlet viscosity
//! range; changing viscosity beyond that range refuses. No affinity scaling or
//! motor heat is invented. Head data are not material properties.

/// Joint pump operating point and common-solid thermal feedback.
pub mod solid;

use fs_exec::Cx;
use fs_material::liquid::ResolvedLiquidState;
use fs_qty::{Length, Pressure, VolumetricFlowRate};
use crate::SourceProvenance;
use super::LiquidResistanceModel;
use super::{LiquidFeedbackBudget, ParallelLiquidNetwork, ParallelLiquidSolution,
    ParallelPressureBudget, ParallelPressureError, checkpoint, positive};

/// One nominal sample at the curve's declared, fixed operating speed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PumpHeadPoint {
    /// Volume flow at the pump inlet, m3/s.
    pub flow: VolumetricFlowRate,
    /// Hydraulic head, m; not a pressure or an efficiency.
    pub head: Length,
}

/// Immutable, bounded, piecewise-linear equipment data with explicit applicability.
#[derive(Clone, Debug, PartialEq)]
pub struct PumpHeadCurve {
    points: Vec<PumpHeadPoint>,
    gravity_m_per_s2: f64,
    kinematic_viscosity_domain_m2_s: (f64, f64),
    source: SourceProvenance,
}
impl PumpHeadCurve {
    /// Admit 2..=1024 increasing flows, nonincreasing nonnegative heads, positive
    /// gravity and a finite positive inclusive inlet-kinematic-viscosity domain.
    /// Samples outside the supplied flow span are never extrapolated. The source
    /// identifies the fixed-speed curve; metadata is not experimental validation.
    pub fn try_new(
        points: Vec<PumpHeadPoint>, gravity_m_per_s2: f64,
        kinematic_viscosity_domain_m2_s: (f64, f64), source: SourceProvenance,
    ) -> Result<Self, PumpOperatingError> {
        let (lo, hi) = kinematic_viscosity_domain_m2_s;
        if !(2..=1024).contains(&points.len()) || !positive(gravity_m_per_s2)
            || !positive(lo) || !positive(hi) || lo > hi
            || source.citation.trim().is_empty() || source.identifier.trim().is_empty()
            || points.iter().any(|p| !p.flow.value().is_finite() || p.flow.value() < 0.0
                || !p.head.value().is_finite() || p.head.value() < 0.0)
            || points.windows(2).any(|p| p[1].flow.value() <= p[0].flow.value()
                || p[1].head.value() > p[0].head.value())
        {
            return Err(PumpOperatingError::Invalid("invalid pump curve, gravity, source or viscosity domain"));
        }
        Ok(Self { points, gravity_m_per_s2, kinematic_viscosity_domain_m2_s, source })
    }
    /// Original samples, including both admitted endpoint flows.
    pub fn points(&self) -> &[PumpHeadPoint] { &self.points }
    /// Original equipment-data citation; not a material receipt.
    pub const fn source(&self) -> &SourceProvenance { &self.source }
    /// Gravity used by the head convention, m/s2.
    pub const fn gravity_m_per_s2(&self) -> f64 { self.gravity_m_per_s2 }
    /// The explicitly supported inlet kinematic-viscosity range, m2/s.
    pub const fn kinematic_viscosity_domain_m2_s(&self) -> (f64, f64) {
        self.kinematic_viscosity_domain_m2_s
    }
    fn admit_inlet(&self, inlet: &ResolvedLiquidState) -> Result<(f64, f64), PumpOperatingError> {
        let density = inlet.density_kg_m3();
        let viscosity = inlet.dynamic_viscosity_pa_s() / density;
        let (lo, hi) = self.kinematic_viscosity_domain_m2_s;
        if !positive(density) || !positive(viscosity) || viscosity < lo || viscosity > hi {
            return Err(PumpOperatingError::OutsideCurve {
                field: "inlet kinematic viscosity", value: viscosity, lower: lo, upper: hi,
            });
        }
        Ok((density, viscosity))
    }
    /// Evaluate head and hydraulic power at the actual inlet material and flow.
    /// Shaft/electrical power and efficiency are deliberately not inferred.
    pub fn evaluate(
        &self, inlet: &ResolvedLiquidState, flow: VolumetricFlowRate,
    ) -> Result<PumpSupplyState, PumpOperatingError> {
        let q = flow.value();
        let (density, viscosity) = self.admit_inlet(inlet)?;
        let lo = self.points[0].flow.value();
        let hi = self.points[self.points.len()-1].flow.value();
        if !q.is_finite() || q < lo || q > hi {
            return Err(PumpOperatingError::OutsideCurve {
                field: "inlet volume flow", value: q, lower: lo, upper: hi,
            });
        }
        let upper = self.points.partition_point(|p| p.flow.value() < q);
        let head = if self.points[upper].flow.value() == q {
            self.points[upper].head.value()
        } else {
            let left = self.points[upper-1];
            let right = self.points[upper];
            let fraction = (q-left.flow.value()) / (right.flow.value()-left.flow.value());
            (1.0-fraction)*left.head.value() + fraction*right.head.value()
        };
        let pressure = density * (self.gravity_m_per_s2 * head);
        let power = pressure * q;
        if !pressure.is_finite() || !power.is_finite()
            || (head > 0.0 && pressure <= 0.0) || (pressure > 0.0 && q > 0.0 && power <= 0.0)
        { return Err(PumpOperatingError::Invalid("pump pressure or hydraulic power is unrepresentable")); }
        Ok(PumpSupplyState { inlet_volume_flow_m3_s: q, inlet_density_kg_m3: density,
            inlet_kinematic_viscosity_m2_s: viscosity, head_m: head,
            available_pressure_pa: pressure, hydraulic_power_w: power })
    }
}

/// Actual pump evaluation; hydraulic output is not shaft power or heat deposition.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PumpSupplyState {
    /// Actual supply volume flow, m3/s.
    pub inlet_volume_flow_m3_s: f64,
    /// Actual common-inlet density, kg/m3.
    pub inlet_density_kg_m3: f64,
    /// Admitted inlet kinematic viscosity, m2/s.
    pub inlet_kinematic_viscosity_m2_s: f64,
    /// Interpolated fixed-speed head, m.
    pub head_m: f64,
    /// Density-scaled pump pressure, Pa.
    pub available_pressure_pa: f64,
    /// Pressure times supply volume flow, W; not shaft/electrical power.
    pub hydraulic_power_w: f64,
}

/// Pressure bracket and absolute/relative residual budget for pump/network matching.
#[derive(Clone, Copy, Debug)]
pub struct PumpOperatingControls {
    /// Positive, increasing network pressure bounds, Pa. No automatic expansion.
    pub pressure_bracket_pa: (f64, f64),
    /// Positive absolute pressure-residual tolerance, Pa.
    pub absolute_pressure_tolerance_pa: f64,
    /// Relative pressure-residual tolerance in [0,1).
    pub relative_pressure_tolerance: f64,
    /// Complete network evaluations, including endpoints, in 2..=128.
    pub max_network_evaluations: usize,
}
impl PumpOperatingControls {
    fn admit(self) -> Result<(), PumpOperatingError> {
        let (lo, hi) = self.pressure_bracket_pa;
        if !positive(lo) || !positive(hi) || lo >= hi
            || !positive(self.absolute_pressure_tolerance_pa)
            || !(0.0..1.0).contains(&self.relative_pressure_tolerance)
            || !(2..=128).contains(&self.max_network_evaluations)
        { return Err(PumpOperatingError::Invalid("invalid pump pressure bracket or numerical controls")); }
        Ok(())
    }
}

/// A failed source or network publishes no candidate state and never falls back.
#[derive(Debug)]
pub enum PumpOperatingError {
    Invalid(&'static str),
    OutsideCurve { field: &'static str, value: f64, lower: f64, upper: f64 },
    Network(ParallelPressureError),
    NotBracketed { lower_residual_pa: f64, upper_residual_pa: f64 },
    NotConverged { evaluations: usize, last_residual_pa: f64 },
    SolidBudget { performed: usize, maximum: usize },
}
impl core::fmt::Display for PumpOperatingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material pump operating point: {self:?}")
    }
}
impl std::error::Error for PumpOperatingError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self { Self::Network(error) => Some(error), _ => None }
    }
}
impl From<ParallelPressureError> for PumpOperatingError {
    fn from(error: ParallelPressureError) -> Self { Self::Network(error) }
}

/// Retains the original source and the complete state actually used for matching.
#[derive(Debug)]
pub struct PumpOperatingPoint<'s, T> {
    /// Original immutable equipment source and applicability.
    pub curve: &'s PumpHeadCurve,
    /// Original accepted network or coupled-solid result.
    pub state: T,
    /// Pump evaluation at that same result's inlet flow and material.
    pub supply: PumpSupplyState,
    /// Common pressure used by the accepted network solve, Pa.
    pub network_pressure_pa: f64,
    /// Signed network pressure minus available pump pressure, Pa.
    pub pressure_residual_pa: f64,
    /// Complete network candidates evaluated, including endpoints.
    pub network_evaluations: usize,
    /// Across ALL pressure candidates, not just the final network result.
    pub total_flow_evaluations: usize,
    /// Common-solid callbacks across every pressure candidate; zero for fixed walls.
    pub total_solid_solves: usize,
}

impl<'p, 'a, 'card, H: LiquidResistanceModel> ParallelLiquidNetwork<'p, 'a, 'card, H> {
    fn admit_pump_inlet(&self, cx: &Cx<'_>, pump: &PumpHeadCurve) -> Result<(), PumpOperatingError> {
        checkpoint(cx).map_err(ParallelPressureError::from)?;
        let path = self.branches[0].path;
        let inlet = path.law.at_temperature(fs_qty::Temperature::new(path.inlet.temperature_k()))
            .map_err(super::LiquidFeedbackError::from).map_err(ParallelPressureError::from)?;
        pump.admit_inlet(&inlet)?;
        Ok(())
    }
    /// Intersect a fixed-speed pump head curve with temperature-dependent branches.
    /// Each trial calls the existing network calculation with one shared work
    /// meter. One branch is the same API for a series circuit. The returned
    /// pressure residual is checked, never inferred from a small pressure bracket.
    /// All branch flow brackets and original omitted-heating gates still apply.
    pub fn solve_with_pump<'s>(
        &self, cx: &Cx<'_>, walls: &[f64], pump: &'s PumpHeadCurve,
        local: LiquidFeedbackBudget, work: ParallelPressureBudget, controls: PumpOperatingControls,
    ) -> Result<PumpOperatingPoint<'s, ParallelLiquidSolution<'a>>, PumpOperatingError> {
        controls.admit()?;
        work.admit()?;
        self.admit_pump_inlet(cx, pump)?;
        let mut performed = 0;
        match_pump(cx, pump, controls, |drive| {
            Ok(self.march_metered(cx, walls, drive, local, work, &mut performed)?)
        })
    }
}

// The pressure matcher keeps the original payload instead of reconstructing a
// thermal result from detached flow/temperature scalars. Shared by both drivers.
trait PumpLoad {
    fn fluid(&self) -> &ParallelLiquidSolution<'_>;
    fn flow_evaluations(&self) -> usize;
    fn solid_solves(&self) -> usize { 0 }
}
impl PumpLoad for ParallelLiquidSolution<'_> {
    fn fluid(&self) -> &ParallelLiquidSolution<'_> { self }
    fn flow_evaluations(&self) -> usize { self.flow_evaluations }
}

fn match_pump<'s, T: PumpLoad>(
    cx: &Cx<'_>, pump: &'s PumpHeadCurve, controls: PumpOperatingControls,
    mut evaluate: impl FnMut(Pressure) -> Result<T, PumpOperatingError>,
) -> Result<PumpOperatingPoint<'s, T>, PumpOperatingError> {
    controls.admit()?;
    let mut total_flow = 0usize;
    let mut total_solid = 0usize;
    let mut probe = |pressure: f64, calls: usize| -> Result<_, PumpOperatingError> {
        checkpoint(cx).map_err(ParallelPressureError::from)?;
        let state = evaluate(Pressure::new(pressure))?;
        let fluid = state.fluid();
        let supply = pump.evaluate(&fluid.inlet_material, VolumetricFlowRate::new(fluid.inlet_volume_flow_m3_per_s))?;
        let residual = pressure-supply.available_pressure_pa;
        if !residual.is_finite() || fluid.driving_pressure_pa != pressure {
            return Err(PumpOperatingError::Invalid("pump matching received an inconsistent or nonfinite network state"));
        }
        total_flow = total_flow.checked_add(state.flow_evaluations())
            .ok_or(PumpOperatingError::Invalid("pump flow-evaluation count overflow"))?;
        total_solid = total_solid.checked_add(state.solid_solves())
            .ok_or(PumpOperatingError::Invalid("pump solid-evaluation count overflow"))?;
        checkpoint(cx).map_err(ParallelPressureError::from)?;
        Ok(PumpOperatingPoint { curve: pump, state, supply, network_pressure_pa: pressure,
            pressure_residual_pa: residual, network_evaluations: calls,
            total_flow_evaluations: total_flow, total_solid_solves: total_solid })
    };
    let closed = |p: &PumpOperatingPoint<'_, T>| p.pressure_residual_pa.abs()
        <= controls.absolute_pressure_tolerance_pa.max(controls.relative_pressure_tolerance
            * p.network_pressure_pa.max(p.supply.available_pressure_pa));
    let (mut low, mut high) = controls.pressure_bracket_pa;
    let lower = probe(low, 1)?;
    if closed(&lower) { return Ok(lower); }
    let mut low_residual = lower.pressure_residual_pa;
    drop(lower);
    let upper = probe(high, 2)?;
    if closed(&upper) { return Ok(upper); }
    let mut last_residual = upper.pressure_residual_pa;
    if last_residual.is_sign_positive() == low_residual.is_sign_positive() {
        return Err(PumpOperatingError::NotBracketed { lower_residual_pa: low_residual, upper_residual_pa: last_residual });
    }
    drop(upper);
    let mut evaluations = 2;
    for calls in 3..=controls.max_network_evaluations {
        let middle = 0.5*low + 0.5*high;
        if middle == low || middle == high { break; }
        let candidate = probe(middle, calls)?;
        evaluations = calls;
        if closed(&candidate) { return Ok(candidate); }
        last_residual = candidate.pressure_residual_pa;
        if last_residual.is_sign_positive() == low_residual.is_sign_positive() {
            low = middle; low_residual = last_residual;
        } else { high = middle; }
    }
    Err(PumpOperatingError::NotConverged { evaluations, last_residual_pa: last_residual })
}
