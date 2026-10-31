//! Pressure-driven material/thermal operating points on an existing liquid path.
//!
//! Re-solve the existing temperature-dependent heat path at each mass-flow
//! trial. Hydraulic rho/mu are evaluated at the arithmetic BULK mean, not the
//! thermal film temperature. Each selected material claim and nonthermal query
//! coordinate is retained. The driving pressure is a differential pressure;
//! the material's thermodynamic pressure coordinate stays at its declared
//! reference. This is a low-compressibility, regional approximation, not an
//! evolving pressure EOS, boiling, momentum-acceleration or pump model.
//! Opt-in `solve_pressure_gravity_driven` adds the hydrostatic head from those
//! same bulk source densities at EVERY trial, including zero-pump circulation.
//!
//! The default pressure path does not include viscous heating. That omission requires
//! an explicit absolute hydraulic-power allowance; excessive omitted power
//! refuses. No global energy certificate is inferred from pressure closure.
//! A user-supplied positive mass-flow bracket is mandatory. All trials must be
//! source/correlation supported; errors never trigger a different material law.
//! The opt-in [`heating`] path instead integrates the original pressure loss
//! into the same fluid enthalpy trajectory and requires zero omitted power.

pub mod solid;
/// Local pressure-loss heating on the original temperature-feedback trajectory.
pub mod heating;
/// Parallel material-driven paths with enthalpy-conserving outlet mixing.
pub mod parallel;
/// Source-bound continuum permeability in the same native liquid solve.
pub mod porous;
use porous::{LiquidResistanceModel, PorousLiquidReceipt};

use fs_exec::Cx;
use fs_material::liquid::ResolvedLiquidState;
use fs_qty::{MassFlowRate, Pressure, Temperature};

use crate::graph::thermal::medium::series::hydraulic::{
    MaterialFlowResistance, ResistanceError, ResistanceEvaluation,
};
use crate::graph::thermal::medium::series::hydraulic::pressure::{
    HydrostaticPressure, SeriesGravity, SeriesPressureError,
};
use super::{LiquidFeedbackBudget, LiquidFeedbackError, LiquidHeatMarch,
    LiquidHeatPath, LiquidHeatSection, LiquidEnthalpyError, checkpoint, positive};

/// Exact pressure-loss assignment in thermal-path order. Geometry coefficients
/// are authored/calibrated assumptions, not another material-property table.
#[derive(Clone, Copy, Debug)]
pub struct LiquidResistanceSection<'a, R = MaterialFlowResistance> {
    /// Must match the corresponding thermal region exactly.
    pub region: &'a str,
    /// Current-material pressure-loss law and explicit Reynolds support.
    pub resistance: R,
}

/// Bounded positive-flow pressure solve, distinct from local thermal tolerances.
#[derive(Clone, Copy, Debug)]
pub struct PressureDrivenBudget {
    /// Strictly increasing, positive mass-flow bounds [kg/s]. No auto expansion.
    pub mass_flow_bracket_kg_per_s: (f64, f64),
    /// Pressure closure absolute floor [Pa], strictly positive.
    pub pressure_absolute_tolerance_pa: f64,
    /// Relative pressure closure tolerance in [0,1), zero disables this term.
    pub pressure_relative_tolerance: f64,
    /// Total thermal-path evaluations, including endpoints, in 2..=128.
    pub max_flow_evaluations: usize,
    /// Explicit maximum omitted hydraulic dissipation [W], finite and >= 0.
    /// Every evaluated trial must fit, not only the eventual operating point.
    pub maximum_neglected_hydraulic_power_w: f64,
}
impl PressureDrivenBudget {
    fn admit(self, drive: Pressure) -> Result<(), PressureDrivenError> {
        if !positive(drive.value()) {
            return Err(PressureDrivenError::Invalid("positive imposed pressure required without gravity"));
        }
        self.admit_with_gravity(drive)
    }

    fn admit_with_gravity(self, drive: Pressure) -> Result<(), PressureDrivenError> {
        let (lo,hi) = self.mass_flow_bracket_kg_per_s;
        if !positive(lo) || !positive(hi) || lo >= hi || !drive.value().is_finite()
            || !positive(self.pressure_absolute_tolerance_pa)
            || !(0.0..1.0).contains(&self.pressure_relative_tolerance)
            || !(2..=128).contains(&self.max_flow_evaluations)
            || !self.maximum_neglected_hydraulic_power_w.is_finite()
            || self.maximum_neglected_hydraulic_power_w < 0.0
        { return Err(PressureDrivenError::Invalid("invalid pressure drive, bracket or budget")); }
        Ok(())
    }
}

/// A refused operating point never contains a partially accepted thermal state.
#[derive(Debug)]
pub enum PressureDrivenError {
    /// Invalid pressure/geometry assignment or nonfinite result.
    Invalid(&'static str),
    /// Existing material/caloric/correlation/thermal refusal, preserved as-is.
    Thermal(LiquidFeedbackError),
    /// Original source-density hydrostatic geometry/arithmetic refusal.
    Gravity(SeriesPressureError),
    /// Existing pressure-law/material refusal, with exact region attribution.
    Resistance { region: String, source: ResistanceError },
    /// Both endpoint residuals have the same sign. Multiple roots are not ruled
    /// out, but this explicitly supplied bracket cannot select one by bisection.
    NotBracketed { lower_residual_pa: f64, upper_residual_pa: f64 },
    /// No pressure-closed state before the evaluation budget/float resolution.
    NotConverged { evaluations: usize, last_residual_pa: f64 },
    /// Shared solid-callback limit exhausted across all flow trials.
    SolidSolveBudget { performed: usize, maximum: usize },
    /// Neglecting viscous/inertial heating exceeds the caller's explicit budget.
    NeglectedDissipation { power_w: f64, allowance_w: f64 },
}
impl core::fmt::Display for PressureDrivenError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "pressure-driven material flow: {self:?}")
    }
}
impl std::error::Error for PressureDrivenError {}
impl From<LiquidFeedbackError> for PressureDrivenError {
    fn from(value: LiquidFeedbackError) -> Self { Self::Thermal(value) }
}
impl From<LiquidEnthalpyError> for PressureDrivenError {
    fn from(value: LiquidEnthalpyError) -> Self { Self::Thermal(value.into()) }
}

/// Source-bound bulk state used by one hydraulic region.
#[derive(Clone, Debug, PartialEq)]
pub struct LiquidHydraulicRegion {
    /// Matching thermal trace.
    pub region: String,
    /// Source-resolved bulk-mean state and all four property-use receipts.
    pub material: ResolvedLiquidState,
    /// Original shape, material snapshot, pressure, flow, tangent and power.
    pub loss: ResistanceEvaluation,
    /// Original porous skeleton and actual pressure faces, when present.
    pub porous: Option<PorousLiquidReceipt>,
}

/// Pressure and dissipation evaluated on the SAME accepted coolant trajectory.
#[derive(Clone, Debug, PartialEq)]
pub struct LiquidPressureState {
    /// Source-bound pressure contributions in stream order.
    pub regions: Vec<LiquidHydraulicRegion>,
    /// Sum of forward pressure drops [Pa].
    pub total_pressure_drop_pa: f64,
    /// Sum of local dp*Q [W], not dp_total times an unrelated reference Q.
    pub total_dissipated_power_w: f64,
}

/// Operating point with its original thermal or coupled-solid payload retained.
#[derive(Debug)]
pub struct PressureDrivenSolution<T> {
    /// Exact thermal solve evaluated at `mass_flow_kg_per_s`.
    pub thermal: T,
    /// Exact bulk-material pressure evaluations for that same solve.
    pub hydraulics: LiquidPressureState,
    /// One conserved mass-flow rate found by the pressure solve [kg/s].
    pub mass_flow_kg_per_s: f64,
    /// Original prescribed pressure difference [Pa].
    pub driving_pressure_pa: f64,
    /// Signed sum(dp_loss)-drive-hydrostatic_gain residual [Pa].
    pub pressure_residual_pa: f64,
    /// Actual full thermal-path evaluations, including endpoints.
    pub flow_evaluations: usize,
    /// Explicit geometry/gravity used in the accepted thermal-pressure solve.
    pub gravity: Option<SeriesGravity>,
    /// Hydrostatic head at that solve's actual bulk densities, when requested.
    pub hydrostatic: Option<HydrostaticPressure>,
}

impl<'a, 'card> LiquidHeatPath<'a, 'card> {
    /// Determine mass flow from differential pressure and actual wall temperatures.
    /// Reuses `try_new` and `march` at every trial, including all source pins,
    /// local coefficient/outlet convergence and sensible-enthalpy transport.
    /// The template's stored mass flow is NOT silently used as the solution.
    ///
    /// The bracket isolates a candidate, not a uniqueness or stability proof.
    /// Returned hydraulic states are bulk-mean evaluations even when heat
    /// transfer uses FilmMean. All references, cards and geometry are immutable.
    ///
    /// # Errors
    /// Invalid configuration/assignment, unsupported trials, cancellation,
    /// excessive neglected dissipation, missing bracket or pressure nonconvergence.
    pub fn solve_pressure_driven(
        &self, cx: &Cx<'_>, walls: &[f64], resistances: &[LiquidResistanceSection<'_, impl LiquidResistanceModel>],
        drive: Pressure, local: LiquidFeedbackBudget, budget: PressureDrivenBudget,
    ) -> Result<PressureDrivenSolution<LiquidHeatMarch<'a>>, PressureDrivenError> {
        self.solve_pressure_driven_inner(cx, walls, resistances, drive, None, local, budget)
    }

    /// Close source-dependent thermal transport and buoyant circulation together.
    /// Each trial marches the original enthalpy path and reevaluates rho(T), mu(T)
    /// before balancing imposed pressure PLUS the connected path's gravity head.
    /// Zero imposed drive is supported; a signed imposed drive can oppose gravity.
    ///
    /// The positive-flow bracket selects one orientation; reverse circulation
    /// requires reversing the authored thermal path and its inlet. This is a
    /// quasi-steady single-phase loop with a prescribed inlet thermal reservoir,
    /// not a closed-inventory transient or a stability/uniqueness proof. Regional
    /// bulk density is an approximation; refine the thermal regions to resolve
    /// a varying profile. Existing hydraulic-heating omission budgets still apply.
    /// Pressure-dependent porous laws are not admitted by this overload because
    /// their local absolute-pressure faces need a separate gravity-aware model.
    pub fn solve_pressure_gravity_driven(
        &self, cx: &Cx<'_>, walls: &[f64], resistances: &[LiquidResistanceSection<'_>],
        drive: Pressure, gravity: &SeriesGravity, local: LiquidFeedbackBudget,
        budget: PressureDrivenBudget,
    ) -> Result<PressureDrivenSolution<LiquidHeatMarch<'a>>, PressureDrivenError> {
        self.solve_pressure_driven_inner(cx, walls, resistances, drive, Some(gravity), local, budget)
    }

    fn solve_pressure_driven_inner(
        &self, cx: &Cx<'_>, walls: &[f64], resistances: &[LiquidResistanceSection<'_, impl LiquidResistanceModel>],
        drive: Pressure, gravity: Option<&SeriesGravity>, local: LiquidFeedbackBudget, budget: PressureDrivenBudget,
    ) -> Result<PressureDrivenSolution<LiquidHeatMarch<'a>>, PressureDrivenError> {
        checkpoint(cx)?;
        if let Some(path) = gravity {
            budget.admit_with_gravity(drive)?;
            if path.section_count() != self.sections.len() {
                return Err(PressureDrivenError::Invalid("one gravity segment per thermal section required"));
            }
        } else { budget.admit(drive)?; }
        local.admit()?;
        self.admit_resistances(resistances)?;
        if walls.len() != self.sections.len() {
            return Err(PressureDrivenError::Invalid("one wall temperature per hydraulic/thermal section is required"));
        }
        for &wall in walls { super::state_at_temperature(self.source.curve(), wall)?; }
        solve_pressure_root_with_gravity(cx, drive, gravity, budget, |mass| {
            let path = self.at_pressure_trial_mass(cx, mass)?;
            let thermal = path.march(cx, walls, local)?;
            let pressure = self.pressure_state(cx, resistances, &thermal, budget)?;
            Ok((thermal, pressure))
        })
    }

    fn admit_resistances(&self, resistances: &[LiquidResistanceSection<'_, impl LiquidResistanceModel>]) -> Result<(), PressureDrivenError> {
        porous::path_outlet_pressure(resistances)?;
        if resistances.len() != self.sections.len() {
            return Err(PressureDrivenError::Invalid("one resistance must be assigned per thermal section"));
        }
        for (hydraulic, thermal) in resistances.iter().zip(&self.sections) {
            if hydraulic.region != thermal.name
                || hydraulic.resistance.flow_area() != thermal.geometry.flow_area
            {
                return Err(PressureDrivenError::Invalid("hydraulic region/flow area differs from the thermal path"));
            }
        }
        Ok(())
    }

    fn at_pressure_trial_mass(&self, cx: &Cx<'_>, mass: f64) -> Result<Self, PressureDrivenError> {
        let sections: Vec<_> = self.sections.iter().map(|s| LiquidHeatSection {
            region: &s.name, geometry: s.geometry, correlation: s.correlation,
        }).collect();
        Ok(Self::try_new(cx, self.law, self.source, MassFlowRate::new(mass),
            Temperature::new(self.inlet.temperature_k()), &sections, self.property_temperature)?)
    }

    fn pressure_state(
        &self, cx: &Cx<'_>, resistances: &[LiquidResistanceSection<'_, impl LiquidResistanceModel>],
        thermal: &LiquidHeatMarch<'a>, budget: PressureDrivenBudget,
    ) -> Result<LiquidPressureState, PressureDrivenError> {
        if thermal.regions.len() != resistances.len() {
            return Err(PressureDrivenError::Invalid("incomplete thermal trajectory"));
        }
        let mut regions = Vec::with_capacity(resistances.len());
        let mut total_pressure = 0.0;
        let mut total_power = 0.0;
        for (section, flow) in resistances.iter().zip(&thermal.regions) {
            checkpoint(cx)?;
            if section.region != flow.region { return Err(PressureDrivenError::Invalid("pressure/heat region mismatch")); }
            let bulk = 0.5*flow.inlet.temperature_k() + 0.5*flow.outlet.temperature_k();
            let material = self.law.at_temperature(Temperature::new(bulk)).map_err(LiquidFeedbackError::from)?;
            let loss = section.resistance.evaluate_material(&material, bulk, MassFlowRate::new(thermal.mass_flow_kg_per_s))
                .map_err(|source| PressureDrivenError::Resistance { region: section.region.to_owned(), source })?;
            total_pressure += loss.pressure_drop_pa;
            total_power += loss.dissipated_power_w;
            if !positive(total_pressure) || !total_power.is_finite() {
                return Err(PressureDrivenError::Invalid("aggregate pressure/power is unrepresentable"));
            }
            regions.push(LiquidHydraulicRegion { region: section.region.to_owned(), material, loss, porous: None });
        }
        if total_power > budget.maximum_neglected_hydraulic_power_w {
            return Err(PressureDrivenError::NeglectedDissipation { power_w: total_power,
                allowance_w: budget.maximum_neglected_hydraulic_power_w });
        }
        self.admit_porous_trajectory(cx, resistances, &thermal.regions, &mut regions)?;
        checkpoint(cx)?;
        Ok(LiquidPressureState { regions, total_pressure_drop_pa: total_pressure, total_dissipated_power_w: total_power })
    }
}

// One bounded operating-point algorithm for fixed-wall and coupled-solid paths.
// Never use the fixed-material tangent as a thermal-feedback derivative.
fn solve_pressure_root<T>(
    cx: &Cx<'_>, drive: Pressure, budget: PressureDrivenBudget,
    evaluate: impl FnMut(f64) -> Result<(T, LiquidPressureState), PressureDrivenError>,
) -> Result<PressureDrivenSolution<T>, PressureDrivenError> {
    solve_pressure_root_with_gravity(cx, drive, None, budget, evaluate)
}

fn solve_pressure_root_with_gravity<T>(
    cx: &Cx<'_>, drive: Pressure, gravity: Option<&SeriesGravity>, budget: PressureDrivenBudget,
    mut evaluate: impl FnMut(f64) -> Result<(T, LiquidPressureState), PressureDrivenError>,
) -> Result<PressureDrivenSolution<T>, PressureDrivenError> {
    if gravity.is_some() { budget.admit_with_gravity(drive)?; }
    else { budget.admit(drive)?; }
    let mut calls = 0;
    let mut probe = |mass| -> Result<PressureDrivenSolution<T>, PressureDrivenError> {
        checkpoint(cx)?;
        let (thermal, hydraulics) = evaluate(mass)?;
        calls += 1;
        let hydrostatic = gravity.map(|path| path.evaluate_densities(
            hydraulics.regions.iter().map(|region| region.loss.transport.density_kg_m3()), mass,
            || cx.checkpoint().map_err(|_| SeriesPressureError::Cancelled),
        )).transpose().map_err(PressureDrivenError::Gravity)?;
        let effective_drive = drive.value()
            + hydrostatic.as_ref().map_or(0.0, HydrostaticPressure::pressure_gain_pa);
        let residual = hydraulics.total_pressure_drop_pa-effective_drive;
        if !residual.is_finite() { return Err(PressureDrivenError::Invalid("pressure residual is nonfinite")); }
        checkpoint(cx)?;
        Ok(PressureDrivenSolution { thermal, hydraulics, mass_flow_kg_per_s: mass,
            driving_pressure_pa: drive.value(), pressure_residual_pa: residual, flow_evaluations: calls,
            gravity: gravity.cloned(), hydrostatic })
    };
    let closed = |value: &PressureDrivenSolution<T>| {
        value.pressure_residual_pa.abs() <= budget.pressure_absolute_tolerance_pa.max(
            budget.pressure_relative_tolerance*(drive.value()
                + value.hydrostatic.as_ref().map_or(0.0, HydrostaticPressure::pressure_gain_pa))
                .abs().max(value.hydraulics.total_pressure_drop_pa.abs()))
    };
    let (mut low,mut high) = budget.mass_flow_bracket_kg_per_s;
    let a = probe(low)?;
    if closed(&a) { return Ok(a); }
    let mut low_residual = a.pressure_residual_pa;
    drop(a);
    let b = probe(high)?;
    if closed(&b) { return Ok(b); }
    let mut last_residual = b.pressure_residual_pa;
    if low_residual.is_sign_positive() == last_residual.is_sign_positive() {
        return Err(PressureDrivenError::NotBracketed { lower_residual_pa: low_residual, upper_residual_pa: last_residual });
    }
    drop(b);
    let mut evaluations = 2;
    for _ in 2..budget.max_flow_evaluations {
        let mass = 0.5*low+0.5*high;
        if mass == low || mass == high { break; }
        let value = probe(mass)?;
        evaluations = value.flow_evaluations;
        if closed(&value) { return Ok(value); }
        last_residual = value.pressure_residual_pa;
        if last_residual.is_sign_positive() == low_residual.is_sign_positive() {
            low = mass; low_residual = last_residual;
        } else { high = mass; }
    }
    Err(PressureDrivenError::NotConverged { evaluations, last_residual_pa: last_residual })
}
