//! Adopt a full solid/air correction only after production-law reassembly.
//! Goal convergence, physical acceptance and failed attempts stay distinct.

use std::collections::BTreeMap;
use fs_conduction::{ConductionError, ConductionSolution, ThermalBoundary};
use fs_conduction::adjoint::LinearGoalSolveConfig;
use crate::conjugate::AirMarch;
use super::{AirPath, ConductionProblem, Cx, LinearAirMaximumSolve, LinearConfig,
    LinearGoalAnalysisConfig, Result, RobinFeedbackAnalysisConfig, SpectralAirMaximumSolve,
    SpectralMaximumControl, ThermalInterfaces, bad, finite, poll, solve, spectral};

/// Explicit physical acceptance, in addition to the original solution's
/// absolute residual threshold. No tolerance is inferred from the goal bound.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PhysicalCoolingGates {
    /// Relative whole-solid energy closure in [0,1].
    pub energy_relative_tolerance: f64,
    /// Largest allowed reference change after using freshly integrated walls, K.
    pub reference_tolerance_k: f64,
    /// Absolute watt floor for each branch and the Robin decomposition check.
    pub balance_tolerance_w: f64,
    /// Relative watt allowance in [0,1), using each branch's OWN heat-rate scale.
    pub balance_relative_tolerance: f64,
}
impl PhysicalCoolingGates {
    fn admit(self) -> Result<()> {
        if !(self.energy_relative_tolerance.is_finite()
            && (0.0..=1.0).contains(&self.energy_relative_tolerance))
            || !(self.reference_tolerance_k.is_finite() && self.reference_tolerance_k > 0.0)
            || !(self.balance_tolerance_w.is_finite() && self.balance_tolerance_w >= 0.0)
            || !(self.balance_relative_tolerance.is_finite()
                && (0.0..1.0).contains(&self.balance_relative_tolerance))
        { return Err(bad("invalid physical cooling energy, reference or watt gates")); }
        Ok(())
    }
    fn watts(self, scale: f64) -> Result<f64> {
        Ok(self.balance_tolerance_w.max(finite(self.balance_relative_tolerance * scale)?))
    }
}

/// Why a computed candidate was not adopted. Original caller data are intact.
#[derive(Debug, Clone, PartialEq)]
pub enum PhysicalCoolingRefusal {
    /// Independent physical solid assembly, residual, material or energy refusal.
    Solid(ConductionError),
    /// Reintegrated walls changed a branch reference beyond its kelvin gate.
    Reference { branch: usize, delta_k: f64, limit_k: f64 },
    /// At least one segment disagreed with its independently integrated Robin flux.
    Interface { branch: usize, imbalance_w: f64, limit_w: f64 },
    /// The branch's total heat did not match its inlet/outlet enthalpy change.
    Enthalpy { branch: usize, imbalance_w: f64, limit_w: f64 },
    /// All Robin parts, INCLUDING non-air boundaries, disagree with the domain total.
    Decomposition { imbalance_w: f64, limit_w: f64 },
}

/// The separately scaled physical checks for one independently supplied inlet.
#[derive(Debug, Clone, PartialEq)]
pub struct PhysicalBranchCheck {
    /// Maximum reference change after fresh wall integration, K.
    pub reference_delta_k: f64,
    /// Maximum per-segment solid/air heat-rate difference, W.
    pub interface_imbalance_w: f64,
    /// Difference between total air heat and inlet/outlet enthalpy change, W.
    pub enthalpy_imbalance_w: f64,
    /// Admitted watt error for THIS branch, never another branch's scale.
    pub watt_limit: f64,
}

/// A physically admitted field, exact Robin boundary and freshly marched air.
#[derive(Debug, Clone, PartialEq)]
pub struct AcceptedLinearCooling {
    /// Physical residual, energy and contact/Robin fluxes belong to this field.
    /// Original iteration histories remain historical; see `correction` for work.
    pub solid: ConductionSolution,
    /// Exact preserved partition with the references used in physical reassembly.
    pub boundary: ThermalBoundary,
    /// Freshly integrated wall means in branch-major stream-wise order.
    pub wall_temperatures_k: Vec<f64>,
    /// Production exponential-law marches from those walls, in input path order.
    pub air: Vec<AirMarch>,
    /// Separate per-branch gates, including each branch's own magnitude scale.
    pub branches: Vec<PhysicalBranchCheck>,
    /// Independently accumulated whole-domain versus per-region Robin difference.
    pub robin_decomposition_imbalance_w: f64,
    /// True only if this accepted field also meets the STORED affine goal bound.
    /// Physical gates do not turn that bound into an assembly/continuum theorem.
    pub stored_goal_met: bool,
    /// Whether the adopted nodal temperature differs from the original field.
    pub temperature_changed: bool,
}

/// Complete correction work survives a physical rejection. Only `accepted`
/// contains a publishable physical solution; `correction` is always a candidate.
#[derive(Debug, Clone, PartialEq)]
pub struct PhysicalAirMaximumPolish {
    /// Candidate, original check, proof preparation and all correction work.
    pub correction: SpectralAirMaximumSolve,
    /// Present only after all physical gates and the final cancellation check.
    pub accepted: Option<AcceptedLinearCooling>,
    /// Explicit rejection of the returned candidate when no state was admitted.
    pub physical_refusal: Option<PhysicalCoolingRefusal>,
    /// Actual physical reassembly attempts, including rejected goal-eligible
    /// fields and the final budget/stagnation candidate. At most goal_checks+1.
    pub physical_checks: usize,
    /// Failed physical attempts; a later accepted field does not erase this work.
    pub physical_rejections: usize,
}

type CandidateAcceptance = std::result::Result<AcceptedLinearCooling, PhysicalCoolingRefusal>;

fn refuse(cx: &Cx<'_>, refusal: PhysicalCoolingRefusal) -> Result<CandidateAcceptance> {
    poll(cx)?;
    Ok(Err(refusal))
}

/// Correct until the stored maximum goal AND physical publication gates pass.
///
/// A goal-eligible field is rebuilt with its PRODUCTION exponential-law air
/// references, then independently reassembled and checked against the original
/// watt-residual threshold, energy/reference tolerances, per-segment heat balance,
/// branch enthalpy and full Robin decomposition. A physical rejection continues
/// the SAME full-feedback FGMRES loop; it does not reset iterations, tighten an
/// unrelated tolerance, or repeat inverse/response preparation. In particular,
/// a loose maximum goal no longer stops before physical residuals are acceptable.
///
/// Only a field passing both requirements terminates with GoalTolerance. At
/// budget/stagnation, the best numerical candidate is physically checked once
/// more and may be returned as accepted with stored_goal_met=false. Rejected
/// attempts are counted even when a later field passes. No caller state changes
/// on refusal or cancellation; a gate callback's provisional output never
/// escapes before the complete driver and final cancellation check succeed.
///
/// Each physical check uses one ordinary assembly/report pass, not a new solve.
/// There are at most correction.solution.solid.goal_checks+1 physical attempts.
/// Source, material, geometry and flow stay fixed; nonlinear conductivity and
/// radiation are unsupported. Stored-affine bounds are not coefficient-lowering,
/// continuum or experimental-validation evidence.
///
/// # Errors
/// Invalid controls/baseline, ordinary solve or air-law failures, nonfinite
/// arithmetic/allocation and cancellation. Final physical refusal is a result
/// value with all correction work retained; it is never numerical goal success.
#[allow(clippy::too_many_arguments)]
pub fn polish_linear_maximum_with_spectral(
    cx: &Cx<'_>, problem: ConductionProblem<'_>, interfaces: Option<&ThermalInterfaces>,
    paths: &[AirPath], linear: LinearConfig, original: &ConductionSolution,
    vertices: &[usize], solid_config: LinearGoalAnalysisConfig,
    feedback_config: RobinFeedbackAnalysisConfig, control: LinearGoalSolveConfig,
    spectral_control: SpectralMaximumControl, gates: PhysicalCoolingGates,
) -> Result<PhysicalAirMaximumPolish> {
    poll(cx)?;
    gates.admit()?;
    let threshold = original.report.residual_threshold;
    if !(threshold.is_finite() && threshold >= 0.0)
        || original.report.elements != problem.mesh.element_count()
        || original.temperature.len() != problem.mesh.vertex_count()
    { return Err(bad("original physical solution has an invalid residual threshold or mesh shape")); }
    let spectral::PreparedMaximum { analyzer, initial_analysis, preparation } = spectral::prepare(
        cx, problem, interfaces, paths, linear, &original.temperature, vertices,
        solid_config, feedback_config, control, spectral_control,
    )?;
    let mut accepted = None;
    let mut physical_refusal = None;
    let mut physical_checks = 0_usize;
    let mut physical_rejections = 0_usize;
    let solid = analyzer.solve_maximum_to_goal_admitted(
        cx, &original.temperature, vertices, control, |temperature, analysis| -> Result<bool> {
            let (_, candidate_air) = solve::air_from_temperature(cx, &analyzer, paths, temperature)?;
            physical_checks = physical_checks.checked_add(1).ok_or_else(|| bad("physical check count overflow"))?;
            match revalidate_candidate(cx, problem, interfaces, paths, original, temperature,
                &candidate_air, analysis.meets_absolute_tolerance(control.absolute_tolerance), gates)?
            {
                Ok(candidate) => {
                    accepted = Some(candidate);
                    physical_refusal = None;
                    Ok(true)
                }
                Err(refusal) => {
                    physical_rejections += 1;
                    physical_refusal = Some(refusal);
                    Ok(false)
                }
            }
        },
    )?;
    let (wall_temperatures_k, air) = solve::air_from_temperature(cx, &analyzer, paths, &solid.temperature)?;
    // Non-success returns the best NUMERICAL candidate. Its physical check may
    // differ from the last rejected trial, so never attach that trial's refusal.
    if accepted.is_none() {
        physical_checks = physical_checks.checked_add(1).ok_or_else(|| bad("physical check count overflow"))?;
        match revalidate_candidate(cx, problem, interfaces, paths, original, &solid.temperature,
            &air, solid.analysis.meets_absolute_tolerance(control.absolute_tolerance), gates)?
        {
            Ok(candidate) => { accepted = Some(candidate); physical_refusal = None; }
            Err(refusal) => { physical_rejections += 1; physical_refusal = Some(refusal); }
        }
    }
    let correction = SpectralAirMaximumSolve {
        solution: LinearAirMaximumSolve { solid, wall_temperatures_k, air },
        initial_analysis, preparation,
    };
    poll(cx)?;
    Ok(PhysicalAirMaximumPolish { correction, accepted, physical_refusal,
        physical_checks, physical_rejections })
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn revalidate_candidate(
    cx: &Cx<'_>, problem: ConductionProblem<'_>, interfaces: Option<&ThermalInterfaces>,
    paths: &[AirPath], original: &ConductionSolution, temperature: &[f64],
    candidate_air: &[AirMarch], goal_met: bool, gates: PhysicalCoolingGates,
) -> Result<CandidateAcceptance> {
    let references: Vec<_> = candidate_air.iter().flat_map(|air| &air.segments)
        .map(|row| (row.region.as_str(), row.reference_temperature_k)).collect();
    let port_count = references.len();
    let physical = original.revalidate_linear_robin_temperature(cx, problem, interfaces,
        temperature, &references, original.report.residual_threshold, gates.energy_relative_tolerance);
    let (solid, boundary) = match physical {
        Ok(physical) => physical,
        Err(error @ ConductionError::Cancelled { .. }) => return Err(error.into()),
        Err(error) => return refuse(cx, PhysicalCoolingRefusal::Solid(error)),
    };
    let mut fluxes = BTreeMap::new();
    let mut robin_total = 0.0;
    for row in &solid.report.robin_fluxes {
        poll(cx)?;
        if fluxes.insert(row.region.as_str(), row).is_some() {
            return Err(bad("physical Robin report repeats a region"));
        }
        robin_total = finite(robin_total + row.heat_rate_w)?;
    }
    let mut walls = Vec::new();
    walls.try_reserve_exact(port_count).map_err(|_| bad("physical wall allocation refused"))?;
    let mut air = Vec::new();
    air.try_reserve_exact(paths.len()).map_err(|_| bad("physical air allocation refused"))?;
    let mut branches = Vec::new();
    branches.try_reserve_exact(paths.len()).map_err(|_| bad("physical branch allocation refused"))?;
    for (branch, path) in paths.iter().enumerate() {
        poll(cx)?;
        let start = walls.len();
        for segment in path.segments() {
            poll(cx)?;
            let flux = fluxes.get(segment.region()).ok_or_else(|| bad("missing physical air-region flux"))?;
            walls.push(finite(flux.mean_wall_temperature_k)?);
        }
        let marched = path.march(&walls[start..])?;
        let mut reference_delta = 0.0_f64;
        let mut imbalance = 0.0_f64;
        let mut scale = 0.0_f64;
        for row in &marched.segments {
            poll(cx)?;
            let flux = fluxes.get(row.region.as_str()).ok_or_else(|| bad("missing marched-region flux"))?;
            reference_delta = reference_delta.max(finite(row.reference_temperature_k - flux.mean_reference_temperature_k)?.abs());
            imbalance = imbalance.max(finite(row.heat_rate_w - flux.heat_rate_w)?.abs());
            scale = scale.max(row.heat_rate_w.abs()).max(flux.heat_rate_w.abs());
        }
        let limit = gates.watts(scale)?;
        let enthalpy = finite(path.capacity_rate_w_per_k()
            * finite(marched.outlet_temperature_k - path.inlet_temperature_k())?)?;
        let enthalpy_error = finite(enthalpy - marched.total_heat_rate_w)?.abs();
        if reference_delta > gates.reference_tolerance_k {
            return refuse(cx, PhysicalCoolingRefusal::Reference {
                branch, delta_k: reference_delta, limit_k: gates.reference_tolerance_k });
        }
        if imbalance > limit {
            return refuse(cx, PhysicalCoolingRefusal::Interface { branch, imbalance_w: imbalance, limit_w: limit });
        }
        if enthalpy_error > limit {
            return refuse(cx, PhysicalCoolingRefusal::Enthalpy { branch, imbalance_w: enthalpy_error, limit_w: limit });
        }
        branches.push(PhysicalBranchCheck { reference_delta_k: reference_delta,
            interface_imbalance_w: imbalance, enthalpy_imbalance_w: enthalpy_error, watt_limit: limit });
        air.push(marched);
    }
    let decomposition = finite(robin_total - solid.report.energy.robin_out_w)?.abs();
    let limit = gates.watts(robin_total.abs().max(solid.report.energy.robin_out_w.abs()))?;
    if decomposition > limit {
        return refuse(cx, PhysicalCoolingRefusal::Decomposition { imbalance_w: decomposition, limit_w: limit });
    }
    let changed = solid.temperature.iter().zip(&original.temperature).any(|(a, b)| a.to_bits() != b.to_bits());
    poll(cx)?;
    Ok(Ok(AcceptedLinearCooling {
        solid, boundary, wall_temperatures_k: walls, air, branches,
        robin_decomposition_imbalance_w: decomposition, stored_goal_met: goal_met,
        temperature_changed: changed,
    }))
}
