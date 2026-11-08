//! Original-source steady flow through a perfectly mixed chamber at a fixed wall.
//!
//! This reusable port owns no fluid or wall inventory. Its pressure, feed and
//! native PH policy remain fixed; every trial retains the original phase
//! calorics and computes source-integrated enthalpy and entropy increments.

use super::{
    MAX_HEAT_EXCHANGER_DUTY_EVALUATIONS, MAX_HEAT_EXCHANGER_PT_EVALUATIONS,
    MixtureHeatExchangerError as Error, MixtureHeatExchangerInlet, Work,
    phase_error, poll, positive,
};
use crate::chemistry::ThermalConductance;
use fs_material::thermochemical_gas::nonideal::flash::equilibrium::{
    PhaseEnthalpyControls, PhaseEnthalpySolution, prepared::PreparedPhaseEnthalpyPath,
};
use fs_qty::{Power, Temperature};

/// Original caloric policy and finite physical/work limits for one flowing port.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixtureWallStreamControl {
    /// Increasing signed heat powers into the stream [W], containing zero.
    /// The passive initial-flux bound further restricts this interval.
    pub duty_bracket_w: [f64; 2],
    /// Original PH domain, phase policy and per-query PT limit.
    pub equilibrium: PhaseEnthalpyControls,
    /// Absolute allowance for source-target and contact closure [W].
    pub absolute_power_tolerance_w: f64,
    /// Relative allowance scaled by actual heat power.
    pub relative_power_tolerance: f64,
    /// Numerical allowance for negative stream-plus-wall entropy [W/K].
    pub entropy_tolerance_w_k: f64,
    /// Duty queries per evaluation, including zero and the nonzero endpoint.
    pub maximum_duty_evaluations: usize,
    /// Whole-call PT allowance for the one-shot function. Prepared evaluations
    /// cap their supplied remaining allowance by this same limit.
    pub maximum_pt_evaluations: usize,
}
impl MixtureWallStreamControl {
    /// Declare the signed power interval and original PH source policy.
    #[must_use]
    pub fn new(duty_bracket_w: [f64; 2], equilibrium: PhaseEnthalpyControls) -> Self {
        Self {
            duty_bracket_w,
            equilibrium,
            absolute_power_tolerance_w: 1e-5,
            relative_power_tolerance: 1e-9,
            entropy_tolerance_w_k: 1e-7,
            maximum_duty_evaluations: MAX_HEAT_EXCHANGER_DUTY_EVALUATIONS,
            maximum_pt_evaluations: MAX_HEAT_EXCHANGER_PT_EVALUATIONS,
        }
    }

    /// Admit the original source policy and bounded numeric controls.
    pub fn validate(self) -> Result<(), Error> {
        self.equilibrium.validate().map_err(|e| phase_error(0, 0, e))?;
        let [lo, hi] = self.duty_bracket_w;
        if !lo.is_finite() || !hi.is_finite() || lo > 0.0 || hi < 0.0 || lo >= hi
            || !positive(self.absolute_power_tolerance_w)
            || !self.relative_power_tolerance.is_finite()
            || !(0.0..1.0).contains(&self.relative_power_tolerance)
            || !self.entropy_tolerance_w_k.is_finite() || self.entropy_tolerance_w_k < 0.0
            || !(1..=MAX_HEAT_EXCHANGER_DUTY_EVALUATIONS).contains(&self.maximum_duty_evaluations)
            || !(1..=MAX_HEAT_EXCHANGER_PT_EVALUATIONS).contains(&self.maximum_pt_evaluations)
        {
            return Err(Error::Invalid("mixed wall-stream controls"));
        }
        Ok(())
    }

    fn tolerance(self, scale: f64) -> f64 {
        self.absolute_power_tolerance_w + self.relative_power_tolerance * scale
    }
}

/// An original phase-equilibrium stream prepared once for repeated wall trials.
///
/// The inlet selection and its work are immutable. An evaluation reports only
/// new work; a surrounding finite-wall solve charges preparation once, then
/// passes its remaining shared allowance into each evaluation.
#[derive(Clone, Debug)]
pub struct PreparedMixtureWallStream<'m> {
    inlet: MixtureHeatExchangerInlet<'m>,
    path: PreparedPhaseEnthalpyPath<'m>,
    controls: MixtureWallStreamControl,
}
impl<'m> PreparedMixtureWallStream<'m> {
    /// Select the original inlet once and reserve source accuracy in watts.
    pub fn new(
        inlet: MixtureHeatExchangerInlet<'m>,
        mut controls: MixtureWallStreamControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(&mut cancelled)?;
        controls.validate()?;
        if !positive(inlet.temperature.value()) || !positive(inlet.pressure.value())
            || !positive(inlet.mass_flow.value())
        {
            return Err(Error::Invalid("positive mixed-stream inlet and flow required"));
        }
        controls.equilibrium.enthalpy.enthalpy_tolerance_j_kg = controls.equilibrium
            .enthalpy.enthalpy_tolerance_j_kg
            .min(0.125 * controls.absolute_power_tolerance_w / inlet.mass_flow.value());
        controls.equilibrium.enthalpy.maximum_flash_evaluations = controls.equilibrium
            .enthalpy.maximum_flash_evaluations.min(controls.maximum_pt_evaluations);
        controls.validate()?;
        let path = inlet.model.prepare_phase_enthalpy_path(
            inlet.temperature, inlet.pressure, inlet.pressure,
            controls.equilibrium, &mut cancelled,
        ).map_err(|e| phase_error(0, 0, e))?;
        Ok(Self { inlet, path, controls })
    }

    /// Original PT selections consumed when preparing the immutable inlet.
    #[must_use]
    pub fn preparation_pt_evaluations(&self) -> usize {
        self.path.preparation_pt_evaluations()
    }

    /// Preparation's native EOS work ceiling, not an observed EOS count.
    #[must_use]
    pub fn preparation_eos_evaluation_upper_bound(&self) -> usize {
        self.path.preparation_eos_evaluation_upper_bound()
    }

    /// Solve a perfectly mixed steady chamber against one fixed temperature.
    ///
    /// Signed heat power satisfies `q = G * (T_wall - T_out(q / mass_flow))`.
    /// Stream enthalpy already includes pressure-flow work. There is no shaft
    /// work, pressure drop, kinetic/potential change, fluid storage, or plug-flow
    /// geometry. The supplied wall is a temperature boundary for this query;
    /// callers may couple it to a separately conserved finite wall.
    ///
    /// Returned counters cover only this evaluation. The effective work limits
    /// are the smaller of the supplied remaining allowances and the prepared
    /// controls. Zero heat reuses the inlet without new PT work. A weak thermal
    /// drive tightens contact closure; if the fixed source policy cannot resolve
    /// it, the final balance refuses the result instead of inventing heat.
    pub fn evaluate(
        &self,
        wall_temperature: Temperature,
        conductance: ThermalConductance,
        remaining_duty_evaluations: usize,
        remaining_pt_evaluations: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MixtureWallStreamStep, Error> {
        poll(&mut cancelled)?;
        let tw = wall_temperature.value();
        let g = conductance.value();
        if !positive(tw) || !g.is_finite() || g < 0.0 {
            return Err(Error::Invalid("positive wall temperature and nonnegative conductance required"));
        }
        let gap = tw - self.inlet.temperature.value();
        let initial_flux = g * gap;
        if !initial_flux.is_finite() || initial_flux == 0.0 && g != 0.0 && gap != 0.0 {
            return Err(Error::Invalid("unrepresentable wall-stream heat power"));
        }
        let mut controls = self.controls;
        if initial_flux != 0.0 {
            controls.absolute_power_tolerance_w =
                controls.absolute_power_tolerance_w.min(0.125 * initial_flux.abs());
        }
        if !positive(controls.absolute_power_tolerance_w) {
            return Err(Error::Invalid("unrepresentable wall-stream power tolerance"));
        }
        let maximum_duty = controls.maximum_duty_evaluations.min(remaining_duty_evaluations);
        let maximum_pt = controls.maximum_pt_evaluations.min(remaining_pt_evaluations);
        let mut work = Work::default();
        let mut evaluate = |q: f64, work: &mut Work| {
            poll(&mut cancelled)?;
            if work.duty == maximum_duty || work.pt == maximum_pt && q != 0.0 {
                return Err(work.exhausted());
            }
            work.duty += 1;
            let increment = q / self.inlet.mass_flow.value();
            if !increment.is_finite() || q != 0.0 && increment == 0.0 {
                return Err(Error::Invalid("unrepresentable stream enthalpy increment"));
            }
            let outlet = self.path.solve_increment(increment, maximum_pt - work.pt, &mut cancelled)
                .map_err(|e| phase_error(0, work.pt, e))?;
            work.pt += outlet.pt_evaluations();
            work.eos += outlet.eos_evaluation_upper_bound();
            let flux = g * (tw - outlet.state().temperature_k());
            if !flux.is_finite() || !(q - flux).is_finite() {
                return Err(Error::Invalid("unrepresentable mixed-stream contact power"));
            }
            Ok(Trial { duty: q, flux, outlet })
        };
        let initial = evaluate(0.0, &mut work)?;
        let accepted = if initial_flux == 0.0 {
            initial
        } else {
            let lo = controls.duty_bracket_w[0].max(initial_flux.min(0.0));
            let hi = controls.duty_bracket_w[1].min(initial_flux.max(0.0));
            if lo >= hi {
                return Err(Error::Invalid("empty passive wall-stream duty bracket"));
            }
            let endpoint = evaluate(if initial_flux > 0.0 { hi } else { lo }, &mut work)?;
            let (mut lower, mut upper) = if initial_flux > 0.0 {
                ([lo, initial.residual()], [hi, endpoint.residual()])
            } else {
                ([lo, endpoint.residual()], [hi, initial.residual()])
            };
            let mut candidate = endpoint;
            loop {
                let tolerance = 0.5 * controls.tolerance(candidate.duty.abs().max(candidate.flux.abs()));
                if !positive(tolerance) {
                    return Err(Error::Invalid("unrepresentable wall-stream closure tolerance"));
                }
                if candidate.residual().abs() <= tolerance { break candidate; }
                if lower[1] > 0.0 || upper[1] < 0.0 {
                    return Err(Error::NotBracketed {
                        duty_bracket_w: [lower[0], upper[0]],
                        residual_w: [lower[1], upper[1]],
                    });
                }
                let fraction = -lower[1] / (upper[1] - lower[1]);
                let fraction = if (0.05..0.95).contains(&fraction) { fraction } else { 0.5 };
                let q = lower[0] + fraction * (upper[0] - lower[0]);
                if q == lower[0] || q == upper[0] { return Err(Error::NoConvergence); }
                candidate = evaluate(q, &mut work)?;
                if candidate.residual() < 0.0 {
                    lower = [q, candidate.residual()];
                } else {
                    upper = [q, candidate.residual()];
                }
            }
        };
        let power = self.inlet.mass_flow.value() * accepted.outlet.accepted_increment_j_kg()
            .ok_or(Error::Invalid("missing source enthalpy increment"))?;
        let entropy = self.inlet.mass_flow.value() * accepted.outlet.entropy_change_j_kg_k()
            .ok_or(Error::Invalid("missing source entropy increment"))?;
        let target_residual = power - accepted.duty;
        let contact_residual = power - accepted.flux;
        let tolerance = controls.tolerance(power.abs().max(accepted.flux.abs()));
        if !positive(tolerance) || [target_residual, contact_residual].iter()
            .any(|r| !r.is_finite() || r.abs() > tolerance)
        {
            return Err(Error::Balance {
                energy_residual_w: target_residual,
                thermal_residual_w: [contact_residual, accepted.flux - accepted.duty],
            });
        }
        // Debit the boundary by the actual source-integrated heat, once.
        let wall_entropy = -power / tw;
        let production = entropy + wall_entropy;
        if !production.is_finite() || production < -controls.entropy_tolerance_w_k {
            return Err(Error::Entropy {
                production_w_k: production, tolerance_w_k: controls.entropy_tolerance_w_k,
            });
        }
        poll(&mut cancelled)?;
        Ok(MixtureWallStreamStep {
            outlet_equilibrium: accepted.outlet,
            target_duty: Power::new(accepted.duty),
            enthalpy_flow_change_w: power,
            enthalpy_target_residual_w: target_residual,
            heat_transfer_residual_w: contact_residual,
            entropy_flow_change_w_k: entropy,
            wall_entropy_flow_change_w_k: wall_entropy,
            entropy_production_w_k: production,
            duty_evaluations: work.duty,
            pt_evaluations: work.pt,
            eos_evaluation_upper_bound: work.eos,
        })
    }
}

struct Trial { duty: f64, flux: f64, outlet: PhaseEnthalpySolution }
impl Trial { fn residual(&self) -> f64 { self.duty - self.flux } }

/// A steady outlet, source-integrated heat and entropy, and actual work receipts.
#[derive(Clone, Debug, PartialEq)]
pub struct MixtureWallStreamStep {
    /// Original PH outlet and receipt, whose work covers only the accepted query.
    pub outlet_equilibrium: PhaseEnthalpySolution,
    /// Requested signed heat into the stream [W].
    pub target_duty: Power,
    /// Actual native `mass_flow * Delta h` [W]; use this for a finite wall debit.
    pub enthalpy_flow_change_w: f64,
    /// Actual native power minus the requested PH power [W].
    pub enthalpy_target_residual_w: f64,
    /// Actual native power minus `G * (T_wall - T_out)` [W].
    pub heat_transfer_residual_w: f64,
    /// Native `mass_flow * Delta s` [W/K].
    pub entropy_flow_change_w_k: f64,
    /// Boundary entropy flow `-actual_heat / T_wall` [W/K].
    pub wall_entropy_flow_change_w_k: f64,
    /// Unclipped stream-plus-boundary entropy production [W/K].
    pub entropy_production_w_k: f64,
    /// Started duty queries across this evaluation.
    pub duty_evaluations: usize,
    /// New original PT calls; the one-shot wrapper also includes preparation.
    pub pt_evaluations: usize,
    /// Native EOS work ceilings with the same scope as `pt_evaluations`.
    pub eos_evaluation_upper_bound: usize,
}

/// Prepare and solve one mixed flow against a fixed wall temperature.
///
/// All work counters include preparation. Reuse `PreparedMixtureWallStream`
/// when an outer finite-wall solve needs several immutable temperature trials.
pub fn exchange_mixture_stream_with_wall(
    inlet: MixtureHeatExchangerInlet<'_>,
    wall_temperature: Temperature,
    conductance: ThermalConductance,
    controls: MixtureWallStreamControl,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MixtureWallStreamStep, Error> {
    let prepared = PreparedMixtureWallStream::new(inlet, controls, &mut cancelled)?;
    let preparation_pt = prepared.preparation_pt_evaluations();
    let mut result = prepared.evaluate(
        wall_temperature, conductance, controls.maximum_duty_evaluations,
        controls.maximum_pt_evaluations - preparation_pt, &mut cancelled,
    ).map_err(|e| match e {
        Error::WorkLimit { duty_evaluations, pt_evaluations } => Error::WorkLimit {
            duty_evaluations, pt_evaluations: preparation_pt + pt_evaluations,
        },
        Error::Equilibrium { side, completed_pt_evaluations, source } => Error::Equilibrium {
            side, completed_pt_evaluations: preparation_pt + completed_pt_evaluations, source,
        },
        other => other,
    })?;
    result.pt_evaluations += preparation_pt;
    result.eos_evaluation_upper_bound += prepared.preparation_eos_evaluation_upper_bound();
    Ok(result)
}
