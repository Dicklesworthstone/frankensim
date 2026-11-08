//! Steady heat exchange between two perfectly mixed, original-law streams.
//!
//! Each chamber has its own fixed pressure and conserved species feed. The
//! original PH increment owner resolves both outlets, including phase changes;
//! neither averaged heat capacity nor an additional latent-energy term is used.

/// Original-source parallel and counterflow plug-flow heat exchange.
pub mod spatial;
/// Original-phase steady mixed flow against a prescribed wall temperature.
pub mod wall_stream;
/// Finite original material storage coupled to quasi-steady mixture streams.
#[cfg(feature = "experimental-euler")]
pub mod wall;
#[cfg(feature = "experimental-euler")]
pub mod spatial_wall;

use crate::chemistry::ThermalConductance;
use fs_material::thermochemical_gas::nonideal::ThermochemicalPengRobinsonModel;
use fs_material::thermochemical_gas::nonideal::flash::{
    ThermochemicalFlashError,
    enthalpy::MAX_ENTHALPY_FLASH_EVALUATIONS,
    equilibrium::{
        PhaseEnthalpyControls, PhaseEnthalpySolution, prepared::PreparedPhaseEnthalpyPath,
    },
};
use fs_qty::{MassFlowRate, Power, Pressure, Temperature};

/// Bounded paired duty queries, including zero and the nonzero bracket endpoint.
pub const MAX_HEAT_EXCHANGER_DUTY_EVALUATIONS: usize = 132;
/// Maximum original PT work across both streams and every duty query.
pub const MAX_HEAT_EXCHANGER_PT_EVALUATIONS: usize =
    2 * MAX_ENTHALPY_FLASH_EVALUATIONS * MAX_HEAT_EXCHANGER_DUTY_EVALUATIONS;

/// An immutable feed into a perfectly mixed, constant-pressure chamber.
#[derive(Clone, Copy, Debug)]
pub struct MixtureHeatExchangerInlet<'m> {
    /// Original species calorics, composition, EOS and reference convention.
    pub model: &'m ThermochemicalPengRobinsonModel,
    /// Positive inlet temperature inside this stream's source bracket.
    pub temperature: Temperature,
    /// Positive pressure, held at both inlet and outlet of this chamber.
    pub pressure: Pressure,
    /// Positive steady mass flow. The heat exchanger does not predict flow.
    pub mass_flow: MassFlowRate,
}

/// Independent source policies and shared physical/work allowances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixtureHeatExchangerControl {
    /// Increasing finite signed powers into the left stream [W], containing
    /// zero. The passive initial-flux bound further limits this interval; both
    /// remaining endpoints must be supported by both original source domains.
    pub duty_bracket_w: [f64; 2],
    /// Left original PH domain, phase policy and per-query PT limit.
    pub left_equilibrium: PhaseEnthalpyControls,
    /// Right original PH domain, phase policy and per-query PT limit.
    pub right_equilibrium: PhaseEnthalpyControls,
    /// Absolute allowance for energy and both heat-transfer equations [W].
    pub absolute_power_tolerance_w: f64,
    /// Relative allowance, scaled by actual duty and heat-transfer power.
    pub relative_power_tolerance: f64,
    /// Numerical allowance for negative combined entropy production [W/K].
    pub entropy_tolerance_w_k: f64,
    /// Whole-call paired duty queries, including bracket endpoints.
    pub maximum_duty_evaluations: usize,
    /// Whole-call original PT queries, including each inlet prepared once.
    pub maximum_pt_evaluations: usize,
}
impl MixtureHeatExchangerControl {
    /// Declare a duty interval and independent original caloric source controls.
    #[must_use]
    pub fn new(
        duty_bracket_w: [f64; 2],
        left_equilibrium: PhaseEnthalpyControls,
        right_equilibrium: PhaseEnthalpyControls,
    ) -> Self {
        Self {
            duty_bracket_w,
            left_equilibrium,
            right_equilibrium,
            absolute_power_tolerance_w: 1e-5,
            relative_power_tolerance: 1e-9,
            entropy_tolerance_w_k: 1e-7,
            maximum_duty_evaluations: MAX_HEAT_EXCHANGER_DUTY_EVALUATIONS,
            maximum_pt_evaluations: MAX_HEAT_EXCHANGER_PT_EVALUATIONS,
        }
    }

    /// Admit both source policies and finite physical/work limits before work.
    pub fn validate(self) -> Result<(), MixtureHeatExchangerError> {
        self.left_equilibrium
            .validate()
            .map_err(|e| phase_error(0, 0, e))?;
        self.right_equilibrium
            .validate()
            .map_err(|e| phase_error(1, 0, e))?;
        let [lo, hi] = self.duty_bracket_w;
        if !lo.is_finite()
            || !hi.is_finite()
            || lo > 0.0
            || hi < 0.0
            || lo >= hi
            || !positive(self.absolute_power_tolerance_w)
            || !self.relative_power_tolerance.is_finite()
            || !(0.0..1.0).contains(&self.relative_power_tolerance)
            || !self.entropy_tolerance_w_k.is_finite()
            || self.entropy_tolerance_w_k < 0.0
            || !(1..=MAX_HEAT_EXCHANGER_DUTY_EVALUATIONS).contains(&self.maximum_duty_evaluations)
            || !(1..=MAX_HEAT_EXCHANGER_PT_EVALUATIONS).contains(&self.maximum_pt_evaluations)
        {
            return Err(MixtureHeatExchangerError::Invalid(
                "heat-exchanger controls",
            ));
        }
        Ok(())
    }

    fn tolerance(self, scale: f64) -> f64 {
        self.absolute_power_tolerance_w + self.relative_power_tolerance * scale
    }
}

/// No partial pair of outlet states is published on refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum MixtureHeatExchangerError {
    /// Invalid input or unrepresentable physical coordinate.
    Invalid(&'static str),
    /// Original source refusal; side 0 is left, side 1 is right.
    Equilibrium {
        /// Stream whose original PH owner refused the query.
        side: usize,
        /// Completed PT calls before this source query began.
        completed_pt_evaluations: usize,
        /// Original refusal, including attempted source work when available.
        source: ThermochemicalFlashError,
    },
    /// The shared duty or PT allowance was exhausted before the next query.
    WorkLimit {
        /// Paired duty queries already started.
        duty_evaluations: usize,
        /// Original PT calls already completed.
        pt_evaluations: usize,
    },
    /// The admitted duty interval does not contain a heat-transfer root.
    NotBracketed {
        /// Actual interval after applying the passive initial-flux bound [W].
        duty_bracket_w: [f64; 2],
        /// Signed heat-transfer residual at each endpoint [W].
        residual_w: [f64; 2],
    },
    /// Floating-point duty resolution was exhausted before physical closure.
    NoConvergence,
    /// Actual source enthalpy flows failed energy or heat-transfer closure.
    Balance {
        /// Sum of the original stream enthalpy-flow changes [W].
        energy_residual_w: f64,
        /// Left/right source power minus their respective interface power [W].
        thermal_residual_w: [f64; 2],
    },
    /// Combined original stream entropy production violated its allowance.
    Entropy {
        /// Unclipped source entropy production [W/K].
        production_w_k: f64,
        /// Declared numerical allowance [W/K].
        tolerance_w_k: f64,
    },
    /// Cancellation before the complete pair was published.
    Cancelled,
}
impl core::fmt::Display for MixtureHeatExchangerError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "equilibrium mixture heat exchanger: {self:?}")
    }
}
impl std::error::Error for MixtureHeatExchangerError {}
type Error = MixtureHeatExchangerError;

/// Accepted outlet pair, physical balances and whole-call source work.
#[derive(Clone, Debug, PartialEq)]
pub struct MixtureHeatExchangerStep {
    /// Original left endpoint and PH receipt; its counts cover only this trial.
    pub left_equilibrium: PhaseEnthalpySolution,
    /// Original right endpoint and PH receipt; its counts cover only this trial.
    pub right_equilibrium: PhaseEnthalpySolution,
    /// Common signed interface duty, positive into the left stream. The right
    /// stream receives its exact negative as its caloric target.
    pub duty_to_left: Power,
    /// Actual original `mass_flow * Delta h` for `[left, right]` [W].
    pub enthalpy_flow_change_w: [f64; 2],
    /// Sum of both actual enthalpy-flow changes [W].
    pub energy_residual_w: f64,
    /// Both actual enthalpy-flow changes minus their interface powers [W].
    pub thermal_residual_w: [f64; 2],
    /// Original `mass_flow * Delta s` for `[left, right]` [W/K].
    pub entropy_flow_change_w_k: [f64; 2],
    /// Sum of original stream entropy changes [W/K], never clipped.
    pub entropy_production_w_k: f64,
    /// Paired duty queries over the entire solve.
    pub duty_evaluations: usize,
    /// Original PT calls across both streams and all duty queries.
    pub pt_evaluations: usize,
    /// Sum of original EOS work ceilings, not an observed total EOS count.
    pub eos_evaluation_upper_bound: usize,
}

/// Resolve two steady, perfectly mixed chambers coupled by a supplied `UA`.
///
/// At each chamber's fixed pressure, `Delta h_left = Q / mass_flow_left` and
/// `Delta h_right = -Q / mass_flow_right`. A bounded safeguarded secant closes
/// `Q = UA * (T_right_out - T_left_out)` using the original PH increment owners.
/// Each species feed remains separate; energy and entropy reference datums
/// cancel before rounding. Stream enthalpy already includes pressure-flow work.
///
/// This model assumes no shaft work, external heat, pressure drop, wall storage,
/// or kinetic/potential energy change. It describes perfectly mixed chambers,
/// not counterflow/parallel-flow plug-flow geometry, and adds no phase kinetics
/// or global phase-stability proof. Conductance and mass flows are supplied.
/// Zero conductance or equal inlet temperatures returns exact zero duty after
/// one original PT selection per stream. Every duty trial reuses those prepared
/// inlets; no accepted endpoint is recomputed.
pub fn exchange_steady_mixture_heat(
    left: MixtureHeatExchangerInlet<'_>,
    right: MixtureHeatExchangerInlet<'_>,
    conductance: ThermalConductance,
    mut controls: MixtureHeatExchangerControl,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MixtureHeatExchangerStep, Error> {
    poll(&mut cancelled)?;
    controls.validate()?;
    if [left, right].iter().any(|s| {
        !positive(s.temperature.value())
            || !positive(s.pressure.value())
            || !positive(s.mass_flow.value())
    }) || !conductance.value().is_finite()
        || conductance.value() < 0.0
    {
        return Err(Error::Invalid(
            "positive inlet states/flows and nonnegative conductance required",
        ));
    }
    let initial_gap = right.temperature.value() - left.temperature.value();
    let initial_flux = conductance.value() * initial_gap;
    if !initial_flux.is_finite()
        || initial_flux == 0.0 && conductance.value() != 0.0 && initial_gap != 0.0
    {
        return Err(Error::Invalid(
            "unrepresentable initial heat-transfer power",
        ));
    }
    if initial_flux != 0.0 {
        controls.absolute_power_tolerance_w = controls
            .absolute_power_tolerance_w
            .min(0.125 * initial_flux.abs());
        // Independent source errors get only an eighth of the absolute power
        // allowance each, including for unequal stream mass flows.
        for (source, flow) in [
            (&mut controls.left_equilibrium, left.mass_flow),
            (&mut controls.right_equilibrium, right.mass_flow),
        ] {
            source.enthalpy.enthalpy_tolerance_j_kg = source
                .enthalpy
                .enthalpy_tolerance_j_kg
                .min(0.125 * controls.absolute_power_tolerance_w / flow.value());
        }
        controls.validate()?;
    }
    let mut work = Work::default();
    let paths = [
        prepare(
            left,
            0,
            controls.left_equilibrium,
            controls,
            &mut work,
            &mut cancelled,
        )?,
        prepare(
            right,
            1,
            controls.right_equilibrium,
            controls,
            &mut work,
            &mut cancelled,
        )?,
    ];
    let evaluate = |q, work: &mut Work, cancel: &mut _| {
        trial(
            [left, right],
            &paths,
            conductance.value(),
            q,
            controls,
            work,
            cancel,
        )
    };
    let initial = evaluate(0.0, &mut work, &mut cancelled)?;
    let accepted = if initial_flux == 0.0 {
        initial
    } else {
        let lo = controls.duty_bracket_w[0].max(initial_flux.min(0.0));
        let hi = controls.duty_bracket_w[1].min(initial_flux.max(0.0));
        if lo >= hi {
            return Err(Error::Invalid("empty passive heat-exchanger duty bracket"));
        }
        let endpoint = evaluate(
            if initial_flux > 0.0 { hi } else { lo },
            &mut work,
            &mut cancelled,
        )?;
        let (mut lower, mut upper) = if initial_flux > 0.0 {
            ([lo, initial.residual()], [hi, endpoint.residual()])
        } else {
            ([lo, endpoint.residual()], [hi, initial.residual()])
        };
        let mut candidate = endpoint;
        loop {
            let tolerance =
                0.5 * controls.tolerance(candidate.duty.abs().max(candidate.flux.abs()));
            if !positive(tolerance) {
                return Err(Error::Invalid("unrepresentable power tolerance"));
            }
            if candidate.residual().abs() <= tolerance {
                break candidate;
            }
            if lower[1] > 0.0 || upper[1] < 0.0 {
                return Err(Error::NotBracketed {
                    duty_bracket_w: [lower[0], upper[0]],
                    residual_w: [lower[1], upper[1]],
                });
            }
            let fraction = -lower[1] / (upper[1] - lower[1]);
            let fraction = if (0.05..0.95).contains(&fraction) {
                fraction
            } else {
                0.5
            };
            let q = lower[0] + fraction * (upper[0] - lower[0]);
            if q == lower[0] || q == upper[0] {
                return Err(Error::NoConvergence);
            }
            candidate = evaluate(q, &mut work, &mut cancelled)?;
            if candidate.residual() < 0.0 {
                lower = [q, candidate.residual()];
            } else {
                upper = [q, candidate.residual()];
            }
        }
    };
    finish(accepted, [left, right], controls, work, &mut cancelled)
}

#[derive(Default)]
struct Work {
    duty: usize,
    pt: usize,
    eos: usize,
}
impl Work {
    fn exhausted(&self) -> Error {
        Error::WorkLimit {
            duty_evaluations: self.duty,
            pt_evaluations: self.pt,
        }
    }
}
struct Trial {
    duty: f64,
    flux: f64,
    outlets: [PhaseEnthalpySolution; 2],
}
impl Trial {
    fn residual(&self) -> f64 {
        self.duty - self.flux
    }
}
fn prepare<'m>(
    inlet: MixtureHeatExchangerInlet<'m>,
    side: usize,
    mut source: PhaseEnthalpyControls,
    controls: MixtureHeatExchangerControl,
    work: &mut Work,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PreparedPhaseEnthalpyPath<'m>, Error> {
    let remaining = controls.maximum_pt_evaluations - work.pt;
    if remaining == 0 {
        return Err(work.exhausted());
    }
    source.enthalpy.maximum_flash_evaluations =
        source.enthalpy.maximum_flash_evaluations.min(remaining);
    let path = inlet
        .model
        .prepare_phase_enthalpy_path(
            inlet.temperature,
            inlet.pressure,
            inlet.pressure,
            source,
            &mut *cancelled,
        )
        .map_err(|e| phase_error(side, work.pt, e))?;
    work.pt += path.preparation_pt_evaluations();
    work.eos += path.preparation_eos_evaluation_upper_bound();
    Ok(path)
}

fn trial(
    inlets: [MixtureHeatExchangerInlet<'_>; 2],
    paths: &[PreparedPhaseEnthalpyPath<'_>; 2],
    ua: f64,
    duty: f64,
    controls: MixtureHeatExchangerControl,
    work: &mut Work,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Trial, Error> {
    poll(cancelled)?;
    if work.duty == controls.maximum_duty_evaluations {
        return Err(work.exhausted());
    }
    work.duty += 1;
    let mut resolve = |side: usize, power: f64| {
        let remaining = controls.maximum_pt_evaluations - work.pt;
        if remaining == 0 && power != 0.0 {
            return Err(work.exhausted());
        }
        let stream = inlets[side];
        let increment = power / stream.mass_flow.value();
        if !increment.is_finite() || power != 0.0 && increment == 0.0 {
            return Err(Error::Invalid("unrepresentable stream enthalpy increment"));
        }
        let outlet = paths[side]
            .solve_increment(increment, remaining, &mut *cancelled)
            .map_err(|e| phase_error(side, work.pt, e))?;
        work.pt += outlet.pt_evaluations();
        work.eos += outlet.eos_evaluation_upper_bound();
        Ok(outlet)
    };
    let outlets = [
        resolve(0, duty)?,
        resolve(1, if duty == 0.0 { 0.0 } else { -duty })?,
    ];
    let flux = ua * (outlets[1].state().temperature_k() - outlets[0].state().temperature_k());
    if !flux.is_finite() || !(duty - flux).is_finite() {
        return Err(Error::Invalid("unrepresentable outlet heat-transfer power"));
    }
    Ok(Trial {
        duty,
        flux,
        outlets,
    })
}

fn finish(
    trial: Trial,
    inlets: [MixtureHeatExchangerInlet<'_>; 2],
    controls: MixtureHeatExchangerControl,
    work: Work,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<MixtureHeatExchangerStep, Error> {
    let mut enthalpy = [0.0; 2];
    let mut entropy = [0.0; 2];
    for side in 0..2 {
        enthalpy[side] = inlets[side].mass_flow.value()
            * trial.outlets[side]
                .accepted_increment_j_kg()
                .ok_or(Error::Invalid("missing source enthalpy increment"))?;
        entropy[side] = inlets[side].mass_flow.value()
            * trial.outlets[side]
                .entropy_change_j_kg_k()
                .ok_or(Error::Invalid("missing source entropy increment"))?;
    }
    let energy = enthalpy[0] + enthalpy[1];
    let thermal = [enthalpy[0] - trial.flux, enthalpy[1] + trial.flux];
    let tolerance = controls.tolerance(
        enthalpy[0]
            .abs()
            .max(enthalpy[1].abs())
            .max(trial.flux.abs()),
    );
    if !positive(tolerance)
        || [energy, thermal[0], thermal[1]]
            .iter()
            .any(|r| !r.is_finite() || r.abs() > tolerance)
    {
        return Err(Error::Balance {
            energy_residual_w: energy,
            thermal_residual_w: thermal,
        });
    }
    let production = entropy[0] + entropy[1];
    if !production.is_finite() || production < -controls.entropy_tolerance_w_k {
        return Err(Error::Entropy {
            production_w_k: production,
            tolerance_w_k: controls.entropy_tolerance_w_k,
        });
    }
    poll(cancelled)?;
    let [left_equilibrium, right_equilibrium] = trial.outlets;
    Ok(MixtureHeatExchangerStep {
        left_equilibrium,
        right_equilibrium,
        duty_to_left: Power::new(trial.duty),
        enthalpy_flow_change_w: enthalpy,
        energy_residual_w: energy,
        thermal_residual_w: thermal,
        entropy_flow_change_w_k: entropy,
        entropy_production_w_k: production,
        duty_evaluations: work.duty,
        pt_evaluations: work.pt,
        eos_evaluation_upper_bound: work.eos,
    })
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn phase_error(
    side: usize,
    completed_pt_evaluations: usize,
    source: ThermochemicalFlashError,
) -> Error {
    if source == ThermochemicalFlashError::Cancelled {
        Error::Cancelled
    } else {
        Error::Equilibrium {
            side,
            completed_pt_evaluations,
            source,
        }
    }
}
