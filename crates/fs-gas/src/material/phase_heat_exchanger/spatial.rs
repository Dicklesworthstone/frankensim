//! One-dimensional parallel/counterflow heat exchange in enthalpy coordinates.
//!
//! At each cumulative heat fraction, both original PH paths recover the local
//! temperatures and phases. Integrating reciprocal co-located temperature gap
//! gives the required conductance; one bounded duty root closes the exchanger.

use super::{
    MAX_HEAT_EXCHANGER_PT_EVALUATIONS, MixtureHeatExchangerControl, MixtureHeatExchangerError,
    MixtureHeatExchangerInlet, MixtureHeatExchangerStep, Trial, Work, finish, phase_error, poll,
    positive, prepare,
};
use crate::chemistry::ThermalConductance;
use crate::material::compensated_sum;
use fs_material::thermochemical_gas::nonideal::flash::equilibrium::{
    PhaseEnthalpyControls, PhaseEnthalpySolution, prepared::PreparedPhaseEnthalpyPath,
};
use std::collections::VecDeque;

/// Maximum original PT queries across preparations, quadrature and duty roots.
pub const MAX_SPATIAL_PT_EVALUATIONS: usize = 1_048_576;
/// Maximum adaptive Simpson panels in one complete duty query.
pub const MAX_SPATIAL_PANELS: usize = 4096;

/// Direction of the right stream relative to the left stream's spatial axis.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FlowArrangement {
    /// Both inlets lie at the same end; their cumulative heat changes oppose.
    Parallel,
    /// The right inlet lies at the left outlet end. Outlet temperatures may
    /// cross; the temperature gap must remain passive at co-located stations.
    Counterflow,
}

/// Original source controls, duty closure and explicit spatial quadrature limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpatialHeatExchangerControl {
    /// Finite increasing powers into the left stream [W], containing zero.
    pub duty_bracket_w: [f64; 2],
    /// Original left PH bracket, tolerances, phase policy and per-query work.
    pub left_equilibrium: PhaseEnthalpyControls,
    /// Original right PH bracket, tolerances, phase policy and per-query work.
    pub right_equilibrium: PhaseEnthalpyControls,
    /// Absolute duty-root and original stream-energy tolerance [W].
    pub absolute_power_tolerance_w: f64,
    /// Relative duty-root and original stream-energy tolerance.
    pub relative_power_tolerance: f64,
    /// Numerical allowance for negative combined original stream entropy [W/K].
    pub entropy_tolerance_w_k: f64,
    /// Whole-call duty queries, including the analytic zero-duty bracket.
    pub maximum_duty_evaluations: usize,
    /// Whole-call PT work, including both immutable inlet preparations.
    pub maximum_pt_evaluations: usize,
    /// Absolute estimated conductance quadrature error allowance [W/K].
    /// This is separate from the duty-root residual and is not a certificate.
    pub absolute_conductance_tolerance_w_k: f64,
    /// Estimated conductance error allowance relative to supplied total UA.
    pub relative_conductance_tolerance: f64,
    /// Maximum adaptive panels per duty query; their native samples are cached.
    pub maximum_panels: usize,
    /// Maximum dyadic refinement depth, including pinch-bound refinement.
    pub maximum_refinement_depth: usize,
}
impl SpatialHeatExchangerControl {
    /// Declare physical duty limits and independent original source policies.
    #[must_use]
    pub fn new(
        duty_bracket_w: [f64; 2],
        left_equilibrium: PhaseEnthalpyControls,
        right_equilibrium: PhaseEnthalpyControls,
    ) -> Self {
        let base =
            MixtureHeatExchangerControl::new(duty_bracket_w, left_equilibrium, right_equilibrium);
        Self {
            duty_bracket_w,
            left_equilibrium,
            right_equilibrium,
            absolute_power_tolerance_w: base.absolute_power_tolerance_w,
            relative_power_tolerance: base.relative_power_tolerance,
            entropy_tolerance_w_k: base.entropy_tolerance_w_k,
            maximum_duty_evaluations: base.maximum_duty_evaluations,
            maximum_pt_evaluations: MAX_SPATIAL_PT_EVALUATIONS,
            absolute_conductance_tolerance_w_k: 1e-6,
            relative_conductance_tolerance: 1e-7,
            maximum_panels: 1024,
            maximum_refinement_depth: 32,
        }
    }

    /// Validate physical/source controls and bounded spatial work before compute.
    pub fn validate(self) -> Result<(), SpatialHeatExchangerError> {
        self.physical().validate()?;
        if !(1..=MAX_SPATIAL_PT_EVALUATIONS).contains(&self.maximum_pt_evaluations)
            || !positive(self.absolute_conductance_tolerance_w_k)
            || !self.relative_conductance_tolerance.is_finite()
            || !(0.0..1.0).contains(&self.relative_conductance_tolerance)
            || !(1..=MAX_SPATIAL_PANELS).contains(&self.maximum_panels)
            || !(1..=48).contains(&self.maximum_refinement_depth)
        {
            return Err(MixtureHeatExchangerError::Invalid(
                "spatial exchanger quadrature/work controls",
            )
            .into());
        }
        Ok(())
    }

    // Reuse the original physical admission and final balance owner. The larger
    // spatial PT budget is enforced by this module, not this projection.
    fn physical(self) -> MixtureHeatExchangerControl {
        MixtureHeatExchangerControl {
            duty_bracket_w: self.duty_bracket_w,
            left_equilibrium: self.left_equilibrium,
            right_equilibrium: self.right_equilibrium,
            absolute_power_tolerance_w: self.absolute_power_tolerance_w,
            relative_power_tolerance: self.relative_power_tolerance,
            entropy_tolerance_w_k: self.entropy_tolerance_w_k,
            maximum_duty_evaluations: self.maximum_duty_evaluations,
            maximum_pt_evaluations: self
                .maximum_pt_evaluations
                .min(MAX_HEAT_EXCHANGER_PT_EVALUATIONS),
        }
    }
}

/// No partial endpoint or spatial profile is returned on refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum SpatialHeatExchangerError {
    /// Original physical admission, source refusal, work limit or balance failure.
    Exchanger(MixtureHeatExchangerError),
    /// Refinement could not meet quadrature accuracy or prove the panel's
    /// positive gap under the native monotone caloric inversion assumption.
    QuadratureLimit {
        /// Adaptive panels allocated in this duty query.
        panels: usize,
        /// Depth of the unresolved panel.
        depth: usize,
        /// Monotone-temperature lower bound on its co-located gap [K].
        gap_lower_bound_k: f64,
        /// Coarse/fine conductance difference on this panel [W/K].
        estimated_error_w_k: f64,
    },
    /// Cancellation before complete publication.
    Cancelled,
}
impl From<MixtureHeatExchangerError> for SpatialHeatExchangerError {
    fn from(error: MixtureHeatExchangerError) -> Self {
        match error {
            MixtureHeatExchangerError::Cancelled => Self::Cancelled,
            other => Self::Exchanger(other),
        }
    }
}
impl core::fmt::Display for SpatialHeatExchangerError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spatial equilibrium mixture heat exchanger: {self:?}")
    }
}
impl std::error::Error for SpatialHeatExchangerError {}
type Error = SpatialHeatExchangerError;

/// Co-located original-source states on the accepted adaptive spatial profile.
#[derive(Clone, Debug, PartialEq)]
pub struct SpatialHeatExchangerStation {
    /// Fraction of total exchanged heat measured along the left-flow direction.
    /// This generally differs from axial position and conductance fraction.
    pub heat_fraction: f64,
    /// Fraction of cumulative integrated UA. It equals axial position fraction
    /// only when conductance per unit length is uniform.
    pub conductance_fraction: f64,
    /// Original left state and increment from the immutable left inlet.
    pub left_equilibrium: PhaseEnthalpySolution,
    /// Original right state and increment from the immutable right inlet.
    pub right_equilibrium: PhaseEnthalpySolution,
}

/// Accepted plug-flow duty, source balances and spatial accuracy evidence.
#[derive(Clone, Debug, PartialEq)]
pub struct SpatialHeatExchangerStep {
    /// Both actual source outlets, energy/entropy audits and all original work.
    /// The thermal residual uses the integrated spatial law, not outlet gap.
    pub balance: MixtureHeatExchangerStep,
    /// Declared relative stream directions.
    pub arrangement: FlowArrangement,
    /// Ordered retained quadrature stations, including both physical ends.
    pub profile: Vec<SpatialHeatExchangerStation>,
    /// Accepted adaptive panels; zero when no heat transfer is required.
    pub panels: usize,
    /// Original PH queries across both streams and every duty trial.
    pub caloric_queries: usize,
    /// Minimum sampled co-located gap in the passive heat-transfer direction [K].
    pub minimum_temperature_gap_k: f64,
    /// Absolute coarse/fine conductance differences summed over accepted panels.
    /// An estimate only: phase boundaries can reduce quadrature convergence order.
    pub quadrature_error_estimate_w_k: f64,
    /// Conditional power uncertainty from the integral error estimate [W],
    /// plus the scalar duty residual. It excludes native PH inversion error
    /// and model/source-data error; the original PH tolerances remain separate.
    pub estimated_duty_error_w: f64,
    /// Integrated required UA minus supplied UA [W/K]; absent at exact rest.
    pub conductance_residual_w_k: Option<f64>,
}

/// Resolve source-exact caloric profiles for steady parallel or counterflow.
///
/// For signed duty Q and heat fraction u, the left increment is u Q / mdot_L.
/// The right increment is -u Q / mdot_R for parallel flow and -(1-u) Q / mdot_R
/// for counterflow. Original PH inversions therefore conserve energy without
/// chaining approximate inlet states. The continuous balance is
/// `UA = Q integral[du / (T_right(u) - T_left(u))]`.
///
/// Nested Simpson panels reuse endpoints/midpoints; native PH work grows
/// linearly with sampled nodes per duty query. Counterflow panel admission uses
/// monotone T(h) endpoint bounds to exclude an unobserved local pinch, refining
/// even when Simpson error is zero. Witnessed nonpositive gaps or a monotone
/// lower bound on required UA may shrink the upper duty bracket without fully
/// resolving a near-pinch trial. Source refusals and unresolved panels remain
/// failures. Estimated quadrature power uncertainty shares the duty allowance.
///
/// No pressure drop, axial conduction, wall storage, external heat, shaft work,
/// species exchange or kinetic/potential energy change is modeled. The phase
/// policy remains the original finite screen, not a global stability proof.
/// Quadrature tolerances control an estimate, separately from the physical
/// duty-root and source balances. No fixed-Cp effectiveness model is used.
pub fn exchange_spatial_mixture_heat(
    left: MixtureHeatExchangerInlet<'_>,
    right: MixtureHeatExchangerInlet<'_>,
    conductance: ThermalConductance,
    arrangement: FlowArrangement,
    mut controls: SpatialHeatExchangerControl,
    mut cancelled: impl FnMut() -> bool,
) -> Result<SpatialHeatExchangerStep, Error> {
    poll(&mut cancelled)?;
    controls.validate()?;
    let inlets = [left, right];
    if inlets.iter().any(|s| {
        !positive(s.temperature.value())
            || !positive(s.pressure.value())
            || !positive(s.mass_flow.value())
    }) || !conductance.value().is_finite()
        || conductance.value() < 0.0
    {
        return Err(MixtureHeatExchangerError::Invalid(
            "positive spatial inlet states/flows and nonnegative UA required",
        )
        .into());
    }
    let gap = right.temperature.value() - left.temperature.value();
    let direction = if gap < 0.0 { -1.0 } else { 1.0 };
    let initial_flux = conductance.value() * gap.abs();
    if !initial_flux.is_finite() || initial_flux == 0.0 && conductance.value() != 0.0 && gap != 0.0
    {
        return Err(
            MixtureHeatExchangerError::Invalid("unrepresentable spatial initial power").into(),
        );
    }
    if initial_flux != 0.0 {
        controls.absolute_power_tolerance_w = controls
            .absolute_power_tolerance_w
            .min(0.125 * initial_flux);
        for (source, stream) in [
            (&mut controls.left_equilibrium, left),
            (&mut controls.right_equilibrium, right),
        ] {
            source.enthalpy.enthalpy_tolerance_j_kg = source
                .enthalpy
                .enthalpy_tolerance_j_kg
                .min(0.125 * controls.absolute_power_tolerance_w / stream.mass_flow.value());
        }
        controls.validate()?;
    }
    let mut work = Work::default();
    let physical = controls.physical();
    let paths = [
        prepare(
            left,
            0,
            controls.left_equilibrium,
            physical,
            &mut work,
            &mut cancelled,
        )?,
        prepare(
            right,
            1,
            controls.right_equilibrium,
            physical,
            &mut work,
            &mut cancelled,
        )?,
    ];
    let context = Context {
        inlets,
        paths: &paths,
        controls,
        ua: conductance.value(),
        direction,
        arrangement,
    };
    let mut queries = 0;
    work.duty = 1; // Analytic zero-duty bracket, after both original preparations.
    if initial_flux == 0.0 {
        let node = context.sample(0.0, 0.0, &mut work, &mut queries, &mut cancelled)?;
        let [l, r] = node.outlets;
        let profile = vec![
            SpatialHeatExchangerStation {
                heat_fraction: 0.0,
                conductance_fraction: 0.0,
                left_equilibrium: l.clone(),
                right_equilibrium: r.clone(),
            },
            SpatialHeatExchangerStation {
                heat_fraction: 1.0,
                conductance_fraction: 1.0,
                left_equilibrium: l.clone(),
                right_equilibrium: r.clone(),
            },
        ];
        let balance = finish(
            Trial {
                duty: 0.0,
                flux: 0.0,
                outlets: [l, r],
            },
            inlets,
            physical,
            work,
            &mut cancelled,
        )?;
        return Ok(SpatialHeatExchangerStep {
            balance,
            arrangement,
            profile,
            panels: 0,
            caloric_queries: queries,
            minimum_temperature_gap_k: gap.abs(),
            quadrature_error_estimate_w_k: 0.0,
            estimated_duty_error_w: 0.0,
            conductance_residual_w_k: None,
        });
    }
    let bound = if direction > 0.0 {
        controls.duty_bracket_w[1]
    } else {
        -controls.duty_bracket_w[0]
    };
    let mut lower = [0.0, -initial_flux];
    let mut upper = [bound.min(initial_flux), f64::NAN];
    if !positive(upper[0]) {
        return Err(
            MixtureHeatExchangerError::Invalid("empty passive spatial duty bracket").into(),
        );
    }
    let mut amplitude = upper[0];
    let accepted = loop {
        poll(&mut cancelled)?;
        if work.duty == controls.maximum_duty_evaluations {
            return Err(work.exhausted().into());
        }
        work.duty += 1;
        match context.integrate(
            direction * amplitude,
            &mut work,
            &mut queries,
            &mut cancelled,
        )? {
            Outcome::Pinched => upper = [amplitude, amplitude],
            Outcome::AboveConductance(residual_bound) => upper = [amplitude, residual_bound],
            Outcome::Resolved(trial) => {
                let flux = context.ua / trial.integral;
                let residual = amplitude - flux;
                let uncertainty = flux * (trial.error / (trial.integral - trial.error));
                let tolerance = 0.5
                    * (controls.absolute_power_tolerance_w
                        + controls.relative_power_tolerance * amplitude.abs().max(flux.abs()));
                if !positive(tolerance) || !residual.is_finite() || !uncertainty.is_finite() {
                    return Err(MixtureHeatExchangerError::Invalid(
                        "unrepresentable spatial root residual",
                    )
                    .into());
                }
                if residual.abs() + uncertainty <= tolerance {
                    break trial;
                }
                if amplitude == upper[0] && residual < 0.0 {
                    let (duty_bracket_w, residual_w) = if direction > 0.0 {
                        ([0.0, amplitude], [lower[1], residual])
                    } else {
                        ([-amplitude, 0.0], [-residual, -lower[1]])
                    };
                    return Err(MixtureHeatExchangerError::NotBracketed {
                        duty_bracket_w,
                        residual_w,
                    }
                    .into());
                }
                if residual < 0.0 {
                    lower = [amplitude, residual];
                } else {
                    upper = [amplitude, residual];
                }
            }
        }
        let fraction = -lower[1] / (upper[1] - lower[1]);
        let fraction = if (0.05..0.95).contains(&fraction) {
            fraction
        } else {
            0.5
        };
        amplitude = lower[0] + fraction * (upper[0] - lower[0]);
        if amplitude == lower[0] || amplitude == upper[0] {
            return Err(MixtureHeatExchangerError::NoConvergence.into());
        }
    };
    publish(accepted, &context, physical, work, queries, &mut cancelled)
}

struct Context<'a, 'm> {
    inlets: [MixtureHeatExchangerInlet<'m>; 2],
    paths: &'a [PreparedPhaseEnthalpyPath<'m>; 2],
    controls: SpatialHeatExchangerControl,
    ua: f64,
    direction: f64,
    arrangement: FlowArrangement,
}
struct Node {
    u: f64,
    gap: f64,
    outlets: [PhaseEnthalpySolution; 2],
}
#[derive(Clone, Copy)]
struct Panel {
    a: usize,
    m: usize,
    b: usize,
    depth: usize,
}
struct Leaf {
    middle: usize,
    end: usize,
    integral: [f64; 2],
    error: f64,
}
struct IntegralTrial {
    duty: f64,
    nodes: Vec<Node>,
    leaves: Vec<Leaf>,
    integral: f64,
    error: f64,
}
enum Outcome {
    Resolved(IntegralTrial),
    Pinched,
    AboveConductance(f64),
}

impl Context<'_, '_> {
    // A partition-wide lower bound on integral(1/gap). Each interval's upper
    // gap follows from monotone source T(h); it remains a valid exclusion if an
    // unsampled internal pinch makes this trial physically inadmissible.
    fn integral_lower_bound(&self, nodes: &[Node], a: usize, b: usize) -> Result<f64, Error> {
        let upper_gap = if self.arrangement == FlowArrangement::Parallel {
            nodes[a].gap
        } else {
            self.direction
                * (nodes[b].outlets[1].state().temperature_k()
                    - nodes[a].outlets[0].state().temperature_k())
        };
        let bound = (nodes[b].u - nodes[a].u) / upper_gap;
        if !positive(upper_gap) || !positive(bound) {
            return Err(
                MixtureHeatExchangerError::Invalid("unrepresentable monotone gap bound").into(),
            );
        }
        Ok(bound)
    }

    fn above_conductance(&self, duty: f64, integral_lower_bound: f64) -> Option<Outcome> {
        let residual_bound = duty.abs() - self.ua / integral_lower_bound;
        let tolerance = self.controls.absolute_power_tolerance_w
            + self.controls.relative_power_tolerance * duty.abs();
        (residual_bound > tolerance).then_some(Outcome::AboveConductance(residual_bound))
    }

    fn sample(
        &self,
        u: f64,
        duty: f64,
        work: &mut Work,
        queries: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Node, Error> {
        poll(cancelled)?;
        let right_fraction = if self.arrangement == FlowArrangement::Parallel {
            u
        } else {
            1.0 - u
        };
        let powers = [u * duty, -right_fraction * duty];
        let mut resolve = |side: usize| {
            let increment = powers[side] / self.inlets[side].mass_flow.value();
            if !increment.is_finite() || powers[side] != 0.0 && increment == 0.0 {
                return Err(MixtureHeatExchangerError::Invalid(
                    "unrepresentable spatial enthalpy increment",
                ));
            }
            let remaining = self.controls.maximum_pt_evaluations - work.pt;
            *queries += 1;
            let result = self.paths[side]
                .solve_increment(increment, remaining, &mut *cancelled)
                .map_err(|e| phase_error(side, work.pt, e))?;
            work.pt += result.pt_evaluations();
            work.eos += result.eos_evaluation_upper_bound();
            Ok(result)
        };
        let outlets = [resolve(0)?, resolve(1)?];
        let gap = self.direction
            * (outlets[1].state().temperature_k() - outlets[0].state().temperature_k());
        if !gap.is_finite() {
            return Err(MixtureHeatExchangerError::Invalid("nonfinite spatial gap").into());
        }
        Ok(Node { u, gap, outlets })
    }

    fn integrate(
        &self,
        duty: f64,
        work: &mut Work,
        queries: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Outcome, Error> {
        let mut nodes = Vec::new();
        for u in [0.0, 0.5, 1.0] {
            let node = self.sample(u, duty, work, queries, cancelled)?;
            if node.gap <= 0.0 {
                return Ok(Outcome::Pinched);
            }
            nodes.push(node);
        }
        let requested_allowance = self.controls.absolute_conductance_tolerance_w_k
            + self.controls.relative_conductance_tolerance * self.ua;
        let power_tolerance = self.controls.absolute_power_tolerance_w
            + self.controls.relative_power_tolerance * duty.abs();
        // At a duty root I=UA/|Q|. Reserve at most an eighth of the physical
        // power allowance for integral uncertainty, regardless of a loose
        // caller quadrature request. The actual conditional power estimate is
        // checked again by the root owner before acceptance.
        let power_allowance = (power_tolerance / duty.abs()) * (self.ua * 0.125);
        let allowance = requested_allowance.min(power_allowance);
        if !positive(allowance) {
            return Err(
                MixtureHeatExchangerError::Invalid("unrepresentable quadrature allowance").into(),
            );
        }
        let mut pending = VecDeque::from([Panel {
            a: 0,
            m: 1,
            b: 2,
            depth: 0,
        }]);
        let mut leaves = Vec::new();
        let mut panels = 1;
        let mut integral_lower_bound = self.integral_lower_bound(&nodes, 0, 2)?;
        if let Some(above) = self.above_conductance(duty, integral_lower_bound) {
            return Ok(above);
        }
        while let Some(panel) = pending.pop_front() {
            poll(cancelled)?;
            let a = &nodes[panel.a];
            let m = &nodes[panel.m];
            let b = &nodes[panel.b];
            let width = b.u - a.u;
            let coarse = width * (1.0 / a.gap + 4.0 / m.gap + 1.0 / b.gap) / 6.0;
            let bound = if self.arrangement == FlowArrangement::Parallel {
                b.gap
            } else {
                self.direction
                    * (a.outlets[1].state().temperature_k() - b.outlets[0].state().temperature_k())
            };
            let positions = [a.u + (m.u - a.u) * 0.5, m.u + (b.u - m.u) * 0.5];
            if !(a.u < positions[0]
                && positions[0] < m.u
                && m.u < positions[1]
                && positions[1] < b.u)
            {
                return Err(MixtureHeatExchangerError::NoConvergence.into());
            }
            let quarter = nodes.len();
            for u in positions {
                let node = self.sample(u, duty, work, queries, cancelled)?;
                if node.gap <= 0.0 {
                    return Ok(Outcome::Pinched);
                }
                nodes.push(node);
            }
            let left = width
                * (1.0 / nodes[panel.a].gap + 4.0 / nodes[quarter].gap + 1.0 / nodes[panel.m].gap)
                / 12.0;
            let right = width
                * (1.0 / nodes[panel.m].gap
                    + 4.0 / nodes[quarter + 1].gap
                    + 1.0 / nodes[panel.b].gap)
                / 12.0;
            let error = ((left + right) - coarse).abs();
            if !positive(left) || !positive(right) || !error.is_finite() {
                return Err(MixtureHeatExchangerError::Invalid(
                    "unrepresentable spatial quadrature",
                )
                .into());
            }
            // No /15 Richardson assumption: phase-boundary kinks may reduce
            // Simpson's smooth-function order. This remains an error estimate.
            if bound > 0.0 && duty.abs() * error <= allowance * width {
                leaves.push(Leaf {
                    middle: panel.m,
                    end: panel.b,
                    integral: [left, right],
                    error,
                });
            } else {
                if panels == self.controls.maximum_panels
                    || panel.depth == self.controls.maximum_refinement_depth
                {
                    return Err(Error::QuadratureLimit {
                        panels,
                        depth: panel.depth,
                        gap_lower_bound_k: bound,
                        estimated_error_w_k: duty.abs() * error,
                    });
                }
                panels += 1;
                integral_lower_bound = compensated_sum(
                    [
                        integral_lower_bound,
                        -self.integral_lower_bound(&nodes, panel.a, panel.b)?,
                        self.integral_lower_bound(&nodes, panel.a, panel.m)?,
                        self.integral_lower_bound(&nodes, panel.m, panel.b)?,
                    ]
                    .into_iter(),
                );
                if let Some(above) = self.above_conductance(duty, integral_lower_bound) {
                    return Ok(above);
                }
                // Breadth-first refinement strengthens the complete partition
                // bound before an arbitrarily thin near-pinch interval can
                // consume the depth limit. Samples remain cached by node index.
                pending.push_back(Panel {
                    a: panel.a,
                    m: quarter,
                    b: panel.m,
                    depth: panel.depth + 1,
                });
                pending.push_back(Panel {
                    a: panel.m,
                    m: quarter + 1,
                    b: panel.b,
                    depth: panel.depth + 1,
                });
            }
        }
        leaves.sort_by(|a, b| nodes[a.end].u.total_cmp(&nodes[b.end].u));
        let integral = compensated_sum(leaves.iter().flat_map(|p| p.integral));
        let error = compensated_sum(leaves.iter().map(|p| p.error));
        if !positive(integral) || !error.is_finite() || error >= integral {
            return Err(
                MixtureHeatExchangerError::Invalid("unresolved reciprocal-gap integral").into(),
            );
        }
        Ok(Outcome::Resolved(IntegralTrial {
            duty,
            nodes,
            leaves,
            integral,
            error,
        }))
    }
}

fn publish(
    trial: IntegralTrial,
    context: &Context<'_, '_>,
    physical: MixtureHeatExchangerControl,
    work: Work,
    queries: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<SpatialHeatExchangerStep, Error> {
    let flux = context.direction * context.ua / trial.integral;
    let uncertainty =
        (context.ua / trial.integral) * (trial.error / (trial.integral - trial.error));
    let estimate = (trial.duty - flux).abs() + uncertainty;
    let conductance_residual = trial.duty.abs() * trial.integral - context.ua;
    if !estimate.is_finite() || !conductance_residual.is_finite() {
        return Err(
            MixtureHeatExchangerError::Invalid("unrepresentable spatial error estimate").into(),
        );
    }
    let right_out = if context.arrangement == FlowArrangement::Counterflow {
        0
    } else {
        2
    };
    let outlets = [
        trial.nodes[2].outlets[0].clone(),
        trial.nodes[right_out].outlets[1].clone(),
    ];
    let minimum_gap = trial
        .nodes
        .iter()
        .map(|n| n.gap)
        .fold(f64::INFINITY, f64::min);
    let panels = trial.leaves.len();
    let mut selected = vec![(0, 0.0)];
    // Match the compensated integral reduction in the same sorted leaf order,
    // so the last normalized coordinate is one without clamping a mismatch.
    let (mut accumulated, mut correction) = (0.0, 0.0);
    let mut accumulate = |value: f64| {
        let adjusted = value - correction;
        let next = accumulated + adjusted;
        correction = (next - accumulated) - adjusted;
        accumulated = next;
        accumulated / trial.integral
    };
    for leaf in &trial.leaves {
        selected.push((leaf.middle, accumulate(leaf.integral[0])));
        selected.push((leaf.end, accumulate(leaf.integral[1])));
    }
    let mut positions = vec![usize::MAX; trial.nodes.len()];
    for (position, &(index, _)) in selected.iter().enumerate() {
        positions[index] = position;
    }
    let mut stations: Vec<_> = (0..selected.len()).map(|_| None).collect();
    for (index, node) in trial.nodes.into_iter().enumerate() {
        let position = positions[index];
        if position != usize::MAX {
            let [left_equilibrium, right_equilibrium] = node.outlets;
            stations[position] = Some(SpatialHeatExchangerStation {
                heat_fraction: node.u,
                conductance_fraction: selected[position].1,
                left_equilibrium,
                right_equilibrium,
            });
        }
    }
    let profile = stations.into_iter().collect::<Option<Vec<_>>>().ok_or(
        MixtureHeatExchangerError::Invalid("incomplete accepted spatial profile"),
    )?;
    let balance = finish(
        Trial {
            duty: trial.duty,
            flux,
            outlets,
        },
        context.inlets,
        physical,
        work,
        cancelled,
    )?;
    Ok(SpatialHeatExchangerStep {
        balance,
        arrangement: context.arrangement,
        profile,
        panels,
        caloric_queries: queries,
        minimum_temperature_gap_k: minimum_gap,
        quadrature_error_estimate_w_k: trial.duty.abs() * trial.error,
        estimated_duty_error_w: estimate,
        conductance_residual_w_k: Some(conductance_residual),
    })
}
