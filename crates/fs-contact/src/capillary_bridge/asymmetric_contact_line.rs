//! Dynamic planar liquid adhesion between two independently sourced surfaces.
//!
//! Both wetted spans are physical state. Their circular menisci share the
//! original fixed liquid volume, while each plate retains its own Young angle,
//! line friction and advancing/receding thresholds. No symmetric source is
//! synthesized. The low-speed line law is zeta_i*u_i = gamma*(cos(theta_Yi)
//! - cos(theta_i)) with a separately sourced pinning resistance; see
//! https://doi.org/10.1038/s41526-021-00186-0, equation (2), and the resisting-line
//! model discussed in https://doi.org/10.1021/acs.langmuir.3c01550.
//!
//! The gap is one discrete-gradient block; the TWO footprints form the other.
//! An average-vector-field integral within the footprint block preserves plate
//! exchange and the original symmetric reduction. An active-set Newton solve
//! admits independently pinned/advancing/receding surfaces. Analytic geometric
//! Hessians give the coupled normal-force tangent, and original surface work
//! closes with each surface's separate nonnegative line losses.
//!
//! The model remains two-dimensional, isothermal, fixed-volume, acute-angle,
//! connected and contact-line dominated. It neglects bulk-flow dissipation,
//! liquid inertia, evaporation and pinch-off; a lost neck is a refusal. Scalar
//! geometry roots and quadrature estimates are bounded numerical calculations,
//! not interval, stability or physical-trajectory certificates.

use super::BridgeLimits;
use super::contact_line::{
    ADVANCING_CONTACT_ANGLE_PROPERTY, CONTACT_LINE_FRICTION_PROPERTY, ContactLineConfig,
    RECEDING_CONTACT_ANGLE_PROPERTY,
};
pub use super::contact_line::{ContactLineError, ContactLineRegime};
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::CancelGate;
use fs_material::liquid::DYNAMIC_VISCOSITY_DIMS;
use fs_material::liquid::surface::{
    ResolvedLiquidSurfaceState, ResolvedLiquidWettingState, require_frozen_coefficient_support,
};
use fs_math::det;
use fs_qty::{Dims, QuantitySpec};

mod geometry;
mod numerics;

/// The same original-source mechanical continuation used by symmetric wetting.
pub type AsymmetricContactLineLoadLine<'a> =
    super::contact_line::load_line::DynamicContactLineLoadLine<
        'a,
        SourceAsymmetricContactLineBridge,
    >;
pub type AsymmetricContactLineLoadLineStep =
    super::contact_line::load_line::DynamicContactLineLoadLineStep<AsymmetricContactLineStep>;

impl super::contact_line::load_line::sealed::Sealed for SourceAsymmetricContactLineBridge {}
impl super::contact_line::load_line::DynamicContactLineModel for SourceAsymmetricContactLineBridge {
    type State = AsymmetricContactLineState;
    type Step = AsymmetricContactLineStep;
    fn state_gap(state: Self::State) -> f64 {
        state.gap_m()
    }
    fn advance(
        &self,
        previous: Self::State,
        gap_m: f64,
        dt_s: f64,
        gate: &CancelGate,
    ) -> Result<Self::Step, ContactLineError> {
        SourceAsymmetricContactLineBridge::advance(self, previous, gap_m, dt_s, gate)
    }
    fn step_view(step: &Self::Step) -> super::contact_line::load_line::DynamicContactLineStepView {
        super::contact_line::load_line::DynamicContactLineStepView {
            gap_m: step.state().gap_m(),
            mean_normal_force_n: step.mean_normal_force_n,
            mean_force_gap_derivative_n_per_m: step.mean_force_gap_derivative_n_per_m,
            estimated_force_error_n: step.estimated_force_error_n,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AsymmetricContactLineConfig {
    /// Common physical/numerical ceilings; both original sources remain distinct.
    pub contact_line: ContactLineConfig,
    pub maximum_backtracks: usize,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AsymmetricContactLineProfile {
    pub gap_m: f64,
    pub wetted_spans_m: [f64; 2],
    pub contact_angles_radians: [f64; 2],
    pub meniscus_radius_m: f64,
    pub neck_width_m: f64,
    pub neck_height_m: f64,
    pub liquid_volume_m3: f64,
    pub liquid_gas_area_m2: f64,
    pub wetted_areas_m2: [f64; 2],
    pub liquid_minus_gas_pressure_pa: f64,
    pub surface_free_energy_j: f64,
    /// Force on the upper plate, along increasing gap.
    pub normal_force_n: f64,
    /// Actual lower/upper resultants in the SAME increasing-gap direction.
    pub plate_normal_forces_n: [f64; 2],
    pub footprint_energy_gradients_n: [f64; 2],
    pub geometric_area_residual_m2: f64,
    pub estimated_volume_error_m3: f64,
    pub estimated_surface_energy_error_j: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AsymmetricContactLineState {
    identity: ContentHash,
    profile: AsymmetricContactLineProfile,
}
impl AsymmetricContactLineState {
    pub const fn gap_m(self) -> f64 {
        self.profile.gap_m
    }
    pub const fn wetted_spans_m(self) -> [f64; 2] {
        self.profile.wetted_spans_m
    }
    pub const fn contact_angles_radians(self) -> [f64; 2] {
        self.profile.contact_angles_radians
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AsymmetricContactLineStep {
    state: AsymmetricContactLineState,
    pub profile: AsymmetricContactLineProfile,
    pub regimes: [ContactLineRegime; 2],
    pub mean_normal_force_n: f64,
    pub mean_force_gap_derivative_n_per_m: f64,
    pub surface_energy_change_j: f64,
    pub contact_line_dissipation_j: f64,
    pub contact_line_dissipation_by_surface_j: [f64; 2],
    pub viscous_line_dissipation_by_surface_j: [f64; 2],
    pub pinning_dissipation_by_surface_j: [f64; 2],
    pub work_on_liquid_j: f64,
    pub energy_balance_residual_j: f64,
    pub footprint_force_residuals_n: [f64; 2],
    pub estimated_force_error_n: f64,
    pub estimated_energy_error_j: f64,
    pub iterations: usize,
    /// Profile/inverse-geometry solves, including all quadrature nodes. Each
    /// retains the separate geometry-iteration ceiling in the configuration.
    pub evaluations: usize,
}
impl AsymmetricContactLineStep {
    pub const fn state(self) -> AsymmetricContactLineState {
        self.state
    }
}

#[derive(Debug)]
pub struct SourceAsymmetricContactLineBridge {
    liquid: ResolvedLiquidSurfaceState,
    wetting: [ResolvedLiquidWettingState; 2],
    volume: f64,
    width: f64,
    area: f64,
    gamma: f64,
    gamma_width: f64,
    young_cosines: [f64; 2],
    zeta: [f64; 2],
    advancing_thresholds: [f64; 2],
    receding_thresholds: [f64; 2],
    limits: BridgeLimits,
    config: AsymmetricContactLineConfig,
    identity: ContentHash,
}

impl SourceAsymmetricContactLineBridge {
    pub fn new(
        liquid: ResolvedLiquidSurfaceState,
        wetting: [ResolvedLiquidWettingState; 2],
        volume_m3: f64,
        invariant_width_m: f64,
        limits: BridgeLimits,
        config: AsymmetricContactLineConfig,
    ) -> Result<Self, ContactLineError> {
        limits.validate()?;
        config.validate()?;
        if liquid.rheology_axes().is_some() {
            return Err(ContactLineError::Invalid(
                "fixed Newtonian liquid source required",
            ));
        }
        let gamma = liquid
            .surface_tension_n_m()
            .ok_or(ContactLineError::Invalid(
                "original active liquid/gas surface tension",
            ))?;
        let area = volume_m3 / invariant_width_m;
        let gamma_width = gamma * invariant_width_m;
        if [
            volume_m3,
            invariant_width_m,
            area,
            gamma_width,
            liquid.density_kg_m3() * volume_m3,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err(ContactLineError::Invalid(
                "positive representable liquid volume and extrusion width",
            ));
        }
        let mut zeta = [0.0; 2];
        let mut young_cosines = [0.0; 2];
        let mut advancing_thresholds = [0.0; 2];
        let mut receding_thresholds = [0.0; 2];
        for i in 0..2 {
            wetting[i].require_liquid_state(&liquid)?;
            let theta = wetting[i].contact_angle_radians();
            young_cosines[i] = wetting[i].young_dupre_law().map_or_else(
                || det::cos(theta),
                |law| -law.wetted_surface_energy_difference_j_m2() / gamma,
            );
            let mut values = [0.0; 3];
            for (j, (name, dims)) in [
                (CONTACT_LINE_FRICTION_PROPERTY, DYNAMIC_VISCOSITY_DIMS),
                (ADVANCING_CONTACT_ANGLE_PROPERTY, Dims::NONE),
                (RECEDING_CONTACT_ANGLE_PROPERTY, Dims::NONE),
            ]
            .into_iter()
            .enumerate()
            {
                let property =
                    wetting[i]
                        .resolved()
                        .property(name)
                        .ok_or(ContactLineError::Invalid(
                            "missing original per-surface contact-line property",
                        ))?;
                if property.requirement().quantity() != QuantitySpec::dimensional(dims) {
                    return Err(ContactLineError::Invalid(
                        "per-surface contact-line source quantity",
                    ));
                }
                require_frozen_coefficient_support(property)?;
                values[j] = property.value_si();
            }
            let c = config.contact_line;
            if values.iter().any(|x| !x.is_finite() || *x <= 0.0)
                || !(c.minimum_angle_radians <= values[2]
                    && values[2] <= theta
                    && theta <= values[1]
                    && values[1] <= c.maximum_angle_radians)
                || !(young_cosines[i] > 0.0 && young_cosines[i] < 1.0)
            {
                return Err(ContactLineError::Invalid(
                    "each surface requires ordered acute receding/Young/advancing angles and positive friction",
                ));
            }
            zeta[i] = values[0];
            advancing_thresholds[i] = if theta == values[1] {
                0.0
            } else {
                gamma_width * (young_cosines[i] - det::cos(values[1]))
            };
            receding_thresholds[i] = if theta == values[2] {
                0.0
            } else {
                gamma_width * (det::cos(values[2]) - young_cosines[i])
            };
            if !advancing_thresholds[i].is_finite()
                || advancing_thresholds[i] < 0.0
                || !receding_thresholds[i].is_finite()
                || receding_thresholds[i] < 0.0
            {
                return Err(ContactLineError::Invalid(
                    "representable per-surface pinning thresholds",
                ));
            }
        }
        let mut hash =
            DomainHasher::new("org.frankensim.fs-contact.asymmetric-contact-line-bridge.v1");
        hash.update(liquid.material().identity().as_bytes());
        for surface in &wetting {
            hash.update(surface.resolved().identity().as_bytes());
        }
        for value in [
            volume_m3,
            invariant_width_m,
            limits.minimum_gap_m,
            limits.maximum_gap_m,
            limits.maximum_wetted_span_m,
            limits.maximum_abs_force_n,
            limits.maximum_abs_surface_energy_j,
            limits.maximum_abs_step_work_j,
            limits.maximum_transition_release_j,
            limits.energy_absolute_tolerance_j,
            limits.energy_relative_tolerance,
            limits.gravity_m_s2,
            limits.maximum_bond_number,
            limits.maximum_capillary_number,
            limits.maximum_weber_number,
        ]
        .into_iter()
        .chain(config.numbers())
        {
            hash.update(&value.to_bits().to_le_bytes());
        }
        for n in [
            config.contact_line.maximum_iterations,
            config.contact_line.maximum_geometry_iterations,
            config.contact_line.maximum_evaluations,
            config.maximum_backtracks,
        ] {
            hash.update(&(n as u64).to_le_bytes());
        }
        Ok(Self {
            liquid,
            wetting,
            volume: volume_m3,
            width: invariant_width_m,
            area,
            gamma,
            gamma_width,
            young_cosines,
            zeta,
            advancing_thresholds,
            receding_thresholds,
            limits,
            config,
            identity: hash.finalize(),
        })
    }

    pub const fn liquid(&self) -> &ResolvedLiquidSurfaceState {
        &self.liquid
    }
    pub const fn wetting(&self) -> &[ResolvedLiquidWettingState; 2] {
        &self.wetting
    }
    pub const fn config(&self) -> AsymmetricContactLineConfig {
        self.config
    }
    pub const fn limits(&self) -> BridgeLimits {
        self.limits
    }
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    pub const fn line_frictions_pa_s(&self) -> [f64; 2] {
        self.zeta
    }

    /// Independent initial dynamic angles; excess surface energy is retained.
    pub fn initialize(
        &self,
        gap_m: f64,
        dynamic_angles_radians: [f64; 2],
    ) -> Result<AsymmetricContactLineState, ContactLineError> {
        let evaluation = self.from_angles(gap_m, dynamic_angles_radians)?;
        Ok(AsymmetricContactLineState {
            identity: self.identity,
            profile: evaluation.profile,
        })
    }
    /// Admit canonical footprints directly, with the same bounded geometry root.
    pub fn initialize_spans(
        &self,
        gap_m: f64,
        wetted_spans_m: [f64; 2],
        gate: &CancelGate,
    ) -> Result<AsymmetricContactLineState, ContactLineError> {
        let evaluation = self.at(gap_m, wetted_spans_m, &mut 0, gate)?;
        Ok(AsymmetricContactLineState {
            identity: self.identity,
            profile: evaluation.profile,
        })
    }
    pub fn evaluate(
        &self,
        state: AsymmetricContactLineState,
    ) -> Result<AsymmetricContactLineProfile, ContactLineError> {
        if state.identity != self.identity {
            return Err(ContactLineError::StateMismatch);
        }
        Ok(state.profile)
    }

    pub fn advance(
        &self,
        previous: AsymmetricContactLineState,
        gap_m: f64,
        dt_s: f64,
        gate: &CancelGate,
    ) -> Result<AsymmetricContactLineStep, ContactLineError> {
        poll(gate)?;
        let p = self.evaluate(previous)?;
        let c = self.config.contact_line;
        if !dt_s.is_finite() || dt_s <= 0.0 || dt_s > c.maximum_step_s {
            return Err(ContactLineError::Domain(
                "positive admitted physical timestep",
            ));
        }
        if !gap_m.is_finite() || (gap_m - p.gap_m).abs() > c.maximum_gap_increment_m {
            return Err(ContactLineError::Domain("physical gap increment"));
        }
        self.check_gap(gap_m)?;
        let b0 = p.wetted_spans_m;
        let friction = self.zeta.map(|z| z * self.width / (2.0 * dt_s));
        if friction.iter().any(|x| !x.is_finite() || *x <= 0.0) {
            return Err(ContactLineError::Invalid("finite positive line resistance"));
        }
        let mut evaluations = 0;
        let mut u = [0.0; 2];
        let mut gradient =
            self.footprint_gradient(p.gap_m, gap_m, b0, b0, &mut evaluations, gate)?;
        let mut selected = self.residual(u, friction, &gradient);
        let mut iterations = 0;
        for iteration in 0..=c.maximum_iterations {
            iterations = iteration;
            let tolerance = self.force_tolerance(
                gradient.value[0]
                    .abs()
                    .max(gradient.value[1].abs())
                    .max((friction[0] * u[0]).abs())
                    .max((friction[1] * u[1]).abs()),
            );
            let pinned_exact =
                (0..2).all(|i| selected.regimes[i] != ContactLineRegime::Pinned || u[i] == 0.0);
            let directions = (0..2).all(|i| match selected.regimes[i] {
                ContactLineRegime::Pinned => true,
                ContactLineRegime::Advancing => u[i] >= 0.0,
                ContactLineRegime::Receding => u[i] <= 0.0,
            });
            let error = gradient.error[0].max(gradient.error[1]);
            if maximum(selected.forces) + error <= tolerance && pinned_exact && directions {
                break;
            }
            if iteration == c.maximum_iterations {
                return Err(ContactLineError::NotConverged);
            }
            let correction = solve_active(
                selected.regimes,
                friction,
                gradient.jacobian,
                selected.forces.map(|x| -x),
            )?;
            let old_merit = maximum(selected.forces);
            let mut fraction = 1.0;
            let mut accepted = None;
            for _ in 0..self.config.maximum_backtracks {
                poll(gate)?;
                let candidate = core::array::from_fn(|i| {
                    if fraction == 1.0 && selected.regimes[i] == ContactLineRegime::Pinned {
                        0.0
                    } else {
                        u[i] + fraction * correction[i]
                    }
                });
                if candidate
                    .iter()
                    .any(|x| !x.is_finite() || x.abs() > c.maximum_span_increment_m)
                {
                    fraction *= 0.5;
                    continue;
                }
                let spans = core::array::from_fn(|i| b0[i] + candidate[i]);
                match self.footprint_gradient(p.gap_m, gap_m, b0, spans, &mut evaluations, gate) {
                    Ok(next) => {
                        let next_selected = self.residual(candidate, friction, &next);
                        let next_merit = maximum(next_selected.forces);
                        if next_merit < old_merit
                            || next_merit + next.error[0].max(next.error[1]) <= tolerance
                        {
                            accepted = Some((candidate, next, next_selected));
                            break;
                        }
                    }
                    Err(ContactLineError::Domain(_)) => {}
                    Err(error) => return Err(error),
                }
                fraction *= 0.5;
            }
            let (next_u, next_gradient, next_selected) =
                accepted.ok_or(ContactLineError::NotConverged)?;
            if next_u == u {
                return Err(ContactLineError::NotConverged);
            }
            u = next_u;
            gradient = next_gradient;
            selected = next_selected;
        }
        let spans = core::array::from_fn(|i| b0[i] + u[i]);
        // Work uses the actual representable canonical coordinates, including
        // any final floating-point addition to an existing footprint.
        let actual_u = core::array::from_fn(|i| spans[i] - b0[i]);
        if actual_u != u {
            u = actual_u;
            gradient =
                self.footprint_gradient(p.gap_m, gap_m, b0, spans, &mut evaluations, gate)?;
            selected = self.residual(u, friction, &gradient);
            let tolerance =
                self.force_tolerance(gradient.value[0].abs().max(gradient.value[1].abs()));
            if maximum(selected.forces) + gradient.error[0].max(gradient.error[1]) > tolerance {
                return Err(ContactLineError::NotConverged);
            }
        }
        let endpoint = self.at(gap_m, spans, &mut evaluations, gate)?;
        let normal = self.gap_gradient(p.gap_m, gap_m, b0, spans, &mut evaluations, gate)?;
        let response = solve_active(
            selected.regimes,
            friction,
            gradient.jacobian,
            core::array::from_fn(|i| {
                if selected.regimes[i] == ContactLineRegime::Pinned {
                    0.0
                } else {
                    -gradient.gap_derivative[i]
                }
            }),
        )?;
        let mean_force = -normal.value;
        let tangent = -(normal.gap_derivative
            + normal.span_derivatives[0] * response[0]
            + normal.span_derivatives[1] * response[1]);
        // The positive viscous diagonal bounds the footprint response to a
        // force perturbation on this convex acute-angle footprint branch.
        let footprint_error = (selected.forces[0].abs() + gradient.error[0])
            .hypot(selected.forces[1].abs() + gradient.error[1])
            / friction[0].min(friction[1]);
        let force_error = normal.error
            + normal.span_derivatives[0].hypot(normal.span_derivatives[1]) * footprint_error;
        let force_tolerance = self.force_tolerance(mean_force.abs());
        if !tangent.is_finite()
            || !mean_force.is_finite()
            || !force_error.is_finite()
            || force_error > force_tolerance
        {
            return Err(ContactLineError::Accuracy {
                estimated_error_n: force_error,
                tolerance_n: force_tolerance,
            });
        }
        let viscous = core::array::from_fn(|i| friction[i] * u[i] * u[i]);
        let pinning = core::array::from_fn(|i| selected.thresholds[i] * u[i]);
        let losses: [f64; 2] = core::array::from_fn(|i| viscous[i] + pinning[i]);
        if viscous
            .iter()
            .chain(pinning.iter())
            .any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err(ContactLineError::Invalid(
                "nonnegative per-surface line work",
            ));
        }
        let dissipation = losses[0] + losses[1];
        let dh = gap_m - p.gap_m;
        let work = normal.value * dh;
        let change = endpoint.profile.surface_free_energy_j - p.surface_free_energy_j;
        let energy_residual = change + dissipation - work;
        let energy_error = p.estimated_surface_energy_error_j
            + endpoint.profile.estimated_surface_energy_error_j
            + normal.error * dh.abs()
            + gradient.error[0] * u[0].abs()
            + gradient.error[1] * u[1].abs();
        let scale = p
            .surface_free_energy_j
            .abs()
            .max(endpoint.profile.surface_free_energy_j.abs())
            .max(work.abs())
            .max(dissipation);
        let tolerance =
            self.limits.energy_absolute_tolerance_j + self.limits.energy_relative_tolerance * scale;
        if !work.is_finite()
            || !energy_error.is_finite()
            || energy_residual.abs() + energy_error > tolerance
        {
            return Err(ContactLineError::Energy {
                residual_j: energy_residual,
                tolerance_j: tolerance,
            });
        }
        super::limit(
            mean_force.abs(),
            self.limits.maximum_abs_force_n,
            "asymmetric wetting mean force",
        )?;
        super::limit(
            work.abs(),
            self.limits.maximum_abs_step_work_j,
            "asymmetric wetting work",
        )?;
        super::limit(
            dissipation,
            self.limits.maximum_abs_step_work_j,
            "asymmetric wetting line loss",
        )?;
        let speeds = u.map(|v| (v / (2.0 * dt_s)).abs());
        if speeds.iter().any(|v| *v > c.maximum_contact_line_speed_m_s) {
            return Err(ContactLineError::Domain("per-surface contact-line speed"));
        }
        let speed = speeds[0].max(speeds[1]).max((dh / dt_s).abs());
        let length = spans
            .into_iter()
            .chain(b0)
            .fold(gap_m.max(p.gap_m), f64::max);
        super::limit(
            self.liquid.dynamic_viscosity_pa_s() * speed / self.gamma,
            self.limits.maximum_capillary_number,
            "asymmetric wetting capillary number",
        )?;
        super::limit(
            self.liquid.density_kg_m3() * speed * speed * length / self.gamma,
            self.limits.maximum_weber_number,
            "asymmetric wetting Weber number",
        )?;
        poll(gate)?;
        Ok(AsymmetricContactLineStep {
            state: AsymmetricContactLineState {
                identity: self.identity,
                profile: endpoint.profile,
            },
            profile: endpoint.profile,
            regimes: selected.regimes,
            mean_normal_force_n: mean_force,
            mean_force_gap_derivative_n_per_m: tangent,
            surface_energy_change_j: change,
            contact_line_dissipation_j: dissipation,
            contact_line_dissipation_by_surface_j: losses,
            viscous_line_dissipation_by_surface_j: viscous,
            pinning_dissipation_by_surface_j: pinning,
            work_on_liquid_j: work,
            energy_balance_residual_j: energy_residual,
            footprint_force_residuals_n: selected.forces,
            estimated_force_error_n: force_error,
            estimated_energy_error_j: energy_error,
            iterations,
            evaluations,
        })
    }

    fn residual(
        &self,
        u: [f64; 2],
        friction: [f64; 2],
        g: &numerics::FootprintGradient,
    ) -> Selected {
        let mut selected = Selected {
            regimes: [ContactLineRegime::Pinned; 2],
            thresholds: [0.0; 2],
            forces: [0.0; 2],
        };
        for i in 0..2 {
            if g.value[i] < -self.advancing_thresholds[i] {
                selected.regimes[i] = ContactLineRegime::Advancing;
                selected.thresholds[i] = self.advancing_thresholds[i];
            } else if g.value[i] > self.receding_thresholds[i] {
                selected.regimes[i] = ContactLineRegime::Receding;
                selected.thresholds[i] = -self.receding_thresholds[i];
            } else {
                selected.thresholds[i] = -g.value[i];
            }
            selected.forces[i] = friction[i] * u[i]
                + if selected.regimes[i] == ContactLineRegime::Pinned {
                    0.0
                } else {
                    g.value[i] + selected.thresholds[i]
                };
        }
        selected
    }
    fn force_tolerance(&self, scale: f64) -> f64 {
        self.config.contact_line.force_absolute_tolerance_n
            + self.config.contact_line.force_relative_tolerance * scale
    }
    fn check_gap(&self, gap: f64) -> Result<(), ContactLineError> {
        if !gap.is_finite() || gap < self.limits.minimum_gap_m || gap > self.limits.maximum_gap_m {
            return Err(ContactLineError::Domain("positive admitted plate gap"));
        }
        Ok(())
    }
}

struct Selected {
    regimes: [ContactLineRegime; 2],
    thresholds: [f64; 2],
    forces: [f64; 2],
}

fn solve_active(
    regimes: [ContactLineRegime; 2],
    friction: [f64; 2],
    hessian: [[f64; 2]; 2],
    rhs: [f64; 2],
) -> Result<[f64; 2], ContactLineError> {
    let mut a = [[0.0; 2]; 2];
    for i in 0..2 {
        if regimes[i] != ContactLineRegime::Pinned {
            a[i] = hessian[i];
        }
        a[i][i] += friction[i];
    }
    let schur = a[1][1] - (a[1][0] / a[0][0]) * a[0][1];
    if !a[0][0].is_finite() || a[0][0] <= 0.0 || !schur.is_finite() || schur <= 0.0 {
        return Err(ContactLineError::Invalid(
            "positive contact-line active Jacobian",
        ));
    }
    let y = (rhs[1] - a[1][0] * (rhs[0] / a[0][0])) / schur;
    let x = (rhs[0] - a[0][1] * y) / a[0][0];
    if !x.is_finite() || !y.is_finite() {
        return Err(ContactLineError::NotConverged);
    }
    Ok([x, y])
}

impl AsymmetricContactLineConfig {
    fn numbers(self) -> [f64; 12] {
        let c = self.contact_line;
        [
            c.minimum_angle_radians,
            c.maximum_angle_radians,
            c.minimum_neck_width_m,
            c.maximum_gap_increment_m,
            c.maximum_span_increment_m,
            c.maximum_step_s,
            c.maximum_contact_line_speed_m_s,
            c.maximum_angle_departure_radians,
            c.geometry_absolute_tolerance_m,
            c.geometry_relative_tolerance,
            c.force_absolute_tolerance_n,
            c.force_relative_tolerance,
        ]
    }
    fn validate(self) -> Result<(), ContactLineError> {
        let c = self.contact_line;
        if self.numbers().iter().any(|v| !v.is_finite() || *v <= 0.0)
            || c.minimum_angle_radians >= c.maximum_angle_radians
            || c.maximum_angle_radians >= core::f64::consts::FRAC_PI_2
            || c.geometry_relative_tolerance >= 1.0
            || c.force_relative_tolerance >= 1.0
            || !(1..=128).contains(&c.maximum_iterations)
            || !(1..=128).contains(&c.maximum_geometry_iterations)
            || c.maximum_evaluations == 0
            || !(1..=64).contains(&self.maximum_backtracks)
        {
            return Err(ContactLineError::Invalid(
                "positive asymmetric wetting physical/numerical budgets",
            ));
        }
        Ok(())
    }
}
fn maximum(v: [f64; 2]) -> f64 {
    v[0].abs().max(v[1].abs())
}
fn poll(gate: &CancelGate) -> Result<(), ContactLineError> {
    if gate.is_requested() {
        Err(ContactLineError::Cancelled)
    } else {
        Ok(())
    }
}
