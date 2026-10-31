//! Inertial circular-capillary segments with a discrete mechanical energy balance.
//!
//! Both incompressible phase columns share one axial velocity. Their effective
//! mass is affine in wet length, with an explicitly supplied additional mass.
//! The continuous inertial force is M(x)*du/dt + M'(x)*u^2/2. Reservoir exchange
//! is separately declared as C*|u|*u, whose loss is C*|u|^3; no hidden entry-loss
//! coefficient is inferred. For one liquid C=rho*A/2 recovers the bidirectional
//! reservoir convention of Fricke et al. (2023), doi:10.1016/j.physd.2023.133895.
//!
//! A velocity-midpoint discrete gradient preserves the actual endpoint kinetic
//! energy increment, Young wall energy and the open-system gravitational
//! account. Bulk, contact-line, hysteresis and reservoir losses stay separate.
//! A velocity zero is a resolved *discrete* stopping event. The caller owns any
//! remaining time, and must restart from rest to pin or select the next branch.
//! Tube-end crossings are refused rather than clipped. This is a lumped
//! inertance closure, not resolved developing-flow profiles or arbitrary
//! contact-line physics.

use super::{ContactLineMotion, InclinedTwoFluidCircularMeniscusInput};
use core::f64::consts::PI;

/// Explicit low-speed contact-line coefficients and their angular domain.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InertialContactLine {
    pub friction_pa_s: f64,
    pub advancing_angle_radians: f64,
    pub receding_angle_radians: f64,
    pub maximum_angle_departure_radians: f64,
}

/// An original two-fluid tube plus physical velocity and explicit end models.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InertialCapillaryInput {
    pub meniscus: InclinedTwoFluidCircularMeniscusInput,
    pub initial_velocity_m_s: f64,
    /// Additional moving mass [kg], not an implicit dry-start regularization.
    pub added_mass_kg: f64,
    /// Advancing reservoir force C*|u|*u; C has units kg/m.
    pub advancing_reservoir_loss_kg_m: f64,
    /// Receding reservoir force C*|u|*u; C has units kg/m.
    pub receding_reservoir_loss_kg_m: f64,
    pub maximum_velocity_m_s: f64,
    pub contact_line: Option<InertialContactLine>,
}

/// Bounded pressure evaluations and independent residual/energy tolerances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InertialCapillarySolve {
    pub maximum_iterations: usize,
    pub pressure_absolute_tolerance_pa: f64,
    pub velocity_absolute_tolerance_m_s: f64,
    pub time_absolute_tolerance_s: f64,
    pub energy_absolute_tolerance_j: f64,
    pub energy_relative_tolerance: f64,
}

/// Why this numerical segment ended.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum InertialCapillaryEnd {
    /// The complete requested duration was consumed.
    TimeLimit,
    /// The endpoint velocity is exactly zero; the caller owns remaining time.
    Stopped,
}

/// State and distinct energy channels for the actually consumed segment.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InertialCapillaryStep {
    pub initial_wet_length_m: f64,
    pub wet_length_m: f64,
    pub initial_velocity_m_s: f64,
    pub velocity_m_s: f64,
    pub elapsed_time_s: f64,
    pub end: InertialCapillaryEnd,
    pub transferred_volume_m3: f64,
    pub initial_effective_mass_kg: f64,
    pub effective_mass_kg: f64,
    pub initial_kinetic_energy_j: f64,
    pub kinetic_energy_j: f64,
    pub kinetic_energy_change_j: f64,
    /// Input pressure times signed transferred volume. Absolute for a sealed owner.
    pub inlet_pressure_work_j: f64,
    /// Downstream midpoint pressure times signed transferred volume.
    pub downstream_pressure_work_j: f64,
    /// Net pressure work, with the pressure difference formed before multiplication.
    pub pressure_work_j: f64,
    pub wall_surface_energy_change_j: f64,
    pub gravitational_storage_change_j: f64,
    pub outlet_potential_energy_j: f64,
    pub liquid_viscous_dissipation_j: f64,
    pub displaced_viscous_dissipation_j: f64,
    pub contact_line_viscous_dissipation_j: f64,
    pub pinning_dissipation_j: f64,
    pub reservoir_dissipation_j: f64,
    pub total_dissipation_j: f64,
    pub motion: ContactLineMotion,
    pub initial_apparent_angle_radians: f64,
    pub final_apparent_angle_radians: f64,
    pub pressure_residual_pa: f64,
    pub energy_balance_residual_j: f64,
    pub iterations: usize,
}

/// No caller state is modified by any refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum InertialCapillaryError {
    Invalid(&'static str),
    Cancelled,
    StepTooLarge(&'static str),
    NotConverged { iterations: usize, residual_pa: f64 },
    EnergyBalance { residual_j: f64, tolerance_j: f64 },
}
impl core::fmt::Display for InertialCapillaryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "inertial capillary: {self:?}")
    }
}
impl std::error::Error for InertialCapillaryError {}

/// Original downstream owner's discrete pressure and analytic sensitivities.
pub(crate) struct Pressure<T> {
    pub downstream_pressure_pa: f64,
    /// Actual downstream pressure at the candidate endpoint, for stopping admission.
    pub endpoint_pressure_pa: f64,
    /// dP_mid/dx_final, holding the initial state and elapsed time fixed.
    pub derivative_wet_length_pa_m: f64,
    /// dP_mid/dt, holding the initial state and final wet length fixed.
    pub derivative_time_pa_s: f64,
    pub payload: T,
}
pub(crate) struct InertialPressureStep<T> {
    pub step: InertialCapillaryStep,
    pub pressure: Pressure<T>,
}
#[derive(Debug, PartialEq)]
pub(crate) enum InertialSolveError<E> {
    Physics(InertialCapillaryError),
    Pressure(E),
}
impl<E> From<InertialCapillaryError> for InertialSolveError<E> {
    fn from(value: InertialCapillaryError) -> Self {
        Self::Physics(value)
    }
}

/// Advance one constant-boundary segment, resolving a velocity-zero event.
///
/// The midpoint map is second order on smooth moving branches. Its stopping
/// time is a numerical event time, not the exact continuous stopping time.
/// Resolve temporal accuracy with subdivision; a successful segment can consume
/// less than the requested time only when it returns `Stopped`.
pub fn step_inertial_circular_meniscus(
    input: InertialCapillaryInput,
    controls: InertialCapillarySolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<InertialCapillaryStep, InertialCapillaryError> {
    let result = solve_with_pressure(
        input,
        controls,
        |_, _, _| {
            Ok::<_, core::convert::Infallible>(Pressure {
                downstream_pressure_pa: 0.0,
                endpoint_pressure_pa: 0.0,
                derivative_wet_length_pa_m: 0.0,
                derivative_time_pa_s: 0.0,
                payload: (),
            })
        },
        &mut cancelled,
    );
    match result {
        Ok(result) => Ok(result.step),
        Err(InertialSolveError::Physics(error)) => Err(error),
        Err(InertialSolveError::Pressure(impossible)) => match impossible {},
    }
}

/// Common mechanics seam. The callback owns downstream storage and heat.
pub(crate) fn solve_with_pressure<T, E, C: FnMut() -> bool>(
    input: InertialCapillaryInput,
    controls: InertialCapillarySolve,
    mut callback: impl FnMut(f64, f64, &mut C) -> Result<Pressure<T>, E>,
    cancelled: &mut C,
) -> Result<InertialPressureStep<T>, InertialSolveError<E>> {
    poll(cancelled)?;
    let model = Model::new(input, controls)?;
    let h = input.meniscus.meniscus.meniscus.duration_s;
    let u0 = input.initial_velocity_m_s;
    let mut budget = Budget {
        used: 0,
        maximum: controls.maximum_iterations,
        residual_pa: 0.0,
    };
    let mut seed = if u0 == 0.0 || h == 0.0 {
        Some(pressure_at(
            model.x0(),
            h,
            &mut callback,
            cancelled,
            &mut budget,
        )?)
    } else {
        None
    };
    let branch = model.branch(seed.as_ref().map(|p| p.downstream_pressure_pa))?;
    if h == 0.0 {
        return finish(
            model,
            branch,
            0.0,
            u0,
            seed.expect("identity pressure"),
            InertialCapillaryEnd::TimeLimit,
            0.0,
            budget.used,
            cancelled,
        );
    }
    if branch.motion == ContactLineMotion::Pinned {
        let pressure = seed.expect("pinned pressure");
        // A midpoint reaction is insufficient if heating has already taken the
        // endpoint outside the static interval. The caller must subdivide.
        let endpoint_drive = model.young_drive(model.x0(), pressure.endpoint_pressure_pa);
        if endpoint_drive > model.advancing_pressure + controls.pressure_absolute_tolerance_pa
            || endpoint_drive < -model.receding_pressure - controls.pressure_absolute_tolerance_pa
        {
            return Err(InertialCapillaryError::StepTooLarge(
                "pinned endpoint crossed a contact-angle threshold",
            )
            .into());
        }
        return finish(
            model,
            branch,
            h,
            0.0,
            pressure,
            InertialCapillaryEnd::TimeLimit,
            0.0,
            budget.used,
            cancelled,
        );
    }
    let sign = branch.sign;
    let base = 0.5 * u0;
    let (_, _) = model.state(h, base)?;
    let initial_pressure = if let Some(p) = seed.take() {
        p
    } else {
        let (x, _) = model.state(h, base)?;
        pressure_at(x, h, &mut callback, cancelled, &mut budget)?
    };
    let (residual, tangent) = model.residual(branch, h, base, &initial_pressure)?;
    budget.residual_pa = residual / model.area;
    if u0 != 0.0 && (sign * residual >= 0.0 || converged_velocity(residual, tangent, model)) {
        return stopping(
            model,
            branch,
            h,
            initial_pressure,
            residual,
            &mut callback,
            cancelled,
            &mut budget,
        );
    }
    if converged_velocity(residual, tangent, model) {
        return finish(
            model,
            branch,
            h,
            base,
            initial_pressure,
            if u0 != 0.0 {
                InertialCapillaryEnd::Stopped
            } else {
                InertialCapillaryEnd::TimeLimit
            },
            residual,
            budget.used,
            cancelled,
        );
    }
    let raw_limit = if sign > 0.0 {
        (0.5 * u0 + 0.5 * input.maximum_velocity_m_s)
            .min(((model.length() - model.x0()) / h).next_down())
    } else {
        (0.5 * u0 - 0.5 * input.maximum_velocity_m_s).max((-model.x0() / h).next_up())
    };
    // Only continue the increasing residual branch connected to the initial
    // state. The known affine derivative is bounded below before downstream
    // stiffness is added; negative density/viscosity contrasts can require dt
    // refinement instead of selecting a remote algebraic root.
    let derivative0 = model.base_tangent(branch, h, base);
    if !derivative0.is_finite() || derivative0 <= 0.0 {
        return Err(InertialCapillaryError::StepTooLarge(
            "cannot continue a stable increasing inertial residual",
        )
        .into());
    }
    let slope = model.derivative_slope(branch, h);
    let mut limit = raw_limit;
    if sign * slope < 0.0 {
        let restricted = base - 0.5 * derivative0 / slope;
        limit = if sign > 0.0 {
            limit.min(restricted)
        } else {
            limit.max(restricted)
        };
    }
    if !limit.is_finite() || sign * (limit - base) <= 0.0 {
        return Err(InertialCapillaryError::StepTooLarge(
            "velocity, tube or monotonicity limit reached; reduce dt",
        )
        .into());
    }
    let mut near = base;
    let mut far = limit;
    let mut current = base;
    let mut residual = residual;
    let mut tangent = tangent;
    loop {
        let newton = current - residual / tangent;
        let low = near.min(far);
        let high = near.max(far);
        let next = if newton.is_finite() && newton > low && newton < high {
            newton
        } else {
            near + 0.5 * (far - near)
        };
        if next == current || next == near || next == far {
            return Err(InertialCapillaryError::StepTooLarge(
                "inertial root exhausted its admitted velocity interval",
            )
            .into());
        }
        let (x, _) = model.state(h, next)?;
        let pressure = pressure_at(x, h, &mut callback, cancelled, &mut budget)?;
        let (r, d) = model.residual(branch, h, next, &pressure)?;
        budget.residual_pa = r / model.area;
        if converged_velocity(r, d, model) {
            return finish(
                model,
                branch,
                h,
                next,
                pressure,
                InertialCapillaryEnd::TimeLimit,
                r,
                budget.used,
                cancelled,
            );
        }
        if sign * r < 0.0 {
            near = next;
        } else {
            far = next;
        }
        current = next;
        residual = r;
        tangent = d;
    }
}

#[derive(Clone, Copy)]
struct Model {
    input: InertialCapillaryInput,
    controls: InertialCapillarySolve,
    area: f64,
    mass0: f64,
    mass_slope: f64,
    young_pressure: f64,
    head0: f64,
    head_slope: f64,
    bulk0: f64,
    bulk_slope: f64,
    line_resistance: f64,
    advancing_pressure: f64,
    receding_pressure: f64,
}
#[derive(Clone, Copy)]
struct Branch {
    sign: f64,
    pin_force: f64,
    reservoir: f64,
    angle_cosine: f64,
    motion: ContactLineMotion,
}
impl Model {
    fn x0(self) -> f64 {
        self.input.meniscus.meniscus.meniscus.wet_length_m
    }
    fn length(self) -> f64 {
        self.input.meniscus.meniscus.meniscus.tube_length_m
    }
    fn young_drive(self, x: f64, downstream: f64) -> f64 {
        let m = self.input.meniscus.meniscus.meniscus;
        ((m.applied_pressure_pa - downstream) + self.young_pressure)
            - (self.head0 + self.head_slope * (x - self.x0()))
    }
    fn new(
        input: InertialCapillaryInput,
        c: InertialCapillarySolve,
    ) -> Result<Self, InertialCapillaryError> {
        let m = input.meniscus.meniscus.meniscus;
        let rho = input.meniscus.density_kg_m3;
        let rho_d = input.meniscus.displaced_density_kg_m3;
        let mu_d = input.meniscus.meniscus.displaced_dynamic_viscosity_pa_s;
        let g = input.meniscus.opposing_gravity_m_s2;
        if ![
            m.radius_m,
            m.tube_length_m,
            rho,
            input.maximum_velocity_m_s,
            c.pressure_absolute_tolerance_pa,
            c.velocity_absolute_tolerance_m_s,
            c.time_absolute_tolerance_s,
            c.energy_absolute_tolerance_j,
        ]
        .into_iter()
        .all(|x| x.is_finite() && x > 0.0)
            || ![
                m.dynamic_viscosity_pa_s,
                mu_d,
                rho_d,
                m.surface_tension_n_m,
                m.duration_s,
                input.added_mass_kg,
                input.advancing_reservoir_loss_kg_m,
                input.receding_reservoir_loss_kg_m,
                c.energy_relative_tolerance,
            ]
            .into_iter()
            .all(|x| x.is_finite() && x >= 0.0)
            || ![
                m.wet_length_m,
                m.applied_pressure_pa,
                g,
                input.initial_velocity_m_s,
            ]
            .into_iter()
            .all(f64::is_finite)
            || !(-1.0..=1.0).contains(&m.contact_angle_cosine)
            || m.wet_length_m < 0.0
            || m.wet_length_m > m.tube_length_m
            || input.initial_velocity_m_s.abs() > input.maximum_velocity_m_s
            || !(1..=256).contains(&c.maximum_iterations)
        {
            return Err(InertialCapillaryError::Invalid(
                "finite tube/state, nonnegative losses, positive bounds and 1..=256 evaluations required",
            ));
        }
        let area = PI * m.radius_m * m.radius_m;
        let mass0 = area * (rho * m.wet_length_m + rho_d * (m.tube_length_m - m.wet_length_m))
            + input.added_mass_kg;
        let mass_slope = area * (rho - rho_d);
        let young_pressure = (2.0 * m.surface_tension_n_m / m.radius_m) * m.contact_angle_cosine;
        let head0 = g * (rho * m.wet_length_m + rho_d * (m.tube_length_m - m.wet_length_m));
        let head_slope = g * (rho - rho_d);
        let bulk0 = 8.0
            * PI
            * (m.dynamic_viscosity_pa_s * m.wet_length_m
                + mu_d * (m.tube_length_m - m.wet_length_m));
        let bulk_slope = 8.0 * PI * (m.dynamic_viscosity_pa_s - mu_d);
        let (line_resistance, advancing_pressure, receding_pressure) =
            if let Some(line) = input.contact_line {
                let young = m.contact_angle_cosine.acos();
                // acos(cos(theta)) can differ by a few ulps even for the same
                // supplied angle. Permit only that representation-scale equality.
                let angle_roundoff = 8.0 * f64::EPSILON * young.abs().max(1.0);
                if !line.friction_pa_s.is_finite()
                    || line.friction_pa_s < 0.0
                    || ![
                        line.advancing_angle_radians,
                        line.receding_angle_radians,
                        line.maximum_angle_departure_radians,
                    ]
                    .into_iter()
                    .all(|x| x.is_finite() && x > 0.0)
                    || line.receding_angle_radians > line.advancing_angle_radians
                    || line.receding_angle_radians > young + angle_roundoff
                    || young > line.advancing_angle_radians + angle_roundoff
                    || line.advancing_angle_radians >= PI
                    || line.maximum_angle_departure_radians > PI
                    || m.surface_tension_n_m <= 0.0
                {
                    return Err(InertialCapillaryError::Invalid(
                        "ordered sourced contact-line coefficients and positive tension required",
                    ));
                }
                let scale = 2.0 * m.surface_tension_n_m / m.radius_m;
                (
                    2.0 * PI * m.radius_m * line.friction_pa_s,
                    scale * (m.contact_angle_cosine - line.advancing_angle_radians.cos()).max(0.0),
                    scale * (line.receding_angle_radians.cos() - m.contact_angle_cosine).max(0.0),
                )
            } else {
                (0.0, 0.0, 0.0)
            };
        if ![area, mass0].into_iter().all(|x| x.is_finite() && x > 0.0)
            || ![
                mass_slope,
                young_pressure,
                head0,
                head_slope,
                bulk0,
                bulk_slope,
                line_resistance,
                advancing_pressure,
                receding_pressure,
            ]
            .into_iter()
            .all(f64::is_finite)
            || advancing_pressure < 0.0
            || receding_pressure < 0.0
        {
            return Err(InertialCapillaryError::Invalid(
                "effective mass or constitutive scales are unrepresentable",
            ));
        }
        Ok(Self {
            input,
            controls: c,
            area,
            mass0,
            mass_slope,
            young_pressure,
            head0,
            head_slope,
            bulk0,
            bulk_slope,
            line_resistance,
            advancing_pressure,
            receding_pressure,
        })
    }
    fn branch(self, at_rest_pressure: Option<f64>) -> Result<Branch, InertialCapillaryError> {
        let m = self.input.meniscus.meniscus.meniscus;
        let u = self.input.initial_velocity_m_s;
        let force = at_rest_pressure.map_or(0.0, |p| self.young_drive(self.x0(), p));
        let motion = if u > 0.0 || (u == 0.0 && force > self.advancing_pressure) {
            ContactLineMotion::Advancing
        } else if u < 0.0 || (u == 0.0 && force < -self.receding_pressure) {
            ContactLineMotion::Receding
        } else {
            ContactLineMotion::Pinned
        };
        let (sign, threshold, reservoir, cosine) = match motion {
            ContactLineMotion::Advancing => (
                1.0,
                self.advancing_pressure,
                self.input.advancing_reservoir_loss_kg_m,
                self.input
                    .contact_line
                    .map_or(m.contact_angle_cosine, |l| l.advancing_angle_radians.cos()),
            ),
            ContactLineMotion::Receding => (
                -1.0,
                -self.receding_pressure,
                self.input.receding_reservoir_loss_kg_m,
                self.input
                    .contact_line
                    .map_or(m.contact_angle_cosine, |l| l.receding_angle_radians.cos()),
            ),
            ContactLineMotion::Pinned => (
                0.0,
                force,
                0.0,
                if self.input.contact_line.is_some() {
                    (m.contact_angle_cosine - force / (2.0 * m.surface_tension_n_m / m.radius_m))
                        .clamp(-1.0, 1.0)
                } else {
                    m.contact_angle_cosine
                },
            ),
        };
        if !force.is_finite() || !threshold.is_finite() || !(self.area * threshold).is_finite() {
            return Err(InertialCapillaryError::Invalid(
                "contact-line branch force is unrepresentable",
            ));
        }
        Ok(Branch {
            sign,
            pin_force: self.area * threshold,
            reservoir,
            angle_cosine: cosine,
            motion,
        })
    }
    fn state(self, h: f64, mean_velocity: f64) -> Result<(f64, f64), InertialCapillaryError> {
        let displacement = h * mean_velocity;
        let x = self.x0() + displacement;
        let u = mean_velocity + (mean_velocity - self.input.initial_velocity_m_s);
        if ![h, mean_velocity, displacement, x, u]
            .into_iter()
            .all(f64::is_finite)
            || (h > 0.0 && mean_velocity != 0.0 && (displacement == 0.0 || x == self.x0()))
        {
            return Err(InertialCapillaryError::Invalid(
                "inertial displacement/velocity is unrepresentable",
            ));
        }
        if x < 0.0
            || x > self.length()
            || (displacement > 0.0 && x == self.length())
            || (displacement < 0.0 && x == 0.0)
        {
            return Err(InertialCapillaryError::StepTooLarge(
                "tube-end event requires its boundary model",
            ));
        }
        if u.abs() > self.input.maximum_velocity_m_s {
            return Err(InertialCapillaryError::StepTooLarge(
                "declared velocity domain exceeded",
            ));
        }
        Ok((x, u))
    }
    fn derivative_slope(self, b: Branch, h: f64) -> f64 {
        4.0 * self.mass_slope + h * self.bulk_slope + 2.0 * b.sign * b.reservoir
    }
    fn base_tangent(self, b: Branch, h: f64, w: f64) -> f64 {
        let u1 = w + (w - self.input.initial_velocity_m_s);
        2.0 * self.mass0 / h
            + 2.0 * self.mass_slope * u1
            + self.bulk0
            + self.line_resistance
            + h * self.bulk_slope * w
            + 2.0 * b.reservoir * w.abs()
            + 0.5 * h * self.area * self.head_slope
    }
    fn residual<T>(
        self,
        b: Branch,
        h: f64,
        w: f64,
        pressure: &Pressure<T>,
    ) -> Result<(f64, f64), InertialCapillaryError> {
        let (_, u1) = self.state(h, w)?;
        let inertia = (2.0 * self.mass0 / h) * (w - self.input.initial_velocity_m_s)
            + 0.5 * self.mass_slope * u1 * u1;
        let xbar = self.x0() + 0.5 * h * w;
        let bulk = self.bulk0 + self.bulk_slope * (xbar - self.x0());
        let force =
            inertia + (bulk + self.line_resistance) * w + b.pin_force + b.reservoir * w.abs() * w
                - self.area * self.young_drive(xbar, pressure.downstream_pressure_pa);
        let derivative =
            self.base_tangent(b, h, w) + self.area * h * pressure.derivative_wet_length_pa_m;
        if !force.is_finite() || !derivative.is_finite() || derivative <= 0.0 || bulk < 0.0 {
            return Err(InertialCapillaryError::StepTooLarge(
                "inertial residual or positive tangent cannot be represented",
            ));
        }
        Ok((force, derivative))
    }
    fn angle(self, b: Branch, u: f64) -> Result<f64, InertialCapillaryError> {
        let m = self.input.meniscus.meniscus.meniscus;
        if let Some(line) = self.input.contact_line {
            let cosine = b.angle_cosine - line.friction_pa_s * u / m.surface_tension_n_m;
            if !cosine.is_finite() || cosine.abs() > 1.0 + 32.0 * f64::EPSILON {
                return Err(InertialCapillaryError::Invalid(
                    "moving apparent contact angle leaves its domain",
                ));
            }
            let angle = cosine.clamp(-1.0, 1.0).acos();
            if (angle - m.contact_angle_cosine.acos()).abs() > line.maximum_angle_departure_radians
            {
                return Err(InertialCapillaryError::Invalid(
                    "apparent-angle departure exceeds its declared bound",
                ));
            }
            Ok(angle)
        } else {
            Ok(m.contact_angle_cosine.acos())
        }
    }
}

struct Budget {
    used: usize,
    maximum: usize,
    residual_pa: f64,
}
fn poll<E>(cancelled: &mut impl FnMut() -> bool) -> Result<(), InertialSolveError<E>> {
    if cancelled() {
        Err(InertialCapillaryError::Cancelled.into())
    } else {
        Ok(())
    }
}
fn pressure_at<T, E, C: FnMut() -> bool>(
    x: f64,
    h: f64,
    callback: &mut impl FnMut(f64, f64, &mut C) -> Result<Pressure<T>, E>,
    cancelled: &mut C,
    budget: &mut Budget,
) -> Result<Pressure<T>, InertialSolveError<E>> {
    poll(cancelled)?;
    if budget.used == budget.maximum {
        return Err(InertialCapillaryError::NotConverged {
            iterations: budget.used,
            residual_pa: budget.residual_pa,
        }
        .into());
    }
    budget.used += 1;
    let p = callback(x, h, cancelled).map_err(InertialSolveError::Pressure)?;
    if ![
        p.downstream_pressure_pa,
        p.endpoint_pressure_pa,
        p.derivative_wet_length_pa_m,
    ]
    .into_iter()
    .all(|v| v.is_finite() && v >= 0.0)
        || !p.derivative_time_pa_s.is_finite()
    {
        return Err(InertialCapillaryError::Invalid(
            "downstream pressure or stable discrete sensitivity is invalid",
        )
        .into());
    }
    Ok(p)
}
fn converged_velocity(residual: f64, tangent: f64, m: Model) -> bool {
    residual.abs() / m.area <= m.controls.pressure_absolute_tolerance_pa
        && (residual == 0.0
            || (residual / tangent).abs() <= m.controls.velocity_absolute_tolerance_m_s)
}
fn stopping<T, E, C: FnMut() -> bool>(
    m: Model,
    b: Branch,
    requested: f64,
    mut pressure: Pressure<T>,
    mut residual: f64,
    callback: &mut impl FnMut(f64, f64, &mut C) -> Result<Pressure<T>, E>,
    cancelled: &mut C,
    budget: &mut Budget,
) -> Result<InertialPressureStep<T>, InertialSolveError<E>> {
    let w = 0.5 * m.input.initial_velocity_m_s;
    if w == 0.0 {
        return Err(InertialCapillaryError::Invalid("stopping velocity is unrepresentable").into());
    }
    let mut lower = 0.0;
    let mut upper = requested;
    let mut time = requested;
    loop {
        let derivative = (m.mass0 * m.input.initial_velocity_m_s / time) / time
            + 0.5 * m.bulk_slope * w * w
            + 0.5 * m.area * m.head_slope * w
            + m.area * (pressure.derivative_wet_length_pa_m * w + pressure.derivative_time_pa_s);
        if !derivative.is_finite() || b.sign * derivative <= 0.0 {
            return Err(InertialCapillaryError::StepTooLarge(
                "cannot continue a stable stopping-time residual",
            )
            .into());
        }
        if residual.abs() / m.area <= m.controls.pressure_absolute_tolerance_pa
            && (residual == 0.0
                || (residual / derivative).abs() <= m.controls.time_absolute_tolerance_s)
        {
            let (x, _) = m.state(time, w)?;
            let endpoint_force =
                m.young_drive(x, pressure.endpoint_pressure_pa) - b.pin_force / m.area;
            if b.sign * endpoint_force >= -m.controls.pressure_absolute_tolerance_pa {
                return Err(InertialCapillaryError::StepTooLarge(
                    "midpoint stopping would replace asymptotic coasting or unopposed motion",
                )
                .into());
            }
            return finish(
                m,
                b,
                time,
                w,
                pressure,
                InertialCapillaryEnd::Stopped,
                residual,
                budget.used,
                cancelled,
            );
        }
        if b.sign * residual >= 0.0 {
            upper = time;
        } else {
            lower = time;
        }
        let newton = time - residual / derivative;
        let next = if newton.is_finite() && newton > lower && newton < upper {
            newton
        } else {
            lower + 0.5 * (upper - lower)
        };
        if next == time || next == lower || next == upper || next <= 0.0 {
            return Err(InertialCapillaryError::NotConverged {
                iterations: budget.used,
                residual_pa: residual / m.area,
            }
            .into());
        }
        let (x, _) = m.state(next, w)?;
        pressure = pressure_at(x, next, callback, cancelled, budget)?;
        let (r, _) = m.residual(b, next, w, &pressure)?;
        residual = r;
        budget.residual_pa = residual / m.area;
        time = next;
    }
}
fn finish<T, E>(
    m: Model,
    b: Branch,
    h: f64,
    w: f64,
    pressure: Pressure<T>,
    end: InertialCapillaryEnd,
    residual_force: f64,
    iterations: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<InertialPressureStep<T>, InertialSolveError<E>> {
    let (x, mut u1) = m.state(h, w)?;
    if end == InertialCapillaryEnd::Stopped {
        u1 = 0.0;
    }
    let a = m.input.meniscus.meniscus.meniscus;
    let displacement = h * w;
    let volume = m.area * displacement;
    let xbar = m.x0() + 0.5 * displacement;
    let mass1 = m.mass0 + m.mass_slope * displacement;
    let massbar = m.mass0 + 0.5 * m.mass_slope * displacement;
    let u0 = m.input.initial_velocity_m_s;
    let kinetic0 = 0.5 * m.mass0 * u0 * u0;
    let kinetic1 = 0.5 * mass1 * u1 * u1;
    let delta_kinetic =
        massbar * w * (u1 - u0) + 0.25 * m.mass_slope * (u1 * u1 + u0 * u0) * displacement;
    let inlet_work = a.applied_pressure_pa * volume;
    let downstream_work = pressure.downstream_pressure_pa * volume;
    let pressure_work = (a.applied_pressure_pa - pressure.downstream_pressure_pa) * volume;
    let wall = -m.young_pressure * volume;
    let gravity = m.head_slope * xbar * volume;
    let outlet = m.input.meniscus.displaced_density_kg_m3
        * m.input.meniscus.opposing_gravity_m_s2
        * m.length()
        * volume;
    let liquid_loss = 8.0 * PI * a.dynamic_viscosity_pa_s * xbar * w * w * h;
    let displaced_loss = 8.0
        * PI
        * m.input.meniscus.meniscus.displaced_dynamic_viscosity_pa_s
        * (m.length() - xbar)
        * w
        * w
        * h;
    let line_loss = m.line_resistance * w * w * h;
    let pinning_loss = b.pin_force * displacement;
    let reservoir_loss = b.reservoir * w.abs() * w * w * h;
    let loss = liquid_loss + displaced_loss + line_loss + pinning_loss + reservoir_loss;
    let residual = pressure_work - wall - gravity - outlet - delta_kinetic - loss;
    let scale = pressure_work.abs()
        + wall.abs()
        + gravity.abs()
        + outlet.abs()
        + delta_kinetic.abs()
        + loss;
    let tolerance =
        m.controls.energy_absolute_tolerance_j + m.controls.energy_relative_tolerance * scale;
    let angle0 = m.angle(b, u0)?;
    let angle1 = m.angle(b, u1)?;
    if ![
        mass1,
        volume,
        kinetic0,
        kinetic1,
        delta_kinetic,
        inlet_work,
        downstream_work,
        pressure_work,
        wall,
        gravity,
        outlet,
        liquid_loss,
        displaced_loss,
        line_loss,
        pinning_loss,
        reservoir_loss,
        loss,
        residual,
        tolerance,
        residual_force / m.area,
    ]
    .into_iter()
    .all(f64::is_finite)
        || mass1 <= 0.0
        || [
            liquid_loss,
            displaced_loss,
            line_loss,
            pinning_loss,
            reservoir_loss,
        ]
        .into_iter()
        .any(|v| v < 0.0)
        || (displacement != 0.0 && volume == 0.0)
        || (u0 != 0.0 && kinetic0 == 0.0)
        || (u1 != 0.0 && kinetic1 == 0.0)
    {
        return Err(InertialCapillaryError::Invalid(
            "state or mechanical accounting is unrepresentable",
        )
        .into());
    }
    if residual.abs() > tolerance {
        return Err(InertialCapillaryError::EnergyBalance {
            residual_j: residual,
            tolerance_j: tolerance,
        }
        .into());
    }
    poll(cancelled)?;
    Ok(InertialPressureStep {
        step: InertialCapillaryStep {
            initial_wet_length_m: m.x0(),
            wet_length_m: x,
            initial_velocity_m_s: u0,
            velocity_m_s: u1,
            elapsed_time_s: h,
            end,
            transferred_volume_m3: volume,
            initial_effective_mass_kg: m.mass0,
            effective_mass_kg: mass1,
            initial_kinetic_energy_j: kinetic0,
            kinetic_energy_j: kinetic1,
            kinetic_energy_change_j: delta_kinetic,
            inlet_pressure_work_j: inlet_work,
            downstream_pressure_work_j: downstream_work,
            pressure_work_j: pressure_work,
            wall_surface_energy_change_j: wall,
            gravitational_storage_change_j: gravity,
            outlet_potential_energy_j: outlet,
            liquid_viscous_dissipation_j: liquid_loss,
            displaced_viscous_dissipation_j: displaced_loss,
            contact_line_viscous_dissipation_j: line_loss,
            pinning_dissipation_j: pinning_loss,
            reservoir_dissipation_j: reservoir_loss,
            total_dissipation_j: loss,
            motion: b.motion,
            initial_apparent_angle_radians: angle0,
            final_apparent_angle_radians: angle1,
            pressure_residual_pa: residual_force / m.area,
            energy_balance_residual_j: residual,
            iterations,
        },
        pressure,
    })
}

#[cfg(test)]
#[path = "inertia_tests.rs"]
mod tests;
