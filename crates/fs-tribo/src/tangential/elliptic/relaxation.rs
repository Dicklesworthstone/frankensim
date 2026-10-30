//! Continuous associated-Coulomb relaxation under held load and body mobility.
//!
//! This is the exact separable velocity flow with numerical scalar quadrature,
//! not a sequence of midpoint impulses. In material axes, write the kinetic
//! capacity as `a*d`, `max(d)=1`, and set `u=d*v`, `A=d*W*d`. The time change
//! `d tau/dt=a/|u|` gives `du/dtau=-A*u`. Diagonalizing this symmetric matrix
//! reduces physical time, net displacement and unsigned travel to integrals on
//! a finite interval. The endpoint `z=0` is the actual finite stopping event.
//! Impulse and traction-slip work have closed forms; an average force times
//! net displacement is generally NOT the work of this curved trajectory.
//!
//! Geometry, inertia, load, coefficients and temperature stay fixed. There is
//! no external tangential drive, pose update, normal constraint or source query.
//! Quadrature errors are adaptive estimates, not interval certificates.

use super::*;

#[derive(Clone, Copy, Debug)]
pub struct EllipticRelaxationInput {
    pub velocity_m_s: [f64; 2],
    /// Full symmetric positive-definite tangent mobility, in the current frame.
    pub mobility_per_kg: [[f64; 2]; 2],
    pub normal_n: f64,
    pub dt_s: f64,
    pub maximum_force_n: f64,
    pub contact_temperature_k: Option<f64>,
    pub nominal_pressure_pa: Option<f64>,
}

#[derive(Clone, Copy, Debug)]
pub struct EllipticRelaxationControls {
    pub time_tolerance_s: f64,
    pub travel_tolerance_m: f64,
    pub maximum_evaluations: usize,
    pub maximum_root_iterations: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct EllipticRelaxationReport {
    pub state: EllipticCoulombState,
    pub velocity_after_m_s: [f64; 2],
    /// Positive traction impulse; the force on body A has the opposite sign.
    pub traction_impulse_ns: [f64; 2],
    pub slip_increment_m: [f64; 2],
    pub distance_m: f64,
    pub dissipated_j: f64,
    /// End of the kinetic phase within this interval. A frictionless direction
    /// can continue moving after this event when one capacity is exactly zero.
    pub kinetic_phase_end_s: Option<f64>,
    pub time_residual_s: f64,
    /// Estimates for [elapsed time, displacement axis 0, axis 1, unsigned path].
    pub quadrature_error: [f64; 4],
    pub evaluations: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum EllipticRelaxationError {
    Constitutive(TriboError),
    Invalid(&'static str),
    Cancelled,
    EvaluationBudget { evaluations: usize },
    UnresolvedQuadrature,
    UnresolvedTime,
}
impl From<TriboError> for EllipticRelaxationError {
    fn from(e: TriboError) -> Self {
        Self::Constitutive(e)
    }
}
impl core::fmt::Display for EllipticRelaxationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "directional relaxation: {self:?}")
    }
}
impl std::error::Error for EllipticRelaxationError {}
type Error = EllipticRelaxationError;

struct Work<'a, C> {
    controls: EllipticRelaxationControls,
    evaluations: usize,
    cancelled: &'a mut C,
}
impl<C: FnMut() -> bool> Work<'_, C> {
    fn poll(&mut self) -> Result<(), Error> {
        if (self.cancelled)() {
            Err(Error::Cancelled)
        } else {
            Ok(())
        }
    }
    fn evaluate(&mut self, f: &impl Fn(f64) -> [f64; 4], x: f64) -> Result<[f64; 4], Error> {
        self.poll()?;
        if self.evaluations >= self.controls.maximum_evaluations {
            return Err(Error::EvaluationBudget {
                evaluations: self.evaluations,
            });
        }
        self.evaluations += 1;
        let value = f(x);
        if value.iter().any(|x| !x.is_finite()) {
            return Err(Error::Invalid("unrepresentable relaxation integrand"));
        }
        Ok(value)
    }
    fn integrate(
        &mut self,
        f: &impl Fn(f64) -> [f64; 4],
        a: f64,
        b: f64,
    ) -> Result<([f64; 4], [f64; 4]), Error> {
        if a == b {
            return Ok(([0.0; 4], [0.0; 4]));
        }
        let values = [
            self.evaluate(f, a)?,
            self.evaluate(f, 0.5 * a + 0.5 * b)?,
            self.evaluate(f, b)?,
        ];
        let tolerance = [
            self.controls.time_tolerance_s / 8.0,
            self.controls.travel_tolerance_m / 8.0,
            self.controls.travel_tolerance_m / 8.0,
            self.controls.travel_tolerance_m / 8.0,
        ];
        self.panel(f, a, b, values, tolerance, 0)
    }
    fn panel(
        &mut self,
        f: &impl Fn(f64) -> [f64; 4],
        a: f64,
        b: f64,
        values: [[f64; 4]; 3],
        tolerance: [f64; 4],
        depth: usize,
    ) -> Result<([f64; 4], [f64; 4]), Error> {
        let mid = 0.5 * a + 0.5 * b;
        let left = self.evaluate(f, 0.5 * a + 0.5 * mid)?;
        let right = self.evaluate(f, 0.5 * mid + 0.5 * b)?;
        let coarse = simpson(a, b, values);
        let fine_left = simpson(a, mid, [values[0], left, values[1]]);
        let fine_right = simpson(mid, b, [values[1], right, values[2]]);
        let fine: [f64; 4] = core::array::from_fn(|i| fine_left[i] + fine_right[i]);
        let error: [f64; 4] = core::array::from_fn(|i| (fine[i] - coarse[i]).abs() / 15.0);
        if error.iter().zip(tolerance).all(|(e, t)| *e <= t) {
            return Ok((
                core::array::from_fn(|i| fine[i] + (fine[i] - coarse[i]) / 15.0),
                error,
            ));
        }
        if depth >= 48 || a == mid || b == mid {
            return Err(Error::UnresolvedQuadrature);
        }
        let half = tolerance.map(|t| t * 0.5);
        let (l, le) = self.panel(f, a, mid, [values[0], left, values[1]], half, depth + 1)?;
        let (r, re) = self.panel(f, mid, b, [values[1], right, values[2]], half, depth + 1)?;
        Ok((
            core::array::from_fn(|i| l[i] + r[i]),
            core::array::from_fn(|i| le[i] + re[i]),
        ))
    }
    fn clock_root(
        &mut self,
        f: &impl Fn(f64) -> [f64; 4],
        target: f64,
        mut x: f64,
        increasing: bool,
    ) -> Result<(f64, [f64; 4], [f64; 4], f64), Error> {
        let mut lower = 0.0;
        let mut upper = 1.0;
        for _ in 0..self.controls.maximum_root_iterations {
            let (value, error) = if increasing {
                self.integrate(f, 0.0, x)?
            } else {
                self.integrate(f, x, 1.0)?
            };
            let residual = value[0] - target;
            if residual.abs() + error[0] <= self.controls.time_tolerance_s {
                return Ok((x, value, error, residual));
            }
            if (residual > 0.0) == increasing {
                upper = x;
            } else {
                lower = x;
            }
            let slope = self.evaluate(f, x)?[0];
            let newton = x - residual / if increasing { slope } else { -slope };
            let next = if newton.is_finite() && newton > lower && newton < upper {
                newton
            } else {
                0.5 * lower + 0.5 * upper
            };
            if next == x || next == lower || next == upper {
                return Err(Error::UnresolvedTime);
            }
            x = next;
        }
        Err(Error::UnresolvedTime)
    }
}
fn simpson(a: f64, b: f64, y: [[f64; 4]; 3]) -> [f64; 4] {
    core::array::from_fn(|i| ((y[0][i] + y[2][i]) / 6.0 + 2.0 * y[1][i] / 3.0) * (b - a))
}
fn checked(x: f64) -> Result<f64, Error> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(Error::Invalid("unrepresentable relaxation response"))
    }
}
fn length(x: [f64; 2]) -> f64 {
    det::hypot(x[0], x[1])
}
fn power(z: f64, p: f64) -> f64 {
    if z == 0.0 {
        0.0
    } else {
        det::exp(p * det::ln(z))
    }
}
fn decay(s: f64, p: f64) -> f64 {
    if s == 0.0 { 1.0 } else { det::exp(-p * s) }
}
fn lost(s: f64, p: f64) -> f64 {
    if s == 0.0 {
        0.0
    } else if s == f64::INFINITY {
        1.0
    } else {
        -det::expm1(-p * s)
    }
}

impl EllipticCoulomb {
    /// Evolve the original native law through curved sliding and finite arrest.
    /// The returned checkpoint retains its static coefficients and texture.
    /// All work is bounded, and no accepted checkpoint is mutated on failure.
    pub fn relax_velocity(
        &self,
        state: &EllipticCoulombState,
        input: &EllipticRelaxationInput,
        controls: EllipticRelaxationControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<EllipticRelaxationReport, Error> {
        self.check(state)?;
        positive_finite(input.dt_s, "relaxation_dt_s")?;
        nonnegative_finite(input.normal_n, "relaxation_normal_n")?;
        nonnegative_finite(input.maximum_force_n, "relaxation_maximum_force_n")?;
        for x in input.velocity_m_s {
            finite(x, "relaxation_velocity_m_s")?;
        }
        if let Some(t) = input.contact_temperature_k {
            positive_finite(t, "contact_temperature_k")?;
        }
        if let Some(p) = input.nominal_pressure_pa {
            nonnegative_finite(p, "nominal_pressure_pa")?;
        }
        if !controls.time_tolerance_s.is_finite()
            || controls.time_tolerance_s <= 0.0
            || !controls.travel_tolerance_m.is_finite()
            || controls.travel_tolerance_m <= 0.0
            || controls.maximum_evaluations == 0
            || controls.maximum_root_iterations == 0
        {
            return Err(Error::Invalid(
                "positive finite relaxation tolerances and work limits required",
            ));
        }
        let mobility = SymmetricTangentialCompliance::new(input.mobility_per_kg)?;
        let scale = input.mobility_per_kg[0][0].max(input.mobility_per_kg[1][1]);
        let (w, determinant) = mobility.material(state.axis)?;
        if determinant <= 0.0 || scale <= 0.0 {
            return Err(Error::Invalid(
                "relaxation requires positive-definite body mobility",
            ));
        }
        let v = to_principal(state.axis, input.velocity_m_s)?;
        let caps = capacities(self.kinetic_mu, input.normal_n)?;
        let maximum = caps[0].max(caps[1]);
        if maximum > input.maximum_force_n {
            return Err(Error::Invalid(
                "kinetic capacity exceeds admitted force budget",
            ));
        }
        let mut work = Work {
            controls,
            evaluations: 0,
            cancelled,
        };
        work.poll()?;
        let mut report = EllipticRelaxationReport {
            state: state.clone(),
            velocity_after_m_s: input.velocity_m_s,
            traction_impulse_ns: [0.0; 2],
            slip_increment_m: [0.0; 2],
            distance_m: 0.0,
            dissipated_j: 0.0,
            kinetic_phase_end_s: None,
            time_residual_s: 0.0,
            quadrature_error: [0.0; 4],
            evaluations: 0,
        };
        if maximum == 0.0 || v == [0.0; 2] {
            report.slip_increment_m = input.velocity_m_s.map(|v| v * input.dt_s);
            report.distance_m = length(input.velocity_m_s) * input.dt_s;
            report.kinetic_phase_end_s = Some(0.0);
        } else if caps[0] == 0.0 || caps[1] == 0.0 {
            let j = usize::from(caps[0] == 0.0);
            let k = 1 - j;
            let diagonal = w[2 * j] * scale;
            let stop = checked(v[j].abs() / diagonal / caps[j])?;
            let active = input.dt_s.min(stop);
            let impulse = checked(caps[j] * active * v[j].signum())?;
            let mut end = v;
            end[j] = if input.dt_s >= stop {
                0.0
            } else {
                v[j] - diagonal * impulse
            };
            end[k] -= scale * w[1] * impulse;
            let travel =
                [0, 1].map(|i| 0.5 * (v[i] + end[i]) * active + end[i] * (input.dt_s - active));
            let mut impulses = [0.0; 2];
            impulses[j] = impulse;
            report.velocity_after_m_s = from_principal(state.axis, end)?;
            report.traction_impulse_ns = from_principal(state.axis, impulses)?;
            report.slip_increment_m = from_principal(state.axis, travel)?;
            report.dissipated_j = checked(0.5 * impulse * (v[j] + end[j]))?;
            let (path, error) = if active == 0.0 {
                ([0.0; 4], [0.0; 4])
            } else {
                work.integrate(
                    &|x| {
                        [
                            0.0,
                            0.0,
                            0.0,
                            active
                                * length([
                                    v[0] * (1.0 - x) + end[0] * x,
                                    v[1] * (1.0 - x) + end[1] * x,
                                ]),
                        ]
                    },
                    0.0,
                    1.0,
                )?
            };
            report.distance_m = path[3] + length(end) * (input.dt_s - active);
            report.quadrature_error = error;
            if input.dt_s >= stop {
                report.kinetic_phase_end_s = Some(stop);
            }
        } else {
            let d = caps.map(|x| x / maximum);
            if d.contains(&0.0) {
                return Err(Error::Invalid("unrepresentable kinetic capacity ratio"));
            }
            let a = [d[0] * d[0] * w[0], d[0] * d[1] * w[1], d[1] * d[1] * w[2]];
            let (lambda, axis) =
                super::compliance::spectrum(a, determinant * d[0] * d[0] * d[1] * d[1])?;
            let lambda = lambda.map(|x| x * scale);
            if lambda.iter().any(|x| !x.is_finite() || *x <= 0.0) {
                return Err(Error::Invalid("unrepresentable relaxation spectrum"));
            }
            let u = [d[0] * v[0], d[1] * v[1]];
            if (0..2).any(|i| v[i] != 0.0 && u[i] == 0.0) {
                return Err(Error::Invalid("unrepresentable weighted slip velocity"));
            }
            let c = to_principal(axis, u)?;
            let ratio = checked(lambda[0] / lambda[1])?;
            let denominator = checked(maximum * lambda[1])?;
            if denominator == 0.0 {
                return Err(Error::Invalid("unrepresentable relaxation timescale"));
            }
            if ratio == 1.0 || c[0] == 0.0 || c[1] == 0.0 {
                let rate = if c[1] == 0.0 { lambda[0] } else { lambda[1] };
                let stop = checked(length(u) / maximum / rate)?;
                let active = input.dt_s.min(stop);
                let fraction = if input.dt_s >= stop {
                    1.0
                } else {
                    active / stop
                };
                let end = v.map(|x| x * (1.0 - fraction));
                let force = [
                    maximum * d[0] * u[0] / length(u),
                    maximum * d[1] * u[1] / length(u),
                ];
                let impulse = force.map(|f| f * active);
                report.velocity_after_m_s = from_principal(state.axis, end)?;
                report.traction_impulse_ns = from_principal(state.axis, impulse)?;
                report.slip_increment_m =
                    from_principal(state.axis, v.map(|x| x * (1.0 - 0.5 * fraction) * active))?;
                report.distance_m = length(report.slip_increment_m);
                report.dissipated_j =
                    checked(0.5 * (impulse[0] * (v[0] + end[0]) + impulse[1] * (v[1] + end[1])))?;
                if input.dt_s >= stop {
                    report.kinetic_phase_end_s = Some(stop);
                }
            } else {
                let integrand = |x: f64| {
                    // z=x^8 regularizes the finite-stop endpoint. In z the
                    // clock contains z^(2*(ratio-1)), whose derivatives are
                    // singular for ordinary nearly equal eigenvalues. The
                    // Jacobian 8*x^7 makes the clock and travel integrands
                    // sufficiently smooth for adaptive Simpson without changing
                    // the physical curve, tolerance or recursion/work budget.
                    let x2 = x * x;
                    let x4 = x2 * x2;
                    let z = x4 * x4;
                    let jacobian = 8.0 * x4 * x2 * x;
                    let h = length([c[0] * power(z, ratio - 1.0), c[1]]);
                    let principal = from_eigen(axis, [c[0] * power(z, ratio), c[1] * z]);
                    let velocity = [principal[0] / d[0], principal[1] / d[1]];
                    let clock = (h / denominator) * jacobian;
                    [
                        clock,
                        clock * velocity[0],
                        clock * velocity[1],
                        clock * length(velocity),
                    ]
                };
                let fast_time = checked(length(u) / maximum / lambda[0])?;
                let (s, integral, error) = if input.dt_s < 0.25 * fast_time {
                    // Short steps never need the full stopping-time integral.
                    // |u(tau)| >= |u0| exp(-lambda_fast*tau) and
                    // 1-exp(-y) >= y/(1+y) provide a conservative bracket
                    // at twice the requested time. Integrating in s=-log(z)
                    // preserves tiny durations that cannot resolve 1-z.
                    let fraction = 2.0 * input.dt_s / fast_time;
                    let upper_s = checked((fraction / ratio) / (1.0 - fraction))?;
                    if upper_s == 0.0 {
                        return Err(Error::UnresolvedTime);
                    }
                    let clock_scale = checked((2.0 * input.dt_s / length(u)) / (1.0 - fraction))?;
                    let short = |x: f64| {
                        let eigen = [
                            c[0] * decay(upper_s * x, ratio),
                            c[1] * decay(upper_s * x, 1.0),
                        ];
                        let principal = from_eigen(axis, eigen);
                        let velocity = [principal[0] / d[0], principal[1] / d[1]];
                        let clock = clock_scale * length(eigen);
                        [
                            clock,
                            clock * velocity[0],
                            clock * velocity[1],
                            clock * length(velocity),
                        ]
                    };
                    let (x, value, error, residual) =
                        work.clock_root(&short, input.dt_s, 0.5, true)?;
                    report.time_residual_s = residual;
                    (upper_s * x, value, error)
                } else {
                    let (whole, whole_error) = work.integrate(&integrand, 0.0, 1.0)?;
                    if input.dt_s > whole[0] + whole_error[0] {
                        report.kinetic_phase_end_s = Some(whole[0]);
                        (f64::INFINITY, whole, whole_error)
                    } else {
                        if input.dt_s >= whole[0] - whole_error[0] {
                            return Err(Error::UnresolvedTime);
                        }
                        let start = power((1.0 - input.dt_s / whole[0]).clamp(0.0, 1.0), 0.125);
                        let (x, value, error, residual) =
                            work.clock_root(&integrand, input.dt_s, start, false)?;
                        report.time_residual_s = residual;
                        (-8.0 * det::ln(x), value, error)
                    }
                };
                let remaining = [decay(s, ratio), decay(s, 1.0)];
                let ue = from_eigen(axis, [c[0] * remaining[0], c[1] * remaining[1]]);
                report.velocity_after_m_s =
                    from_principal(state.axis, [ue[0] / d[0], ue[1] / d[1]])?;
                let integrated_u = from_eigen(
                    axis,
                    [
                        c[0] * lost(s, ratio) / lambda[0],
                        c[1] * lost(s, 1.0) / lambda[1],
                    ],
                );
                report.traction_impulse_ns =
                    from_principal(state.axis, [d[0] * integrated_u[0], d[1] * integrated_u[1]])?;
                report.slip_increment_m = from_principal(state.axis, [integral[1], integral[2]])?;
                report.distance_m = integral[3];
                let weighted = [c[0] / det::sqrt(lambda[0]), c[1] / det::sqrt(lambda[1])];
                report.dissipated_j = checked(
                    0.5 * weighted[0] * weighted[0] * lost(s, 2.0 * ratio)
                        + 0.5 * weighted[1] * weighted[1] * lost(s, 2.0),
                )?;
                report.quadrature_error = [
                    error[0],
                    state.axis[0].abs() * error[1] + state.axis[1].abs() * error[2],
                    state.axis[1].abs() * error[1] + state.axis[0].abs() * error[2],
                    error[3],
                ];
            }
        }
        for value in report
            .velocity_after_m_s
            .into_iter()
            .chain(report.traction_impulse_ns)
            .chain(report.slip_increment_m)
            .chain([report.distance_m, report.dissipated_j])
        {
            checked(value)?;
        }
        if report.distance_m < 0.0 || report.dissipated_j < 0.0 {
            return Err(Error::Invalid("negative physical travel or dissipation"));
        }
        report.state.contact_temperature_k = input.contact_temperature_k;
        for i in 0..2 {
            report.state.slip_m[i] = checked(state.slip_m[i] + report.slip_increment_m[i])?;
        }
        report.evaluations = work.evaluations;
        work.poll()?;
        Ok(report)
    }
}
fn from_eigen(axis: [f64; 2], x: [f64; 2]) -> [f64; 2] {
    [
        axis[0] * x[0] - axis[1] * x[1],
        axis[1] * x[0] + axis[0] * x[1],
    ]
}
