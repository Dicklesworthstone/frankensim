//! Adaptive P1 moments of the existing Green field, not another acoustic solve.
//!
//! Keep original-panel barycentric coordinates through every subdivision. A
//! P1 trace is evaluated at the quadrature point, never at the panel centre.
//! The same three single-layer moments sum to its P0 moment by partition of
//! unity. Kernel signs, analytic spatial derivatives, geometric separation,
//! phase resolution and quadrature rules all belong to the parent owner.
use super::*;

type Moments = [Integral; 3];

/// Generic open-screen double-layer evaluation. Geometry admission belongs to
/// the screen owner; this shares the existing adaptive P1 moments and kernels.
pub(crate) fn open_pressure(
    k: f64,
    triangles: &[Triangle],
    normals: &[Point],
    jumps: &[[C64; 3]],
    points: &[Point],
    separation: f64,
    options: Options,
    mut cancelled: impl FnMut() -> bool,
) -> Result<Vec<C64>, HelmholtzError> {
    options.validate()?;
    if !k.is_finite()
        || k <= 0.0
        || !separation.is_finite()
        || separation <= 0.0
        || normals.len() != triangles.len()
        || jumps.len() != triangles.len()
        || points.iter().flatten().any(|v| !v.is_finite())
        || jumps
            .iter()
            .flatten()
            .any(|v| !v.re.is_finite() || !v.im.is_finite())
    {
        return Err(bad("invalid open-screen receiver request"));
    }
    if points
        .len()
        .checked_mul(triangles.len())
        .and_then(|v| v.checked_mul(80))
        .is_none_or(|v| v > options.maximum_kernel_evaluations)
    {
        return Err(bad("near-field kernel-evaluation budget exhausted"));
    }
    let mut work = Work {
        options,
        evaluations: 0,
    };
    let mut out = Vec::with_capacity(points.len());
    for &x in points {
        let mut p = C64::ZERO;
        for ((&t, &normal), jump) in triangles.iter().zip(normals).zip(jumps) {
            if distance(x, t) <= separation {
                return Err(bad("open-screen receiver touches surface"));
            }
            let moments = integrate_linear(
                k,
                x,
                t,
                normal,
                [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
                0,
                &mut work,
                false,
                &mut cancelled,
            )?;
            for a in 0..3 {
                p = p + moments[a].d * jump[a];
            }
        }
        out.push(p);
    }
    if cancelled() {
        return Err(HelmholtzError::Cancelled);
    }
    Ok(out)
}

fn rule_linear(
    k: f64,
    x: Point,
    t: Triangle,
    normal: Point,
    basis: Triangle,
    nodes: &[f64],
    weights: &[f64],
    work: &mut Work,
    with_velocity: bool,
) -> Result<Moments, HelmholtzError> {
    work.evaluations = work
        .evaluations
        .checked_add(nodes.len() * nodes.len())
        .filter(|n| *n <= work.options.maximum_kernel_evaluations)
        .ok_or_else(|| bad("near-field kernel-evaluation budget exhausted"))?;
    let u = sub(t[1], t[0]);
    let v = sub(t[2], t[1]);
    let area2 = norm(cross(u, sub(t[2], t[0])));
    let mut out = [Integral::default(); 3];
    for (&a, &wa) in nodes.iter().zip(weights) {
        for (&b, &wb) in nodes.iter().zip(weights) {
            let y = std::array::from_fn(|c| t[0][c] + a * u[c] + a * b * v[c]);
            let sample = kernel(k, sub(x, y), normal)?;
            let local = [1. - a, a * (1. - b), a * b];
            let jac = wa * wb * area2 * a;
            for node in 0..3 {
                let shape = (0..3).map(|c| local[c] * basis[c][node]).sum::<f64>();
                out[node].accumulate(&sample, jac * shape, with_velocity)?;
            }
        }
    }
    if out.iter().any(|c| {
        ![c.s.re, c.s.im, c.d.re, c.d.im, c.abs_s, c.abs_d]
            .iter()
            .all(|v| v.is_finite())
    }) {
        return Err(bad("nonfinite linear near-field quadrature"));
    }
    Ok(out)
}

#[allow(clippy::too_many_arguments)]
fn integrate_linear(
    k: f64,
    x: Point,
    t: Triangle,
    normal: Point,
    basis: Triangle,
    depth: usize,
    work: &mut Work,
    with_velocity: bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Moments, HelmholtzError> {
    if cancelled() {
        return Err(HelmholtzError::Cancelled);
    }
    let coarse = rule_linear(k, x, t, normal, basis, &X4, &W4, work, with_velocity)?;
    let mut fine = rule_linear(k, x, t, normal, basis, &X8, &W8, work, with_velocity)?;
    let size = [
        norm(sub(t[1], t[0])),
        norm(sub(t[2], t[1])),
        norm(sub(t[0], t[2])),
    ]
    .into_iter()
    .fold(0.0_f64, f64::max);
    let mut accepted = size <= 2. * distance(x, t) && k * size <= 3.;
    let tol = work.options.relative_tolerance;
    for (f, c) in fine.iter_mut().zip(&coarse) {
        f.error_s = (f.s - c.s).abs();
        f.error_d = (f.d - c.d).abs();
        let gradient_ok = !with_velocity || f.gradient.accept(&c.gradient, tol);
        accepted &= f.error_s <= tol * f.abs_s && f.error_d <= tol * f.abs_d && gradient_ok;
    }
    if accepted {
        return Ok(fine);
    }
    if depth == work.options.maximum_depth {
        return Err(bad("linear near-field triangle refinement depth exhausted"));
    }
    let mut sum = [Integral::default(); 3];
    for (child, shape) in children(t).into_iter().zip(children(basis)) {
        let values = integrate_linear(
            k,
            x,
            child,
            normal,
            shape,
            depth + 1,
            work,
            with_velocity,
            cancelled,
        )?;
        for (total, value) in sum.iter_mut().zip(values) {
            total.add(value);
        }
    }
    Ok(sum)
}

/// Internal rows; the reciprocal solver is the public source-identity boundary.
pub(crate) struct Rows {
    rows: Vec<Vec<Moments>>,
    omega_rho: f64,
    pub(crate) with_velocity: bool,
    pub(crate) kernel_evaluations: usize,
}
impl Geometry<'_> {
    pub(crate) fn prepare_linear(
        &self,
        k: f64,
        medium: Medium,
        options: Options,
        with_velocity: bool,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Rows, HelmholtzError> {
        if cancelled() {
            return Err(HelmholtzError::Cancelled);
        }
        options.validate()?;
        let omega_rho = k * medium.sound_speed * medium.density;
        if !k.is_finite()
            || k <= 0.
            || !medium.density.is_finite()
            || medium.density <= 0.
            || !medium.sound_speed.is_finite()
            || medium.sound_speed <= 0.
            || !omega_rho.is_finite()
            || omega_rho <= 0.
        {
            return Err(bad(
                "linear near-field evaluation needs positive finite k and medium",
            ));
        }
        let triangles = self
            .surface
            .triangles()
            .ok_or_else(|| bad("missing retained receiver triangles"))?;
        // Every panel consumes at least 16+64 samples. Preflight before retaining
        // rows, including requests whose geometry could never fit the work cap.
        if self
            .points
            .len()
            .checked_mul(triangles.len())
            .and_then(|n| n.checked_mul(80))
            .is_none_or(|n| n > options.maximum_kernel_evaluations)
        {
            return Err(bad("near-field kernel-evaluation budget exhausted"));
        }
        let mut work = Work {
            options,
            evaluations: 0,
        };
        let mut rows = Vec::with_capacity(self.points.len());
        for &point in &self.points {
            let mut row = Vec::new();
            row.try_reserve_exact(triangles.len()).map_err(|_| {
                HelmholtzError::AllocationFailed {
                    what: "linear near-field moments",
                }
            })?;
            for (i, &t) in triangles.iter().enumerate() {
                row.push(integrate_linear(
                    k,
                    point,
                    t,
                    self.surface.normals()[i],
                    [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
                    0,
                    &mut work,
                    with_velocity,
                    &mut cancelled,
                )?);
            }
            rows.push(row);
        }
        if cancelled() {
            return Err(HelmholtzError::Cancelled);
        }
        Ok(Rows {
            rows,
            omega_rho,
            with_velocity,
            kernel_evaluations: work.evaluations,
        })
    }
}
impl Rows {
    /// Evaluate D[p_P1]+S[lambda_P0] and optionally its analytic spatial gradient.
    /// The source callback addresses ORIGINAL panels/vertices, not leaf triangles.
    pub(crate) fn evaluate(
        &self,
        source: impl Fn(usize) -> ([C64; 3], C64),
        velocity: bool,
    ) -> Result<VectorObservation, HelmholtzError> {
        if velocity && !self.with_velocity {
            return Err(bad(
                "linear particle velocity requires derivative preparation",
            ));
        }
        let mut pressure = Vec::with_capacity(self.rows.len());
        let mut errors = Vec::with_capacity(self.rows.len());
        let mut vectors = Vec::with_capacity(self.rows.len());
        let mut vector_errors = Vec::with_capacity(self.rows.len());
        for row in &self.rows {
            let (mut re, mut im, mut cr, mut ci, mut error) =
                ([0.; 4], [0.; 4], [0.; 4], [0.; 4], [0.; 4]);
            for (i, moments) in row.iter().enumerate() {
                let (p, lambda) = source(i);
                if p.iter()
                    .chain(std::iter::once(&lambda))
                    .any(|v| !v.re.is_finite() || !v.im.is_finite())
                {
                    return Err(bad("nonfinite linear source coefficients"));
                }
                for (cell, p) in moments.iter().zip(p) {
                    let term = cell.d * p + cell.s * lambda;
                    compensated(&mut re[0], &mut cr[0], term.re);
                    compensated(&mut im[0], &mut ci[0], term.im);
                    error[0] += p.abs() * cell.error_d + lambda.abs() * cell.error_s;
                    if velocity {
                        for c in 0..3 {
                            let term = cell.gradient.d[c] * p + cell.gradient.s[c] * lambda;
                            compensated(&mut re[c + 1], &mut cr[c + 1], term.re);
                            compensated(&mut im[c + 1], &mut ci[c + 1], term.im);
                            error[c + 1] += p.abs() * cell.gradient.error_d[c]
                                + lambda.abs() * cell.gradient.error_s[c];
                        }
                    }
                }
            }
            if re.iter().chain(&im).chain(&error).any(|v| !v.is_finite()) {
                return Err(bad("nonfinite linear pressure, gradient or discrepancy"));
            }
            pressure.push(C64::new(re[0], im[0]));
            errors.push(error[0]);
            if velocity {
                // Euler momentum for exp(-i omega t), not normal surface motion.
                let inverse = C64::new(0., -1. / self.omega_rho);
                let values = std::array::from_fn(|c| C64::new(re[c + 1], im[c + 1]) * inverse);
                let estimates = std::array::from_fn(|c| error[c + 1] / self.omega_rho);
                if values
                    .iter()
                    .any(|v: &C64| !v.re.is_finite() || !v.im.is_finite())
                    || estimates.iter().any(|v: &f64| !v.is_finite())
                {
                    return Err(bad("linear particle velocity overflow"));
                }
                vectors.push(values);
                vector_errors.push(estimates);
            }
        }
        Ok(VectorObservation {
            scalar: Observation {
                pressure,
                quadrature_error_estimate_pa: errors,
            },
            particle_velocity_m_s: vectors,
            quadrature_error_estimate_m_s: vector_errors,
        })
    }
}
