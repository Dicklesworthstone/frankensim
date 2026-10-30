//! Damped continuation of the parent Burton–Miller operator for CQ references.
//! With the Laplace transform exp(-s t), k=i*s/c and q=-rho*s*v. The parent
//! point kernels, equivalent-disc self terms and static row closures are shared;
//! this module introduces no boundary discretization or passivity repair.
use super::*;

/// Two-sided complex-frequency continuation and analytic wavenumber derivatives
/// for exterior resonance pencils. Unlike CQ, this explicitly admits Im(k)<0.
#[path = "helmholtz_complex.rs"]
pub mod complex;

#[derive(Clone, Copy)]
pub(crate) struct DampedWaveNumber {
    pub(crate) k: C64,
    pub(crate) phase_origin_m: f64,
}
impl WaveNumber for DampedWaveNumber {
    fn value(self) -> C64 {
        self.k
    }
    fn outgoing(self, r: f64) -> C64 {
        let distance = r - self.phase_origin_m;
        expik(self.k.re * distance).scale(det::exp(-self.k.im * distance))
    }
    fn radial_factor(self, r: f64) -> C64 {
        C64::new(-1.0 / r - self.k.im, self.k.re)
    }
    fn disc_self(self, area: f64) -> (C64, C64) {
        let a = det::sqrt(area / core::f64::consts::PI);
        let ik = C64::new(-self.k.im, self.k.re);
        let z = ik.scale(a);
        // exp(z)-1 and exp(z)-1-z must retain their small terms independently.
        let (em1, remainder) = if z.abs() < 0.25 {
            let mut term = z;
            let mut tail = C64::ZERO;
            for j in 2..=20 {
                term = (term * z).scale(1.0 / f64::from(j));
                tail = tail + term;
            }
            (z + tail, tail)
        } else {
            let em1 = expik(z.im).scale(det::exp(z.re)) - C64::ONE;
            (em1, em1 - z)
        };
        (em1 * ik.scale(2.0).recip(), remainder.scale(-0.5 / a))
    }
    fn cbie_real(self) -> f64 {
        // Only nonzero-alpha Burton–Miller is admitted below. Its assembly
        // never enters the parent's real-frequency triangle-CBIE branches.
        debug_assert_eq!(self.k.im, 0.0);
        self.k.re
    }
}

/// One boundary solution at an explicit Laplace point. Fields are immutable so
/// they cannot be cross-wired with a different geometry, medium or frequency.
pub struct LaplaceSolution {
    s: C64,
    medium: Medium,
    fingerprint: u64,
    pressure: Vec<C64>,
    velocity: Vec<C64>,
}
impl LaplaceSolution {
    /// Laplace frequency [1/s], Re s >= 0, excluding zero.
    pub fn frequency_per_s(&self) -> C64 {
        self.s
    }
    /// Constant-panel surface pressure [Pa].
    pub fn pressure_pa(&self) -> &[C64] {
        &self.pressure
    }
}

fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}
fn bad(what: &'static str) -> HelmholtzError {
    HelmholtzError::BadParameter { what }
}

/// Solve the SAME dense Burton–Miller system at k=i*s/c. A batch shares its
/// factorization. The parent dense-panel, resolution and 256-source caps apply
/// to |k|, including the decay length; complex frequencies do not bypass them.
/// # Errors
/// Nonfinite/left-half-plane s, invalid fields or medium, parent work/resolution
/// refusals, singular factorization, or a nonfinite pressure.
pub fn solve_radiation_laplace_batch(
    surface: &SpherePanels,
    s: C64,
    medium: Medium,
    fields: &[&[C64]],
) -> Result<Vec<LaplaceSolution>, HelmholtzError> {
    if !finite(s)
        || s.re < 0.0
        || s.abs() == 0.0
        || !medium.density.is_finite()
        || medium.density <= 0.0
        || !medium.sound_speed.is_finite()
        || medium.sound_speed <= 0.0
    {
        return Err(bad(
            "Laplace BEM requires nonzero Re(s)>=0 and positive finite medium",
        ));
    }
    let n = surface.centroids().len();
    surface
        .validate()
        .map_err(|_| bad("invalid Laplace source surface"))?;
    if fields.is_empty()
        || fields.len() > MAX_RADIATION_FIELDS_PER_BATCH
        || fields
            .iter()
            .any(|f| f.len() != n || f.iter().any(|&v| !finite(v)))
    {
        return Err(bad("invalid Laplace boundary-field batch"));
    }
    let k = C64::new(-s.im / medium.sound_speed, s.re / medium.sound_speed);
    admit_resolution(surface, k.abs())?;
    let wave = DampedWaveNumber {
        k,
        phase_origin_m: 0.0,
    };
    // Same Burton-Miller coupling sign as the real-frequency owner (`alpha_for`).
    let alpha = C64::new(0.0, 1.0) * k.recip();
    let (a, b) = assemble_dense_with_static(surface, wave, alpha, None);
    if a.iter().chain(&b).any(|&v| !finite(v)) {
        return Err(bad("nonfinite damped BEM matrix"));
    }
    let lu = lu_complex(&a, n).map_err(|_| HelmholtzError::Singular)?;
    let q_scale = s.scale(-medium.density);
    let mut solutions = Vec::with_capacity(fields.len());
    for &field in fields {
        let q: Vec<_> = field.iter().map(|&v| q_scale * v).collect();
        let mut pressure = vec![C64::ZERO; n];
        for i in 0..n {
            for j in 0..n {
                pressure[i] = pressure[i] + b[i * n + j] * q[j];
            }
        }
        lu.solve(&mut pressure);
        if pressure.iter().any(|&p| !finite(p)) {
            return Err(bad("nonfinite damped BEM solution"));
        }
        solutions.push(LaplaceSolution {
            s,
            medium,
            fingerprint: surface_fingerprint(surface),
            pressure,
            velocity: field.to_vec(),
        });
    }
    Ok(solutions)
}

/// Pressure [Pa] and its Cartesian gradient [Pa/m] at one finite receiver.
#[derive(Debug, Clone, Copy)]
pub struct PointField {
    pub pressure_pa: C64,
    pub gradient_pa_per_m: [C64; 3],
}

// This is the centroid Green representation already exposed by the parent.
// Both Fourier bakes and damped CQ use these same analytic normal derivatives.
// No near-field quadrature accuracy is claimed for a centroid receiver row.
pub(super) fn fields_from_traces(
    surface: &SpherePanels,
    pressure: &[C64],
    velocity: &[C64],
    wave: DampedWaveNumber,
    q_scale: C64,
    points: &[[f64; 3]],
) -> Result<Vec<PointField>, HelmholtzError> {
    let mut fields = Vec::with_capacity(points.len());
    for &point in points {
        let mut p = C64::ZERO;
        let mut gradient = [C64::ZERO; 3];
        for j in 0..pressure.len() {
            let y = surface.centroids()[j];
            let ny = surface.normals()[j];
            let area = surface.areas()[j];
            let q = q_scale * velocity[j];
            let (g, d, _, _) = kernels(wave, point, [0.0; 3], y, ny);
            p = p + (d * pressure[j] - g * q).scale(area);
            for c in 0..3 {
                let mut axis = [0.0; 3];
                axis[c] = 1.0;
                let (_, _, ds, dd) = kernels(wave, point, axis, y, ny);
                gradient[c] = gradient[c] + (dd * pressure[j] - ds * q).scale(area);
            }
        }
        if !finite(p) || gradient.iter().any(|&v| !finite(v)) {
            return Err(bad("nonfinite exterior pressure/gradient"));
        }
        fields.push(PointField {
            pressure_pa: p,
            gradient_pa_per_m: gradient,
        });
    }
    Ok(fields)
}

/// Existing centroid observation with analytic Cartesian derivatives, exp(-iwt).
/// Geometry/clearance and solution identity are checked by existing owners.
/// # Errors
/// Invalid/interior/near receiver, cross-wired solution or nonfinite field.
pub fn observe_frequency(
    surface: &SpherePanels,
    solution: &RadiationSolution,
    points: &[[f64; 3]],
    minimum_clearance_m: f64,
) -> Result<Vec<PointField>, HelmholtzError> {
    crate::near_field::Geometry::new(surface, points, minimum_clearance_m)?;
    exterior_pressure_at_points(surface, solution, solution.medium, &[])?;
    fields_from_traces(
        surface,
        &solution.pressure,
        &solution.velocity,
        DampedWaveNumber {
            k: C64::from_re(solution.k),
            phase_origin_m: 0.0,
        },
        C64::new(
            0.0,
            solution.k * solution.medium.sound_speed * solution.medium.density,
        ),
        points,
    )
}

/// Damped observation; optional known delay extraction is evaluated INSIDE each
/// Green exponential, avoiding underflow followed by an overflowing exp(s*tau).
/// Each advance must be <= the actual point-to-surface clearance/c. Output is
/// exp(s*advance)*H(s), including the same factor on the pressure gradient.
/// # Errors
/// Cross-wiring, receiver admission or an excessive/nonfinite delay advance.
pub fn observe_laplace(
    surface: &SpherePanels,
    solution: &LaplaceSolution,
    points: &[[f64; 3]],
    minimum_clearance_m: f64,
    advances_s: &[f64],
) -> Result<Vec<PointField>, HelmholtzError> {
    let geometry = crate::near_field::Geometry::new(surface, points, minimum_clearance_m)?;
    if solution.fingerprint != surface_fingerprint(surface) || advances_s.len() != points.len() {
        return Err(bad("cross-wired Laplace solution or delay shape"));
    }
    let k = C64::new(
        -solution.s.im / solution.medium.sound_speed,
        solution.s.re / solution.medium.sound_speed,
    );
    let mut fields = Vec::with_capacity(points.len());
    for ((&point, &advance), &clearance) in
        points.iter().zip(advances_s).zip(geometry.clearances_m())
    {
        if !advance.is_finite()
            || advance < 0.0
            || advance > clearance / solution.medium.sound_speed
        {
            return Err(bad("Laplace advance exceeds geometric propagation delay"));
        }
        let wave = DampedWaveNumber {
            k,
            phase_origin_m: advance * solution.medium.sound_speed,
        };
        fields.extend(fields_from_traces(
            surface,
            &solution.pressure,
            &solution.velocity,
            wave,
            solution.s.scale(-solution.medium.density),
            &[point],
        )?);
    }
    Ok(fields)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn g0_damped_kernels_and_real_axis_share_normal_derivatives() {
        let x = [0.4, 0.2, 0.7];
        let y = [-0.1, 0.0, 0.2];
        let nx = [1.0, 0.0, 0.0];
        let ny = [0.0, 0.0, 1.0];
        let wave = DampedWaveNumber {
            k: C64::new(2.3, 0.4),
            phase_origin_m: 0.0,
        };
        let (_, _, ds, dd) = kernels(wave, x, nx, y, ny);
        let eps = 1e-5;
        let plus = kernels(wave, [x[0] + eps, x[1], x[2]], nx, y, ny);
        let minus = kernels(wave, [x[0] - eps, x[1], x[2]], nx, y, ny);
        assert!((ds - (plus.0 - minus.0).scale(0.5 / eps)).abs() < 1e-8);
        assert!((dd - (plus.1 - minus.1).scale(0.5 / eps)).abs() < 1e-8);
        let a = kernels(2.3, x, nx, y, ny);
        let b = kernels(
            DampedWaveNumber {
                k: C64::from_re(2.3),
                phase_origin_m: 0.0,
            },
            x,
            nx,
            y,
            ny,
        );
        for (a, b) in [a.0, a.1, a.2, a.3].into_iter().zip([b.0, b.1, b.2, b.3]) {
            assert!((a - b).abs() < 2e-15);
        }
    }
}

