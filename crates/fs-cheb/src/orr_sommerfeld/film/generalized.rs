//! Generalized-Newtonian free-surface spectrum of an original resolved base.
//!
//! On the parent's thickness/mean-speed scales, let mu be apparent viscosity
//! and M=d(mu*gamma)/d(gamma), divided by the SAME reference viscosity. The
//! viscous operator is (D²+a²)[M(D²+a²)phi]-4a²D(mu Dphi), not a Newtonian
//! operator with a fitted viscosity. Convection retains the actual U and U''.
//! At the surface, eta=-(D²+a²)phi/U''_s and the viscous normal stress is
//! D[M(D²+a²)phi]-4a²mu Dphi, including M'_s. See the primary derivation cited
//! by the parent module, equations (20), (22)-(29), with wall-outward y.
//!
//! This module assembles only the physical coefficients. The existing native
//! basis, equilibration, LU/QR eigensolve, paired-order matching and original
//! pencil SVD residual remain the sole numerical solver and acceptance path.

use super::{
    FilmEigenmode, FilmParameters, FilmSpectrumError, Grid, Pencil, PreparedFilmSpectrum, poll,
    positive, validate_parameters,
};
use fs_math::c64::C64;

/// One original base point on y/h in [0,1], on common h, Ubar and mu_ref scales.
///
/// These are the actual constitutive/base fields, not independent fit parameters.
/// The caller supplies a fully developed, no-slip, shear-free base with unit
/// depth-averaged velocity, nonnegative shear and positive viscosities.
/// Coefficients must have finite one-sided derivatives through second order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FilmProfilePoint {
    pub velocity: f64,
    /// d(U/Ubar)/d(y/h).
    pub velocity_gradient: f64,
    /// d²(U/Ubar)/d(y/h)², from original momentum and stress tangent.
    pub velocity_curvature: f64,
    /// mu/mu_ref, multiplying normal strain perturbations.
    pub apparent_viscosity: f64,
    /// d(mu*gamma)/d(gamma)/mu_ref, multiplying shear perturbations.
    pub stress_tangent: f64,
}

/// The same immutable source at both native orders, prepared once per section.
///
/// Native differentiation resolves viscosity variations once, so wavelength
/// queries reuse the full base and never resample or refit the material law.
pub struct PreparedGeneralizedFilm<'a> {
    solver: &'a PreparedFilmSpectrum,
    coarse: Profile,
    fine: Profile,
    surface_stress_tangent_gradient: f64,
}

impl PreparedFilmSpectrum {
    /// Required sample coordinates, descending from free surface (1) to wall (0).
    /// The two slices refer to the original cached coarse/fine coefficient grids.
    pub fn profile_coordinates(&self) -> (&[f64], &[f64]) {
        (&self.coarse.nodes, &self.fine.nodes)
    }

    /// Prepare coefficients from the SAME physical source at both grid orders.
    ///
    /// The exact one-sided d(M/mu_ref)/d(y/h) at the free surface is supplied
    /// separately; it enters normal stress and is not silently set to zero.
    /// For a smooth nonconstant Carreau law with transition exponent a>=2 it
    /// is exactly zero. Native nodal differentiation supplies interior mu', M',
    /// M''. Paired-order eigenvalue admission also resolves these coefficients.
    pub fn prepare_profile(
        &self,
        coarse: &[FilmProfilePoint],
        fine: &[FilmProfilePoint],
        surface_stress_tangent_gradient: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedGeneralizedFilm<'_>, FilmSpectrumError> {
        poll(&mut cancelled)?;
        if !surface_stress_tangent_gradient.is_finite() {
            return Err(FilmSpectrumError::Invalid(
                "unrepresentable surface differential-viscosity gradient",
            ));
        }
        let coarse = Profile::new(coarse, self.coarse.nodes.len(), &mut cancelled)?;
        let fine = Profile::new(fine, self.fine.nodes.len(), &mut cancelled)?;
        if coarse.points[0] != fine.points[0] || coarse.points.last() != fine.points.last() {
            return Err(FilmSpectrumError::Invalid(
                "both spectral grids must retain the same physical boundary state",
            ));
        }
        poll(&mut cancelled)?;
        Ok(PreparedGeneralizedFilm {
            solver: self,
            coarse,
            fine,
            surface_stress_tangent_gradient,
        })
    }
}

impl PreparedGeneralizedFilm<'_> {
    /// Rightmost temporal modes, with exactly the parent's residual/resolution
    /// admission. Re uses the SAME reference viscosity as every profile point.
    /// No spatial amplification or universal turbulence criterion is inferred.
    pub fn solve(
        &self,
        parameters: FilmParameters,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<FilmEigenmode>, FilmSpectrumError> {
        poll(&mut cancelled)?;
        validate_parameters(parameters)?;
        let coarse = self.solver.coarse.generalized_pencil(
            parameters,
            &self.coarse,
            self.surface_stress_tangent_gradient,
            &mut cancelled,
        )?;
        let fine = self.solver.fine.generalized_pencil(
            parameters,
            &self.fine,
            self.surface_stress_tangent_gradient,
            &mut cancelled,
        )?;
        self.solver.solve_pencils(coarse, fine, cancelled)
    }
}

struct Profile {
    points: Vec<FilmProfilePoint>,
    mu_gradient: Vec<f64>,
    tangent_gradient: Vec<f64>,
    tangent_curvature: Vec<f64>,
}

impl Profile {
    fn new(
        points: &[FilmProfilePoint],
        count: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, FilmSpectrumError> {
        if points.len() != count {
            return Err(FilmSpectrumError::Invalid("generalized film profile shape"));
        }
        for point in points {
            poll(cancelled)?;
            if !point.velocity.is_finite()
                || point.velocity < 0.0
                || !point.velocity_gradient.is_finite()
                || point.velocity_gradient < 0.0
                || !point.velocity_curvature.is_finite()
                || point.velocity_curvature >= 0.0
                || !positive(point.apparent_viscosity)
                || !positive(point.stress_tangent)
            {
                return Err(FilmSpectrumError::Invalid(
                    "invalid original generalized film base field",
                ));
            }
        }
        let surface = points[0];
        let wall = points[count - 1];
        if !positive(surface.velocity)
            || surface.velocity_gradient != 0.0
            || wall.velocity != 0.0
            || !positive(wall.velocity_gradient)
            || (surface.apparent_viscosity - surface.stress_tangent).abs()
                > 64.0 * f64::EPSILON * surface.apparent_viscosity
        {
            return Err(FilmSpectrumError::Invalid(
                "generalized film needs original no-slip/shear-free boundary values",
            ));
        }
        // Original Chebyshev derivative, transformed from [-1,1] to [0,1].
        let derivative = crate::diff_matrix(count - 1);
        let mu: Vec<_> = points
            .iter()
            .map(|point| point.apparent_viscosity)
            .collect();
        let tangent: Vec<_> = points.iter().map(|point| point.stress_tangent).collect();
        let mu_gradient = differentiate(&derivative, &mu, cancelled)?;
        let tangent_gradient = differentiate(&derivative, &tangent, cancelled)?;
        let tangent_curvature = differentiate(&derivative, &tangent_gradient, cancelled)?;
        Ok(Self {
            points: points.to_vec(),
            mu_gradient,
            tangent_gradient,
            tangent_curvature,
        })
    }
}

fn differentiate(
    derivative: &[f64],
    values: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>, FilmSpectrumError> {
    let n = values.len();
    let mut result = vec![0.0; n];
    for i in 0..n {
        poll(cancelled)?;
        let mut sum = 0.0;
        for j in 0..n {
            // Row-difference form annihilates a constant exactly, preserving the
            // exact Newtonian coefficient limit without subtracting large sums.
            if i != j {
                sum = derivative[i * n + j].mul_add(values[j] - values[i], sum);
            }
        }
        result[i] = 2.0 * sum;
        if !result[i].is_finite() {
            return Err(FilmSpectrumError::Invalid(
                "generalized viscosity derivative is not representable",
            ));
        }
    }
    Ok(result)
}

impl Grid {
    fn generalized_pencil(
        &self,
        data: FilmParameters,
        profile: &Profile,
        surface_tangent_gradient: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Pencil, FilmSpectrumError> {
        let n = self.nodes.len();
        let mut a = vec![C64::ZERO; n * n];
        let mut b = vec![C64::ZERO; n * n];
        let [d0, d1, d2, d3, d4] = &self.derivatives;
        let alpha = data.kh;
        let alpha2 = alpha * alpha;
        let alpha4 = alpha2 * alpha2;
        for i in 1..(n - 1) {
            poll(cancelled)?;
            let point = profile.points[i];
            let mu = point.apparent_viscosity;
            let m = point.stress_tangent;
            let dm = profile.tangent_gradient[i];
            let d2m = profile.tangent_curvature[i];
            let dmu = profile.mu_gradient[i];
            for j in 0..n {
                let source = i * n + j;
                let target = (i - 1) * n + j;
                let laplace = d2[source] - alpha2 * d0[source];
                // Expanded product rule retains analytic derivatives of the
                // wall-factored basis. Differentiating sampled products instead
                // would alias its two highest-degree columns (degree N+2).
                let viscous = m * d4[source]
                    + 2.0 * dm * d3[source]
                    + (d2m + 2.0 * alpha2 * m - 4.0 * alpha2 * mu) * d2[source]
                    + (2.0 * alpha2 * dm - 4.0 * alpha2 * dmu) * d1[source]
                    + (alpha2 * d2m + alpha4 * m) * d0[source];
                a[target] = C64::new(
                    viscous / data.reynolds,
                    -alpha * (point.velocity * laplace - point.velocity_curvature * d0[source]),
                );
                b[target] = C64::from_re(laplace);
            }
        }
        let surface = profile.points[0];
        let restoring = data.normal_gravity + data.capillarity * alpha2;
        for j in 0..n {
            let eta = -(d2[j] + alpha2 * d0[j]) / surface.velocity_curvature;
            let viscous = surface.stress_tangent * d3[j]
                + surface_tangent_gradient * (d2[j] + alpha2 * d0[j])
                + alpha2 * (surface.stress_tangent - 4.0 * surface.apparent_viscosity) * d1[j];
            let normal = (n - 2) * n + j;
            a[normal] = C64::new(
                viscous / data.reynolds,
                -alpha * (surface.velocity * d1[j] + restoring * eta),
            );
            b[normal] = C64::from_re(d1[j]);
            let kinematic = (n - 1) * n + j;
            a[kinematic] = C64::new(0.0, -alpha * (d0[j] + surface.velocity * eta));
            b[kinematic] = C64::from_re(eta);
        }
        Pencil::equilibrated(a, b, n, cancelled)
    }
}

#[cfg(test)]
mod tests;
