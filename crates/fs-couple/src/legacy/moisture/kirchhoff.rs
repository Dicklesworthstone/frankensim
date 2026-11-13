//! Homogeneous nonlinear diffusion through the original source's potential.
//!
//! For Phi'(M)=D(M), the internal face flux is
//! rho_dry*(Phi(M_left)-Phi(M_right))/dx. A symmetric secant of Phi therefore
//! gives the original conservative M-matrix without freezing a dry D(0)=0
//! endpoint into an impermeable barrier. Bounded Picard iteration still checks
//! the complete endpoint equation; this is not an explicit wetting-front tracker.
//!
//! Preparation and exact piecewise-linear source integration belong to
//! fs-material. One chart covers the maximum-principle moisture interval for
//! the complete trial. No material callback or source search occurs per face.

use super::*;

/// The scalar surface root is monotone and bounded independently of Picard work.
const MAX_SURFACE_ITERATIONS: usize = 128;

impl MoistureSlab {
    pub(super) fn kirchhoff_profile(
        &self,
        surface: [Option<f64>; 2],
        prior_bounds: Option<[f64; 2]>,
        cx: &Cx<'_>,
    ) -> Result<KirchhoffMoistureDiffusivity, MoistureError> {
        poll(cx)?;
        let mut lower = f64::INFINITY;
        let mut upper: f64 = 0.0;
        for (i, m) in self
            .moisture
            .iter()
            .copied()
            .chain(surface.into_iter().flatten())
            .chain(prior_bounds.into_iter().flatten())
            .enumerate()
        {
            if i % 128 == 0 {
                poll(cx)?;
            }
            lower = lower.min(m);
            upper = upper.max(m);
        }
        let result = KirchhoffMoistureDiffusivity::from_binding(
            &self.card,
            &self.transport,
            &self.diffusivity_property,
            [lower, upper],
        )?;
        poll(cx)?;
        Ok(result)
    }
}

/// An original zero integral is admissible; underflow of positive transport is
/// not silently reinterpreted as a source-defined impermeable interface.
pub(super) fn scaled_diffusivity(scale: f64, diffusivity: f64) -> Result<f64, MoistureError> {
    if diffusivity == 0.0 {
        return Ok(0.0);
    }
    let coefficient = scale * diffusivity;
    if !positive(coefficient) {
        return Err(MoistureError::Invalid(
            "Kirchhoff face coefficient overflow or underflow",
        ));
    }
    Ok(coefficient)
}

fn relatively_equal(a: f64, b: f64) -> bool {
    a == b || (a - b).abs() / a.abs().max(b.abs()) <= 128.0 * f64::EPSILON
}

/// Half-cell potential difference in series with the existing surface law.
/// The returned fraction reconstructs the actual surface content from the cell
/// and reservoir values, just like the original constant-coefficient boundary.
#[allow(clippy::too_many_arguments)]
pub(super) fn surface_coupling(
    profile: &KirchhoffMoistureDiffusivity,
    face: &ResolvedMoistureFace,
    dt: f64,
    dx: f64,
    cell: f64,
    scale: f64,
    cx: &Cx<'_>,
) -> Result<(f64, f64), MoistureError> {
    poll(cx)?;
    let Some(equilibrium) = face.equilibrium else {
        return Ok((0.0, 0.0));
    };
    if face.coefficient_m_s == Some(0.0) {
        return Ok((0.0, 0.0));
    }
    let mean = profile.mean_diffusivity(cell, equilibrium)?;
    if face.coefficient_m_s.is_none() {
        let coefficient = 2.0 * scaled_diffusivity(scale, mean)?;
        if !coefficient.is_finite() {
            return Err(MoistureError::Invalid(
                "Kirchhoff boundary coefficient overflow",
            ));
        }
        return Ok((coefficient, 1.0));
    }
    if cell == equilibrium || mean == 0.0 {
        return face.coupling(dt, dx, mean, scale);
    }

    // Set M_surface=M_cell+f*(M_eq-M_cell), with f in [0,1]. The
    // scalar balance is 2*mean_D(M_cell,M_surface)*f/dx=h*(1-f).
    // Its unnormalised residual is monotone because D>=0 and h>0,
    // including desorption. Using the stable series fraction avoids products
    // of large transfer speeds and loss of a small flux through 1-f subtraction.
    let mut bracket: [f64; 2] = [0.0, 1.0];
    let mut iterations = 0;
    for iteration in 1..=MAX_SURFACE_ITERATIONS {
        iterations = iteration;
        poll(cx)?;
        let fraction = bracket[0] + 0.5 * (bracket[1] - bracket[0]);
        let surface = fraction.mul_add(equilibrium - cell, cell);
        let local_mean = profile.mean_diffusivity(cell, surface)?;
        let coupling = face.coupling(dt, dx, local_mean, scale)?;

        // Recheck the actual reported surface against the SAME source integral.
        // This also retains a very small representable flux when the fraction
        // rounds to one: no subtraction of two nearly equal reservoir values
        // is used to manufacture the coefficient.
        let reported_surface = coupling.1.mul_add(equilibrium - cell, cell);
        let checked_mean = profile.mean_diffusivity(cell, reported_surface)?;
        let checked = face.coupling(dt, dx, checked_mean, scale)?;
        if relatively_equal(coupling.0, checked.0) && relatively_equal(coupling.1, checked.1) {
            return Ok(checked);
        }
        if fraction == bracket[0] || fraction == bracket[1] {
            break;
        }
        if fraction > coupling.1 {
            bracket[1] = fraction;
        } else {
            bracket[0] = fraction;
        }
    }
    Err(MoistureError::KirchhoffSurfaceNotConverged {
        iterations,
        fraction_bracket: bracket,
    })
}
