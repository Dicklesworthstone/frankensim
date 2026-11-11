//! Driven frequency derivatives with the original molecular dispersion.
//!
//! Differentiate A(k)*p=B(k)*q with q=i*rho*omega*v. One original complex-k
//! assembly and one LU serve the values and tangents of every imposed drive.

use super::*;

#[path = "frequency/observation.rs"]
mod observation;
pub use observation::MolecularPointFrequencyResponse;

#[path = "frequency/aperture.rs"]
mod aperture;
pub use aperture::MolecularApertureFrequencyResponse;

/// Original radiation and its analytic tangent per rad/s, with geometry,
/// equilibrium gas and reservoir parameters held. The derivative includes the
/// original molecular phase AND attenuation slopes and any supplied drive slope.
#[derive(Clone, Debug)]
pub struct MolecularRadiationFrequencyResponse {
    pub response: MolecularRadiation,
    /// dk/domega [s/m] from the same original source as the boundary response.
    pub wavenumber_derivative_s_m: C64,
    /// Boundary pressure derivative [Pa/(rad/s)], in original panel order.
    pub pressure_derivative_pa_per_rad_s: Vec<C64>,
    /// Imposed outward velocity derivative [(m/s)/(rad/s)].
    pub velocity_derivative_m_s_per_rad_s: Vec<C64>,
    /// Signed derivative of raw discrete mean boundary work [W/(rad/s)].
    /// This includes work eventually absorbed in the gas, not only far flux.
    pub boundary_work_derivative_w_per_rad_s: f64,
}

impl MolecularExterior {
    /// Solve prescribed outward velocities and their frequency derivatives.
    /// `None` holds all velocity phasors fixed; supplied slopes have one finite
    /// complete panel field per drive and are derivatives with respect to
    /// angular frequency, not Hz. This is the prescribed-velocity boundary
    /// problem of [`Self::solve_batch`]; no surface admittance is inferred.
    ///
    /// The exact differentiated equation is
    /// `A*p' = B*q' + k'*(B_k*q - A_k*p)`, with
    /// `q'=i*rho*(v+omega*v')`. Replacing k' by 1/c or q by i*rho*c*k*v
    /// would omit material dispersion or change the physical momentum law.
    ///
    /// The original complex assembly admits four matrices (64*N*N bytes)
    /// under `ComplexOptions`; its documented budget excludes LU and returned
    /// fields. One factorization and two back-solves per drive are used. The
    /// ordinary values-only solve retains its existing two-matrix allocation.
    /// Cancellation is checked before admission, at every assembly/RHS row,
    /// inside the original LU/solves, between drives and before publication.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_batch_with_frequency_derivative(
        &self,
        surface: &SpherePanels,
        omega: f64,
        fields: &[&[C64]],
        field_derivatives: Option<&[&[C64]]>,
        formulation: Formulation,
        options: &ComplexOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<MolecularRadiationFrequencyResponse>, MolecularError> {
        poll(&mut cancelled)?;
        if !omega.is_finite() || omega < self.band_rad_s[0] || omega > self.band_rad_s[1] {
            return Err(bad(
                "molecular exterior frequency is outside its declared band",
            ));
        }
        let n = surface.areas().len();
        let valid_field = |v: &[C64]| v.len() == n && v.iter().all(|&z| finite(z));
        if fields.is_empty()
            || fields.len() > MAX_RADIATION_FIELDS_PER_BATCH
            || fields.iter().any(|v| !valid_field(v))
            || field_derivatives.is_some_and(|slopes| {
                slopes.len() != fields.len() || slopes.iter().any(|v| !valid_field(v))
            })
        {
            return Err(bad(
                "molecular frequency solve needs bounded complete finite drives and slopes",
            ));
        }
        let rho = self.source.gas().density;
        let q_scale = C64::new(0.0, rho * omega);
        let q_scale_derivative = C64::new(0.0, rho);
        if !finite(q_scale) {
            return Err(bad(
                "molecular exterior boundary scales are unrepresentable",
            ));
        }
        let material = self.source.relaxing_frequency_response(omega)?;
        let k = material.wavenumber_per_m;
        let dk = material.wavenumber_derivative_s_m;
        let op = complex::assemble_with_cancel(surface, k, formulation, options, &mut cancelled)?;
        poll(&mut cancelled)?;
        let lu = lu_complex_with_cancel(&op.a, n, &mut cancelled)?;
        // The original A is no longer needed; the factor owns its coefficients.
        drop(op.a);
        let zero = vec![C64::ZERO; n];
        let mut results = Vec::with_capacity(fields.len());
        let fingerprint = surface_fingerprint(surface);
        for (case, &velocity) in fields.iter().enumerate() {
            poll(&mut cancelled)?;
            let slope = field_derivatives.map_or(zero.as_slice(), |values| values[case]);
            let mut pressure = vec![C64::ZERO; n];
            for i in 0..n {
                poll(&mut cancelled)?;
                for (j, &v) in velocity.iter().enumerate() {
                    pressure[i] = pressure[i] + op.b[i * n + j] * (q_scale * v);
                }
            }
            lu.solve_with_cancel(&mut pressure, &mut cancelled)?;
            let mut derivative = vec![C64::ZERO; n];
            for i in 0..n {
                poll(&mut cancelled)?;
                let mut geometry = C64::ZERO;
                let mut source = C64::ZERO;
                for j in 0..n {
                    let index = i * n + j;
                    let q = q_scale * velocity[j];
                    let dq = q_scale_derivative * velocity[j] + q_scale * slope[j];
                    source = source + op.b[index] * dq;
                    geometry = geometry + op.db_dk[index] * q - op.da_dk[index] * pressure[j];
                }
                derivative[i] = source + dk * geometry;
            }
            lu.solve_with_cancel(&mut derivative, &mut cancelled)?;
            let mut boundary_work_w = 0.0;
            let mut boundary_work_derivative_w_per_rad_s = 0.0;
            for i in 0..n {
                if !finite(pressure[i]) || !finite(derivative[i]) {
                    return Err(bad("nonfinite molecular pressure or frequency derivative"));
                }
                let area = surface.areas()[i];
                boundary_work_w += 0.5 * (pressure[i] * velocity[i].conj()).re * area;
                boundary_work_derivative_w_per_rad_s += 0.5
                    * (derivative[i] * velocity[i].conj() + pressure[i] * slope[i].conj()).re
                    * area;
            }
            if !boundary_work_w.is_finite() || !boundary_work_derivative_w_per_rad_s.is_finite() {
                return Err(bad("nonfinite molecular boundary work or its derivative"));
            }
            results.push(MolecularRadiationFrequencyResponse {
                response: MolecularRadiation {
                    source: Arc::clone(&self.source),
                    omega,
                    k,
                    fingerprint,
                    pressure,
                    velocity: velocity.to_vec(),
                    boundary_work_w,
                },
                wavenumber_derivative_s_m: dk,
                pressure_derivative_pa_per_rad_s: derivative,
                velocity_derivative_m_s_per_rad_s: slope.to_vec(),
                boundary_work_derivative_w_per_rad_s,
            });
        }
        poll(&mut cancelled)?;
        Ok(results)
    }
}

#[cfg(test)]
#[path = "frequency/tests.rs"]
mod tests;
