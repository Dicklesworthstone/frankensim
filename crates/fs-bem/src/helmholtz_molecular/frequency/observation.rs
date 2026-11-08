//! Analytic frequency slopes of the original finite receiver representation.
use super::*;

/// Pressure, Cartesian gradient and their angular-frequency derivatives at a
/// fixed exterior receiver. Values use the original centroid Green quadrature;
/// this does not certify near-surface quadrature or spatial convergence.
#[derive(Clone, Copy, Debug)]
pub struct MolecularPointFrequencyResponse {
    pub field: PointField,
    /// d(p)/d(omega) [Pa/(rad/s)], including the imposed drive slope.
    pub pressure_derivative_pa_per_rad_s: C64,
    /// d(grad(p))/d(omega) [(Pa/m)/(rad/s)] in Cartesian order.
    pub gradient_derivative_pa_per_m_per_rad_s: [C64; 3],
}

impl MolecularPointFrequencyResponse {
    /// p'/p [s]: real part is d ln|p|/domega; imaginary part is group delay
    /// d arg(p)/domega for the original exp(-i omega t) convention. This
    /// includes source phase slopes and is not a wavefront or energy velocity.
    /// A zero pressure has no defined phase; nonfinite data or ratios refuse.
    pub fn pressure_log_frequency_derivative_s(&self) -> Result<C64, MolecularError> {
        if self.field.pressure_pa == C64::ZERO
            || !finite(self.field.pressure_pa)
            || !finite(self.pressure_derivative_pa_per_rad_s)
        {
            return Err(bad("pressure frequency slope requires a nonzero finite response and finite derivative"));
        }
        let slope = self.pressure_derivative_pa_per_rad_s / self.field.pressure_pa;
        if !finite(slope) {
            return Err(bad("pressure logarithmic frequency slope is unrepresentable"));
        }
        Ok(slope)
    }

    /// Group delay [s] of this pressure response, including any supplied drive
    /// phase slope. Positive delay has positive phase slope for exp(-i omega t).
    pub fn pressure_group_delay_s(&self) -> Result<f64, MolecularError> {
        self.pressure_log_frequency_derivative_s().map(|slope| slope.im)
    }
}

impl MolecularRadiationFrequencyResponse {
    /// Observe the same solved material, traces and their frequency tangent.
    /// The chain rule differentiates BOTH boundary traces and complex Green
    /// kernels; q'=i*rho*(v+omega*v') retains the original momentum law.
    /// Geometry and receiver caps/clearance are the original observer's. Three
    /// original point-kernel evaluations per panel produce pressure, all three
    /// gradient components and their slopes in one traversal, with cancellation
    /// between receivers and bounded panel tiles. No boundary solve or frequency
    /// difference occurs.
    pub fn observe(
        &self,
        surface: &SpherePanels,
        points: &[[f64; 3]],
        minimum_clearance_m: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<MolecularPointFrequencyResponse>, MolecularError> {
        observe_traces(
            &self.response,
            &self.pressure_derivative_pa_per_rad_s,
            &self.velocity_derivative_m_s_per_rad_s,
            self.wavenumber_derivative_s_m,
            surface,
            points,
            minimum_clearance_m,
            cancelled,
        )
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn observe_traces(
    response: &MolecularRadiation,
    pressure_derivative: &[C64],
    velocity_derivative: &[C64],
    dk: C64,
    surface: &SpherePanels,
    points: &[[f64; 3]],
    minimum_clearance_m: f64,
    mut cancelled: impl FnMut() -> bool,
) -> Result<Vec<MolecularPointFrequencyResponse>, MolecularError> {
    poll(&mut cancelled)?;
    let n = surface.areas().len();
    if response.fingerprint != surface_fingerprint(surface)
        || response.pressure.len() != n
        || response.velocity.len() != n
        || pressure_derivative.len() != n
        || velocity_derivative.len() != n
        || pressure_derivative.iter().chain(velocity_derivative).any(|&z| !finite(z))
        || !finite(dk)
    {
        return Err(bad("molecular receiver needs the original surface and complete finite tangents"));
    }
    crate::near_field::Geometry::new(surface, points, minimum_clearance_m)?;
    let rho = response.source.gas().density;
    let q_scale = C64::new(0.0, rho * response.omega);
    let q_derivative_scale = C64::new(0.0, rho);
    let mut output = Vec::with_capacity(points.len());
    for &point in points {
        poll(&mut cancelled)?;
        let mut values = [C64::ZERO; 4];
        let mut derivatives = [C64::ZERO; 4];
        for j in 0..n {
            if j % 128 == 0 {
                poll(&mut cancelled)?;
            }
            let p = response.pressure[j];
            let dp = pressure_derivative[j];
            let q = q_scale * response.velocity[j];
            let dq = q_derivative_scale * response.velocity[j] + q_scale * velocity_derivative[j];
            let area = surface.areas()[j];
            for axis in 0..3 {
                let mut normal = [0.0; 3];
                normal[axis] = 1.0;
                let (value, derivative) = complex::point_with_derivative::<true>(
                    response.k, point, normal, surface.centroids()[j], surface.normals()[j],
                );
                let mut contract = |component: usize, single: usize, double: usize| {
                    values[component] = values[component]
                        + (value[double] * p - value[single] * q).scale(area);
                    derivatives[component] = derivatives[component]
                        + (value[double] * dp - value[single] * dq
                            + dk * (derivative[double] * p - derivative[single] * q)).scale(area);
                };
                // The source-normal kernels are independent of the target
                // direction, so the first Cartesian call also supplies pressure.
                if axis == 0 {
                    contract(0, 0, 1);
                }
                contract(axis + 1, 2, 3);
            }
        }
        if values.iter().chain(&derivatives).any(|&z| !finite(z)) {
            return Err(bad("nonfinite molecular receiver field or frequency derivative"));
        }
        output.push(MolecularPointFrequencyResponse {
            field: PointField {
                pressure_pa: values[0],
                gradient_pa_per_m: [values[1], values[2], values[3]],
            },
            pressure_derivative_pa_per_rad_s: derivatives[0],
            gradient_derivative_pa_per_m_per_rad_s: [derivatives[1], derivatives[2], derivatives[3]],
        });
    }
    poll(&mut cancelled)?;
    Ok(output)
}

#[cfg(test)]
#[path = "observation/tests.rs"]
mod tests;
