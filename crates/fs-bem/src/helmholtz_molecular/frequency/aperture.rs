//! Frequency-dependent self and mutual aperture loads and finite receivers.
use super::*;
use super::super::aperture as original_aperture;

/// Original aperture reaction with exact frequency slopes from the same source
/// solves. Shapes and reference planes stay fixed; each basis has unit imposed
/// generalized velocity, so its velocity derivative is zero.
#[derive(Clone, Debug)]
pub struct MolecularApertureFrequencyResponse {
    response: MolecularApertureResponse,
    impedance_derivative: Vec<C64>,
    pressure_derivatives: Vec<Vec<C64>>,
    dk: C64,
}

impl MolecularApertureFrequencyResponse {
    /// Original source-bound matrix, diagnostics and boundary solutions.
    pub fn response(&self) -> &MolecularApertureResponse {
        &self.response
    }

    /// Row-major dZ/domega [(N s/m)/(rad/s)] in the original aperture order.
    /// A frequency derivative is not a passive impedance and has no positivity
    /// constraint; the original response retains its own passivity diagnostic.
    pub fn impedance_derivative(&self) -> &[C64] {
        &self.impedance_derivative
    }

    /// Plane-port dZ/domega [(Pa s/m^3)/(rad/s)], with the same fixed volume
    /// gains as `response().impedance().acoustic(0)`. A zero-net-flow shape has
    /// no plane-port meaning and refuses; generalized impedance remains valid.
    pub fn acoustic_derivative(&self) -> Result<Vec<C64>, MolecularError> {
        let matrix = self.response.impedance();
        if matrix.volume_gains.iter().any(|q| !q.is_finite() || *q == 0.0) {
            return Err(bad("a zero-flow aperture has no plane-port impedance derivative"));
        }
        let n = matrix.apertures;
        let mut output = Vec::with_capacity(n * n);
        for (i, &value) in self.impedance_derivative.iter().enumerate() {
            let normalized = value.scale(1.0 / (matrix.volume_gains[i / n] * matrix.volume_gains[i % n]));
            if !finite(normalized) {
                return Err(bad("aperture acoustic derivative is unrepresentable"));
            }
            output.push(normalized);
        }
        Ok(output)
    }

    /// Observe the reaction's own outward volume flows and optional dU/domega.
    /// Both U*p' and U'*p enter the pressure tangent; the original fixed shape
    /// gives v'=sum shape*U'/gain. Traces are combined before a single receiver
    /// traversal, retaining mutual scattering and molecular propagation once.
    #[allow(clippy::too_many_arguments)]
    pub fn observe_volume_flows(
        &self,
        surface: &SpherePanels,
        volume_flows: &[C64],
        flow_derivatives: Option<&[C64]>,
        points: &[[f64; 3]],
        minimum_clearance_m: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<MolecularPointFrequencyResponse>, MolecularError> {
        poll(&mut cancelled)?;
        if flow_derivatives.is_some_and(|slopes| {
            slopes.len() != volume_flows.len() || slopes.iter().any(|&z| !finite(z))
        }) {
            return Err(bad("aperture receiver needs one finite frequency slope per volume flow"));
        }
        let combined = self.response.combined_volume_flows(volume_flows, &mut cancelled)?;
        let mut pressure_derivative = vec![C64::ZERO; combined.pressure.len()];
        let mut velocity_derivative = vec![C64::ZERO; combined.velocity.len()];
        for (b, solution) in self.response.solutions.iter().enumerate() {
            poll(&mut cancelled)?;
            let gain = self.response.matrix.volume_gains[b];
            let amplitude = volume_flows[b].scale(1.0 / gain);
            let slope = flow_derivatives.map_or(C64::ZERO, |values| values[b]).scale(1.0 / gain);
            if !finite(slope) {
                return Err(bad("aperture flow frequency slope is unrepresentable"));
            }
            for j in 0..pressure_derivative.len() {
                if j % 128 == 0 {
                    poll(&mut cancelled)?;
                }
                pressure_derivative[j] = pressure_derivative[j]
                    + self.pressure_derivatives[b][j] * amplitude + solution.pressure[j] * slope;
                velocity_derivative[j] = velocity_derivative[j] + solution.velocity[j] * slope;
            }
        }
        super::observation::observe_traces(
            &combined, &pressure_derivative, &velocity_derivative, self.dk,
            surface, points, minimum_clearance_m, cancelled,
        )
    }
}

impl MolecularExterior {
    /// Compute all aperture reaction values and frequency slopes with one
    /// original complex assembly and LU. Original shape/reference admission,
    /// work-conjugate force projection, diagnostics and acoustic normalization
    /// are retained. The derivative is never treated as a passive load itself.
    /// The complex four-matrix budget and existing source-count cap apply.
    pub fn aperture_response_with_frequency_derivative(
        &self,
        surface: &SpherePanels,
        apertures: &[Aperture],
        omega: f64,
        formulation: Formulation,
        options: &ComplexOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularApertureFrequencyResponse, MolecularError> {
        poll(&mut cancelled)?;
        original_aperture::admit_apertures(surface, apertures)?;
        let fields: Vec<Vec<_>> = apertures.iter()
            .map(|a| a.shape.iter().map(|&v| C64::from_re(v)).collect()).collect();
        let refs: Vec<_> = fields.iter().map(Vec::as_slice).collect();
        let solved = self.solve_batch_with_frequency_derivative(
            surface, omega, &refs, None, formulation, options, &mut cancelled,
        )?;
        let dk = solved[0].wavenumber_derivative_s_m;
        let mut solutions = Vec::with_capacity(solved.len());
        let mut pressure_derivatives = Vec::with_capacity(solved.len());
        for result in solved {
            solutions.push(result.response);
            pressure_derivatives.push(result.pressure_derivative_pa_per_rad_s);
        }
        poll(&mut cancelled)?;
        let derivative_refs: Vec<_> = pressure_derivatives.iter().map(Vec::as_slice).collect();
        let impedance_derivative = original_aperture::project_pressures(surface, apertures, &derivative_refs)?;
        let response = self.aperture_response_from_solutions(
            surface, apertures, omega, formulation, solutions, &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(MolecularApertureFrequencyResponse {
            response, impedance_derivative, pressure_derivatives, dk,
        })
    }
}

#[cfg(test)]
#[path = "aperture/tests.rs"]
mod tests;
