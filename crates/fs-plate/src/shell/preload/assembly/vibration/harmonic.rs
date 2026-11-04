//! Harmonic perturbations of the actual accepted loaded shell assembly.
//!
//! The same complete projected M/K as the transient owner and all glue C
//! entries enter fs-modal's coupled solve. Physical forces and observations
//! retain the original support/tie map. The static preload is never reapplied.

use super::{AssemblyVibration, checkpoint};
use fs_exec::Cx;
use fs_math::c64::C64;
use fs_modal::harmonic::{
    HarmonicFrequencyControls, HarmonicFrequencyResponse, HarmonicResponse,
    HarmonicResponseControls, PreparedHarmonicResponse,
};

/// Frequency-dependent loads at original work-conjugate physical attachments.
pub mod ports;

/// Cold projection limits, independent of each frequency factor's controls.
#[derive(Clone, Copy, Debug)]
pub struct AssemblyHarmonicBudget {
    /// Complete retained mode count; no mode is dropped to meet this limit.
    pub maximum_modes: usize,
    /// Upper bound on n*(nnz(M)+nnz(K)) + 2*n*n*free multiply-add terms.
    pub maximum_projection_terms: u64,
    /// Two projected matrices and two free-coordinate scratch vectors [bytes].
    /// Excludes borrowed modes, allocator metadata and the shared factor owner.
    pub maximum_projection_bytes: usize,
}

/// One immutable physical projection, reusable throughout a frequency sweep.
///
/// This is response in the supplied certified windows, not a full-spectrum
/// truncation certificate. Material/radiation loss, nonlinear motion and
/// changes of a connection's active branch remain outside this model.
pub struct AssemblyHarmonicModel<'v, 'body> {
    vibration: &'v AssemblyVibration<'body>,
    mass: Vec<f64>,
    stiffness: Vec<f64>,
}

/// A single frequency factor, reusable for multiple source and receiver cases.
pub struct PreparedAssemblyHarmonicResponse<'m, 'v, 'body> {
    model: &'m AssemblyHarmonicModel<'v, 'body>,
    response: PreparedHarmonicResponse,
    relative_tolerance: f64,
}

/// Peak complex amplitudes with the shared convention Re(amplitude exp(-iωt)).
#[derive(Debug)]
pub struct AssemblyHarmonicResponse {
    /// Full physical component perturbations [m/rad], without static preload.
    pub perturbation_displacements: Vec<Vec<C64>>,
    /// Full physical component velocity amplitudes [m/s, rad/s].
    pub velocities: Vec<Vec<C64>>,
    /// Original retained-coordinate solve, residual and cycle-mean work/loss.
    /// Its viscosity loss is exactly the supplied glue loss, counted once.
    pub modal: HarmonicResponse,
    /// Independent cycle-mean input power from physical component loads [W].
    /// At DC this is zero; a static load is not treated as an oscillation.
    pub average_component_power_w: f64,
}

/// The same physical assembled response and its exact frequency tangent.
/// Geometry, equilibrium, support/tie maps and the retained basis stay fixed.
#[derive(Debug)]
pub struct AssemblyHarmonicFrequencyResponse {
    /// Complete retained response/tangent, residuals and physical power slopes.
    pub modal: HarmonicFrequencyResponse,
    /// Physical displacement perturbation amplitudes, excluding static preload.
    pub perturbation_displacements: Vec<Vec<C64>>,
    /// Physical velocity amplitudes.
    pub velocities: Vec<Vec<C64>>,
    /// Physical displacement derivatives per rad/s, with the original tie map.
    pub displacement_derivatives: Vec<Vec<C64>>,
    /// Physical velocity derivatives per rad/s, including the kinematic term.
    pub velocity_derivatives: Vec<Vec<C64>>,
    /// Independently accumulated original physical source work [W].
    pub average_component_power_w: f64,
    /// Independently accumulated source-work derivative [W per rad/s].
    pub average_component_power_derivative: f64,
}

impl<'v, 'body> AssemblyHarmonicModel<'v, 'body> {
    /// Project once from the unchanged accepted equilibrium and its true modes.
    /// No exact-time propagator, inverse mass or extra eigenproblem is prepared.
    pub fn new(
        vibration: &'v AssemblyVibration<'body>,
        budget: AssemblyHarmonicBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, String> {
        checkpoint(cx)?;
        let n = vibration.mode_count();
        if n == 0 || n > budget.maximum_modes {
            return Err("assembled harmonic response requires complete mode capacity".into());
        }
        let (mass, stiffness) = vibration.projected_pencil(
            budget.maximum_projection_terms,
            budget.maximum_projection_bytes,
            cx,
        )?;
        if vibration.glue_damping().len() != mass.len() {
            return Err("assembled harmonic damping and source dimensions disagree".into());
        }
        checkpoint(cx)?;
        Ok(Self {
            vibration,
            mass,
            stiffness,
        })
    }

    /// Original Phi^T M Phi, retained without replacing it by identity.
    pub fn mass_matrix(&self) -> &[f64] {
        &self.mass
    }

    /// Original Phi^T K_t Phi, including finite off-diagonal projection terms.
    pub fn stiffness_matrix(&self) -> &[f64] {
        &self.stiffness
    }

    /// Complete Phi^T C_glue Phi, including off-diagonal coupling.
    pub fn damping_matrix(&self) -> &[f64] {
        self.vibration.glue_damping()
    }

    /// Factor the full K - ω²M - iωC through the existing shared owner.
    /// Undamped singular frequencies refuse; no damping/jitter is manufactured.
    pub fn prepare(
        &self,
        angular_frequency_rad_s: f64,
        controls: HarmonicResponseControls,
        cx: &Cx<'_>,
    ) -> Result<PreparedAssemblyHarmonicResponse<'_, 'v, 'body>, String> {
        let response = PreparedHarmonicResponse::try_new(
            &self.mass,
            self.damping_matrix(),
            &self.stiffness,
            self.vibration.mode_count(),
            angular_frequency_rad_s,
            controls,
            cx,
        )
        .map_err(|e| e.to_string())?;
        checkpoint(cx)?;
        Ok(PreparedAssemblyHarmonicResponse {
            model: self,
            response,
            relative_tolerance: controls.maximum_relative_residual,
        })
    }
}

impl PreparedAssemblyHarmonicResponse<'_, '_, '_> {
    /// Apply arbitrary phased component forces/moments in original coordinates.
    /// Fixed and eliminated slave coordinates use the same T^T load map and
    /// T observation map as equilibrium and transient response. The solver and
    /// loaded equilibrium are immutable, including on cancellation or refusal.
    pub fn solve(
        &self,
        component_loads: &[&[C64]],
        cx: &Cx<'_>,
    ) -> Result<AssemblyHarmonicResponse, String> {
        let force = self.project_loads(component_loads, cx)?;
        let modal = self.response.solve(&force, cx).map_err(|e| e.to_string())?;
        let perturbation_displacements = self.expand(&modal.displacement, cx)?;
        let velocities = self.expand(&modal.velocity, cx)?;
        let (power, scale) = component_power(component_loads, &velocities, cx)?;
        self.check_power(power, modal.average_input_power_w, scale)?;
        checkpoint(cx)?;
        Ok(AssemblyHarmonicResponse {
            perturbation_displacements,
            velocities,
            modal,
            average_component_power_w: power,
        })
    }

    /// Differentiate this frequency's physical forced response using the same
    /// cached full M/K/C factor. Optional source slopes use the same original
    /// component force/moment coordinates; None holds their phasors fixed.
    /// The full nonproportional glue loss and support/tie projection remain.
    /// This supplies actual physical receivers to resonance refinement without
    /// a frequency finite difference, extra factorization or modal decoupling.
    pub fn solve_with_frequency_derivative(
        &self,
        component_loads: &[&[C64]],
        component_load_derivatives: Option<&[&[C64]]>,
        controls: HarmonicFrequencyControls,
        cx: &Cx<'_>,
    ) -> Result<AssemblyHarmonicFrequencyResponse, String> {
        let force = self.project_loads(component_loads, cx)?;
        let force_derivative = component_load_derivatives
            .map(|loads| self.project_loads(loads, cx))
            .transpose()?;
        let modal = self
            .response
            .solve_with_frequency_derivative(
                &force,
                force_derivative.as_deref(),
                None,
                controls,
                cx,
            )
            .map_err(|e| e.to_string())?;
        let perturbation_displacements = self.expand(&modal.response.displacement, cx)?;
        let velocities = self.expand(&modal.response.velocity, cx)?;
        let displacement_derivatives = self.expand(&modal.displacement_derivative, cx)?;
        let velocity_derivatives = self.expand(&modal.velocity_derivative, cx)?;
        let (power, scale) = component_power(component_loads, &velocities, cx)?;
        self.check_power(power, modal.response.average_input_power_w, scale)?;
        let (mut slope, mut slope_scale) =
            component_power(component_loads, &velocity_derivatives, cx)?;
        if let Some(loads) = component_load_derivatives {
            let (source, source_scale) = component_power(loads, &velocities, cx)?;
            slope += source;
            slope_scale += source_scale;
        }
        self.check_power(slope, modal.average_input_power_derivative, slope_scale)?;
        checkpoint(cx)?;
        Ok(AssemblyHarmonicFrequencyResponse {
            modal,
            perturbation_displacements,
            velocities,
            displacement_derivatives,
            velocity_derivatives,
            average_component_power_w: power,
            average_component_power_derivative: slope,
        })
    }

    fn check_power(&self, physical: f64, modal: f64, scale: f64) -> Result<(), String> {
        if !physical.is_finite()
            || !scale.is_finite()
            || (physical - modal).abs() > self.relative_tolerance * scale.max(f64::MIN_POSITIVE)
        {
            return Err("assembled harmonic physical/modal cycle power disagrees".into());
        }
        Ok(())
    }

    fn project_loads(&self, component_loads: &[&[C64]], cx: &Cx<'_>) -> Result<Vec<C64>, String> {
        checkpoint(cx)?;
        let components = &self.model.vibration.equilibrium.displacements;
        if component_loads.len() != components.len()
            || component_loads
                .iter()
                .zip(components)
                .any(|(force, shape)| force.len() != shape.len())
        {
            return Err("assembled harmonic loads must match every physical component".into());
        }
        if component_loads
            .iter()
            .flat_map(|v| v.iter())
            .any(|v| !v.re.is_finite() || !v.im.is_finite())
        {
            return Err("assembled harmonic loads must be finite complex amplitudes".into());
        }
        let project = |imaginary: bool| {
            let loads: Vec<Vec<_>> = component_loads
                .iter()
                .map(|v| {
                    v.iter()
                        .map(|z| if imaginary { z.im } else { z.re })
                        .collect()
                })
                .collect();
            let loads: Vec<_> = loads.iter().map(Vec::as_slice).collect();
            self.model.vibration.project_loads(&loads, cx)
        };
        let real = project(false)?;
        let imaginary = project(true)?;
        let force: Vec<_> = real
            .into_iter()
            .zip(imaginary)
            .map(|(re, im)| C64::new(re, im))
            .collect();
        checkpoint(cx)?;
        Ok(force)
    }

    fn expand(&self, coordinates: &[C64], cx: &Cx<'_>) -> Result<Vec<Vec<C64>>, String> {
        let real: Vec<_> = coordinates.iter().map(|z| z.re).collect();
        let imaginary: Vec<_> = coordinates.iter().map(|z| z.im).collect();
        let real = self.model.vibration.perturbation_displacements(&real, cx)?;
        let imaginary = self
            .model
            .vibration
            .perturbation_displacements(&imaginary, cx)?;
        let result: Vec<Vec<_>> = real
            .into_iter()
            .zip(imaginary)
            .map(|(re, im)| {
                re.into_iter()
                    .zip(im)
                    .map(|(re, im)| C64::new(re, im))
                    .collect()
            })
            .collect();
        if result
            .iter()
            .flatten()
            .any(|z| !z.re.is_finite() || !z.im.is_finite())
        {
            return Err("assembled harmonic physical output overflow".into());
        }
        checkpoint(cx)?;
        Ok(result)
    }
}

fn component_power(
    loads: &[&[C64]],
    velocities: &[Vec<C64>],
    cx: &Cx<'_>,
) -> Result<(f64, f64), String> {
    let (mut power, mut scale) = (0.0_f64, 0.0_f64);
    for (force, velocity) in loads.iter().zip(velocities) {
        checkpoint(cx)?;
        for (&f, &v) in force.iter().zip(velocity) {
            power += 0.5 * (f.conj() * v).re;
            scale += 0.5 * f.abs() * v.abs();
        }
    }
    Ok((power, scale))
}

#[cfg(test)]
mod tests;
