//! Element-local pressure, particle velocity and acoustic energy flow.
//!
//! Gradients use the original hierarchical basis and physical element map.
//! A probe retains one coherent local material sample from the caller's
//! original frozen source, rather than using the geometry's uniform medium.
//! At a material interface, the explicitly selected tetrahedron determines
//! the one-sided trace; gradients are never averaged across the interface.

use super::time::HeterogeneousTimeState;
use super::{
    AcousticFeError, HeterogeneousError, HeterogeneousHelmholtzFe, Medium, finite, invalid,
    medium_coefficients,
};
use fs_math::c64::C64;

/// Complex peak fields in `exp(+i omega t)` and time-averaged energy flow.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicAcousticSample {
    /// Pressure [Pa], coherently interpolated before forming intensity.
    pub pressure_pa: C64,
    /// Physical, element-local gradient of pressure [Pa/m].
    pub pressure_gradient_pa_m: [C64; 3],
    /// Particle velocity `i grad(p)/(omega rho)` [m/s].
    pub particle_velocity_m_s: [C64; 3],
    /// Transported energy flux `Re(p conj(v))/2` [W/m^2].
    pub active_intensity_w_m2: [f64; 3],
    /// Reactive energy flux `Im(p conj(v))/2` [W/m^2]. This is not bulk heat.
    pub reactive_intensity_w_m2: [f64; 3],
}

/// Actual instantaneous fields of the pressure-potential formulation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransientAcousticSample {
    /// Pressure `p = psi_dot` [Pa].
    pub pressure_pa: f64,
    /// Pressure potential `psi` [Pa s]; its constant offset is a gauge.
    pub pressure_potential_pa_s: f64,
    /// Physical, element-local potential gradient [Pa s/m].
    pub potential_gradient_pa_s_m: [f64; 3],
    /// Particle velocity `-grad(psi)/rho` [m/s].
    pub particle_velocity_m_s: [f64; 3],
    /// Instantaneous transported energy flux `p v` [W/m^2].
    pub instantaneous_intensity_w_m2: [f64; 3],
}

struct ProbeShape {
    dof: usize,
    value: f64,
    gradient_m_inv: [f64; 3],
}

/// Reusable observation at one point in the original spatial/material owner.
///
/// Preparation retains the local high-order value/gradient stencil and one
/// original material sample. Subsequent fields require O(local DOFs) work and
/// no allocation, geometry reconstruction, material query or differencing.
/// A probe is valid only for that same frozen equilibrium; callers retain
/// responsibility for passing the original source when preparing it.
pub struct HeterogeneousProbe<'m, 'f, 'c> {
    model: &'m HeterogeneousHelmholtzFe<'f, 'c>,
    element: usize,
    position_m: [f64; 3],
    medium: Medium,
    inverse_density: f64,
    shapes: Vec<ProbeShape>,
}

impl<'f, 'c> HeterogeneousHelmholtzFe<'f, 'c> {
    /// Prepare a one-sided local pressure/velocity/intensity observation.
    ///
    /// `material_at` has the same source contract as [`Self::new`]: evaluate
    /// the ORIGINAL frozen material at the supplied original element and SI
    /// position. It is called exactly once after location validation, and
    /// its error propagates unchanged. Density and speed undergo the same
    /// checks as assembly; neither is interpolated from quadrature samples.
    /// This lower FE layer cannot certify caller-owned material identity.
    ///
    /// `barycentric` follows the selected tetrahedron's stored vertex order.
    /// On a shared face, select each neighboring element separately to
    /// observe its trace. FE gradient/flux error still requires h/p and
    /// material-quadrature convergence; continuity is not imposed afterward.
    pub fn probe<'m, E: From<AcousticFeError>>(
        &'m self,
        element: usize,
        barycentric: [f64; 4],
        material_at: &mut impl FnMut(usize, [f64; 3]) -> Result<Medium, E>,
    ) -> Result<HeterogeneousProbe<'m, 'f, 'c>, E> {
        self.validate_sample(self.ndof(), element, barycentric)?;
        let mut position_m = [0.0; 3];
        for (&vertex, &bary) in self.geometry.tetrahedra()[element].iter().zip(&barycentric) {
            for (axis, position) in position_m.iter_mut().enumerate() {
                *position += bary * self.geometry.positions[vertex as usize][axis];
            }
        }
        if position_m.iter().any(|x| !x.is_finite()) {
            return Err(invalid("unrepresentable heterogeneous observation point").into());
        }
        let medium = material_at(element, position_m)?;
        let (inverse_density, _) = medium_coefficients(medium)?;
        let dofs = self.geometry.space.element_dofs(element);
        let mut shapes = Vec::with_capacity(dofs.len());
        for (dof, shape) in dofs {
            let value = shape.eval(barycentric, self.geometry.order);
            let mut gradient_m_inv = [0.0; 3];
            for (a, derivative) in shape
                .d_lambda(barycentric, self.geometry.order)
                .into_iter()
                .enumerate()
            {
                for (axis, gradient) in gradient_m_inv.iter_mut().enumerate() {
                    *gradient += derivative * self.barycentric_gradients[element][a][axis];
                }
            }
            if !value.is_finite() || gradient_m_inv.iter().any(|g| !g.is_finite()) {
                return Err(invalid("unrepresentable heterogeneous observation stencil").into());
            }
            shapes.push(ProbeShape {
                dof,
                value,
                gradient_m_inv,
            });
        }
        Ok(HeterogeneousProbe {
            model: self,
            element,
            position_m,
            medium,
            inverse_density,
            shapes,
        })
    }
}

impl HeterogeneousProbe<'_, '_, '_> {
    /// Original tetrahedron defining this one-sided observation.
    #[must_use]
    pub const fn element(&self) -> usize {
        self.element
    }

    /// Physical observation position [m].
    #[must_use]
    pub const fn position_m(&self) -> [f64; 3] {
        self.position_m
    }

    /// Exact original local material sample retained during preparation.
    #[must_use]
    pub const fn medium(&self) -> Medium {
        self.medium
    }

    /// Observe original hierarchical complex peak pressure coefficients at
    /// their physical frequency. Use the same frequency as the harmonic solve.
    /// All contributions are superposed before forming quadratic intensity,
    /// preserving interference and evanescent reactive flow.
    ///
    /// The coefficient vector must belong to this original model. Its length
    /// and locally contributing entries are checked; a bare slice carries no
    /// independently verifiable model or frequency identity.
    pub fn harmonic(
        &self,
        pressure_pa: &[C64],
        frequency_hz: f64,
    ) -> Result<HarmonicAcousticSample, HeterogeneousError> {
        let omega = core::f64::consts::TAU * frequency_hz;
        let velocity_factor = self.inverse_density / omega;
        if pressure_pa.len() != self.model.ndof()
            || !omega.is_finite()
            || omega <= 0.0
            || !velocity_factor.is_finite()
            || velocity_factor <= 0.0
        {
            return Err(invalid(
                "harmonic observation requires original coefficients and representable positive frequency",
            )
            .into());
        }
        let mut pressure = C64::ZERO;
        let mut gradient = [C64::ZERO; 3];
        for shape in &self.shapes {
            let coefficient = pressure_pa[shape.dof];
            pressure = pressure + coefficient.scale(shape.value);
            for (axis, component) in gradient.iter_mut().enumerate() {
                *component = *component + coefficient.scale(shape.gradient_m_inv[axis]);
            }
        }
        let velocity = gradient.map(|g| C64::new(-g.im, g.re).scale(velocity_factor));
        let intensity = velocity.map(|v| pressure.scale(0.5) * v.conj());
        if !finite(pressure)
            || gradient.iter().chain(&velocity).chain(&intensity).any(|&v| !finite(v))
        {
            return Err(invalid("nonfinite harmonic acoustic observation").into());
        }
        Ok(HarmonicAcousticSample {
            pressure_pa: pressure,
            pressure_gradient_pa_m: gradient,
            particle_velocity_m_s: velocity,
            active_intensity_w_m2: intensity.map(|i| i.re),
            reactive_intensity_w_m2: intensity.map(|i| i.im),
        })
    }

    /// Observe a transient state belonging to this exact spatial/material
    /// operator. Velocity comes from potential, not from pressure or its rate;
    /// the instantaneous flux has no harmonic factor of one half.
    pub fn transient(
        &self,
        state: &HeterogeneousTimeState<'_, '_, '_>,
    ) -> Result<TransientAcousticSample, HeterogeneousError> {
        if !core::ptr::eq(self.model, state.model) {
            return Err(invalid("transient observation belongs to another spatial/material operator").into());
        }
        let mut pressure = 0.0;
        let mut potential = 0.0;
        let mut gradient = [0.0; 3];
        for shape in &self.shapes {
            pressure += state.pressure_pa()[shape.dof] * shape.value;
            let coefficient = state.pressure_potential_pa_s()[shape.dof];
            potential += coefficient * shape.value;
            for (axis, component) in gradient.iter_mut().enumerate() {
                *component += coefficient * shape.gradient_m_inv[axis];
            }
        }
        let velocity = gradient.map(|g| -g * self.inverse_density);
        let intensity = velocity.map(|v| pressure * v);
        if !pressure.is_finite()
            || !potential.is_finite()
            || gradient.iter().chain(&velocity).chain(&intensity).any(|v| !v.is_finite())
        {
            return Err(invalid("nonfinite transient acoustic observation").into());
        }
        Ok(TransientAcousticSample {
            pressure_pa: pressure,
            pressure_potential_pa_s: potential,
            potential_gradient_pa_s_m: gradient,
            particle_velocity_m_s: velocity,
            instantaneous_intensity_w_m2: intensity,
        })
    }
}

#[cfg(test)]
mod tests;
