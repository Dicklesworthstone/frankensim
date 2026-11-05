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

mod surface;
pub use surface::{
    ElementFace, HarmonicSurfacePower, HeterogeneousSurfaceProbe, SurfaceProbeLimits,
};

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

/// Bounded inverse-map preparation for an observation specified in SI space.
///
/// The caller selects the original tetrahedron/material side. A successful
/// location lies in its closed reference simplex, and its mapped position
/// differs from the requested point by at most `absolute_tolerance_m` on each
/// Cartesian axis. This is a residual criterion, not a global injectivity or
/// point-containment certificate for the mesh.
#[derive(Clone, Copy, Debug)]
pub struct ProbeLocationLimits {
    /// Maximum Newton updates, from 1 through 128.
    pub max_iterations: usize,
    /// Maximum successively halved trial steps per update, from 1 through 64.
    pub max_backtracks: usize,
    /// Positive finite physical-position tolerance [m]. Choose it independently
    /// of the field discretization error; no unitless tolerance replaces it.
    pub absolute_tolerance_m: f64,
}

// Euclidean projection of a Newton TRIAL onto the four-coordinate simplex.
// Truncating the Newton direction at its first face can stall on curved
// boundaries; projected trials can continue tangentially. Every trial is
// evaluated against the unchanged physical target before it is accepted.
fn reference_trial(mut lambda: [f64; 4]) -> Option<[f64; 4]> {
    if lambda.iter().any(|x| !x.is_finite()) {
        return None;
    }
    if lambda.iter().any(|x| !(0.0..=1.0).contains(x)) {
        let mut sorted = lambda;
        sorted.sort_by(|a, b| b.total_cmp(a));
        let mut sum = 0.0;
        let mut threshold = None;
        for (i, &value) in sorted.iter().enumerate() {
            sum += value;
            let candidate = (sum - 1.0) / (i + 1) as f64;
            if candidate.is_finite() && value > candidate {
                threshold = Some(candidate);
            }
        }
        let threshold = threshold?;
        lambda = lambda.map(|x| (x - threshold).max(0.0));
    }
    let mut anchor = 0;
    for i in 1..4 {
        if lambda[i] > lambda[anchor] {
            anchor = i;
        }
    }
    lambda[anchor] = 1.0
        - (0..4)
            .filter(|&i| i != anchor)
            .map(|i| lambda[i])
            .sum::<f64>();
    lambda
        .iter()
        .all(|x| x.is_finite() && (0.0..=1.0).contains(x))
        .then_some(lambda)
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
    reference_barycentric: [f64; 4],
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
    /// observe its trace. The original map supplies the physical position and
    /// pointwise inverse-transpose Jacobian, including on curved elements.
    /// FE gradient/flux error still requires h/p, geometry and material-quadrature
    /// convergence; continuity is not imposed afterward.
    pub fn probe<'m, E: From<AcousticFeError>>(
        &'m self,
        element: usize,
        barycentric: [f64; 4],
        material_at: &mut impl FnMut(usize, [f64; 3]) -> Result<Medium, E>,
    ) -> Result<HeterogeneousProbe<'m, 'f, 'c>, E> {
        self.validate_sample(self.ndof(), element, barycentric)?;
        let (position_m, barycentric_gradients) = if self.geometry.is_curved() {
            let point = self.geometry.geometry_point(element, barycentric)?;
            (point.position_m, point.barycentric_gradients)
        } else {
            let vertices = self.geometry.tetrahedra()[element];
            let mut anchor = 0;
            for local in 1..4 {
                if barycentric[local] > barycentric[anchor] {
                    anchor = local;
                }
            }
            let origin = self.geometry.positions[vertices[anchor] as usize];
            let mut position_m = origin;
            // Retain the affine map about a contributing vertex: it preserves
            // vertices and constant face coordinates exactly, so a rounded
            // boundary sample does not leave the original source's domain.
            for (local, (&vertex, &bary)) in vertices.iter().zip(&barycentric).enumerate() {
                if local == anchor || bary == 0.0 {
                    continue;
                }
                for (axis, position) in position_m.iter_mut().enumerate() {
                    let offset = self.geometry.positions[vertex as usize][axis] - origin[axis];
                    *position = bary.mul_add(offset, *position);
                }
            }
            (position_m, self.barycentric_gradients[element])
        };
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
                    *gradient += derivative * barycentric_gradients[a][axis];
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
            reference_barycentric: barycentric,
            position_m,
            medium,
            inverse_density,
            shapes,
        })
    }

    /// Prepare an observation at a physical point in a caller-selected element.
    ///
    /// Projected, damped Newton uses the same admitted geometric map and
    /// pointwise physical gradients as assembly. Its bounded line search stays
    /// within the closed reference simplex and only accepts a smaller residual
    /// against the unchanged physical target.
    /// The source is sampled exactly once, at the recovered position, after
    /// successful location. Outside/unresolved points and exhausted location
    /// limits make no source call. Cancellation is checked throughout location
    /// and around final preparation; no partial probe is published. No
    /// nearest-element fallback or silent change to the physical target changes
    /// the material side.
    pub fn probe_at_position<'m, E: From<AcousticFeError>>(
        &'m self,
        element: usize,
        position_m: [f64; 3],
        limits: ProbeLocationLimits,
        material_at: &mut impl FnMut(usize, [f64; 3]) -> Result<Medium, E>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<HeterogeneousProbe<'m, 'f, 'c>, E> {
        if element >= self.geometry.tetrahedra().len()
            || position_m.iter().any(|x| !x.is_finite())
            || !(1..=128).contains(&limits.max_iterations)
            || !(1..=64).contains(&limits.max_backtracks)
            || !limits.absolute_tolerance_m.is_finite()
            || limits.absolute_tolerance_m <= 0.0
        {
            return Err(invalid(
                "physical acoustic probe location or limits outside admitted range",
            )
            .into());
        }
        let checkpoint = |cancelled: &mut dyn FnMut() -> bool| {
            if cancelled() {
                Err(AcousticFeError::Cancelled)
            } else {
                Ok(())
            }
        };
        let residual = |mapped: [f64; 3]| -> ([f64; 3], f64) {
            let difference = core::array::from_fn(|axis| position_m[axis] - mapped[axis]);
            let norm = difference.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
            (difference, norm)
        };
        let mut barycentric = [0.25; 4];
        for iteration in 0..=limits.max_iterations {
            checkpoint(cancelled)?;
            let point = self.geometry.geometry_point(element, barycentric)?;
            let (difference, norm) = residual(point.position_m);
            if !norm.is_finite() {
                return Err(invalid("unrepresentable physical acoustic probe residual").into());
            }
            if norm <= limits.absolute_tolerance_m {
                checkpoint(cancelled)?;
                let probe = self.probe(element, barycentric, material_at)?;
                checkpoint(cancelled)?;
                return Ok(probe);
            }
            if iteration == limits.max_iterations {
                break;
            }
            let mut direction = point.barycentric_gradients.map(|gradient| {
                gradient
                    .iter()
                    .zip(difference)
                    .map(|(g, r)| g * r)
                    .sum::<f64>()
            });
            // Three independent reference coordinates keep the unit sum;
            // never normalize a candidate after evaluating its residual.
            direction[0] = -(direction[1] + direction[2] + direction[3]);
            if direction.iter().any(|x| !x.is_finite()) {
                return Err(invalid("unrepresentable physical acoustic probe Newton step").into());
            }
            let mut step = 1.0_f64;
            let mut accepted = false;
            for _ in 0..limits.max_backtracks {
                checkpoint(cancelled)?;
                if step <= 0.0 || !step.is_finite() {
                    break;
                }
                let mut candidate = core::array::from_fn(|i| barycentric[i] + step * direction[i]);
                candidate[0] = 1.0 - (candidate[1] + candidate[2] + candidate[3]);
                if let Some(candidate) = reference_trial(candidate) {
                    let trial = self.geometry.geometry_point(element, candidate)?;
                    let (_, trial_norm) = residual(trial.position_m);
                    if trial_norm.is_finite() && trial_norm < norm {
                        barycentric = candidate;
                        accepted = true;
                        break;
                    }
                }
                step *= 0.5;
            }
            if !accepted {
                return Err(
                    invalid("physical acoustic probe inverse map did not make progress").into(),
                );
            }
        }
        Err(invalid("physical acoustic probe exceeded its inverse-map iteration budget").into())
    }
}

impl HeterogeneousProbe<'_, '_, '_> {
    /// Original tetrahedron defining this one-sided observation.
    #[must_use]
    pub const fn element(&self) -> usize {
        self.element
    }

    /// Recovered reference point in the selected tetrahedron's vertex order.
    #[must_use]
    pub const fn reference_barycentric(&self) -> [f64; 4] {
        self.reference_barycentric
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
            || gradient
                .iter()
                .chain(&velocity)
                .chain(&intensity)
                .any(|&v| !finite(v))
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
            return Err(invalid(
                "transient observation belongs to another spatial/material operator",
            )
            .into());
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
            || gradient
                .iter()
                .chain(&velocity)
                .chain(&intensity)
                .any(|v| !v.is_finite())
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

#[cfg(test)]
mod curved_tests;
