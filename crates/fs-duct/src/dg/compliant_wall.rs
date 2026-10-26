//! Reciprocal local wall mechanics in the original DG time operator.
//!
//! With circumference P, area A and gas bulk modulus B, the equations are
//! (A/B) p_t + q_x = -P v, sigma v_t = p-K w-r v, and w_t=v.
//! Uniform coefficients make all three mass matrices proportional to the same
//! exact polynomial Gram matrix. In energy coordinates the pressure/velocity
//! coupling is g=sqrt(P B/(A sigma)); the wall oscillator frequency is
//! sqrt(K/sigma). Both couplings are skew. Only -(r/sigma) v_tilde dissipates.
//! This is the same local WallPin law as the frequency-domain duct, including
//! radial inertia, rather than a retuned wave speed or a separate time step.

use super::*;
use fs_phs::WallPin;

#[derive(Debug, Clone)]
pub(super) struct CompliantWall {
    pin: WallPin,
    start: usize,
    displacement_scale: f64,
    velocity_scale: f64,
}

/// Physical wall fields on the original DG nodes, read without advancing time.
#[derive(Debug, Clone, PartialEq)]
pub struct DgCompliantWallState {
    /// Outward radial displacement [m], relative to the reference cylinder.
    pub displacement_m: Vec<f64>,
    /// Outward radial velocity [m/s].
    pub velocity_m_s: Vec<f64>,
    /// Integral of P (sigma v² + K w²)/2 [J], excluding gas storage.
    pub stored_energy_j: f64,
    /// Integral of P r v² [W], excluding viscothermal and terminal losses.
    pub dissipated_power_w: f64,
}

fn invalid(what: &'static str) -> DuctError {
    DuctError::BadParameter { what }
}

impl DgDuct {
    /// Exact common radius from the original cylindrical geometry, if present.
    /// A sampled taper or profile is never inferred to be a uniform cylinder.
    #[must_use]
    pub const fn uniform_circular_radius_m(&self) -> Option<f64> {
        self.uniform_radius_m
    }

    /// Attach physical radial wall inertia, stiffness and viscosity before
    /// constructing a state or propagator. The same exact DG exponential then
    /// advances gas, wall and any previously attached viscothermal memories.
    /// Every pressure port sees the wall reaction and includes wall storage.
    ///
    /// The wall covers the full lateral surface of identical circular segments
    /// in a uniform gas. Its coefficients are the existing WallPin law; a
    /// material adapter must derive them from geometry and actual source data.
    /// This is a frozen linear local membrane, without axial shell waves,
    /// prestress, external wall radiation or an independently moving mesh.
    /// Matched ends retain the original GAS characteristic impedance; they are
    /// a resistive termination, not a reflectionless elastic-tube extension.
    ///
    /// Attach viscothermal memory FIRST if required. Wall mechanics are then
    /// appended without replacing any gas-memory coordinate. Initial acoustic
    /// states start with the wall at rest; reusing an old state is refused by
    /// the original operator identity. Source validity and input bandwidth
    /// remain the caller's responsibility, as for the constant WallPin owner.
    ///
    /// # Errors
    /// Invalid coefficients, duplicate attachment, nonuniform gas or radius,
    /// or total state beyond the existing dense-propagator cap.
    pub fn with_compliant_wall(mut self, pin: WallPin) -> Result<Self, DuctError> {
        if self.compliant_wall.is_some() {
            return Err(invalid("DG compliant wall is already attached"));
        }
        let radius = self
            .uniform_radius_m
            .ok_or_else(|| invalid("DG compliant wall requires identical circular cylinders"))?;
        if self
            .gas_field
            .as_ref()
            .is_some_and(|field| field.uniform_state().is_none())
        {
            return Err(invalid("DG compliant wall requires uniform material gas"));
        }
        fs_phs::wall_specific_impedance(&pin, 1.0).map_err(|_| {
            invalid("DG compliant wall requires passive finite mass, stiffness and resistance")
        })?;
        let bulk = self.gas.density * self.gas.sound_speed * self.gas.sound_speed;
        // P/A=2/r. No geometric averaging or nodal mass lumping is used.
        let pressure_per_displacement = 2.0 * bulk / radius;
        let coupling = det::sqrt(pressure_per_displacement / pin.surface_density);
        let omega = det::sqrt(pin.stiffness_per_area / pin.surface_density);
        let decay = pin.resistance / pin.surface_density;
        let displacement_scale =
            1.0 / det::sqrt(pressure_per_displacement * pin.stiffness_per_area);
        let velocity_scale = 1.0 / det::sqrt(pressure_per_displacement * pin.surface_density);
        if [
            bulk,
            pressure_per_displacement,
            coupling,
            omega,
            displacement_scale,
            velocity_scale,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.0)
            || !decay.is_finite()
            || (pin.resistance > 0.0 && decay == 0.0)
        {
            return Err(invalid(
                "DG compliant wall coefficients are unrepresentable",
            ));
        }
        let m = self.field_dofs;
        let original = self.total_dofs();
        let total = original
            .checked_add(2 * m)
            .filter(|n| *n <= MAX_DG_PROPAGATOR_DOFS)
            .ok_or_else(|| {
                invalid("DG acoustic plus compliant wall state exceeds the dense DOF cap")
            })?;
        let mut operator = vec![0.0; total * total];
        for i in 0..original {
            operator[i * total..i * total + original]
                .copy_from_slice(&self.l_tilde[i * original..(i + 1) * original]);
        }
        for i in 0..m {
            let displacement = original + i;
            let velocity = original + m + i;
            operator[i * total + velocity] = -coupling;
            operator[velocity * total + i] = coupling;
            operator[displacement * total + velocity] = omega;
            operator[velocity * total + displacement] = -omega;
            operator[velocity * total + velocity] = -decay;
        }
        for (_, port) in &mut self.ports {
            port.resize(total, 0.0);
        }
        self.compliant_wall = Some(CompliantWall {
            pin,
            start: original,
            displacement_scale,
            velocity_scale,
        });
        self.memory_dofs += 2 * m;
        self.l_tilde = operator;
        self.g_block = None;
        self.operator_identity = hash_operator(
            self.mesh,
            self.inlet,
            self.outlet,
            &self.elements,
            &self.l_tilde,
            &self.ports,
        );
        Ok(self)
    }

    /// The original physical wall coefficients, if attached.
    #[must_use]
    pub fn compliant_wall(&self) -> Option<WallPin> {
        self.compliant_wall.as_ref().map(|wall| wall.pin)
    }

    /// Read actual wall motion and storage from the native acoustic state.
    /// This is an observation only; no pressure is replayed through a second
    /// oscillator. Nodal ordering is the original DG node_positions ordering.
    ///
    /// # Errors
    /// Missing wall, foreign/nonfinite state, or unrepresentable observation.
    pub fn compliant_wall_state(&self, state: &DgState) -> Result<DgCompliantWallState, DuctError> {
        let wall = self
            .compliant_wall
            .as_ref()
            .ok_or_else(|| invalid("DG duct has no compliant wall"))?;
        if state.operator_identity != self.operator_identity
            || state.tilde.len() != self.total_dofs()
            || state.tilde.iter().any(|value| !value.is_finite())
        {
            return Err(invalid(
                "DG wall observation requires a finite state of this exact operator",
            ));
        }
        let m = self.field_dofs;
        let local = self.mesh.order + 1;
        let mut displacement_m = vec![0.0; m];
        let mut velocity_m_s = vec![0.0; m];
        for (e, element) in self.elements.iter().enumerate() {
            let base = e * local;
            for i in 0..local {
                for j in 0..local {
                    let map = element.mp_inv_sqrt[i * local + j];
                    displacement_m[base + i] += map * state.tilde[wall.start + base + j];
                    velocity_m_s[base + i] += map * state.tilde[wall.start + m + base + j];
                }
                displacement_m[base + i] *= wall.displacement_scale;
                velocity_m_s[base + i] *= wall.velocity_scale;
            }
        }
        let stored_energy_j = 0.5
            * state.tilde[wall.start..wall.start + 2 * m]
                .iter()
                .map(|value| value * value)
                .sum::<f64>();
        let dissipated_power_w = self.compliant_wall_loss(state);
        if displacement_m
            .iter()
            .chain(&velocity_m_s)
            .any(|value| !value.is_finite())
            || !stored_energy_j.is_finite()
            || !dissipated_power_w.is_finite()
        {
            return Err(invalid("DG wall observation is unrepresentable"));
        }
        Ok(DgCompliantWallState {
            displacement_m,
            velocity_m_s,
            stored_energy_j,
            dissipated_power_w,
        })
    }

    pub(super) fn compliant_wall_loss(&self, state: &DgState) -> f64 {
        self.compliant_wall.as_ref().map_or(0.0, |wall| {
            wall.pin.resistance / wall.pin.surface_density
                * state.tilde[wall.start + self.field_dofs..wall.start + 2 * self.field_dofs]
                    .iter()
                    .map(|value| value * value)
                    .sum::<f64>()
        })
    }
}

#[cfg(test)]
mod tests;
