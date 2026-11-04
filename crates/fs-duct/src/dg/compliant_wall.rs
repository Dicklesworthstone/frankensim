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
    relaxation: Vec<WallRelaxationBranch>,
}

#[derive(Debug, Clone)]
struct WallRelaxationBranch {
    start: usize,
    rate_per_s: f64,
    pressure_scale: f64,
}

/// Physical wall fields on the original DG nodes, read without advancing time.
#[derive(Debug, Clone, PartialEq)]
pub struct DgCompliantWallState {
    /// Outward radial displacement [m], relative to the reference cylinder.
    pub displacement_m: Vec<f64>,
    /// Outward radial velocity [m/s].
    pub velocity_m_s: Vec<f64>,
    /// Each retained Maxwell arm's radial reaction pressure [Pa] at DG nodes.
    /// Positive pressure opposes outward wall motion. Exact zero arms are absent.
    pub relaxation_pressure_pa: Vec<Vec<f64>>,
    /// Complete wall energy [J]: inertia, equilibrium spring and Maxwell springs.
    /// Gas storage is excluded; relaxation_stored_energy_j is already included.
    pub stored_energy_j: f64,
    /// Energy [J] in the retained Maxwell springs alone.
    pub relaxation_stored_energy_j: f64,
    /// Kelvin and Maxwell dissipation [W], excluding gas and terminal losses.
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
            relaxation: Vec::new(),
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

    /// Attach the existing generalized-Maxwell material spectrum to the same
    /// radial wall and gas operator. The base pin contains equilibrium stiffness
    /// K_inf; each source modulus E_j gives K_j = (K_inf/E_inf) E_j. No fitted
    /// damping or independently stepped constitutive history is introduced.
    ///
    /// Each radial reaction pressure obeys z_j' = K_j v - z_j/tau_j and acts
    /// as -z_j on the wall. Its energy is integral P z_j^2/(2 K_j), and its
    /// loss is integral P z_j^2/(K_j tau_j). Energy-normalized coupling to
    /// wall velocity is skew; only -1/tau_j dissipates. All states advance in
    /// the original exact-ZOH exponential, so midpoint-specific clock limits
    /// do not apply. Input bandwidth and dense exponential accuracy still do.
    ///
    /// The source's exact zero arms are omitted; an empty/zero spectrum is
    /// bit-identical to with_compliant_wall. New memories start relaxed.
    ///
    /// # Errors
    /// Invalid source spectrum, unrepresentable positive coefficients, or the
    /// original wall/geometry/DOF refusals. No positive arm may be dropped.
    pub fn with_compliant_wall_relaxation(
        self,
        pin: WallPin,
        material: &fs_material::visco::GeneralizedMaxwell,
    ) -> Result<Self, DuctError> {
        if material.terms.len() > MAX_DG_PROPAGATOR_DOFS {
            return Err(invalid(
                "DG wall source spectrum exceeds the bounded state capacity",
            ));
        }
        let material =
            fs_material::visco::GeneralizedMaxwell::new(material.e_inf, material.terms.clone())
                .map_err(|_| invalid("DG wall requires a passive generalized-Maxwell spectrum"))?;
        let count = material
            .terms
            .iter()
            .filter(|(modulus, _)| *modulus > 0.0)
            .count();
        let mut duct = self.with_compliant_wall(pin)?;
        if count == 0 {
            return Ok(duct);
        }
        let m = duct.field_dofs;
        let original = duct.total_dofs();
        let extra = count
            .checked_mul(m)
            .ok_or_else(|| invalid("DG wall memory DOF overflow"))?;
        let total = original
            .checked_add(extra)
            .filter(|n| *n <= MAX_DG_PROPAGATOR_DOFS)
            .ok_or_else(|| {
                invalid("DG acoustic plus Maxwell wall state exceeds the dense DOF cap")
            })?;
        let wall = duct
            .compliant_wall
            .as_ref()
            .expect("wall was attached above");
        let velocity_start = wall.start + m;
        let radius = duct.uniform_radius_m.expect("wall requires a cylinder");
        let bulk = duct.gas.density * duct.gas.sound_speed * duct.gas.sound_speed;
        let beta = 2.0 * bulk / radius;
        let factor = pin.stiffness_per_area / material.e_inf;
        let mut branches = Vec::with_capacity(count);
        let mut couplings = Vec::with_capacity(count);
        for (modulus, tau) in material.terms {
            if modulus == 0.0 {
                continue;
            }
            let stiffness = factor * modulus;
            let coupling = det::sqrt(stiffness / pin.surface_density);
            let rate_per_s = 1.0 / tau;
            let pressure_scale = det::sqrt(stiffness / beta);
            if [stiffness, coupling, rate_per_s, pressure_scale]
                .iter()
                .any(|value| !value.is_finite() || *value <= 0.0)
            {
                return Err(invalid(
                    "DG Maxwell wall reduction lost a positive coefficient",
                ));
            }
            branches.push(WallRelaxationBranch {
                start: original + branches.len() * m,
                rate_per_s,
                pressure_scale,
            });
            couplings.push(coupling);
        }
        let mut operator = vec![0.0; total * total];
        for i in 0..original {
            operator[i * total..i * total + original]
                .copy_from_slice(&duct.l_tilde[i * original..(i + 1) * original]);
        }
        for (branch, coupling) in branches.iter().zip(couplings) {
            for i in 0..m {
                let velocity = velocity_start + i;
                let memory = branch.start + i;
                operator[velocity * total + memory] = -coupling;
                operator[memory * total + velocity] = coupling;
                operator[memory * total + memory] = -branch.rate_per_s;
            }
        }
        for (_, port) in &mut duct.ports {
            port.resize(total, 0.0);
        }
        duct.compliant_wall
            .as_mut()
            .expect("wall was attached above")
            .relaxation = branches;
        duct.memory_dofs += extra;
        duct.l_tilde = operator;
        duct.operator_identity = hash_operator(
            duct.mesh,
            duct.inlet,
            duct.outlet,
            &duct.elements,
            &duct.l_tilde,
            &duct.ports,
        );
        Ok(duct)
    }

    /// The original equilibrium/inertial/Kelvin coefficients, if attached.
    /// Maxwell branch storage and pressure are returned by compliant_wall_state.
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
        let mut relaxation_pressure_pa = Vec::with_capacity(wall.relaxation.len());
        let mut relaxation_stored_energy_j = 0.0;
        for branch in &wall.relaxation {
            let mut pressure = vec![0.0; m];
            for (e, element) in self.elements.iter().enumerate() {
                let base = e * local;
                for i in 0..local {
                    for j in 0..local {
                        pressure[base + i] += element.mp_inv_sqrt[i * local + j]
                            * state.tilde[branch.start + base + j];
                    }
                    pressure[base + i] *= branch.pressure_scale;
                }
            }
            relaxation_stored_energy_j += 0.5
                * state.tilde[branch.start..branch.start + m]
                    .iter()
                    .map(|value| value * value)
                    .sum::<f64>();
            relaxation_pressure_pa.push(pressure);
        }
        let stored_energy_j = 0.5
            * state.tilde[wall.start..wall.start + 2 * m]
                .iter()
                .map(|value| value * value)
                .sum::<f64>()
            + relaxation_stored_energy_j;
        let dissipated_power_w = self.compliant_wall_loss(state);
        if displacement_m
            .iter()
            .chain(&velocity_m_s)
            .chain(relaxation_pressure_pa.iter().flatten())
            .any(|value| !value.is_finite())
            || !stored_energy_j.is_finite()
            || !dissipated_power_w.is_finite()
        {
            return Err(invalid("DG wall observation is unrepresentable"));
        }
        Ok(DgCompliantWallState {
            displacement_m,
            velocity_m_s,
            relaxation_pressure_pa,
            stored_energy_j,
            relaxation_stored_energy_j,
            dissipated_power_w,
        })
    }

    pub(super) fn compliant_wall_loss(&self, state: &DgState) -> f64 {
        self.compliant_wall.as_ref().map_or(0.0, |wall| {
            let mut loss = wall.pin.resistance / wall.pin.surface_density
                * state.tilde[wall.start + self.field_dofs..wall.start + 2 * self.field_dofs]
                    .iter()
                    .map(|value| value * value)
                    .sum::<f64>();
            for branch in &wall.relaxation {
                loss += branch.rate_per_s
                    * state.tilde[branch.start..branch.start + self.field_dofs]
                        .iter()
                        .map(|value| value * value)
                        .sum::<f64>();
            }
            loss
        })
    }
}

#[cfg(test)]
mod tests;
