//! Original Maxwell constitutive histories on their actual DG wall support.
//! There is one reciprocal energy-coordinate assembly for uniform and regional
//! spectra. No common relaxation fit, ghost arm or constitutive time step exists.
use super::*;
use fs_material::visco::GeneralizedMaxwell;

struct RelaxationSource {
    field_start: usize,
    dofs: usize,
    pin: WallPin,
    material: GeneralizedMaxwell,
}

fn admit_spectrum(material: &GeneralizedMaxwell) -> Result<GeneralizedMaxwell, DuctError> {
    if material.terms.len() > MAX_DG_PROPAGATOR_DOFS {
        return Err(invalid("DG wall source spectrum exceeds the bounded state capacity"));
    }
    GeneralizedMaxwell::new(material.e_inf, material.terms.clone())
        .map_err(|_| invalid("DG wall requires a passive generalized-Maxwell spectrum"))
}

impl DgDuct {
    /// Attach the existing generalized-Maxwell material spectrum to the same
    /// radial wall and gas operator. The base pin contains equilibrium stiffness
    /// K_inf; each source modulus E_j gives K_j = (K_inf/E_inf) E_j.
    ///
    /// Each radial reaction obeys z_j' = K_j v - z_j/tau_j and acts as -z_j
    /// on the wall. Its energy is integral P z_j^2/(2 K_j), and its loss is
    /// integral P z_j^2/(K_j tau_j). Energy-normalized coupling is skew;
    /// only -1/tau_j dissipates. All states advance in the original exact-ZOH
    /// exponential. Input bandwidth and dense exponential accuracy still apply.
    ///
    /// Exact zero arms are omitted; an empty/zero spectrum is bit-identical to
    /// with_compliant_wall. New memories start relaxed.
    ///
    /// # Errors
    /// Invalid spectrum, unrepresentable positive coefficients, or the existing
    /// wall/geometry/DOF refusals. No positive source arm may be dropped.
    pub fn with_compliant_wall_relaxation(
        self,
        pin: WallPin,
        material: &GeneralizedMaxwell,
    ) -> Result<Self, DuctError> {
        self.with_compliant_wall_relaxation_controlled(pin, material, &mut || false)
    }

    /// Cancellable preparation retaining every positive Maxwell source arm.
    pub fn with_compliant_wall_relaxation_controlled(
        self,
        pin: WallPin,
        material: &GeneralizedMaxwell,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, DuctError> {
        check_cancelled(cancelled)?;
        let material = admit_spectrum(material)?;
        let dofs = self.field_dofs;
        self.with_compliant_wall_controlled(pin, cancelled)?
            .append_wall_relaxation(
                &[RelaxationSource { field_start: 0, dofs, pin, material }],
                cancelled,
            )
    }

    /// Retain one independently declared Maxwell spectrum per original DG
    /// element. None supplies no hereditary arms on that element; its original
    /// elastic/Kelvin pin remains active. Source spectra may have different
    /// moduli, relaxation times and positive-arm counts.
    ///
    /// Only local nodes receive each arm's state and reaction. Every positive
    /// arm retains its own energy and decay in the same reciprocal gas/wall
    /// operator. No spectrum is averaged and no arm outside its support is
    /// allocated. Exactly homogeneous assignments take the original uniform
    /// path, including state order and operator identity.
    ///
    /// Material boundaries must coincide with existing element boundaries.
    /// Original cylindrical, frozen, local-membrane and cold-binding limits apply.
    ///
    /// # Errors
    /// Missing/extra element assignments, invalid sources, positive-coefficient
    /// underflow/overflow, cancellation or the complete-state DOF cap.
    pub fn with_compliant_wall_relaxation_regions(
        self,
        pins: &[WallPin],
        spectra: &[Option<&GeneralizedMaxwell>],
    ) -> Result<Self, DuctError> {
        self.with_compliant_wall_relaxation_regions_controlled(pins, spectra, &mut || false)
    }

    /// Cancellable preparation; a failure returns no partially coupled duct.
    pub fn with_compliant_wall_relaxation_regions_controlled(
        self,
        pins: &[WallPin],
        spectra: &[Option<&GeneralizedMaxwell>],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, DuctError> {
        check_cancelled(cancelled)?;
        if pins.is_empty() || pins.len() != self.elements.len() || spectra.len() != pins.len() {
            return Err(invalid("DG regional wall needs one pin and optional spectrum per element"));
        }
        if pins.iter().all(|pin| *pin == pins[0]) && spectra.iter().all(|law| *law == spectra[0]) {
            return match spectra[0] {
                Some(material) => self.with_compliant_wall_relaxation_controlled(
                    pins[0], material, cancelled,
                ),
                None => self.with_compliant_wall_controlled(pins[0], cancelled),
            };
        }
        let local = self.mesh.order + 1;
        let mut sources = Vec::new();
        for (element, (&pin, material)) in pins.iter().zip(spectra).enumerate() {
            check_cancelled(cancelled)?;
            if let Some(material) = material {
                sources.push(RelaxationSource {
                    field_start: element * local,
                    dofs: local,
                    pin,
                    material: admit_spectrum(material)?,
                });
            }
        }
        self.with_compliant_wall_regions_controlled(pins, cancelled)?
            .append_wall_relaxation(&sources, cancelled)
    }

    fn append_wall_relaxation(
        mut self,
        sources: &[RelaxationSource],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, DuctError> {
        let m = self.field_dofs;
        let original = self.total_dofs();
        let wall = self.compliant_wall.as_ref().expect("wall attached before its source histories");
        let velocity_start = wall.start + m;
        let beta = if wall.physical_velocity_maps.is_some() {
            None
        } else {
            let radius = self.uniform_radius_m.expect("wall requires a cylinder");
            let bulk = self.gas.density * self.gas.sound_speed * self.gas.sound_speed;
            Some(2.0 * bulk / radius)
        };
        let mut branches = Vec::new();
        let mut couplings = Vec::new();
        let mut total = original;
        for source in sources {
            let factor = source.pin.stiffness_per_area / source.material.e_inf;
            for (source_arm, &(modulus, tau)) in source.material.terms.iter().enumerate() {
                check_cancelled(cancelled)?;
                if modulus == 0.0 {
                    continue;
                }
                let start = total;
                total = total.checked_add(source.dofs)
                    .filter(|n| *n <= MAX_DG_PROPAGATOR_DOFS)
                    .ok_or_else(|| invalid("DG acoustic plus Maxwell wall state exceeds the dense DOF cap"))?;
                let stiffness = factor * modulus;
                let coupling = det::sqrt(stiffness / source.pin.surface_density);
                let rate_per_s = 1.0 / tau;
                let pressure_scale = match beta {
                    Some(beta) => det::sqrt(stiffness / beta),
                    // The geometric observation map retains the reference wall's
                    // mass even when the actual local inertia is different.
                    None => det::sqrt(stiffness * wall.pin.surface_density),
                };
                if [stiffness, coupling, rate_per_s, pressure_scale].iter()
                    .any(|value| !value.is_finite() || *value <= 0.0)
                {
                    return Err(invalid("DG Maxwell wall reduction lost a positive coefficient"));
                }
                branches.push(WallRelaxationBranch {
                    start,
                    field_start: source.field_start,
                    dofs: source.dofs,
                    source_arm,
                    stiffness_per_area: stiffness,
                    rate_per_s,
                    pressure_scale,
                });
                couplings.push(coupling);
            }
        }
        if branches.is_empty() {
            return Ok(self);
        }
        let mut operator = vec![0.0; total * total];
        for i in 0..original {
            check_cancelled(cancelled)?;
            operator[i * total..i * total + original]
                .copy_from_slice(&self.l_tilde[i * original..(i + 1) * original]);
        }
        for (branch, coupling) in branches.iter().zip(couplings) {
            for i in 0..branch.dofs {
                check_cancelled(cancelled)?;
                let velocity = velocity_start + branch.field_start + i;
                let memory = branch.start + i;
                operator[velocity * total + memory] = -coupling;
                operator[memory * total + velocity] = coupling;
                operator[memory * total + memory] = -branch.rate_per_s;
            }
        }
        check_cancelled(cancelled)?;
        for (_, port) in &mut self.ports {
            port.resize(total, 0.0);
        }
        self.compliant_wall.as_mut().expect("wall attached above").relaxation = branches;
        self.memory_dofs += total - original;
        self.l_tilde = operator;
        self.operator_identity = hash_operator(
            self.mesh, self.inlet, self.outlet, &self.elements, &self.l_tilde, &self.ports,
        );
        check_cancelled(cancelled)?;
        Ok(self)
    }
}
