//! Opt-in energy correction on admitted pressure/velocity-equilibrium contacts.
//!
//! Species, represented momenta and physical face fluxes retain their original
//! owners. Only energy is reconstructed, with the original mixture calorics.
//! The resulting numerical energy defect is explicit; it is neither heat nor
//! mechanical work and does not establish an entropy or shock-capturing theorem.

use super::{Cell, Error, MixtureGas, allocated, poll, positive};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

/// Authored contact admission and whole-attempt numerical energy budget.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContactEquilibriumControl {
    /// Relative pressure tolerance and velocity tolerance in sound-speed units.
    /// Must be finite, positive and at most 1e-6; bulk speed never sets the scale.
    pub relative_tolerance: f64,
    /// Maximum sum of absolute RK-weighted energy corrections [J].
    /// Finite and nonnegative. Opposing corrections cannot cancel this budget.
    pub maximum_energy_defect_j: f64,
}
impl ContactEquilibriumControl {
    /// Admit explicit tolerances without changing a workspace or field.
    pub fn validate(self) -> Result<(), Error> {
        if !positive(self.relative_tolerance)
            || self.relative_tolerance > 1e-6
            || !self.maximum_energy_defect_j.is_finite()
            || self.maximum_energy_defect_j < 0.0
        {
            return Err(Error::Invalid("contact-equilibrium controls"));
        }
        Ok(())
    }
}

/// Actual represented energy changes from accepted correction maps.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ContactCorrection {
    /// Signed RK-weighted numerical energy change [J], not physical heating.
    pub signed_energy_j: f64,
    /// Sum of absolute RK-weighted energy changes [J].
    pub absolute_energy_j: f64,
    /// Unique cells with a nonzero represented correction in this attempt.
    pub corrected_cells: usize,
}

/// The native transport report and its separately exposed energy defect.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContactEquilibriumStepReport<R> {
    /// Native field, boundary, clock and conservative-balance report.
    pub transport: R,
    /// Numerical correction, separate from physical boundary/source work.
    pub correction: ContactCorrection,
    /// Native energy residual minus the signed numerical correction [J].
    pub corrected_energy_residual_j: f64,
}

#[derive(Clone, Copy, Default)]
struct Reference {
    pressure: f64,
    velocity: [f64; 3],
    sound_speed: f64,
    eligible: bool,
    corrected: bool,
}
impl Reference {
    fn matches(self, pressure: f64, velocity: [f64; 3], tolerance: f64) -> bool {
        (pressure - self.pressure).abs() <= tolerance * pressure.max(self.pressure)
            && velocity
                .into_iter()
                .zip(self.velocity)
                .all(|(a, b)| (a - b).abs() <= tolerance * self.sound_speed)
    }
}

/// Caller-owned contact references and compensated energy accounting.
///
/// Preparation allocates one flat reference array. Attempts and projections
/// allocate nothing. Owners must inspect every actual incident face, excluding
/// unsupported physical boundaries, and publish the receipt only with the
/// complete accepted field. Call `begin_attempt` before each halved retry.
pub struct ContactWorkspace {
    references: Vec<Reference>,
    control: Option<ContactEquilibriumControl>,
    correction: ContactCorrection,
    signed_roundoff: f64,
    absolute_roundoff: f64,
}
impl ContactWorkspace {
    /// Required numerical heap storage; fixed-size owner fields are excluded.
    pub fn required_workspace_bytes(cells: usize) -> Result<usize, Error> {
        if cells == 0 {
            return Err(Error::Invalid("empty contact-equilibrium workspace"));
        }
        cells
            .checked_mul(size_of::<Reference>())
            .ok_or(Error::Invalid(
                "contact-equilibrium workspace size overflow",
            ))
    }

    /// Reserve the reusable cell references within the caller's byte budget.
    pub fn new(
        cells: usize,
        maximum_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(&mut cancelled)?;
        if Self::required_workspace_bytes(cells)? > maximum_bytes {
            return Err(Error::Invalid("contact-equilibrium workspace byte budget"));
        }
        let references = allocated(cells, Reference::default())?;
        if references.capacity() * size_of::<Reference>() > maximum_bytes {
            return Err(Error::Invalid("contact-equilibrium workspace byte budget"));
        }
        poll(&mut cancelled)?;
        Ok(Self {
            references,
            control: None,
            correction: ContactCorrection::default(),
            signed_roundoff: 0.0,
            absolute_roundoff: 0.0,
        })
    }

    /// Actual retained numerical heap capacity in bytes.
    pub fn workspace_bytes(&self) -> usize {
        self.references.capacity() * size_of::<Reference>()
    }

    /// Check the original cell count and authored controls without mutation.
    pub fn admit(&self, cells: usize, control: ContactEquilibriumControl) -> Result<(), Error> {
        control.validate()?;
        if cells != self.references.len() {
            return Err(Error::Invalid("contact-equilibrium field dimensions"));
        }
        Ok(())
    }

    /// Reset private scratch and accounting for one complete RK attempt.
    pub fn begin_attempt(&mut self, control: ContactEquilibriumControl) -> Result<(), Error> {
        control.validate()?;
        self.references.fill(Reference::default());
        self.control = Some(control);
        self.correction = ContactCorrection::default();
        self.signed_roundoff = 0.0;
        self.absolute_roundoff = 0.0;
        Ok(())
    }

    /// Retain an original cell's admitted pressure and full velocity. Unused
    /// velocity dimensions are zero. Later stages intersect eligibility with
    /// this original reference; they never replace it or re-enable a cell.
    pub fn set_reference(
        &mut self,
        index: usize,
        pressure: f64,
        velocity: [f64; 3],
        sound_speed: f64,
        first_stage: bool,
    ) -> Result<(), Error> {
        let control = self.control()?;
        admit_trace(pressure, velocity)?;
        if !positive(sound_speed) {
            return Err(Error::Invalid("contact-equilibrium sound speed"));
        }
        let reference = self
            .references
            .get_mut(index)
            .ok_or(Error::Invalid("contact-equilibrium cell index"))?;
        if first_stage {
            *reference = Reference {
                pressure,
                velocity,
                sound_speed,
                eligible: true,
                corrected: false,
            };
        } else {
            if !positive(reference.pressure) {
                return Err(Error::Invalid(
                    "contact-equilibrium reference not initialized",
                ));
            }
            reference.eligible &= reference.matches(pressure, velocity, control.relative_tolerance);
        }
        Ok(())
    }

    /// Compare BOTH reconstructed face traces with the recipient cell's
    /// reference. Matching face traces alone do not exclude a pressure gradient.
    /// The conservative HLL fallback disqualifies this recipient.
    #[allow(clippy::too_many_arguments)]
    pub fn observe_pair(
        &mut self,
        index: usize,
        pressure_left: f64,
        velocity_left: [f64; 3],
        pressure_right: f64,
        velocity_right: [f64; 3],
        hll_fallback: bool,
    ) -> Result<(), Error> {
        let control = self.control()?;
        let reference = self.reference(index)?;
        admit_trace(pressure_left, velocity_left)?;
        admit_trace(pressure_right, velocity_right)?;
        self.references[index].eligible &= !hll_fallback
            && reference.matches(pressure_left, velocity_left, control.relative_tolerance)
            && reference.matches(pressure_right, velocity_right, control.relative_tolerance);
        Ok(())
    }

    /// Exclude a physical boundary or other unsupported stencil for this attempt.
    pub fn exclude(&mut self, index: usize) -> Result<(), Error> {
        self.reference(index)?;
        self.references[index].eligible = false;
        Ok(())
    }

    /// Reconstruct only energy from the original calorics and represented
    /// momenta. Axial/planar/spatial momenta have lengths 1/2/3. Excluded cells
    /// retain `energy` exactly. Source refusal, cancellation or an exceeded
    /// budget changes neither the receipt nor the contact mask.
    ///
    /// Supply weight 1/2 for each forward-Euler projection and 1 for the final
    /// SSP average. Volume and weight must be finite and positive; weight is
    /// at most one. Only actual represented changes enter the energy budget.
    #[allow(clippy::too_many_arguments)]
    pub fn project_energy(
        &mut self,
        gas: &MixtureGas,
        index: usize,
        partial_densities: &[f64],
        momenta: &[f64],
        energy: f64,
        volume: f64,
        weight: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64, Error> {
        poll(&mut cancelled)?;
        let control = self.control()?;
        let reference = self.reference(index)?;
        if partial_densities.len() != gas.species_count()
            || !(1..=3).contains(&momenta.len())
            || !momenta.iter().all(|x| x.is_finite())
            || !energy.is_finite()
            || !positive(volume)
            || !positive(weight)
            || weight > 1.0
        {
            return Err(Error::Invalid(
                "contact-equilibrium state, volume or RK weight",
            ));
        }
        if !reference.eligible {
            poll(&mut cancelled)?;
            return Ok(energy);
        }
        let (density, concentration) = gas.admit_densities(partial_densities)?;
        let mut velocity = [0.0; 3];
        let mut kinetic = [0.0; 3];
        for (axis, &momentum) in momenta.iter().enumerate() {
            velocity[axis] = momentum / density;
            kinetic[axis] = 0.5 * momentum * velocity[axis];
        }
        if !velocity.into_iter().chain(kinetic).all(f64::is_finite) {
            return Err(Error::Invalid(
                "unrepresentable contact-equilibrium kinetic energy",
            ));
        }
        if !reference.matches(reference.pressure, velocity, control.relative_tolerance) {
            poll(&mut cancelled)?;
            self.references[index].eligible = false;
            return Ok(energy);
        }
        let temperature = (reference.pressure / concentration) / R;
        let (thermal, _) = gas.caloric(partial_densities, temperature)?;
        let mut projected = thermal;
        for value in &kinetic[..momenta.len()] {
            projected += value;
        }
        // Match the existing normal/spatial reduction order and retain the
        // original precision gate after every kinetic component is removed.
        let mut reduced = projected;
        for value in &kinetic[1..momenta.len()] {
            reduced -= value;
        }
        let recovered = reduced - kinetic[0];
        let tolerance = density * gas.domain().absolute_energy_tolerance
            + gas.domain().relative_energy_tolerance * thermal.abs();
        if !gas.admissible_cell(Cell {
            partial_densities,
            momentum: momenta[0],
            energy: reduced,
        }) || !tolerance.is_finite()
            || !recovered.is_finite()
            || (recovered - thermal).abs() > tolerance
        {
            return Err(Error::Invalid(
                "contact-equilibrium energy loses source precision",
            ));
        }
        let delta = projected - energy;
        if delta == 0.0 {
            poll(&mut cancelled)?;
            return Ok(energy);
        }
        let extensive = (delta * weight) * volume;
        if !extensive.is_finite() || (delta != 0.0 && extensive == 0.0) {
            return Err(Error::Invalid(
                "unrepresentable contact-equilibrium energy defect",
            ));
        }
        let (signed, signed_roundoff) = add(
            self.correction.signed_energy_j,
            self.signed_roundoff,
            extensive,
        );
        let (absolute, absolute_roundoff) = add(
            self.correction.absolute_energy_j,
            self.absolute_roundoff,
            extensive.abs(),
        );
        if ![signed, signed_roundoff, absolute, absolute_roundoff]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(Error::Invalid(
                "contact-equilibrium energy accounting overflow",
            ));
        }
        if absolute > control.maximum_energy_defect_j {
            return Err(Error::ContactEnergyBudget {
                observed_j: absolute,
                maximum_j: control.maximum_energy_defect_j,
            });
        }
        poll(&mut cancelled)?;
        self.correction.signed_energy_j = signed;
        self.correction.absolute_energy_j = absolute;
        self.signed_roundoff = signed_roundoff;
        self.absolute_roundoff = absolute_roundoff;
        if delta != 0.0 && !reference.corrected {
            self.correction.corrected_cells += 1;
            self.references[index].corrected = true;
        }
        Ok(projected)
    }

    /// Current attempt's accumulated corrections; owners publish it atomically.
    pub const fn receipt(&self) -> ContactCorrection {
        self.correction
    }

    fn control(&self) -> Result<ContactEquilibriumControl, Error> {
        self.control.ok_or(Error::Invalid(
            "contact-equilibrium attempt not initialized",
        ))
    }
    fn reference(&self, index: usize) -> Result<Reference, Error> {
        self.references
            .get(index)
            .copied()
            .filter(|r| positive(r.pressure))
            .ok_or(Error::Invalid(
                "contact-equilibrium reference not initialized",
            ))
    }
}

fn admit_trace(pressure: f64, velocity: [f64; 3]) -> Result<(), Error> {
    if !positive(pressure) || !velocity.into_iter().all(f64::is_finite) {
        return Err(Error::Invalid("contact-equilibrium pressure or velocity"));
    }
    Ok(())
}
fn add(total: f64, correction: f64, value: f64) -> (f64, f64) {
    let y = value - correction;
    let next = total + y;
    (next, (next - total) - y)
}

#[cfg(test)]
mod tests;
