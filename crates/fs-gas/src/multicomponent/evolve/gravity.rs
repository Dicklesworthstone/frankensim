//! Well-balanced isothermal gravity on the original evolving-species Euler step.
//!
//! The reference has uniform composition and temperature. Perturbations retain
//! every species equation and the original composition-dependent calorics.
//! Only a stationary constant-area mesh, reflecting ends and first-order space
//! are admitted. This supplies no molecular diffusive equilibrium, chemistry,
//! moving boundary, reservoir potential convention or higher-order reconstruction.
//!
//! Momentum subtracts the same original reference face tractions before
//! differencing, then adds g*(rho-rho_ref) with rho=sum(rho_species). Energy
//! work uses the SAME Riemann species mass fluxes, never a separate rho*u.
//! With phi=-g*(x-x0), internal exchanges cancel in sum V*(E+rho*phi) at both
//! original SSP-RK2 stages. E retains the native fixed species-energy gauge.

use super::*;
use crate::dynamics::IsothermalAtmosphere;
use fs_math::det;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

const CLOSED: Boundaries = Boundaries::Ends {
    left: Boundary::Reflecting,
    right: Boundary::Reflecting,
};

/// Original gas and geometry with an immutable, source-admitted atmosphere.
///
/// Reference mass fractions must sum to exactly one in the native compensated
/// reduction; they are not normalized, clipped or seeded. Zero species remain
/// exactly absent. Preparation stores physical face states and analytic
/// exponential cell averages evaluated to floating precision. Successful time
/// steps allocate no additional scratch.
#[derive(Clone, Debug)]
pub struct IsothermalMixtureGravity {
    gas: MixtureGas,
    mesh: Mesh1d,
    reference: IsothermalAtmosphere,
    fractions: Vec<f64>,
    cells: Vec<f64>,
    cell_primitive: Vec<Primitive>,
    faces: Vec<f64>,
    face_primitive: Vec<Primitive>,
    traction: Vec<f64>,
}

/// One native accepted step and its independent gravitational-energy account.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixtureGravityStep {
    /// Original Euler transaction and its species, force and energy accounts.
    pub transport: StepReport,
    /// Discrete cell-center sum phi(x_i)*M_i before the step [J], not the
    /// within-cell integral of the continuous exponential atmosphere.
    pub potential_energy_before_j: f64,
    /// The same discrete potential-energy quadrature after the step [J].
    pub potential_energy_after_j: f64,
    /// Reduced from local mass changes, preserving small exchanges [J].
    pub potential_energy_change_j: f64,
    /// Native gas-energy residual plus gravitational work plus potential change.
    /// Species retain their original formation references in the native report.
    pub total_energy_residual_j: f64,
    /// Actual stage species mass integrated against physical g [kg m/s].
    pub physical_gravity_impulse_kg_m_s: f64,
    /// Integrated reference-traction force minus g times reference cell mass.
    /// This exposes the represented reference-balance correction separately
    /// from physical gravity; it is never presented as an additional material force.
    pub reference_traction_impulse_defect_kg_m_s: f64,
}

#[derive(Clone, Copy, Default)]
pub(super) struct Balance {
    before: f64,
    after: f64,
    change: f64,
    residual: f64,
    pub(super) physical_impulse: f64,
    pub(super) reference_impulse_defect: f64,
}

impl IsothermalMixtureGravity {
    /// Prepare analytic exponential averages with the original species bank.
    ///
    /// Face areas must be identical and each volume must equal area times
    /// coordinate width to roundoff. The pressure datum and temperature use
    /// the shared isothermal-atmosphere convention, without a constant gamma.
    pub fn try_new(
        gas: MixtureGas,
        mesh: Mesh1d,
        reference: IsothermalAtmosphere,
        mass_fractions: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(&mut cancelled)?;
        let n = gas.species_count();
        let count = mesh.cell_count();
        let width = n + 2;
        let (fraction_sum, concentration) = gas.admit_densities(mass_fractions)?;
        let temperature = reference.temperature.value();
        let pressure = reference.reference_pressure.value();
        let g = reference.acceleration.value();
        let datum = reference.reference_position.value();
        let rt = R * concentration * temperature;
        if fraction_sum != 1.0
            || !g.is_finite()
            || !datum.is_finite()
            || !positive(pressure)
            || !positive(rt)
            || !(gas.domain().minimum_temperature..=gas.domain().maximum_temperature)
                .contains(&temperature)
        {
            return Err(Error::Invalid(
                "mixture atmosphere reference or normalized mass fractions",
            ));
        }
        let area = mesh.uniform_cross_section().ok_or(Error::Invalid(
            "mixture gravity requires a constant-area passage",
        ))?;
        for i in 0..count {
            poll(&mut cancelled)?;
            let volume = area * (mesh.positions()[i + 1] - mesh.positions()[i]);
            if !positive(volume)
                || (volume - mesh.volumes()[i]).abs()
                    > 64.0 * f64::EPSILON * volume.max(mesh.volumes()[i])
            {
                return Err(Error::Invalid(
                    "mixture gravity volume disagrees with physical geometry",
                ));
            }
        }
        let mut fractions = allocated(n, 0.0)?;
        fractions.copy_from_slice(mass_fractions);
        let mut cells = allocated(count.checked_mul(width).ok_or(Error::Allocation)?, 0.0)?;
        let mut faces = allocated(
            (count + 1).checked_mul(width).ok_or(Error::Allocation)?,
            0.0,
        )?;
        let mut cell_primitive = allocated(count, Primitive::default())?;
        let mut face_primitive = allocated(count + 1, Primitive::default())?;
        let mut traction = allocated(count + 1, 0.0)?;
        let beta = g / rt;
        if !beta.is_finite() {
            return Err(Error::Invalid(
                "unrepresentable mixture atmosphere scale height",
            ));
        }
        let mut densities = [0.0; MAX_SPECIES];
        let mut flux = [0.0; MAX_COMPONENTS];
        let state_at =
            |p: f64, out: &mut [f64], densities: &mut [f64]| -> Result<Primitive, Error> {
                for (rho, &fraction) in densities.iter_mut().zip(mass_fractions) {
                    *rho = (p / rt) * fraction;
                    if fraction > 0.0 && *rho == 0.0 {
                        return Err(Error::Invalid("mixture atmosphere lost a positive species"));
                    }
                }
                gas.state_at(densities, 0.0, p, out)
            };
        for i in 0..=count {
            poll(&mut cancelled)?;
            let offset = mesh.positions()[i] - datum;
            let exponent = beta * offset;
            let p = pressure * det::exp(exponent);
            if !exponent.is_finite() || !(g * offset).is_finite() || !positive(p) {
                return Err(Error::Invalid(
                    "unrepresentable mixture atmosphere pressure or potential",
                ));
            }
            let q = &mut faces[i * width..(i + 1) * width];
            let primitive = state_at(p, q, &mut densities[..n])?;
            // Cache exactly the traction used by the original boundary/face
            // kernel, including the source's pressure representation.
            if i == 0 || i == count {
                gas.wall_flux(
                    cell(q, n),
                    primitive,
                    if i == 0 { -1.0 } else { 1.0 },
                    &mut flux[..width],
                    &mut cancelled,
                )?;
            } else {
                gas.interface(
                    q,
                    primitive,
                    q,
                    primitive,
                    &mut flux[..width],
                    &mut cancelled,
                )?;
            }
            face_primitive[i] = primitive;
            traction[i] = flux[n];
        }
        for i in 0..count {
            poll(&mut cancelled)?;
            let left = mesh.positions()[i];
            let right = mesh.positions()[i + 1];
            let center = left + 0.5 * (right - left);
            let z = 0.5 * beta * (right - left);
            let average = if z.abs() < 1e-3 {
                let z2 = z * z;
                1.0 + z2 * (1.0 / 6.0 + z2 * (1.0 / 120.0 + z2 / 5040.0))
            } else {
                (det::exp(z) - det::exp(-z)) / (2.0 * z)
            };
            let p = pressure * det::exp(beta * (center - datum)) * average;
            cell_primitive[i] = state_at(
                p,
                &mut cells[i * width..(i + 1) * width],
                &mut densities[..n],
            )?;
        }
        let result = Self {
            gas,
            mesh,
            reference,
            fractions,
            cells,
            cell_primitive,
            faces,
            face_primitive,
            traction,
        };
        result.potential(&result.cells, &mut cancelled)?;
        totals(&result.cells, result.mesh.volumes(), n, &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Original species bank, caloric domain and formation references.
    pub const fn gas(&self) -> &MixtureGas {
        &self.gas
    }
    /// Original constant-area physical mesh.
    pub const fn mesh(&self) -> &Mesh1d {
        &self.mesh
    }
    /// Authored temperature, pressure datum and signed acceleration.
    pub const fn atmosphere(&self) -> IsothermalAtmosphere {
        self.reference
    }
    /// Unmodified reference composition in the source's canonical order.
    pub fn reference_mass_fractions(&self) -> &[f64] {
        &self.fractions
    }

    /// Original source-admitted cell average and its primitive diagnostics.
    pub fn reference_cell(&self, index: usize) -> Option<(Cell<'_>, Primitive)> {
        let width = self.gas.species_count() + 2;
        let values = self
            .cells
            .get(index.checked_mul(width)?..index.checked_add(1)?.checked_mul(width)?)?;
        Some((cell(values, width - 2), *self.cell_primitive.get(index)?))
    }

    /// Original source-admitted physical face and its primitive diagnostics.
    pub fn reference_face(&self, index: usize) -> Option<(Cell<'_>, Primitive)> {
        let width = self.gas.species_count() + 2;
        let values = self
            .faces
            .get(index.checked_mul(width)?..index.checked_add(1)?.checked_mul(width)?)?;
        Some((cell(values, width - 2), *self.face_primitive.get(index)?))
    }

    /// Create the native flow at exact reference rest with reflecting ends.
    pub fn initial_flow(
        &self,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MixtureEuler1d, Error> {
        poll(&mut cancelled)?;
        let mut states = allocated(self.cells.len(), 0.0)?;
        states.copy_from_slice(&self.cells);
        MixtureEuler1d::from_flat(
            self.gas.clone(),
            self.mesh.try_copy()?,
            states,
            CLOSED,
            0.0,
            &mut cancelled,
        )
    }

    /// Gravitational potential energy of an admitted flow on this exact mesh.
    pub fn potential_energy(
        &self,
        flow: &MixtureEuler1d,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64, Error> {
        self.admit(flow, &mut cancelled)?;
        self.potential(&flow.states, &mut cancelled)
    }

    fn admit(
        &self,
        flow: &MixtureEuler1d,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        poll(cancelled)?;
        if flow.gas != self.gas || flow.mesh != self.mesh || flow.boundaries != CLOSED {
            return Err(Error::Invalid(
                "mixture gravity source, mesh or reflecting ends differ",
            ));
        }
        Ok(())
    }

    pub(super) fn primitive(
        &self,
        i: usize,
        q: &[f64],
        guess: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Primitive, Error> {
        let width = self.gas.species_count() + 2;
        if q == &self.cells[i * width..(i + 1) * width] {
            Ok(self.cell_primitive[i])
        } else {
            self.gas.decode(q, Some(guess), cancelled)
        }
    }

    pub(super) fn reconstruct(
        &self,
        i: usize,
        face: usize,
        q: &[f64],
        primitive: Primitive,
        out: &mut [f64],
    ) -> Result<Primitive, Error> {
        let n = self.gas.species_count();
        let width = n + 2;
        if q == &self.cells[i * width..(i + 1) * width] {
            out.copy_from_slice(&self.faces[face * width..(face + 1) * width]);
            return Ok(self.face_primitive[face]);
        }
        // Positive relative-density reconstruction preserves the current
        // composition, including exact zeros, while recovering the same
        // reference face from either adjacent reference cell.
        let ratio = self.face_primitive[face].density / self.cell_primitive[i].density;
        let pressure = primitive.pressure
            * (self.face_primitive[face].pressure / self.cell_primitive[i].pressure);
        let mut densities = [0.0; MAX_SPECIES];
        for k in 0..n {
            densities[k] = q[k] * ratio;
            if q[k] > 0.0 && densities[k] == 0.0 {
                return Err(Error::Invalid(
                    "mixture gravity face lost a positive species",
                ));
            }
        }
        self.gas
            .state_at(&densities[..n], primitive.velocity, pressure, out)
    }

    pub(super) fn correct_forward(
        &self,
        i: usize,
        q: &[f64],
        dt: f64,
        faces: &[f64],
        out: &mut [f64],
        exchange: &mut Exchange,
    ) -> Result<(), Error> {
        let n = self.gas.species_count();
        let width = n + 2;
        let volume = self.mesh.volumes()[i];
        let area = self.mesh.face_areas()[i];
        let g = self.reference.acceleration.value();
        let rho = sum(q[..n].iter().copied());
        let perturbation_force = g * (rho - self.cell_primitive[i].density) * volume;
        let left = &faces[i * width..(i + 1) * width];
        let right = &faces[(i + 1) * width..(i + 2) * width];
        let center = self.center(i);
        let power = g
            * area
            * ((center - self.mesh.positions()[i]) * sum(left[..n].iter().copied())
                + (self.mesh.positions()[i + 1] - center) * sum(right[..n].iter().copied()));
        out[n] = q[n]
            + (dt / volume)
                * (area * (left[n] - self.traction[i]) - area * (right[n] - self.traction[i + 1])
                    + perturbation_force);
        let advected_energy = out[n + 1];
        out[n + 1] += (dt / volume) * power;
        let reference_force = area * self.traction[i + 1] - area * self.traction[i];
        let force = reference_force + perturbation_force;
        let physical_force = g * rho * volume;
        let reference_defect = reference_force - g * self.cell_primitive[i].density * volume;
        let source_energy_defect = 0.5 * ((out[n + 1] - advected_energy) * volume - dt * power);
        exchange.add([
            0.5 * dt * force,
            0.5 * dt * power,
            0.5 * dt * physical_force,
            0.5 * dt * reference_defect,
            source_energy_defect,
        ])?;
        Ok(())
    }

    fn center(&self, i: usize) -> f64 {
        self.mesh.positions()[i] + 0.5 * (self.mesh.positions()[i + 1] - self.mesh.positions()[i])
    }

    fn phi(&self, i: usize) -> f64 {
        -self.reference.acceleration.value()
            * (self.center(i) - self.reference.reference_position.value())
    }

    fn potential(
        &self,
        states: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, Error> {
        let n = self.gas.species_count();
        let mut result = Exchange::default();
        for (i, q) in states.chunks_exact(n + 2).enumerate() {
            poll(cancelled)?;
            result.add([
                sum(q[..n].iter().copied()) * self.mesh.volumes()[i] * self.phi(i),
                0.0,
                0.0,
                0.0,
                0.0,
            ])?;
        }
        Ok(result.values[0])
    }

    pub(super) fn balance(
        &self,
        initial: &[f64],
        candidate: &[f64],
        gas_residual: f64,
        source: ForceChange,
        pending: Balance,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Balance, Error> {
        let n = self.gas.species_count();
        let before = self.potential(initial, cancelled)?;
        let after = self.potential(candidate, cancelled)?;
        let mut change = Exchange::default();
        for (i, (old, new)) in initial
            .chunks_exact(n + 2)
            .zip(candidate.chunks_exact(n + 2))
            .enumerate()
        {
            poll(cancelled)?;
            let delta = sum(new[..n].iter().zip(&old[..n]).map(|(a, b)| a - b));
            change.add([
                delta * self.mesh.volumes()[i] * self.phi(i),
                0.0,
                0.0,
                0.0,
                0.0,
            ])?;
        }
        let change = change.values[0];
        let residual = sum([gas_residual, source.work, change]);
        if !residual.is_finite() {
            return Err(Error::Invalid(
                "unrepresentable gas plus gravitational energy balance",
            ));
        }
        Ok(Balance {
            before,
            after,
            change,
            residual,
            ..pending
        })
    }
}

#[derive(Default)]
pub(super) struct Exchange {
    values: [f64; 5],
    corrections: [f64; 5],
}
impl Exchange {
    fn add(&mut self, values: [f64; 5]) -> Result<(), Error> {
        for (k, value) in values.into_iter().enumerate() {
            let y = value - self.corrections[k];
            let next = self.values[k] + y;
            self.corrections[k] = (next - self.values[k]) - y;
            self.values[k] = next;
            if !value.is_finite() || !next.is_finite() || !self.corrections[k].is_finite() {
                return Err(Error::Invalid(
                    "unrepresentable mixture gravitational exchange",
                ));
            }
        }
        Ok(())
    }
    pub(super) fn physical_impulse(&self) -> f64 {
        self.values[2]
    }
    pub(super) fn reference_impulse_defect(&self) -> f64 {
        self.values[3]
    }
    pub(super) fn finish(self) -> ForceChange {
        ForceChange {
            impulse: self.values[0],
            work: self.values[1],
            energy_residual: self.values[4],
        }
    }
}

impl MixtureEuler1d {
    /// Advance actual species, momentum and energy under the prepared gravity.
    ///
    /// Requires the same source and stationary geometry, reflecting ends and
    /// explicit FirstOrder reconstruction. Ordinary stepping and held forcing
    /// are separate entry points; no persistent gravity mode is inferred.
    /// Native EOS, CFL and atomic retry/cancellation gates remain authoritative.
    pub fn step_with_isothermal_gravity(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        gravity: &IsothermalMixtureGravity,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MixtureGravityStep, Error> {
        gravity.admit(self, &mut cancelled)?;
        if controls.reconstruction != Reconstruction::FirstOrder {
            return Err(Error::Invalid(
                "mixture gravity currently requires explicit FirstOrder reconstruction",
            ));
        }
        let mut balance = Balance::default();
        let transport = if gravity.reference.acceleration.value() == 0.0 {
            let step = self.step(maximum_duration, controls, &mut cancelled)?;
            balance.residual = step.conservation_residual.energy;
            step
        } else {
            self.advance_inner_with_gravity::<false>(
                maximum_duration,
                controls,
                [None, None],
                None,
                None,
                None,
                Some((gravity, &mut balance)),
                &mut cancelled,
            )?
        };
        Ok(MixtureGravityStep {
            transport,
            potential_energy_before_j: balance.before,
            potential_energy_after_j: balance.after,
            potential_energy_change_j: balance.change,
            total_energy_residual_j: balance.residual,
            physical_gravity_impulse_kg_m_s: balance.physical_impulse,
            reference_traction_impulse_defect_kg_m_s: balance.reference_impulse_defect,
        })
    }
}
