//! Pressure-only outflow with the incident trace's original species calorics.
//!
//! One borrowed caloric law supplies the same shock/fan and sonic sampling used
//! by the frozen-source owner. Composition and tangential velocities cross the
//! outgoing contact unchanged; their mass, momentum and energy share one flux.

use super::{LocalCaloric, characteristic_error};
use crate::dynamics::characteristic::shared::pressure_outlet_state;
use crate::multicomponent::{
    Cell, Error, FaceFlux, MAX_NORMAL_COMPONENTS, MAX_SPATIAL_COMPONENTS, MixtureGas, NormalCell,
    NormalPrimitive, Primitive, SpatialCell, SpatialPrimitive, poll, positive,
};

impl MixtureGas {
    /// Outgoing species/normal-momentum/energy flux at a nonnegative absolute
    /// backpressure [Pa]. The original local composition supplies compression
    /// shocks, isentropic fans and sonic sampling. Departing supersonic waves
    /// leave the interior trace unchanged. A fixed pressure is not a
    /// nonreflecting boundary.
    ///
    /// The outward normal is -1 or +1; flux is oriented along the increasing
    /// coordinate. Backflow refuses because pressure alone supplies neither
    /// entering composition nor thermodynamics. Output has `species_count()+2`
    /// entries and remains unchanged on every refusal or cancellation.
    /// Successful evaluation allocates no heap workspace or new caloric owner.
    /// Zero backpressure denotes exterior vacuum; the positive-density sonic
    /// face must be supported by the incident species' original caloric laws.
    /// Exterior expansion, internal vacuum and rarefied flow are not resolved.
    pub fn pressure_outlet_flux(
        &self,
        interior: Cell<'_>,
        outward_normal: f64,
        backpressure: f64,
        out: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        let primitive = self.primitive(interior, &mut cancelled)?;
        self.pressure_outlet_admitted_flux(
            interior.partial_densities,
            primitive,
            outward_normal,
            backpressure,
            out,
            &mut cancelled,
        )
    }

    /// Planar pressure outlet carrying the incident tangential velocity and its
    /// full kinetic energy. Output has `species_count()+3` entries: species,
    /// normal momentum, transverse momentum and total referenced energy.
    /// The backflow and atomic-output policy is that of [`Self::pressure_outlet_flux`].
    pub fn pressure_outlet_normal_flux(
        &self,
        interior: NormalCell<'_>,
        outward_normal: f64,
        backpressure: f64,
        out: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        let primitive = self.normal_primitive(interior, &mut cancelled)?;
        self.pressure_outlet_components(
            interior.partial_densities,
            primitive,
            [interior.transverse_momentum / primitive.density],
            outward_normal,
            backpressure,
            out,
            &mut cancelled,
        )
    }

    /// Reconstructed-trace path to [`Self::pressure_outlet_normal_flux`], using
    /// original source and full kinetic-energy admission without caloric inversion.
    pub fn pressure_outlet_normal_flux_from_primitive(
        &self,
        interior: NormalPrimitive<'_>,
        outward_normal: f64,
        backpressure: f64,
        out: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        poll(&mut cancelled)?;
        let mut state = [0.0; MAX_NORMAL_COMPONENTS];
        let primitive =
            self.normal_conserved_into(interior, &mut state[..self.species_count() + 3])?;
        self.pressure_outlet_components(
            interior.partial_densities,
            primitive,
            [interior.transverse_velocity],
            outward_normal,
            backpressure,
            out,
            &mut cancelled,
        )
    }

    /// Three-dimensional pressure outlet carrying both incident tangential
    /// velocities and their kinetic energy. Output has `species_count()+4`
    /// entries: species, normal momentum, two transverse momenta and full
    /// referenced energy. Refusal and cancellation leave the entire output intact.
    pub fn pressure_outlet_spatial_flux(
        &self,
        interior: SpatialCell<'_>,
        outward_normal: f64,
        backpressure: f64,
        out: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        let primitive = self.spatial_primitive(interior, &mut cancelled)?;
        self.pressure_outlet_components(
            interior.partial_densities,
            primitive,
            interior.transverse_momentum.map(|m| m / primitive.density),
            outward_normal,
            backpressure,
            out,
            &mut cancelled,
        )
    }

    /// Reconstructed-trace path to [`Self::pressure_outlet_spatial_flux`],
    /// retaining full source precision without a second temperature inversion.
    pub fn pressure_outlet_spatial_flux_from_primitive(
        &self,
        interior: SpatialPrimitive<'_>,
        outward_normal: f64,
        backpressure: f64,
        out: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        poll(&mut cancelled)?;
        let mut state = [0.0; MAX_SPATIAL_COMPONENTS];
        let primitive =
            self.spatial_conserved_into(interior, &mut state[..self.species_count() + 4])?;
        self.pressure_outlet_components(
            interior.partial_densities,
            primitive,
            interior.transverse_velocity,
            outward_normal,
            backpressure,
            out,
            &mut cancelled,
        )
    }

    /// Internal reconstruction path: the caller has already source-admitted
    /// this primitive and its corresponding canonical partial densities.
    pub(in crate::multicomponent) fn pressure_outlet_admitted_flux(
        &self,
        partial_densities: &[f64],
        primitive: Primitive,
        outward_normal: f64,
        backpressure: f64,
        out: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        self.pressure_outlet_components(
            partial_densities,
            primitive,
            [],
            outward_normal,
            backpressure,
            out,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn pressure_outlet_components<const TANGENTS: usize>(
        &self,
        partial_densities: &[f64],
        inside: Primitive,
        transverse_velocity: [f64; TANGENTS],
        outward_normal: f64,
        backpressure: f64,
        out: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        poll(cancelled)?;
        let n = self.species_count();
        let energy_index = n + TANGENTS + 1;
        if TANGENTS > 2
            || out.len() != energy_index + 1
            || !matches!(outward_normal, -1.0 | 1.0)
            || !transverse_velocity.iter().all(|v| v.is_finite())
        {
            return Err(Error::Invalid(
                "mixture pressure outlet axis, normal or shear",
            ));
        }
        let law = LocalCaloric::new(self, partial_densities)?;
        let port = pressure_outlet_state(
            &law,
            Primitive {
                velocity: outward_normal * inside.velocity,
                ..inside
            },
            backpressure,
            cancelled,
        )
        .map_err(characteristic_error)?;
        self.outgoing_characteristic_flux(
            &law,
            port.primitive,
            transverse_velocity,
            outward_normal,
            port.signal_speed,
            out,
            cancelled,
        )
    }

    /// One conservative donor flux for pressure and acoustic outlet closures.
    /// The borrowed local law retains the current species fractions and their
    /// original energy references, including exact-zero absent species.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn outgoing_characteristic_flux<const TANGENTS: usize>(
        &self,
        law: &LocalCaloric<'_>,
        face: Primitive,
        transverse_velocity: [f64; TANGENTS],
        outward_normal: f64,
        signal_speed: f64,
        out: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        let n = self.species_count();
        let energy_index = n + TANGENTS + 1;
        if TANGENTS > 2
            || out.len() != energy_index + 1
            || !matches!(outward_normal, -1.0 | 1.0)
            || !transverse_velocity.iter().all(|v| v.is_finite())
        {
            return Err(Error::Invalid(
                "mixture pressure outlet axis, normal or shear",
            ));
        }
        let velocity = outward_normal * face.velocity;
        let mut candidate = [0.0; MAX_SPATIAL_COMPONENTS];
        let admitted = law.state_at(
            face.density,
            velocity,
            face.pressure,
            &mut candidate[..n + 2],
        )?;
        candidate[energy_index] = candidate[n + 1];
        for (axis, transverse) in transverse_velocity.into_iter().enumerate() {
            let momentum = admitted.density * transverse;
            candidate[n + 1 + axis] = momentum;
            candidate[energy_index] += 0.5 * momentum * transverse;
        }
        // The full conservative state must retain the original sensible energy,
        // including after transverse kinetic energy is subtracted again. A
        // finite total alone would hide a lost thermal state at extreme shear.
        if TANGENTS != 0 {
            let mut reduced_energy = candidate[energy_index];
            for momentum in &candidate[n + 1..energy_index] {
                reduced_energy -= 0.5 * momentum * (momentum / admitted.density);
            }
            let reduced = Cell {
                partial_densities: &candidate[..n],
                momentum: candidate[n],
                energy: reduced_energy,
            };
            let recovered =
                reduced_energy - 0.5 * reduced.momentum * (reduced.momentum / admitted.density);
            let (thermal, _) = self.caloric(reduced.partial_densities, admitted.temperature)?;
            let tolerance = admitted.density * self.domain().absolute_energy_tolerance
                + self.domain().relative_energy_tolerance * thermal.abs();
            if !self.admissible_cell(reduced)
                || !recovered.is_finite()
                || (recovered - thermal).abs() > tolerance
            {
                return Err(Error::Invalid(
                    "mixture pressure outlet loses declared caloric precision",
                ));
            }
        }
        let energy_flux = (candidate[energy_index] + face.pressure) * velocity;
        for value in &mut candidate[..energy_index] {
            *value *= velocity;
        }
        candidate[n] += face.pressure;
        candidate[energy_index] = energy_flux;
        let signal_speed = signal_speed.max(velocity.abs() + admitted.sound_speed);
        if !candidate[..=energy_index].iter().all(|v| v.is_finite())
            || !positive(candidate[n])
            || !positive(signal_speed)
        {
            return Err(Error::Invalid(
                "unrepresentable mixture pressure outlet flux",
            ));
        }
        poll(cancelled)?;
        out.copy_from_slice(&candidate[..=energy_index]);
        Ok(FaceFlux {
            signal_speed,
            used_hll_fallback: false,
        })
    }
}

#[cfg(test)]
mod tests;
