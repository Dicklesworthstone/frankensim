//! Finite non-isothermal absorbing layers coupled through two diffuse ports.
//!
//! The conduction slab owner preserves direction through the ENTIRE stack;
//! intermediate material interfaces do not redistribute radiation. Exterior
//! port irradiation and outgoing power use the enclosure's explicit diffuse
//! closure. This does not resolve angular illumination or claim Lambertian
//! emergence from an uncoated nonscattering slab.
//!
//! The original slab owns its matched-index or sourced Fresnel-interface
//! response, including trapped internal rays. Its no-scattering and source
//! support assumptions remain in force. Port net flux is not another
//! thermal source: only the separate layer absorption-minus-emission heats
//! the material. Layer temperatures need not be equal.

use fs_conduction::radiation::slab::AbsorbingSlab;

use super::sheets::{BoundaryCoefficients, coefficients};
use super::{
    Layer, LayerAbsorption, MAX_RADIATION_SURFACES, ParticipatingEnclosure, Power, RadiationError,
    RadiationReport, RadiationSolve, RadiationWorkspace, SurfaceEmissivity, Temperature, TwoZone,
    checkpoint, finite,
};

/// Bound on ALL material layers in one enclosure, independent of surface count.
pub const MAX_SLAB_LAYERS: usize = 64;

/// Equal-area exterior ports of one original prepared absorbing stack.
#[derive(Clone, Debug, PartialEq)]
pub struct SlabPair {
    /// Original enclosure face indices: first-layer side, then last-layer side.
    pub surfaces: [usize; 2],
    /// Frozen, source-bearing optical depth and angular-response owner.
    pub optics: AbsorbingSlab,
}

pub(super) fn admit_axis(
    n: usize,
    slabs: &[SlabPair],
) -> Result<([bool; MAX_RADIATION_SURFACES], usize), RadiationError> {
    if n > MAX_RADIATION_SURFACES || slabs.len() > n / 2 {
        return Err(RadiationError::Invalid("invalid slab surface axis"));
    }
    let mut ports = [false; MAX_RADIATION_SURFACES];
    let mut layers = 0;
    for slab in slabs {
        let [a, b] = slab.surfaces;
        if a >= n || b >= n || a == b || ports[a] || ports[b] {
            return Err(RadiationError::Invalid(
                "slab ports must be distinct, in range and unshared",
            ));
        }
        let count = slab.optics.layer_count();
        if count == 0 || count > MAX_SLAB_LAYERS - layers {
            return Err(RadiationError::Invalid(
                "at most 64 total absorbing slab layers are admitted",
            ));
        }
        layers += count;
        ports[a] = true;
        ports[b] = true;
    }
    Ok((ports, layers))
}

pub(super) fn layer_powers(slabs: &[SlabPair], values: &[f64]) -> Vec<Vec<Power>> {
    let mut offset = 0;
    slabs
        .iter()
        .map(|slab| {
            let end = offset + slab.optics.layer_count();
            let result = values[offset..end]
                .iter()
                .copied()
                .map(Power::new)
                .collect();
            offset = end;
            result
        })
        .collect()
}

impl ParticipatingEnclosure {
    /// Evaluate original gray gas/opaque boundaries and non-isothermal layers.
    ///
    /// Gray cards cover ONLY non-port faces in increasing enclosure order.
    /// Layer temperatures follow slab-then-layer order. Each port's temperature
    /// must equal its adjacent first/last layer, but the two ports may differ.
    /// No heat is booked to ports in addition to the material layer powers.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_with_slabs(
        &self,
        surface_temperatures: &[Temperature],
        emissivities: &[SurfaceEmissivity],
        slabs: &[SlabPair],
        layer_temperatures: &[Vec<Temperature>],
        gas_temperatures: [Temperature; 2],
        absorption: &[LayerAbsorption; 2],
        controls: RadiationSolve,
        work: &mut RadiationWorkspace,
        cancelled: impl FnMut() -> bool,
    ) -> Result<RadiationReport, RadiationError> {
        if slabs.is_empty() {
            return Err(RadiationError::Invalid(
                "provide at least one absorbing slab",
            ));
        }
        self.evaluate_gray(
            surface_temperatures,
            emissivities,
            &[],
            slabs,
            layer_temperatures,
            gas_temperatures,
            absorption,
            controls,
            work,
            cancelled,
        )
    }

    /// Use original gas temperatures and refuse paths through an absent layer.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_room_with_slabs(
        &self,
        room: &TwoZone,
        surface_temperatures: &[Temperature],
        emissivities: &[SurfaceEmissivity],
        slabs: &[SlabPair],
        layer_temperatures: &[Vec<Temperature>],
        absorption: &[LayerAbsorption; 2],
        controls: RadiationSolve,
        work: &mut RadiationWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RadiationReport, RadiationError> {
        for path in &self.paths {
            checkpoint(&mut cancelled)?;
            for segment in &path.segments {
                if segment.length.value() > 0.0 && room.layer(segment.layer).volume().value() == 0.0
                {
                    return Err(RadiationError::Invalid(
                        "optical path crosses an absent gas layer",
                    ));
                }
            }
        }
        self.evaluate_with_slabs(
            surface_temperatures,
            emissivities,
            slabs,
            layer_temperatures,
            [
                room.layer(Layer::Lower).temperature(),
                room.layer(Layer::Upper).temperature(),
            ],
            absorption,
            controls,
            work,
            cancelled,
        )
    }

    pub(super) fn admit_slab_state(
        &self,
        surface_temperatures: &[Temperature],
        slabs: &[SlabPair],
        layer_temperatures: &[Vec<Temperature>],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[bool; MAX_RADIATION_SURFACES], RadiationError> {
        let (ports, _) = admit_axis(self.geometry.len(), slabs)?;
        if layer_temperatures.len() != slabs.len() {
            return Err(RadiationError::Invalid(
                "slab temperature and source axes differ",
            ));
        }
        for (slab, temperatures) in slabs.iter().zip(layer_temperatures) {
            checkpoint(cancelled)?;
            let [a, b] = slab.surfaces;
            let count = slab.optics.layer_count();
            if temperatures.len() != count
                || self.geometry.areas_m2()[a] != self.geometry.areas_m2()[b]
            {
                return Err(RadiationError::Invalid(
                    "slab ports need equal area and one temperature per layer",
                ));
            }
            let mut values = [0.0; MAX_SLAB_LAYERS];
            for (value, t) in values.iter_mut().zip(temperatures) {
                *value = t.value();
            }
            slab.optics
                .validate_temperatures(&values[..count])
                .map_err(|_| {
                    RadiationError::Invalid("slab temperature is outside original optical support")
                })?;
            if surface_temperatures[a] != temperatures[0]
                || surface_temperatures[b] != temperatures[count - 1]
            {
                return Err(RadiationError::Invalid(
                    "slab port temperatures must match their adjacent layers",
                ));
            }
        }
        Ok(ports)
    }

    pub(super) fn gray_slab_emissivities(
        &self,
        temperatures: &[Temperature],
        emissivities: &[SurfaceEmissivity],
        slabs: &[SlabPair],
        layer_temperatures: &[Vec<Temperature>],
        values: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), RadiationError> {
        let ports = self.admit_slab_state(temperatures, slabs, layer_temperatures, cancelled)?;
        if emissivities.len() != self.geometry.len() - 2 * slabs.len() {
            return Err(RadiationError::Invalid(
                "supply one original gray source for each non-port face",
            ));
        }
        let mut sources = emissivities.iter();
        for (i, value) in values.iter_mut().enumerate() {
            *value = 0.0;
            if !ports[i] {
                let source = sources
                    .next()
                    .ok_or(RadiationError::Invalid("missing opaque source"))?;
                if source
                    .temperature_validity_k()
                    .is_some_and(|(lo, hi)| !(lo..=hi).contains(&temperatures[i].value()))
                {
                    return Err(RadiationError::Invalid(
                        "surface temperature is outside original emissivity validity",
                    ));
                }
                *value = source.value();
            }
        }
        Ok(())
    }

    pub(super) fn prepare_slab_optics(
        &self,
        emissivities: &[f64],
        slabs: &[SlabPair],
        absorption: [f64; 2],
        work: &mut RadiationWorkspace,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), RadiationError> {
        self.prepare_exchange(absorption, work, cancelled)?;
        self.factor_boundary_optics(&slab_coefficients(emissivities, slabs), work, cancelled)
    }

    pub(super) fn solve_slab_emission(
        &self,
        emissivities: &[f64],
        slabs: &[SlabPair],
        controls: RadiationSolve,
        work: &mut RadiationWorkspace,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(f64, f64), RadiationError> {
        let n = self.geometry.len();
        let coefficients = slab_coefficients(emissivities, slabs);
        let mut source = [0.0; MAX_RADIATION_SURFACES];
        for i in 0..n {
            source[i] = finite(emissivities[i] * work.emission[i])?;
        }
        let mut offset = 0;
        for slab in slabs {
            checkpoint(cancelled)?;
            let end = offset + slab.optics.layer_count();
            let emission = slab
                .optics
                .emission_to_faces(&work.layer_emission[offset..end])
                .map_err(|_| RadiationError::Invalid("original slab emission response refused"))?;
            source[slab.surfaces[0]] = emission[0];
            source[slab.surfaces[1]] = emission[1];
            offset = end;
        }
        let irradiation =
            self.solve_boundary_field(&coefficients, &source[..n], work, cancelled)?;
        let mut maximum = 0.0_f64;
        for i in 0..n {
            let residual = if coefficients.partner[i] == i {
                emissivities[i] * (work.radiosity[i] - work.emission[i])
                    - coefficients.reflection[i] * work.powers[i] / self.geometry.areas_m2()[i]
            } else {
                work.radiosity[i]
                    - source[i]
                    - coefficients.reflection[i] * irradiation[i]
                    - coefficients.transmission[i] * irradiation[coefficients.partner[i]]
            };
            maximum = maximum.max(finite(residual)?.abs());
        }
        offset = 0;
        for slab in slabs {
            checkpoint(cancelled)?;
            let [a, b] = slab.surfaces;
            let end = offset + slab.optics.layer_count();
            slab.optics
                .layer_net_fluxes(
                    [irradiation[a], irradiation[b]],
                    &work.layer_emission[offset..end],
                    &mut work.layer_powers[offset..end],
                )
                .map_err(|_| {
                    RadiationError::Invalid("original slab layer heat response refused")
                })?;
            for power in &mut work.layer_powers[offset..end] {
                *power = finite(self.geometry.areas_m2()[a] * *power)?;
            }
            // Exterior ports are exchange surfaces, not extra heat inventories.
            work.powers[a] = 0.0;
            work.powers[b] = 0.0;
            offset = end;
        }
        let balance = finite(
            work.powers
                .iter()
                .chain(&work.layer_powers[..offset])
                .copied()
                .sum(),
        )?;
        if maximum > controls.surface_residual_w_m2 || balance.abs() > controls.energy_residual_w {
            return Err(RadiationError::Residual {
                surface_w_m2: maximum,
                energy_w: balance,
            });
        }
        checkpoint(cancelled)?;
        Ok((maximum, balance))
    }
}

fn slab_coefficients(emissivities: &[f64], slabs: &[SlabPair]) -> BoundaryCoefficients {
    let mut result = coefficients(emissivities, &[]);
    for slab in slabs {
        let [a, b] = slab.surfaces;
        for (side, (i, j)) in [(a, b), (b, a)].into_iter().enumerate() {
            result.reflection[i] = slab.optics.reflectance();
            result.transmission[i] = slab.optics.transmittance();
            result.partner[i] = j;
            result.diagonal[i] = 1.0;
            result.thermal_anchor[i] = slab.optics.face_coupling().iter().any(|x| x[side] > 0.0);
        }
    }
    result
}

/// Implicit layer-resolved radiation heating of original thermal cells.
pub mod thermal;

#[cfg(test)]
mod tests;
