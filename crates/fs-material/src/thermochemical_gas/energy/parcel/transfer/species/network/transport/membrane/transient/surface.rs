//! Ordered source-backed exterior gas/polymer transfer resistance.
//!
//! The existing interface permeance multiplies gas-equivalent pressure drop.
//! This is an explicitly authored independent-species surface law, not an
//! inferred mixed-gas film correlation or an additional sorbed storage model.

use super::*;
use crate::thermochemical_gas::ThermochemicalGasModel;

type Surfaces<'a> = [Option<InterfaceBinding<'a>>; 2];
type Coefficients = [Option<MolarPermeance>; 2];

pub(super) fn resolve<'a>(
    species: &SpeciesId,
    temperature_k: f64,
    layers: &[MembraneDiffusionLayerSource<'a>],
    sources: &[Option<MembraneInterfaceMaterial<'a>>; 2],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Surfaces<'a>, Coefficients), Error> {
    if sources.iter().all(Option::is_none) {
        return Ok(([None, None], [None, None]));
    }
    let requirement = ScalarPropertyRequirement::try_new(
        MOLAR_INTERFACE_PERMEANCE_PROPERTY,
        MOLAR_PERMEANCE_DIMS,
        ScalarAdmissibility::NonNegative,
    )
    .map_err(source_error)?;
    let mut bound = [None, None];
    let mut coefficients = [None, None];
    for side in 0..2 {
        poll(cancelled)?;
        let Some(binding) = &sources[side] else {
            continue;
        };
        let (polymer, adjacent) = if side == 0 {
            (
                &binding.interface.surface_b().material,
                layers[0].material.material.id(),
            )
        } else {
            (
                &binding.interface.surface_a().material,
                layers[layers.len() - 1].material.material.id(),
            )
        };
        if &binding.species != species || polymer != adjacent {
            return Err(Error::Invalid(
                "surface source must bind the same species and original ordered adjacent polymer",
            ));
        }
        let point = selection_point(binding.context, temperature_k)?;
        let original = resolve_interface_state_point(
            binding.interface,
            &point,
            &[requirement.clone()],
            binding.selection.clone(),
        )
        .map_err(source_error)?;
        let property = original
            .property(MOLAR_INTERFACE_PERMEANCE_PROPERTY)
            .ok_or(Error::Invalid(
                "resolved molar surface permeance is missing",
            ))?;
        admit_henry(property)?;
        coefficients[side] =
            Some(MolarPermeance::try_new(property.value_si()).map_err(GasEnergyError::from)?);
        bound[side] = Some(InterfaceBinding {
            binding: binding.clone(),
            original,
        });
    }
    Ok((bound, coefficients))
}

impl Bindings<'_> {
    /// Source identity only. The interface claim must itself qualify the actual
    /// independent-species mixture; pure-gas measurements gain no new domain.
    pub(super) fn admit_gas_surfaces(&self, catalog: &ThermochemicalGasModel) -> Result<(), Error> {
        for side in 0..2 {
            let Some(surface) = &self.surfaces[side] else {
                continue;
            };
            let exterior = if side == 0 {
                &surface.binding.interface.surface_a().material
            } else {
                &surface.binding.interface.surface_b().material
            };
            let material = catalog
                .transport_materials()
                .iter()
                .find(|law| law.species == self.species)
                .ok_or(Error::Invalid(
                    "surface species has no original gas transport material",
                ))?;
            if exterior != material.material.id() {
                return Err(Error::Invalid(
                    "surface exterior must identify the original catalog species transport material",
                ));
            }
        }
        Ok(())
    }

    pub(super) fn admit_surfaces<const N: usize>(
        &self,
        samples: [(&MembraneDiffusionState, &[f64], &[[Option<f64>; 2]]); N],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[Option<[ResolvedInterfaceStatePoint; 2]>; 2], Error> {
        let mut result = [None, None];
        for side in 0..2 {
            poll(cancelled)?;
            let Some(binding) = &self.surfaces[side] else {
                continue;
            };
            let (mut lo, mut hi) = (f64::INFINITY, 0.0_f64);
            for (state, _, faces) in samples {
                let index = if side == 0 {
                    0
                } else {
                    state.model().cell_count()
                };
                for trace in faces[index] {
                    let pressure = trace.ok_or(Error::Invalid(
                        "coincident surface and bulk seals leave a trace unknown; cannot qualify the surface source"))?;
                    lo = lo.min(pressure);
                    hi = hi.max(pressure);
                }
            }
            let low = binding
                .original
                .at_query_point(
                    binding.binding.interface,
                    &updated_point(binding.binding.context, self.temperature_k, lo)?,
                )
                .map_err(source_error)?;
            poll(cancelled)?;
            let high = if lo == hi {
                low.clone()
            } else {
                binding
                    .original
                    .at_query_point(
                        binding.binding.interface,
                        &updated_point(binding.binding.context, self.temperature_k, hi)?,
                    )
                    .map_err(source_error)?
            };
            result[side] = Some([low, high]);
        }
        Ok(result)
    }
}
