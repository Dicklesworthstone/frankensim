//! Original-source admission and atomic endpoint refresh, independent of heat.
use super::*;
use fs_em::plane_wave::DielectricResponse;
use fs_gas::chemistry::energy::ElectronHeavyCollisionSources;
use fs_material::ionized_gas::{IonizedGasError, IonizedGasMomentumSource};
use fs_qty::{
    QuantitySpec, Temperature,
    semantic::{QuantityKind, ValueForm},
};
use fs_thermochem::coulomb::{CoulombErrorV1, CoulombModelV1};

pub(super) struct BoundLayer<'m> {
    pub background: &'m ResolvedMaterialStatePoint,
    pub collisions: ElectronHeavyCollisionSources,
    pub electron_w: f64,
    pub heavy_w: f64,
    pub total_w: f64,
}

pub(super) fn prepare<'m>(
    wave: &PreparedIonizedGasLayerIllumination<'m>,
    input: IonizedGasLayerInput<'_>,
    area: Area,
    i: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<BoundLayer<'m>, Error> {
    let film = wave
        .layer_source(i)
        .ok_or_else(|| invalid(Some(i), "original finite film missing"))?;
    let volume = area.value() * film.thickness_m;
    if !volume.is_finite()
        || volume <= 0.0
        || (volume / input.initial.volume().value() - 1.0).abs() > 64.0 * f64::EPSILON
    {
        return Err(invalid(
            Some(i),
            "gas volume must equal positive physical area times film thickness",
        ));
    }
    let electron = input
        .initial
        .electron_index()
        .ok_or_else(|| invalid(Some(i), "structural electron source required"))?;
    if input.sources.is_empty()
        || input.sources.iter().any(|source| {
            let atom = source.source().channel();
            atom.thermodynamics() != input.initial.thermodynamics()
                || atom.electron_index() != electron
        })
        || input.optically_thin_bremsstrahlung.is_some_and(|source| {
            source.thermodynamics() != input.initial.thermodynamics()
                || source.electron_index() != electron
        })
    {
        return Err(invalid(
            Some(i),
            "atomic radiation sources must bind this exact original gas bank",
        ));
    }
    // The prepared dielectric law includes carrier friction, not direct
    // bound-free absorption of the coherent drive. Check every original edge,
    // including absent stages that this interval may create. h is exact SI.
    let drive_photon_j = (6.626_070_15e-34 / core::f64::consts::TAU) * wave.omega_rad_s();
    if input
        .sources
        .iter()
        .any(|source| drive_photon_j >= source.source().threshold().value())
    {
        return Err(invalid(
            Some(i),
            "coherent drive must lie below every original photoionization threshold",
        ));
    }
    let bound = input
        .transport_bank
        .evaluate(input.initial, &mut *cancelled)
        .map_err(|error| transport_error(i, error))?;
    let background = match film.response {
        DielectricResponse::IonizedGas(plasma) => {
            let source = plasma.transport().source();
            let actual = bound
                .as_ref()
                .ok_or_else(|| invalid(Some(i), "wave carriers cannot disappear at admission"))?;
            if source.state() != input.initial
                || actual.carriers() != source.carriers()
                || actual.limits() != source.limits()
                || actual.coulomb_model() != source.coulomb_model()
            {
                return Err(invalid(
                    Some(i),
                    "wave gas, carrier laws, Coulomb source and transport limits must match independently",
                ));
            }
            plasma.material()
        }
        DielectricResponse::Lossless(background) if bound.is_none() => background,
        _ => {
            return Err(invalid(
                Some(i),
                "all films must be original ionized gas or exactly neutral lossless gas",
            ));
        }
    };
    if background.card_identity() != input.original_background_card.content_hash() {
        return Err(invalid(
            Some(i),
            "bound-charge refresh requires the original immutable material card",
        ));
    }
    admit_background(input.initial, background, i)?;
    let id = input.initial.thermodynamics().models()[electron].species();
    let law = input
        .transport_bank
        .sources()
        .find(|(species, _)| *species == id)
        .map(|(_, law)| law)
        .ok_or_else(|| invalid(Some(i), "complete original electron momentum law required"))?;
    let IonizedGasMomentumSource::ElectronNeutral {
        cross_sections,
        controls,
        maximum_target_thermal_speed_ratio,
    } = law
    else {
        return Err(invalid(
            Some(i),
            "thermal exchange requires original target-resolved electron momentum sources",
        ));
    };
    // A neutral bank has no evaluated transport, but must retain the original
    // Coulomb policy for electrons and ionic stages born during this interval.
    let coulomb = match (bound.as_ref(), input.transport_bank.coulomb_source()) {
        (Some(transport), _) => transport.coulomb_model().cloned(),
        (None, None) => None,
        (None, Some(source)) => Some(
            CoulombModelV1::new_with_cutoff(
                input.initial,
                source.screening(),
                source.cutoff(),
                source.limits(),
                &mut *cancelled,
            )
            .map_err(|error| match error {
                CoulombErrorV1::Cancelled => Error::Cancelled,
                error => transport_error(
                    i,
                    IonizedGasRefreshError::Transport(IonizedGasError::Coulomb(error)),
                ),
            })?,
        ),
    };
    let absorbed = wave
        .layer_absorbed_power_w_m2(i)
        .ok_or_else(|| invalid(Some(i), "prepared physical film absorption missing"))?;
    let (electron_w_m2, heavy_w_m2) = match wave.layer_carrier_heating(i) {
        Some(heat) => (heat.electron_w_m2, heat.heavy_w_m2),
        None if bound.is_none() && absorbed == 0.0 => (0.0, 0.0),
        _ => return Err(invalid(Some(i), "original carrier heat partition missing")),
    };
    let powers = [electron_w_m2, heavy_w_m2, absorbed].map(|power| power * area.value());
    if powers.iter().any(|power| !power.is_finite())
        || absorbed < 0.0
        || [electron_w_m2, heavy_w_m2, absorbed]
            .iter()
            .zip(powers)
            .any(|(&before, after)| before != 0.0 && after == 0.0)
    {
        return Err(invalid(Some(i), "physical film power is unrepresentable"));
    }
    Ok(BoundLayer {
        background,
        collisions: ElectronHeavyCollisionSources {
            momentum_sources: cross_sections.clone(),
            coulomb,
            maximum_charged_mole_fraction: input
                .transport_bank
                .limits()
                .maximum_charged_mole_fraction,
            maximum_target_thermal_speed_ratio: *maximum_target_thermal_speed_ratio,
            rate_controls: *controls,
        },
        electron_w: powers[0],
        heavy_w: powers[1],
        total_w: powers[2],
    })
}

fn admit_background(
    state: &FrozenIonizedGasStateV1,
    background: &ResolvedMaterialStatePoint,
    i: usize,
) -> Result<(), Error> {
    for (axis, value) in background.query_point() {
        let value = *value;
        let (expected, kind) = match axis.as_str() {
            "T" | "T_h" => (
                state.heavy_temperature().value(),
                QuantityKind::AbsoluteTemperature,
            ),
            "T_e" => (
                state.electron_temperature().value(),
                QuantityKind::AbsoluteTemperature,
            ),
            "p" => (state.pressure().value(), QuantityKind::Pressure),
            _ => return Err(invalid(Some(i), "unmapped original background condition")),
        };
        let quantity = background.axis_quantities().get(axis).copied();
        let typed = quantity
            .and_then(QuantitySpec::semantic_type)
            .is_some_and(|ty| ty.kind() == kind && ty.form() == ValueForm::Static);
        let legacy_t = axis == "T"
            && quantity.is_none_or(|q| q == QuantitySpec::dimensional(Temperature::DIMS));
        if !(typed || legacy_t)
            || !value.is_finite()
            || value <= 0.0
            || (value != expected
                && (axis != "p" || (value / expected - 1.0).abs() > 512.0 * f64::EPSILON))
        {
            return Err(invalid(
                Some(i),
                "original background units and conditions must match the gas",
            ));
        }
    }
    Ok(())
}

pub(super) fn refresh(
    input: IonizedGasLayerInput<'_>,
    layer: &BoundLayer<'_>,
    next: &FrozenIonizedGasStateV1,
    controls: IonizedGasReactingHeatControls,
    i: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<RefreshedIonizedGasLayer, Error> {
    let old = input.initial;
    let heavy_k = (next.heavy_temperature().value() - old.heavy_temperature().value()).abs();
    let electron_k = if old
        .electron_index()
        .is_some_and(|e| old.amounts()[e].value() > 0.0)
    {
        (next.electron_temperature().value() - old.electron_temperature().value()).abs()
    } else {
        0.0
    };
    if heavy_k.max(electron_k) > controls.maximum_temperature_increment_k {
        return Err(Error::TemperatureIncrement {
            layer: i,
            heavy_k,
            electron_k,
            maximum_k: controls.maximum_temperature_increment_k,
        });
    }
    for (s, (new, old_amount)) in next.amounts().iter().zip(old.amounts()).enumerate() {
        poll(cancelled)?;
        let change_mol = (new.value() - old_amount.value()).abs();
        let maximum_mol = controls.maximum_absolute_population_increment_mol
            + controls.maximum_relative_population_increment * old_amount.value();
        if !maximum_mol.is_finite() {
            return Err(invalid(
                Some(i),
                "species increment allowance is unrepresentable",
            ));
        }
        if change_mol > maximum_mol {
            return Err(Error::PopulationIncrement {
                layer: i,
                species: old.thermodynamics().models()[s].species().clone(),
                change_mol,
                maximum_mol,
            });
        }
    }
    let transport = input
        .transport_bank
        .evaluate(next, &mut *cancelled)
        .map_err(|error| transport_error(i, error))?;
    let mut query = layer.background.query();
    for (axis, _) in layer.background.query_point() {
        let value = match axis.as_str() {
            "T" | "T_h" => next.heavy_temperature().value(),
            "T_e" => next.electron_temperature().value(),
            "p" => next.pressure().value(),
            _ => return Err(invalid(Some(i), "unmapped original background condition")),
        };
        query = query
            .with(axis.clone(), value)
            .map_err(|_| invalid(Some(i), "refreshed background coordinate is inadmissible"))?;
    }
    let background = layer
        .background
        .at_query_point(input.original_background_card, &query)
        .map_err(|error| Error::Background { layer: i, error })?;
    Ok(RefreshedIonizedGasLayer {
        transport,
        background,
    })
}
