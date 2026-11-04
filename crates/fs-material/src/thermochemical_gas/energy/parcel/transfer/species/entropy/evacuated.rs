//! Closed inventory and entropy change when a rigid empty volume first fills.
//!
//! The original donor supplies every conserved species and its caloric/entropy
//! reference. Vacuum supplies no state or temperature. The existing zero-amount
//! entropy branch evaluates the receiver's new inventory, including its mixing
//! and volume contribution; no logarithm of zero is evaluated.

use super::{
    Accumulator, ClosedGasChanges, ClosedGasParcel, Error, Reference, Temperature, amount, check,
    contributions, same_source, source_error,
};

fn final_amount(parcel: &ClosedGasParcel, reference: &Reference<'_>) -> Result<f64, Error> {
    let components = parcel.model().thermodynamics().components();
    match components.binary_search_by(|model| model.species().cmp(reference.model.species())) {
        Ok(index) => amount(parcel, index),
        Err(_) => Ok(0.0),
    }
}

/// Actual closed change from one finite donor to that donor and a newly filled
/// fixed-volume receiver. Every source must come from a positive original donor
/// inventory. Source constants cancel before floating-point evaluation.
///
/// The receiver's original mass, energy and entropy are exactly zero. Its final
/// volume is the authored empty vessel's volume; the physical owner must retain
/// that geometry when reconstructing the new gas. The initial temperature used
/// as the caloric reference belongs ONLY to the original donor.
///
/// Each pass is bounded by the original gas owner's species limit. The returned
/// entropy may be negative: the enclosing transport owner must apply its own
/// absolute entropy budget and refuse an inadmissible time step. This endpoint
/// account does not establish an unconditional time-discretization theorem.
pub(in crate::thermochemical_gas::energy::parcel::transfer::species) fn newly_filled_pair_change(
    original: &ClosedGasParcel,
    donor: (u64, &ClosedGasParcel),
    receiver: (u64, &ClosedGasParcel),
    mass_tolerance_kg: f64,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ClosedGasChanges, Error> {
    check(&mut cancelled)?;
    if !mass_tolerance_kg.is_finite()
        || mass_tolerance_kg <= 0.0
        || donor.0 == receiver.0
        || donor.1.volume_m3() != original.volume_m3()
    {
        return Err(Error::Invalid(
            "newly filled pair requires distinct identities, fixed donor volume and positive mass tolerance",
        ));
    }

    let before = original.model().thermodynamics().components();
    // Validate the complete endpoints first, including any original species
    // that was removed from the donor by a selective transfer. A source cannot
    // enter a closed pair merely because its name or coefficients are similar.
    for (_, parcel) in [donor, receiver] {
        for (index, model) in parcel
            .model()
            .thermodynamics()
            .components()
            .iter()
            .enumerate()
        {
            check(&mut cancelled)?;
            let original_index = before
                .binary_search_by(|source| source.species().cmp(model.species()))
                .map_err(|_| Error::Invalid("newly filled pair introduced an unknown species"))?;
            if !same_source(&before[original_index], model) {
                return Err(Error::Invalid(
                    "newly filled pair changed an original species source",
                ));
            }
            if amount(original, original_index)? == 0.0 && amount(parcel, index)? != 0.0 {
                return Err(Error::Invalid(
                    "newly filled pair created an absent original inventory",
                ));
            }
        }
    }

    let mut mass = Accumulator::default();
    mass.add(donor.1.mass_kg())?;
    mass.add(receiver.1.mass_kg())?;
    mass.add(-original.mass_kg())?;
    let total_mass_residual_kg = mass.total()?;
    if total_mass_residual_kg.abs() > mass_tolerance_kg {
        return Err(Error::Invalid(
            "newly filled pair did not conserve total gas mass",
        ));
    }

    let initial_temperature = original.state().state().temperature;
    let mut entropy = Accumulator::default();
    let mut energy = Accumulator::default();
    let mut maximum_species_mass_residual_kg = 0.0_f64;
    for (index, model) in before.iter().enumerate() {
        check(&mut cancelled)?;
        let n0 = amount(original, index)?;
        if n0 == 0.0 {
            continue;
        }
        let mut reference = Reference {
            model,
            temperature_k: initial_temperature,
            volume_m3: original.volume_m3(),
            amount_mol: n0,
            mass_change: Accumulator::default(),
        };
        let remaining = final_amount(donor.1, &reference)?;
        let received = final_amount(receiver.1, &reference)?;
        for (endpoint, old_amount, next_amount) in
            [(donor, n0, remaining), (receiver, 0.0, received)]
        {
            check(&mut cancelled)?;
            if old_amount == 0.0 && next_amount == 0.0 {
                // An absent species has no temperature or partial pressure.
                continue;
            }
            let temperature = endpoint.1.state().state().temperature;
            for term in contributions(
                &reference,
                old_amount,
                next_amount,
                initial_temperature,
                temperature,
                endpoint.1.volume_m3(),
                endpoint.0,
            )? {
                entropy.add(term)?;
            }
            if next_amount > 0.0 {
                let du = model
                    .molar_internal_energy_change(
                        Temperature::new(initial_temperature),
                        Temperature::new(temperature),
                    )
                    .map_err(|error| Error::Parcel {
                        id: endpoint.0,
                        error: source_error(error),
                    })?
                    .value();
                energy.add(next_amount * du)?;
            }
            reference
                .mass_change
                .add((next_amount - old_amount) * model.molar_mass().value())?;
        }
        let residual = reference.mass_change.total()?;
        maximum_species_mass_residual_kg = maximum_species_mass_residual_kg.max(residual.abs());
        if residual.abs() > mass_tolerance_kg {
            return Err(Error::Invalid(
                "newly filled pair did not conserve each source-qualified species",
            ));
        }
    }
    check(&mut cancelled)?;
    Ok(ClosedGasChanges {
        total_mass_residual_kg,
        maximum_species_mass_residual_kg,
        energy_change_j: energy.total()?,
        entropy_change_j_k: entropy.total()?,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::thermochemical_gas::{
        ThermochemicalGasModel,
        tests::{model, point, species_with_energy_reference, transport_in_domain},
    };
    use fs_thermochem::{
        Composition, CompositionBasis, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
        mixture::FrozenIdealGasMixtureModelV1,
    };

    const VOLUME: f64 = 0.001;
    const MASS_TOLERANCE: f64 = 1e-13;
    const MOLAR_MASS: f64 = 0.028;

    fn pure(name: &str, cp_r: f64, datum: f64) -> ThermochemicalGasModel {
        ThermochemicalGasModel::new(
            FrozenIdealGasMixtureModelV1::new(
                vec![species_with_energy_reference(
                    name, MOLAR_MASS, cp_r, 0.0, datum,
                )],
                Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
            )
            .unwrap(),
            vec![transport_in_domain(
                name,
                1.7e-5,
                0.025,
                [200.0, 2000.0],
                [1.0, 1e7],
            )],
            None,
        )
        .unwrap()
    }

    fn parcel(model: ThermochemicalGasModel, t: f64, p: f64, volume: f64) -> ClosedGasParcel {
        ClosedGasParcel::new(model, point(t, p), volume).unwrap()
    }

    #[test]
    fn g1_empty_receiver_recovers_ideal_mixture_free_expansion_entropy() {
        let source = model(vec![0.75, 0.25], 0.001);
        let original = parcel(source.clone(), 400.0, 100000.0, VOLUME);
        // One original volume expands into three total volumes at unchanged
        // temperature and composition. The empty receiver owns two volumes.
        let donor = parcel(source.clone(), 400.0, 100000.0 / 3.0, VOLUME);
        let receiver = parcel(source, 400.0, 100000.0 / 3.0, 2.0 * VOLUME);
        let change = newly_filled_pair_change(
            &original,
            (5, &donor),
            (8, &receiver),
            MASS_TOLERANCE,
            || false,
        )
        .unwrap();
        let molar_mass = original
            .state()
            .thermodynamics()
            .properties()
            .molar_mass()
            .value();
        let expected = original.mass_kg() / molar_mass * R * 3.0_f64.ln();
        assert!((change.entropy_change_j_k - expected).abs() < 1e-12);
        assert_eq!(change.energy_change_j, 0.0);
        assert!(change.total_mass_residual_kg.abs() <= MASS_TOLERANCE);
        assert!(change.maximum_species_mass_residual_kg <= MASS_TOLERANCE);
    }

    #[test]
    fn g3_energy_conserving_first_fill_exposes_a_negative_entropy_overshoot() {
        let fraction = 0.1;
        let cv_r = 1.5;
        let temperature_ratio = 1.0 + 0.5 / cv_r;
        let donor_t = 400.0 * (1.0 - fraction * temperature_ratio) / (1.0 - fraction);
        let receiver_t = 400.0 * temperature_ratio;
        let receiver_volume = VOLUME * 0.01;
        let mut entropy_changes = Vec::new();
        for datum in [0.0, 1e12] {
            let source = pure("A", cv_r + 1.0, datum);
            let original = parcel(source.clone(), 400.0, 100000.0, VOLUME);
            let donor = parcel(
                source.clone(),
                donor_t,
                100000.0 * (1.0 - fraction * temperature_ratio),
                VOLUME,
            );
            let receiver = parcel(
                source,
                receiver_t,
                100000.0 * fraction * temperature_ratio * (VOLUME / receiver_volume),
                receiver_volume,
            );
            let change = newly_filled_pair_change(
                &original,
                (5, &donor),
                (8, &receiver),
                MASS_TOLERANCE,
                || false,
            )
            .unwrap();
            // Independent ideal-gas entropy at fixed volumes: n Cv ln(T/T0)
            // minus n R ln((n/V)/(n0/V0)). This oversized molecular deposit
            // conserves energy but cannot be accepted as a passive closed step.
            let expected = original.mass_kg() / MOLAR_MASS
                * R
                * ((1.0 - fraction) * (cv_r * (donor_t / 400.0).ln() - (1.0 - fraction).ln())
                    + fraction
                        * (cv_r * temperature_ratio.ln()
                            - (fraction * VOLUME / receiver_volume).ln()));
            assert!((change.entropy_change_j_k - expected).abs() < 1e-12);
            assert!(change.entropy_change_j_k < -1e-3);
            assert!(change.energy_change_j.abs() < 1e-10);
            assert!(change.total_mass_residual_kg.abs() <= MASS_TOLERANCE);
            entropy_changes.push(change.entropy_change_j_k);
        }
        assert_eq!(entropy_changes[0], entropy_changes[1]);
    }

    #[test]
    fn g4_empty_receiver_audit_refuses_source_geometry_mass_and_cancellation_errors() {
        let source = pure("A", 2.5, 0.0);
        let original = parcel(source.clone(), 400.0, 100000.0, VOLUME);
        let donor = parcel(source.clone(), 400.0, 50000.0, VOLUME);
        let receiver = donor.clone();
        for model in [pure("B", 2.5, 0.0), pure("A", 3.5, 0.0)] {
            let changed = parcel(model, 400.0, 50000.0, VOLUME);
            assert!(matches!(
                newly_filled_pair_change(
                    &original,
                    (5, &donor),
                    (8, &changed),
                    MASS_TOLERANCE,
                    || false
                ),
                Err(Error::Invalid(_))
            ));
        }
        let changed_volume = parcel(source.clone(), 400.0, 25000.0, 2.0 * VOLUME);
        assert!(matches!(
            newly_filled_pair_change(
                &original,
                (5, &changed_volume),
                (8, &receiver),
                MASS_TOLERANCE,
                || false
            ),
            Err(Error::Invalid(_))
        ));
        let added_mass = parcel(source, 400.0, 55000.0, VOLUME);
        assert!(matches!(
            newly_filled_pair_change(
                &original,
                (5, &donor),
                (8, &added_mass),
                MASS_TOLERANCE,
                || false
            ),
            Err(Error::Invalid(_))
        ));
        let mut polls = 0;
        assert!(matches!(
            newly_filled_pair_change(
                &original,
                (5, &donor),
                (8, &receiver),
                MASS_TOLERANCE,
                || {
                    polls += 1;
                    polls == 5
                }
            ),
            Err(Error::Cancelled)
        ));
        assert_eq!(polls, 5);
    }
}
