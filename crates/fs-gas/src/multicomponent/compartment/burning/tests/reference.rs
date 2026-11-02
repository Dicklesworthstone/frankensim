//! Conserved elemental-reference shifts must change neither power nor heating.
use super::*;

fn shifted_mechanism(heat: f64, cp_shift: f64, reference: f64) -> IdealGasKineticsV1 {
    let original = fixture_mechanism(heat, cp_shift, false);
    let formation = [0.0, -0.5 * heat, -0.25 * heat, 0.0, 0.0];
    // Molecular mass is this fixture's exact conserved linear combination of
    // C/H/N/O counts. Shift conserved elemental references, not species freely.
    let models = (0..5)
        .map(|i| {
            source(
                NAMES[i],
                MASSES[i],
                3.5 + if i == 1 { cp_shift } else { 0.0 },
                0.02,
                MASSES[i] * reference + formation[i],
            )
        })
        .collect();
    let bank = IdealGasEquilibriumV1::new(
        original.thermodynamics().elemental_matrix().clone(),
        original.thermodynamics().charge_vector().clone(),
        models,
    )
    .unwrap();
    IdealGasKineticsV1::new(
        bank,
        "shifted-global-reaction-reference",
        1,
        original.reactions().to_vec(),
        0.0,
    )
    .unwrap()
}

#[test]
fn g3_power_conversion_and_actual_chemical_heat_preserve_elemental_reference_invariance() {
    for heat in [40.0, 80.0] {
        for reference in [0.0, 2.0_f64.powi(45), -2.0_f64.powi(45)] {
            let mechanism = shifted_mechanism(heat, 1.5, reference);
            let zone = room(&mechanism, FEED, true);
            let owner = adapter(&zone, &mechanism);
            for temperature in [2.0, 8.0] {
                let delta_h = R * (1.5 * temperature - heat);
                let requested = owner
                    .extent_for_power(
                        Power::new(20.0),
                        Time::new(0.3),
                        Temperature::new(temperature),
                        || false,
                    )
                    .unwrap();
                close(requested.value(), -6.0 / delta_h, 1e-16);
                let report = owner
                    .evaluate(
                        &zone,
                        Layer::Upper,
                        requested,
                        Time::new(0.3),
                        None,
                        Temperature::new(temperature),
                        || false,
                    )
                    .unwrap();
                close(report.reaction_enthalpy.value(), delta_h, 2e-12);
                close(report.actual_extent.value(), requested.value(), 1e-16);
                close(report.released_power.value(), 20.0, 2e-13);
                // The stored chemical correction belongs to u(Tmin), not the
                // independently selected enthalpy-report temperature.
                let reference_u = R * (1.5 * zone.gas().domain().minimum_temperature - heat);
                close(
                    report.chemical_rates.referenced_energy_w,
                    -report.actual_extent.value() * reference_u / 0.3,
                    2e-12,
                );
                assert!(report.chemical_reference_roundoff_w.is_finite());
            }
        }
    }
}

#[test]
fn g1_reference_shifted_burning_matches_actual_population_first_law_and_pressure_work() {
    for heat in [40.0, 80.0] {
        let mut baseline: Option<(f64, f64)> = None;
        for reference in [0.0, 2.0_f64.powi(45), -2.0_f64.powi(45)] {
            let mechanism = shifted_mechanism(heat, 0.0, reference);
            let mut zone = room(&mechanism, FEED, true);
            let owner = adapter(&zone, &mechanism);
            let before = zone.layer(Layer::Upper).species_masses().to_vec();
            let before_n = sum(before.iter().zip(MASSES).map(|(m, mass)| m / mass));
            let result = owner
                .apply_split(
                    &mut zone,
                    Layer::Upper,
                    Amount::new(0.1),
                    Time::new(0.3),
                    Temperature::new(8.0),
                    ZoneSolve::default(),
                    || false,
                )
                .unwrap();
            let after = zone.layer(Layer::Upper).species_masses();
            let after_n = sum(after.iter().zip(MASSES).map(|(m, mass)| m / mass));
            // Products were initially absent. This oracle uses the actual
            // published populations, independently of the requested extent and
            // stored caloric residual. Phi(T)=u_sensible/R plus its zero offset.
            let reacted = 0.5 * after[1] / MASSES[1] + 0.25 * after[2] / MASSES[2];
            let phi = |t: f64| 2.5 * t + 0.01 * t * t;
            let c = (before_n / after_n - 1.0) * phi(2.0) + heat * reacted / after_n;
            let expected = 2.0 + 2.0 * c / (2.54 + (2.54_f64.powi(2) + 0.04 * c).sqrt());
            let temperature = zone.layer(Layer::Upper).temperature().value();
            close(temperature, expected, 2e-9);
            close(
                R * (after_n * phi(temperature) - before_n * phi(2.0) - heat * reacted),
                0.0,
                2e-8,
            );
            close(
                zone.pressure().value(),
                after_n * R * temperature / zone.volume().value(),
                2e-9,
            );
            assert_eq!(zone.time().value(), 0.0); // Split source cannot advance transport time twice.
            assert_eq!(result.relaxation.interface_work_j, [0.0, -0.0]);
            assert!(result.burning.chemical_reference_roundoff_w.is_finite());
            if let Some((t, p)) = baseline {
                close(temperature, t, 2e-9);
                close(zone.pressure().value(), p, 2e-9);
            } else {
                baseline = Some((temperature, zone.pressure().value()));
            }
        }
    }
}

#[test]
fn g3_moving_interface_burning_preserves_its_original_mechanical_work_under_reference_shifts() {
    let mut baseline: Option<[f64; 4]> = None;
    for reference in [0.0, 2.0_f64.powi(45), -2.0_f64.powi(45)] {
        let mechanism = shifted_mechanism(40.0, 0.0, reference);
        let mut zone = room(&mechanism, FEED, false);
        let owner = adapter(&zone, &mechanism);
        let before = zone.layers;
        let old_pressure = zone.pressure().value();
        let report = owner
            .apply_split(
                &mut zone,
                Layer::Upper,
                Amount::new(0.1),
                Time::new(0.3),
                Temperature::new(2.0),
                ZoneSolve::default(),
                || false,
            )
            .unwrap();
        let values = [
            zone.layer(Layer::Lower).temperature().value(),
            zone.layer(Layer::Upper).temperature().value(),
            zone.pressure().value(),
            zone.interface_height().value(),
        ];
        assert!(values[0] > 2.0 && values[1] > values[0]);
        let work = 0.5
            * (old_pressure + values[2])
            * (zone.layer(Layer::Lower).volume().value() - before[0].volume().value());
        close(report.relaxation.interface_work_j[0], work, 2e-12);
        close(report.relaxation.interface_work_j[1], -work, 2e-12);
        // The upper layer expands and compresses the lower layer. Source heat
        // belongs only to upper chemistry; equal/opposite p dV stays internal.
        assert!(work < 0.0);
        let stored = sum(zone.layers.iter().zip(before).map(|(after, before)| {
            after.sensible_internal_energy().value() - before.sensible_internal_energy().value()
        }));
        close(
            stored,
            report.burning.actual_extent.value() * 40.0 * R,
            2e-8,
        );
        if let Some(expected) = baseline {
            for (actual, expected) in values.into_iter().zip(expected) {
                close(actual, expected, 2e-9);
            }
        } else {
            baseline = Some(values);
        }
    }
}

#[test]
fn g1_recombination_keeps_pressure_flow_work_out_of_closed_room_chemical_heat() {
    for reference in [0.0, 2.0_f64.powi(45), -2.0_f64.powi(45)] {
        // Manufactured 2 A -> A2: Delta nu=-1, Delta Cp/R=-1.5,
        // Delta Cv/R=-0.5 and retained formation change -40 R.
        let ids = vec![id("A"), id("A2")];
        let bank = IdealGasEquilibriumV1::new(
            ElementalMatrix::new(
                vec![ElementId::new("X").unwrap()],
                ids.clone(),
                vec![vec![1, 2]],
            )
            .unwrap(),
            ChargeVector::new(ids, vec![0, 0]).unwrap(),
            vec![
                source("A", 1.0, 2.5, 0.0, reference),
                source("A2", 2.0, 3.5, 0.0, 2.0 * reference - 40.0),
            ],
        )
        .unwrap();
        let law = ArrheniusRateV1::new(
            Temperature::new(2.0),
            ReactionProgressRateQuantityV1::new(1.0),
            0.0,
            MolarEnergyQuantityV1::new(0.0),
        )
        .unwrap();
        let reaction =
            ElementaryReactionV1::new(reaction_id(), vec![(id("A"), 2)], vec![(id("A2"), 1)], law)
                .unwrap();
        let mechanism =
            IdealGasKineticsV1::new(bank, "recombination-enthalpy", 1, vec![reaction], 0.0)
                .unwrap();
        let gas = MixtureGas::new(
            mechanism.thermodynamics().clone(),
            CaloricDomain::new(Temperature::new(1.0), Temperature::new(50.0)),
            || false,
        )
        .unwrap();
        let mut zone = TwoZone::new(
            gas,
            Area::new(0.5),
            Length::new(2.0),
            Length::new(0.0),
            Pressure::new(10.0 * R),
            [Temperature::new(2.0); 2],
            [vec![1.0, 0.0], vec![1.0, 0.0]],
            || false,
        )
        .unwrap();
        let owner = adapter(&zone, &mechanism);
        let power_extent = owner
            .extent_for_power(
                Power::new(20.0),
                Time::new(0.3),
                Temperature::new(8.0),
                || false,
            )
            .unwrap();
        close(power_extent.value(), 6.0 / (52.0 * R), 1e-16);
        let before_atoms = zone.layer(Layer::Upper).species_masses()[0];
        let report = owner
            .apply_split(
                &mut zone,
                Layer::Upper,
                Amount::new(0.5),
                Time::new(0.3),
                Temperature::new(1.0),
                ZoneSolve::default(),
                || false,
            )
            .unwrap();
        // At this same reference T, Delta h=-41.5R but Delta u=-40.5R.
        // The closed rigid room receives the latter sensible-gauge correction.
        close(report.burning.reaction_enthalpy.value(), -41.5 * R, 2e-12);
        close(
            report.burning.chemical_rates.referenced_energy_w,
            report.burning.actual_extent.value() / 0.3 * 40.5 * R,
            2e-12,
        );
        let after = zone.layer(Layer::Upper).species_masses();
        let (atoms, dimers) = (after[0], after[1] / 2.0);
        let expected_t = (1.5 * before_atoms * 2.0 + 40.0 * dimers) / (1.5 * atoms + 2.5 * dimers);
        let temperature = zone.layer(Layer::Upper).temperature().value();
        close(temperature, expected_t, 2e-9);
        close(atoms + 2.0 * dimers, before_atoms, 2e-15);
        close(
            R * ((1.5 * atoms + 2.5 * dimers) * temperature
                - 40.0 * dimers
                - 1.5 * before_atoms * 2.0),
            0.0,
            2e-8,
        );
        close(
            zone.pressure().value(),
            (atoms + dimers) * R * temperature,
            2e-9,
        );
    }
}
