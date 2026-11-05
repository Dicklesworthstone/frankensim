use super::*;

fn close(actual: f64, expected: f64, relative: f64) {
    assert!((actual - expected).abs() <= relative * expected.abs(),
        "actual={actual:.17e}, expected={expected:.17e}");
}

fn at_kev(reaction: FusionReactionV1, kev: f64) -> FusionReactivityV1 {
    reaction.reactivity(Temperature::new(kev * KELVIN_PER_KEV)).unwrap()
}

#[test]
fn g2_primary_table_viii_reactivities_and_si_conversion() {
    // Independent, four-significant-figure entries printed in table VIII,
    // p.625, in cm³/s. Channel order is DT, DDn, DDp, DHe3. The tolerance
    // covers decimal table rounding, not experimental/model uncertainty.
    let anchors = [
        (1.0, [6.857e-21, 9.933e-23, 1.017e-22, 3.057e-26]),
        (10.0, [1.136e-16, 6.023e-19, 5.781e-19, 2.126e-19]),
        (20.0, [4.330e-16, 2.603e-18, 2.399e-18, 3.482e-18]),
        (50.0, [8.649e-16, 1.133e-17, 9.838e-18, 5.554e-17]),
    ];
    for (kev, expected) in anchors {
        for (reaction, cm3_s) in FusionReactionV1::ALL.into_iter().zip(expected) {
            close(at_kev(reaction, kev).value_m3_s(), cm3_s * 1.0e-6, 5.0e-4);
        }
    }
    // The conventional 10 keV DT value is ~1.136e-22 m³/s, not cm³/s.
    close(at_kev(FusionReactionV1::DeuteriumTritium, 10.0).value_m3_s(),
        1.1361654705836233e-22, 3.0e-14);
}

#[test]
fn g3_analytic_temperature_slope_matches_independent_difference() {
    for reaction in FusionReactionV1::ALL {
        for kev in [1.0, 10.0, 30.0, 50.0, 90.0] {
            let t = kev * KELVIN_PER_KEV;
            let h = t * 1.0e-5;
            let rate = reaction.reactivity(Temperature::new(t)).unwrap();
            let plus = reaction.reactivity(Temperature::new(t + h)).unwrap().value_m3_s();
            let minus = reaction.reactivity(Temperature::new(t - h)).unwrap().value_m3_s();
            close(rate.derivative_m3_s_k(), (plus - minus) / (2.0 * h), 2.0e-8);
            close(rate.logarithmic_temperature_derivative(),
                t * rate.derivative_m3_s_k() / rate.value_m3_s(), 8.0 * f64::EPSILON);
        }
    }
    // DT reactivity turns over; a forced positive derivative is unphysical.
    assert!(at_kev(FusionReactionV1::DeuteriumTritium, 100.0).derivative_m3_s_k() < 0.0);
}

#[test]
fn g0_source_temperature_limits_are_inclusive_and_extrapolation_is_refused() {
    for reaction in FusionReactionV1::ALL {
        let [lower, upper] = reaction.temperature_range();
        assert!(reaction.reactivity(lower).is_ok());
        assert!(reaction.reactivity(upper).is_ok());
        for t in [lower.value().next_down(), upper.value().next_up()] {
            assert!(matches!(reaction.reactivity(Temperature::new(t)),
                Err(FusionErrorV1::OutsideTemperatureDomain { .. })));
        }
        for t in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert_eq!(reaction.reactivity(Temperature::new(t)),
                Err(FusionErrorV1::InvalidTemperature));
        }
    }
    assert_eq!(FusionReactionV1::DeuteriumHelium3.temperature_range_kev(), [0.5, 190.0]);
    assert!(FusionReactionV1::DeuteriumHelium3
        .reactivity(Temperature::new(0.3 * KELVIN_PER_KEV)).is_err());
}

#[test]
fn g0_nuclear_stoichiometry_conserves_baryons_and_charge() {
    for reaction in FusionReactionV1::ALL {
        let mut baryons = 0_i16;
        let mut charge = 0_i16;
        let mut mass_defect = 0.0;
        for species in NuclearSpeciesV1::ALL {
            let coefficient = reaction.stoichiometric_coefficient(species);
            baryons += i16::from(coefficient) * i16::from(species.baryon_number());
            charge += i16::from(coefficient) * i16::from(species.charge_number());
            mass_defect -= f64::from(coefficient) * species.mass_kg();
        }
        assert_eq!(baryons, 0);
        assert_eq!(charge, 0);
        assert!(mass_defect > 0.0);
        close(reaction.q_j(), mass_defect * SPEED_OF_LIGHT_M_S * SPEED_OF_LIGHT_M_S,
            2.0e-12);
    }
}

#[test]
fn g0_identical_pair_counting_preserves_two_deuteron_sink_per_event() {
    for reaction in [FusionReactionV1::DeuteriumDeuteriumNeutron,
        FusionReactionV1::DeuteriumDeuteriumProton] {
        assert_eq!(reaction.reactant_pair_factor(), 0.5);
        assert_eq!(reaction.stoichiometric_coefficient(NuclearSpeciesV1::Deuterium), -2);
        assert_eq!(reaction.reactants(), [NuclearSpeciesV1::Deuterium; 2]);
    }
    assert_eq!(FusionReactionV1::DeuteriumTritium.reactant_pair_factor(), 1.0);
    let ratio = |t| at_kev(FusionReactionV1::DeuteriumDeuteriumNeutron, t).value_m3_s()
        / at_kev(FusionReactionV1::DeuteriumDeuteriumProton, t).value_m3_s();
    // Neither an assumed equal branch split nor one shared fit reproduces this.
    assert!(ratio(1.0) < 1.0);
    assert!(ratio(50.0) > 1.1);
}

#[test]
fn g0_q_energy_partition_and_reaction_weighted_thermal_energy() {
    let q_mev = [17.5892479, 3.2688570, 4.0326639, 18.3530549];
    for (reaction, expected) in FusionReactionV1::ALL.into_iter().zip(q_mev) {
        close(reaction.q_j() / (1.0e6 * ELECTRONVOLT_J), expected, 2.0e-8);
        close(reaction.charged_q_j() + reaction.neutron_q_j(), reaction.q_j(), 4.0*f64::EPSILON);
        let fractions = reaction.product_relative_energy_fractions();
        assert_eq!(fractions[0] + fractions[1], 1.0);
        let products = reaction.products();
        close(fractions[0] * products[0].mass_kg(),
            fractions[1] * products[1].mass_kg(), 4.0*f64::EPSILON);
        assert_eq!(reaction.neutron_q_j() > 0.0, products[1] == NuclearSpeciesV1::Neutron);
        for kev in [1.0, 10.0, 90.0] {
            let rate = at_kev(reaction, kev);
            assert!(rate.mean_reacting_relative_energy_j().is_finite());
            assert!(rate.mean_reacting_relative_energy_j() > 0.0);
        }
    }
    // Reacting pairs preferentially sample the energetic tail at 10 keV;
    // debiting two unweighted average ion energies would undercount this.
    let rate = at_kev(FusionReactionV1::DeuteriumTritium, 10.0);
    assert!(rate.mean_reacting_relative_energy_j() >
        1.5 * BOLTZMANN_J_PER_K * rate.temperature().value());
}

#[test]
fn g5_repeated_evaluations_have_identical_bits() {
    for reaction in FusionReactionV1::ALL {
        let a = at_kev(reaction, 15.0);
        let b = at_kev(reaction, 15.0);
        assert_eq!(a.value_m3_s().to_bits(), b.value_m3_s().to_bits());
        assert_eq!(a.derivative_m3_s_k().to_bits(), b.derivative_m3_s_k().to_bits());
    }
}
