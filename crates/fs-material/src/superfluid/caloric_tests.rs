use super::*;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "{actual:.17e} != {expected:.17e}"
    );
}

fn reference() -> SuperfluidCaloricReference {
    SuperfluidCaloricReference {
        pressure_pa: 100_000.0,
        density_at_min_temperature_kg_m3: 145.0,
        specific_enthalpy_at_min_temperature_j_kg: -20.0,
        specific_entropy_at_min_temperature_j_kg_k: 100.0,
    }
}

// Synthetic authored data test the constitutive mathematics; these numbers are
// deliberately not presented as measured helium properties or as a helium fit.
fn knot(t: f64, cp: f64, alpha: f64) -> SuperfluidCaloricKnot {
    SuperfluidCaloricKnot {
        temperature_k: t,
        specific_heat_cp_j_kg_k: cp,
        superfluid_fraction: 0.95 - 0.2 * t,
        adiabatic_sound_speed_m_s: 240.0 - 5.0 * t,
        thermal_expansion_per_k: alpha,
        normal_viscosity_pa_s: 1.2e-6 + 0.1e-6 * t,
    }
}

fn law(knots: Vec<SuperfluidCaloricKnot>) -> IsobaricSuperfluidLaw {
    IsobaricSuperfluidLaw::try_new(ContentHash([0x43; 32]), reference(), knots).unwrap()
}

#[test]
fn integrated_linear_cp_agrees_with_independent_analytic_enthalpy_and_entropy() {
    // Cp(T)=600 T: h-h0=300(T²-1), s-s0=600(T-1).
    let law = law(vec![knot(1.0, 600.0, -0.01), knot(2.0, 1200.0, -0.01)]);
    for t in [1.0, 1.01, 1.25, 1.5, 1.99, 2.0] {
        let state = law.state_at_temperature_k(t).unwrap();
        close(
            state.specific_enthalpy_j_kg(),
            -20.0 + 300.0 * (t * t - 1.0),
            2e-14,
        );
        close(
            state.specific_entropy_j_kg_k(),
            100.0 + 600.0 * (t - 1.0),
            2e-14,
        );
        close(
            state.two_fluid().density_kg_m3(),
            145.0 * (0.01 * (t - 1.0)).exp(),
            2e-14,
        );
        close(
            state.two_fluid().specific_heat_cp_j_kg_k(),
            600.0 * t,
            2e-14,
        );
        close(
            law.temperature_at_specific_enthalpy(state.specific_enthalpy_j_kg())
                .unwrap(),
            t,
            2e-14,
        );
    }
    assert_eq!(law.pressure_pa(), 100_000.0);
    assert_eq!(law.temperature_bounds_k(), [1.0, 2.0]);
    assert_eq!(law.specific_enthalpy_bounds_j_kg(), [-20.0, 880.0]);
}

#[test]
fn fixed_density_subset_has_the_constant_cp_logarithmic_entropy() {
    let law = law(vec![knot(1.0, 700.0, 0.0), knot(2.0, 700.0, 0.0)]);
    let a = law.state_at_temperature_k(1.0).unwrap();
    let b = law.state_at_temperature_k(2.0).unwrap();
    assert_eq!(a.two_fluid().density_kg_m3(), b.two_fluid().density_kg_m3());
    close(
        b.specific_entropy_j_kg_k(),
        100.0 + 700.0 * 2.0_f64.ln(),
        2e-14,
    );
    close(
        law.specific_entropy_change_j_kg_k(-20.0, 680.0).unwrap(),
        700.0 * 2.0_f64.ln(),
        2e-14,
    );
    close(b.two_fluid().specific_heat_cv_j_kg_k(), 700.0, 2e-14);
    assert_ne!(a.two_fluid().sound_modes(), b.two_fluid().sound_modes());
}

#[test]
fn differential_thermodynamic_identities_hold_between_knots() {
    let law = law(vec![
        knot(1.0, 600.0, -0.01),
        knot(1.8, 1100.0, -0.003),
        knot(2.4, 1800.0, 0.004),
    ]);
    for t in [1.2, 1.6, 2.0, 2.3] {
        let eps = 1e-5;
        let left = law.state_at_temperature_k(t - eps).unwrap();
        let state = law.state_at_temperature_k(t).unwrap();
        let right = law.state_at_temperature_k(t + eps).unwrap();
        let cp = state.two_fluid().specific_heat_cp_j_kg_k();
        close(
            (right.specific_enthalpy_j_kg() - left.specific_enthalpy_j_kg()) / (2.0 * eps),
            cp,
            1e-9,
        );
        close(
            (right.specific_entropy_j_kg_k() - left.specific_entropy_j_kg_k()) / (2.0 * eps),
            cp / t,
            1e-9,
        );
        close(
            (right.two_fluid().density_kg_m3().ln() - left.two_fluid().density_kg_m3().ln())
                / (2.0 * eps),
            -state.two_fluid().thermal_expansion_per_k(),
            1e-9,
        );
        let p = state.two_fluid();
        close(
            p.isothermal_compressibility_per_pa() - p.adiabatic_compressibility_per_pa(),
            t * p.thermal_expansion_per_k().powi(2) / (p.density_kg_m3() * cp),
            1e-20,
        );
        assert!(p.sound_modes().slow_m_s > 0.0);
        assert!(p.sound_modes().fast_m_s >= p.sound_modes().slow_m_s);
    }
}

#[test]
fn sign_changing_expansion_has_its_interior_density_maximum() {
    let law = law(vec![knot(1.0, 700.0, -0.01), knot(2.0, 700.0, 0.01)]);
    let middle = law.state_at_temperature_k(1.5).unwrap();
    close(
        middle.two_fluid().density_kg_m3(),
        145.0 * 0.0025_f64.exp(),
        2e-14,
    );
    close(middle.two_fluid().thermal_expansion_per_k(), 0.0, 2e-14);
    close(
        law.state_at_temperature_k(2.0)
            .unwrap()
            .two_fluid()
            .density_kg_m3(),
        145.0,
        2e-14,
    );
}

#[test]
fn inverse_is_monotone_for_increasing_decreasing_and_nearly_constant_capacity() {
    for capacities in [
        [500.0, 1500.0, 800.0],
        [1500.0, 500.0, 700.0],
        [700.0, 700.0 + 1e-10, 700.0 - 1e-10],
    ] {
        let law = law(vec![
            knot(1.0, capacities[0], -0.01),
            knot(1.6, capacities[1], -0.005),
            knot(2.4, capacities[2], 0.0),
        ]);
        let mut previous_h = f64::NEG_INFINITY;
        for i in 0..=140 {
            let t = if i == 140 {
                2.4
            } else {
                1.0 + 0.01 * f64::from(i)
            };
            let h = law.specific_enthalpy_at_temperature_k(t).unwrap();
            assert!(h > previous_h);
            previous_h = h;
            close(law.temperature_at_specific_enthalpy(h).unwrap(), t, 2e-14);
            let state = law.state_at_specific_enthalpy(h).unwrap();
            assert_eq!(state.specific_enthalpy_j_kg().to_bits(), h.to_bits());
        }
        let min_cp = capacities.into_iter().fold(f64::INFINITY, f64::min);
        assert!(law.maximum_temperature_slope() >= 1.0 / min_cp);
        close(law.maximum_temperature_slope(), 1.0 / min_cp, 1e-16);
    }
}

#[test]
fn entropy_is_additive_antisymmetric_and_retains_sub_temperature_ulp_heat() {
    let law = law(vec![
        knot(1.0, 600.0, -0.01),
        knot(1.6, 960.0, -0.005),
        knot(2.4, 1440.0, 0.0),
    ]);
    let h = |t| law.specific_enthalpy_at_temperature_k(t).unwrap();
    let ds = |a, b| law.specific_entropy_change_j_kg_k(a, b).unwrap();
    close(ds(h(1.2), h(2.2)), 600.0, 2e-14);
    close(
        ds(h(1.2), h(2.2)),
        ds(h(1.2), h(1.6)) + ds(h(1.6), h(2.2)),
        2e-14,
    );
    assert_eq!(ds(h(2.2), h(1.2)), -ds(h(1.2), h(2.2)));
    assert_eq!(ds(h(1.6), h(1.6)), 0.0);
    let mut r = reference();
    r.specific_enthalpy_at_min_temperature_j_kg = 0.0;
    let tiny = IsobaricSuperfluidLaw::try_new(
        ContentHash([0x43; 32]),
        r,
        vec![knot(1.0, 600.0, 0.0), knot(2.0, 600.0, 0.0)],
    )
    .unwrap();
    assert_eq!(tiny.temperature_at_specific_enthalpy(1e-20).unwrap(), 1.0);
    assert_eq!(
        tiny.specific_entropy_change_j_kg_k(0.0, 1e-20).unwrap(),
        1e-20
    );
}

#[test]
fn native_thermal_interface_owns_the_same_enthalpy_entropy_and_slope() {
    let law = law(vec![knot(1.0, 600.0, -0.01), knot(2.0, 1200.0, -0.01)]);
    let thermal: &dyn ThermalLaw = &law;
    assert_eq!(thermal.storage_basis(), ThermalStorageBasis::Enthalpy);
    assert_eq!(thermal.caloric_identity(), Some(*law.identity().as_bytes()));
    close(thermal.temperature_k(355.0).unwrap(), 1.5, 2e-14);
    assert_eq!(
        thermal.maximum_temperature_slope(),
        law.maximum_temperature_slope()
    );
    close(
        thermal
            .specific_entropy_change_j_kg_k(-20.0, 880.0)
            .unwrap()
            .unwrap(),
        600.0,
        2e-14,
    );
    assert_eq!(
        thermal.temperature_k(881.0),
        Err(DepositionError::MaterialDomain)
    );
    assert_eq!(
        thermal.specific_entropy_change_j_kg_k(-21.0, 880.0),
        Err(DepositionError::MaterialDomain)
    );
}

#[test]
fn caloric_identity_binds_datum_pressure_reference_and_every_authored_property() {
    let data = vec![knot(1.0, 600.0, -0.01), knot(2.0, 1200.0, -0.01)];
    let original = law(data.clone());
    assert_eq!(law(data.clone()).identity(), original.identity());
    let mut pressure = reference();
    pressure.pressure_pa += 1.0;
    assert_ne!(
        IsobaricSuperfluidLaw::try_new(ContentHash([0x43; 32]), pressure, data.clone())
            .unwrap()
            .identity(),
        original.identity()
    );
    let mut datum = reference();
    datum.specific_enthalpy_at_min_temperature_j_kg += 1.0;
    let changed =
        IsobaricSuperfluidLaw::try_new(ContentHash([0x43; 32]), datum, data.clone()).unwrap();
    assert_ne!(changed.identity(), original.identity());
    close(
        changed
            .specific_entropy_change_j_kg_k(-19.0, 881.0)
            .unwrap(),
        original
            .specific_entropy_change_j_kg_k(-20.0, 880.0)
            .unwrap(),
        2e-14,
    );
    let mut knots = data.clone();
    knots[1].normal_viscosity_pa_s *= 1.01;
    assert_ne!(law(knots).identity(), original.identity());
    assert_ne!(
        IsobaricSuperfluidLaw::try_new(ContentHash([0x44; 32]), reference(), data)
            .unwrap()
            .identity(),
        original.identity()
    );
}

#[test]
fn invalid_domains_endpoints_and_nonfinite_calorics_are_refused() {
    let data = vec![knot(1.0, 600.0, -0.01), knot(2.0, 1200.0, -0.01)];
    let original = law(data.clone());
    for t in [0.0, 0.999, 2.001, f64::INFINITY, f64::NAN] {
        assert!(original.state_at_temperature_k(t).is_err());
    }
    for h in [-21.0, 881.0, f64::INFINITY, f64::NAN] {
        assert!(original.temperature_at_specific_enthalpy(h).is_err());
    }
    for cp in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        let mut bad = data.clone();
        bad[0].specific_heat_cp_j_kg_k = cp;
        assert!(IsobaricSuperfluidLaw::try_new(ContentHash([0x43; 32]), reference(), bad).is_err());
    }
    for fraction in [0.0, 1.0, -0.1, 1.1, f64::NAN] {
        let mut bad = data.clone();
        bad[0].superfluid_fraction = fraction;
        assert!(IsobaricSuperfluidLaw::try_new(ContentHash([0x43; 32]), reference(), bad).is_err());
    }
    let mut bad = data.clone();
    bad[1].temperature_k = 1.0;
    assert!(IsobaricSuperfluidLaw::try_new(ContentHash([0x43; 32]), reference(), bad).is_err());
    let mut invalid_reference = reference();
    invalid_reference.specific_enthalpy_at_min_temperature_j_kg = 1e100;
    assert!(
        IsobaricSuperfluidLaw::try_new(ContentHash([0x43; 32]), invalid_reference, data.clone())
            .is_err()
    );
    assert!(
        IsobaricSuperfluidLaw::try_new(ContentHash([0x43; 32]), reference(), vec![data[0]])
            .is_err()
    );
    assert!(
        IsobaricSuperfluidLaw::try_new(
            ContentHash([0x43; 32]),
            reference(),
            vec![data[0]; MAX_SUPERFLUID_CALORIC_KNOTS + 1]
        )
        .is_err()
    );
}

#[test]
fn response_is_continuous_at_knots_and_interior_density_failure_is_not_hidden() {
    let law = law(vec![
        knot(1.0, 600.0, -0.01),
        knot(1.6, 1000.0, 0.003),
        knot(2.4, 1400.0, 0.005),
    ]);
    let left = law.state_at_temperature_k(1.6 - 1e-9).unwrap();
    let middle = law.state_at_temperature_k(1.6).unwrap();
    let right = law.state_at_temperature_k(1.6 + 1e-9).unwrap();
    for state in [left, right] {
        close(
            state.specific_enthalpy_j_kg(),
            middle.specific_enthalpy_j_kg(),
            1e-8,
        );
        close(
            state.specific_entropy_j_kg_k(),
            middle.specific_entropy_j_kg_k(),
            1e-8,
        );
        close(
            state.two_fluid().density_kg_m3(),
            middle.two_fluid().density_kg_m3(),
            1e-8,
        );
        close(
            state.two_fluid().sound_modes().slow_m_s,
            middle.two_fluid().sound_modes().slow_m_s,
            1e-8,
        );
    }
    // Both endpoint expansion integrals are zero overall; the interior maximum
    // still overflows and must be caught during admission.
    assert!(
        IsobaricSuperfluidLaw::try_new(
            ContentHash([0x43; 32]),
            reference(),
            vec![knot(1.0, 600.0, -3000.0), knot(2.0, 1200.0, 3000.0)]
        )
        .is_err()
    );
}
