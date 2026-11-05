//! G0/G1 six-source isotope identity, thermal support and plasma admission.

use super::*;

// Deliberately neither nuclear order nor canonical lexical order.
const NAMES: [&str; 6] = [
    "z_proton",
    "b_deuterium",
    "x_tritium",
    "a_helium3",
    "m_helium4",
    "electron",
];

fn sources() -> Vec<Nasa9StandardStateModelV1> {
    [
        NuclearSpeciesV1::Proton.mass_kg(),
        NuclearSpeciesV1::Deuterium.mass_kg(),
        NuclearSpeciesV1::Tritium.mass_kg(),
        NuclearSpeciesV1::Helium3.mass_kg(),
        NuclearSpeciesV1::Helium4.mass_kg(),
        ELECTRON_MASS_KG,
    ]
    .into_iter()
    .zip(NAMES)
    .map(|(mass, name)| source(name, mass))
    .collect()
}

fn bank(models: Vec<Nasa9StandardStateModelV1>, charges: Vec<i128>) -> IdealGasEquilibriumV1 {
    let names: Vec<_> = NAMES.map(id).into();
    IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("H").unwrap(), ElementId::new("He").unwrap()],
            names.clone(),
            vec![vec![1, 1, 1, 0, 0, 0], vec![0, 0, 0, 1, 1, 0]],
        )
        .unwrap(),
        ChargeVector::new(names, charges).unwrap(),
        models,
    )
    .unwrap()
}

fn bind(bank: IdealGasEquilibriumV1) -> Result<LightIonFusionPlasmaV1, DtFusionPlasmaErrorV1> {
    LightIonFusionPlasmaV1::new(
        bank,
        id(NAMES[0]),
        id(NAMES[1]),
        id(NAMES[2]),
        id(NAMES[3]),
        id(NAMES[4]),
    )
}

fn fixture() -> LightIonFusionPlasmaV1 {
    bind(bank(sources(), vec![1, 1, 1, 2, 2, -1])).unwrap()
}

fn amounts(plasma: &LightIonFusionPlasmaV1, heavy: [f64; 5]) -> [Amount; 6] {
    let mut result = [Amount::new(0.0); 6];
    for (index, amount) in plasma.indices().into_iter().zip(heavy) {
        result[index] = Amount::new(amount);
    }
    result[plasma.electron_index()] =
        Amount::new(heavy[0] + heavy[1] + heavy[2] + 2.0 * heavy[3] + 2.0 * heavy[4]);
    result
}

fn state(
    plasma: &LightIonFusionPlasmaV1,
    heavy: [f64; 5],
    th: f64,
    te: f64,
) -> FrozenIonizedGasStateV1 {
    FrozenIonizedGasStateV1::new(
        plasma.thermodynamics(),
        &amounts(plasma, heavy),
        Volume::new(2.0),
        Temperature::new(th),
        Temperature::new(te),
        1e-12,
        || false,
    )
    .unwrap()
}

#[test]
fn g0_six_original_sources_retain_explicit_isotopes_and_canonical_order() {
    let original = bank(sources(), vec![1, 1, 1, 2, 2, -1]);
    let plasma = bind(original.clone()).unwrap();
    assert_eq!(plasma.indices(), [5, 1, 4, 0, 3]);
    assert_eq!(plasma.electron_index(), 2);
    assert!(core::ptr::eq(
        original.models().as_ptr(),
        plasma.thermodynamics().models().as_ptr()
    ));
    let heavy = [0.0005, 0.001, 0.002, 0.0003, 0.0007];
    let point = state(&plasma, heavy, 1e8, 2e7);
    let before = point.clone();
    plasma.admit_state(&point, || false).unwrap();
    plasma
        .admit_state_with_limits(&point, limits(), || false)
        .unwrap();
    assert_eq!(point, before);
    let r = crate::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
    let nh: f64 = heavy.iter().sum();
    let ne = heavy[0] + heavy[1] + heavy[2] + 2.0 * (heavy[3] + heavy[4]);
    assert!((point.pressure().value() / (r * (nh * 1e8 + ne * 2e7) / 2.0) - 1.0).abs() < 1e-14);
    assert!((point.heavy_internal_energy().value() / (1.5 * r * nh * 1e8) - 1.0).abs() < 1e-14);
    assert!((point.electron_internal_energy().value() / (1.5 * r * ne * 2e7) - 1.0).abs() < 1e-14);
}

#[test]
fn g0_rejects_helium_isotope_swaps_aliases_atomic_masses_and_incomplete_banks() {
    let original = bank(sources(), vec![1, 1, 1, 2, 2, -1]);
    for mapping in [
        [NAMES[0], NAMES[1], NAMES[2], NAMES[4], NAMES[3]],
        [NAMES[1], NAMES[0], NAMES[2], NAMES[3], NAMES[4]],
        [NAMES[0], NAMES[1], NAMES[1], NAMES[3], NAMES[4]],
        [NAMES[0], NAMES[1], NAMES[2], "missing", NAMES[4]],
    ] {
        let [p, d, t, he3, he4] = mapping.map(id);
        assert!(LightIonFusionPlasmaV1::new(original.clone(), p, d, t, he3, he4).is_err());
    }
    let mut models = sources();
    models[3] = source(
        NAMES[3],
        NuclearSpeciesV1::Helium3.mass_kg() + 2.0 * ELECTRON_MASS_KG,
    );
    assert!(bind(bank(models, vec![1, 1, 1, 2, 2, -1])).is_err());
    assert!(bind(bank(sources(), vec![1, 1, 1, 1, 2, -1])).is_err());
    // The original bank already owns structural-electron admission; an
    // invalid positive, nucleus-free carrier cannot reach the fusion binding.
    let ids = NAMES.map(id).to_vec();
    assert!(
        IdealGasEquilibriumV1::new(
            ElementalMatrix::new(
                vec![ElementId::new("H").unwrap(), ElementId::new("He").unwrap()],
                ids.clone(),
                vec![vec![1, 1, 1, 0, 0, 0], vec![0, 0, 0, 1, 1, 0]],
            )
            .unwrap(),
            ChargeVector::new(ids, vec![1, 1, 1, 2, 2, 1]).unwrap(),
            sources(),
        )
        .is_err()
    );
    let dt = super::fixture();
    assert!(bind(dt.thermodynamics().clone()).is_err());
}

#[test]
fn g0_secondary_products_require_original_support_and_sensible_references() {
    // Every absent possible product is checked, even for initially pure D+.
    for product in [0, 2, 3, 4] {
        let mut models = sources();
        models[product] = changed_source(&models[product], |card| {
            card.validity = ValidityDomain::unconstrained().with("T", 1e4, 1e7)
        });
        let plasma = bind(bank(models, vec![1, 1, 1, 2, 2, -1])).unwrap();
        let initial = state(&plasma, [0.0, 0.001, 0.0, 0.0, 0.0], 1e8, 2e7);
        assert!(matches!(plasma.admit_state(&initial, || false),
            Err(DtFusionPlasmaErrorV1::Source { species, .. }) if species == id(NAMES[product])));
    }
    for (parameter, value) in [("a7", 1.0), ("a2", 3.5), ("a3", 1e-8)] {
        let mut models = sources();
        models[0] = changed_source(&models[0], |card| {
            card.parameters.get_mut(parameter).unwrap().value = value
        });
        assert!(bind(bank(models, vec![1, 1, 1, 2, 2, -1])).is_err());
    }
    let mut models = sources();
    models[3] = models[3]
        .with_caloric_continuity(
            Temperature::new(1e8),
            MolarEnergyQuantityV1::new(0.0),
            MolarThermalQuantityV1::new(0.0),
        )
        .unwrap();
    assert!(bind(bank(models, vec![1, 1, 1, 2, 2, -1])).is_err());
}

#[test]
fn g1_six_species_scales_match_direct_formulas_including_helium3_charge_square() {
    let plasma = fixture();
    for heavy in [
        [0.0005, 0.001, 0.002, 0.0003, 0.0007],
        [0.0, 0.0, 0.0, 0.001, 0.0],
    ] {
        let initial = state(&plasma, heavy, 1e7, 1e8);
        let actual = plasma
            .admit_state_with_limits(&initial, limits(), || false)
            .unwrap();
        let densities: Vec<_> = (0..6)
            .map(|i| initial.number_density_per_m3(i).unwrap())
            .collect();
        let total_density: f64 = densities.iter().sum();
        let pi = core::f64::consts::PI;
        let radius = (3.0 / (4.0 * pi * total_density)).cbrt();
        let mut expected = [0.0_f64; 3];
        for (i, n) in densities.into_iter().enumerate() {
            if n == 0.0 {
                continue;
            }
            let mass = plasma.thermodynamics().models()[i].molar_mass().value() / AVOGADRO_PER_MOL;
            let t = initial.species_temperature(i).unwrap().value();
            let z = plasma.thermodynamics().charge_vector().get(i).unwrap() as f64;
            expected[0] = expected[0].max(
                z * z * ELEMENTARY_CHARGE_C.powi(2)
                    / (4.0 * pi * VACUUM_PERMITTIVITY_F_M * radius * BOLTZMANN_J_PER_K * t),
            );
            expected[1] = expected[1]
                .max(n * (PLANCK_J_S / (2.0 * pi * mass * BOLTZMANN_J_PER_K * t).sqrt()).powi(3));
            expected[2] =
                expected[2].max(BOLTZMANN_J_PER_K * t / (mass * SPEED_OF_LIGHT_M_S.powi(2)));
        }
        for (got, want) in [
            actual.maximum_coulomb_coupling,
            actual.maximum_thermal_phase_space_density,
            actual.maximum_thermal_energy_rest_energy_ratio,
        ]
        .into_iter()
        .zip(expected)
        {
            assert!((got / want - 1.0).abs() < 1e-12, "{got} != {want}");
        }
        let mut restrictive = limits();
        restrictive.maximum_coulomb_coupling = 0.5 * expected[0];
        assert_eq!(
            plasma.admit_state_with_limits(&initial, restrictive, || false),
            Err(DtFusionPlasmaErrorV1::Domain("Coulomb coupling"))
        );
    }
}

#[test]
fn g0_network_refuses_changed_sources_unbalanced_inventory_and_cancellation() {
    let plasma = fixture();
    let initial = state(&plasma, [0.0, 0.001, 0.002, 0.0, 0.0], 1e8, 2e7);
    let mut models = sources();
    models[3] = changed_source(&models[3], |card| {
        card.provenance.source.push_str(" alternate source")
    });
    let other = bind(bank(models, vec![1, 1, 1, 2, 2, -1])).unwrap();
    assert!(other.admit_state(&initial, || false).is_err());
    let original = amounts(&plasma, [0.0001, 0.001, 0.002, 0.0002, 0.0003]);
    let admit = |amounts: &[Amount], volume, physical_limits| {
        plasma.admit_point(
            Temperature::new(1e8),
            Temperature::new(2e7),
            amounts,
            Volume::new(volume),
            physical_limits,
            || false,
        )
    };
    admit(&original, 2.0, limits()).unwrap();
    for bad in [-1.0, f64::INFINITY, f64::NAN] {
        let mut changed = original;
        changed[plasma.indices()[3]] = Amount::new(bad);
        assert!(admit(&changed, 2.0, limits()).is_err());
    }
    assert!(admit(&original[..5], 2.0, limits()).is_err());
    assert!(admit(&[Amount::new(0.0); 6], 2.0, limits()).is_err());
    assert!(admit(&original, 0.0, limits()).is_err());
    let mut changed = original;
    changed[plasma.electron_index()] = Amount::new(0.0);
    assert!(admit(&changed, 2.0, limits()).is_err());
    let mut changed = original;
    changed[plasma.indices()[3]] = Amount::new(2.0 * changed[plasma.indices()[3]].value());
    assert!(admit(&changed, 2.0, limits()).is_err());
    assert_eq!(
        plasma.admit_state_with_limits(&initial, limits(), || true),
        Err(DtFusionPlasmaErrorV1::Cancelled)
    );
    let before = initial.clone();
    let mut polls = 0;
    assert_eq!(
        plasma.admit_state_with_limits(&initial, limits(), || {
            polls += 1;
            polls == 4
        }),
        Err(DtFusionPlasmaErrorV1::Cancelled)
    );
    assert_eq!(initial, before);
    // Full fuel exhaustion can leave proton/He4 ash without artificial seeds.
    admit(
        &amounts(&plasma, [0.001, 0.0, 0.0, 0.0, 0.002]),
        2.0,
        limits(),
    )
    .unwrap();
}

#[test]
fn g0_network_temperature_support_checks_secondary_products_and_original_gaps() {
    for (second_start, entropy_jump) in [(1e7, false), (2e7, false), (1e7, true)] {
        let mut models = sources();
        let mut lower = models[3].regions()[0].card().clone();
        lower.validity = ValidityDomain::unconstrained().with("T", 1e4, 1e7);
        let mut upper = lower.clone();
        upper.validity = ValidityDomain::unconstrained().with("T", second_start, 2e9);
        if entropy_jump {
            upper.parameters.get_mut("a8").unwrap().value = 1.0;
        }
        models[3] = Nasa9StandardStateModelV1::new(
            models[3].species().clone(),
            models[3].molar_mass(),
            models[3].convention().clone(),
            vec![
                Nasa9RegionV1::from_card(lower).unwrap(),
                Nasa9RegionV1::from_card(upper).unwrap(),
            ],
        )
        .unwrap();
        let result = bind(bank(models, vec![1, 1, 1, 2, 2, -1]));
        if entropy_jump {
            assert!(result.is_err());
            continue;
        }
        let plasma = result.unwrap();
        assert_eq!(
            plasma
                .admit_temperature_interval(
                    Temperature::new(1e6),
                    Temperature::new(1e8),
                    Temperature::new(1e5),
                    Temperature::new(1e7),
                    || false,
                )
                .is_ok(),
            second_start == 1e7
        );
    }
}
