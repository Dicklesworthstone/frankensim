//! G1/G3/G4 material-source and pressure/enthalpy conformance.
use super::*;
use crate::thermochemical_gas::energy::equilibrium::{
    GasEquilibriumChamberSolve,
    tests::{fixture, isomers},
};
use crate::thermochemical_gas::tests::{
    point, species, species_with_energy_reference, transport, transport_in_domain,
};
use fs_evidence::ValidityDomain;
use fs_thermochem::equilibrium::composition::EquilibriumSolveV1;
use fs_thermochem::{
    Nasa9RegionV1, Nasa9StandardStateModelV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
};

fn controls() -> GasEquilibriumEnthalpySolve {
    let mut equilibrium = EquilibriumSolveV1::default();
    equilibrium.relative_balance_tolerance = 1e-13;
    equilibrium.dimensionless_potential_tolerance = 1e-12;
    GasEquilibriumEnthalpySolve {
        equilibrium: EquilibriumEnthalpySolveV1 {
            equilibrium,
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 1500.0,
            enthalpy_tolerance_j: 1e-7,
            maximum_temperature_iterations: 128,
        },
        mass_tolerance_kg: 1e-12,
        enthalpy_tolerance_j: 1e-6,
    }
}
fn feed() -> [(SpeciesId, Amount); 1] {
    [(SpeciesId::new("A").unwrap(), Amount::new(2.0))]
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}
fn expected_temperature(heat: f64) -> f64 {
    let y = |t: f64| 1.0 / (1.0 + (-24000.0 / (R * t)).exp());
    let (mut low, mut high) = (200.0, 1500.0);
    for _ in 0..80 {
        let t = 0.5 * (low + high);
        let h = 2.0 * (4.0 * R * (t - 900.0) - 24000.0 * (y(t) - y(900.0)));
        if h < heat {
            low = t;
        } else {
            high = t;
        }
    }
    0.5 * (low + high)
}

#[test]
fn g1_reacting_heating_and_cooling_preserve_actual_material_sources_and_mass() {
    let (bank, chemistry) = isomers();
    let original = bank.clone();
    let context = point(900.0, 1e5).with("sample_id", 7.0).unwrap();
    let mut path =
        EquilibriumGasEnthalpyPath::new(&bank, &chemistry, &context, &feed(), controls(), || false)
            .unwrap();
    let reference = path.reference_state().clone();
    for (heat, pressure) in [(4000.0, 1e5), (8000.0, 2e5), (-6000.0, 8e4)] {
        let state = path
            .sample(Pressure::new(pressure), Energy::new(heat), || false)
            .unwrap();
        let gas = state.parcel();
        let t = expected_temperature(heat);
        close(gas.state().state().temperature, t, 2e-8);
        assert_eq!(
            gas.mass_kg().to_bits(),
            reference.parcel().mass_kg().to_bits()
        );
        close(gas.volume_m3(), 2.0 * R * t / pressure, 1e-12);
        close(
            state.balance().enthalpy_change_j,
            heat,
            controls().enthalpy_tolerance_j,
        );
        assert_eq!(gas.point().axes()["sample_id"], 7.0);
        assert_eq!(
            gas.model().thermodynamics().components(),
            chemistry.models()
        );
        assert_eq!(
            gas.model().transport_materials(),
            bank.transport_materials()
        );
        let frozen_cp = gas.mass_kg() * gas.state().state().specific_heat_cp;
        assert!(state.native_sample().response.heat_capacity_p_j_k > frozen_cp);
        assert!(state.equilibrium().equilibrium().is_bound_to(&chemistry));
        assert_eq!(state.equilibrium().equilibrium().target_enthalpy(), None);
    }
    assert_eq!(path.reference_state(), &reference);
    assert_eq!(bank, original);
}

#[test]
fn g3_small_material_enthalpy_increment_survives_large_formation_datums() {
    for offset in [-1e13, 0.0, 1e13] {
        let (bank, chemistry) = fixture(
            vec![species_with_energy_reference("A", 0.01, 4.0, 0.0, offset)],
            vec![transport("A", 1.7e-5, 0.025)],
            &["H"],
            vec![vec![1]],
        );
        let mut c = controls();
        c.equilibrium.enthalpy_tolerance_j = 1e-8;
        c.enthalpy_tolerance_j = 1e-8;
        c.equilibrium.maximum_temperature_iterations = 1;
        let mut path = EquilibriumGasEnthalpyPath::new(
            &bank,
            &chemistry,
            &point(900.0, 1e5),
            &feed(),
            c,
            || false,
        )
        .unwrap();
        let result = path
            .sample(Pressure::new(30000.0), Energy::new(0.001), || false)
            .unwrap();
        close(
            result.parcel().state().state().temperature,
            900.0 + 0.001 / (8.0 * R),
            1e-12,
        );
        close(result.balance().enthalpy_change_j, 0.001, 1e-8);
        assert_eq!(result.parcel().mass_kg(), 0.02);
        assert_eq!(result.native_sample().temperature_iterations, 1);
    }
}

#[test]
fn g4_material_refusal_and_final_cancellation_invalidate_only_the_current_query() {
    let (_, chemistry) = isomers();
    let (bank, _) = fixture(
        chemistry.models().to_vec(),
        vec![
            transport_in_domain("A", 1.7e-5, 0.025, [200.0, 950.0], [1000.0, 1e7]),
            transport_in_domain("B", 2.0e-5, 0.15, [200.0, 950.0], [1000.0, 1e7]),
        ],
        &["H"],
        vec![vec![1, 1]],
    );
    let make = || {
        EquilibriumGasEnthalpyPath::new(
            &bank,
            &chemistry,
            &point(900.0, 1e5),
            &feed(),
            controls(),
            || false,
        )
        .unwrap()
    };
    let mut path = make();
    let reference = path.reference_state().clone();
    assert!(matches!(
        path.sample(Pressure::new(1e5), Energy::new(10000.0), || false),
        Err(GasEquilibriumEnthalpyError::Material(_))
    ));
    assert!(path.current().is_none());
    assert_eq!(path.reference_state(), &reference);
    assert!(
        path.sample(Pressure::new(1e5), Energy::new(1000.0), || false)
            .is_ok()
    );
    let mut complete = make();
    let mut count = 0;
    complete
        .sample(Pressure::new(1e5), Energy::new(1000.0), || {
            count += 1;
            false
        })
        .unwrap();
    let mut path = make();
    let mut polls = 0;
    assert!(matches!(
        path.sample(Pressure::new(1e5), Energy::new(1000.0), || {
            polls += 1;
            polls == count
        }),
        Err(GasEquilibriumEnthalpyError::Cancelled)
    ));
    assert_eq!(polls, count);
    assert!(path.current().is_none());
    assert_eq!(path.reference_state(), &reference);
    assert!(
        path.sample(Pressure::new(1e5), Energy::new(1000.0), || false)
            .is_ok()
    );
}

#[test]
fn g3_material_hp_cannot_bypass_the_original_inlet_source_path() {
    for gap in [0.0, 100.0] {
        let source = species("A", 0.01, 4.0, 0.0);
        let mut lower = source.regions()[0].card().clone();
        lower.validity = ValidityDomain::unconstrained().with("T", 200.0, 1000.0);
        let mut upper = lower.clone();
        upper.validity = ValidityDomain::unconstrained().with("T", 1000.0 + gap, 2000.0);
        upper.parameters.get_mut("a7").unwrap().value -= 4000.0;
        let model = Nasa9StandardStateModelV1::new(
            source.species().clone(),
            source.molar_mass(),
            source.convention().clone(),
            vec![
                Nasa9RegionV1::from_card(lower).unwrap(),
                Nasa9RegionV1::from_card(upper).unwrap(),
            ],
        )
        .unwrap();
        let (bank, chemistry) = fixture(
            vec![model],
            vec![transport("A", 1.7e-5, 0.025)],
            &["H"],
            vec![vec![1]],
        );
        let mut c = controls();
        c.equilibrium.maximum_temperature_k = 800.0;
        let c = GasEquilibriumChamberSolve {
            equilibrium: c.equilibrium,
            mass_tolerance_kg: c.mass_tolerance_kg,
        };
        let context = point(300.0, 1e5);
        let target = chemistry
            .feed_enthalpy(Temperature::new(1500.0), &feed(), || false)
            .unwrap();
        // The explicit target is supported entirely on its declared lower
        // chart; the single-temperature inlet may not jump between raw charts.
        let explicit = bank
            .equilibrium_chamber_hp_enthalpy(&chemistry, &context, &feed(), target, c, || false)
            .unwrap();
        close(explicit.parcel().state().state().temperature, 500.0, 1e-10);
        assert!(matches!(
            bank.equilibrium_chamber_hp(
                &chemistry,
                &context,
                &feed(),
                Temperature::new(1500.0),
                c,
                || false,
            ),
            Err(GasEquilibriumChamberError::Chemistry(
                CompositionEquilibriumErrorV1::Invalid(_)
            ))
        ));
    }
}
