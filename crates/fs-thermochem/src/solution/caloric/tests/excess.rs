//! Whole-owner caloric checks for signed heat capacities of mixing.

use super::*;

fn pair(a: f64, b: f64, c: f64, reference: f64) -> BinaryRegularSolutionV1 {
    BinaryRegularSolutionV1::new(
        [SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()],
        a,
        b,
        Temperature::new(100.0),
        Temperature::new(1000.0),
    )
    .unwrap()
    .with_excess_heat_capacity(c, Temperature::new(reference))
    .unwrap()
}

fn caloric(a: f64, b: f64, c: f64, reference: f64) -> BinaryLiquidCaloricModelV1 {
    BinaryLiquidCaloricModelV1::new(
        pair(a, b, c, reference),
        [source("A", 100.0), source("B", 160.0)],
    )
    .unwrap()
}

fn equilibrium(model: &BinaryLiquidCaloricModelV1, t: f64) -> BinaryCaloricStateV1 {
    model
        .equilibrate(
            Temperature::new(t),
            [2.0, 3.0],
            BinaryCoexistenceOptionsV1 {
                relative_tolerance: 4.0 * f64::EPSILON,
                maximum_iterations: 256,
            },
            || false,
        )
        .unwrap()
}

#[test]
fn g1_source_reduced_entropy_retains_signed_excess_and_accepted_coexistence() {
    for c in [-40.0, 40.0] {
        let model = caloric(0.0, 100.0, c, 400.0);
        let pure = model
            .equilibrate(
                Temperature::new(100.0),
                [1.0, 0.0],
                controls().coexistence,
                || false,
            )
            .unwrap();
        assert_eq!(pure.reduced_entropy_j_per_k().unwrap(), 0.0);
        let ideal = -5.0 * R * (0.4_f64 * 0.4_f64.ln() + 0.6_f64 * 0.6_f64.ln());
        for t in [350.0_f64, 450.0] {
            let state = equilibrium(&model, t);
            assert_eq!(state.equilibrium().unwrap().phases().len(), 1);
            let excess = 5.0 * 0.4 * 0.6 * c * (t / 400.0).ln();
            near(
                state.reduced_entropy_j_per_k().unwrap(),
                680.0 * (t / 100.0).ln() + ideal + excess,
                2e-12,
            );
            assert_eq!(
                excess.is_sign_positive(),
                (c * (t - 400.0)).is_sign_positive()
            );
        }
        let cold = equilibrium(&model, 350.0);
        let hot = equilibrium(&model, 450.0);
        near(
            hot.reduced_entropy_j_per_k().unwrap() - cold.reduced_entropy_j_per_k().unwrap(),
            model.entropy_change_j_per_k(&cold, &hot).unwrap(),
            2e-12,
        );

        // At T=T0=400 K, chi=1200/400=3 independently of signed C.
        // Use the independently solved dilute binodal, not the feed x_A=.4.
        let separated = caloric(0.0, 1200.0, c, 400.0);
        let state = equilibrium(&separated, 400.0);
        assert_eq!(state.equilibrium().unwrap().phases().len(), 2);
        let x = 0.070_720_181_679_944_8_f64;
        let phase = -5.0 * R * (x * x.ln() + (1.0 - x) * (1.0 - x).ln());
        near(
            state.reduced_entropy_j_per_k().unwrap(),
            680.0 * 4.0_f64.ln() + phase,
            2e-11,
        );
        assert!((phase - ideal).abs() > 10.0);
    }
}

#[test]
fn g1_signed_homogeneous_mixing_cp_closes_exact_enthalpy_and_entropy() {
    for c in [-80.0, 80.0] {
        let model = caloric(0.0, 0.0, c, 400.0);
        let cold = equilibrium(&model, 350.0);
        let hot = equilibrium(&model, 450.0);
        assert_eq!(cold.equilibrium().unwrap().phases().len(), 1);
        assert_eq!(hot.equilibrium().unwrap().phases().len(), 1);
        let cp = 2.0 * 100.0 + 3.0 * 160.0 + 5.0 * c * 0.4 * 0.6;
        near(cold.isobaric_heat_capacity_j_per_k().unwrap(), cp, 2e-13);
        near(
            model.enthalpy_change_j(&cold, &hot).unwrap(),
            cp * 100.0,
            2e-11,
        );
        near(
            model.entropy_change_j_per_k(&cold, &hot).unwrap(),
            cp * (450.0_f64 / 350.0).ln(),
            2e-13,
        );
        assert_eq!(
            model.enthalpy_change_j(&cold, &hot).unwrap(),
            -model.enthalpy_change_j(&hot, &cold).unwrap()
        );
        assert_eq!(
            model.entropy_change_j_per_k(&cold, &hot).unwrap(),
            -model.entropy_change_j_per_k(&hot, &cold).unwrap()
        );
        let balance = model
            .add_heat(&cold, cp * 100.0, controls(), || false)
            .unwrap();
        near(balance.final_state().temperature().value(), 450.0, 1e-9);
        assert!(balance.energy_residual_j().abs() < 1e-7);
    }
}

#[test]
fn g1_equilibrium_cp_includes_signed_mixing_and_phase_redistribution() {
    for c in [-20.0, 80.0] {
        let model = caloric(0.0, 1000.0, c, 400.0);
        for t in [350.0, 400.0, 450.0] {
            let state = equilibrium(&model, t);
            assert_eq!(state.equilibrium().unwrap().phases().len(), 2);
            let step = 0.002;
            let low = equilibrium(&model, t - step);
            let high = equilibrium(&model, t + step);
            let cp = state.isobaric_heat_capacity_j_per_k().unwrap();
            near(
                cp,
                model.enthalpy_change_j(&low, &high).unwrap() / (2.0 * step),
                5e-6,
            );
            near(
                cp,
                t * model.entropy_change_j_per_k(&low, &high).unwrap() / (2.0 * step),
                5e-6,
            );
            let x = state.equilibrium().unwrap().phases()[0]
                .state()
                .mole_fractions()[0];
            let frozen = 680.0 + 5.0 * c * x * (1.0 - x);
            assert!(cp > frozen, "equilibrium redistribution must be included");
        }
        let cold = equilibrium(&model, 350.0);
        let hot = equilibrium(&model, 650.0);
        assert_eq!(hot.equilibrium().unwrap().phases().len(), 1);
        let heat = model.enthalpy_change_j(&cold, &hot).unwrap();
        let forward = model.add_heat(&cold, heat, controls(), || false).unwrap();
        near(forward.final_state().temperature().value(), 650.0, 1e-9);
        let reverse = model
            .add_heat(forward.final_state(), -heat, controls(), || false)
            .unwrap();
        near(reverse.final_state().temperature().value(), 350.0, 1e-9);
        near(
            forward.entropy_change_j_per_k(),
            -reverse.entropy_change_j_per_k(),
            2e-9,
        );
    }
}

#[test]
fn g1_mixing_cp_retains_small_heat_against_large_interaction_references() {
    let model = caloric(-1e18, 1e18, 80.0, 300.0);
    let a = 350.0_f64;
    let b = a.next_up();
    let low = model.homogeneous(Temperature::new(a), [2.0, 3.0]).unwrap();
    let high = model.homogeneous(Temperature::new(b), [2.0, 3.0]).unwrap();
    // Absolute source enthalpies cannot resolve this heat; local increments
    // still recover the original pure Cp plus the authored excess Cp exactly.
    assert_eq!(low.enthalpy_j(), high.enthalpy_j());
    let cp = 680.0 + 5.0 * 80.0 * 0.4 * 0.6;
    near(
        model.enthalpy_change_j(&low, &high).unwrap(),
        cp * (b - a),
        1e-25,
    );
    near(
        model.entropy_change_j_per_k(&low, &high).unwrap(),
        cp * ((b - a) / a),
        1e-27,
    );
    assert!(model.entropy_change_j_per_k(&low, &high).unwrap() > 0.0);
}

#[test]
fn g1_critical_response_keeps_mixing_cp_and_declares_its_one_sided_tangent() {
    let model = BinaryLiquidCaloricModelV1::new(
        pair(0.0, 1000.0, 80.0, 500.0),
        [source("A", 100.0), source("B", 100.0)],
    )
    .unwrap();
    let critical = model
        .equilibrate(
            Temperature::new(500.0),
            [1.0, 1.0],
            controls().coexistence,
            || false,
        )
        .unwrap();
    assert!(critical.isobaric_heat_capacity_j_per_k().is_err());
    let response = critical.isobaric_heat_capacity_response().unwrap();
    assert_eq!(response.heat_capacity_j_per_k(), 240.0);
    assert_eq!(
        response.qualification(),
        BinaryCaloricHeatCapacityQualificationV1::HomogeneousPhaseBoundaryTangent
    );
    let below = model
        .equilibrate(
            Temperature::new(500.0_f64.next_down()),
            [1.0, 1.0],
            controls().coexistence,
            || false,
        )
        .unwrap();
    near(
        below.isobaric_heat_capacity_j_per_k().unwrap(),
        240.0 + 3.0 * R,
        3e-12,
    );
    // At T0 the excess enthalpy vanishes when B=0: a critical composition
    // then has a unique Cp despite the nonzero excess-Cp coefficient itself.
    let stationary = BinaryLiquidCaloricModelV1::new(
        pair(2.0, 0.0, 80.0, 500.0),
        [source("A", 100.0), source("B", 100.0)],
    )
    .unwrap();
    let critical = stationary
        .equilibrate(
            Temperature::new(500.0),
            [1.0, 1.0],
            controls().coexistence,
            || false,
        )
        .unwrap();
    assert_eq!(critical.isobaric_heat_capacity_j_per_k().unwrap(), 240.0);
}

#[test]
fn g3_negative_mixing_cp_admission_checks_composition_and_internal_source_knots() {
    // Exact composition-wide condition: |C| < (sqrt(Cp_A)+sqrt(Cp_B))².
    // This valid asymmetric pair would fail the crude |C| < 4*min(Cp) test.
    let valid = BinaryLiquidCaloricModelV1::new(
        pair(0.0, 0.0, -120.0, 400.0),
        [source("A", 1.0), source("B", 100.0)],
    )
    .unwrap();
    assert!(
        valid
            .homogeneous(Temperature::new(400.0), [10.0, 1.0])
            .unwrap()
            .isobaric_heat_capacity_j_per_k()
            .unwrap()
            > 0.0
    );
    for c in [-121.0, -122.0, -f64::MAX] {
        assert!(
            BinaryLiquidCaloricModelV1::new(
                pair(0.0, 0.0, c, 400.0),
                [source("A", 1.0), source("B", 100.0)]
            )
            .is_err()
        );
    }
    let dipped = |species| {
        PureLiquidHeatCapacityV1::new(
            SpeciesId::new(species).unwrap(),
            Temperature::new(100.0),
            [(100.0, 100.0), (400.0, 1.0), (1000.0, 100.0)]
                .into_iter()
                .map(|(t, cp)| (Temperature::new(t), MolarThermalQuantityV1::new(cp)))
                .collect(),
        )
        .unwrap()
    };
    assert!(
        BinaryLiquidCaloricModelV1::new(pair(0.0, 0.0, -10.0, 400.0), [dipped("A"), dipped("B")])
            .is_err()
    );
}

#[test]
fn g1_storage_bound_covers_negative_mixing_cp_and_actual_equilibrium_response() {
    let domain = [Temperature::new(300.0), Temperature::new(700.0)];
    for c in [-80.0, 80.0] {
        let model = caloric(0.0, 1000.0, c, 400.0);
        let bound = model
            .isobaric_heat_capacity_lower_bound_j_per_k([2.0, 3.0], domain)
            .unwrap();
        let expected = 680.0 + 1.2 * c.min(0.0);
        assert!(bound <= expected && bound > expected * (1.0 - 3e-15));
        for t in [300.0, 350.0, 400.0, 450.0, 600.0, 700.0] {
            assert!(
                equilibrium(&model, t)
                    .isobaric_heat_capacity_response()
                    .unwrap()
                    .heat_capacity_j_per_k()
                    >= bound
            );
        }
        let pure = model
            .isobaric_heat_capacity_lower_bound_j_per_k([1.0, 0.0], domain)
            .unwrap();
        assert!(pure <= 100.0 && pure > 100.0 * (1.0 - 1e-15));
        assert!(
            model
                .isobaric_heat_capacity_lower_bound_j_per_k([0.0, 0.0], domain)
                .is_err()
        );
        assert!(
            model
                .isobaric_heat_capacity_lower_bound_j_per_k(
                    [2.0, 3.0],
                    [Temperature::new(99.0), domain[1]]
                )
                .is_err()
        );
    }
}
