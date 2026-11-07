use super::*;
use crate::solution::{BinaryCoexistenceOptionsV1, BinaryRegularSolutionV1};
use crate::{MolarThermalQuantityV1, SpeciesId};

fn model(a: f64, b: f64, excess_cp: f64) -> BinaryLiquidCaloricModelV1 {
    let pair = BinaryRegularSolutionV1::new(
        [SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()],
        a,
        b,
        Temperature::new(100.0),
        Temperature::new(1000.0),
    )
    .unwrap()
    .with_excess_heat_capacity(excess_cp, Temperature::new(400.0))
    .unwrap();
    let capacities = [
        ("A", [75.0, 115.0, 95.0, 125.0]),
        ("B", [150.0, 170.0, 210.0, 180.0]),
    ]
    .map(|(species, cp)| {
        PureLiquidHeatCapacityV1::new(
            SpeciesId::new(species).unwrap(),
            Temperature::new(100.0),
            [100.0, 350.0, 650.0, 1000.0]
                .into_iter()
                .zip(cp)
                .map(|(t, c)| (Temperature::new(t), MolarThermalQuantityV1::new(c)))
                .collect(),
        )
        .unwrap()
    });
    BinaryLiquidCaloricModelV1::new(pair, capacities).unwrap()
}

fn point(t: f64, x: f64) -> BinaryLiquidTransportPointV1 {
    BinaryLiquidTransportPointV1 {
        temperature: Temperature::new(t),
        mole_fractions: [x, 1.0 - x],
    }
}

fn coefficients(q: f64) -> BinaryLiquidTransportCoefficientsV1 {
    BinaryLiquidTransportCoefficientsV1 {
        diffusivity_m2_s: 3.0e-9,
        thermal_conductivity_w_m_k: 0.25,
        heat_of_transport_j_per_mol: q,
    }
}

fn controls() -> BinaryLiquidTransportSolveV1 {
    BinaryLiquidTransportSolveV1 {
        temperature_bracket: [Temperature::new(200.0), Temperature::new(800.0)],
        temperature_tolerance_k: 1e-9,
        energy_tolerance_j: 1e-7,
        maximum_iterations: 96,
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:.17e}, expected {expected:.17e}, tolerance {tolerance:e}"
    );
}

fn face(
    model: &BinaryLiquidCaloricModelV1,
    points: [BinaryLiquidTransportPointV1; 2],
    coefficients: BinaryLiquidTransportCoefficientsV1,
) -> BinaryLiquidTransportFaceV1 {
    model
        .transport_face(
            points,
            AmountConcentration::new(40_000.0),
            2e-4,
            coefficients,
            || false,
        )
        .unwrap()
}

// Independent endpoint chemical potential, using the original pure H/S and
// the thermodynamic model's chemical potentials; no face-helper arithmetic.
fn chemical_potential_over_t(
    model: &BinaryLiquidCaloricModelV1,
    point: BinaryLiquidTransportPointV1,
) -> f64 {
    let mut pure = [0.0; 2];
    for (i, source) in model.heat_capacities().iter().enumerate() {
        pure[i] = source
            .enthalpy_change_j_per_mol(source.reference_temperature(), point.temperature)
            .unwrap()
            / point.temperature.value()
            - source
                .entropy_change_j_per_mol_k(source.reference_temperature(), point.temperature)
                .unwrap();
    }
    let state = model
        .thermodynamics()
        .evaluate_composition(point.temperature, point.mole_fractions)
        .unwrap();
    let reduced = state.reduced_chemical_potentials();
    pure[0] - pure[1] + R * (reduced[0].unwrap() - reduced[1].unwrap())
}

#[test]
fn g1_uniform_composition_has_no_spurious_caloric_thermal_diffusion() {
    // Both heat capacities vary across multiple source segments. Signed excess
    // Cp changes chi(T), partial enthalpies, and the transported energy chart.
    for excess_cp in [-40.0, 40.0] {
        let model = model(0.25, 120.0, excess_cp);
        for temperatures in [[230.0, 720.0], [350.0, 350.0_f64.next_up()]] {
            let points = temperatures.map(|t| point(t, 0.3));
            let result = face(&model, points, coefficients(0.0));
            assert_eq!(result.chemical_effort_j_per_mol_k, 0.0);
            assert_eq!(result.molar_flux_mol_m2_s, [0.0, 0.0]);
            let fourier = -0.25 * ((temperatures[1] - temperatures[0]) / 2e-4);
            assert_eq!(result.energy_flux_w_m2, fourier);
            assert_eq!(result.reduced_heat_flux_w_m2, fourier);
            assert!(result.entropy_production_w_m2_k > 0.0);
        }
    }
}

#[test]
fn g1_original_calorics_close_the_nonisothermal_onsager_entropy_identity() {
    for excess_cp in [-40.0, 40.0] {
        let model = model(0.25, 120.0, excess_cp);
        let points = [point(260.0, 0.23), point(730.0, 0.61)];
        for q in [-2500.0, 0.0, 2500.0] {
            let result = face(&model, points, coefficients(q));
            let j = result.molar_flux_mol_m2_s[0];
            let delta_beta = 1.0 / 730.0 - 1.0 / 260.0;
            let full_effort = chemical_potential_over_t(&model, points[1])
                - chemical_potential_over_t(&model, points[0]);
            near(
                result.chemical_effort_j_per_mol_k
                    + result.exchange_enthalpy_j_per_mol * delta_beta,
                full_effort,
                8e-14,
            );
            near(
                result.entropy_production_w_m2_k,
                -j * full_effort + result.energy_flux_w_m2 * delta_beta,
                3e-12,
            );
            assert_eq!(j, -result.molar_flux_mol_m2_s[1]);
            assert!(result.entropy_production_w_m2_k > 0.0);
            let reverse = face(&model, [points[1], points[0]], coefficients(q));
            assert_eq!(reverse.molar_flux_mol_m2_s, [-j, j]);
            assert_eq!(reverse.energy_flux_w_m2, -result.energy_flux_w_m2);
            assert_eq!(
                reverse.exchange_enthalpy_j_per_mol,
                result.exchange_enthalpy_j_per_mol
            );
            assert_eq!(
                reverse.entropy_production_w_m2_k,
                result.entropy_production_w_m2_k
            );
        }
    }
}

#[test]
fn g1_soret_and_dufour_use_the_same_authored_heat_of_transport() {
    let model = model(0.25, 120.0, 40.0);
    let mut kinetic = coefficients(2400.0);
    kinetic.thermal_conductivity_w_m_k = 0.0;
    let thermal = face(&model, [point(300.0, 0.3), point(600.0, 0.3)], kinetic);
    let expected_j = thermal.mobility_mol2_k_j_m_s
        * kinetic.heat_of_transport_j_per_mol
        * (1.0 / 600.0 - 1.0 / 300.0)
        / 2e-4;
    near(
        thermal.molar_flux_mol_m2_s[0],
        expected_j,
        4.0 * f64::EPSILON * expected_j.abs(),
    );
    assert!(thermal.molar_flux_mol_m2_s[0] < 0.0);
    for points in [
        [point(300.0, 0.3), point(600.0, 0.3)],
        [point(450.0, 0.3), point(450.0, 0.7)],
    ] {
        let result = face(&model, points, kinetic);
        let j = result.molar_flux_mol_m2_s[0];
        near(
            result.reduced_heat_flux_w_m2,
            kinetic.heat_of_transport_j_per_mol * j,
            1e-15,
        );
        let expected_energy =
            (result.exchange_enthalpy_j_per_mol + kinetic.heat_of_transport_j_per_mol) * j;
        near(
            result.energy_flux_w_m2,
            expected_energy,
            8.0 * f64::EPSILON * expected_energy.abs(),
        );
        assert!(result.entropy_production_w_m2_k > 0.0);
    }
}

#[test]
fn g1_ideal_fick_limit_and_critical_cubic_response_survive_small_gradients() {
    let ideal = model(0.0, 0.0, 0.0);
    let points = [point(450.0, 0.3), point(450.0, 0.300001)];
    let result = face(&ideal, points, coefficients(0.0));
    let expected =
        -40_000.0 * 3e-9 * (points[1].mole_fractions[0] - points[0].mole_fractions[0]) / 2e-4;
    near(
        result.molar_flux_mol_m2_s[0],
        expected,
        expected.abs() * 1e-10,
    );

    let critical = model(2.0, 0.0, 0.0);
    let delta = 2.0_f64.powi(-25);
    let result = face(
        &critical,
        [point(450.0, 0.5 - delta), point(450.0, 0.5 + delta)],
        coefficients(0.0),
    );
    // ln[(1/2+d)/(1/2-d)] - 4d = 16d^3/3 + O(d^5),
    // independently evaluated at the two symmetric endpoints.
    let expected = R * (32.0 / 3.0) * delta.powi(3);
    near(
        result.chemical_effort_j_per_mol_k,
        expected,
        expected * 2e-14,
    );
    assert!(result.chemical_effort_j_per_mol_k > 0.0);
    assert!(result.molar_flux_mol_m2_s[0] < 0.0);
    assert!(result.entropy_production_w_m2_k > 0.0);
}

#[test]
fn g1_homogeneous_inventory_changes_carry_original_species_and_mixing_calorics() {
    for excess_cp in [-40.0, 40.0] {
        let model = model(0.25, 120.0, excess_cp);
        let initial = model
            .homogeneous(Temperature::new(270.0), [2.0, 3.0])
            .unwrap();
        // Deliberately also change total moles: the helper must carry the
        // extensive mixing term, not silently assume equal-volume exchange.
        let final_state = model
            .homogeneous(Temperature::new(740.0), [3.4, 2.1])
            .unwrap();
        let change = model
            .homogeneous_transport_change(&initial, &final_state)
            .unwrap();
        near(
            change.enthalpy_change_j,
            final_state.enthalpy_j() - initial.enthalpy_j(),
            1e-10,
        );
        near(
            change.entropy_change_j_per_k,
            final_state.reduced_entropy_j_per_k().unwrap()
                - initial.reduced_entropy_j_per_k().unwrap(),
            5e-13,
        );
        let reverse = model
            .homogeneous_transport_change(&final_state, &initial)
            .unwrap();
        near(change.enthalpy_change_j, -reverse.enthalpy_change_j, 1e-10);
        near(
            change.entropy_change_j_per_k,
            -reverse.entropy_change_j_per_k,
            5e-13,
        );
    }
    // Large source references erase this thermal increment from absolute H/S.
    // The original local integration must nevertheless retain its sign/size.
    let large = model(-1e18, 1e18, 80.0);
    let t = 350.0_f64;
    let initial = large.homogeneous(Temperature::new(t), [2.0, 3.0]).unwrap();
    let final_state = large
        .homogeneous(Temperature::new(t.next_up()), [2.0, 3.0])
        .unwrap();
    assert_eq!(initial.enthalpy_j(), final_state.enthalpy_j());
    let change = large
        .homogeneous_transport_change(&initial, &final_state)
        .unwrap();
    let cp = 2.0 * 115.0 + 3.0 * 170.0 + 5.0 * 80.0 * 0.4 * 0.6;
    near(change.enthalpy_change_j, cp * (t.next_up() - t), 1e-25);
    near(
        change.entropy_change_j_per_k,
        cp * ((t.next_up() - t) / t),
        1e-27,
    );
}

#[test]
fn g1_inventory_enthalpy_inverse_recovers_temperature_without_a_phase_flash() {
    let model = model(0.25, 120.0, 40.0);
    let initial = model
        .homogeneous(Temperature::new(310.0), [4.0, 7.0])
        .unwrap();
    let target = model
        .homogeneous(Temperature::new(523.456), [4.3, 6.7])
        .unwrap();
    let heat = target.enthalpy_j() - initial.enthalpy_j();
    let result = model
        .homogeneous_transport_candidate(&initial, [4.3, 6.7], heat, controls(), || false)
        .unwrap();
    assert_eq!(result.final_state().component_moles(), [4.3, 6.7]);
    assert!(result.final_state().equilibrium().is_none());
    near(result.final_state().temperature().value(), 523.456, 1e-9);
    near(result.energy_residual_j(), 0.0, 1e-7);
    near(
        result.change().entropy_change_j_per_k,
        target.reduced_entropy_j_per_k().unwrap() - initial.reduced_entropy_j_per_k().unwrap(),
        2e-10,
    );
    let reverse = model
        .homogeneous_transport_candidate(
            result.final_state(),
            [4.0, 7.0],
            -heat,
            controls(),
            || false,
        )
        .unwrap();
    near(reverse.final_state().temperature().value(), 310.0, 2e-9);
    assert_eq!(initial.component_moles(), [4.0, 7.0]);
    assert_eq!(initial.temperature().value(), 310.0);
}

#[test]
fn g3_convex_admission_checks_the_interior_excess_capacity_maximum() {
    let bracket = [Temperature::new(250.0), Temperature::new(650.0)];
    let unsafe_model = model(1.9, 80.0, 40.0);
    assert!(
        bracket
            .into_iter()
            .all(|t| unsafe_model.thermodynamics().chi(t).unwrap() < 2.0)
    );
    let interior = Temperature::new(400.0 - R * 80.0 / 40.0);
    assert!(unsafe_model.thermodynamics().chi(interior).unwrap() > 2.0);
    assert!(
        unsafe_model
            .admit_transport_temperature_interval(bracket)
            .is_err()
    );
    let admitted = model(1.79, 80.0, 40.0);
    near(
        admitted
            .admit_transport_temperature_interval(bracket)
            .unwrap(),
        admitted.thermodynamics().chi(interior).unwrap(),
        2e-15,
    );
    assert_eq!(
        model(2.0, 0.0, 0.0)
            .admit_transport_temperature_interval(bracket)
            .unwrap(),
        2.0
    );
}

#[test]
fn g3_transport_refusals_are_atomic_and_zero_kinetics_still_exposes_thermodynamics() {
    let model = model(0.25, 120.0, 40.0);
    let points = [point(300.0, 0.2), point(600.0, 0.6)];
    let zero = BinaryLiquidTransportCoefficientsV1 {
        diffusivity_m2_s: 0.0,
        thermal_conductivity_w_m_k: 0.0,
        heat_of_transport_j_per_mol: 0.0,
    };
    let result = face(&model, points, zero);
    assert_ne!(result.chemical_effort_j_per_mol_k, 0.0);
    assert_ne!(result.exchange_enthalpy_j_per_mol, 0.0);
    assert_eq!(result.molar_flux_mol_m2_s, [0.0; 2]);
    assert_eq!(result.energy_flux_w_m2, 0.0);
    assert_eq!(result.entropy_production_w_m2_k, 0.0);
    for invalid in [
        BinaryLiquidTransportCoefficientsV1 {
            diffusivity_m2_s: -1.0,
            ..zero
        },
        BinaryLiquidTransportCoefficientsV1 {
            thermal_conductivity_w_m_k: -1.0,
            ..zero
        },
        BinaryLiquidTransportCoefficientsV1 {
            heat_of_transport_j_per_mol: f64::NAN,
            ..zero
        },
    ] {
        assert!(
            model
                .transport_face(
                    points,
                    AmountConcentration::new(40_000.0),
                    2e-4,
                    invalid,
                    || false
                )
                .is_err()
        );
    }
    for invalid in [point(300.0, 0.0), point(300.0, 1.0), point(99.0, 0.2)] {
        assert!(
            model
                .transport_face(
                    [invalid, points[1]],
                    AmountConcentration::new(40_000.0),
                    2e-4,
                    zero,
                    || false
                )
                .is_err()
        );
    }
    let initial = model
        .homogeneous(Temperature::new(300.0), [2.0, 3.0])
        .unwrap();
    let before = initial.clone();
    let mut polls = 0;
    let cancelled =
        model.homogeneous_transport_candidate(&initial, [2.1, 2.9], 400.0, controls(), || {
            polls += 1;
            polls == 3
        });
    assert_eq!(cancelled.unwrap_err(), BinaryCaloricErrorV1::Cancelled);
    assert!(matches!(
        model.homogeneous_transport_candidate(&initial, [2.1, 2.9], 1e8, controls(), || false),
        Err(BinaryCaloricErrorV1::HeatOutsideBracket { .. })
    ));
    let equilibrium = model
        .equilibrate(
            Temperature::new(300.0),
            [2.0, 3.0],
            BinaryCoexistenceOptionsV1::default(),
            || false,
        )
        .unwrap();
    assert!(
        model
            .homogeneous_transport_change(&initial, &equilibrium)
            .is_err()
    );
    assert_eq!(initial, before);
}

#[test]
fn g1_trace_mobility_keeps_independent_positive_component_support() {
    let model = model(0.0, 0.0, 0.0);
    let points = [1e-300, 2e-300].map(|trace| BinaryLiquidTransportPointV1 {
        temperature: Temperature::new(400.0),
        mole_fractions: [trace, 1.0],
    });
    let result = model
        .transport_face(
            points,
            AmountConcentration::new(1e300),
            1.0,
            BinaryLiquidTransportCoefficientsV1 {
                diffusivity_m2_s: 1e-9,
                thermal_conductivity_w_m_k: 0.0,
                heat_of_transport_j_per_mol: 0.0,
            },
            || false,
        )
        .unwrap();
    let expected = -1.5e-9 * 2.0_f64.ln();
    near(
        result.molar_flux_mol_m2_s[0],
        expected,
        expected.abs() * 2e-15,
    );
    assert!(result.entropy_production_w_m2_k > 0.0);
}
