//! Original equilibrium calorics drive the native nonlinear thermal owner.

use super::*;
use fs_material::liquid::solution::caloric::storage::conduction::BinaryLiquidConductionLaw;
use fs_material::thermal::conduction::ConductivityCurve;
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionChain, ImplicitConductionControls, ImplicitConductionLaw,
    ImplicitConductionWorkspace,
};
use fs_thermal::enthalpy::conduction::{ConductionObserver, ConductiveLink, ConductivityLaw};

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}"
    );
}

fn conductivity() -> ConductivityCurve {
    ConductivityCurve::new(vec![[290.0, 4.0], [400.0, 9.5], [600.0, 8.0]]).unwrap()
}

#[test]
fn equilibrium_specific_heat_is_intensive_and_matches_original_enthalpy_derivative() {
    let card = card(CompositionBasis::MoleFraction, 1.0);
    let liquid = liquid(&card, 900.0);
    let thermal = BinaryLiquidCaloricLaw::try_new(&liquid, sources(), "pressure").unwrap();
    let storage = thermal
        .thermal_storage([2.0, 3.0], controls(), || false)
        .unwrap();
    let doubled = thermal
        .thermal_storage([4.0, 6.0], controls(), || false)
        .unwrap();
    let sensible_cp = 500.0 / 0.336;
    // Independent high-precision solution of the original binodal equation
    // and its analytic derivative, rounded only for these f64 assertions.
    for (t, expected) in [
        (300.0, 1577.156_735_939_303_2),
        (350.0, 1608.580_639_961_647_8),
        (430.0, 1660.515_075_290_422_4),
        (500.0, sensible_cp),
    ] {
        let cp = storage
            .specific_heat_cp_j_kg_k(Temperature::new(t), || false)
            .unwrap();
        close(cp, expected, 2e-8);
        close(
            cp,
            doubled
                .specific_heat_cp_j_kg_k(Temperature::new(t), || false)
                .unwrap(),
            2e-11,
        );
        let dt = 0.005;
        let enthalpy = |t| {
            storage
                .specific_enthalpy_j_kg(Temperature::new(t), || false)
                .unwrap()
        };
        close(cp, (enthalpy(t + dt) - enthalpy(t - dt)) / (2.0 * dt), 2e-5);
        if t < 450.0 {
            assert!(cp > sensible_cp);
        } else {
            close(cp, sensible_cp, 2e-12);
        }
    }
    assert!(matches!(
        storage.specific_heat_cp_j_kg_k(Temperature::new(350.0), || true),
        Err(BinaryLiquidCaloricError::Cancelled)
    ));
    for t in [289.0, 601.0, f64::NAN] {
        assert!(
            storage
                .specific_heat_cp_j_kg_k(Temperature::new(t), || false)
                .is_err()
        );
    }
}

#[test]
fn equilibrium_conduction_shares_one_caloric_owner_and_analytic_local_tangents() {
    let card = card(CompositionBasis::MoleFraction, 1.0);
    let liquid = liquid(&card, 900.0);
    let thermal = BinaryLiquidCaloricLaw::try_new(&liquid, sources(), "pressure").unwrap();
    let storage = thermal
        .thermal_storage([2.0, 3.0], controls(), || false)
        .unwrap();
    let law = BinaryLiquidConductionLaw::try_new(&storage, conductivity()).unwrap();
    assert!(core::ptr::eq(law.caloric(), &storage));
    assert_eq!(
        law.enthalpy_bounds_j_kg(),
        storage.specific_enthalpy_domain_j_kg()
    );
    assert_eq!(law.conductivity(), &conductivity());
    for t in [350.0, 430.0, 500.0] {
        let h = storage
            .specific_enthalpy_j_kg(Temperature::new(t), || false)
            .unwrap();
        let state = law.state_at_specific_enthalpy(h).unwrap();
        let cp = storage
            .specific_heat_cp_j_kg_k(Temperature::new(t), || false)
            .unwrap();
        let (k, dk_dt) = if t < 400.0 {
            (4.0 + 0.05 * (t - 290.0), 0.05)
        } else {
            (9.5 - 0.0075 * (t - 400.0), -0.0075)
        };
        close(state.temperature_k, t, 2e-9);
        close(state.temperature_slope_kg_k_j, 1.0 / cp, 2e-14);
        close(state.conductivity_w_m_k, k, 2e-10);
        close(state.conductivity_slope_kg_w_m_k_j, dk_dt / cp, 2e-15);
        close(law.conductivity_w_m_k(h).unwrap(), k, 2e-10);
        let delta = 0.1;
        let lower = law.state_at_specific_enthalpy(h - delta).unwrap();
        let upper = law.state_at_specific_enthalpy(h + delta).unwrap();
        close(
            state.temperature_slope_kg_k_j,
            (upper.temperature_k - lower.temperature_k) / (2.0 * delta),
            2e-9,
        );
        close(
            state.conductivity_slope_kg_w_m_k_j,
            (upper.conductivity_w_m_k - lower.conductivity_w_m_k) / (2.0 * delta),
            2e-10,
        );
        assert_eq!(law.time_observables(h).unwrap().phase_mass_fraction, None);
    }
    let curve = conductivity();
    close(
        curve
            .conductivity_and_temperature_derivative(290.0)
            .unwrap()
            .1,
        0.05,
        1e-17,
    );
    close(
        curve
            .conductivity_and_temperature_derivative(400.0)
            .unwrap()
            .1,
        -0.0075,
        1e-17,
    );
    close(
        curve
            .conductivity_and_temperature_derivative(600.0)
            .unwrap()
            .1,
        -0.0075,
        1e-17,
    );
    assert!(
        BinaryLiquidConductionLaw::try_new(
            &storage,
            ConductivityCurve::new(vec![[300.0, 1.0], [600.0, 1.0]]).unwrap()
        )
        .is_err()
    );
    for h in [
        -1.0,
        storage.specific_enthalpy_domain_j_kg()[1].next_up(),
        f64::NAN,
    ] {
        assert!(law.state_at_specific_enthalpy(h).is_err());
    }
}

#[test]
fn critical_corner_uses_explicit_one_sided_tangent_without_claiming_unique_cp() {
    let card = card(CompositionBasis::MoleFraction, 1.0);
    let liquid = liquid(&card, 900.0);
    let thermal = BinaryLiquidCaloricLaw::try_new(&liquid, sources(), "pressure").unwrap();
    let storage = thermal
        .thermal_storage(
            [1.0, 1.0],
            BinaryCaloricSolveV1 {
                temperature_bracket: [Temperature::new(450.0), Temperature::new(600.0)],
                ..controls()
            },
            || false,
        )
        .unwrap();
    assert!(
        storage
            .specific_heat_cp_j_kg_k(Temperature::new(450.0), || false)
            .is_err()
    );
    let law = BinaryLiquidConductionLaw::try_new(&storage, conductivity()).unwrap();
    let state = law.state_at_specific_enthalpy(0.0).unwrap();
    close(state.temperature_k, 450.0, 0.0);
    close(state.temperature_slope_kg_k_j, 0.118 / 200.0, 1e-18);
}

#[test]
fn native_implicit_conduction_crosses_miscibility_and_matches_independent_heat_balance() {
    let card = card(CompositionBasis::MoleFraction, 1.0);
    let liquid = liquid(&card, 900.0);
    let thermal = BinaryLiquidCaloricLaw::try_new(&liquid, sources(), "pressure").unwrap();
    let storage = thermal
        .thermal_storage([2.0, 3.0], controls(), || false)
        .unwrap();
    let law = BinaryLiquidConductionLaw::try_new(
        &storage,
        ConductivityCurve::new(vec![[290.0, 10.0], [600.0, 10.0]]).unwrap(),
    )
    .unwrap();
    let h = [350.0, 590.0].map(|t| {
        storage
            .specific_enthalpy_j_kg(Temperature::new(t), || false)
            .unwrap()
    });
    let initial = DepositionSimulation::new(
        vec![
            birth(1, 0.336, h[0], 0.0),
            birth(2, 0.672, h[1], 0.0),
            TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: 1,
                    a: CellId(1),
                    b: CellId(2),
                    area_m2: 0.1,
                    conductance_w_k: 0.0,
                }),
            },
        ],
        &[(MaterialId(1), &storage)],
        budget(),
    )
    .unwrap();
    let (initial, _) = initial.advance_to(0.0, config(1.0), || Ok(())).unwrap();
    let conduction = ConductionObserver::new(
        &initial,
        &[(MaterialId(1), &law as &dyn ConductivityLaw)],
        &[ConductiveLink {
            link_id: 1,
            half_distances_m: [0.1, 0.1],
            contact_resistance_m2_k_w: 0.0,
        }],
        &[],
        || Ok(()),
    )
    .unwrap();
    let implicit = ImplicitConductionChain::new(
        &initial,
        &conduction,
        &[(MaterialId(1), &law as &dyn ImplicitConductionLaw)],
        || Ok(()),
    )
    .unwrap();
    // dt=500 s exceeds the cold cell's 80 s conservative explicit limit,
    // including the original timeline's 0.8 safety factor.
    // At the accepted endpoint both cells are homogeneous, so the independent
    // BE balance is linear in transferred heat Q, but includes the original
    // cold cell's miscibility enthalpy. The binodal below was computed at
    // 75-digit precision independently of the production equilibrium solver.
    let dt = 500.0;
    let dilute = 0.129_488_816_883_917_17;
    let mixing_heat = 5.0 * R * 900.0 * (0.4 * 0.6 - dilute * (1.0 - dilute));
    let conductance = 5.0;
    let q = dt * conductance * (240.0 + mixing_heat / 500.0)
        / (1.0 + dt * conductance * (1.0 / 500.0 + 1.0 / 1000.0));
    let expected = [350.0 + (q - mixing_heat) / 500.0, 590.0 - q / 1000.0];
    assert!(expected[0] > 450.0);
    let before = initial.cells().clone();
    let (next, report) = implicit
        .advance_to(
            &initial,
            dt,
            ImplicitConductionControls {
                maximum_dt_s: dt,
                cell_tolerance_j: 2e-6,
                total_tolerance_j: 1e-5,
                ..ImplicitConductionControls::default()
            },
            &mut ImplicitConductionWorkspace::new(2).unwrap(),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(report.energy.substeps, 1);
    assert_eq!(report.rejected_intervals, 0);
    close(next.temperature_k(CellId(1)).unwrap(), expected[0], 2e-7);
    close(next.temperature_k(CellId(2)).unwrap(), expected[1], 2e-7);
    close(
        0.336 * (next.cells()[&CellId(1)].specific_h_j_kg - h[0]),
        q,
        1e-5,
    );
    close(next.energy_j(), initial.energy_j(), 1e-5);
    assert_eq!(report.energy.external_heat_j, 0.0);
    assert!(report.energy.residual_j.abs() < 1e-5);
    assert!(report.entropy.thermal_storage_change_j_k > 0.0);
    assert!(report.entropy.conduction_generation_j_k > 0.0);
    assert_eq!(next.mass_kg(), initial.mass_kg());
    assert_eq!(initial.cells(), &before);
    assert_eq!(initial.time_s(), 0.0);
}
