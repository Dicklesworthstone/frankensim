//! Source-bound dwell recovery through actual thermal and free-sided consumers.

use super::*;
use fs_material::fiber::Uniaxial;
use fs_material::plastic::chaboche::recovery::STATIC_RECOVERY_RATE_DIMS;

const RATES: [&str; 2] = ["branch_fast_static_rate", "branch_slow_static_rate"];
fn rate_names() -> [String; 2] {
    RATES.map(str::to_owned)
}
fn recovery_properties(varying: bool) -> Vec<PropertyClaim> {
    let mut data = if varying {
        rate_properties()
    } else {
        properties([0.0; 2])
    };
    data.extend([
        scalar(RATES[0], STATIC_RECOVERY_RATE_DIMS, 1.0),
        scalar(RATES[1], STATIC_RECOVERY_RATE_DIMS, 0.2),
    ]);
    if varying {
        curve(
            data.iter_mut().find(|p| p.key.name() == RATES[0]).unwrap(),
            vec![(280.0, 0.1), (300.0, 1.0), (310.0, 4.0), (320.0, 6.0)],
        );
        curve(
            data.iter_mut().find(|p| p.key.name() == RATES[1]).unwrap(),
            vec![(280.0, 0.0), (300.0, 0.2), (320.0, 0.3)],
        );
    }
    data
}
fn selected(source: ThermalChabocheSource<'_, 2>) -> ThermalChabocheSource<'_, 2> {
    source
        .with_static_recovery(MaterialPropertySelection::SingleClaimOnly, rate_names())
        .unwrap()
}
fn expected_rates(t: f64) -> [f64; 2] {
    let fast = if t <= 300.0 {
        0.1 + 0.9 * (t - 280.0) / 20.0
    } else if t <= 310.0 {
        1.0 + 0.3 * (t - 300.0)
    } else {
        4.0 + 0.2 * (t - 310.0)
    };
    let slow = if t <= 300.0 {
        0.01 * (t - 280.0)
    } else {
        0.2 + 0.005 * (t - 300.0)
    };
    [fast, slow]
}

#[test]
fn sourced_hot_dwell_and_reversal_retain_all_branches_original_energy_and_axial_tangent() {
    let fixture = Fixture::new(recovery_properties(true));
    for source in [fixture.binding(), fixture.rate(false), fixture.rate(true)].map(selected) {
        let mut law = source.at_uniaxial_temperatures(300.0, 300.0, 0.01).unwrap();
        let mut state = law.initial_state();
        let mut old_t = 300.0;
        let mut old_strain = 0.0;
        for (t, strain, dt) in [
            (300.0, 0.012, 0.01),
            (315.0, 0.012, 0.15),
            (315.0, 0.012, 0.4),
            (290.0, -0.009, 0.12),
        ] {
            let interval = source.at_uniaxial_temperatures(old_t, t, dt).unwrap();
            law.admit_successor(&interval).unwrap();
            law = interval;
            let retained = state.clone();
            let old_tensor = law.tensor_strain(old_strain, &state).unwrap();
            let old_storage = law.stored_energy(old_strain, &state).unwrap().total_j_m3();
            let next = law
                .advance_with_energy(old_strain, &state, strain, tol())
                .unwrap();
            let borrowed = source
                .advance_with_energy_for_step(
                    h(old_t),
                    h(t),
                    &old_tensor,
                    state.mechanical.tensor_history(),
                    &next.strain,
                    tol(),
                    dt,
                )
                .unwrap();
            assert_eq!(state, retained);
            assert_eq!(next.state.mechanical.tensor_history(), &borrowed.state);
            assert_eq!(next.thermal.mechanical, borrowed.thermal.mechanical);
            assert_eq!(next.state.temperature_k(), t);
            let d = next.state.mechanical.equivalent_plastic_strain()
                - state.mechanical.equivalent_plastic_strain();
            let dp = next.state.mechanical.plastic_strain() - state.mechanical.plastic_strain();
            let x0 = state.mechanical.backstress_pa();
            let x1 = next.state.mechanical.backstress_pa();
            let rates = expected_rates(t);
            let prepared = law.tensor_law().static_recovery_law().unwrap();
            assert_eq!(prepared.duration_s(), dt);
            let mut physical_recovery = 0.0;
            let mut coefficient = 0.0;
            let mut storage = 0.0;
            for i in 0..2 {
                near(prepared.rates_per_s()[i], rates[i], 1e-14);
                let c0 = C[i] + DC[i] * (old_t - 300.0);
                let c1 = C[i] + DC[i] * (t - 300.0);
                near(
                    x1[i],
                    (x0[i] * c1 / c0 + c1 * dp) / (1.0 + GAMMA[i] * d + rates[i] * dt),
                    2e-6,
                );
                physical_recovery += (GAMMA[i] * d + rates[i] * dt) * x1[i] * x1[i] / c1;
                coefficient += 0.5 * (c1 - c0) * (x0[i] / c0).powi(2);
                storage += 0.5 * x1[i] * x1[i] / c1;
                let receipt = law
                    .material()
                    .property(RATES[i])
                    .unwrap()
                    .answer()
                    .receipt
                    .selected;
                let original = source
                    .base()
                    .reference()
                    .material()
                    .property(RATES[i])
                    .unwrap()
                    .answer()
                    .receipt
                    .selected;
                assert_eq!(receipt, original);
            }
            assert!(physical_recovery > 0.0);
            near(next.recovery_dissipation_j_m3, physical_recovery, 1e-7);
            near(
                next.thermal.mechanical.energy.physical_dissipation_j_m3,
                SY * d + physical_recovery + next.overstress_dissipation_j_m3,
                1e-7,
            );
            near(
                next.thermal.coefficient_storage_change_j_m3,
                coefficient,
                1e-7,
            );
            near(
                next.thermal
                    .mechanical
                    .energy
                    .final_state
                    .kinematic_hardening_j_m3,
                0.5 * HK * next.state.mechanical.plastic_strain().powi(2) + storage,
                1e-7,
            );
            near(
                next.thermal.stored_change_j_m3(),
                law.stored_energy(strain, &next.state).unwrap().total_j_m3() - old_storage,
                1e-6,
            );
            let above = law.response(strain + 1e-7, &state).unwrap().stress_pa;
            let below = law.response(strain - 1e-7, &state).unwrap().stress_pa;
            near(next.tangent_pa, (above - below) / 2e-7, 10.0);
            law.validate_history(strain, &next.state, tol().yield_pa)
                .unwrap();
            state = next.state;
            old_t = t;
            old_strain = strain;
        }
    }
}

#[test]
fn source_rates_require_complete_support_distinct_roles_actual_time_and_original_family() {
    let data = recovery_properties(true);
    let fixture = Fixture::new(data.clone());
    let source = selected(fixture.binding());
    assert_eq!(source.static_recovery_properties(), Some(&rate_names()));
    assert!(source.resolve_enthalpy(0.0, h(315.0)).is_err());
    assert!(source.at_uniaxial_temperatures(300.0, 315.0, 0.0).is_err());
    assert!(source.at_uniaxial_temperatures(300.0, 325.0, 0.1).is_err());
    assert!(
        source
            .clone()
            .with_static_recovery(MaterialPropertySelection::SingleClaimOnly, rate_names(),)
            .is_err()
    );
    for names in [
        [RATES[0].to_owned(), RATES[0].to_owned()],
        [schema()[0].recovery.clone(), RATES[1].to_owned()],
        ["missing_time_recovery".to_owned(), RATES[1].to_owned()],
    ] {
        assert!(
            fixture
                .binding()
                .with_static_recovery(MaterialPropertySelection::SingleClaimOnly, names,)
                .is_err()
        );
    }
    for case in 0..3 {
        let mut changed = data.clone();
        let property = changed
            .iter_mut()
            .find(|p| p.key.name() == RATES[0])
            .unwrap();
        match case {
            0 => *property = scalar(RATES[0], Pressure::DIMS, 1.0),
            1 => curve(
                property,
                vec![(280.0, 0.1), (300.0, 1.0), (310.0, -0.1), (320.0, 6.0)],
            ),
            _ => curve(property, vec![(295.0, 0.1), (300.0, 1.0), (320.0, 6.0)]),
        }
        let invalid = Fixture::new(changed);
        assert!(
            invalid
                .binding()
                .with_static_recovery(MaterialPropertySelection::SingleClaimOnly, rate_names(),)
                .is_err()
        );
    }
    let accepted = source.at_uniaxial_temperatures(300.0, 315.0, 0.2).unwrap();
    let without = fixture
        .binding()
        .at_uniaxial_temperatures(315.0, 315.0, 0.3)
        .unwrap();
    assert!(accepted.admit_successor(&without).is_err());
    let wrong_clock = source.at_uniaxial_temperatures(300.0, 315.0, 0.3).unwrap();
    assert!(accepted.admit_successor(&wrong_clock).is_err());
    let next = source.at_uniaxial_temperatures(315.0, 305.0, 0.3).unwrap();
    accepted.admit_successor(&next).unwrap();

    let mut altered = data;
    altered
        .iter_mut()
        .find(|p| p.key.name() == RATES[0])
        .unwrap()
        .provenance
        .source = "different independent static-recovery calibration".into();
    let foreign = Fixture::new(altered);
    let foreign = selected(foreign.binding())
        .at_uniaxial_temperatures(315.0, 305.0, 0.3)
        .unwrap();
    assert!(accepted.admit_successor(&foreign).is_err());
}

#[test]
fn uniform_adiabatic_dwell_heats_from_recovery_without_consuming_plastic_time_on_retry() {
    let fixture = Fixture::new(recovery_properties(false));
    let source = selected(fixture.binding());
    let point = CaloricPoint::try_new(fixture.query(), 0.0, 1e-6, 1e-3)
        .unwrap()
        .with_voce_hardening(MaterialPropertySelection::SingleClaimOnly)
        .unwrap()
        .with_chaboche_branches(MaterialPropertySelection::SingleClaimOnly, schema())
        .unwrap()
        .with_static_recovery(MaterialPropertySelection::SingleClaimOnly, rate_names())
        .unwrap();
    let reservoir = *point.caloric_state();
    let load_strain = [0.012, -0.006, -0.006, 0.0, 0.0, 0.0];
    let loaded = point
        .advance(
            load_strain,
            0.001,
            feedback_config(),
            |input| finite_reservoir_proposal(input, reservoir, 0.0),
            || false,
        )
        .unwrap()
        .next;
    let strain = loaded.history().plastic.plastic_strain;
    let unloaded = loaded
        .advance(
            strain,
            0.001,
            feedback_config(),
            |input| finite_reservoir_proposal(input, reservoir, 0.0),
            || false,
        )
        .unwrap()
        .next;
    assert_eq!(*unloaded.stress_pa(), [0.0; 6]);
    let retained = (
        unloaded.history().clone(),
        *unloaded.caloric_state(),
        unloaded.time_s(),
    );
    assert!(matches!(
        unloaded.advance(
            strain,
            0.3,
            feedback_config(),
            |input| finite_reservoir_proposal(input, reservoir, 0.0),
            || true,
        ),
        Err(FeedbackError::Cancelled)
    ));
    assert_eq!(
        retained,
        (
            unloaded.history().clone(),
            *unloaded.caloric_state(),
            unloaded.time_s()
        )
    );
    let next = unloaded
        .advance(
            strain,
            0.3,
            feedback_config(),
            |input| finite_reservoir_proposal(input, reservoir, 0.0),
            || false,
        )
        .unwrap();
    assert_eq!(next.next.history().plastic, unloaded.history().plastic);
    assert_eq!(*next.next.stress_pa(), [0.0; 6]);
    assert_eq!(next.mechanical_work_j, 0.0);
    assert_eq!(next.transport.external_heat_j, 0.0);
    assert!(next.physical_heat_j > 0.0);
    assert!(next.next.caloric_state().temperature_k() > unloaded.caloric_state().temperature_k());
    assert_eq!(next.reversible_heat_j, 0.0);
    let mut recovery = 0.0;
    for (i, r) in [1.0, 0.2].into_iter().enumerate() {
        let x0 = unloaded.history().backstress_pa[i];
        let x1 = next.next.history().backstress_pa[i];
        for j in 0..6 {
            near(x1[j], x0[j] / (1.0 + r * 0.3), 1e-8);
        }
        recovery += 1.5 * r * 0.3 * inner(&x1, &x1) / C[i];
    }
    near(next.physical_heat_j, unloaded.volume_m3() * recovery, 1e-12);
    let change = unloaded.mass_kg()
        * (next.next.caloric_state().specific_enthalpy_j_kg()
            - unloaded.caloric_state().specific_enthalpy_j_kg());
    near(change, next.physical_heat_j, 1e-8);
    near(
        change
            + unloaded.volume_m3()
                * (next.mechanical.energy.stored_change_j_m3()
                    + next.mechanical.energy.numerical_dissipation_j_m3),
        0.0,
        1e-8,
    );
    let native = source
        .advance_with_energy_for_step(
            unloaded.caloric_state().specific_enthalpy_j_kg(),
            next.next.caloric_state().specific_enthalpy_j_kg(),
            unloaded.strain(),
            unloaded.history(),
            &strain,
            tol(),
            0.3,
        )
        .unwrap();
    assert_eq!(*next.next.history(), native.state);
    assert_eq!(next.mechanical, native.thermal.mechanical);
    assert!(
        unloaded
            .clone()
            .with_static_recovery(MaterialPropertySelection::SingleClaimOnly, rate_names(),)
            .is_err()
    );
}
