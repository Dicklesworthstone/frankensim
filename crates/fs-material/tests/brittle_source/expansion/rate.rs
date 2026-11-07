//! Timed thermal deformation retains native mechanics, work and endpoint rates.
use super::super::{card, near, properties, scalar, settings};
use super::{query, resolve};
use fs_evidence::ValidityDomain;
use fs_matdb::{InterpolationPolicy, MaterialCard, PropertyClaim, PropertyValue};
use fs_material::brittle::source::expansion::ResolvedThermalBrittle;
use fs_material::brittle::source::expansion::rate::resolve_thermal_viscous_brittle;
use fs_material::brittle::source::rate::{
    BRITTLE_NORMAL_VISCOSITY_PROPERTY, BRITTLE_TANGENTIAL_VISCOSITY_PROPERTY,
    BRITTLE_VISCOSITY_DIMS,
};
use fs_material::state_point::{
    INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
    MaterialPropertySelection, YOUNG_MODULUS_PROPERTY,
};
use fs_material::tensor::contract;
use fs_material::{SmallStrainLaw, Voigt};
use fs_qty::{Pressure, QuantitySpec, Temperature};

fn claims() -> Vec<PropertyClaim> {
    let mut claims = properties();
    let mut alpha = scalar(
        LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
        INVERSE_TEMPERATURE_DIMS,
        1e-5,
    );
    alpha.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(300.0, 1e-5), (400.0, 3e-5)],
        dims: INVERSE_TEMPERATURE_DIMS,
    };
    alpha.validity = ValidityDomain::unconstrained().with_quantity(
        "T",
        QuantitySpec::dimensional(Temperature::DIMS),
        300.0,
        400.0,
    );
    alpha.interpolation = InterpolationPolicy::LinearInside;
    claims.extend([
        alpha,
        scalar(
            BRITTLE_NORMAL_VISCOSITY_PROPERTY,
            BRITTLE_VISCOSITY_DIMS,
            1e8,
        ),
        scalar(
            BRITTLE_TANGENTIAL_VISCOSITY_PROPERTY,
            BRITTLE_VISCOSITY_DIMS,
            2e8,
        ),
    ]);
    claims
}

fn viscous(source: &MaterialCard, reference: f64, current: f64) -> ResolvedThermalBrittle {
    resolve_thermal_viscous_brittle(
        source,
        &query(reference),
        &query(current),
        MaterialPropertySelection::SingleClaimOnly,
        settings(0.01),
        1.0,
    )
    .unwrap()
}

fn near_tensor(actual: &Voigt, expected: &Voigt) {
    for i in 0..6 {
        near(actual[i], expected[i]);
    }
}

#[test]
fn free_thermal_motion_has_no_dashpot_stress_while_restraint_drives_timed_tension() {
    let source = card(claims());
    let start = viscous(&source, 400.0, 400.0);
    let cold = start.at_temperature(&source, 300.0).unwrap();
    let virgin = start.band().law().initial_state();
    let saved = virgin;
    let prepared = start.transfer_coefficients(&virgin, &cold, 1e-6).unwrap();
    assert_eq!(prepared.parameter_work_j_m3, 0.0);
    assert_eq!(prepared.damage_dissipated_j_m3, 0.0);
    let change = cold.free_linear_strain() - start.free_linear_strain();
    near(change, -0.002);
    let free_total = [change, change, change, 0.0, 0.0, 0.0];
    for dt in [0.02, 1e-9] {
        // Even a fast thermal/free displacement trajectory has zero mechanical
        // rate. The adhesive jump-rate bound does not constrain matched expansion.
        let path = cold
            .trial_total_over_time_with_mean_stress(&free_total, &prepared.state, dt)
            .unwrap();
        assert_eq!(path.endpoint.stress_pa, [0.0; 6]);
        assert_eq!(path.mean_stress_pa, [0.0; 6]);
        assert_eq!(path.endpoint.equilibrium.state.strain(), &[0.0; 6]);
        assert_eq!(path.endpoint.equilibrium.state.damage(), 0.0);
        assert_eq!(path.endpoint.viscous_dissipated_j_m3, 0.0);
        assert_eq!(path.path_work_j_m3, 0.0);
        let free_rate = change / dt;
        let total_rate = [free_rate, free_rate, free_rate, 0.0, 0.0, 0.0];
        assert_eq!(
            cold.mechanical_strain_rate(&total_rate, free_rate).unwrap(),
            [0.0; 6]
        );
        assert_eq!(
            cold.stress_at_rate(&path.endpoint.equilibrium.state, &total_rate, free_rate)
                .unwrap(),
            [0.0; 6]
        );
    }

    let dt = 0.02;
    let restrained = cold
        .trial_total_over_time_with_mean_stress(&[0.0; 6], &prepared.state, dt)
        .unwrap();
    let endpoint = &restrained.endpoint;
    assert!(endpoint.stress_pa[0] > endpoint.equilibrium.stress_pa[0]);
    assert!(endpoint.equilibrium.state.damage() > 0.0);
    assert!(endpoint.equilibrium.energy.dissipated > 0.0);
    assert!(endpoint.viscous_dissipated_j_m3 > 0.0);
    let eigenstrain_work = -change * restrained.mean_stress_pa[..3].iter().sum::<f64>();
    near(eigenstrain_work, restrained.path_work_j_m3);
    near(
        eigenstrain_work,
        endpoint.equilibrium.energy.stored_change
            + endpoint.equilibrium.energy.dissipated
            + endpoint.viscous_dissipated_j_m3,
    );
    near_tensor(
        &cold.total_strain(&endpoint.equilibrium.state).unwrap(),
        &[0.0; 6],
    );
    near_tensor(
        &cold
            .stress_at_rate(&endpoint.equilibrium.state, &[0.0; 6], change / dt)
            .unwrap(),
        &endpoint.stress_pa,
    );
    assert!(
        cold.trial_total_over_time(&[0.0; 6], &prepared.state, 1e-9)
            .is_err()
    );
    // Holding total velocity at zero is a different physical observation from
    // following the free thermal velocity, even at zero mechanical strain.
    assert!(
        cold.stress_at_rate(&prepared.state, &[0.0; 6], change / dt)
            .unwrap()[0]
            > 0.0
    );
    assert!(cold.trial_total(&[0.0; 6], &prepared.state).is_err());
    assert_eq!(virgin, saved);
    assert_eq!(prepared.state.damage(), 0.0);
}

#[test]
fn total_strain_jacobians_and_actual_endpoint_rates_keep_the_native_thermal_chain_rule() {
    let source = card(claims());
    let start = viscous(&source, 300.0, 300.0);
    let next = start.at_temperature(&source, 301.0).unwrap();
    let dt = 0.02;
    let old_total = [0.00025, 0.0001, -0.00005, 0.00003, 0.00002, 0.00001];
    let total = [0.0004, 0.00016, -0.00002, 0.00006, 0.00004, -0.00003];
    let old = start
        .trial_total_over_time(&old_total, &start.band().law().initial_state(), dt)
        .unwrap()
        .equilibrium
        .state;
    let prepared = start.transfer_coefficients(&old, &next, 1e-6).unwrap();
    let saved = prepared.state;
    let endpoint_law = next.small_strain_law_over_time(dt).unwrap();
    let mean_law = next.mean_small_strain_law_over_time(dt).unwrap();
    assert_eq!(endpoint_law.duration_s(), dt);
    assert_eq!(mean_law.duration_s(), dt);
    let (path, mean_tangent) = mean_law.trial(&total, &prepared.state).unwrap();
    assert_eq!(path.endpoint.equilibrium.state.damage(), 0.0);
    let mean_tangent = mean_tangent.expect("strictly tensile intact path is smooth");
    let endpoint_tangent = path.endpoint.tangent_pa.expect("intact endpoint is smooth");
    let native = next
        .band()
        .law()
        .trial_over_time_with_mean_stress_tangent(
            &next.mechanical_strain(&total).unwrap(),
            &prepared.state,
            dt,
        )
        .unwrap();
    assert_eq!(path, native.0);
    assert_eq!(Some(mean_tangent), native.1);
    assert_eq!(
        endpoint_law.trial(&total, &prepared.state).unwrap(),
        path.endpoint
    );
    assert_eq!(
        endpoint_law.stress(&total, &prepared.state),
        path.endpoint.stress_pa
    );
    assert_eq!(
        endpoint_law.tangent(&total, &prepared.state),
        endpoint_tangent
    );
    assert_eq!(
        mean_law.stress(&total, &prepared.state),
        path.mean_stress_pa
    );
    assert_eq!(mean_law.tangent(&total, &prepared.state), mean_tangent);
    assert_eq!(
        endpoint_law.update_state(&total, &prepared.state),
        path.endpoint.equilibrium.state
    );
    assert_eq!(
        mean_law.update_state(&total, &prepared.state),
        path.endpoint.equilibrium.state
    );
    let stiffness = next.band().law().elastic().stiffness();
    for j in 0..6 {
        let h = 1e-8;
        let mut plus = total;
        let mut minus = total;
        plus[j] += h;
        minus[j] -= h;
        let plus = next
            .trial_total_over_time_with_mean_stress(&plus, &prepared.state, dt)
            .unwrap();
        let minus = next
            .trial_total_over_time_with_mean_stress(&minus, &prepared.state, dt)
            .unwrap();
        for i in 0..6 {
            // The mean elastic part is half the endpoint stiffness; the native
            // eta/dt part is NOT halved. Temperature is fixed in this derivative.
            assert!(
                (endpoint_tangent[i][j] - mean_tangent[i][j] - 0.5 * stiffness[i][j]).abs() < 1e-3
            );
            for (analytic, numerical) in [
                (
                    endpoint_tangent[i][j],
                    (plus.endpoint.stress_pa[i] - minus.endpoint.stress_pa[i]) / (2.0 * h),
                ),
                (
                    mean_tangent[i][j],
                    (plus.mean_stress_pa[i] - minus.mean_stress_pa[i]) / (2.0 * h),
                ),
            ] {
                assert!(
                    (analytic - numerical).abs()
                        <= 3e-6 * analytic.abs().max(numerical.abs()).max(3e3),
                    "{analytic:e} != {numerical:e}"
                );
            }
        }
    }

    let change: Voigt = core::array::from_fn(|i| total[i] - old_total[i]);
    let free_change = next.free_linear_strain() - start.free_linear_strain();
    let total_rate = change.map(|value| 2.0 * value / dt);
    let free_rate = free_change / dt;
    let mechanical_rate = next.mechanical_strain_rate(&total_rate, free_rate).unwrap();
    for i in 0..6 {
        near(
            mechanical_rate[i],
            total_rate[i] - if i < 3 { free_rate } else { 0.0 },
        );
    }
    let accepted = path.endpoint.equilibrium.state;
    let observed = next
        .stress_at_rate(&accepted, &total_rate, free_rate)
        .unwrap();
    assert_eq!(
        observed,
        next.band()
            .law()
            .stress_at_rate(&accepted, &mechanical_rate)
            .unwrap()
    );
    let missing_thermal_rate = next
        .band()
        .law()
        .stress_at_rate(&accepted, &total_rate)
        .unwrap();
    assert!((observed[0] - missing_thermal_rate[0]).abs() > 1.0);
    assert!((observed[0] - path.endpoint.stress_pa[0]).abs() > 1.0);
    near(
        contract(&path.mean_stress_pa, &change)
            - free_change * path.mean_stress_pa[..3].iter().sum::<f64>(),
        path.path_work_j_m3,
    );
    assert_eq!(prepared.state, saved);
}

#[test]
fn coefficient_displacement_eigenstrain_fracture_and_viscous_accounts_close_separately() {
    let mut claims = claims();
    let modulus = claims
        .iter_mut()
        .find(|p| p.key.name() == YOUNG_MODULUS_PROPERTY)
        .unwrap();
    modulus.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(300.0, 1e9), (400.0, 0.5e9)],
        dims: Pressure::DIMS,
    };
    modulus.validity = ValidityDomain::unconstrained().with_quantity(
        "T",
        QuantitySpec::dimensional(Temperature::DIMS),
        300.0,
        400.0,
    );
    modulus.interpolation = InterpolationPolicy::LinearInside;
    let source = card(claims);
    let start = viscous(&source, 300.0, 300.0);
    let next = start.at_temperature(&source, 400.0).unwrap();
    let old_total = [0.0003, 0.0002, -0.0001, 0.0001, 0.0002, -0.00015];
    let old = start
        .trial_total_over_time(&old_total, &start.band().law().initial_state(), 0.02)
        .unwrap()
        .equilibrium
        .state;
    assert_eq!(old.damage(), 0.0);
    let saved = old;
    let coefficient = start.transfer_coefficients(&old, &next, 1e-6).unwrap();
    near(
        coefficient.parameter_work_j_m3,
        -0.5 * old.stored_energy_j_m3(),
    );
    assert_eq!(coefficient.state.strain(), old.strain());
    let total = [0.006, 0.0001, -0.0001, 0.001, 0.0005, 0.0];
    let path = next
        .trial_total_over_time_with_mean_stress(&total, &coefficient.state, 0.04)
        .unwrap();
    let endpoint = &path.endpoint;
    assert!(endpoint.equilibrium.state.damage() > 0.0);
    assert!(endpoint.equilibrium.energy.dissipated > 0.0);
    assert!(endpoint.viscous_dissipated_j_m3 > 0.0);
    let change: Voigt = core::array::from_fn(|i| total[i] - old_total[i]);
    let displacement_work = contract(&path.mean_stress_pa, &change);
    let eigenstrain_work = -(next.free_linear_strain() - start.free_linear_strain())
        * path.mean_stress_pa[..3].iter().sum::<f64>();
    near(displacement_work + eigenstrain_work, path.path_work_j_m3);
    let fracture = coefficient.damage_dissipated_j_m3 + endpoint.equilibrium.energy.dissipated;
    let stored_change = endpoint.equilibrium.state.stored_energy_j_m3() - old.stored_energy_j_m3();
    near(
        coefficient.parameter_work_j_m3 + displacement_work + eigenstrain_work,
        stored_change + fracture + endpoint.viscous_dissipated_j_m3,
    );
    near(
        endpoint.equilibrium.state.dissipated_j_m3() - old.dissipated_j_m3(),
        fracture,
    );
    near(
        endpoint.total_work_j_m3(),
        endpoint.equilibrium.energy.work + endpoint.viscous_dissipated_j_m3,
    );
    assert_eq!(old, saved);
}

#[test]
fn original_no_viscosity_response_and_timed_owner_domain_refusals_are_retained() {
    let source = card(claims());
    let original = resolve(&source, 300.0, 300.0);
    let next = original.at_temperature(&source, 301.0).unwrap();
    assert!(next.band().law().viscosity().is_none());
    let state = original
        .transfer_coefficients(&original.band().law().initial_state(), &next, 1e-6)
        .unwrap()
        .state;
    let total = [0.004, 0.0, 0.0, 0.001, 0.0, 0.0];
    let untimed = next.mean_small_strain_law().trial(&total, &state).unwrap();
    let timed = next
        .trial_total_over_time_with_mean_stress_tangent(&total, &state, 0.02)
        .unwrap();
    assert_eq!(timed.0.endpoint.equilibrium, untimed.0.endpoint);
    assert_eq!(timed.0.endpoint.stress_pa, untimed.0.endpoint.stress_pa);
    assert_eq!(timed.0.mean_stress_pa, untimed.0.mean_stress_pa);
    assert_eq!(timed.0.path_work_j_m3, untimed.0.path_work_j_m3);
    assert_eq!(timed.1, untimed.1);
    assert_eq!(timed.0.endpoint.viscous_dissipated_j_m3, 0.0);
    // An inactive dashpot must not reject finite rates merely because their
    // otherwise irrelevant mechanical-rate difference is unrepresentable.
    let huge_rate = [f64::MAX; 6];
    assert!(next.mechanical_strain_rate(&huge_rate, -f64::MAX).is_err());
    assert_eq!(
        next.stress_at_rate(&timed.0.endpoint.equilibrium.state, &huge_rate, -f64::MAX)
            .unwrap(),
        untimed.0.endpoint.stress_pa
    );
    assert!(
        next.stress_at_rate(
            &timed.0.endpoint.equilibrium.state,
            &huge_rate,
            f64::INFINITY
        )
        .is_err()
    );

    let active = viscous(&source, 300.0, 300.0);
    let virgin = active.band().law().initial_state();
    for dt in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(active.trial_total_over_time(&total, &virgin, dt).is_err());
        assert!(
            active
                .trial_total_over_time_with_mean_stress(&total, &virgin, dt)
                .is_err()
        );
        assert!(active.small_strain_law_over_time(dt).is_err());
        assert!(active.mean_small_strain_law_over_time(dt).is_err());
    }
    assert!(
        active
            .mechanical_strain_rate(&[0.0; 6], f64::INFINITY)
            .is_err()
    );
    assert!(active.mechanical_strain_rate(&[f64::NAN; 6], 0.0).is_err());
    assert!(active.stress_at_rate(&virgin, &[0.0; 6], 200.0).is_err());
    assert!(
        active
            .trial_total_over_time(&[0.051, 0.0, 0.0, 0.0, 0.0, 0.0], &virgin, 0.02)
            .is_err()
    );
    assert!(active.trial_total(&total, &virgin).is_err());
    let identical = active.clone();
    let foreign = identical.band().law().initial_state();
    assert!(
        active
            .trial_total_over_time(&total, &foreign, 0.02)
            .is_err()
    );
    assert!(active.stress_at_rate(&foreign, &[0.0; 6], 0.0).is_err());
    let endpoint = active.small_strain_law_over_time(0.02).unwrap();
    assert!(
        endpoint
            .stress(&total, &foreign)
            .iter()
            .all(|value| value.is_nan())
    );
    assert!(!endpoint.update_state(&total, &foreign).is_valid());
    let mean = active.mean_small_strain_law_over_time(0.02).unwrap();
    let (held, tangent) = mean.trial(&[0.0; 6], &virgin).unwrap();
    assert!(tangent.is_none());
    assert!(
        mean.tangent(&[0.0; 6], &virgin)
            .iter()
            .flatten()
            .all(|value| value.is_nan())
    );
    assert_eq!(held.path_work_j_m3, 0.0);
    let failed = active
        .trial_total_over_time(&[0.03, 0.0, 0.0, 0.0, 0.0, 0.0], &virgin, 0.02)
        .unwrap();
    assert!(failed.equilibrium.state.failed());
    assert_eq!(
        active
            .stress_at_rate(&failed.equilibrium.state, &huge_rate, -f64::MAX)
            .unwrap(),
        failed.equilibrium.stress_pa
    );
    assert!(
        active
            .stress_at_rate(&failed.equilibrium.state, &huge_rate, f64::NAN)
            .is_err()
    );
    assert_eq!(virgin, active.band().law().initial_state());
}
