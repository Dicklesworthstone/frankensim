//! G0/G1: owned thermal intervals retain the original source and native work.

use super::*;
use fs_material::fiber::Uniaxial;
use fs_material::plastic::thermal::tabulated::axial::UniaxialThermalTabulatedState;

fn axial_fields(
    temperature: f64,
    total: f64,
    state: &UniaxialThermalTabulatedState,
    expansion: f64,
) -> (f64, f64, f64) {
    let e = E + DE * (temperature - 300.0);
    let hk = HK + DHK * (temperature - 300.0);
    let p = state.mechanical.plastic_strain();
    let alpha = state.mechanical.equivalent_plastic_strain();
    let stress = e * (total - expansion * (temperature - 300.0) - p);
    let elastic = stress * stress / (2.0 * e);
    let storage = elastic + 0.5 * hk * p * p + primitive(alpha, ordinates(temperature));
    let entropy =
        expansion * stress - DE / e * elastic - 0.5 * DHK * p * p - primitive(alpha, SLOPES);
    (storage, entropy, stress)
}

#[test]
fn heated_constant_force_crosses_measured_knot_with_exact_tensor_work_and_retained_extension() {
    let expansion = 2e-5;
    let fixture = Fixture::new(properties(expansion, true));
    let source = fixture.source();
    let cold = source.at_uniaxial_temperatures(300.0, 300.0, 0.1).unwrap();
    let initial = cold.initial_state();
    let alpha = 0.0059;
    let force_stress = 4e6 + HK * alpha;
    let cold_strain = alpha + force_stress / E;
    let loaded = cold
        .advance_with_energy(0.0, &initial, cold_strain, tolerance())
        .unwrap();
    near(loaded.stress_pa, force_stress, 1e-6);
    near(
        loaded.state.mechanical.equivalent_plastic_strain(),
        alpha,
        2e-14,
    );
    assert_eq!(loaded.state.temperature_k(), 300.0);

    let hot = source.at_uniaxial_temperatures(300.0, 320.0, 0.2).unwrap();
    cold.admit_successor(&hot).unwrap();
    assert!(std::ptr::eq(
        cold.hardening_source(),
        hot.hardening_source()
    ));
    let e_hot = E + 20.0 * DE;
    let hk_hot = HK + 20.0 * DHK;
    let slope_hot = 4e8;
    let alpha_hot = (force_stress - 3.6e6 + slope_hot * 0.006) / (hk_hot + slope_hot);
    let hot_strain = expansion * 20.0 + alpha_hot + force_stress / e_hot;
    let heated = hot
        .advance_with_energy(cold_strain, &loaded.state, hot_strain, tolerance())
        .unwrap();
    near(heated.stress_pa, force_stress, 1e-6);
    near(
        heated.state.mechanical.equivalent_plastic_strain(),
        alpha_hot,
        2e-14,
    );
    assert!(heated.state.mechanical.equivalent_plastic_strain() > 0.006);
    assert_eq!(heated.state.temperature_k(), 320.0);
    near(
        heated.tangent_pa,
        e_hot * (hk_hot + slope_hot) / (e_hot + hk_hot + slope_hot),
        1e-6,
    );
    let perturbation = 1e-8;
    let plus = hot
        .response(hot_strain + perturbation, &loaded.state)
        .unwrap();
    let minus = hot
        .response(hot_strain - perturbation, &loaded.state)
        .unwrap();
    near(
        (plus.stress_pa - minus.stress_pa) / (2.0 * perturbation),
        heated.tangent_pa,
        0.5,
    );
    let native = hot
        .tensor_law()
        .response(&heated.strain, &loaded.state.mechanical.tensor_history())
        .unwrap();
    near(native.stress_pa[0], heated.stress_pa, 1e-6);
    for stress in &native.stress_pa[1..] {
        near(*stress, 0.0, 1e-6);
    }
    near(native.state.alpha, alpha_hot, 2e-14);

    let old_total = cold.tensor_strain(cold_strain, &loaded.state).unwrap();
    let borrowed = source
        .advance_with_energy_for_step(
            0.0,
            h(320.0),
            &old_total,
            &loaded.state.mechanical.tensor_history(),
            &heated.strain,
            tolerance(),
            0.2,
        )
        .unwrap();
    assert_eq!(heated.thermal.mechanical.state, borrowed.mechanical.state);
    near(
        heated.thermal.reversible_heat_j_m3,
        borrowed.reversible_heat_j_m3,
        0.0,
    );
    near(
        heated.thermal.thermoelastic_storage_change_j_m3,
        borrowed.thermoelastic_storage_change_j_m3,
        0.0,
    );
    let (psi0, entropy0, _) = axial_fields(300.0, cold_strain, &loaded.state, expansion);
    let (psi1, entropy1, _) = axial_fields(320.0, hot_strain, &heated.state, expansion);
    near(
        hot.stored_energy(cold_strain, &loaded.state)
            .unwrap()
            .total_j_m3(),
        psi0,
        1e-7,
    );
    near(
        hot.stored_energy(hot_strain, &heated.state)
            .unwrap()
            .total_j_m3(),
        psi1,
        1e-7,
    );
    near(heated.thermal.stored_change_j_m3(), psi1 - psi0, 1e-7);
    near(
        heated.thermal.thermoelastic_storage_change_j_m3,
        320.0 * entropy1 - 300.0 * entropy0,
        1e-7,
    );
    let work = force_stress * (hot_strain - cold_strain);
    near(heated.thermal.mechanical_work_j_m3, work, 1e-7);
    near(
        heated.thermal.mechanical.energy.physical_dissipation_j_m3,
        ordinates(320.0)[0] * (alpha_hot - alpha),
        1e-7,
    );
    near(
        psi1 - psi0
            + heated.thermal.thermoelastic_storage_change_j_m3
            + heated.thermal.caloric_heat_j_m3()
            + heated.thermal.mechanical.energy.numerical_dissipation_j_m3,
        work,
        1e-5,
    );

    let cooling = source.at_uniaxial_temperatures(320.0, 300.0, 0.35).unwrap();
    hot.admit_successor(&cooling).unwrap();
    let residual = heated.state.mechanical.plastic_strain();
    let released = cooling
        .advance_with_energy(hot_strain, &heated.state, residual, tolerance())
        .unwrap();
    near(released.stress_pa, 0.0, 1e-6);
    near(released.state.mechanical.plastic_strain(), residual, 2e-14);
    near(
        released.state.mechanical.equivalent_plastic_strain(),
        alpha_hot,
        2e-14,
    );
    near(
        released.thermal.mechanical.energy.physical_dissipation_j_m3,
        0.0,
        1e-7,
    );
    assert!(residual > alpha);
    assert_eq!(initial.mechanical.equivalent_plastic_strain(), 0.0);
}

#[test]
fn linear_and_power_thermal_holds_use_actual_duration_and_keep_viscous_heat_separate() {
    let expansion = 2e-5;
    let fixture = Fixture::new(properties(expansion, false));
    let strain = 0.012;
    for power in [false, true] {
        let base = if power {
            fixture
                .base()
                .with_power_law_overstress(MaterialPropertySelection::SingleClaimOnly)
        } else {
            fixture
                .base()
                .with_linear_overstress(MaterialPropertySelection::SingleClaimOnly)
        }
        .unwrap();
        let source = base.with_tabulated_hardening(SINGLE, SAMPLE_LIMIT).unwrap();
        let cold = source.at_uniaxial_temperatures(300.0, 300.0, 0.02).unwrap();
        let initial = cold.initial_state();
        let loaded = cold
            .advance_with_energy(0.0, &initial, strain, tolerance())
            .unwrap();
        let old = loaded.state.clone();
        let old_alpha = old.mechanical.equivalent_plastic_strain();
        assert!(old_alpha > 0.002 && old_alpha < 0.006);
        let mut flow = Vec::new();
        for duration in [0.02, 0.4] {
            let hot = source
                .at_uniaxial_temperatures(300.0, 315.0, duration)
                .unwrap();
            cold.admit_successor(&hot).unwrap();
            assert_eq!(hot.duration_s(), duration);
            let held = hot
                .advance_with_energy(strain, &old, strain, tolerance())
                .unwrap();
            let rate = hot.tensor_law().rate_law().unwrap();
            let viscosity = 9.25e6;
            near(rate.viscosity_pa_s(), viscosity, 1e-7);
            // The 315 K segment of the table between its 0.006 and 0.012 knots.
            let hot_yield = ordinates(315.0);
            let slope = (hot_yield[3] - hot_yield[2]) / (KNOTS[3].0 - KNOTS[2].0);
            let stiffness = E + HK + slope;
            let excess = E * (strain - 15.0 * expansion)
                - stiffness * old_alpha
                - (hot_yield[2] - slope * KNOTS[2].0);
            // The endpoint lies beyond .006, so its original affine segment
            // gives an independent linear or quadratic equation in sqrt(delta).
            let delta = if power {
                let k = (2e6 * viscosity / duration).sqrt();
                let root = 2.0 * excess / (k + (k * k + 4.0 * stiffness * excess).sqrt());
                root * root
            } else {
                excess / (stiffness + viscosity / duration)
            };
            let expected = old_alpha + delta;
            assert!(expected > 0.006 && expected < 0.012);
            near(
                held.state.mechanical.equivalent_plastic_strain(),
                expected,
                3e-13,
            );
            assert_eq!(held.state.temperature_k(), 315.0);
            assert!(held.stress_pa < loaded.stress_pa);
            let overstress = if power {
                2e6 * (viscosity * delta / (2e6 * duration)).sqrt()
            } else {
                viscosity * delta / duration
            };
            near(
                held.thermal.mechanical.energy.physical_dissipation_j_m3,
                (ordinates(315.0)[0] + overstress) * delta,
                1e-5,
            );
            assert!(held.thermal.mechanical.energy.numerical_dissipation_j_m3 > 0.0);
            near(held.thermal.mechanical_work_j_m3, 0.0, 0.0);
            hot.validate_history(strain, &held.state, tolerance().yield_pa)
                .unwrap();
            let again = source
                .at_uniaxial_temperatures(315.0, 315.0, duration)
                .unwrap();
            hot.admit_successor(&again).unwrap();
            let relaxed = again
                .advance_with_energy(strain, &held.state, strain, tolerance())
                .unwrap();
            assert!(relaxed.state.mechanical.equivalent_plastic_strain() > expected);
            assert!(relaxed.stress_pa < held.stress_pa);
            flow.push(expected);
        }
        assert!(flow[1] > flow[0]);
        assert_eq!(old, loaded.state);
        assert_eq!(initial.mechanical.equivalent_plastic_strain(), 0.0);
    }
}

#[test]
fn original_source_and_caloric_identity_are_required_before_interval_rebinding() {
    let alternate_cp = "alternate-specific-heat-capacity";
    let mut claims = properties(0.0, false);
    claims.push(property(alternate_cp, CP_DIMS, CP, 0.0));
    let fixture = Fixture::new(claims);
    let source = fixture.source();
    let cold = source.at_uniaxial_temperatures(300.0, 300.0, 0.1).unwrap();
    assert_eq!(cold.density_kg_m3(), RHO);
    assert_eq!(
        cold.caloric_curve_identity(),
        fixture.caloric.curve().identity()
    );
    assert_eq!(cold.material().card_identity(), fixture.card.content_hash());
    let cloned = cold.clone();
    assert!(std::ptr::eq(cold.tensor_law(), cloned.tensor_law()));
    let hot = source.at_uniaxial_temperatures(300.0, 320.0, 0.5).unwrap();
    cold.admit_successor(&hot).unwrap();
    assert!(
        cold.admit_successor(&source.at_uniaxial_temperatures(301.0, 320.0, 0.5).unwrap())
            .is_err()
    );
    assert!(
        source
            .at_uniaxial_temperatures(300.0, 320.001, 0.1)
            .is_err()
    );
    assert!(source.at_uniaxial_temperatures(300.0, 310.0, 0.0).is_err());

    // The original card and numerical Cp/rho chart are the same. Only the
    // selected caloric claim differs, which must still refuse history transfer.
    let alternative = ResolvedHeatCapacityCurve::try_new(
        &fixture.card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new(alternate_cp, CP_DIMS),
            density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            temperature_axis: TEMPERATURE,
            lower: &point(280.0),
            upper: &point(320.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        -20.0 * CP,
        1e-6,
    )
    .unwrap();
    assert_eq!(alternative.curve().knots(), fixture.caloric.curve().knots());
    assert_ne!(
        alternative.curve().identity(),
        fixture.caloric.curve().identity()
    );
    let mut foreign_query = fixture.query();
    foreign_query.caloric = &alternative;
    let foreign = ThermoplasticSource::try_new(foreign_query, 0.0)
        .unwrap()
        .with_tabulated_hardening(SINGLE, SAMPLE_LIMIT)
        .unwrap()
        .at_uniaxial_temperatures(300.0, 320.0, 0.5)
        .unwrap();
    assert!(cold.admit_successor(&foreign).is_err());
    let changed_kinetics = fixture
        .base()
        .with_linear_overstress(MaterialPropertySelection::SingleClaimOnly)
        .unwrap()
        .with_tabulated_hardening(SINGLE, SAMPLE_LIMIT)
        .unwrap()
        .at_uniaxial_temperatures(300.0, 320.0, 0.5)
        .unwrap();
    assert!(cold.admit_successor(&changed_kinetics).is_err());

    let initial = cold.initial_state();
    let loaded = cold
        .advance_with_energy(0.0, &initial, 0.01, tolerance())
        .unwrap();
    let old = loaded.state.clone();
    assert!(
        hot.advance_with_energy(0.01, &old, 0.1, tolerance())
            .is_err()
    );
    assert_eq!(old, loaded.state);
    let accepted = hot
        .advance_with_energy(0.01, &old, 0.011, tolerance())
        .unwrap();
    assert!(hot.response(0.011, &accepted.state).is_err());
    assert!(cold.tensor_strain(0.011, &accepted.state).is_err());
    hot.validate_history(0.011, &accepted.state, tolerance().yield_pa)
        .unwrap();
    assert_eq!(old, loaded.state);
}

#[test]
fn stress_free_thermal_motion_reconstructs_lateral_expansion_and_obeys_full_tensor_domain() {
    let expansion = 2e-5;
    let fixture = Fixture::new(properties(expansion, true));
    let source = fixture.source();
    let law = source.at_uniaxial_temperatures(300.0, 320.0, 0.25).unwrap();
    let initial = law.initial_state();
    let free = 20.0 * expansion;
    let next = law
        .advance_with_energy(0.0, &initial, free, tolerance())
        .unwrap();
    for value in &next.strain[..3] {
        near(*value, free, 0.0);
    }
    for value in &next.strain[3..] {
        near(*value, 0.0, 0.0);
    }
    near(next.stress_pa, 0.0, 0.0);
    near(next.state.mechanical.equivalent_plastic_strain(), 0.0, 0.0);
    near(next.thermal.caloric_heat_j_m3(), 0.0, 0.0);
    near(
        law.stored_energy(0.0, &initial).unwrap().total_j_m3(),
        0.0,
        0.0,
    );
    near(
        law.stored_energy(free, &next.state).unwrap().total_j_m3(),
        0.0,
        0.0,
    );
    assert_eq!(law.tensor_strain(free, &next.state).unwrap(), next.strain);

    // A permissive scalar axial value cannot hide unsupported lateral thermal
    // strain. Low E keeps the original thermoelastic capacity bound admissible.
    let mut large_expansion = properties(0.00125, false);
    large_expansion[0] = property(YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 1e7, 0.0);
    let bounded_fixture = Fixture::new(large_expansion);
    let bounded = bounded_fixture
        .source()
        .at_uniaxial_temperatures(300.0, 320.0, 0.25)
        .unwrap();
    let virgin = bounded.initial_state();
    let axial = 0.025;
    assert!(
        axial
            < bounded
                .tensor_law()
                .tabulated_law()
                .backbone()
                .elastic
                .strain_limit
    );
    assert!(bounded.response(axial, &virgin).is_err());
    assert_eq!(virgin.mechanical.equivalent_plastic_strain(), 0.0);
}
