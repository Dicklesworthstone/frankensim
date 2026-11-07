//! G0/G1/G3 checks of original thermal cyclic histories on free-sided fibers.
//! Manufactured source data only; these are not measured alloy calibration.

use super::*;
use fs_material::fiber::Uniaxial;

fn branch_entropy(temperature: f64, axial_backstress: [f64; 2]) -> f64 {
    (0..2)
        .map(|i| {
            let modulus = C[i] + DC[i] * (temperature - 300.0);
            -0.5 * DC[i] * (axial_backstress[i] / modulus).powi(2)
        })
        .sum()
}

#[test]
fn owned_cyclic_intervals_retain_every_branch_and_native_energy_through_hot_reversal() {
    let fixture = Fixture::new(rate_properties());
    for source in [fixture.binding(), fixture.rate(false), fixture.rate(true)] {
        let mut interval = source.at_uniaxial_temperatures(300.0, 300.0, 4e-4).unwrap();
        let mut state = interval.initial_state();
        let mut old_temperature = 300.0;
        let mut old_strain = 0.0;
        for (temperature, strain, duration) in [
            (300.0, 0.012, 4e-4),
            (315.0, 0.012, 8e-4),
            (305.0, -0.009, 6e-4),
            (290.0, -0.007, 9e-4),
        ] {
            let next_interval = source
                .at_uniaxial_temperatures(old_temperature, temperature, duration)
                .unwrap();
            interval.admit_successor(&next_interval).unwrap();
            interval = next_interval;
            let saved = state.clone();
            let old_tensor = interval.tensor_strain(old_strain, &state).unwrap();
            let old_storage = interval.stored_energy(old_strain, &state).unwrap().total_j_m3();
            let candidate = interval
                .advance_with_energy(old_strain, &state, strain, tol())
                .unwrap();
            let borrowed = source
                .advance_with_energy_for_step(
                    h(old_temperature),
                    h(temperature),
                    &old_tensor,
                    state.mechanical.tensor_history(),
                    &candidate.strain,
                    tol(),
                    duration,
                )
                .unwrap();
            assert_eq!(state, saved, "proposing an interval must not mutate history");
            assert_eq!(candidate.state.mechanical.tensor_history(), &borrowed.state);
            assert_eq!(candidate.thermal.mechanical, borrowed.thermal.mechanical);
            near(
                candidate.thermal.caloric_heat_j_m3(),
                borrowed.thermal.caloric_heat_j_m3(),
                tol().energy_j_m3,
            );
            assert_eq!(candidate.state.temperature_k(), temperature);
            assert_ne!(
                candidate.state.mechanical.tensor_history().backstress_pa[0],
                candidate.state.mechanical.tensor_history().backstress_pa[1],
            );
            let old_a = state.mechanical.backstress_pa();
            let new_a = candidate.state.mechanical.backstress_pa();
            let alpha = candidate.state.mechanical.equivalent_plastic_strain();
            let delta = alpha - state.mechanical.equivalent_plastic_strain();
            let mut coefficient = 0.0;
            let mut branch_storage = 0.0;
            let mut recovery = 0.0;
            for i in 0..2 {
                let c0 = C[i] + DC[i] * (old_temperature - 300.0);
                let c1 = C[i] + DC[i] * (temperature - 300.0);
                let dp = candidate.state.mechanical.plastic_strain()
                    - state.mechanical.plastic_strain();
                near(
                    new_a[i],
                    ((old_a[i] / c0) * c1 + c1 * dp) / (1.0 + GAMMA[i] * delta),
                    2e-6,
                );
                coefficient += 0.5 * (c1 - c0) * (old_a[i] / c0).powi(2);
                branch_storage += 0.5 * new_a[i] * new_a[i] / c1;
                recovery += GAMMA[i] * delta * new_a[i] * new_a[i] / c1;
            }
            near(candidate.recovery_dissipation_j_m3, recovery, 1e-7);
            near(
                candidate.thermal.mechanical.energy.final_state.kinematic_hardening_j_m3,
                0.5 * HK * candidate.state.mechanical.plastic_strain().powi(2) + branch_storage,
                1e-7,
            );
            near(candidate.thermal.coefficient_storage_change_j_m3, coefficient, 1e-7);
            near(
                candidate.thermal.thermoelastic_storage_change_j_m3,
                temperature * branch_entropy(temperature, new_a)
                    - old_temperature * branch_entropy(old_temperature, old_a),
                1e-6,
            );
            near(
                candidate.thermal.mechanical.energy.physical_dissipation_j_m3,
                SY * delta
                    + candidate.recovery_dissipation_j_m3
                    + candidate.overstress_dissipation_j_m3,
                1e-7,
            );
            let new_storage = interval.stored_energy(strain, &candidate.state).unwrap().total_j_m3();
            near(
                candidate.thermal.stored_change_j_m3(),
                new_storage - old_storage,
                1e-6,
            );
            interval
                .validate_history(strain, &candidate.state, tol().yield_pa)
                .unwrap();
            state = candidate.state;
            old_strain = strain;
            old_temperature = temperature;
        }
    }
}

#[test]
fn axial_native_perzyna_matches_independent_flow_tangent_and_physical_heat() {
    let fixture = Fixture::new(rate_properties());
    for power in [false, true] {
        let source = fixture.rate(power);
        for duration in [0.002, 0.013] {
            let interval = source
                .at_uniaxial_temperatures(300.0, 300.0, duration)
                .unwrap();
            for alpha in [1e-7, 0.001, 0.004] {
                let viscosity = 2e6;
                let exponent = if power { 3.0 } else { 1.0 };
                let viscous = if power {
                    1e6 * (viscosity * alpha / (duration * 1e6)).powf(1.0 / exponent)
                } else {
                    viscosity * alpha / duration
                };
                let saturation = Q0 * -(-B0 * alpha).exp_m1();
                let branches: f64 = (0..2)
                    .map(|i| C[i] * alpha / (1.0 + GAMMA[i] * alpha))
                    .sum();
                let stress = SY + H * alpha + saturation + HK * alpha + branches + viscous;
                let strain = alpha + stress / E;
                let old = interval.initial_state();
                let response = interval.response(strain, &old).unwrap();
                near(response.state.mechanical.equivalent_plastic_strain(), alpha, 1e-14);
                near(response.state.mechanical.plastic_strain(), alpha, 1e-14);
                near(response.stress_pa, stress, 1e-4);
                let hardening = H + HK + Q0 * B0 * (-B0 * alpha).exp()
                    + (0..2)
                        .map(|i| C[i] / (1.0 + GAMMA[i] * alpha).powi(2))
                        .sum::<f64>()
                    + viscous / (exponent * alpha);
                near(response.tangent_pa, E * (hardening / (E + hardening)), 1e-3);
                near(
                    response.strain[1],
                    -NU * (strain - alpha) - 0.5 * alpha,
                    1e-14,
                );
                let checked = interval.advance_with_energy(0.0, &old, strain, tol()).unwrap();
                near(checked.overstress_dissipation_j_m3, viscous * alpha, 1e-6);
                let retained = checked.state.clone();
                interval.validate_history(strain, &retained, tol().yield_pa).unwrap();
                assert_eq!(retained, checked.state, "inspection must not consume another dt");
                assert!(interval.response(strain, &retained).unwrap().state.mechanical
                    .equivalent_plastic_strain() > retained.mechanical.equivalent_plastic_strain());
            }
        }
    }

    // Reverse yielding is shifted by each retained branch and Prager memory.
    let source = fixture.binding();
    let interval = source.at_uniaxial_temperatures(300.0, 300.0, 0.01).unwrap();
    let loaded = interval.response(0.012, &interval.initial_state()).unwrap();
    let p = loaded.state.mechanical.plastic_strain();
    let alpha = loaded.state.mechanical.equivalent_plastic_strain();
    let radius = SY + H * alpha + Q0 * -(-B0 * alpha).exp_m1();
    let translation = HK * p + loaded.state.mechanical.backstress_pa().iter().sum::<f64>();
    let reverse_stress = translation - radius;
    assert!(translation > 0.0 && reverse_stress > -radius);
    let reverse_strain = p + reverse_stress / E;
    let elastic = interval.response(reverse_strain + 1e-6, &loaded.state).unwrap();
    near(elastic.state.mechanical.equivalent_plastic_strain(), alpha, 1e-15);
    let reversed = interval.response(reverse_strain - 1e-6, &loaded.state).unwrap();
    assert!(reversed.state.mechanical.equivalent_plastic_strain() > alpha);
    assert!(reversed.state.mechanical.plastic_strain() < p);
}

#[test]
fn heating_unloaded_fiber_transports_branch_coordinates_without_creating_flow_or_heat() {
    let mut data = properties(DC);
    data.push(scalar(
        LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
        INVERSE_TEMPERATURE_DIMS,
        1e-5,
    ));
    let fixture = Fixture::new(data);
    let source = fixture.binding();
    let cold = source.at_uniaxial_temperatures(300.0, 300.0, 0.1).unwrap();
    let loaded = cold.advance_with_energy(0.0, &cold.initial_state(), 0.012, tol()).unwrap();
    let p = loaded.state.mechanical.plastic_strain();
    let unloaded = cold.advance_with_energy(0.012, &loaded.state, p, tol()).unwrap();
    near(unloaded.stress_pa, 0.0, tol().yield_pa);
    let hot = source.at_uniaxial_temperatures(300.0, 310.0, 0.2).unwrap();
    cold.admit_successor(&hot).unwrap();
    let free = 1e-4;
    let heated = hot.advance_with_energy(p, &unloaded.state, p + free, tol()).unwrap();
    near(heated.stress_pa, 0.0, tol().yield_pa);
    near(
        heated.state.mechanical.equivalent_plastic_strain(),
        unloaded.state.mechanical.equivalent_plastic_strain(),
        1e-15,
    );
    near(heated.state.mechanical.plastic_strain(), p, 1e-15);
    let old_a = unloaded.state.mechanical.backstress_pa();
    let new_a = heated.state.mechanical.backstress_pa();
    for i in 0..2 {
        near(new_a[i] / (C[i] + 10.0 * DC[i]), old_a[i] / C[i], 1e-15);
    }
    near(heated.strain[1], -0.5 * p + free, 1e-14);
    near(heated.thermal.caloric_heat_j_m3(), 0.0, tol().energy_j_m3);
    near(heated.recovery_dissipation_j_m3, 0.0, tol().energy_j_m3);
    near(heated.overstress_dissipation_j_m3, 0.0, tol().energy_j_m3);
    near(heated.thermal.mechanical.energy.numerical_dissipation_j_m3, 0.0, tol().energy_j_m3);
}

#[test]
fn thermal_cyclic_intervals_refuse_foreign_sources_clocks_and_domains_atomically() {
    let fixture = Fixture::new(rate_properties());
    let source = fixture.rate(true);
    let cold = source.at_uniaxial_temperatures(300.0, 300.0, 0.01).unwrap();
    let loaded = cold.advance_with_energy(0.0, &cold.initial_state(), 0.01, tol()).unwrap();
    let original = loaded.state.clone();
    let hot = source.at_uniaxial_temperatures(300.0, 315.0, 0.02).unwrap();
    cold.admit_successor(&hot).unwrap();
    let continued = hot.advance_with_energy(0.01, &original, 0.01, tol()).unwrap();
    assert!(hot.response(0.01, &continued.state).is_err());
    assert!(cold.admit_successor(
        &source.at_uniaxial_temperatures(301.0, 315.0, 0.02).unwrap()
    ).is_err());
    for dt in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(source.at_uniaxial_temperatures(300.0, 310.0, dt).is_err());
    }
    assert!(source.at_uniaxial_temperatures(300.0, 321.0, 0.01).is_err());
    assert!(hot.advance_with_energy(0.01, &original, 0.2, tol()).is_err());
    assert!(cold.validate_history(0.01, &original, f64::NAN).is_err());
    assert_eq!(loaded.state, original);

    let mut foreign_data = rate_properties();
    foreign_data.iter_mut().find(|p| p.key.name() == CP_NAME).unwrap()
        .provenance.source = "different original caloric measurement".into();
    let foreign = Fixture::new(foreign_data);
    let foreign_interval = foreign.rate(true)
        .at_uniaxial_temperatures(300.0, 315.0, 0.02).unwrap();
    assert_ne!(cold.caloric_curve_identity(), foreign_interval.caloric_curve_identity());
    assert!(cold.admit_successor(&foreign_interval).is_err());
    assert_eq!(cold.material().card_identity(), source.base().reference().material().card_identity());
    assert_eq!(continued.state.mechanical.tensor_history().backstress_pa.len(), 2);
}

#[test]
fn explicitly_absent_cyclic_branches_still_return_against_the_hot_yield_surface() {
    let mut data = properties([0.0, 0.0]);
    for name in [Q.to_owned(), schema()[0].modulus.clone(), schema()[1].modulus.clone()] {
        let property = data.iter_mut().find(|p| p.key.name() == name).unwrap();
        property.value = PropertyValue::Scalar { value: 0.0, dims: Pressure::DIMS };
    }
    affine(&mut data, YIELD_STRESS_PROPERTY, SY, -80_000.0);
    let fixture = Fixture::new(data);
    let source = fixture.binding();
    let cold = source.at_uniaxial_temperatures(300.0, 300.0, 0.1).unwrap();
    let initial = cold.initial_state();
    let loaded = cold.advance_with_energy(0.0, &initial, 0.012, tol()).unwrap();
    let hot = source.at_uniaxial_temperatures(300.0, 315.0, 0.1).unwrap();
    cold.admit_successor(&hot).unwrap();
    let heated = hot.advance_with_energy(0.012, &loaded.state, 0.012, tol()).unwrap();
    assert!(heated.state.mechanical.equivalent_plastic_strain()
        > loaded.state.mechanical.equivalent_plastic_strain());
    assert!(heated.stress_pa < loaded.stress_pa);
    assert_eq!(heated.state.mechanical.backstress_pa(), [0.0; 2]);
    assert_eq!(heated.recovery_dissipation_j_m3, 0.0);
    near(heated.thermal.mechanical.energy.residual_j_m3, 0.0, tol().energy_j_m3);
}
