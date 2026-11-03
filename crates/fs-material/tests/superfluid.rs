//! G0/G1/G3: two-fluid thermodynamics, currents and physical source balances.
//! The round-number states are manufactured numerical fixtures, not helium data.

use fs_material::fluid_acoustics::FluidAcoustics;
use fs_material::superfluid::{TwoFluidProperties, TwoFluidState, TwoFluidVelocity};

fn properties() -> TwoFluidProperties {
    TwoFluidProperties {
        temperature_k: 2.0,
        density_kg_m3: 4.0,
        superfluid_fraction: 0.75,
        entropy_j_kg_k: 4.0,
        specific_heat_cp_j_kg_k: 16.0,
        adiabatic_sound_speed_m_s: 10.0,
        thermal_expansion_per_k: 0.0,
        normal_viscosity_pa_s: 0.03,
    }
}

fn state() -> TwoFluidState {
    TwoFluidState::try_new(properties()).unwrap()
}

fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 4e-13 * scale.max(expected.abs()),
        "{actual:e} != {expected:e}"
    );
}

fn close_vector(actual: [f64; 3], expected: [f64; 3], scale: f64) {
    for i in 0..3 {
        close(actual[i], expected[i], scale);
    }
}

#[test]
fn nonzero_expansion_closes_thermodynamics_and_two_distinct_modes() {
    let mut p = properties();
    p.thermal_expansion_per_k = 0.4;
    let material = TwoFluidState::try_new(p).unwrap();
    let cp = material.specific_heat_cp_j_kg_k();
    let cv = material.specific_heat_cv_j_kg_k();
    let kt = material.isothermal_compressibility_per_pa();
    close(
        cp - cv,
        p.temperature_k * p.thermal_expansion_per_k.powi(2) / (p.density_kg_m3 * kt),
        cp,
    );
    close(
        material.adiabatic_compressibility_per_pa() / kt,
        cv / cp,
        1.0,
    );

    // Independent energy-coordinate two-fluid wave matrix G G^T. Its
    // eigenvalue trace/determinant must equal the returned squared speeds.
    let c = p.adiabatic_sound_speed_m_s;
    let c20 = p.entropy_j_kg_k * (3.0 * p.temperature_k / cp).sqrt();
    let coupling = c * p.thermal_expansion_per_k * (p.temperature_k / cp).sqrt() * c20;
    let modes = material.sound_modes();
    let fast2 = modes.fast_m_s.powi(2);
    let slow2 = modes.slow_m_s.powi(2);
    close(
        fast2 + slow2,
        c * c + coupling * coupling + c20 * c20,
        fast2,
    );
    close(fast2 * slow2, c * c * c20 * c20, fast2 * slow2);
    assert!(modes.fast_m_s > c);
    assert!(modes.slow_m_s < c20);

    p.thermal_expansion_per_k *= -1.0;
    let opposite = TwoFluidState::try_new(p).unwrap();
    assert_eq!(
        opposite.sound_modes(),
        modes,
        "expansion sign changes polarization, not speeds"
    );
}

#[test]
fn zero_expansion_and_small_entropy_mode_are_resolved_without_cancellation() {
    let material = state();
    assert_eq!(material.specific_heat_cv_j_kg_k(), 16.0);
    close(material.sound_modes().fast_m_s, 10.0, 10.0);
    close(material.sound_modes().slow_m_s, 6.0_f64.sqrt(), 3.0);
    close(material.normal_kinematic_viscosity_m2_s(), 0.03, 0.03);

    let mut p = properties();
    p.entropy_j_kg_k = 1e-12;
    let weak = TwoFluidState::try_new(p).unwrap();
    close(
        weak.sound_modes().slow_m_s,
        1e-12 * (3.0_f64 / 8.0).sqrt(),
        1e-12,
    );
    assert!(weak.sound_modes().slow_m_s > 0.0);

    // The uncoupled density mode need not be the faster mode for authored data.
    p.entropy_j_kg_k = 100.0;
    let reversed = TwoFluidState::try_new(p).unwrap();
    close(reversed.sound_modes().slow_m_s, 10.0, 10.0);
    assert!(reversed.sound_modes().fast_m_s > 50.0);
}

#[test]
fn counterflow_carries_heat_and_entropy_with_no_total_mass_flow() {
    let material = state();
    let relative = [0.4, -0.2, 0.1];
    let velocity = material
        .velocities_from_barycentric([0.0; 3], relative)
        .unwrap();
    let currents = material.currents(velocity).unwrap();
    close_vector(currents.mass_kg_m2_s, [0.0; 3], 1.0);
    close_vector(currents.heat_w_m2, relative.map(|w| 24.0 * w), 10.0);
    close_vector(
        currents.entropy_w_m2_k.map(|s| 2.0 * s),
        currents.heat_w_m2,
        10.0,
    );
    let energy = material.kinetic_energy_density(velocity).unwrap();
    close(energy.counterflow_j_m3, 0.5 * (3.0 / 4.0) * 0.21, 1.0);
    assert_eq!(energy.barycentric_j_m3, 0.0);

    let boost = [0.7, -0.4, 0.2];
    let moved = material
        .velocities_from_barycentric(boost, relative)
        .unwrap();
    let moved_currents = material.currents(moved).unwrap();
    close_vector(moved_currents.heat_w_m2, currents.heat_w_m2, 10.0);
    close_vector(moved_currents.mass_kg_m2_s, boost.map(|v| 4.0 * v), 4.0);
    // Direct component kinetic storage, independently of the mean/relative split.
    let component_energy = (0..3)
        .map(|i| 0.5 * moved.normal_m_s[i].powi(2) + 1.5 * moved.superfluid_m_s[i].powi(2))
        .sum::<f64>();
    close(
        material.kinetic_energy_density(moved).unwrap().total_j_m3,
        component_energy,
        4.0,
    );
}

#[test]
fn temperature_gradient_drives_opposed_components_and_fountain_balance() {
    let material = state();
    let dt = [0.1, -0.3, 0.0];
    let driven = material
        .thermomechanical_acceleration([0.0; 3], dt)
        .unwrap();
    close_vector(driven.normal_m_s2, dt.map(|g| -12.0 * g), 4.0);
    close_vector(driven.superfluid_m_s2, dt.map(|g| 4.0 * g), 4.0);
    for i in 0..3 {
        close(
            driven.normal_m_s2[i] + 3.0 * driven.superfluid_m_s2[i],
            0.0,
            4.0,
        );
    }
    // grad(p) = rho*s*grad(T) makes the superfluid chemical potential uniform.
    let fountain = material
        .thermomechanical_acceleration(dt.map(|g| 16.0 * g), dt)
        .unwrap();
    close_vector(fountain.superfluid_m_s2, [0.0; 3], 4.0);
    close_vector(fountain.normal_m_s2, dt.map(|g| -16.0 * g), 5.0);
}

#[test]
fn linear_friction_preserves_momentum_and_transfers_relative_energy_once() {
    let material = state();
    let mean = [0.7, -0.4, 0.2];
    let relative = [0.4, -0.2, 0.1];
    let initial = material
        .velocities_from_barycentric(mean, relative)
        .unwrap();
    let gamma = material.reduced_density_kg_m3() * 2.0;
    let whole = material.relax_linear_friction(initial, gamma, 0.3).unwrap();
    close_vector(
        material.barycentric_velocity_m_s(whole.velocity).unwrap(),
        mean,
        1.0,
    );
    close_vector(
        material.counterflow_velocity_m_s(whole.velocity).unwrap(),
        relative.map(|w| w * (-0.6_f64).exp()),
        1.0,
    );
    assert_eq!(
        whole.superfluid_impulse_kg_m2_s,
        whole.normal_impulse_kg_m2_s.map(|x| -x)
    );
    let before = material.kinetic_energy_density(initial).unwrap();
    let after = material.kinetic_energy_density(whole.velocity).unwrap();
    close(
        before.total_j_m3 - after.total_j_m3,
        whole.dissipated_energy_j_m3,
        4.0,
    );
    close(
        whole.dissipated_energy_j_m3,
        before.counterflow_j_m3 * (1.0 - (-1.2_f64).exp()),
        1.0,
    );
    let a = material.relax_linear_friction(initial, gamma, 0.1).unwrap();
    let b = material
        .relax_linear_friction(a.velocity, gamma, 0.2)
        .unwrap();
    close_vector(b.velocity.normal_m_s, whole.velocity.normal_m_s, 1.0);
    close_vector(
        b.velocity.superfluid_m_s,
        whole.velocity.superfluid_m_s,
        1.0,
    );
    close(
        a.dissipated_energy_j_m3 + b.dissipated_energy_j_m3,
        whole.dissipated_energy_j_m3,
        1.0,
    );

    let tiny = material
        .relax_linear_friction(initial, gamma, 1e-20)
        .unwrap();
    assert!(
        tiny.dissipated_energy_j_m3 > 0.0,
        "expm1 retains a sub-ulp loss fraction"
    );
    let stiff = material
        .relax_linear_friction(initial, gamma, 1000.0)
        .unwrap();
    close_vector(stiff.velocity.normal_m_s, mean, 1.0);
    close_vector(stiff.velocity.superfluid_m_s, mean, 1.0);
    close(stiff.dissipated_energy_j_m3, before.counterflow_j_m3, 1.0);
}

#[test]
fn gorter_mellink_has_cubic_drag_and_exact_bounded_heat() {
    let material = state();
    let initial = material
        .velocities_from_barycentric([0.0; 3], [2.0, 0.0, 0.0])
        .unwrap();
    let coefficient = 0.25;
    let step = material
        .relax_gorter_mellink(initial, coefficient, 0.5)
        .unwrap();
    let expected_w = 2.0 / 5.0_f64.sqrt();
    close(
        material.counterflow_velocity_m_s(step.velocity).unwrap()[0],
        expected_w,
        2.0,
    );
    close(step.dissipated_energy_j_m3, 1.5 * 0.8, 1.5);
    close_vector(
        material.currents(step.velocity).unwrap().mass_kg_m2_s,
        [0.0; 3],
        1.0,
    );
    let after = material
        .kinetic_energy_density(step.velocity)
        .unwrap()
        .total_j_m3;
    close(after + step.dissipated_energy_j_m3, 1.5, 1.5);

    let short = material
        .relax_gorter_mellink(initial, coefficient, 1e-8)
        .unwrap();
    // F_n = -A*rho_n*rho_s*|w|²*w independently predicts initial normal force.
    let observed_force = short.normal_impulse_kg_m2_s[0] / 1e-8;
    assert!((observed_force / -6.0 - 1.0).abs() < 1e-7);
    let a = material
        .relax_gorter_mellink(initial, coefficient, 0.2)
        .unwrap();
    let b = material
        .relax_gorter_mellink(a.velocity, coefficient, 0.3)
        .unwrap();
    close_vector(b.velocity.normal_m_s, step.velocity.normal_m_s, 2.0);
    close(
        a.dissipated_energy_j_m3 + b.dissipated_energy_j_m3,
        step.dissipated_energy_j_m3,
        1.5,
    );
}

#[test]
fn scalar_acoustics_requires_an_explicit_uncoupled_projection() {
    let material = state();
    let first = material.first_sound_projection().unwrap();
    assert_eq!(first.source(), &material);
    let acoustic = first.acoustic_state().unwrap();
    assert_eq!(acoustic.density_kg_m3(), 4.0);
    assert_eq!(acoustic.sound_speed_m_s(), 10.0);
    assert_eq!(acoustic.adiabatic_bulk_modulus_pa(), 400.0);
    let mut p = properties();
    p.thermal_expansion_per_k = 1e-12;
    assert!(
        TwoFluidState::try_new(p)
            .unwrap()
            .first_sound_projection()
            .is_err()
    );
}

#[test]
fn singular_phases_bad_states_and_negative_dissipation_are_refused() {
    for fraction in [0.0, 1.0, -0.1, 1.1, f64::NAN] {
        let mut p = properties();
        p.superfluid_fraction = fraction;
        assert!(TwoFluidState::try_new(p).is_err());
    }
    for value in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        let mut p = properties();
        p.temperature_k = value;
        assert!(TwoFluidState::try_new(p).is_err());
    }
    let mut p = properties();
    p.normal_viscosity_pa_s = -1.0;
    assert!(TwoFluidState::try_new(p).is_err());
    p = properties();
    p.thermal_expansion_per_k = f64::INFINITY;
    assert!(TwoFluidState::try_new(p).is_err());
    let material = state();
    let initial = TwoFluidVelocity {
        normal_m_s: [1.0; 3],
        superfluid_m_s: [0.0; 3],
    };
    for bad in [-1.0, f64::INFINITY, f64::NAN] {
        assert!(material.relax_linear_friction(initial, bad, 1.0).is_err());
        assert!(material.relax_gorter_mellink(initial, bad, 1.0).is_err());
        assert!(material.relax_linear_friction(initial, 1.0, bad).is_err());
    }
    assert_eq!(
        material
            .relax_linear_friction(initial, 0.0, 1.0)
            .unwrap()
            .velocity,
        initial
    );
    assert_eq!(
        material
            .relax_gorter_mellink(initial, 1.0, 0.0)
            .unwrap()
            .velocity,
        initial
    );
    let bad = TwoFluidVelocity {
        normal_m_s: [f64::NAN; 3],
        ..initial
    };
    assert!(material.currents(bad).is_err());
    assert!(material.relax_linear_friction(bad, 0.0, 0.0).is_err());
    // Finite inputs whose rate product overflows must not masquerade as an
    // exactly fully damped state with an invented heat deposit.
    assert!(
        material
            .relax_linear_friction(initial, f64::MAX, 10.0)
            .is_err()
    );
    assert!(
        material
            .relax_gorter_mellink(initial, f64::MAX, 10.0)
            .is_err()
    );
}
