//! Independent linear-overstress references, tangent and physical heat checks.
use fs_material::elastic::IsotropicElastic;
use fs_material::plastic::J2Plasticity;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::rate::LinearOverstressJ2;
use fs_material::{SmallStrainLaw, Voigt};

fn law() -> J2Plasticity {
    J2Plasticity::new(IsotropicElastic::new(1.2e9, 0.2, 0.05).unwrap(), 3e6, 1e8).unwrap()
}
fn limits() -> EnergyTolerances {
    EnergyTolerances {
        energy_j_m3: 1e-7,
        yield_pa: 1e-6,
    }
}
fn close(a: f64, b: f64, allowance: f64) {
    assert!(
        (a - b).abs() <= allowance,
        "{a} != {b}, allowance {allowance}"
    );
}

#[test]
fn g1_radial_return_has_independent_rate_and_heat_with_energetic_hardening() {
    let backbone = law();
    let (_, mu) = backbone.elastic.lame();
    let eta = 6e8;
    let dt = 0.25;
    let rate = LinearOverstressJ2::new(backbone.clone(), eta, dt).unwrap();
    let strain = [0.01, -0.005, -0.005, 0.0, 0.0, 0.0];
    let old = backbone.initial_state();
    let next = rate
        .advance_with_energy(&[0.0; 6], &old, &strain, limits())
        .unwrap();
    // Pure deviatoric axial loading gives q_trial=3 mu eps_x and p_x=alpha.
    let delta =
        (3.0 * mu * strain[0] - backbone.yield_stress) / (3.0 * mu + backbone.hardening + eta / dt);
    let elastic_energy = 1.5 * mu * (strain[0] - delta).powi(2);
    let hardening_energy = 0.5 * backbone.hardening * delta * delta;
    let heat = backbone.yield_stress * delta + eta / dt * delta * delta;
    close(next.state.alpha, delta, 1e-17);
    close(next.state.plastic_strain[0], delta, 1e-17);
    close(next.stress_pa[0], 2.0 * mu * (strain[0] - delta), 1e-8);
    close(
        backbone.yield_function(&next.stress_pa, next.state.alpha),
        eta / dt * delta,
        1e-8,
    );
    close(next.energy.final_state.elastic_j_m3, elastic_energy, 1e-9);
    close(
        next.energy.final_state.hardening_j_m3,
        hardening_energy,
        1e-9,
    );
    close(next.energy.physical_dissipation_j_m3, heat, 1e-9);
    close(
        next.energy.numerical_dissipation_j_m3,
        elastic_energy + hardening_energy,
        1e-9,
    );
    assert!(next.energy.residual_j_m3.abs() < limits().energy_j_m3);
    assert!(
        backbone
            .validate_committed_history(&strain, &next.state, 1e-6)
            .is_err()
    );
    assert!(
        backbone
            .validate_overstress_history(&strain, &next.state)
            .is_ok()
    );
    assert_eq!(old, backbone.initial_state());
}

#[test]
fn g1_fixed_strain_relaxes_on_actual_time_and_physical_loss_heats() {
    let backbone = law();
    let (_, mu) = backbone.elastic.lame();
    let eta = 8e8;
    let strain = [0.01, -0.005, -0.005, 0.0, 0.0, 0.0];
    let total_time = 0.5;
    let f0 = 3.0 * mu * strain[0] - backbone.yield_stress;
    let stiffness = 3.0 * mu + backbone.hardening;
    let exact = f0 / stiffness * (1.0 - (-stiffness * total_time / eta).exp());
    let mut errors = Vec::new();
    for steps in [8, 16, 32] {
        let dt = total_time / steps as f64;
        let rate = LinearOverstressJ2::new(backbone.clone(), eta, dt).unwrap();
        // An instantaneous elastic preparation can carry positive overstress.
        let mut history = backbone.initial_state();
        let mut heat = 0.0;
        let mut numerical = 0.0;
        for _ in 0..steps {
            let next = rate
                .advance_with_energy(&strain, &history, &strain, limits())
                .unwrap();
            close(next.energy.mechanical_work_j_m3, 0.0, 0.0);
            assert!(next.energy.physical_dissipation_j_m3 > 0.0);
            heat += next.energy.physical_dissipation_j_m3;
            numerical += next.energy.numerical_dissipation_j_m3;
            history = next.state;
        }
        let discrete = f0 / stiffness * (1.0 - (1.0 + stiffness * dt / eta).powi(-steps));
        close(history.alpha, discrete, 2e-17);
        let stored = 1.5 * mu * (strain[0] - history.alpha).powi(2)
            + 0.5 * backbone.hardening * history.alpha * history.alpha;
        close(
            stored + heat + numerical,
            1.5 * mu * strain[0] * strain[0],
            2e-8,
        );
        errors.push((history.alpha - exact).abs());
    }
    assert!(errors[1] < 0.55 * errors[0]);
    assert!(errors[2] < 0.55 * errors[1]);
}

#[test]
fn g0_discrete_tangent_and_material_time_scaling_match_native_stress() {
    let backbone = law();
    let rate = LinearOverstressJ2::new(backbone.clone(), 2e8, 0.2).unwrap();
    let state = rate.update_state(
        &[0.004, -0.001, -0.001, 0.002, 0.0, 0.0],
        &rate.initial_state(),
    );
    let strain: Voigt = [0.009, -0.003, -0.002, 0.004, 0.0005, -0.001];
    let tangent = rate.tangent(&strain, &state);
    for j in 0..6 {
        let mut plus = strain;
        plus[j] += 1e-8;
        let mut minus = strain;
        minus[j] -= 1e-8;
        let sp = rate.stress(&plus, &state);
        let sm = rate.stress(&minus, &state);
        for i in 0..6 {
            close((sp[i] - sm[i]) / 2e-8, tangent[i][j], 1.0);
        }
    }
    let equivalent = LinearOverstressJ2::new(backbone.clone(), 4e8, 0.4).unwrap();
    assert_eq!(
        equivalent.stress(&strain, &state),
        rate.stress(&strain, &state)
    );
    assert_eq!(
        equivalent.update_state(&strain, &state),
        rate.update_state(&strain, &state)
    );
    let slow = LinearOverstressJ2::new(backbone.clone(), 2e8, 0.02).unwrap();
    assert!(slow.update_state(&strain, &state).alpha < rate.update_state(&strain, &state).alpha);
    for (eta, dt) in [
        (0.0, 1.0),
        (-1.0, 1.0),
        (1.0, 0.0),
        (f64::MAX, f64::MIN_POSITIVE),
    ] {
        assert!(LinearOverstressJ2::new(backbone.clone(), eta, dt).is_err());
    }
}
