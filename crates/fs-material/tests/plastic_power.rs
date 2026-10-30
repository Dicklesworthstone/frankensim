//! Independent manufactured power-law creep/relaxation and native tensor checks.
use fs_material::elastic::IsotropicElastic;
use fs_material::plastic::J2Plasticity;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::rate::{LinearOverstressJ2, PowerLawOverstressJ2};
use fs_material::{SmallStrainLaw, Voigt};

fn backbone() -> J2Plasticity {
    J2Plasticity::new_mixed(IsotropicElastic::new(1.2e9, 0.2, 0.05).unwrap(),
        3e6, 1e8, 7e7).unwrap()
}
fn limits() -> EnergyTolerances {
    EnergyTolerances { energy_j_m3: 2e-7, yield_pa: 2e-6 }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance, "{a} != {b}, allowance {tolerance}");
}
fn axial(e: f64) -> Voigt { [e, -0.5*e, -0.5*e, 0.0, 0.0, 0.0] }

#[test]
fn g0_n_one_preserves_the_original_linear_history_tangent_and_energy_exactly() {
    let original = LinearOverstressJ2::new(backbone(), 6e8, 0.25).unwrap();
    let strain = [0.01, -0.003, -0.002, 0.002, -0.001, 0.0004];
    let history = original.update_state(&axial(0.004), &original.initial_state());
    for scale in [1e-6, 2e6, 1e24] {
        let nonlinear = PowerLawOverstressJ2::new(backbone(), 6e8, scale, 1.0, 0.25).unwrap();
        assert_eq!(original.stress(&strain, &history), nonlinear.stress(&strain, &history));
        assert_eq!(original.tangent(&strain, &history), nonlinear.tangent(&strain, &history));
        assert_eq!(original.update_state(&strain, &history), nonlinear.update_state(&strain, &history));
        assert_eq!(original.advance_with_energy(&axial(0.004), &history, &strain, limits()),
            nonlinear.advance_with_energy(&axial(0.004), &history, &strain, limits()));
    }
}

#[test]
fn g1_quadratic_closed_form_checks_mixed_hardening_and_physical_heat() {
    let material = backbone();
    let (_, mu) = material.elastic.lame();
    let (eta, scale, dt) = (6e8, 4e6, 0.25);
    let law = PowerLawOverstressJ2::new(material.clone(), eta, scale, 2.0, dt).unwrap();
    let strain = axial(0.01);
    let next = law.advance_with_energy(&[0.0; 6], &law.initial_state(), &strain, limits()).unwrap();
    let a = 3.0*mu + material.hardening + material.kinematic_hardening;
    let f_trial = 3.0*mu*strain[0] - material.yield_stress;
    // Solve f + A*dt*f^2/(eta*s0) = f_trial algebraically, with no root iteration.
    let f = 2.0*f_trial / (1.0 + (1.0 + 4.0*a*dt*f_trial/(eta*scale)).sqrt());
    let delta = dt*f*f/(eta*scale);
    close(next.state.alpha, delta, 3e-17);
    close(next.state.plastic_strain[0], delta, 3e-17);
    close(next.stress_pa[0], 2.0*mu*(strain[0]-delta), 1e-7);
    close(material.yield_function_with_state(&next.stress_pa, &next.state), f, 1e-7);
    close(next.energy.physical_dissipation_j_m3, (material.yield_stress+f)*delta, 1e-8);
    close(next.energy.final_state.hardening_j_m3, 0.5*material.hardening*delta*delta, 1e-9);
    close(next.energy.final_state.kinematic_hardening_j_m3,
        0.5*material.kinematic_hardening*delta*delta, 1e-9);
    assert!(next.energy.residual_j_m3.abs() < limits().energy_j_m3);
}

#[test]
fn g1_constant_stress_creep_retains_actual_time_and_accumulated_history() {
    let material = J2Plasticity::new(backbone().elastic, 3e6, 0.0).unwrap();
    let (_, mu) = material.elastic.lame();
    let (eta, scale, f, dt) = (1e9, 5e6, 4e6, 0.05);
    for n in [1.4, 3.2, 8.0] {
        let law = PowerLawOverstressJ2::new(material.clone(), eta, scale, n, dt).unwrap();
        let rate = scale/eta * f64::powf(f/scale, n);
        let q = material.yield_stress+f;
        let mut history = law.initial_state();
        let mut old = axial(q/(3.0*mu));
        for step in 1..=12 {
            // Prescribe the independently known strain history for constant q.
            let strain = axial(q/(3.0*mu)+rate*dt*step as f64);
            let next = law.advance_with_energy(&old, &history, &strain, limits()).unwrap();
            close(next.state.alpha, rate*dt*step as f64, 8e-17);
            close(next.stress_pa[0], 2.0*q/3.0, 2e-7);
            close(next.energy.physical_dissipation_j_m3, q*rate*dt, 1e-8);
            history = next.state;
            old = strain;
        }
    }
}

#[test]
fn g1_fixed_strain_nonlinear_relaxation_converges_to_analytical_solution() {
    let material = backbone();
    let (_, mu) = material.elastic.lame();
    let (eta, scale, time) = (8e8, 4e6, 0.5);
    let strain = axial(0.01);
    let a = 3.0*mu+material.hardening+material.kinematic_hardening;
    let f0 = 3.0*mu*strain[0]-material.yield_stress;
    let exact = (f0-f0/(1.0+a*time*f0/(eta*scale)))/a;
    let mut errors = Vec::new();
    for steps in [8, 16, 32] {
        let dt = time/steps as f64;
        let law = PowerLawOverstressJ2::new(material.clone(), eta, scale, 2.0, dt).unwrap();
        let mut history = law.initial_state();
        let mut heat = 0.0;
        let mut numerical = 0.0;
        for _ in 0..steps {
            let next = law.advance_with_energy(&strain, &history, &strain, limits()).unwrap();
            close(next.energy.mechanical_work_j_m3, 0.0, 0.0);
            assert!(next.energy.physical_dissipation_j_m3 > 0.0);
            heat += next.energy.physical_dissipation_j_m3;
            numerical += next.energy.numerical_dissipation_j_m3;
            history = next.state;
        }
        let storage = 1.5*mu*(strain[0]-history.alpha).powi(2)
            +0.5*(material.hardening+material.kinematic_hardening)*history.alpha.powi(2);
        close(storage+heat+numerical, 1.5*mu*strain[0].powi(2), 1e-6);
        errors.push((history.alpha-exact).abs());
    }
    assert!(errors[1] < 0.6*errors[0]);
    assert!(errors[2] < 0.6*errors[1]);
}

#[test]
fn g0_consistent_tensor_tangent_includes_power_derivative_and_shear_weights() {
    for exponent in [1.25, 2.0, 4.5] {
        let law = PowerLawOverstressJ2::new(backbone(), 2e8, 4e6, exponent, 0.2).unwrap();
        let state = law.update_state(&[0.004, -0.001, -0.001, 0.002, 0.0, 0.0],
            &law.initial_state());
        let strain = [0.009, -0.003, -0.002, 0.004, 0.0005, -0.001];
        let tangent = law.tangent(&strain, &state);
        for j in 0..6 {
            let mut plus = strain; plus[j] += 1e-8;
            let mut minus = strain; minus[j] -= 1e-8;
            let sp = law.stress(&plus, &state);
            let sm = law.stress(&minus, &state);
            for i in 0..6 { close((sp[i]-sm[i])/2e-8, tangent[i][j], 2.0); }
        }
        // Doubling eta and physical dt leaves the complete flow law invariant.
        let scaled = PowerLawOverstressJ2::new(backbone(), 4e8, 4e6, exponent, 0.4).unwrap();
        for i in 0..6 { close(law.stress(&strain, &state)[i], scaled.stress(&strain, &state)[i], 1e-7); }
    }
}

#[test]
fn g0_elastic_domain_invalid_sources_and_unrepresentable_return_are_explicit() {
    let law = PowerLawOverstressJ2::new(backbone(), 2e8, 4e6, 2.5, 0.2).unwrap();
    for strain in [[0.0; 6], axial(0.001), [0.002, 0.002, 0.002, 0.0, 0.0, 0.0]] {
        let next = law.advance_with_energy(&[0.0; 6], &law.initial_state(), &strain, limits()).unwrap();
        assert_eq!(next.state, law.initial_state());
        assert_eq!(next.energy.physical_dissipation_j_m3, 0.0);
        assert_eq!(law.tangent(&strain, &law.initial_state()), backbone().elastic.stiffness());
    }
    for (scale, exponent) in [(0.0, 2.0), (-1.0, 2.0), (1.0, 0.5), (1.0, f64::INFINITY)] {
        assert!(PowerLawOverstressJ2::new(backbone(), 2e8, scale, exponent, 0.2).is_err());
    }
    let enormous = PowerLawOverstressJ2::new(backbone(), 2e8, 1e8, 1e308, 0.2).unwrap();
    assert!(enormous.advance_with_energy(&[0.0; 6], &enormous.initial_state(),
        &axial(0.01), limits()).is_err());
}
