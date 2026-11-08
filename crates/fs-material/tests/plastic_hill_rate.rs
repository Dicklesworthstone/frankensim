//! Synthetic kinetic, energetic and reduction regressions. Explicit coefficients
//! exercise the model equations; they do not calibrate a measured alloy.

use fs_material::IsotropicElastic;
use fs_material::plastic::J2Plasticity;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::hill::oriented::OrientedHill48;
use fs_material::plastic::hill::rate::LinearOverstressHill48;
use fs_material::plastic::hill::{Hill48, Hill48Coefficients};
use fs_material::plastic::rate::VoceOverstressJ2;
use fs_material::plastic::reduction::{PlaneStressReduction, ReductionConfig};
use fs_material::plastic::voce::VoceJ2;
use fs_material::tensor::{contract, rotate, trace};
use fs_material::{SmallStrainLaw, Voigt};

const E: f64 = 1200.0;
const NU: f64 = 0.2;
const Y: f64 = 3.0;
const H: f64 = 40.0;
const C: f64 = 100.0;
const Q: f64 = 5.0;
const B: f64 = 180.0;
const ETA: f64 = 8.0;
const DT: f64 = 0.04;
const TOL: EnergyTolerances = EnergyTolerances {
    energy_j_m3: 1e-10,
    yield_pa: 1e-9,
};

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.16e}, expected={expected:.16e}, allowance={tolerance:.3e}"
    );
}

fn coefficients() -> Hill48Coefficients {
    Hill48Coefficients {
        f: 0.25,
        g: 0.9,
        h: 0.45,
        l: 1.4,
        m: 2.1,
        n: 1.8,
    }
}

fn native(coefficients: Hill48Coefficients, saturation: f64) -> Hill48 {
    Hill48::new(
        IsotropicElastic::new(E, NU, 0.05).unwrap(),
        Y,
        H,
        coefficients,
    )
    .unwrap()
    .with_linear_kinematic_hardening(C)
    .unwrap()
    .with_voce_isotropic_hardening(saturation, B)
    .unwrap()
}

#[test]
fn held_anisotropic_shear_follows_independent_discrete_relaxation_and_heat_balance() {
    let coefficients = coefficients();
    for (axis, coefficient) in [
        (3, coefficients.n),
        (4, coefficients.l),
        (5, coefficients.m),
    ] {
        let backbone = native(coefficients, 0.0);
        let law = LinearOverstressHill48::new(backbone.clone(), ETA, DT).unwrap();
        let mut strain = [0.0; 6];
        strain[axis] = 0.012;
        let first = law
            .advance_with_energy(&[0.0; 6], &law.initial_state(), &strain, TOL)
            .unwrap();
        assert!(first.response.yield_residual_pa > TOL.yield_pa);
        assert!(
            backbone
                .validate_committed_history(&strain, &first.response.state, TOL.yield_pa)
                .is_err()
        );
        law.validate_history(&strain, &first.response.state)
            .unwrap();
        let modulus = ETA / DT;
        let reduction = H + (E / (1.0 + NU) + C) * coefficient;
        let mut state = first.response.state;
        let mut excess = first.response.yield_residual_pa;
        let mut stress = first.response.stress_pa[axis];
        for _ in 0..6 {
            let retained = state.clone();
            let step = law
                .advance_with_energy(&strain, &state, &strain, TOL)
                .unwrap();
            assert_eq!(state, retained);
            assert_eq!(
                step,
                law.advance_with_energy(&strain, &state, &strain, TOL)
                    .unwrap()
            );
            let delta = excess / (reduction + modulus);
            near(step.response.plastic_increment, delta, 3e-15);
            near(
                step.response.yield_residual_pa,
                excess * modulus / (reduction + modulus),
                3e-11,
            );
            assert!(step.response.stress_pa[axis] < stress);
            near(step.energy.mechanical_work_j_m3, 0.0, 0.0);
            near(
                step.energy.physical_dissipation_j_m3,
                (Y + modulus * delta) * delta,
                3e-13,
            );
            near(
                -step.energy.stored_change_j_m3(),
                step.energy.physical_dissipation_j_m3 + step.energy.numerical_dissipation_j_m3,
                3e-13,
            );
            assert!(step.energy.physical_dissipation_j_m3 > 0.0);
            assert!(step.energy.numerical_dissipation_j_m3 >= 0.0);
            near(trace(&step.response.state.plastic_strain), 0.0, 3e-17);
            excess = step.response.yield_residual_pa;
            stress = step.response.stress_pa[axis];
            state = step.response.state;
        }
    }
}

#[test]
fn anisotropic_shear_matches_prescribed_multiplier_and_exact_viscous_voce_slope() {
    let coefficients = coefficients();
    for duration in [0.004, 0.04, 0.4] {
        let law = LinearOverstressHill48::new(native(coefficients, Q), ETA, duration).unwrap();
        for (axis, coefficient) in [
            (3, coefficients.n),
            (4, coefficients.l),
            (5, coefficients.m),
        ] {
            // Construct a target from the independent associated-flow equations.
            let alpha = 0.002;
            let modulus = ETA / duration;
            let radius = Y + H * alpha + Q * -(-B * alpha).exp_m1() + modulus * alpha;
            let plastic = alpha * (coefficient / 2.0).sqrt();
            let stress = radius / (2.0 * coefficient).sqrt() + C * plastic;
            let twice_mu = E / (1.0 + NU);
            let mut strain = [0.0; 6];
            strain[axis] = plastic + stress / twice_mu;
            let step = law
                .advance_with_energy(&[0.0; 6], &law.initial_state(), &strain, TOL)
                .unwrap();
            near(step.response.state.alpha, alpha, 2e-15);
            near(step.response.state.plastic_strain[axis], plastic, 2e-15);
            near(step.response.stress_pa[axis], stress, 3e-12);
            near(step.response.yield_residual_pa, modulus * alpha, 3e-12);
            let endpoint_modulus = C + (H + Q * B * (-B * alpha).exp() + modulus) / coefficient;
            near(
                step.response.tangent_pa[axis][axis],
                twice_mu * endpoint_modulus / (twice_mu + endpoint_modulus),
                3e-9,
            );
            let original = law
                .local_law()
                .stored_energy(&strain, &step.response.state)
                .unwrap();
            assert_eq!(step.energy.final_state, original);
            near(
                step.energy.physical_dissipation_j_m3,
                Y * alpha + modulus * alpha * alpha,
                3e-13,
            );
        }
    }
}

#[test]
fn isotropic_limit_matches_joint_voce_j2_history_tangent_and_energy_during_cycles_and_holds() {
    let backbone = native(Hill48Coefficients::isotropic(), Q);
    let law = LinearOverstressHill48::new(backbone.clone(), ETA, DT).unwrap();
    let j2 = VoceOverstressJ2::new_linear(
        VoceJ2::new(
            J2Plasticity::new_mixed(backbone.elastic().clone(), Y, H, 1.5 * C).unwrap(),
            Q,
            B,
        )
        .unwrap(),
        ETA,
        DT,
    )
    .unwrap();
    let path = [
        [0.0001, -0.00005, -0.00005, 0.0002, 0.0, 0.0],
        [0.008, -0.003, -0.005, 0.002, 0.0, 0.001],
        [0.008, -0.003, -0.005, 0.002, 0.0, 0.001],
        [-0.006, 0.002, 0.004, -0.003, 0.002, 0.001],
        [0.01, -0.004, -0.006, 0.005, -0.002, 0.0],
        [0.0; 6],
    ];
    let mut hs = law.initial_state();
    let mut js = j2.initial_state();
    let mut old = [0.0; 6];
    for strain in path {
        let h = law.advance_with_energy(&old, &hs, &strain, TOL).unwrap();
        let j = j2.advance_with_energy(&old, &js, &strain, TOL).unwrap();
        let tangent = j2.tangent(&strain, &js);
        for i in 0..6 {
            near(h.response.stress_pa[i], j.stress_pa[i], 3e-11);
            near(
                h.response.state.plastic_strain[i],
                j.state.plastic_strain[i],
                3e-15,
            );
            for k in 0..6 {
                near(h.response.tangent_pa[i][k], tangent[i][k], 3e-8);
            }
        }
        near(h.response.state.alpha, j.state.alpha, 3e-15);
        near(
            h.energy.final_state.total_j_m3(),
            j.energy.final_state.total_j_m3(),
            3e-12,
        );
        near(
            h.energy.physical_dissipation_j_m3,
            j.energy.physical_dissipation_j_m3,
            3e-12,
        );
        near(
            h.energy.numerical_dissipation_j_m3,
            j.energy.numerical_dissipation_j_m3,
            3e-12,
        );
        hs = h.response.state;
        js = j.state;
        old = strain;
    }
}

#[test]
fn rotated_nonproportional_history_has_the_full_consistent_tangent_and_objective_work() {
    let rotation = [[0.8, -0.48, 0.36], [0.6, 0.64, -0.48], [0.0, 0.6, 0.8]];
    let backbone = native(coefficients(), Q);
    let local = LinearOverstressHill48::new(backbone.clone(), ETA, DT).unwrap();
    let spatial = LinearOverstressHill48::new_oriented(
        OrientedHill48::new(backbone, rotation).unwrap(),
        ETA,
        DT,
    )
    .unwrap();
    let old = [0.008, -0.003, -0.005, 0.002, 0.001, 0.0];
    let state = local.response(&old, &local.initial_state()).unwrap().state;
    let strain = [-0.006, 0.004, 0.002, -0.003, 0.002, 0.001];
    let expected = local
        .advance_with_energy(&old, &state, &strain, TOL)
        .unwrap();
    let spatial_old = rotate(&old, &rotation);
    let spatial_strain = rotate(&strain, &rotation);
    let step = spatial
        .advance_with_energy(&spatial_old, &state, &spatial_strain, TOL)
        .unwrap();
    let expected_stress = rotate(&expected.response.stress_pa, &rotation);
    for i in 0..6 {
        near(step.response.stress_pa[i], expected_stress[i], 3e-11);
        near(
            step.response.state.plastic_strain[i],
            expected.response.state.plastic_strain[i],
            3e-15,
        );
    }
    near(
        step.energy.mechanical_work_j_m3,
        expected.energy.mechanical_work_j_m3,
        3e-12,
    );
    near(
        step.energy.physical_dissipation_j_m3,
        expected.energy.physical_dissipation_j_m3,
        3e-12,
    );
    let perturbation = 1e-7;
    let weights = [1.0, 1.0, 1.0, 2.0, 2.0, 2.0];
    for j in 0..6 {
        let mut plus = spatial_strain;
        let mut minus = spatial_strain;
        plus[j] += perturbation;
        minus[j] -= perturbation;
        let p = spatial.stress(&plus, &state);
        let m = spatial.stress(&minus, &state);
        for i in 0..6 {
            near(
                step.response.tangent_pa[i][j],
                (p[i] - m[i]) / (2.0 * perturbation),
                3e-6,
            );
            near(
                weights[i] * step.response.tangent_pa[i][j],
                weights[j] * step.response.tangent_pa[j][i],
                3e-9,
            );
        }
    }
}

#[test]
fn plane_stress_reduction_keeps_viscous_history_and_condensed_tangent() {
    let rotation = [[0.8, -0.48, 0.36], [0.6, 0.64, -0.48], [0.0, 0.6, 0.8]];
    let law = LinearOverstressHill48::new_oriented(
        OrientedHill48::new(native(coefficients(), Q), rotation).unwrap(),
        ETA,
        DT,
    )
    .unwrap();
    let reduced = PlaneStressReduction::new(
        law,
        ReductionConfig {
            max_iterations: 32,
            stress_tolerance_pa: 1e-10,
            strain_limit: 0.05,
        },
    )
    .unwrap();
    let initial = reduced.initial_state();
    let strain = [0.009, -0.002, 0.003];
    let loaded = reduced.evaluate(strain, &initial).unwrap();
    let held = reduced.evaluate(strain, &loaded.state).unwrap();
    assert!(held.state.history.alpha > loaded.state.history.alpha);
    assert!(held.stress_pa[0] < loaded.stress_pa[0]);
    for component in [2, 4, 5] {
        near(held.full_stress_pa[component], 0.0, 1e-10);
    }
    let energy = reduced
        .law()
        .advance_with_energy(
            &loaded.state.strain,
            &loaded.state.history,
            &held.state.strain,
            TOL,
        )
        .unwrap();
    near(energy.energy.mechanical_work_j_m3, 0.0, 1e-11);
    assert!(energy.energy.physical_dissipation_j_m3 > 0.0);
    let perturbation = 1e-7;
    for j in 0..3 {
        let mut plus = strain;
        let mut minus = strain;
        plus[j] += perturbation;
        minus[j] -= perturbation;
        let p = reduced.evaluate(plus, &loaded.state).unwrap();
        let m = reduced.evaluate(minus, &loaded.state).unwrap();
        for i in 0..3 {
            near(
                held.tangent_pa[i][j],
                (p.stress_pa[i] - m.stress_pa[i]) / (2.0 * perturbation),
                3e-6,
            );
        }
    }
    let de: Voigt = core::array::from_fn(|i| held.state.strain[i] - loaded.state.strain[i]);
    near(contract(&held.full_stress_pa, &de), 0.0, 1e-11);
}

#[test]
fn stress_scaling_preserves_history_and_invalid_physical_intervals_refuse() {
    let reference = LinearOverstressHill48::new(native(coefficients(), Q), ETA, DT).unwrap();
    let strain = [0.008, -0.003, -0.005, 0.002, 0.0, 0.001];
    let state = reference.initial_state();
    let expected = reference.response(&strain, &state).unwrap();
    for scale in [1e-180, 1e180] {
        let backbone = Hill48::new(
            IsotropicElastic::new(E * scale, NU, 0.05).unwrap(),
            Y * scale,
            H * scale,
            coefficients(),
        )
        .unwrap()
        .with_linear_kinematic_hardening(C * scale)
        .unwrap()
        .with_voce_isotropic_hardening(Q * scale, B)
        .unwrap();
        let law = LinearOverstressHill48::new(backbone, ETA * scale, DT).unwrap();
        let response = law.response(&strain, &state).unwrap();
        near(response.state.alpha, expected.state.alpha, 3e-15);
        for i in 0..6 {
            near(response.stress_pa[i] / scale, expected.stress_pa[i], 3e-11);
            for j in 0..6 {
                near(
                    response.tangent_pa[i][j] / scale,
                    expected.tangent_pa[i][j],
                    3e-8,
                );
            }
        }
    }
    for (viscosity, duration) in [
        (0.0, DT),
        (-ETA, DT),
        (ETA, 0.0),
        (ETA, -DT),
        (f64::NAN, DT),
        (ETA, f64::INFINITY),
        (f64::MAX, f64::MIN_POSITIVE),
    ] {
        assert!(
            LinearOverstressHill48::new(native(coefficients(), Q), viscosity, duration).is_err()
        );
    }
    assert!(reference.endpoint_overstress_pa(-1.0).is_err());
    assert!(reference.endpoint_overstress_pa(f64::NAN).is_err());
    let mut invalid = state;
    invalid.plastic_strain = [0.001; 6];
    assert!(reference.validate_history(&strain, &invalid).is_err());
}

#[path = "plastic_hill_rate/power.rs"]
mod power;
