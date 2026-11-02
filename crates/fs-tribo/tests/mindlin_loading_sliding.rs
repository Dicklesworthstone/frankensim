//! Integrated load-distance for the original simultaneous Hertz--Mindlin path.
//! Synthetic elastic/friction parameters verify equations, not wear calibration.
use fs_tribo::mindlin::loading::path::MindlinPathError;
use fs_tribo::mindlin::loading::{HertzMindlinContact, HertzMindlinParameters, HertzMindlinState};
use fs_tribo::mindlin::{MindlinDeresiewiczPoint, MindlinPointParameters};
use fs_tribo::{TangentialLaw, TangentialSolveConfig, TangentialStepInput};

fn model(mu: f64) -> HertzMindlinContact {
    HertzMindlinContact::new(HertzMindlinParameters {
        effective_radius_m: 0.02,
        shear_moduli_pa: [30e9, 20e9],
        poisson_ratios: [0.25, 0.3],
        friction_coefficient: mu,
        maximum_history_intervals: 128,
    })
    .unwrap()
}
fn input(normal_n: f64, du: f64, compliance: f64) -> TangentialStepInput {
    TangentialStepInput {
        normal_n,
        slip_increment_m: [du, 0.0],
        dt_s: 1e-3,
        contact_temperature_k: None,
        nominal_pressure_pa: None,
        compliance_m_per_n: compliance,
        solve: TangentialSolveConfig {
            max_iterations: 64,
            maximum_force_n: 100.0,
            force_absolute_tolerance_n: 1e-13,
            force_relative_tolerance: 1e-13,
        },
    }
}
fn fixed(law: &HertzMindlinContact, state: &HertzMindlinState) -> MindlinDeresiewiczPoint {
    let p = law.parameters();
    MindlinDeresiewiczPoint::new(MindlinPointParameters {
        normal_force_n: state.normal_force_n(),
        contact_radius_m: state.contact_radius_m(),
        shear_moduli_pa: p.shear_moduli_pa,
        poisson_ratios: p.poisson_ratios,
        friction_coefficient: p.friction_coefficient,
        maximum_history_intervals: p.maximum_history_intervals,
    })
    .unwrap()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

#[test]
fn frictionless_paths_integrate_the_nonlinear_hertz_load_through_contact_and_liftoff() {
    let law = model(0.0);
    for (n0, n1) in [
        (0.0, 0.0),
        (0.0, 10.0),
        (10.0, 0.0),
        (10.0, 80.0),
        (80.0, 10.0),
        (10.0, 10.0),
        (10.0, 10.0 * (1.0 + 1e-12)),
    ] {
        let before = law.normal_step(&law.initial_state(), n0).unwrap().state;
        for du in [-2e-7, 2e-7] {
            let receipt = law
                .advance_coupled_path_with_sliding(&before, &input(n1, du, 3.0), &mut || false)
                .unwrap();
            let step = receipt.path_step();
            let mean = if n0 == n1 {
                n0
            } else if n0.min(n1) == 0.0 {
                0.4 * n0.max(n1)
            } else if (n1 - n0).abs() < 1e-8 * n0 {
                0.5 * (n0 + n1)
            } else {
                // Independent Hertz energy difference divided by indentation change.
                let d0 = before.indentation_m();
                let d1 = step.tangential.state.indentation_m();
                0.4 * (n1 * d1 - n0 * d0) / (d1 - d0)
            };
            close(
                receipt.pressure_weighted_sliding_j(),
                mean * du.abs(),
                3e-14 * n0.max(n1) * du.abs(),
            );
            assert_eq!(step.tangential.slip_increment_m[0], du);
            assert_eq!(step.tangential.dissipated_j, 0.0);
        }
    }
    // Refinement must follow affine indentation, not affine normal force.
    let before = law.normal_step(&law.initial_state(), 10.0).unwrap().state;
    let direct = law
        .advance_coupled_path_with_sliding(&before, &input(80.0, 2e-7, 0.0), &mut || false)
        .unwrap();
    let d0 = before.indentation_m();
    let d1 = direct.path_step().tangential.state.indentation_m();
    let mut state = before;
    let mut total = 0.0;
    for i in 1..=16 {
        let depth = d0 + (d1 - d0) * i as f64 / 16.0;
        let load = if i == 16 {
            80.0
        } else {
            10.0 * (depth / d0).powf(1.5)
        };
        let part = law
            .advance_coupled_path_with_sliding(&state, &input(load, 2e-7 / 16.0, 0.0), &mut || {
                false
            })
            .unwrap();
        total += part.pressure_weighted_sliding_j();
        state = part.into_path_step().tangential.state;
    }
    close(total, direct.pressure_weighted_sliding_j(), 2e-12 * total);
    assert!((total - 0.5 * (10.0 + 80.0) * 2e-7).abs() > 0.05 * total);
}

#[test]
fn normal_unloading_drives_microslip_at_zero_macroscopic_shear_and_recontact_is_clean() {
    let law = model(0.6);
    let normal = law.normal_step(&law.initial_state(), 10.0).unwrap().state;
    let d = fixed(&law, &normal).gross_slip_displacement_m();
    let loaded = law
        .advance(&normal, &input(10.0, 2.0 * d, 0.0))
        .unwrap()
        .state;
    for load in [1.25_f64, 0.0] {
        let receipt = law
            .advance_coupled_path_with_sliding(&loaded, &input(load, 0.0, 0.0), &mut || false)
            .unwrap();
        // Saturated shear storage scales with a^5, hence N^(5/3).
        let expected = 0.4 * 10.0 * d * (1.0 - (load / 10.0).powf(5.0 / 3.0));
        close(
            receipt.pressure_weighted_sliding_j(),
            expected,
            3e-13 * expected,
        );
        assert_eq!(receipt.path_step().tangential.slip_increment_m[0], 0.0);
        assert!(receipt.pressure_weighted_sliding_j() > 0.0);
        assert!(receipt.path_step().normal_work_j < 0.0);
    }
    let separated = law
        .advance_coupled_path_with_sliding(&loaded, &input(0.0, 0.0, 0.0), &mut || false)
        .unwrap()
        .into_path_step()
        .tangential
        .state;
    assert_eq!(separated.history_intervals(), 0);
    let regrown = law
        .advance_coupled_path_with_sliding(&separated, &input(10.0, 0.4 * d, 0.0), &mut || false)
        .unwrap();
    assert_eq!(regrown.pressure_weighted_sliding_j(), 0.0);
    assert!(law.stored_energy_j(&regrown.path_step().tangential.state) > 0.0);
}

// Independent time/space return map: integrate q_normal*|plastic travel|
// directly, never native heat or a difference of total work/storage. q_normal
// is the conjugate MDR pressure density, E*max(delta-x^2/R,0), on each half-line.
fn grid_measure(
    law: &HertzMindlinContact,
    depth: f64,
    shear: f64,
    path: &[(f64, f64)],
    n: usize,
) -> f64 {
    let p = law.parameters();
    let e = 1.0
        / ((1.0 - p.poisson_ratios[0]) / (2.0 * p.shear_moduli_pa[0])
            + (1.0 - p.poisson_ratios[1]) / (2.0 * p.shear_moduli_pa[1]));
    let g = 4.0
        / ((2.0 - p.poisson_ratios[0]) / p.shear_moduli_pa[0]
            + (2.0 - p.poisson_ratios[1]) / p.shear_moduli_pa[1]);
    let max_depth = path.iter().map(|x| x.0).fold(0.0, f64::max) * depth;
    let dx = (p.effective_radius_m * max_depth).sqrt() / n as f64;
    let mut strains = vec![0.0; n];
    let mut measure = 0.0;
    let steps = 2 * n;
    for segment in path.windows(2) {
        for time in 0..steps {
            let d0 =
                depth * (segment[0].0 + (segment[1].0 - segment[0].0) * time as f64 / steps as f64);
            let d1 = depth
                * (segment[0].0 + (segment[1].0 - segment[0].0) * (time + 1) as f64 / steps as f64);
            let du = (segment[1].1 - segment[0].1) * shear / steps as f64;
            for (index, strain) in strains.iter_mut().enumerate() {
                let x = (index as f64 + 0.5) * dx;
                let profile = x * x / p.effective_radius_m;
                if profile >= d0.max(d1) {
                    *strain = 0.0;
                    continue;
                }
                let enter = if d1 > d0 {
                    ((profile - d0) / (d1 - d0)).clamp(0.0, 1.0)
                } else {
                    0.0
                };
                let exit = if d1 < d0 {
                    ((d0 - profile) / (d0 - d1)).clamp(0.0, 1.0)
                } else {
                    1.0
                };
                let q0 = e * (d0 + (d1 - d0) * enter - profile).max(0.0);
                let q1 = e * (d0 + (d1 - d0) * exit - profile).max(0.0);
                let old = if enter > 0.0 { 0.0 } else { *strain };
                let travel = du * (exit - enter);
                let bound = p.friction_coefficient * q1 / g;
                let next = (old + travel).clamp(-bound, bound);
                measure += 2.0 * dx * 0.5 * (q0 + q1) * (travel - (next - old)).abs();
                *strain = if exit < 1.0 { 0.0 } else { next };
            }
        }
    }
    measure
}

#[test]
fn independent_pressure_and_plastic_slip_history_refines_to_the_native_path_measure() {
    let trajectory: [(f64, f64); 7] = [
        (0.0, 0.0),
        (1.0, 0.7),
        (0.2, 1.3),
        (1.4, -0.8),
        (0.0, -0.8),
        (0.8, -0.3),
        (0.0, -0.3),
    ];
    for mu in [0.0, 0.6] {
        let law = model(mu);
        let reference = law.normal_step(&law.initial_state(), 10.0).unwrap().state;
        let depth = reference.indentation_m();
        let shear = 2e-7;
        let mut state = law.initial_state();
        let mut measure = 0.0;
        for segment in trajectory.windows(2) {
            let load = 10.0_f64 * segment[1].0.powf(1.5);
            let receipt = law
                .advance_coupled_path_with_sliding(
                    &state,
                    &input(load, (segment[1].1 - segment[0].1) * shear, 0.0),
                    &mut || false,
                )
                .unwrap();
            measure += receipt.pressure_weighted_sliding_j();
            state = receipt.into_path_step().tangential.state;
        }
        let errors = [64, 128, 256]
            .map(|n| (grid_measure(&law, depth, shear, &trajectory, n) - measure).abs() / measure);
        // Entry/yield locations need not align with successive midpoint grids;
        // require decreasing errors and a material reduction over the ladder.
        assert!(
            errors[1] < errors[0] && errors[2] < errors[1] && errors[2] < 0.2 * errors[0],
            "mu={mu}: {errors:?}"
        );
        assert!(errors[2] < 1e-3, "mu={mu}: {errors:?}");
        assert_eq!(state.normal_force_n(), 0.0);
    }
}

#[test]
fn exact_fixed_load_and_compliant_paths_keep_original_budgets_and_atomicity() {
    let law = model(0.6);
    let normal = law.normal_step(&law.initial_state(), 10.0).unwrap().state;
    let fixed = fixed(&law, &normal);
    let d = fixed.gross_slip_displacement_m();
    let before = law
        .advance(&normal, &input(10.0, 0.6 * d, 0.0))
        .unwrap()
        .state;
    let x = input(10.0, -0.2 * d, 0.5 / fixed.initial_stiffness_n_per_m());
    let original = fixed
        .advance_with_surface_sliding(before.fixed_load_state().unwrap(), &x)
        .unwrap();
    let wrapped = law
        .advance_coupled_path_with_sliding(&before, &x, &mut || false)
        .unwrap();
    assert_eq!(
        wrapped.pressure_weighted_sliding_j(),
        original.pressure_weighted_sliding_j()
    );
    let x = input(5.0, -0.2 * d, 0.5 / fixed.initial_stiffness_n_per_m());
    let mut original_polls = 0;
    let native = law
        .advance_coupled_path(&before, &x, &mut || {
            original_polls += 1;
            false
        })
        .unwrap();
    let mut wrapped_polls = 0;
    let candidate = law
        .advance_coupled_path_with_sliding(&before, &x, &mut || {
            wrapped_polls += 1;
            false
        })
        .unwrap();
    assert_eq!(candidate.path_step(), &native);
    assert_eq!(original_polls, wrapped_polls);
    assert_ne!(
        candidate.path_step().tangential.slip_increment_m[0],
        x.slip_increment_m[0]
    );
    let saved = before.clone();
    assert!(matches!(
        law.advance_coupled_path_with_sliding(&before, &x, &mut || true),
        Err(MindlinPathError::Cancelled)
    ));
    assert!(
        model(0.5)
            .advance_coupled_path_with_sliding(&before, &x, &mut || false)
            .is_err()
    );
    let mut invalid = x;
    invalid.solve.maximum_force_n = 0.001;
    assert!(
        law.advance_coupled_path_with_sliding(&before, &invalid, &mut || false)
            .is_err()
    );
    assert_eq!(before, saved);
}
