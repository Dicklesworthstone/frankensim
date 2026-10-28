//! Independent Hertz/MDR checks for explicit varying-load shear histories.
use fs_math::det;
use fs_tribo::mindlin::loading::{HertzMindlinContact, HertzMindlinParameters, HertzMindlinState};
use fs_tribo::mindlin::{MindlinDeresiewiczPoint, MindlinPointParameters};
use fs_tribo::{TangentialLaw, TangentialSolveConfig, TangentialStepInput};

fn model(cap: usize) -> HertzMindlinContact {
    HertzMindlinContact::new(HertzMindlinParameters {
        effective_radius_m: 0.02,
        shear_moduli_pa: [30e9, 20e9],
        poisson_ratios: [0.25, 0.3],
        friction_coefficient: 0.6,
        maximum_history_intervals: cap,
    })
    .unwrap()
}

fn input(normal: f64, du: f64, compliance: f64) -> TangentialStepInput {
    TangentialStepInput {
        normal_n: normal,
        slip_increment_m: [du, 0.0],
        dt_s: 1e-3,
        contact_temperature_k: Some(300.0),
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

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
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

#[test]
fn hertz_normal_work_and_fixed_load_shear_match_independent_elastic_relations() {
    let law = model(32);
    let p = law.parameters();
    let e = 1.0
        / ((1.0 - p.poisson_ratios[0]) / (2.0 * p.shear_moduli_pa[0])
            + (1.0 - p.poisson_ratios[1]) / (2.0 * p.shear_moduli_pa[1]));
    let normal = 12.0;
    let loaded = law.normal_step(&law.initial_state(), normal).unwrap();
    let a = loaded.state.contact_radius_m();
    let indentation = loaded.state.indentation_m();
    close(
        normal,
        4.0 * e * a.powi(3) / (3.0 * p.effective_radius_m),
        2e-12,
    );
    close(indentation, a * a / p.effective_radius_m, 1e-20);
    close(loaded.normal_work_j, 0.4 * normal * indentation, 1e-19);
    assert_eq!(loaded.dissipated_j, 0.0);
    assert_eq!(law.stored_energy_j(&loaded.state), 0.0);
    let original = fixed(&law, &loaded.state);
    let critical = original.gross_slip_displacement_m();
    let mut evolving = loaded.state;
    let mut constant = original.initial_state();
    for ratio in [0.2, 0.6, -0.3, -0.7, 0.1, 0.4, 1.2] {
        let x = input(
            normal,
            ratio * critical,
            0.25 / original.initial_stiffness_n_per_m(),
        );
        let a = law.advance(&evolving, &x).unwrap();
        let b = original.advance(&constant, &x).unwrap();
        assert_eq!(a.traction_n, b.traction_n);
        assert_eq!(a.dissipated_j, b.dissipated_j);
        assert_eq!(a.state.fixed_load_state(), Some(&b.state));
        evolving = a.state;
        constant = b.state;
    }
}

#[test]
fn normal_growth_adds_unstrained_material_without_inventing_shear_energy() {
    let law = model(32);
    let unloaded = law.normal_step(&law.initial_state(), 10.0).unwrap().state;
    let native = fixed(&law, &unloaded);
    let loaded = law
        .advance(
            &unloaded,
            &input(10.0, 0.6 * native.gross_slip_displacement_m(), 0.0),
        )
        .unwrap()
        .state;
    let grown = law.normal_step(&loaded, 80.0).unwrap();
    close(
        grown.state.contact_radius_m(),
        2.0 * loaded.contact_radius_m(),
        1e-17,
    );
    close(
        grown.state.endpoint_force_n(),
        loaded.endpoint_force_n(),
        1e-12,
    );
    assert_eq!(
        law.stored_energy_j(&grown.state),
        law.stored_energy_j(&loaded)
    );
    assert_eq!(grown.dissipated_j, 0.0);
    assert_eq!(grown.tangential_stored_change_j, 0.0);
    // Newly contacting material is initially relaxed, then joins the true
    // enlarged 8aG* stiffness under an additional small shear increment.
    let enlarged = fixed(&law, &grown.state);
    let du = 1e-7 * enlarged.gross_slip_displacement_m();
    let shear = law.advance(&grown.state, &input(80.0, du, 0.0)).unwrap();
    close(
        (shear.state.endpoint_force_n() - grown.state.endpoint_force_n()) / du,
        enlarged.initial_stiffness_n_per_m(),
        2e-6 * enlarged.initial_stiffness_n_per_m(),
    );
}

#[test]
fn unloading_a_saturated_patch_exports_the_analytic_a_fifth_energy_release() {
    let law = model(32);
    let normal = law.normal_step(&law.initial_state(), 16.0).unwrap().state;
    let native = fixed(&law, &normal);
    let saturated = law
        .advance(
            &normal,
            &input(16.0, 2.0 * native.gross_slip_displacement_m(), 0.0),
        )
        .unwrap()
        .state;
    let old_energy = law.stored_energy_j(&saturated);
    let p = law.parameters();
    let g_mdr = 4.0
        / ((2.0 - p.poisson_ratios[0]) / p.shear_moduli_pa[0]
            + (2.0 - p.poisson_ratios[1]) / p.shear_moduli_pa[1]);
    let beta = p.friction_coefficient * law.normal_modulus_pa() / (g_mdr * p.effective_radius_m);
    let independent = 8.0 / 15.0 * g_mdr * beta.powi(2) * saturated.contact_radius_m().powi(5);
    close(old_energy, independent, 2e-14 * independent);
    let reduced = law.normal_step(&saturated, 2.0).unwrap();
    close(
        law.stored_energy_j(&reduced.state),
        old_energy / 32.0,
        3e-14 * old_energy,
    );
    close(
        reduced.dissipated_j,
        31.0 / 32.0 * old_energy,
        3e-14 * old_energy,
    );
    close(reduced.state.endpoint_force_n(), 0.6 * 2.0, 1e-13);
    close(
        reduced.tangential_stored_change_j + reduced.dissipated_j,
        0.0,
        1e-20,
    );
    assert!(reduced.normal_work_j < 0.0);
    // A very small decrement still has positive, resolved release rather than
    // a cancellation-derived zero or negative heat.
    let tiny = law.normal_step(&saturated, 16.0 * (1.0 - 1e-8)).unwrap();
    assert!(tiny.dissipated_j > 0.0);
    close(tiny.dissipated_j / old_energy, 5.0 / 3.0 * 1e-8, 2e-15);
}

#[test]
fn separation_pays_for_old_memory_and_recontact_is_stress_free() {
    let law = model(32);
    let normal = law.normal_step(&law.initial_state(), 10.0).unwrap().state;
    let native = fixed(&law, &normal);
    let loaded = law
        .advance(
            &normal,
            &input(10.0, 0.6 * native.gross_slip_displacement_m(), 0.0),
        )
        .unwrap()
        .state;
    let released = law.normal_step(&loaded, 0.0).unwrap();
    assert_eq!(released.dissipated_j, law.stored_energy_j(&loaded));
    assert_eq!(released.normal_work_j, -loaded.normal_energy_j());
    assert_eq!(released.state.endpoint_force_n(), 0.0);
    assert_eq!(released.state.history_intervals(), 0);
    let free = law
        .advance(&released.state, &input(0.0, 1e-3, 1.0))
        .unwrap();
    assert_eq!(free.traction_n, [0.0; 2]);
    assert_eq!(free.dissipated_j, 0.0);
    let recontact = law.normal_step(&free.state, 10.0).unwrap().state;
    assert_eq!(recontact.endpoint_force_n(), 0.0);
    assert_eq!(law.stored_energy_j(&recontact), 0.0);
    let du = 0.2 * native.gross_slip_displacement_m();
    let actual = law.advance(&recontact, &input(10.0, du, 0.0)).unwrap();
    let virgin = native
        .advance(&native.initial_state(), &input(10.0, du, 0.0))
        .unwrap();
    assert_eq!(actual.traction_n, virgin.traction_n);
    assert_eq!(actual.dissipated_j, virgin.dissipated_j);
    close(actual.state.slip_m(), loaded.slip_m() + 1e-3 + du, 1e-18);
}

#[test]
fn normal_subdivision_preserves_memory_and_load_shear_order_is_physical() {
    let law = model(64);
    let zero = law.normal_step(&law.initial_state(), 10.0).unwrap().state;
    let native = fixed(&law, &zero);
    let du = 0.6 * native.gross_slip_displacement_m();
    let loaded = law.advance(&zero, &input(10.0, du, 0.0)).unwrap().state;
    let direct = law.normal_step(&loaded, 2.0).unwrap();
    let mut subdivided = loaded.clone();
    let mut heat = 0.0;
    for k in 1..=32 {
        let step = law
            .normal_step(&subdivided, 10.0 - 8.0 * k as f64 / 32.0)
            .unwrap();
        heat += step.dissipated_j;
        subdivided = step.state;
    }
    close(
        subdivided.endpoint_force_n(),
        direct.state.endpoint_force_n(),
        2e-13,
    );
    close(
        law.stored_energy_j(&subdivided),
        law.stored_energy_j(&direct.state),
        2e-13 * law.stored_energy_j(&loaded),
    );
    close(
        heat,
        direct.dissipated_j,
        2e-13 * law.stored_energy_j(&loaded),
    );
    let normal_first = law.normal_step(&zero, 80.0).unwrap().state;
    let normal_first = law
        .advance(&normal_first, &input(80.0, du, 0.0))
        .unwrap()
        .state;
    let shear_first = law.normal_step(&loaded, 80.0).unwrap().state;
    assert!(normal_first.endpoint_force_n() > 1.5 * shear_first.endpoint_force_n());
    let x = input(2.0, -du, 0.3 / native.initial_stiffness_n_per_m());
    let a = law.advance(&direct.state, &x).unwrap();
    let b = law.advance(&subdivided, &x).unwrap();
    close(a.traction_n[0], b.traction_n[0], 1e-12);
    close(
        a.dissipated_j,
        b.dissipated_j,
        3e-13 * law.stored_energy_j(&loaded),
    );
}

#[test]
fn invalid_paths_memory_budget_and_foreign_checkpoints_refuse_atomically() {
    let law = model(2);
    let normal = law.normal_step(&law.initial_state(), 10.0).unwrap().state;
    let native = fixed(&law, &normal);
    let loaded = law
        .advance(
            &normal,
            &input(10.0, 0.4 * native.gross_slip_displacement_m(), 0.0),
        )
        .unwrap()
        .state;
    let saved = loaded.clone();
    assert!(law.normal_step(&loaded, 80.0).is_err());
    for value in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(law.normal_step(&loaded, value).is_err());
    }
    assert!(law.advance(&loaded, &input(11.0, 0.0, 0.0)).is_err());
    let mut transverse = input(10.0, 0.0, 0.0);
    transverse.slip_increment_m[1] = 1e-6;
    assert!(law.advance(&loaded, &transverse).is_err());
    let mut other = law.parameters();
    other.effective_radius_m *= 2.0;
    let other = HertzMindlinContact::new(other).unwrap();
    assert!(other.normal_step(&loaded, 10.0).is_err());
    assert!(other.stored_energy_j(&loaded).is_nan());
    assert_eq!(loaded, saved);
    let x = input(10.0, 2.0 * native.gross_slip_displacement_m(), 0.0);
    assert_eq!(
        law.advance(&loaded, &x).unwrap(),
        law.advance(&saved, &x).unwrap()
    );
    let mut frictionless = law.parameters();
    frictionless.friction_coefficient = 0.0;
    let frictionless = HertzMindlinContact::new(frictionless).unwrap();
    let first = frictionless
        .normal_step(&frictionless.initial_state(), 12.0)
        .unwrap()
        .state;
    let slipped = frictionless
        .advance(&first, &input(12.0, 1e-3, 1.0))
        .unwrap()
        .state;
    let changed = frictionless.normal_step(&slipped, 2.0).unwrap();
    assert_eq!(changed.dissipated_j, 0.0);
    assert_eq!(changed.state.endpoint_force_n(), 0.0);
}

// Independent midpoint quadrature of physical MDR springs in a fixed x grid.
// It has no analytical history intervals or normalized-r remapping.
fn grid_path(cells: usize, law: &HertzMindlinContact, scale: f64) -> (f64, f64, f64) {
    let p = law.parameters();
    let e = law.normal_modulus_pa();
    let g = 4.0
        / ((2.0 - p.poisson_ratios[0]) / p.shear_moduli_pa[0]
            + (2.0 - p.poisson_ratios[1]) / p.shear_moduli_pa[1]);
    let amax = det::pow(0.75 * 30.0 * p.effective_radius_m / e, 1.0 / 3.0);
    let dx = amax / cells as f64;
    let mut strain = vec![0.0_f64; cells];
    let mut heat = 0.0;
    let mut radius = 0.0;
    for (load, shift) in path(scale) {
        if shift == 0.0 {
            radius = det::pow(0.75 * load * p.effective_radius_m / e, 1.0 / 3.0);
            for (i, s) in strain.iter_mut().enumerate() {
                let x = (i as f64 + 0.5) * dx;
                let bound = (p.friction_coefficient * e / (g * p.effective_radius_m)
                    * (radius * radius - x * x))
                    .max(0.0);
                let next = (*s).clamp(-bound, bound);
                heat += g * dx * (*s - next) * (*s + next);
                *s = next;
            }
        } else {
            for (i, s) in strain.iter_mut().enumerate() {
                let x = (i as f64 + 0.5) * dx;
                if x >= radius {
                    continue;
                }
                let bound = p.friction_coefficient * e / (g * p.effective_radius_m)
                    * (radius * radius - x * x);
                let next = (*s + shift).clamp(-bound, bound);
                heat += 2.0 * g * dx * bound * (shift - (next - *s)).abs();
                *s = next;
            }
        }
    }
    let force = 2.0 * g * dx * strain.iter().sum::<f64>();
    let energy = g * dx * strain.iter().map(|s| s * s).sum::<f64>();
    (force, energy, heat)
}

fn path(scale: f64) -> [(f64, f64); 12] {
    [
        (10.0, 0.0),
        (10.0, 0.7 * scale),
        (30.0, 0.0),
        (30.0, -0.2 * scale),
        (5.0, 0.0),
        (5.0, -0.9 * scale),
        (18.0, 0.0),
        (18.0, 0.5 * scale),
        (2.0, 0.0),
        (2.0, 0.6 * scale),
        (12.0, 0.0),
        (12.0, -0.3 * scale),
    ]
}

#[test]
fn mixed_normal_and_reversing_shear_history_matches_independent_spatial_refinement() {
    let law = model(64);
    let normal = law.normal_step(&law.initial_state(), 10.0).unwrap().state;
    let scale = fixed(&law, &normal).gross_slip_displacement_m();
    let mut state = law.initial_state();
    let mut heat = 0.0;
    for (load, shift) in path(scale) {
        if shift == 0.0 {
            let step = law.normal_step(&state, load).unwrap();
            heat += step.dissipated_j;
            state = step.state;
        } else {
            let step = law.advance(&state, &input(load, shift, 0.0)).unwrap();
            heat += step.dissipated_j;
            state = step.state;
        }
    }
    let exact = (state.endpoint_force_n(), law.stored_energy_j(&state), heat);
    let mut errors = Vec::new();
    for count in [512, 1024, 2048] {
        let numeric = grid_path(count, &law, scale);
        let error = (numeric.0 - exact.0).abs() / 10.0
            + (numeric.1 - exact.1).abs() / (10.0 * scale)
            + (numeric.2 - exact.2).abs() / (10.0 * scale);
        errors.push(error);
        assert!(
            error < 5.0 / (count * count) as f64,
            "n={count}, error={error:e}"
        );
    }
    assert!(errors[2] < errors[0] / 4.0, "refinement errors={errors:?}");
}
