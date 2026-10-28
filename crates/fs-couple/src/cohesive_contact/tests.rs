use super::*;
use fs_material::cohesive::{BilinearCohesive, CohesiveViscosity, mixed_mode::MixedModeParameters};

fn law() -> BilinearCohesive {
    BilinearCohesive::new_mixed_mode(
        MixedModeParameters {
            normal_stiffness_pa_m: 1e9,
            tangential_stiffness_pa_m: 4e8,
            normal_strength_pa: 2e5,
            tangential_strength_pa: 1.6e5,
            mode_i_fracture_energy_j_m2: 100.0,
            mode_ii_fracture_energy_j_m2: 400.0,
            bk_exponent: 2.0,
        },
        2e9,
        0.02,
    )
    .unwrap()
}

fn config() -> CohesiveContactConfig {
    CohesiveContactConfig {
        max_iterations: 24,
        max_material_evaluations: 400,
        max_backtracks: 12,
        jump_absolute_tolerance_m: 1e-14,
        jump_relative_tolerance: 1e-11,
        maximum_jump_increment_m: 0.02,
        maximum_force_n: 100.0,
        material_energy_tolerance_j_m2: 1e-8,
    }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}
fn mv(a: [[f64; 3]; 3], x: [f64; 3]) -> [f64; 3] {
    a.map(|row| dot(row, x))
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}
fn close_vector(a: [f64; 3], b: [f64; 3], tolerance: f64) {
    for i in 0..3 {
        close(a[i], b[i], tolerance);
    }
}

#[test]
fn full_off_diagonal_compliance_recovers_a_manufactured_elastic_solution() {
    let law = law();
    let old = [5e-5, -4e-5, 3e-5];
    let expected = [1e-4, 8e-5, -6e-5];
    let point = law.bonded_point().trial(old, 1e-8).unwrap().point;
    let area = 2e-5;
    let stiffness = [1e9, 4e8, 4e8];
    let force = core::array::from_fn(|i| area * stiffness[i] * 0.5 * (old[i] + expected[i]));
    // Both symmetric mechanical compliance and the explicitly retained
    // nonsymmetric caller map are tested; no diagonal surrogate reproduces g.
    for compliance in [
        [
            [2e-4, 4e-5, -3e-5],
            [4e-5, 3e-4, 5e-5],
            [-3e-5, 5e-5, 1.5e-4],
        ],
        [
            [2e-4, 7e-5, -3e-5],
            [4e-5, 3e-4, 5e-5],
            [-2e-5, 5e-5, 1.5e-4],
        ],
    ] {
        let response = mv(compliance, force);
        let free = core::array::from_fn(|i| expected[i] + response[i]);
        let mut owner = CohesiveLoadLine::new(area, compliance, config()).unwrap();
        let step = owner.trial(point, free, 0.01, || true).unwrap();
        close_vector(step.material.equilibrium.point.jump_m(), expected, 2e-14);
        close_vector(step.mean_force_n, force, 3e-10);
        assert_eq!(step.material.equilibrium.point.damage(), 0.0);
        assert_eq!(owner.compliance_m_per_n(), compliance);
        assert!(norm_inf(step.compatibility_residual_m) <= step.compatibility_tolerance_m);
        close(
            dot(
                step.mean_force_n,
                core::array::from_fn(|i| expected[i] - old[i]),
            ),
            area * step.material.total_work_j_m2(),
            2e-12,
        );
    }
}

#[test]
fn free_two_body_three_axis_dynamics_conserve_momentum_and_total_fracture_viscous_energy() {
    for law in [
        law(),
        law().with_viscosity(CohesiveViscosity::new(2e6, 7e5, 2.0).unwrap()),
    ] {
        let (ma, mb, area, dt) = (2.0, 3.0, 1e-6, 2e-4);
        let compliance = 0.5 * dt * dt * (1.0 / ma + 1.0 / mb);
        let mut owner = CohesiveLoadLine::new(
            area,
            core::array::from_fn(|i| {
                core::array::from_fn(|j| if i == j { compliance } else { 0.0 })
            }),
            config(),
        )
        .unwrap();
        let center = [0.012, -0.006, 0.003];
        let relative = [0.05, 0.04, -0.025];
        let mut va = core::array::from_fn::<_, 3, _>(|i| center[i] - mb / (ma + mb) * relative[i]);
        let mut vb = core::array::from_fn::<_, 3, _>(|i| center[i] + ma / (ma + mb) * relative[i]);
        let momentum = core::array::from_fn(|i| ma * va[i] + mb * vb[i]);
        let initial = 0.5 * ma * dot(va, va) + 0.5 * mb * dot(vb, vb);
        let (mut xa, mut xb) = ([0.0; 3], [0.0; 3]);
        let mut point = law.bonded_point();
        let (mut fracture, mut viscous) = (0.0, 0.0);
        let mut partial = false;
        for _ in 0..500 {
            let free = core::array::from_fn(|i| xb[i] - xa[i] + dt * (vb[i] - va[i]));
            let step = owner.trial(point, free, dt, || true).unwrap();
            let force = step.mean_force_n;
            for i in 0..3 {
                xa[i] += dt * va[i] + 0.5 * dt * dt * force[i] / ma;
                xb[i] += dt * vb[i] - 0.5 * dt * dt * force[i] / mb;
                va[i] += dt * force[i] / ma;
                vb[i] -= dt * force[i] / mb;
            }
            fracture += area * step.material.equilibrium.energy.dissipated;
            viscous += area * step.material.viscous_dissipated_j_m2;
            point = step.material.equilibrium.point;
            partial |= point.damage() > 0.0 && point.damage() < 1.0;
            close_vector(
                point.jump_m(),
                core::array::from_fn(|i| xb[i] - xa[i]),
                2e-13,
            );
            close_vector(
                core::array::from_fn(|i| ma * va[i] + mb * vb[i]),
                momentum,
                3e-15,
            );
            let kinetic = 0.5 * ma * dot(va, va) + 0.5 * mb * dot(vb, vb);
            close(
                kinetic + area * point.stored_energy_j_m2() + fracture + viscous,
                initial,
                3e-12,
            );
        }
        assert!(partial && point.failed());
        close(fracture, area * point.dissipated_j_m2(), 1e-14);
        assert!(fracture > 0.0);
        assert_eq!(viscous > 0.0, law.viscosity().is_some());
        assert_eq!(point.traction_pa(), [0.0; 3]);
    }
}

#[test]
fn opening_shear_unloading_and_reclosure_preserve_irreversible_mixed_mode_failure() {
    let law = law();
    let mut owner = CohesiveLoadLine::new(1e-6, [[0.0; 3]; 3], config()).unwrap();
    let initial = owner
        .trial(law.bonded_point(), [4e-4, 3e-4, -1e-4], 0.02, || true)
        .unwrap();
    let mut point = initial.material.equilibrium.point;
    let history = point.mixed_mode_history();
    assert!(point.damage() > 0.0 && !point.failed());
    let retained = point.damage();
    let unloaded = owner
        .trial(point, [5e-5, -4e-5, 2e-5], 0.02, || true)
        .unwrap();
    point = unloaded.material.equilibrium.point;
    close(point.damage(), retained, 1e-14);
    assert_eq!(unloaded.material.equilibrium.energy.dissipated, 0.0);
    let broken = owner
        .trial(point, [1e-3, 0.012, -0.003], 0.2, || true)
        .unwrap();
    point = broken.material.equilibrium.point;
    assert!(point.failed());
    assert_eq!(point.mixed_mode_history(), history);
    assert!(norm_inf(broken.mean_force_n) > 0.0);
    assert_eq!(broken.material.traction_pa, [0.0; 3]);
    let closed = owner
        .trial(point, [-5e-4, 0.006, 0.0], 0.2, || true)
        .unwrap();
    assert!(closed.material.equilibrium.point.failed());
    assert_eq!(closed.material.equilibrium.energy.dissipated, 0.0);
    close_vector(closed.mean_force_n, [-1.0 / 6.0, 0.0, 0.0], 2e-12);
    close_vector(closed.material.traction_pa, [-1e6, 0.0, 0.0], 1e-7);
}

#[test]
fn rotating_both_shear_frame_and_full_compliance_rotates_the_same_physical_solution() {
    let law = law();
    let rotation = [[1.0, 0.0, 0.0], [0.0, 0.6, -0.8], [0.0, 0.8, 0.6]];
    let compliance = [
        [2e-6, 4e-7, -3e-7],
        [4e-7, 3e-6, 5e-7],
        [-3e-7, 5e-7, 1.5e-6],
    ];
    let rotated = core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            (0..3)
                .flat_map(|k| {
                    (0..3).map(move |l| rotation[i][k] * compliance[k][l] * rotation[j][l])
                })
                .sum()
        })
    });
    let start = [5e-5, -4e-5, 3e-5];
    let free = [5e-4, 1e-3, -3e-4];
    let a = law.bonded_point().trial(start, 1e-8).unwrap().point;
    let b = law
        .bonded_point()
        .trial(mv(rotation, start), 1e-8)
        .unwrap()
        .point;
    let x = CohesiveLoadLine::new(2e-5, compliance, config())
        .unwrap()
        .trial(a, free, 0.02, || true)
        .unwrap();
    let y = CohesiveLoadLine::new(2e-5, rotated, config())
        .unwrap()
        .trial(b, mv(rotation, free), 0.02, || true)
        .unwrap();
    close_vector(
        mv(rotation, x.material.equilibrium.point.jump_m()),
        y.material.equilibrium.point.jump_m(),
        4e-13,
    );
    close_vector(mv(rotation, x.mean_force_n), y.mean_force_n, 2e-9);
    close(
        x.material.total_work_j_m2(),
        y.material.total_work_j_m2(),
        2e-8,
    );
    close(
        x.material.equilibrium.point.damage(),
        y.material.equilibrium.point.damage(),
        2e-10,
    );
}

#[test]
fn cancellation_force_range_and_evaluation_budget_leave_original_history_unchanged() {
    let law = law();
    let point = law
        .bonded_point()
        .trial([3e-4, 2e-4, 1e-4], 1e-8)
        .unwrap()
        .point;
    let saved = point;
    let free = [6e-4, 7e-4, -2e-4];
    let mut owner = CohesiveLoadLine::new(1e-6, [[0.0; 3]; 3], config()).unwrap();
    assert_eq!(
        owner.trial(point, free, 0.01, || false),
        Err(CohesiveContactError::Cancelled)
    );
    let mut polls = 0;
    assert_eq!(
        owner.trial(point, free, 0.01, || {
            polls += 1;
            polls < 3
        }),
        Err(CohesiveContactError::Cancelled)
    );
    let mut short = config();
    short.material_energy_tolerance_j_m2 = 0.0;
    assert!(CohesiveLoadLine::new(1e-6, [[0.0; 3]; 3], short).is_err());
    short = config();
    short.max_material_evaluations = 1;
    let mut limited = CohesiveLoadLine::new(1e-6, [[0.0; 3]; 3], short).unwrap();
    assert!(matches!(
        limited.trial(point, free, 0.01, || true),
        Err(CohesiveContactError::Budget(_))
    ));
    short = config();
    short.maximum_force_n = 1e-12;
    let mut limited = CohesiveLoadLine::new(1e-6, [[0.0; 3]; 3], short).unwrap();
    assert!(matches!(
        limited.trial(point, free, 0.01, || true),
        Err(CohesiveContactError::Range(_))
    ));
    assert_eq!(point, saved);
    assert_eq!(
        owner.trial(point, free, 0.01, || true),
        owner.trial(saved, free, 0.01, || true)
    );
}

#[test]
fn a_valid_inward_step_at_the_material_range_boundary_uses_one_sided_derivatives() {
    let law = law();
    let boundary = [0.012, 0.0096, 0.0128];
    let inward = boundary.map(|g| 0.95 * g);
    let point = law.bonded_point().trial(boundary, 1e-8).unwrap().point;
    let mut owner = CohesiveLoadLine::new(1e-6, [[0.0; 3]; 3], config()).unwrap();
    let step = owner.trial(point, inward, 0.01, || true).unwrap();
    close_vector(step.material.equilibrium.point.jump_m(), inward, 2e-14);
    assert_eq!(step.mean_force_n, [0.0; 3]);
}
