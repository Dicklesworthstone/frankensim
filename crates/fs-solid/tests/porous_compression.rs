//! An independent hydrostatic reference through the native force-controlled
//! tetrahedral solve. Authored coefficients exercise pressure-sensitive plastic
//! volume change; they do not calibrate a particular foam or finite crushing.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::IsotropicElastic;
use fs_material::plastic::porous::PorousPlasticity;
use fs_solid::{IncrementSettings, SmallStrainLawField, SmallStrainTetProblem, TetAssemblyBudget};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 5,
                kernel_id: 83,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}",
    );
}

#[test]
fn hydrostatic_force_path_retains_plastic_volume_after_unloading_and_reloading() {
    const YOUNG: f64 = 1.0e9;
    const POISSON: f64 = 0.22;
    const YIELD: f64 = 1.0e6;
    const HARDENING: f64 = 1.0e8;
    const SENSITIVITY: f64 = 1.5;
    const PEAK_PRESSURE: f64 = 2.0e6;

    let law = PorousPlasticity::new(
        IsotropicElastic::new(YOUNG, POISSON, 0.05).unwrap(),
        YIELD,
        HARDENING,
        SENSITIVITY,
    )
    .unwrap();
    let nodes = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let tetrahedra = [[0, 1, 2, 3]];
    // Remove only the six rigid-body modes. The normal displacement on each
    // coordinate axis remains free: the solver must find the volume change.
    let prescribed = [(0, 0.0), (1, 0.0), (2, 0.0), (4, 0.0), (5, 0.0), (8, 0.0)];
    // V=1/6 and grad(N_1,N_2,N_3) are the three coordinate unit vectors.
    // These are the balanced nodal forces of uniform stress -p I, including
    // the opposite face's force at node zero.
    let face_force = PEAK_PRESSURE / 6.0;
    let forces = [
        (0, face_force),
        (1, face_force),
        (2, face_force),
        (3, -face_force),
        (7, -face_force),
        (11, -face_force),
    ];

    let bulk = YOUNG / (3.0 * (1.0 - 2.0 * POISSON));
    let normalization = (1.0 + SENSITIVITY * SENSITIVITY / 9.0).sqrt();
    // At hydrostatic yield, phi=a*p/sqrt(c)=Y+H*alpha. Associated flow gives
    // tr(eps_p)=-a*alpha/sqrt(c). These are independent closed-form oracles.
    let expected_alpha = (SENSITIVITY * PEAK_PRESSURE / normalization - YIELD) / HARDENING;
    let permanent_diagonal = -SENSITIVITY * expected_alpha / (3.0 * normalization);
    assert!(expected_alpha > 0.0);

    with_cx(|cx| {
        let mut problem = SmallStrainTetProblem {
            nodes_m: &nodes,
            tetrahedra: &tetrahedra,
            laws: SmallStrainLawField::Uniform(&law),
            prescribed_m: &prescribed,
            nodal_forces_n: &forces,
            load_path: &[0.2, 0.6, 1.0],
            settings: IncrementSettings::default(),
            budget: TetAssemblyBudget::standard(),
        };
        let peak = problem.solve(cx).unwrap();
        close(peak.element_states[0].alpha, expected_alpha, 2.0e-10);
        for axis in 0..3 {
            close(peak.element_stress[0][axis], -PEAK_PRESSURE, 2.0e-3);
            close(
                peak.element_states[0].plastic_strain[axis],
                permanent_diagonal,
                2.0e-10,
            );
            close(
                peak.element_strain[0][axis],
                permanent_diagonal - PEAK_PRESSURE / (3.0 * bulk),
                2.0e-10,
            );
        }

        problem.load_path = &[0.5, 0.0];
        let released = problem
            .solve_from(peak.displacement_m.as_flattened(), &peak.element_states, cx)
            .unwrap();
        assert_eq!(released.element_states, peak.element_states);
        for axis in 0..3 {
            close(released.element_stress[0][axis], 0.0, 2.0e-3);
            close(
                released.element_strain[0][axis],
                permanent_diagonal,
                2.0e-10,
            );
            close(
                released.displacement_m[axis + 1][axis],
                permanent_diagonal,
                2.0e-10,
            );
        }
        let residual_volume_strain: f64 = released.element_strain[0][..3].iter().sum();
        close(residual_volume_strain, 3.0 * permanent_diagonal, 6.0e-10);
        assert!(residual_volume_strain < -0.01);

        // This pressure exceeds virgin yield but stays below the previous
        // maximum. Reloading must follow the elastic slope with retained crush.
        problem.load_path = &[0.5];
        let reloaded = problem
            .solve_from(
                released.displacement_m.as_flattened(),
                &released.element_states,
                cx,
            )
            .unwrap();
        assert_eq!(reloaded.element_states, released.element_states);
        for axis in 0..3 {
            close(
                reloaded.element_stress[0][axis],
                -0.5 * PEAK_PRESSURE,
                2.0e-3,
            );
            close(
                reloaded.element_strain[0][axis],
                permanent_diagonal - 0.5 * PEAK_PRESSURE / (3.0 * bulk),
                2.0e-10,
            );
        }
    });
}
