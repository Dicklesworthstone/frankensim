//! Pointwise stress seeds use mesh-frame Mandel [xx,yy,zz,xy,yz,zx].
//! In particular, shear stresses are sqrt(2) times physical tensor shears.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::linear3::*;
use fs_solid::loads3::transient::adjoint::TransientAdjointError3;
use fs_solid::loads3::transient::adjoint::material::MaterialGradientBudget3;
use fs_solid::loads3::transient::adjoint::support::stress::{StressPullback3, von_mises_limit};
use fs_solid::loads3::{TetLoadGeometry3, VolumeLoad3};

const NODES: [[f64; 3]; 4] = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const ZERO: [[f64; 3]; 4] = [[0.; 3]; 4];
const FRAME: fs_blake3::ContentHash = fs_blake3::ContentHash([153; 32]);

fn with_cx<R>(cancelled: bool, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 153,
                kernel_id: 153,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * (1.0 + a.abs().max(b.abs())),
        "{a:e} != {b:e}"
    );
}
fn thermal(material: &TetElasticMaterial) -> TetThermalStrainState {
    TetThermalStrainState::try_new(
        301.,
        300.,
        [0.0007, 0., 0., 0., 0., 0.],
        material.material_state_identity,
        FRAME,
        0.01,
    )
    .unwrap()
}

#[test]
fn g1_fixed_thermal_stress_has_full_motion_and_direct_material_seeds_without_volume() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6., 60., 0., [153; 32]).unwrap();
        let fixed: Vec<_> = (0..12).collect();
        let problem = TetLinearElasticProblem {
            nodes_m: &NODES,
            tetrahedra: &TETS,
            materials: TetMaterialField::Uniform(&material),
            fixed_dofs: &fixed,
            budget: Default::default(),
        };
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, problem.budget, 4, cx).unwrap();
        let load = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let free = thermal(&material);
        let mut q = ZERO;
        q[1][0] = 0.001;
        q[2][0] = 0.002;
        let run = load
            .begin_transient_with_supports(
                &problem,
                Some(TetThermalStrainField::Uniform(&free)),
                0.,
                &q,
                &ZERO,
                &ZERO,
                Default::default(),
                cx,
            )
            .unwrap();
        let before = run.nodal_motion(cx).unwrap();
        let response = run.recover_response(FRAME, 0.01, cx).unwrap();
        let sigma = response.stress.elements()[0].stress_mandel_pa();
        close(sigma[0], 0.018, 1e-14);
        close(sigma[3], 0.12 / core::f64::consts::SQRT_2, 1e-14);
        let rows = [(0, [2., 0., 0., 3., 0., 0.])];
        let pb = run
            .stress_pullback(&rows, 0.01, Default::default(), cx)
            .unwrap();
        let shear = 180. / core::f64::consts::SQRT_2;
        let expected = [
            [-120. - shear, -shear, 0.],
            [120., shear, 0.],
            [shear, 0., 0.],
            [0., 0., 0.],
        ];
        for (actual, expected) in pb.motion.q.iter().zip(expected) {
            for axis in 0..3 {
                close(actual[axis], expected[axis], 1e-13);
            }
        }
        // Point stress has no 1/6 tetrahedron-volume factor; all these DOFs
        // are constrained, so reducing the seed would discard the objective.
        close(pb.motion.q[1][0], 120., 1e-13);
        assert_eq!(pb.motion.v, ZERO);
        assert_eq!(pb.motion.a, ZERO);
        close(
            pb.material.element_scales[0][0],
            2. * sigma[0] + 3. * sigma[3],
            1e-13,
        );
        assert_eq!(pb.material.element_scales[0][1], 0.);
        assert_eq!(pb.material.rayleigh, [0.; 2]);
        assert_eq!(run.nodal_motion(cx).unwrap(), before);
        let initial = run
            .prescribed_adjoint(Default::default(), Default::default(), cx)
            .unwrap()
            .pullback_initialization(&pb.motion, cx)
            .unwrap();
        assert_eq!(initial.initial.q, pb.motion.q);
        assert!(initial.force.is_empty());
    });
}

// Independent physical perturbations of all full nodal q entries and C's
// scalar multiplier; the forward oracle is the published stress recovery.
fn anisotropic_objective(
    q: [[f64; 3]; 4],
    scale: f64,
    length: f64,
    cx: &Cx<'_>,
) -> (f64, StressPullback3) {
    let mut tensor = [[0.; 6]; 6];
    for (i, diagonal) in [80., 90., 110., 45., 65., 75.].into_iter().enumerate() {
        tensor[i][i] = diagonal;
    }
    for (i, j, value) in [(0, 1, 12.), (0, 3, 5.), (1, 5, -4.), (2, 4, 3.)] {
        tensor[i][j] = value;
        tensor[j][i] = value;
    }
    for row in &mut tensor {
        for entry in row {
            *entry *= scale;
        }
    }
    let material = TetElasticMaterial::try_new_mandel(6., tensor, [154; 32]).unwrap();
    let nodes = NODES.map(|v| v.map(|x| x * length));
    let fixed: Vec<_> = (0..12).collect();
    let problem = TetLinearElasticProblem {
        nodes_m: &nodes,
        tetrahedra: &TETS,
        materials: TetMaterialField::Uniform(&material),
        fixed_dofs: &fixed,
        budget: Default::default(),
    };
    let geometry = TetLoadGeometry3::new(&nodes, &TETS, problem.budget, 4, cx).unwrap();
    let load = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
    let free = thermal(&material);
    let run = load
        .begin_transient_with_supports(
            &problem,
            Some(TetThermalStrainField::Uniform(&free)),
            0.,
            &q,
            &ZERO,
            &ZERO,
            Default::default(),
            cx,
        )
        .unwrap();
    let weights = [0.8, -0.4, 0.3, 1.1, -0.7, 0.6];
    let response = run.recover_response(FRAME, 0.01, cx).unwrap();
    let sigma = response.stress.elements()[0].stress_mandel_pa();
    let value = (0..6).map(|i| weights[i] * sigma[i]).sum();
    let pb = run
        .stress_pullback(&[(0, weights)], 0.01, Default::default(), cx)
        .unwrap();
    (value, pb)
}

#[test]
fn g2_all_mandel_couplings_match_full_q_and_stiffness_differences_at_two_mesh_scales() {
    with_cx(false, |cx| {
        let q = [
            [0.0001, -0.0002, 0.0003],
            [0.001, 0.0004, -0.0002],
            [0.0005, -0.0008, 0.0002],
            [-0.0003, 0.0007, 0.0011],
        ];
        for (scale, length) in [(0.8, 0.7), (1.3, 2.1)] {
            let (_, pb) = anisotropic_objective(q, scale, length, cx);
            for dof in 0..12 {
                let h = 1e-7;
                let mut plus = q;
                let mut minus = q;
                plus[dof / 3][dof % 3] += h;
                minus[dof / 3][dof % 3] -= h;
                let fd = (anisotropic_objective(plus, scale, length, cx).0
                    - anisotropic_objective(minus, scale, length, cx).0)
                    / (2. * h);
                close(pb.motion.q[dof / 3][dof % 3], fd, 2e-9);
            }
            let h = 1e-6;
            let fd = (anisotropic_objective(q, scale + h, length, cx).0
                - anisotropic_objective(q, scale - h, length, cx).0)
                / (2. * h);
            close(pb.material.element_scales[0][0] / scale, fd, 2e-9);
        }
    });
}

#[test]
fn g1_von_mises_limit_uses_mandel_shears_and_has_smooth_hydrostatic_origin() {
    let limit = 4.;
    let (value, seed) = von_mises_limit([limit, 0., 0., 0., 0., 0.], limit).unwrap();
    close(value, 0., 1e-14);
    assert_eq!(seed, [0.5, -0.25, -0.25, 0., 0., 0.]);
    for pressure in [0., 2.5, -3.] {
        assert_eq!(
            von_mises_limit([pressure, pressure, pressure, 0., 0., 0.], limit).unwrap(),
            (-1., [0.; 6])
        );
    }
    // Physical shear tau has Mandel entry sqrt(2)*tau and sigma_vm^2=3*tau^2.
    let tau = 1.7;
    let shear = [0., 0., 0., core::f64::consts::SQRT_2 * tau, 0., 0.];
    close(
        von_mises_limit(shear, limit).unwrap().0,
        3. * tau * tau / (limit * limit) - 1.,
        1e-14,
    );
    let sigma = [2.3, -1.2, 0.4, 1.7, -0.9, 0.2];
    let (value, seed) = von_mises_limit(sigma, limit).unwrap();
    let mut translated = sigma;
    for component in &mut translated[..3] {
        *component += 9.;
    }
    let (translated_value, translated_seed) = von_mises_limit(translated, limit).unwrap();
    close(value, translated_value, 1e-14);
    for i in 0..6 {
        close(seed[i], translated_seed[i], 1e-14);
        let h = 1e-6;
        let mut plus = sigma;
        let mut minus = sigma;
        plus[i] += h;
        minus[i] -= h;
        let fd = (von_mises_limit(plus, limit).unwrap().0
            - von_mises_limit(minus, limit).unwrap().0)
            / (2. * h);
        close(seed[i], fd, 2e-9);
    }
    for invalid in [0., -1., f64::NAN, f64::INFINITY] {
        assert!(von_mises_limit(sigma, invalid).is_err());
    }
    assert!(von_mises_limit([f64::NAN; 6], limit).is_err());
    assert!(von_mises_limit([f64::MAX; 6], f64::MIN_POSITIVE).is_err());
}

#[test]
fn g3_stress_refusals_and_cancelled_retries_leave_accepted_motion_unchanged() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6., 60., 0., [155; 32]).unwrap();
        let fixed: Vec<_> = (0..12).collect();
        let problem = TetLinearElasticProblem {
            nodes_m: &NODES,
            tetrahedra: &TETS,
            materials: TetMaterialField::Uniform(&material),
            fixed_dofs: &fixed,
            budget: Default::default(),
        };
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, problem.budget, 4, cx).unwrap();
        let load = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let free = thermal(&material);
        let mut q = ZERO;
        q[1][0] = 0.001;
        let run = load
            .begin_transient_with_supports(
                &problem,
                Some(TetThermalStrainField::Uniform(&free)),
                0.,
                &q,
                &ZERO,
                &ZERO,
                Default::default(),
                cx,
            )
            .unwrap();
        let before = run.nodal_motion(cx).unwrap();
        let seed = [1.; 6];
        let valid = [(0, seed)];
        let expected = run
            .stress_pullback(&valid, 0.01, Default::default(), cx)
            .unwrap();
        for rows in [
            vec![(1, seed)],
            vec![(0, seed), (0, seed)],
            vec![(0, [f64::NAN; 6])],
        ] {
            assert!(
                run.stress_pullback(&rows, 0.01, Default::default(), cx)
                    .is_err()
            );
        }
        for ceiling in [0., -1., f64::NAN, f64::INFINITY, 0.0005] {
            // Empty objective rows still cannot bypass the full strain gate.
            assert!(
                run.stress_pullback(&[], ceiling, Default::default(), cx)
                    .is_err()
            );
        }
        for budget in [
            MaterialGradientBudget3 {
                maximum_full_dofs: 11,
                ..Default::default()
            },
            MaterialGradientBudget3 {
                maximum_elements: 0,
                ..Default::default()
            },
        ] {
            assert_eq!(
                run.stress_pullback(&valid, 0.01, budget, cx),
                Err(TransientAdjointError3::Budget)
            );
        }
        with_cx(true, |cancelled| {
            assert_eq!(
                run.stress_pullback(&valid, 0.01, Default::default(), cancelled),
                Err(TransientAdjointError3::Cancelled),
            )
        });
        assert_eq!(run.nodal_motion(cx).unwrap(), before);
        assert_eq!(
            run.stress_pullback(&valid, 0.01, Default::default(), cx)
                .unwrap(),
            expected
        );
        // Total strain is now zero, but the frozen free strain produces an
        // elastic strain of -0.0007, which must also obey the same ceiling.
        let thermal_only = load
            .begin_transient_with_supports(
                &problem,
                Some(TetThermalStrainField::Uniform(&free)),
                0.,
                &ZERO,
                &ZERO,
                &ZERO,
                Default::default(),
                cx,
            )
            .unwrap();
        assert!(
            thermal_only
                .stress_pullback(&[], 0.0005, Default::default(), cx)
                .is_err()
        );
        let empty = run
            .stress_pullback(&[], 0.01, Default::default(), cx)
            .unwrap();
        assert_eq!(empty.motion.q, ZERO);
        assert_eq!(empty.material.element_scales, vec![[0.; 2]]);
    });
}
