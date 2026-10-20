//! Physical material sensitivities with actual driven-support kinematics.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::linear3::*;
use fs_solid::loads3::transient::adjoint::TransientAdjointError3;
use fs_solid::loads3::transient::adjoint::material::{
    MaterialGradientBudget3, MaterialScaleGradient3,
};
use fs_solid::loads3::transient::adjoint::support::NodalMotionCotangent3;
use fs_solid::loads3::transient::{TetTransientRun3, TransientConfig3};
use fs_solid::loads3::{TetLoadGeometry3, VolumeLoad3};

const NODES: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const ZERO: [[f64; 3]; 4] = [[0.0; 3]; 4];

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
                seed: 141,
                kernel_id: 141,
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
fn config(eta: f64, zeta: f64) -> TransientConfig3 {
    TransientConfig3 {
        mass_damping_per_s: eta,
        stiffness_damping_s: zeta,
        ..Default::default()
    }
}
fn thermal(material: &TetElasticMaterial, strain: [f64; 6]) -> TetThermalStrainState {
    TetThermalStrainState::try_new(
        301.0,
        300.0,
        strain,
        material.material_state_identity,
        fs_blake3::ContentHash([143; 32]),
        0.01,
    )
    .unwrap()
}
fn seed(nodes: usize) -> NodalMotionCotangent3 {
    NodalMotionCotangent3 {
        q: (0..nodes)
            .map(|n| std::array::from_fn(|c| 0.2 + 0.017 * (3 * n + c) as f64))
            .collect(),
        v: (0..nodes)
            .map(|n| std::array::from_fn(|c| -0.05 + 0.013 * (3 * n + c) as f64))
            .collect(),
        a: (0..nodes)
            .map(|n| std::array::from_fn(|c| 0.006 - 0.001 * (3 * n + c) as f64))
            .collect(),
    }
}
fn objective(
    run: &TetTransientRun3<'_, '_, '_, '_>,
    seed: &NodalMotionCotangent3,
    cx: &Cx<'_>,
) -> f64 {
    let motion = run.nodal_motion(cx).unwrap();
    [&seed.q, &seed.v, &seed.a]
        .into_iter()
        .zip(&motion)
        .flat_map(|(s, m)| s.iter().flatten().zip(m.iter().flatten()))
        .map(|(s, m)| s * m)
        .sum()
}
fn add(target: &mut MaterialScaleGradient3, value: MaterialScaleGradient3) {
    for (sum, entry) in target.element_scales.iter_mut().zip(value.element_scales) {
        for c in 0..2 {
            sum[c] += entry[c];
        }
    }
    for c in 0..2 {
        target.rayleigh[c] += value.rayleigh[c];
    }
}

#[test]
fn g1_initial_material_gradient_keeps_prescribed_mass_stiffness_and_damping() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6.0, 60.0, 0.0, [141; 32]).unwrap();
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
        let problem = TetLinearElasticProblem {
            nodes_m: &NODES,
            tetrahedra: &TETS,
            materials: TetMaterialField::Uniform(&material),
            fixed_dofs: &fixed,
            budget: Default::default(),
        };
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, problem.budget, 4, cx).unwrap();
        let load = geometry
            .assemble(&[], VolumeLoad3::None, &[(3, 0.017)], cx)
            .unwrap();
        let thermal = thermal(&material, [0.0007, 0.0, 0.0, 0.0, 0.0, 0.0]);
        let mut q = ZERO;
        q[0][0] = 0.0002;
        q[1][0] = 0.001;
        let mut v = ZERO;
        v[0][0] = -0.003;
        v[1][0] = 0.002;
        let mut a = ZERO;
        a[0][0] = 0.004;
        let eta = 0.3;
        let zeta = 0.004;
        let run = load
            .begin_transient_with_supports(
                &problem,
                Some(TetThermalStrainField::Uniform(&thermal)),
                0.0,
                &q,
                &v,
                &a,
                config(eta, zeta),
                cx,
            )
            .unwrap();
        let adjoint = run
            .prescribed_adjoint(Default::default(), Default::default(), cx)
            .unwrap();
        let mut s = seed(4);
        s.a[1][0] = 0.04;
        let pb = adjoint.pullback_initialization(&s, cx).unwrap();
        let gradient = run
            .material_scale_pullback(&pb.force, Default::default(), cx)
            .unwrap();
        // Unit P1 tet: M_ff=.1, M_fc=.05, K_ff=10, K_fc=-10.
        let mv = 0.1 * v[1][0] + 0.05 * v[0][0];
        let kv = 10.0 * (v[1][0] - v[0][0]);
        let ku = 10.0 * (q[1][0] - q[0][0]);
        let ft = 10.0 * 0.0007;
        let af = (0.017 + ft - ku - eta * mv - zeta * kv - 0.05 * a[0][0]) / 0.1;
        let lambda = 0.04 / 0.1;
        close(run.state().acceleration[0], af, 1e-13);
        close(pb.force[0], lambda, 1e-13);
        close(
            gradient.element_scales[0][0],
            -lambda * (ku + zeta * kv - ft),
            1e-13,
        );
        close(
            gradient.element_scales[0][1],
            -lambda * (0.1 * af + 0.05 * a[0][0] + eta * mv),
            1e-13,
        );
        close(gradient.rayleigh[0], -lambda * mv, 1e-13);
        close(gradient.rayleigh[1], -lambda * kv, 1e-13);
        // Deleting the constrained terms cannot satisfy these checks.
        let zero_padded_mass = -lambda * (0.1 * af + eta * 0.1 * v[1][0]);
        let zero_padded_stiffness = -lambda * 10.0 * (q[1][0] + zeta * v[1][0] - 0.0007);
        assert!((gradient.element_scales[0][1] - zero_padded_mass).abs() > 1e-5);
        assert!((gradient.element_scales[0][0] - zero_padded_stiffness).abs() > 1e-4);
    });
}

fn multi_step(parameters: [f64; 6], cx: &Cx<'_>) -> (f64, [f64; 6]) {
    let nodes = [
        [0., 0., 0.],
        [1., 0., 0.],
        [0., 1., 0.],
        [0., 0., 1.],
        [1., 1., 1.],
    ];
    let tets = [[0, 1, 2, 3], [1, 2, 3, 4]];
    let fixed = [0, 1, 2, 4, 5, 8];
    let materials: Vec<_> = (0..2)
        .map(|e| {
            // A full SPD Mandel tensor, with independent anisotropy in each tet.
            let c = std::array::from_fn(|i| {
                std::array::from_fn(|j| {
                    parameters[2 * e]
                        * ((if i == j { 40.0 + 10.0 * i as f64 } else { 0.0 })
                            + (i + 1) as f64 * (j + 1) as f64 * (e + 1) as f64)
                })
            });
            TetElasticMaterial::try_new_mandel(
                (6.0 + 3.0 * e as f64) * parameters[2 * e + 1],
                c,
                [141 + e as u8; 32],
            )
            .unwrap()
        })
        .collect();
    let thermal: Vec<_> = materials
        .iter()
        .enumerate()
        .map(|(e, m)| {
            thermal(
                m,
                std::array::from_fn(|i| 0.00002 * (e + 1) as f64 * (i + 1) as f64),
            )
        })
        .collect();
    let problem = TetLinearElasticProblem {
        nodes_m: &nodes,
        tetrahedra: &tets,
        materials: TetMaterialField::PerElement(&materials),
        fixed_dofs: &fixed,
        budget: Default::default(),
    };
    let geometry = TetLoadGeometry3::new(&nodes, &tets, problem.budget, 8, cx).unwrap();
    let mut q = [[0.0; 3]; 5];
    let mut v = q;
    let mut a = q;
    for d in 0..15 {
        q[d / 3][d % 3] = 0.00004 * (d + 1) as f64;
        v[d / 3][d % 3] = -0.0001 * (d + 1) as f64;
        if fixed.contains(&d) {
            a[d / 3][d % 3] = 0.0008 * (d + 1) as f64;
        }
    }
    let load = geometry
        .assemble(
            &[],
            VolumeLoad3::Uniform([0.02, -0.03, 0.01]),
            &[(3, 0.017)],
            cx,
        )
        .unwrap();
    let mut run = load
        .begin_transient_with_supports(
            &problem,
            Some(TetThermalStrainField::PerElement(&thermal)),
            0.0,
            &q,
            &v,
            &a,
            config(parameters[4], parameters[5]),
            cx,
        )
        .unwrap();
    let adjoint = run
        .prescribed_adjoint(Default::default(), Default::default(), cx)
        .unwrap();
    let mut accepted = vec![run.clone()];
    let intervals = [0.02, 0.03, 0.025];
    let mut time = 0.0;
    for (i, h) in intervals.into_iter().enumerate() {
        time += h;
        let supports: Vec<_> = fixed
            .iter()
            .map(|&d| {
                (
                    d,
                    q[d / 3][d % 3] + time * v[d / 3][d % 3] + 0.5 * time * time * a[d / 3][d % 3],
                )
            })
            .collect();
        let load = geometry
            .assemble(
                &[],
                VolumeLoad3::Uniform([0.02, -0.03, 0.01]),
                &[(3, 0.013 - 0.01 * i as f64)],
                cx,
            )
            .unwrap();
        run.advance_with_supports(&load, &supports, h, cx).unwrap();
        accepted.push(run.clone());
    }
    let mut cotangent = seed(nodes.len());
    let value = objective(&run, &cotangent, cx);
    let mut gradient = MaterialScaleGradient3 {
        element_scales: vec![[0.0; 2]; 2],
        rayleigh: [0.0; 2],
    };
    for i in (0..intervals.len()).rev() {
        let pb = adjoint.pullback_step(intervals[i], &cotangent, cx).unwrap();
        add(
            &mut gradient,
            accepted[i + 1]
                .material_scale_pullback(&pb.force, Default::default(), cx)
                .unwrap(),
        );
        cotangent = pb.previous;
    }
    let initial = adjoint.pullback_initialization(&cotangent, cx).unwrap();
    add(
        &mut gradient,
        accepted[0]
            .material_scale_pullback(&initial.force, Default::default(), cx)
            .unwrap(),
    );
    (
        value,
        [
            gradient.element_scales[0][0] / parameters[0],
            gradient.element_scales[0][1] / parameters[1],
            gradient.element_scales[1][0] / parameters[2],
            gradient.element_scales[1][1] / parameters[3],
            gradient.rayleigh[0],
            gradient.rayleigh[1],
        ],
    )
}

#[test]
fn g3_heterogeneous_anisotropic_driven_trajectory_matches_six_material_directions() {
    with_cx(false, |cx| {
        for parameters in [
            [1.0, 1.0, 1.0, 1.0, 0.3, 0.004],
            [1.7, 0.6, 0.8, 1.4, 0.7, 0.009],
        ] {
            let (_, gradient) = multi_step(parameters, cx);
            for c in 0..6 {
                let epsilon = if c == 5 { 1e-6 } else { 1e-5 * parameters[c] };
                let mut plus = parameters;
                let mut minus = parameters;
                plus[c] += epsilon;
                minus[c] -= epsilon;
                let fd = (multi_step(plus, cx).0 - multi_step(minus, cx).0) / (2.0 * epsilon);
                close(gradient[c], fd, 2e-7);
            }
        }
    });
}

#[test]
fn g0_g4_fully_prescribed_material_gradient_is_zero_and_refusals_are_reusable() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6.0, 60.0, 0.25, [141; 32]).unwrap();
        let fixed: Vec<_> = (0..12).collect();
        let problem = TetLinearElasticProblem {
            nodes_m: &NODES,
            tetrahedra: &TETS,
            materials: TetMaterialField::Uniform(&material),
            fixed_dofs: &fixed,
            budget: Default::default(),
        };
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, problem.budget, 4, cx).unwrap();
        let load = geometry
            .assemble(&[], VolumeLoad3::None, &[(3, 0.017)], cx)
            .unwrap();
        let thermal = thermal(&material, [0.0007, 0.0, 0.0, 0.0, 0.0, 0.0]);
        let q = [[0.0001, -0.0002, 0.0003]; 4];
        let v = [[0.001, -0.002, 0.003]; 4];
        let a = [[0.002, -0.003, 0.004]; 4];
        let mut run = load
            .begin_transient_with_supports(
                &problem,
                Some(TetThermalStrainField::Uniform(&thermal)),
                0.0,
                &q,
                &v,
                &a,
                config(0.3, 0.004),
                cx,
            )
            .unwrap();
        let h = 0.02;
        let supports: Vec<_> = fixed
            .iter()
            .map(|&d| {
                (
                    d,
                    q[d / 3][d % 3] + h * v[d / 3][d % 3] + 0.5 * h * h * a[d / 3][d % 3],
                )
            })
            .collect();
        run.advance_with_supports(&load, &supports, h, cx).unwrap();
        let before = run.state().clone();
        let full_before = run.nodal_motion(cx).unwrap();
        let adjoint = run
            .prescribed_adjoint(Default::default(), Default::default(), cx)
            .unwrap();
        let pb = adjoint.pullback_step(h, &seed(4), cx).unwrap();
        assert!(pb.force.is_empty());
        let expected = MaterialScaleGradient3 {
            element_scales: vec![[0.0; 2]],
            rayleigh: [0.0; 2],
        };
        assert_eq!(
            run.material_scale_pullback(&pb.force, Default::default(), cx)
                .unwrap(),
            expected
        );
        assert!(
            run.material_scale_pullback(&[f64::NAN], Default::default(), cx)
                .is_err()
        );
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
                run.material_scale_pullback(&[], budget, cx),
                Err(TransientAdjointError3::Budget)
            );
        }
        with_cx(true, |cancelled| {
            assert_eq!(
                run.material_scale_pullback(&[], Default::default(), cancelled),
                Err(TransientAdjointError3::Cancelled)
            )
        });
        assert_eq!(run.state(), &before);
        assert_eq!(run.nodal_motion(cx).unwrap(), full_before);
        assert_eq!(
            run.material_scale_pullback(&[], Default::default(), cx)
                .unwrap(),
            expected
        );
    });
}
