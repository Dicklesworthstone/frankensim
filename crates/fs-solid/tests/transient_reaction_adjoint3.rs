//! Reactions retain direct force/material terms as well as motion dependence.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::linear3::*;
use fs_solid::loads3::transient::TransientConfig3;
use fs_solid::loads3::transient::adjoint::TransientAdjointError3;
use fs_solid::loads3::transient::adjoint::material::{
    MaterialGradientBudget3, MaterialScaleGradient3,
};
use fs_solid::loads3::{TetLoadGeometry3, VolumeLoad3};

const NODES: [[f64; 3]; 4] = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const ZERO: [[f64; 3]; 4] = [[0.; 3]; 4];
const FRAME: fs_blake3::ContentHash = fs_blake3::ContentHash([149; 32]);

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
                seed: 149,
                kernel_id: 149,
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
fn add(target: &mut MaterialScaleGradient3, value: MaterialScaleGradient3) {
    for (sum, entry) in target.element_scales.iter_mut().zip(value.element_scales) {
        for axis in 0..2 {
            sum[axis] += entry[axis];
        }
    }
    for axis in 0..2 {
        target.rayleigh[axis] += value.rayleigh[axis];
    }
}

#[test]
fn g1_initial_reaction_pullback_includes_direct_and_solved_acceleration_terms() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6., 60., 0., [149; 32]).unwrap();
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
            .assemble(&[], VolumeLoad3::None, &[(0, 0.011), (3, 0.017)], cx)
            .unwrap();
        let thermal = thermal(&material);
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
        let weight = 1.3;
        let run = load
            .begin_transient_with_supports(
                &problem,
                Some(TetThermalStrainField::Uniform(&thermal)),
                0.,
                &q,
                &v,
                &a,
                config(eta, zeta),
                cx,
            )
            .unwrap();
        let before = run.nodal_motion(cx).unwrap();
        let response = run.recover_response(FRAME, 0.01, cx).unwrap();
        let af = run.state().acceleration[0];
        // The fixed corner has K_00=20; its coupling to the one free x DOF
        // is -10. Consistent mass contributes M_00=.1 and M_03=.05.
        let ku = 20. * q[0][0] - 10. * q[1][0];
        let kv = 20. * v[0][0] - 10. * v[1][0];
        let mv = 0.1 * v[0][0] + 0.05 * v[1][0];
        let ma = 0.1 * a[0][0] + 0.05 * af;
        close(
            response.support_reactions_n[0][0],
            ma + eta * mv + zeta * kv + ku + 0.007 - 0.011,
            1e-13,
        );
        let pb = run
            .reaction_pullback(&[(0, weight)], Default::default(), cx)
            .unwrap();
        close(pb.motion.q[0][0], 20. * weight, 1e-13);
        close(pb.motion.q[1][0], -10. * weight, 1e-13);
        close(pb.motion.v[0][0], (0.1 * eta + 20. * zeta) * weight, 1e-13);
        close(pb.motion.v[1][0], (0.05 * eta - 10. * zeta) * weight, 1e-13);
        close(pb.motion.a[0][0], 0.1 * weight, 1e-13);
        close(pb.motion.a[1][0], 0.05 * weight, 1e-13);
        assert_eq!(pb.force, vec![(0, -weight)]);
        close(
            pb.material.element_scales[0][0],
            weight * (ku + zeta * kv + 0.007),
            1e-13,
        );
        close(
            pb.material.element_scales[0][1],
            weight * (ma + eta * mv),
            1e-13,
        );
        close(pb.material.rayleigh[0], weight * mv, 1e-13);
        close(pb.material.rayleigh[1], weight * kv, 1e-13);
        let adjoint = run
            .prescribed_adjoint(Default::default(), Default::default(), cx)
            .unwrap();
        assert_eq!(
            pb.motion,
            adjoint
                .reaction_motion_pullback(&[(0, weight)], cx)
                .unwrap()
        );
        let stationary = load
            .begin_transient(
                &problem,
                Some(TetThermalStrainField::Uniform(&thermal)),
                0.,
                &ZERO,
                &ZERO,
                config(eta, zeta),
                cx,
            )
            .unwrap();
        assert_eq!(
            pb.motion,
            stationary
                .reaction_pullback(&[(0, weight)], Default::default(), cx)
                .unwrap()
                .motion
        );
        let initial = adjoint.pullback_initialization(&pb.motion, cx).unwrap();
        close(initial.force[0], 0.5 * weight, 1e-13);
        let mut total = pb.material;
        add(
            &mut total,
            run.material_scale_pullback(&initial.force, Default::default(), cx)
                .unwrap(),
        );
        close(
            total.element_scales[0][0],
            weight
                * (25. * q[0][0] - 15. * q[1][0] + zeta * (25. * v[0][0] - 15. * v[1][0]) + 0.0105),
            1e-13,
        );
        close(
            total.element_scales[0][1],
            weight * 0.075 * (a[0][0] + eta * v[0][0]),
            1e-13,
        );
        close(total.rayleigh[0], weight * 0.075 * v[0][0], 1e-13);
        close(
            total.rayleigh[1],
            weight * (25. * v[0][0] - 15. * v[1][0]),
            1e-13,
        );
        for rows in [
            vec![(3, 1.)],
            vec![(12, 1.)],
            vec![(0, 1.), (0, 2.)],
            vec![(1, 1.), (0, 1.)],
            vec![(0, f64::NAN)],
        ] {
            assert!(
                run.reaction_pullback(&rows, Default::default(), cx)
                    .is_err()
            );
            assert!(adjoint.reaction_motion_pullback(&rows, cx).is_err());
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
                run.reaction_pullback(&[(0, weight)], budget, cx),
                Err(TransientAdjointError3::Budget)
            );
        }
        with_cx(true, |cancelled| {
            assert_eq!(
                run.reaction_pullback(&[(0, weight)], Default::default(), cancelled),
                Err(TransientAdjointError3::Cancelled)
            )
        });
        assert_eq!(run.nodal_motion(cx).unwrap(), before);
        assert!(
            run.reaction_pullback(&[(0, weight)], Default::default(), cx)
                .is_ok()
        );
    });
}

const STEPS: [f64; 3] = [0.025, 0.04, 0.02];
const SUPPORT: [f64; 4] = [0.0001, 0.00012, 0.00020, 0.00022];
const FREE_FORCE: [f64; 4] = [0.0005, 0.002, -0.001, 0.0015];
const FIXED_FORCE: [f64; 4] = [0.001, -0.002, 0.0015, 0.003];

// Stiffness, density, eta, zeta, support amplitude, fixed-force amplitude,
// free-force amplitude. A fixed-node force changes reactions without motion.
fn trajectory(parameters: [f64; 7], cx: &Cx<'_>) -> (f64, [f64; 7]) {
    let material =
        TetElasticMaterial::try_new(6. * parameters[1], 60. * parameters[0], 0.25, [150; 32])
            .unwrap();
    let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
    let problem = TetLinearElasticProblem {
        nodes_m: &NODES,
        tetrahedra: &TETS,
        materials: TetMaterialField::Uniform(&material),
        fixed_dofs: &fixed,
        budget: Default::default(),
    };
    let thermal = thermal(&material);
    let geometry = TetLoadGeometry3::new(&NODES, &TETS, problem.budget, 4, cx).unwrap();
    let load = geometry
        .assemble(
            &[],
            VolumeLoad3::None,
            &[
                (0, FIXED_FORCE[0] * parameters[5]),
                (3, FREE_FORCE[0] * parameters[6]),
            ],
            cx,
        )
        .unwrap();
    let mut q = ZERO;
    q[0][0] = SUPPORT[0] * parameters[4];
    q[1][0] = 0.00015;
    let mut v = ZERO;
    v[0][0] = -0.0005 * parameters[4];
    v[1][0] = 0.001;
    let mut a = ZERO;
    a[0][0] = 0.002 * parameters[4];
    let mut run = load
        .begin_transient_with_supports(
            &problem,
            Some(TetThermalStrainField::Uniform(&thermal)),
            0.,
            &q,
            &v,
            &a,
            config(parameters[2], parameters[3]),
            cx,
        )
        .unwrap();
    let adjoint = run
        .prescribed_adjoint(Default::default(), Default::default(), cx)
        .unwrap();
    let mut accepted = vec![run.clone()];
    for (i, h) in STEPS.into_iter().enumerate() {
        let load = geometry
            .assemble(
                &[],
                VolumeLoad3::None,
                &[
                    (0, FIXED_FORCE[i + 1] * parameters[5]),
                    (3, FREE_FORCE[i + 1] * parameters[6]),
                ],
                cx,
            )
            .unwrap();
        run.advance_with_supports(&load, &[(0, SUPPORT[i + 1] * parameters[4])], h, cx)
            .unwrap();
        accepted.push(run.clone());
    }
    let reaction = run
        .recover_response(FRAME, 0.01, cx)
        .unwrap()
        .support_reactions_n[0][0];
    let residual = (reaction - 0.005) / 0.02;
    let direct = run
        .reaction_pullback(&[(0, residual / 0.02)], Default::default(), cx)
        .unwrap();
    let fixed_force_gradient = direct.force[0].1 * FIXED_FORCE[3];
    let mut material_gradient = direct.material;
    let mut cotangent = direct.motion;
    let mut support_gradient = 0.;
    let mut free_force_gradient = 0.;
    let support_index = adjoint.fixed_dofs().binary_search(&0).unwrap();
    for i in (0..STEPS.len()).rev() {
        let pb = adjoint.pullback_step(STEPS[i], &cotangent, cx).unwrap();
        add(
            &mut material_gradient,
            accepted[i + 1]
                .material_scale_pullback(&pb.force, Default::default(), cx)
                .unwrap(),
        );
        support_gradient += pb.prescribed_displacement[support_index] * SUPPORT[i + 1];
        free_force_gradient += pb.force[0] * FREE_FORCE[i + 1];
        cotangent = pb.previous;
    }
    let initial = adjoint.pullback_initialization(&cotangent, cx).unwrap();
    add(
        &mut material_gradient,
        accepted[0]
            .material_scale_pullback(&initial.force, Default::default(), cx)
            .unwrap(),
    );
    support_gradient += initial.initial.q[0][0] * SUPPORT[0] - 0.0005 * initial.initial.v[0][0]
        + 0.002 * initial.initial.a[0][0];
    free_force_gradient += initial.force[0] * FREE_FORCE[0];
    (
        0.5 * residual * residual,
        [
            material_gradient.element_scales[0][0] / parameters[0],
            material_gradient.element_scales[0][1] / parameters[1],
            material_gradient.rayleigh[0],
            material_gradient.rayleigh[1],
            support_gradient,
            fixed_force_gradient,
            free_force_gradient,
        ],
    )
}

#[test]
fn g3_reaction_target_gradient_matches_material_and_both_load_control_paths() {
    with_cx(false, |cx| {
        for parameters in [
            [1., 1., 0.4, 0.006, 1., 1., 1.],
            [1.6, 0.7, 0.8, 0.009, 0.8, 1.2, 0.9],
        ] {
            let (_, gradient) = trajectory(parameters, cx);
            for axis in 0..7 {
                let h = if axis == 3 {
                    1e-6
                } else {
                    1e-5 * parameters[axis]
                };
                let mut plus = parameters;
                let mut minus = parameters;
                plus[axis] += h;
                minus[axis] -= h;
                let fd = (trajectory(plus, cx).0 - trajectory(minus, cx).0) / (2. * h);
                close(gradient[axis], fd, 3e-6);
            }
        }
    });
}

#[test]
fn g1_fully_prescribed_reactions_keep_direct_thermal_mass_and_damping_derivatives() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6., 60., 0., [151; 32]).unwrap();
        let fixed: Vec<_> = (0..12).collect();
        let problem = TetLinearElasticProblem {
            nodes_m: &NODES,
            tetrahedra: &TETS,
            materials: TetMaterialField::Uniform(&material),
            fixed_dofs: &fixed,
            budget: Default::default(),
        };
        let thermal = thermal(&material);
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, problem.budget, 4, cx).unwrap();
        let load = geometry
            .assemble(&[], VolumeLoad3::None, &[(0, 0.001)], cx)
            .unwrap();
        let run = load
            .begin_transient_with_supports(
                &problem,
                Some(TetThermalStrainField::Uniform(&thermal)),
                0.,
                &ZERO,
                &[[0.003, 0., 0.]; 4],
                &[[0.02, 0., 0.]; 4],
                config(0.3, 0.004),
                cx,
            )
            .unwrap();
        let reaction = run
            .recover_response(FRAME, 0.01, cx)
            .unwrap()
            .support_reactions_n[0][0];
        close(reaction, 0.005 + 0.000225 + 0.007 - 0.001, 1e-13);
        let pb = run
            .reaction_pullback(&[(0, 1.)], Default::default(), cx)
            .unwrap();
        close(pb.material.element_scales[0][0], 0.007, 1e-13);
        close(pb.material.element_scales[0][1], 0.005225, 1e-13);
        close(pb.material.rayleigh[0], 0.00075, 1e-13);
        close(pb.material.rayleigh[1], 0., 1e-13);
        let adjoint = run
            .prescribed_adjoint(Default::default(), Default::default(), cx)
            .unwrap();
        let initial = adjoint.pullback_initialization(&pb.motion, cx).unwrap();
        assert!(initial.force.is_empty());
        assert_eq!(initial.iterations, 0);
        let implicit = run
            .material_scale_pullback(&initial.force, Default::default(), cx)
            .unwrap();
        assert_eq!(implicit.element_scales, vec![[0.; 2]]);
        assert_eq!(implicit.rayleigh, [0.; 2]);
        assert_eq!(pb.force, vec![(0, -1.)]);
    });
}
