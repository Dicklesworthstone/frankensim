//! Shared spatial plate ports, dynamics, admission, and calibration for two bodies.
#[path = "../examples/calibrate_plate_impact.rs"]
#[allow(dead_code)]
mod calibration;
#[path = "../examples/material_plate_impact.rs"]
#[allow(dead_code)]
mod sourced;

use calibration::{DrivingLoad, adjoint_config, driving_load, method};
use fs_dcontact::Obstacle;
use fs_plate::{
    PlateMesh, PlateModel,
    loading::{PlateLoadBudget, PlatePatchBudget, PlatePatchPoint, PlatePatchStencil},
};
use fs_plate_transient::{
    NoPlateLoad, PlateDynamics, PlateDynamicsBudget, PlateDynamicsParameters, PlateImpact,
    PlateImpactBudget, PlateImpactError, PlateImpactSet, PlateImpactSetBudget, PlateLoad,
};
use fs_solver::FlexiblePreconditioner;
use fs_sparse::Coo;
use fs_time::galpha::{
    SecondOrderProblem, SecondOrderState, TimeSolveError,
    initialization::{InitialSolveConfig, initial_workspace_components, second_order_acceleration},
    second_order_adjoint::{SecondOrderVjp, trajectory::StructuralTrajectoryModel},
};

// Common plate/load prefix, then [mass, K, chi, gap] for each body.
const PARAMETERS: [f64; 13] = [
    1.1, 0.9, 0.04, 0.02, 0.7, 0.2, 700.0, 3.0, 0.005, 0.35, 1100.0, 1.2, 0.008,
];
const POINTS: [[f64; 2]; 2] = [[0.2, 0.2], [0.6, 0.2]];

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * (1.0 + b.abs()),
        "{a:.16e} != {b:.16e}"
    );
}
fn close_vectors(a: &[f64], b: &[f64], tolerance: f64) {
    assert_eq!(a.len(), b.len());
    for (&a, &b) in a.iter().zip(b) {
        close(a, b, tolerance);
    }
}
fn budget() -> PlateImpactSetBudget {
    PlateImpactSetBudget {
        max_bodies: 4,
        max_free_dofs: 128,
        max_contact_points: 16,
        max_parameters: 32,
    }
}

// The same analytic spatial plate pencil used by the single-impact regressions.
fn fixture(stiffness: f64) -> (PlateMesh, PlateModel) {
    let mesh =
        PlateMesh::from_unstructured(vec![(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)], vec![[0, 1, 2]])
            .unwrap();
    let mut k = Coo::new(9, 9);
    let mut m = Coo::new(9, 9);
    for i in 0..9 {
        k.push(i, i, stiffness * (1.0 + i as f64 / 10.0));
        m.push(i, i, 2.0 + i as f64 / 10.0);
    }
    (
        mesh,
        PlateModel {
            k: k.assemble(),
            m: m.assemble(),
            dof_map: (0..9).map(Some).collect(),
            free: 9,
        },
    )
}

fn body<'a, L: PlateLoad + ?Sized>(
    mesh: &PlateMesh,
    plate: &'a PlateModel,
    load: &'a L,
    p: [f64; 4],
    contact: [f64; 4],
    point: [f64; 2],
) -> PlateImpact<'a, L> {
    let dynamics = PlateDynamics::new(
        plate,
        PlateDynamicsParameters {
            stiffness_scale: p[0],
            mass_scale: p[1],
            mass_damping_per_s: p[2],
            stiffness_damping_s: p[3],
        },
        load,
        PlateDynamicsBudget {
            max_dofs: plate.free,
            max_nonzeros: plate.k.nnz() + plate.m.nnz(),
            max_full_dofs: plate.dof_map.len(),
            max_load_parameters: load.parameter_count(),
        },
        &mut || false,
    )
    .unwrap();
    let patch = PlatePatchStencil::locate(
        mesh,
        plate,
        &[PlatePatchPoint {
            point,
            area_m2: 0.004,
        }],
        PlatePatchBudget {
            max_points: 1,
            location: PlateLoadBudget {
                max_nodes: mesh.nodes.len(),
                max_triangles: mesh.tris.len(),
            },
        },
    )
    .unwrap();
    let law = Obstacle::new(
        vec![1.0],
        1,
        1,
        vec![contact[3]],
        vec![1.0],
        contact[1],
        1.5,
        "authored multiple-impact regression".into(),
    )
    .unwrap()
    .with_internal_loss(contact[2])
    .unwrap();
    PlateImpact::new(
        dynamics,
        patch,
        contact[0],
        law,
        PlateImpactBudget {
            max_free_dofs: plate.free + 1,
            max_contact_points: 1,
        },
    )
    .unwrap()
}

fn pair<'a>(
    mesh: &PlateMesh,
    plate: &'a PlateModel,
    load: &'a DrivingLoad,
    p: [f64; 13],
    reversed: bool,
) -> PlateImpactSet<'a, DrivingLoad> {
    let mut bodies: Vec<_> = (0..2)
        .map(|i| {
            let offset = 5 + 4 * i;
            body(
                mesh,
                plate,
                load,
                p[..4].try_into().unwrap(),
                p[offset..offset + 4].try_into().unwrap(),
                POINTS[i],
            )
            .with_contact_calibration()
            .unwrap()
        })
        .collect();
    if reversed {
        bodies.reverse();
    }
    PlateImpactSet::new(bodies, budget()).unwrap()
}
fn state_vectors() -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut q: Vec<_> = (0..11).map(|i| 0.001 * (i as f64 - 2.0)).collect();
    let mut v: Vec<_> = (0..11).map(|i| 0.003 * (i as f64 - 4.0)).collect();
    let a = (0..11).map(|i| 0.02 * (i as f64 - 3.0)).collect();
    q[9] = 0.04;
    q[10] = 0.05;
    v[9] = 0.15;
    v[10] = -0.03;
    (q, v, a)
}
fn swap_bodies(values: &[f64]) -> Vec<f64> {
    let mut result = values.to_vec();
    let n = result.len();
    result.swap(n - 2, n - 1);
    result
}
fn actions<M: SecondOrderProblem>(model: &M, q: &[f64], v: &[f64], d: &[f64]) -> Vec<Vec<f64>> {
    let mut result = vec![vec![f64::NAN; model.dimension()]; 7];
    model.mass_apply(d, &mut result[0]);
    model.damping_apply(d, &mut result[1]);
    model.internal_force(q, &mut result[2]);
    model.tangent_apply(q, d, &mut result[3]);
    model.resistance(q, v, &mut result[4]);
    model.resistance_displacement_apply(q, v, d, &mut result[5]);
    model.resistance_velocity_apply(q, v, d, &mut result[6]);
    result
}
fn consistent_state<L: PlateLoad + ?Sized>(
    model: &PlateImpactSet<'_, L>,
    q: &[f64],
    v: &[f64],
) -> SecondOrderState {
    let n = model.dimension();
    let mut force = vec![0.0; n];
    model.forcing(0.0, &mut force).unwrap();
    let config = InitialSolveConfig {
        restart: n,
        max_cycles: 4,
        tolerance: 1e-12,
    };
    let a = second_order_acceleration(
        model,
        q,
        v,
        &force,
        &model.mass_preconditioner(),
        config,
        initial_workspace_components(n, model.parameter_count(), config).unwrap(),
        &mut || false,
    )
    .unwrap();
    SecondOrderState::new(0.0, q, v, &a.value)
}

#[test]
fn singleton_preserves_original_actions_parameters_probes_and_step() {
    let (mesh, plate) = fixture(4.0);
    let load = driving_load(&mesh, &plate, PARAMETERS[4]).unwrap();
    let make_body = || {
        body(
            &mesh,
            &plate,
            &load,
            PARAMETERS[..4].try_into().unwrap(),
            PARAMETERS[5..9].try_into().unwrap(),
            POINTS[0],
        )
        .with_contact_calibration()
        .unwrap()
    };
    let original = make_body();
    let set = PlateImpactSet::new(vec![make_body()], budget()).unwrap();
    assert_eq!(
        (set.plate_dofs(), set.body_count(), set.dimension()),
        (9, 1, 10)
    );
    assert_eq!(set.parameter_count(), original.parameter_count());
    assert_eq!(set.contact_parameter_range(0), Some(5..9));
    assert!(set.body(1).is_none() && set.contact_parameter_range(1).is_none());
    assert!(std::ptr::eq(set.plate().model(), &plate));
    let (mut q, mut v, mut a) = state_vectors();
    q.truncate(10);
    v.truncate(10);
    a.truncate(10);
    for (actual, expected) in actions(&set, &q, &v, &a)
        .iter()
        .zip(actions(&original, &q, &v, &a))
    {
        close_vectors(actual, &expected, 1e-13);
    }
    let (mut got, mut expected) = (vec![0.0; 9], vec![0.0; 9]);
    set.residual_parameter_vjp(&q, &v, &a, &v, &mut got)
        .unwrap();
    original
        .residual_parameter_vjp(&q, &v, &a, &v, &mut expected)
        .unwrap();
    close_vectors(&got, &expected, 1e-13);
    let (mut got, mut expected) = (vec![0.0; 10], vec![0.0; 10]);
    set.mass_preconditioner().apply(0, &a, &mut got);
    original.mass_preconditioner().apply(0, &a, &mut expected);
    close_vectors(&got, &expected, 1e-13);
    let probe = set.probe(&q, &v).unwrap();
    let one = original.probe(&q, &v).unwrap();
    close(probe.total_energy_j, one.total_energy_j, 1e-13);
    close(
        probe.plate_dissipated_power_w,
        one.plate_dissipated_power_w,
        1e-13,
    );
    close(probe.bodies[0].normal_force_n, one.normal_force_n, 1e-13);
    let step = method(10, 0.001);
    set.forcing(step.forcing_time(0.0).unwrap(), &mut got)
        .unwrap();
    original
        .forcing(step.forcing_time(0.0).unwrap(), &mut expected)
        .unwrap();
    close_vectors(&got, &expected, 1e-13);
    let mut combined = SecondOrderState::new(0.0, &q, &v, &a);
    let mut single = combined.clone();
    step.step(&mut combined, &set, &got).unwrap();
    step.step(&mut single, &original, &expected).unwrap();
    close_vectors(&combined.q, &single.q, 1e-11);
    close_vectors(&combined.v, &single.v, 1e-11);
    close_vectors(&combined.a, &single.a, 1e-9);
}

#[test]
fn contacts_exchange_momentum_through_one_plate_and_influence_each_other() {
    let (mesh, plate) = fixture(0.0);
    let load = driving_load(&mesh, &plate, 0.0).unwrap();
    let mut p = PARAMETERS;
    p[..5].copy_from_slice(&[1.0, 1.0, 0.0, 0.0, 0.0]);
    p[6] = 2000.0;
    let coupled = pair(&mesh, &plate, &load, p, false);
    let mut q = vec![0.0; 11];
    let mut v = q.clone();
    q[9] = 0.02;
    q[10] = 0.02;
    v[9] = 0.8;
    let probe = coupled.probe(&q, &v).unwrap();
    assert!(probe.bodies.iter().all(|body| body.normal_force_n > 0.0));
    let mut resistance = vec![0.0; 11];
    coupled.resistance(&q, &v, &mut resistance);
    close(resistance.iter().sum(), 0.0, 1e-12);
    close(
        dot(&resistance, &v),
        probe
            .bodies
            .iter()
            .map(|body| body.normal_force_n * body.relative_velocity_m_per_s)
            .sum(),
        1e-12,
    );
    let run = |p: [f64; 13]| {
        let model = pair(&mesh, &plate, &load, p, false);
        let mut state = consistent_state(&model, &q, &v);
        let initial_momentum = p[5] * v[9] + p[9] * v[10];
        for _ in 0..40 {
            method(11, 0.001)
                .step(&mut state, &model, &[0.0; 11])
                .unwrap();
            let momentum = p[5] * state.v[9]
                + p[9] * state.v[10]
                + [0, 3, 6]
                    .into_iter()
                    .map(|i| (2.0 + i as f64 / 10.0) * state.v[i])
                    .sum::<f64>();
            close(momentum, initial_momentum, 2e-9);
        }
        let probe = model.probe(&state.q, &state.v).unwrap();
        let plate_energy: f64 = (0..9)
            .map(|i| 0.5 * (2.0 + i as f64 / 10.0) * state.v[i].powi(2))
            .sum();
        close(probe.plate_kinetic_energy_j, plate_energy, 1e-12);
        close(
            probe.total_energy_j,
            probe.plate_kinetic_energy_j
                + probe.plate_elastic_energy_j
                + probe.indenter_kinetic_energy_j
                + probe.contact_energy_j,
            1e-12,
        );
        state
    };
    let touching = run(p);
    p[6] = 0.0; // Remove only the first body's contact, keeping the second body unchanged.
    let first_contact_absent = run(p);
    assert!(
        (touching.v[10] - first_contact_absent.v[10]).abs() > 1e-5,
        "the first body's impulse must change the second body's trajectory through the shared plate"
    );
}

#[test]
fn body_permutation_preserves_physics_and_permutates_contact_parameter_blocks() {
    let (mesh, plate) = fixture(4.0);
    let load = driving_load(&mesh, &plate, PARAMETERS[4]).unwrap();
    let ab = pair(&mesh, &plate, &load, PARAMETERS, false);
    let ba = pair(&mesh, &plate, &load, PARAMETERS, true);
    assert_eq!(ab.contact_parameter_range(0), Some(5..9));
    assert_eq!(ab.contact_parameter_range(1), Some(9..13));
    let (q, v, a) = state_vectors();
    let (pq, pv, pa) = (swap_bodies(&q), swap_bodies(&v), swap_bodies(&a));
    for (actual, permuted) in actions(&ab, &q, &v, &a)
        .iter()
        .zip(actions(&ba, &pq, &pv, &pa))
    {
        close_vectors(actual, &swap_bodies(&permuted), 1e-12);
    }
    let (first, second) = (ab.probe(&q, &v).unwrap(), ba.probe(&pq, &pv).unwrap());
    assert_eq!(first.bodies[0], second.bodies[1]);
    assert_eq!(first.bodies[1], second.bodies[0]);
    close(first.total_energy_j, second.total_energy_j, 1e-13);
    let (mut ga, mut gb) = (vec![0.0; 13], vec![0.0; 13]);
    ab.residual_parameter_vjp(&q, &v, &a, &v, &mut ga).unwrap();
    ba.residual_parameter_vjp(&pq, &pv, &pa, &pv, &mut gb)
        .unwrap();
    close_vectors(&ga[..5], &gb[..5], 1e-12);
    close_vectors(&ga[5..9], &gb[9..13], 1e-12);
    close_vectors(&ga[9..13], &gb[5..9], 1e-12);
    let mut applied = vec![f64::NAN; 13];
    let mut common = vec![f64::NAN; 5];
    ab.forcing_vjp(0.013, &v, &mut applied).unwrap();
    ab.plate().forcing_vjp(0.013, &v[..9], &mut common).unwrap();
    close_vectors(&applied[..5], &common, 1e-13);
    assert_eq!(&applied[5..], &[0.0; 8]);
}

#[test]
fn mixed_contact_actions_and_all_step_parameter_pullbacks_match_differences() {
    let (mesh, plate) = fixture(4.0);
    let load = driving_load(&mesh, &plate, PARAMETERS[4]).unwrap();
    let set = pair(&mesh, &plate, &load, PARAMETERS, false);
    let (q, v, a) = state_vectors();
    let d: Vec<_> = (0..11).map(|i| 0.04 * (i as f64 - 3.0)).collect();
    let seed: Vec<_> = (0..11).map(|i| 0.03 * (5.0 - i as f64)).collect();
    let mut mixed = q.clone();
    mixed[10] = -0.01;
    let probe = set.probe(&mixed, &v).unwrap();
    assert_eq!(probe.bodies[0].active_contact_points, 1);
    assert_eq!(probe.bodies[1].active_contact_points, 0);
    for velocity in [false, true] {
        let (mut action, mut transpose) = (vec![0.0; 11], vec![0.0; 11]);
        if velocity {
            set.resistance_velocity_apply(&mixed, &v, &d, &mut action);
            set.resistance_velocity_transpose_apply(&mixed, &v, &seed, &mut transpose)
                .unwrap();
        } else {
            set.resistance_displacement_apply(&mixed, &v, &d, &mut action);
            set.resistance_displacement_transpose_apply(&mixed, &v, &seed, &mut transpose)
                .unwrap();
        }
        close(dot(&seed, &action), dot(&d, &transpose), 1e-12);
        let (mut qp, mut qm, mut vp, mut vm) = (mixed.clone(), mixed.clone(), v.clone(), v.clone());
        let (plus, minus) = if velocity {
            (&mut vp, &mut vm)
        } else {
            (&mut qp, &mut qm)
        };
        for i in 0..11 {
            plus[i] += 1e-7 * d[i];
            minus[i] -= 1e-7 * d[i];
        }
        let (mut fp, mut fm) = (vec![0.0; 11], vec![0.0; 11]);
        set.resistance(&qp, &vp, &mut fp);
        set.resistance(&qm, &vm, &mut fm);
        for i in 0..11 {
            close(action[i], (fp[i] - fm[i]) / 2e-7, 3e-7);
        }
    }
    // Both contacts are active in this step, so every body's physical suffix
    // is exercised through the actual coupled Newton solve and transpose solve.
    let initial = SecondOrderState::new(0.0, &q, &v, &a);
    let step = method(11, 0.001);
    let time = step.forcing_time(0.0).unwrap();
    let mut force = vec![0.0; 11];
    set.forcing(time, &mut force).unwrap();
    let config = adjoint_config(11);
    let pulled = step
        .step_vjp(
            &initial,
            &set,
            &force,
            (&seed, &d, &a),
            &set.mass_preconditioner(),
            config,
            step.adjoint_workspace_components(13, config).unwrap(),
            &mut || false,
        )
        .unwrap();
    let mut gradient = pulled.parameters;
    let mut applied = vec![0.0; 13];
    set.forcing_vjp(time, &pulled.forcing, &mut applied)
        .unwrap();
    for (g, f) in gradient.iter_mut().zip(applied) {
        *g += f;
    }
    let loss = |p: [f64; 13]| {
        let load = driving_load(&mesh, &plate, p[4]).unwrap();
        let model = pair(&mesh, &plate, &load, p, false);
        let mut end = initial.clone();
        let mut force = vec![0.0; 11];
        model.forcing(time, &mut force).unwrap();
        step.step(&mut end, &model, &force).unwrap();
        dot(&seed, &end.q) + dot(&d, &end.v) + dot(&a, &end.a)
    };
    for i in 0..13 {
        let h = 1e-5 * PARAMETERS[i].abs().max(1e-3);
        let (mut plus, mut minus) = (PARAMETERS, PARAMETERS);
        plus[i] += h;
        minus[i] -= h;
        let fd = (loss(plus) - loss(minus)) / (2.0 * h);
        assert!(
            (gradient[i] - fd).abs() <= 3e-8 + 3e-5 * gradient[i].abs().max(fd.abs()),
            "parameter {i}: {} != {fd}",
            gradient[i]
        );
    }
}

#[test]
fn one_source_endpoint_refusal_rolls_back_all_bodies_and_existing_history() {
    let source = sourced::synthetic_source("multiple-impact source", 1e7, 900.0, 0.4).unwrap();
    let h = source.envelope().maximum_step_s;
    let (mesh, plate) = sourced::build_plate(&source).unwrap();
    let load = NoPlateLoad;
    let authored = body(
        &mesh,
        &plate,
        &load,
        [1.0, 1.0, 0.0, 0.0],
        [0.2, 700.0, 3.0, 0.005],
        [0.1, 0.1],
    )
    .with_contact_calibration()
    .unwrap();
    let bounded = sourced::build_impact(&mesh, &plate, source, &load).unwrap();
    let set = PlateImpactSet::new(vec![authored, bounded], budget()).unwrap();
    let n = set.dimension();
    let offset = set.plate_dofs();
    assert_eq!(set.contact_parameter_range(0), Some(4..8));
    assert_eq!(set.contact_parameter_range(1), Some(8..8));
    let zero = vec![0.0; n];
    let mut state = SecondOrderState::new(0.0, &zero, &zero, &zero);
    sourced::method(n, sourced::STEP_S)
        .step(&mut state, &set, &zero)
        .unwrap();
    assert_eq!(state.history.len(), 1);
    let before = state.clone();
    let mut stage_q = zero.clone();
    let mut stage_v = zero.clone();
    stage_q[offset + 1] = 4.5e-6;
    stage_v[offset + 1] = 0.09;
    set.validate_state(&stage_q, &stage_v).unwrap();
    let mut force = zero;
    force[offset] = 0.04; // This body's accepted state must also remain unchanged.
    force[offset + 1] = 1.8;
    let error = sourced::method(n, h)
        .step(&mut state, &set, &force)
        .unwrap_err();
    assert!(
        matches!(&error, TimeSolveError::ModelRefused { phase: "accepted endpoint", reason }
        if reason.contains("source point absolute rate")),
        "{error}"
    );
    assert_eq!(state, before);
}

#[test]
fn admission_rejects_distinct_owners_mismatched_plate_parameters_and_all_caps() {
    let (mesh, plate) = fixture(4.0);
    let (_, other_plate) = fixture(4.0);
    let load = driving_load(&mesh, &plate, PARAMETERS[4]).unwrap();
    let other_load = driving_load(&mesh, &plate, PARAMETERS[4]).unwrap();
    let common: [f64; 4] = PARAMETERS[..4].try_into().unwrap();
    let contact: [f64; 4] = PARAMETERS[5..9].try_into().unwrap();
    let make_body = |model, load, p| {
        body(&mesh, model, load, p, contact, POINTS[0])
            .with_contact_calibration()
            .unwrap()
    };
    assert!(matches!(
        PlateImpactSet::<DrivingLoad>::new(vec![], budget()),
        Err(PlateImpactError::Invalid(_))
    ));
    for bodies in [
        vec![
            make_body(&plate, &load, common),
            make_body(&other_plate, &load, common),
        ],
        vec![
            make_body(&plate, &load, common),
            make_body(&plate, &other_load, common),
        ],
    ] {
        assert!(matches!(
            PlateImpactSet::new(bodies, budget()),
            Err(PlateImpactError::Invalid(_))
        ));
    }
    for i in 0..4 {
        let mut different = common;
        different[i] += 0.01;
        assert!(
            PlateImpactSet::new(
                vec![
                    make_body(&plate, &load, common),
                    make_body(&plate, &load, different)
                ],
                budget(),
            )
            .is_err()
        );
    }
    for cap in 0..4 {
        let mut limit = budget();
        match cap {
            0 => limit.max_bodies = 1,
            1 => limit.max_free_dofs = 10,
            2 => limit.max_contact_points = 1,
            _ => limit.max_parameters = 12,
        }
        assert!(matches!(
            PlateImpactSet::new(
                vec![
                    make_body(&plate, &load, common),
                    make_body(&plate, &load, common)
                ],
                limit,
            ),
            Err(PlateImpactError::Budget(_))
        ));
    }
}
