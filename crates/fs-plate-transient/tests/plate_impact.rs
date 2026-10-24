//! Reciprocal contact ports and coupled structural solves on an analytic plate
//! pencil. The runnable plate_impact example additionally assembles DKT elements.
use fs_dcontact::Obstacle;
use fs_plate::{
    PlateMesh, PlateModel,
    loading::{PlateLoadBudget, PlatePatchBudget, PlatePatchPoint, PlatePatchStencil},
};
use fs_plate_transient::{
    NoPlateLoad, PlateDynamics, PlateDynamicsBudget, PlateDynamicsParameters, PlateImpact,
    PlateImpactBudget,
};
use fs_solver::{FlexiblePreconditioner, NewtonKrylovConfig};
use fs_sparse::Coo;
use fs_time::galpha::{
    ImplicitSolveConfig, OperatorGeneralizedAlpha, SecondOrderProblem, SecondOrderState,
    initialization::{InitialSolveConfig, initial_workspace_components, second_order_acceleration},
    second_order_adjoint::{SecondOrderAdjointConfig, SecondOrderVjp},
};

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * (1.0 + b.abs()),
        "{a:.16e} != {b:.16e}"
    );
}
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
fn patch(mesh: &PlateMesh, model: &PlateModel) -> PlatePatchStencil {
    PlatePatchStencil::locate(
        mesh,
        model,
        &[
            PlatePatchPoint {
                point: [0.2, 0.2],
                area_m2: 0.001,
            },
            PlatePatchPoint {
                point: [0.6, 0.2],
                area_m2: 0.003,
            },
        ],
        PlatePatchBudget {
            max_points: 2,
            location: PlateLoadBudget {
                max_nodes: 3,
                max_triangles: 1,
            },
        },
    )
    .unwrap()
}
fn law(stiffness: f64, alpha: f64, chi: f64) -> Obstacle {
    Obstacle::new(
        vec![1.0],
        1,
        1,
        vec![0.005],
        vec![2.0],
        stiffness,
        alpha,
        "explicit analytic contact fixture".into(),
    )
    .unwrap()
    .with_internal_loss(chi)
    .unwrap()
}
fn dynamics<'a>(model: &'a PlateModel, p: [f64; 4]) -> PlateDynamics<'a, NoPlateLoad> {
    PlateDynamics::new(
        model,
        PlateDynamicsParameters {
            stiffness_scale: p[0],
            mass_scale: p[1],
            mass_damping_per_s: p[2],
            stiffness_damping_s: p[3],
        },
        &NoPlateLoad,
        PlateDynamicsBudget {
            max_dofs: 9,
            max_nonzeros: model.k.nnz() + model.m.nnz(),
            max_full_dofs: 27,
            max_load_parameters: 0,
        },
        &mut || false,
    )
    .unwrap()
}
fn impact<'a>(
    mesh: &PlateMesh,
    model: &'a PlateModel,
    p: [f64; 4],
    law: Obstacle,
) -> PlateImpact<'a, NoPlateLoad> {
    PlateImpact::new(
        dynamics(model, p),
        patch(mesh, model),
        0.2,
        law,
        PlateImpactBudget {
            max_free_dofs: 10,
            max_contact_points: 1,
        },
    )
    .unwrap()
}
fn weights() -> Vec<f64> {
    // Independently integrated barycentric weights: .25*(.6,.2,.2)
    // + .75*(.2,.6,.2). Footprint area must not rescale the contact law.
    vec![0.3, 0.0, 0.0, 0.5, 0.0, 0.0, 0.2, 0.0, 0.0]
}
struct Identity;
impl FlexiblePreconditioner for Identity {
    fn apply(&self, _: usize, input: &[f64], output: &mut [f64]) {
        output.copy_from_slice(input);
    }
}
fn method(step: f64) -> OperatorGeneralizedAlpha {
    OperatorGeneralizedAlpha::new(
        10,
        step,
        1.0,
        ImplicitSolveConfig {
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1e-12,
                relative_tolerance: 1e-12,
                linear_restart: 10,
                max_linear_cycles: 4,
                forcing_minimum: 1e-13,
                forcing_maximum: 1e-10,
                ..NewtonKrylovConfig::default()
            },
            max_newton_iterations: 16,
        },
    )
}

#[test]
fn contact_uses_reciprocal_work_and_original_law_without_double_area() {
    let (mesh, model) = fixture(4.0);
    let p = [1.1, 0.9, 0.04, 0.02];
    let coupled = impact(&mesh, &model, p, law(700.0, 1.5, 3.0));
    let bare = dynamics(&model, p);
    let b = weights();
    for (closure, rate) in [(0.03, 0.2), (0.03, -1.0), (-0.01, 0.2)] {
        let mut q: Vec<_> = (0..9).map(|i| 0.001 * (i as f64 - 4.0)).collect();
        let mut v: Vec<_> = (0..9).map(|i| 0.01 * (i as f64 - 2.0)).collect();
        q.push(closure + dot(&b, &q));
        v.push(rate + dot(&b, &v));
        let penetration = (closure - 0.005_f64).max(0.0);
        let elastic = 2.0 * 700.0 * penetration.powf(1.5);
        let force = elastic * (1.0 + 3.0 * rate).max(0.0);
        let loss = elastic * (3.0 * rate).max(-1.0);
        let mut resistance = vec![0.0; 10];
        let mut base = vec![0.0; 9];
        coupled.resistance(&q, &v, &mut resistance);
        bare.resistance(&q[..9], &v[..9], &mut base);
        let contact: Vec<_> = resistance[..9]
            .iter()
            .zip(&base)
            .map(|(r, b)| r - b)
            .collect();
        for i in 0..9 {
            close(contact[i], -b[i] * force, 2e-12);
        }
        close(resistance[9], force, 2e-12);
        close(contact.iter().sum::<f64>() + resistance[9], 0.0, 2e-12);
        close(
            dot(&contact, &v[..9]) + resistance[9] * v[9],
            elastic * rate + loss * rate,
            2e-12,
        );
        let probe = coupled.probe(&q, &v).unwrap();
        close(probe.closure_m, closure, 1e-14);
        close(probe.relative_velocity_m_per_s, rate, 1e-14);
        close(probe.normal_force_n, force, 2e-12);
        close(
            probe.contact_energy_j,
            2.0 * 700.0 / 2.5 * penetration.powf(2.5),
            2e-12,
        );
        close(probe.contact_dissipated_power_w, loss * rate, 2e-12);
        assert!(probe.normal_force_n >= 0.0);
        assert!(probe.contact_dissipated_power_w >= 0.0);
    }
}

#[test]
fn contact_partials_and_transposes_match_finite_differences_on_each_smooth_branch() {
    let (mesh, model) = fixture(4.0);
    let coupled = impact(&mesh, &model, [1.1, 0.9, 0.04, 0.02], law(700.0, 1.5, 3.0));
    let direction: Vec<_> = (0..10).map(|i| 0.08 * (i as f64 - 3.0)).collect();
    let seed: Vec<_> = (0..10).map(|i| 0.07 * (5.0 - i as f64)).collect();
    for (closure, rate) in [(0.03, 0.2), (0.03, -1.0), (-0.01, 0.2)] {
        let mut q = vec![0.0; 10];
        q[9] = closure;
        let mut v = vec![0.0; 10];
        v[9] = rate;
        for velocity_partial in [false, true] {
            let mut action = vec![0.0; 10];
            let mut transpose = vec![0.0; 10];
            if velocity_partial {
                coupled.resistance_velocity_apply(&q, &v, &direction, &mut action);
                coupled
                    .resistance_velocity_transpose_apply(&q, &v, &seed, &mut transpose)
                    .unwrap();
            } else {
                coupled.resistance_displacement_apply(&q, &v, &direction, &mut action);
                coupled
                    .resistance_displacement_transpose_apply(&q, &v, &seed, &mut transpose)
                    .unwrap();
            }
            close(dot(&seed, &action), dot(&direction, &transpose), 2e-12);
            let epsilon = 1e-7;
            let (mut qp, mut qm, mut vp, mut vm) = (q.clone(), q.clone(), v.clone(), v.clone());
            let (plus, minus) = if velocity_partial {
                (&mut vp, &mut vm)
            } else {
                (&mut qp, &mut qm)
            };
            for i in 0..10 {
                plus[i] += epsilon * direction[i];
                minus[i] -= epsilon * direction[i];
            }
            let (mut fp, mut fm) = (vec![0.0; 10], vec![0.0; 10]);
            coupled.resistance(&qp, &vp, &mut fp);
            coupled.resistance(&qm, &vm, &mut fm);
            for i in 0..10 {
                close(action[i], (fp[i] - fm[i]) / (2.0 * epsilon), 2e-7);
            }
        }
    }
}

#[test]
fn free_plate_impact_transfers_momentum_and_rebounds_after_separation() {
    let (mesh, model) = fixture(0.0);
    let coupled = impact(&mesh, &model, [1.0, 1.0, 0.0, 0.0], law(2000.0, 1.0, 0.2));
    let q = vec![0.0; 10];
    let mut v = q.clone();
    v[9] = 0.8;
    let force = vec![0.0; 10];
    let config = InitialSolveConfig {
        restart: 10,
        max_cycles: 4,
        tolerance: 1e-12,
    };
    let initial = second_order_acceleration(
        &coupled,
        &q,
        &v,
        &force,
        &coupled.mass_preconditioner(),
        config,
        initial_workspace_components(10, coupled.parameter_count(), config).unwrap(),
        &mut || false,
    )
    .unwrap();
    let mut state = SecondOrderState::new(0.0, &q, &v, &initial.value);
    let initial_energy = coupled.probe(&q, &v).unwrap().total_energy_j;
    let mut touched = false;
    let mut peak_force = 0.0_f64;
    let method = method(0.001);
    for _ in 0..180 {
        method.step(&mut state, &coupled, &force).unwrap();
        let probe = coupled.probe(&state.q, &state.v).unwrap();
        touched |= probe.active_contact_points > 0;
        peak_force = peak_force.max(probe.normal_force_n);
        let momentum = 0.2 * state.v[9]
            + [0, 3, 6]
                .into_iter()
                .map(|i| (2.0 + i as f64 / 10.0) * state.v[i])
                .sum::<f64>();
        close(momentum, 0.2 * 0.8, 2e-9);
    }
    let final_probe = coupled.probe(&state.q, &state.v).unwrap();
    assert!(touched && peak_force > 1.0);
    assert!(state.v[9] < 0.0, "finite-mass body must rebound");
    assert!(
        final_probe.patch.velocity > 0.0,
        "plate must receive momentum"
    );
    assert_eq!(final_probe.active_contact_points, 0);
    assert!(final_probe.total_energy_j < initial_energy);
}

#[test]
fn coupled_step_adjoint_includes_plate_parameters_and_initial_contact_motion() {
    let (mesh, model) = fixture(4.0);
    let p = [1.1, 0.9, 0.04, 0.02];
    let coupled = impact(&mesh, &model, p, law(700.0, 1.5, 3.0));
    let mut q: Vec<_> = (0..10).map(|i| 0.001 * i as f64).collect();
    q[9] = 0.04;
    let mut v: Vec<_> = (0..10).map(|i| 0.002 * i as f64).collect();
    v[9] = 0.15;
    let initial = SecondOrderState::new(0.0, &q, &v, &[0.0; 10]);
    let forcing = [0.0; 10];
    let (q_seed, v_seed, a_seed) = ([0.2; 10], [-0.3; 10], [0.04; 10]);
    let method = method(0.01);
    let adjoint = SecondOrderAdjointConfig {
        restart: 10,
        max_cycles: 4,
        tolerance: 1e-12,
    };
    let gradient = method
        .step_vjp(
            &initial,
            &coupled,
            &forcing,
            (&q_seed, &v_seed, &a_seed),
            &Identity,
            adjoint,
            method.adjoint_workspace_components(4, adjoint).unwrap(),
            &mut || false,
        )
        .unwrap();
    let evaluate = |parameters: [f64; 4], initial: &SecondOrderState| {
        let coupled = impact(&mesh, &model, parameters, law(700.0, 1.5, 3.0));
        let mut end = initial.clone();
        method.step(&mut end, &coupled, &forcing).unwrap();
        dot(&q_seed, &end.q) + dot(&v_seed, &end.v) + dot(&a_seed, &end.a)
    };
    for i in 0..4 {
        let epsilon = 1e-5;
        let (mut plus, mut minus) = (p, p);
        plus[i] += epsilon;
        minus[i] -= epsilon;
        let fd = (evaluate(plus, &initial) - evaluate(minus, &initial)) / (2.0 * epsilon);
        close(gradient.parameters[i], fd, 3e-6);
    }
    let direction: Vec<_> = (0..10).map(|i| 0.03 * (i as f64 - 3.0)).collect();
    for velocity_seed in [false, true] {
        let epsilon = 1e-7;
        let (mut plus, mut minus) = (initial.clone(), initial.clone());
        let (a, b) = if velocity_seed {
            (&mut plus.v, &mut minus.v)
        } else {
            (&mut plus.q, &mut minus.q)
        };
        for i in 0..10 {
            a[i] += epsilon * direction[i];
            b[i] -= epsilon * direction[i];
        }
        let fd = (evaluate(p, &plus) - evaluate(p, &minus)) / (2.0 * epsilon);
        let bar = if velocity_seed {
            &gradient.initial_v
        } else {
            &gradient.initial_q
        };
        close(dot(bar, &direction), fd, 3e-6);
    }
}

#[test]
fn patch_binding_rejects_changed_support_indexing_and_bounded_admission() {
    let (mesh, mut model) = fixture(4.0);
    let original = patch(&mesh, &model);
    original.validate_model(&model).unwrap();
    model.dof_map.swap(0, 3); // Same size, different physical node assignment.
    assert!(original.validate_model(&model).is_err());
    assert!(
        PlateImpact::new(
            dynamics(&model, [1.0, 1.0, 0.0, 0.0]),
            original,
            0.2,
            law(700.0, 1.5, 3.0),
            PlateImpactBudget {
                max_free_dofs: 10,
                max_contact_points: 1
            }
        )
        .is_err()
    );
    for budget in [
        PlateImpactBudget {
            max_free_dofs: 9,
            max_contact_points: 1,
        },
        PlateImpactBudget {
            max_free_dofs: 10,
            max_contact_points: 0,
        },
    ] {
        assert!(
            PlateImpact::new(
                dynamics(&model, [1.0, 1.0, 0.0, 0.0]),
                patch(&mesh, &model),
                0.2,
                law(700.0, 1.5, 3.0),
                budget
            )
            .is_err()
        );
    }
    for mass in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            PlateImpact::new(
                dynamics(&model, [1.0, 1.0, 0.0, 0.0]),
                patch(&mesh, &model),
                mass,
                law(700.0, 1.5, 3.0),
                PlateImpactBudget {
                    max_free_dofs: 10,
                    max_contact_points: 1
                }
            )
            .is_err()
        );
    }
}
