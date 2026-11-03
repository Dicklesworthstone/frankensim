//! G0/G1/G3/G4: reciprocal impulses on two independent native plate trajectories.

use super::*;
use fs_contact::multicontact::{MultiImpactConfig, MultiImpactDeterminacy};
use fs_couple::plate_impact::multiple::{
    PlateImpactSite, PlateRigidImpactSetBudget, PlateRigidImpactSetControls,
};
use fs_couple::plate_impact::plates::{
    PlateImpactSurface, PlatePlateImpact, PlatePlateImpactBudget, PlatePlateImpactReport,
};
use fs_plate_transient::PlateLoad;

fn pair_controls() -> PlateRigidImpactSetControls {
    let single = controls(0.0);
    let mut contact = MultiImpactConfig::default();
    contact.solver.gradient_tolerance = 1e-12;
    let mut impulse = single.impulse;
    impulse.mass_solve.restart = 18;
    PlateRigidImpactSetControls {
        impulse,
        contact,
        point_tolerance_m: single.point_tolerance_m,
        speed_tolerance_m_s: single.speed_tolerance_m_s,
        linear_impulse_tolerance_ns: single.linear_impulse_tolerance_ns,
        angular_impulse_tolerance_nms: single.angular_impulse_tolerance_nms,
        energy_tolerance_j: single.energy_tolerance_j,
    }
}
fn pair_budget() -> PlatePlateImpactBudget {
    let owner = PlateRigidImpactSetBudget {
        plate: PlateRigidImpactBudget {
            location: PlateLoadBudget {
                max_nodes: 6,
                max_triangles: 4,
            },
            max_full_mass_nonzeros: 1000,
            max_workspace_components: 20000,
        },
        max_contacts: 8,
        max_mobility_products: 10000,
    };
    PlatePlateImpactBudget {
        positive: owner,
        negative: owner,
        max_workspace_components: 100000,
    }
}
fn sites(points: &[[f64; 2]], restitution: f64) -> Vec<PlateImpactSite> {
    points
        .iter()
        .map(|&point_xy| PlateImpactSite {
            point_xy,
            activation_tolerance_m: 1e-9,
            restitution,
        })
        .collect()
}
fn surface<'a>(
    mesh: &'a PlateMesh,
    model: &'a PlateModel,
    mass: &'a Csr,
    scale: f64,
    damping: f64,
    height: f64,
) -> PlateImpactSurface<'a, NoPlateLoad> {
    PlateImpactSurface {
        dynamics: PlateDynamics::new(
            model,
            PlateDynamicsParameters {
                stiffness_scale: 1.0,
                mass_scale: scale,
                mass_damping_per_s: damping,
                stiffness_damping_s: 0.0,
            },
            &NoPlateLoad,
            PlateDynamicsBudget {
                max_dofs: 18,
                max_nonzeros: 1000,
                max_full_dofs: 18,
                max_load_parameters: 0,
            },
            &mut || false,
        )
        .unwrap(),
        mesh,
        full_mass: mass,
        reference_height_m: height,
    }
}
fn state(n: usize) -> SecondOrderState {
    SecondOrderState::new(0.7, &vec![0.0; n], &vec![0.0; n], &vec![0.0; n])
}
fn audit(report: &PlatePlateImpactReport) {
    let balance = report.balance.as_ref().unwrap();
    close_vec(balance.linear_balance_residual, Vec3::ZERO);
    close_vec(balance.angular_balance_residual, Vec3::ZERO);
    close(balance.work_energy_residual_j, 0.0);
    close(balance.energy_residual_j, 0.0);
    close(
        balance.kinetic_energy_change_j,
        -report.contact.dissipated_energy,
    );
    for (&a, &b) in report
        .relative_normal_speeds_after_m_s
        .iter()
        .zip(&report.contact.relative_normal_speeds_after)
    {
        close(a, b);
    }
}

// Independent exact P1 area integrals on a four-triangle refinement. No cached
// production point response or structural mass solver constructs this oracle.
fn fine_fixture(shift: [f64; 2]) -> (PlateMesh, PlateModel, Csr) {
    let mut mesh = PlateMesh {
        nodes: vec![
            (0.0, 0.0),
            (1.0, 0.0),
            (0.0, 1.0),
            (0.5, 0.0),
            (0.5, 0.5),
            (0.0, 0.5),
        ],
        tris: vec![[0, 3, 5], [3, 1, 4], [5, 4, 2], [3, 4, 5]],
    };
    for (x, y) in &mut mesh.nodes {
        *x += shift[0];
        *y += shift[1];
    }
    let mut mass = Coo::new(18, 18);
    for &[a, b, c] in &mesh.tris {
        let [p, q, r] = [mesh.nodes[a], mesh.nodes[b], mesh.nodes[c]];
        let area = 0.5 * ((q.0 - p.0) * (r.1 - p.1) - (q.1 - p.1) * (r.0 - p.0));
        for i in [a, b, c] {
            for j in [a, b, c] {
                for component in 0..3 {
                    let rotary = if component == 0 { 1.0 } else { 1.0 / 12.0 };
                    mass.push(
                        3 * i + component,
                        3 * j + component,
                        0.5 * area * rotary * if i == j { 2.0 } else { 1.0 },
                    );
                }
            }
        }
    }
    let full = mass.assemble();
    let model = PlateModel {
        m: full.clone(),
        k: Coo::new(18, 18).assemble(),
        dof_map: (0..18).map(Some).collect(),
        free: 18,
    };
    (mesh, model, full)
}
fn coarse_row(point: [f64; 2]) -> Vec<f64> {
    vec![
        1.0 - point[0] - point[1],
        0.0,
        0.0,
        point[0],
        0.0,
        0.0,
        point[1],
        0.0,
        0.0,
    ]
}
fn fine_rows() -> Vec<Vec<f64>> {
    [
        [(0, 0.25), (3, 0.5), (5, 0.25)],
        [(3, 0.25), (1, 0.5), (4, 0.25)],
        [(5, 0.25), (4, 0.25), (2, 0.5)],
    ]
    .iter()
    .map(|entries| {
        let mut row = vec![0.0; 18];
        for &(node, weight) in entries {
            row[3 * node] = weight;
        }
        row
    })
    .collect()
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

#[test]
fn independently_meshed_plate_pair_matches_dense_mass_response_and_continues_both_owners() {
    let points = [[0.25, 0.125], [0.75, 0.125], [0.125, 0.75]];
    let positive_rows: Vec<_> = points.map(coarse_row).into_iter().collect();
    let negative_rows = fine_rows();
    for restitution in [0.0, 0.4, 1.0] {
        let mut reference = None;
        for scene in 0..2 {
            let shift = if scene == 0 { [0.0, 0.0] } else { [3.5, -2.25] };
            let order = if scene == 0 { [0, 1, 2] } else { [2, 0, 1] };
            let (mut mesh, model, full) = fixture(false);
            mesh.nodes = vec![
                (shift[0], shift[1]),
                (1.0 + shift[0], shift[1]),
                (shift[0], 1.0 + shift[1]),
            ];
            let (other_mesh, other_model, other_full) = fine_fixture(shift);
            let positive_dynamics = surface(&mesh, &model, &full, 1.3, 0.07, 0.0).dynamics;
            let negative_dynamics =
                surface(&other_mesh, &other_model, &other_full, 0.8, 0.19, 0.0).dynamics;
            let ordered_points = order.map(|i| [points[i][0] + shift[0], points[i][1] + shift[1]]);
            let prepared = PlatePlateImpact::prepare(
                surface(&mesh, &model, &full, 1.3, 0.07, 0.0),
                surface(&other_mesh, &other_model, &other_full, 0.8, 0.19, 0.0),
                &sites(&ordered_points, restitution),
                pair_controls(),
                pair_budget(),
                &mut || false,
            )
            .unwrap();
            let positive_mass: Vec<_> = full.to_dense().iter().map(|v| 1.3 * v).collect();
            let negative_mass: Vec<_> = other_full.to_dense().iter().map(|v| 0.8 * v).collect();
            let positive_factor = lu(&positive_mass, 9).unwrap();
            let negative_factor = lu(&negative_mass, 18).unwrap();
            let mut positive_response = positive_rows.clone();
            let mut negative_response = negative_rows.clone();
            for row in &mut positive_response {
                positive_factor.solve(row);
            }
            for row in &mut negative_response {
                negative_factor.solve(row);
            }
            let mut matrix = vec![0.0; 9];
            for i in 0..3 {
                for j in 0..3 {
                    matrix[3 * i + j] = dot(&positive_rows[i], &positive_response[j])
                        + dot(&negative_rows[i], &negative_response[j]);
                }
            }
            for i in 0..3 {
                for j in 0..3 {
                    close(
                        prepared.contact_mobility()[3 * i + j],
                        matrix[3 * order[i] + order[j]],
                    );
                }
            }
            // Construct a strictly compressive independent simultaneous impulse.
            let impulse = [0.2; 3];
            let before: Vec<_> = (0..3)
                .map(|i| -dot(&matrix[3 * i..3 * i + 3], &impulse) / (1.0 + restitution))
                .collect();
            assert!(before.iter().all(|v| *v < 0.0));
            let interpolation: Vec<_> = positive_rows
                .iter()
                .flat_map(|r| [r[0], r[3], r[6]])
                .collect();
            let mut nodal = before.clone();
            lu(&interpolation, 3).unwrap().solve(&mut nodal);
            let mut positive = state(9);
            let mut negative = state(18);
            for i in 0..3 {
                positive.v[3 * i] = nodal[i];
                positive.a[3 * i] = -0.07 * nodal[i];
            }
            let before_positive = positive.clone();
            let before_negative = negative.clone();
            let report = prepared
                .apply(&mut positive, &mut negative, &mut || false)
                .unwrap();
            audit(&report);
            for &j in &report.contact.normal_impulses {
                close(j, 0.2);
            }
            for r in 0..9 {
                close(
                    positive.v[r],
                    before_positive.v[r]
                        + (0..3).map(|i| 0.2 * positive_response[i][r]).sum::<f64>(),
                );
                close(positive.a[r], -0.07 * positive.v[r]);
            }
            for r in 0..18 {
                close(
                    negative.v[r],
                    -(0..3).map(|i| 0.2 * negative_response[i][r]).sum::<f64>(),
                );
                close(negative.a[r], -0.19 * negative.v[r]);
            }
            close(
                report.contact.dissipated_energy,
                -0.5 * (1.0 - restitution) * dot(&impulse, &before),
            );
            assert_eq!(
                (&positive.q, positive.t, positive.steps, &positive.history),
                (
                    &before_positive.q,
                    before_positive.t,
                    before_positive.steps,
                    &before_positive.history
                )
            );
            assert_eq!(
                (&negative.q, negative.t, negative.steps, &negative.history),
                (
                    &before_negative.q,
                    before_negative.t,
                    before_negative.steps,
                    &before_negative.history
                )
            );
            if let Some((pv, nv)) = &reference {
                for (&a, &b) in positive.v.iter().zip(pv) {
                    close(a, b);
                }
                for (&a, &b) in negative.v.iter().zip(nv) {
                    close(a, b);
                }
            } else {
                reference = Some((positive.v.clone(), negative.v.clone()));
            }
            for (state, dynamics, dimension) in [
                (&mut positive, positive_dynamics, 9),
                (&mut negative, negative_dynamics, 18),
            ] {
                let energy_before = dot(&state.v, &state.v);
                OperatorGeneralizedAlpha::new(dimension, 1e-3, 1.0, ImplicitSolveConfig::default())
                    .step(state, &dynamics, &vec![0.0; dimension])
                    .unwrap();
                assert_eq!(state.steps, 1);
                assert_eq!(state.history.len(), 1);
                assert!(state.t > 0.7 && dot(&state.v, &state.v) < energy_before);
            }
        }
    }
}

#[test]
fn initially_separating_and_inactive_rows_keep_the_shared_response_and_defeat_pair_sweeps() {
    let (mut mesh, model, full) = fixture(false);
    mesh.nodes = vec![(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)];
    let points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
    let make = |points: &[[f64; 2]]| {
        PlatePlateImpact::prepare(
            surface(&mesh, &model, &full, 1.0, 0.0, 0.0),
            surface(&mesh, &model, &full, 1.0, 0.0, 0.0),
            &sites(points, 0.0),
            pair_controls(),
            pair_budget(),
            &mut || false,
        )
        .unwrap()
    };
    let prepared = make(&points);
    let mut positive = state(9);
    positive.q[6] = 0.01;
    positive.v[0] = -2.0;
    positive.v[3] = 0.25;
    let mut negative = state(9);
    let mut sequential = (positive.clone(), negative.clone());
    let report = prepared
        .apply(&mut positive, &mut negative, &mut || false)
        .unwrap();
    audit(&report);
    // The active two-row system is [[6,-2],[-2,6]].
    close(report.contact.normal_impulses[0], 11.5 / 32.0);
    close(report.contact.normal_impulses[1], 2.5 / 32.0);
    close(report.contact.normal_impulses[2], 0.0);
    close(report.relative_normal_speeds_after_m_s[0], 0.0);
    close(report.relative_normal_speeds_after_m_s[1], 0.0);
    close(report.relative_normal_speeds_after_m_s[2], -0.875);
    for point in &points[..2] {
        make(&[*point])
            .apply(&mut sequential.0, &mut sequential.1, &mut || false)
            .unwrap();
    }
    assert!((sequential.0.v[0] - positive.v[0]).abs() > 1e-3);
    assert!(sequential.0.v[0] - sequential.1.v[0] < -1e-3);
}

#[test]
fn original_support_cross_inertia_and_redundant_contacts_preserve_total_motion() {
    let (mesh, model, full) = fixture(false);
    let (other_mesh, other_model, original_full) = fixture(true);
    let mut modified = Coo::new(9, 9);
    for i in 0..9 {
        let (cols, values) = original_full.row(i);
        for (&j, &v) in cols.iter().zip(values) {
            modified.push(i, j, v);
        }
    }
    // A complete fixed-slope/free-translation inertia block. Its transpose is
    // retained in the original full mass and vanishes from the reduced block.
    modified.push(1, 3, 0.01);
    modified.push(3, 1, 0.01);
    let other_full = modified.assemble();
    let mut reference = None;
    for duplicate in [false, true] {
        let points = if duplicate {
            vec![[1.2, -0.7], [1.2, -0.7], [2.2, -0.7]]
        } else {
            vec![[1.2, -0.7], [2.2, -0.7]]
        };
        let prepared = PlatePlateImpact::prepare(
            surface(&mesh, &model, &full, 1.0, 0.0, 0.0),
            surface(&other_mesh, &other_model, &other_full, 1.7, 0.0, 0.0),
            &sites(&points, 0.3),
            pair_controls(),
            pair_budget(),
            &mut || false,
        )
        .unwrap();
        let mut positive = state(9);
        positive.v[0] = -1.0;
        positive.v[3] = -0.7;
        let mut negative = state(6);
        let report = prepared
            .apply(&mut positive, &mut negative, &mut || false)
            .unwrap();
        audit(&report);
        if duplicate {
            assert!(matches!(
                report.contact.determinacy,
                MultiImpactDeterminacy::UniqueVelocity(_)
            ));
        }
        let balance = report.balance.as_ref().unwrap();
        assert!(balance.support_impulses[0].is_empty());
        let mut full_change = [0.0; 9];
        for (i, &r) in other_model.dof_map.iter().enumerate() {
            if let Some(r) = r {
                full_change[i] = negative.v[r];
            }
        }
        let dense = other_full.to_dense();
        let mut cross_inertia = 0.0_f64;
        for &(i, reaction) in &balance.support_impulses[1] {
            let inertia: f64 = (0..9)
                .map(|j| 1.7 * dense[9 * i + j] * full_change[j])
                .sum();
            let applied: f64 = points
                .iter()
                .zip(&report.contact.normal_impulses)
                .map(|(&p, &j)| -coarse_row([p[0] - 1.2, p[1] + 0.7])[i] * j)
                .sum();
            close(reaction, inertia - applied);
            cross_inertia = cross_inertia.max(inertia.abs());
            if i == 1 {
                assert!(reaction.abs() > 1e-4);
            }
        }
        assert!(cross_inertia > 1e-3);
        if let Some((pv, nv, loss)) = &reference {
            for (&a, &b) in positive.v.iter().zip(pv) {
                close(a, b);
            }
            for (&a, &b) in negative.v.iter().zip(nv) {
                close(a, b);
            }
            close(report.contact.dissipated_energy, *loss);
        } else {
            reference = Some((
                positive.v.clone(),
                negative.v.clone(),
                report.contact.dissipated_energy,
            ));
        }
    }
}

#[test]
fn collective_work_retains_small_loss_during_large_common_translation() {
    let (mut mesh, model, full) = fixture(false);
    mesh.nodes = vec![(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)];
    let prepared = PlatePlateImpact::prepare(
        surface(&mesh, &model, &full, 1.0, 0.0, 1e12),
        surface(&mesh, &model, &full, 1.0, 0.0, 1e12),
        &sites(&[[0.0, 0.0]], 0.0),
        pair_controls(),
        pair_budget(),
        &mut || false,
    )
    .unwrap();
    let common = 67_108_864.0;
    let unit = 1.0 / 1024.0;
    let mut positive = state(9);
    let mut negative = state(9);
    for i in [0, 3, 6] {
        positive.v[i] = common;
        negative.v[i] = common;
    }
    positive.v[0] -= 6.0 * unit;
    let before = (positive.clone(), negative.clone());
    let report = prepared
        .apply(&mut positive, &mut negative, &mut || false)
        .unwrap();
    audit(&report);
    close(report.contact.normal_impulses[0], unit);
    close(report.contact.dissipated_energy, 3.0 * unit * unit);
    let energy = |v: &[f64]| {
        let mut mv = vec![0.0; 9];
        full.spmv(v, &mut mv);
        0.5 * dot(v, &mv)
    };
    let naive =
        energy(&positive.v) + energy(&negative.v) - energy(&before.0.v) - energy(&before.1.v);
    assert!(
        (naive + report.contact.dissipated_energy).abs() > 0.5 * report.contact.dissipated_energy
    );
}

#[test]
fn pair_limits_synchronized_times_noops_and_late_failures_publish_neither_state() {
    let (mesh, model, full) = fixture(false);
    let point = [[2.2, -0.7]];
    let make = |budget| {
        PlatePlateImpact::prepare(
            surface(&mesh, &model, &full, 1.0, 0.0, 0.0),
            surface(&mesh, &model, &full, 1.0, 0.0, 0.0),
            &sites(&point, 0.4),
            pair_controls(),
            budget,
            &mut || false,
        )
    };
    let mut too_small = pair_budget();
    too_small.max_workspace_components = 1;
    assert!(matches!(
        make(too_small),
        Err(PlateRigidImpactError::Budget(_))
    ));
    let mut overflow = pair_budget();
    overflow.positive.plate.max_workspace_components = usize::MAX;
    assert!(matches!(
        make(overflow),
        Err(PlateRigidImpactError::Budget(_))
    ));
    let mut no_contacts = pair_budget();
    no_contacts.negative.max_contacts = 0;
    assert!(matches!(
        make(no_contacts),
        Err(PlateRigidImpactError::Budget(_))
    ));
    let prepared = make(pair_budget()).unwrap();
    let mut positive = state(9);
    positive.v[3] = 1.0;
    let mut negative = state(9);
    let unchanged = (positive.clone(), negative.clone());
    let zero = prepared
        .apply(&mut positive, &mut negative, &mut || false)
        .unwrap();
    assert!(zero.plates.iter().all(Option::is_none) && zero.balance.is_none());
    assert_eq!((positive.clone(), negative.clone()), unchanged);
    negative.t += 0.1;
    let unchanged = (positive.clone(), negative.clone());
    assert!(matches!(
        prepared.apply(&mut positive, &mut negative, &mut || false),
        Err(PlateRigidImpactError::Invalid(_))
    ));
    assert_eq!((positive, negative), unchanged);
    let mut initial = (state(9), state(9));
    initial.0.v[3] = -1.0;
    let mut accepted = initial.clone();
    let mut polls = 0;
    prepared
        .apply(&mut accepted.0, &mut accepted.1, &mut || {
            polls += 1;
            false
        })
        .unwrap();
    for limit in [0, polls / 2, polls - 1] {
        let mut candidate = initial.clone();
        let mut count = 0;
        let outcome = prepared.apply(&mut candidate.0, &mut candidate.1, &mut || {
            let stop = count == limit;
            count += 1;
            stop
        });
        assert!(outcome.is_err());
        assert_eq!(candidate, initial);
    }
    let mut retry = initial.clone();
    prepared
        .apply(&mut retry.0, &mut retry.1, &mut || false)
        .unwrap();
    assert_eq!(accepted, retry);
    // A small interpolated closing speed can coexist with nodal backgrounds
    // whose ulp exceeds the needed impulse update. The native owner must refuse.
    let prepared = PlatePlateImpact::prepare(
        surface(&mesh, &model, &full, 1.0, 0.0, 0.0),
        surface(&mesh, &model, &full, 1.0, 0.0, 0.0),
        &sites(&[[2.2 - 1e-4, -0.7]], 0.0),
        pair_controls(),
        pair_budget(),
        &mut || false,
    )
    .unwrap();
    let common = 1_152_921_504_606_846_976.0;
    let mut unresolved = (state(9), state(9));
    for i in [0, 3, 6] {
        unresolved.0.v[i] = common;
        unresolved.1.v[i] = common;
    }
    unresolved.0.v[0] -= 1024.0;
    let before = unresolved.clone();
    assert!(
        prepared
            .apply(&mut unresolved.0, &mut unresolved.1, &mut || false)
            .is_err()
    );
    assert_eq!(unresolved, before);
}

#[test]
fn distinct_native_load_types_restart_acceleration_and_second_owner_refusal_is_atomic() {
    use std::cell::Cell;
    struct HeldLoad {
        refuse: Cell<bool>,
    }
    impl PlateLoad for HeldLoad {
        fn parameter_count(&self) -> usize {
            0
        }
        fn forcing(&self, _: f64, output: &mut [f64]) -> Result<(), String> {
            if self.refuse.get() {
                return Err("second original load refused".into());
            }
            output.fill(0.0);
            output[0] = 1.5;
            Ok(())
        }
        fn forcing_vjp(&self, _: f64, _: &[f64], output: &mut [f64]) -> Result<(), String> {
            output.fill(0.0);
            Ok(())
        }
    }
    let (mesh, model, full) = fixture(false);
    let load = HeldLoad {
        refuse: Cell::new(true),
    };
    let negative_dynamics = PlateDynamics::new(
        &model,
        PlateDynamicsParameters {
            stiffness_scale: 1.0,
            mass_scale: 1.7,
            mass_damping_per_s: 0.2,
            stiffness_damping_s: 0.0,
        },
        &load,
        PlateDynamicsBudget {
            max_dofs: 9,
            max_nonzeros: 100,
            max_full_dofs: 9,
            max_load_parameters: 0,
        },
        &mut || false,
    )
    .unwrap();
    let prepared = PlatePlateImpact::prepare(
        surface(&mesh, &model, &full, 1.0, 0.1, 0.0),
        PlateImpactSurface {
            dynamics: negative_dynamics,
            mesh: &mesh,
            full_mass: &full,
            reference_height_m: 0.0,
        },
        &sites(&[[1.2, -0.7]], 0.3),
        pair_controls(),
        pair_budget(),
        &mut || false,
    )
    .unwrap();
    let mut positive = state(9);
    positive.v[0] = -1.0;
    let mut negative = state(9);
    let before = (positive.clone(), negative.clone());
    assert!(matches!(
        prepared.apply(&mut positive, &mut negative, &mut || false),
        Err(PlateRigidImpactError::Load(_))
    ));
    assert_eq!((positive.clone(), negative.clone()), before);
    load.refuse.set(false);
    let report = prepared
        .apply(&mut positive, &mut negative, &mut || false)
        .unwrap();
    audit(&report);
    let dense: Vec<_> = full.to_dense().iter().map(|v| 1.7 * v).collect();
    let mut acceleration = vec![0.0; 9];
    acceleration[0] = 1.5;
    lu(&dense, 9).unwrap().solve(&mut acceleration);
    for i in 0..9 {
        close(negative.a[i], acceleration[i] - 0.2 * negative.v[i]);
        close(positive.a[i], -0.1 * positive.v[i]);
    }
    let mut replay = before;
    prepared
        .apply(&mut replay.0, &mut replay.1, &mut || false)
        .unwrap();
    assert_eq!((positive, negative), replay);
}
