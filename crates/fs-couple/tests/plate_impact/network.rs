//! G0/G1/G3/G4: simultaneous momentum exchange between original plate owners.

use super::plates::{
    coarse_row, dot, fine_fixture, fine_rows, pair_budget, pair_controls, state, surface,
};
use super::*;
use fs_contact::multicontact::{MultiImpactDeterminacy, MultiImpactError};
use fs_couple::plate_impact::multiple::PlateImpactSite;
use fs_couple::plate_impact::network::{
    PlateImpactNetwork, PlateImpactNetworkBudget, PlateImpactNetworkReport, PlateNetworkContact,
};
use fs_couple::plate_impact::plates::{PlateImpactSurface, PlatePlateImpact};
use fs_plate_transient::PlateLoad;

fn network_budget() -> PlateImpactNetworkBudget {
    PlateImpactNetworkBudget {
        max_plates: 4,
        max_contacts: 8,
        max_total_full_dofs: 72,
        max_workspace_components: 200000,
    }
}
fn edge(
    positive_plate: usize,
    negative_plate: usize,
    point_xy: [f64; 2],
    restitution: f64,
) -> PlateNetworkContact {
    PlateNetworkContact {
        positive_plate,
        negative_plate,
        site: PlateImpactSite {
            point_xy,
            restitution,
            activation_tolerance_m: 1e-9,
        },
    }
}
fn audit(report: &PlateImpactNetworkReport) {
    let b = report.balance.as_ref().unwrap();
    close_vec(b.linear_balance_residual, Vec3::ZERO);
    close_vec(b.angular_balance_residual, Vec3::ZERO);
    close(b.work_energy_residual_j, 0.0);
    close(b.energy_residual_j, 0.0);
    close(b.kinetic_energy_change_j, -report.contact.dissipated_energy);
    close(b.relative_impulse_work_j, report.contact.relative_work);
    for (&a, &b) in report
        .relative_normal_speeds_after_m_s
        .iter()
        .zip(&report.contact.relative_normal_speeds_after)
    {
        close(a, b);
    }
}
fn same_bits(a: &SecondOrderState, b: &SecondOrderState) {
    assert_eq!(a, b);
    assert_eq!(a.t.to_bits(), b.t.to_bits());
    for (a, b) in [&a.q, &a.v, &a.a].into_iter().zip([&b.q, &b.v, &b.a]) {
        assert_eq!(
            a.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
            b.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
        );
    }
}

#[test]
fn analytic_three_plate_chain_is_global_permutation_covariant_and_continues_native_owners() {
    let (mesh, model, full) = fixture(false);
    let point = [1.2 + 1.0 / 3.0, -0.7 + 1.0 / 3.0];
    let scales = [1.0, 2.0, 0.7, 1.0];
    let damping = [0.07, 0.17, 0.23, 0.0];
    let velocity = [-2.0, 0.0, 1.0];
    let masses = [3.0, 6.0, 2.1];
    let center = dot(&masses, &velocity) / masses.iter().sum::<f64>();
    // At the centroid the independently integrated consistent P1 mass gives
    // uniform translation, so these are three finite masses, not a modal fit.
    let mobility = [
        1.0 / masses[0] + 1.0 / masses[1],
        -1.0 / masses[1],
        -1.0 / masses[1],
        1.0 / masses[1] + 1.0 / masses[2],
    ];
    for restitution in [0.0, 0.4, 1.0] {
        let expected = velocity.map(|v| (1.0 + restitution) * center - restitution * v);
        let expected_impulse = [
            masses[0] * (expected[0] - velocity[0]),
            masses[2] * (velocity[2] - expected[2]),
        ];
        for permutation in [false, true] {
            let owners = if permutation {
                [2, 0, 3, 1]
            } else {
                [0, 1, 2, 3]
            };
            let order = if permutation { [1, 0] } else { [0, 1] };
            let index = |original| owners.iter().position(|&i| i == original).unwrap();
            let surfaces: Vec<_> = owners
                .iter()
                .map(|&i| surface(&mesh, &model, &full, scales[i], damping[i], 2.0))
                .collect();
            let contacts = order.map(|i| edge(index(i), index(i + 1), point, restitution));
            let prepared = PlateImpactNetwork::prepare(
                &surfaces,
                &contacts,
                pair_controls(),
                &[pair_budget().positive; 4],
                network_budget(),
                &mut || false,
            )
            .unwrap();
            for i in 0..2 {
                for j in 0..2 {
                    close(
                        prepared.contact_mobility()[2 * i + j],
                        mobility[2 * order[i] + order[j]],
                    );
                }
            }
            let mut states: Vec<_> = owners
                .iter()
                .map(|&i| {
                    if i == 3 {
                        // An unreferenced owner is neither clock-synchronized nor
                        // admitted as if it belonged to this event.
                        return SecondOrderState::new(19.0, &[-0.0], &[-0.0], &[7.0]);
                    }
                    let mut s = state(9);
                    for r in [0, 3, 6] {
                        s.v[r] = velocity[i];
                        s.a[r] = -damping[i] * velocity[i];
                    }
                    s
                })
                .collect();
            let original = states.clone();
            let report = prepared.apply(&mut states, &mut || false).unwrap();
            audit(&report);
            for i in 0..2 {
                close(
                    report.contact.normal_impulses[i],
                    expected_impulse[order[i]],
                );
            }
            for i in 0..3 {
                let s = &states[index(i)];
                for r in 0..9 {
                    close(s.v[r], if r % 3 == 0 { expected[i] } else { 0.0 });
                    close(s.a[r], -damping[i] * s.v[r]);
                }
                let before = &original[index(i)];
                assert_eq!(
                    (&s.q, s.t, s.steps, &s.history),
                    (&before.q, before.t, before.steps, &before.history)
                );
                assert!(report.plates[index(i)].is_some());
            }
            same_bits(&states[index(3)], &original[index(3)]);
            assert!(report.plates[index(3)].is_none());
            let loss: f64 = (0..3)
                .map(|i| 0.5 * masses[i] * (velocity[i] * velocity[i] - expected[i] * expected[i]))
                .sum();
            close(report.contact.dissipated_energy, loss);
            if restitution == 0.0 && !permutation {
                let mut sequential = original.clone();
                for i in 0..2 {
                    let pair = PlatePlateImpact::prepare(
                        surface(&mesh, &model, &full, scales[i], damping[i], 2.0),
                        surface(&mesh, &model, &full, scales[i + 1], damping[i + 1], 2.0),
                        &[contacts[i].site],
                        pair_controls(),
                        pair_budget(),
                        &mut || false,
                    )
                    .unwrap();
                    let (left, right) = sequential.split_at_mut(i + 1);
                    pair.apply(&mut left[i], &mut right[0], &mut || false)
                        .unwrap();
                }
                assert!((sequential[0].v[0] - states[0].v[0]).abs() > 0.1);
                assert!(sequential[0].v[0] - sequential[1].v[0] < -0.1);
            }
            for i in 0..3 {
                let s = &mut states[index(i)];
                let before_energy = dot(&s.v, &s.v);
                OperatorGeneralizedAlpha::new(9, 1e-3, 1.0, ImplicitSolveConfig::default())
                    .step(s, &surfaces[index(i)].dynamics, &[0.0; 9])
                    .unwrap();
                assert_eq!((s.steps, s.history.len()), (1, 1));
                assert!(s.t > 0.7 && dot(&s.v, &s.v) < before_energy);
            }
        }
    }
}

#[test]
fn three_nonmatching_plates_use_full_signed_consistent_mass_and_explicit_face_heights() {
    let (mut mesh, model, full) = fixture(false);
    mesh.nodes = vec![(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)];
    let (fine_mesh, fine_model, fine_full) = fine_fixture([0.0, 0.0]);
    let surfaces = [
        surface(&mesh, &model, &full, 1.3, 0.07, 0.04),
        surface(&fine_mesh, &fine_model, &fine_full, 0.8, 0.19, 0.02),
        surface(&mesh, &model, &full, 2.1, 0.31, 0.0),
    ];
    let points = [[0.25, 0.125], [0.75, 0.125], [0.125, 0.75]];
    let contacts = [
        edge(0, 1, points[0], 0.4),
        edge(1, 2, points[1], 0.4),
        edge(0, 2, points[2], 0.4),
    ];
    let prepared = PlateImpactNetwork::prepare(
        &surfaces,
        &contacts,
        pair_controls(),
        &[pair_budget().positive; 3],
        network_budget(),
        &mut || false,
    )
    .unwrap();
    let coarse: Vec<_> = points.map(coarse_row).into_iter().collect();
    let fine = fine_rows();
    let mut rows = vec![
        vec![vec![0.0; 9]; 3],
        vec![vec![0.0; 18]; 3],
        vec![vec![0.0; 9]; 3],
    ];
    for (i, contact) in contacts.iter().enumerate() {
        for (owner, sign) in [
            (contact.positive_plate, 1.0),
            (contact.negative_plate, -1.0),
        ] {
            rows[owner][i] = (if owner == 1 { &fine[i] } else { &coarse[i] })
                .iter()
                .map(|v| sign * v)
                .collect();
        }
    }
    let mut responses = rows.clone();
    for owner in 0..3 {
        let dynamics = surfaces[owner].dynamics;
        let dense: Vec<_> = dynamics
            .model()
            .m
            .to_dense()
            .iter()
            .map(|m| m * dynamics.parameters().mass_scale)
            .collect();
        let factor = lu(&dense, dynamics.model().free).unwrap();
        for row in &mut responses[owner] {
            factor.solve(row);
        }
    }
    let mut matrix = vec![0.0; 9];
    for i in 0..3 {
        for j in 0..3 {
            matrix[3 * i + j] = (0..3).map(|p| dot(&rows[p][i], &responses[p][j])).sum();
            close(matrix[3 * i + j], prepared.contact_mobility()[3 * i + j]);
        }
    }
    let impulse = [0.2; 3];
    let speeds: Vec<_> = (0..3)
        .map(|i| -dot(&matrix[3 * i..3 * i + 3], &impulse) / 1.4)
        .collect();
    assert!(speeds.iter().all(|v| *v < 0.0));
    let mut states: Vec<_> = surfaces
        .iter()
        .enumerate()
        .map(|(p, surface)| {
            let mut s = state(surface.dynamics.model().free);
            for r in 0..s.v.len() {
                s.v[r] = -(0..3).map(|i| impulse[i] * responses[p][i][r]).sum::<f64>() / 1.4;
                s.a[r] = -surface.dynamics.parameters().mass_damping_per_s * s.v[r];
            }
            // All three actual faces coincide, despite different declared reference
            // heights. The two independent original interpolation spaces remain.
            for r in (0..s.q.len()).step_by(3) {
                s.q[r] = -surface.reference_height_m;
            }
            s
        })
        .collect();
    let before = states.clone();
    let report = prepared.apply(&mut states, &mut || false).unwrap();
    audit(&report);
    for i in 0..3 {
        close(report.gaps_m[i], 0.0);
        close(report.relative_normal_speeds_before_m_s[i], speeds[i]);
        close(report.contact.normal_impulses[i], impulse[i]);
        for r in 0..states[i].v.len() {
            close(
                states[i].v[r],
                before[i].v[r] + (0..3).map(|j| impulse[j] * responses[i][j][r]).sum::<f64>(),
            );
            close(
                states[i].a[r],
                -surfaces[i].dynamics.parameters().mass_damping_per_s * states[i].v[r],
            );
        }
    }
}

#[test]
fn support_cross_inertia_redundancy_and_two_plate_specialization_keep_original_physics() {
    let (mesh, model, full) = fixture(false);
    let (supported_mesh, supported_model, original_full) = fixture(true);
    let mut coupled = Coo::new(9, 9);
    for r in 0..9 {
        let (cols, values) = original_full.row(r);
        for (&c, &v) in cols.iter().zip(values) {
            coupled.push(r, c, v);
        }
    }
    coupled.push(1, 3, 0.01);
    coupled.push(3, 1, 0.01);
    let supported_full = coupled.assemble();
    let surfaces = [
        surface(&mesh, &model, &full, 1.0, 0.0, 0.0),
        surface(
            &supported_mesh,
            &supported_model,
            &supported_full,
            1.7,
            0.0,
            0.0,
        ),
        surface(&mesh, &model, &full, 0.9, 0.0, 0.0),
    ];
    let base = [edge(0, 1, [1.2, -0.7], 0.3), edge(1, 2, [2.2, -0.7], 0.3)];
    let mut reference: Option<(Vec<SecondOrderState>, f64)> = None;
    for duplicate in [false, true] {
        let mut contacts = base.to_vec();
        if duplicate {
            contacts.push(base[0]);
        }
        let prepared = PlateImpactNetwork::prepare(
            &surfaces,
            &contacts,
            pair_controls(),
            &[pair_budget().positive; 3],
            network_budget(),
            &mut || false,
        )
        .unwrap();
        let mut states = vec![state(9), state(6), state(9)];
        states[0].v[0] = -1.0;
        states[1].v[0] = -0.7;
        let before = states.clone();
        let report = prepared.apply(&mut states, &mut || false).unwrap();
        audit(&report);
        if duplicate {
            assert!(matches!(
                report.contact.determinacy,
                MultiImpactDeterminacy::UniqueVelocity(_)
            ));
        }
        let middle = report.balance.as_ref().unwrap().plates[1].as_ref().unwrap();
        let mut full_delta = vec![0.0; 9];
        for (r, &free) in supported_model.dof_map.iter().enumerate() {
            if let Some(free) = free {
                full_delta[r] = states[1].v[free] - before[1].v[free];
            }
        }
        let mut actual = vec![0.0; 9];
        supported_full.spmv(&full_delta, &mut actual);
        for &(r, reaction) in &middle.support_impulses {
            let requested: f64 = contacts
                .iter()
                .zip(&report.contact.normal_impulses)
                .map(|(c, &j)| {
                    let sign = if c.positive_plate == 1 { 1.0 } else { -1.0 };
                    sign * j * coarse_row([c.site.point_xy[0] - 1.2, c.site.point_xy[1] + 0.7])[r]
                })
                .sum();
            close(reaction, 1.7 * actual[r] - requested);
            if r == 1 {
                assert!(reaction.abs() > 1e-4);
            }
        }
        if let Some((original, loss)) = &reference {
            for (s, reference) in states.iter().zip(original) {
                for (&a, &b) in s.v.iter().zip(&reference.v) {
                    close(a, b);
                }
            }
            close(report.contact.dissipated_energy, *loss);
        } else {
            reference = Some((states, report.contact.dissipated_energy));
        }
    }
    let contacts = [base[0], edge(0, 1, [2.2, -0.7], 0.3)];
    let network = PlateImpactNetwork::prepare(
        &surfaces[..2],
        &contacts,
        pair_controls(),
        &[pair_budget().positive; 2],
        network_budget(),
        &mut || false,
    )
    .unwrap();
    let pair = PlatePlateImpact::prepare(
        surface(&mesh, &model, &full, 1.0, 0.0, 0.0),
        surface(
            &supported_mesh,
            &supported_model,
            &supported_full,
            1.7,
            0.0,
            0.0,
        ),
        &contacts.map(|c| c.site),
        pair_controls(),
        pair_budget(),
        &mut || false,
    )
    .unwrap();
    for (&a, &b) in network
        .contact_mobility()
        .iter()
        .zip(pair.contact_mobility())
    {
        close(a, b);
    }
    let mut states = [state(9), state(6)];
    states[0].v[0] = -1.0;
    states[0].v[3] = -0.7;
    let mut pair_states = states.clone();
    let report = network.apply(&mut states, &mut || false).unwrap();
    let [positive, negative] = &mut pair_states;
    let pair_report = pair.apply(positive, negative, &mut || false).unwrap();
    for (s, p) in states.iter().zip(&pair_states) {
        for (&a, &b) in s.v.iter().zip(&p.v) {
            close(a, b);
        }
    }
    close(
        report.contact.dissipated_energy,
        pair_report.contact.dissipated_energy,
    );
}

#[test]
fn separating_inactive_and_zero_rows_keep_global_response_and_reject_active_energy_gain() {
    let (mesh, model, full) = fixture(false);
    let surfaces = [
        surface(&mesh, &model, &full, 3.0, 0.0, 0.0),
        surface(&mesh, &model, &full, 2.0, 0.0, 0.0),
        surface(&mesh, &model, &full, 3.0, 0.0, 0.0),
    ];
    let mut contacts = [
        edge(0, 1, [1.2, -0.7], 1.0),
        edge(1, 2, [1.2, -0.7], 1.0),
        edge(0, 2, [1.2, 0.3], 1.0),
    ];
    let make = |contacts: &[PlateNetworkContact]| {
        PlateImpactNetwork::prepare(
            &surfaces,
            contacts,
            pair_controls(),
            &[pair_budget().positive; 3],
            network_budget(),
            &mut || false,
        )
        .unwrap()
    };
    let mut states = [state(9), state(9), state(9)];
    states[0].v[0] = -1.0;
    states[2].v[0] = -0.5;
    states[0].q[6] = 0.01;
    let original = states.clone();
    let active = make(&contacts);
    assert!(matches!(
        active.apply(&mut states, &mut || false),
        Err(PlateRigidImpactError::MultipleContact(
            MultiImpactError::EnergyIncrease { .. }
        ))
    ));
    assert_eq!(states, original);
    for c in &mut contacts {
        c.site.restitution = 0.0;
    }
    let passive = make(&contacts);
    let report = passive.apply(&mut states, &mut || false).unwrap();
    audit(&report);
    // The two active interfaces have W=[[2.5,-1.5],[-1.5,2.5]].
    close(report.contact.normal_impulses[0], 0.4375);
    close(report.contact.normal_impulses[1], 0.0625);
    close(report.contact.normal_impulses[2], 0.0);
    close(report.relative_normal_speeds_after_m_s[0], 0.0);
    close(report.relative_normal_speeds_after_m_s[1], 0.0);
    close(report.relative_normal_speeds_after_m_s[2], -1.0 / 6.0);
    let mut zero = [state(9), state(9), state(9)];
    zero[0].v[4] = -0.0;
    zero[1].a.fill(7.0);
    let before = zero.clone();
    let report = passive.apply(&mut zero, &mut || false).unwrap();
    assert!(report.plates.iter().all(Option::is_none) && report.balance.is_none());
    for (s, b) in zero.iter().zip(&before) {
        same_bits(s, b);
    }
}

#[test]
fn graph_budgets_clocks_cancellation_and_heterogeneous_load_failure_are_atomic() {
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
                return Err("middle original load refused".into());
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
    let held = HeldLoad {
        refuse: Cell::new(false),
    };
    let zero = NoPlateLoad;
    let loads: [&dyn PlateLoad; 3] = [&zero, &held, &zero];
    let surfaces: Vec<_> = loads
        .iter()
        .map(|&load| PlateImpactSurface {
            dynamics: PlateDynamics::new(
                &model,
                PlateDynamicsParameters {
                    stiffness_scale: 1.0,
                    mass_scale: 1.0,
                    mass_damping_per_s: 0.2,
                    stiffness_damping_s: 0.0,
                },
                load,
                PlateDynamicsBudget {
                    max_dofs: 9,
                    max_nonzeros: 100,
                    max_full_dofs: 9,
                    max_load_parameters: 0,
                },
                &mut || false,
            )
            .unwrap(),
            mesh: &mesh,
            full_mass: &full,
            reference_height_m: 0.0,
        })
        .collect();
    let contacts = [edge(0, 1, [1.2, -0.7], 0.4), edge(1, 2, [1.2, -0.7], 0.4)];
    let budgets = [pair_budget().positive; 3];
    let make = |contacts: &[PlateNetworkContact], budgets: &[_], budget| {
        PlateImpactNetwork::prepare(
            &surfaces,
            contacts,
            pair_controls(),
            budgets,
            budget,
            &mut || false,
        )
    };
    for kind in 0..4 {
        let mut budget = network_budget();
        match kind {
            0 => budget.max_plates = 2,
            1 => budget.max_contacts = 1,
            2 => budget.max_total_full_dofs = 26,
            _ => budget.max_workspace_components = 1,
        }
        assert!(matches!(
            make(&contacts, &budgets, budget),
            Err(PlateRigidImpactError::Budget(_))
        ));
    }
    let mut overflow = budgets;
    overflow[1].plate.max_workspace_components = usize::MAX;
    assert!(matches!(
        make(&contacts, &overflow, network_budget()),
        Err(PlateRigidImpactError::Budget(_))
    ));
    let mut limited = budgets;
    limited[1].max_contacts = 1;
    assert!(matches!(
        make(&contacts, &limited, network_budget()),
        Err(PlateRigidImpactError::Budget(_))
    ));
    for malformed in [edge(0, 3, [1.2, -0.7], 0.0), edge(1, 1, [1.2, -0.7], 0.0)] {
        assert!(matches!(
            make(&[malformed], &budgets, network_budget()),
            Err(PlateRigidImpactError::Invalid(_))
        ));
    }
    let prepared = make(&contacts, &budgets, network_budget()).unwrap();
    let mut initial = [state(9), state(9), state(9)];
    initial[0].v[0] = -2.0;
    initial[2].v[0] = 1.0;
    let mut mismatched = initial.clone();
    mismatched[2].t += 0.1;
    let before = mismatched.clone();
    assert!(matches!(
        prepared.apply(&mut mismatched, &mut || false),
        Err(PlateRigidImpactError::Invalid(_))
    ));
    assert_eq!(mismatched, before);
    held.refuse.set(true);
    let mut refused = initial.clone();
    assert!(matches!(
        prepared.apply(&mut refused, &mut || false),
        Err(PlateRigidImpactError::Load(_))
    ));
    assert_eq!(refused, initial);
    held.refuse.set(false);
    let mut accepted = initial.clone();
    let mut polls = 0;
    prepared
        .apply(&mut accepted, &mut || {
            polls += 1;
            false
        })
        .unwrap();
    for limit in [0, polls / 2, polls - 1] {
        let mut candidate = initial.clone();
        let mut count = 0;
        assert!(
            prepared
                .apply(&mut candidate, &mut || {
                    let stop = count == limit;
                    count += 1;
                    stop
                })
                .is_err()
        );
        assert_eq!(candidate, initial);
    }
    let mut retry = initial.clone();
    prepared.apply(&mut retry, &mut || false).unwrap();
    assert_eq!(retry, accepted);
    // Restart uses the original heterogeneous load and damping. Independently
    // reconstruct M*a + C*v = f on the middle plate after publication.
    let balance: Vec<_> = accepted[1]
        .a
        .iter()
        .zip(&accepted[1].v)
        .map(|(a, v)| a + 0.2 * v)
        .collect();
    let mut force = vec![0.0; 9];
    full.spmv(&balance, &mut force);
    for (i, &f) in force.iter().enumerate() {
        close(f, if i == 0 { 1.5 } else { 0.0 });
    }
}
