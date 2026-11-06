//! Simultaneous original junctions on actual native modal and memory states.
use super::*;
use core::f64::consts::PI;
use fs_contact::interface_binding::adhesive_shear::SourceAdhesiveShear;
use fs_contact::normal_patch::dmt::SourceDmtSphereContact;
use fs_couple::render::plate::impact::radiation::{Model, Pole};
use fs_couple::render::schedule::force::coupled::contact::adhesive_dmt::multiple::{
    AdhesiveDmtModalContact, MultipleAdhesiveDmtModalConfig, MultipleAdhesiveDmtModalError,
    MultipleAdhesiveDmtModalSystem,
};
use fs_phs::RelaxationBranch;

mod temperature;

fn gate() -> CancelGate {
    CancelGate::new_clock_free()
}
fn controls() -> MultipleAdhesiveDmtModalConfig {
    MultipleAdhesiveDmtModalConfig {
        maximum_contacts: 4,
        maximum_sweeps: 96,
        maximum_setup_terms: 262_144,
        maximum_source_evaluations: 24_000,
        contact_work_absolute_tolerance_j: 1e-26,
        contact_work_relative_tolerance: 2e-8,
    }
}
fn contact<'a>(
    n: &'a SourceDmtSphereContact,
    s: &'a SourceAdhesiveShear,
    initial: [f64; 3],
    columns: [Vec<f64>; 3],
) -> AdhesiveDmtModalContact<'a> {
    let m = n.normal_state().elastic().interface().surface_materials();
    let mut limits = config();
    limits.minimum_indentation_m = -100.0 * Z;
    AdhesiveDmtModalContact {
        normal: n,
        shear: s,
        surface_materials: [&m[0], &m[1]],
        initial_position_m: initial,
        columns,
        config: limits,
    }
}
fn second_columns() -> [Vec<f64>; 3] {
    // A second physical normal/tangent frame sharing the SAME original bodies.
    let rotation = [[0.8, 0.6, 0.0], [-0.36, 0.48, 0.8], [0.48, -0.64, 0.6]];
    let original = columns();
    core::array::from_fn(|i| {
        (0..6)
            .map(|j| (0..3).map(|k| rotation[i][k] * original[k][j]).sum())
            .collect()
    })
}
fn native_states(n: &CoupledModalSystem) -> Vec<ModalAcousticState> {
    n.components()
        .iter()
        .flat_map(|m| m.states().iter().copied())
        .collect()
}
fn bits(n: &CoupledModalSystem) -> Vec<u64> {
    native_states(n)
        .into_iter()
        .flat_map(|s| {
            [
                s.displacement_m_sqrt_kg.to_bits(),
                s.velocity_m_sqrt_kg_per_s.to_bits(),
            ]
        })
        .chain(n.material_memory().iter().map(|x| x.to_bits()))
        .chain(n.radiation_state().iter().map(|x| x.to_bits()))
        .chain([n.samples_rendered(), n.total_energy_j().unwrap().to_bits()])
        .collect()
}
fn rich_network() -> CoupledModalSystem {
    network([0.0; 3], true, RATE)
        .with_material_relaxation(
            vec![RelaxationBranch {
                projection: columns()[0].iter().flat_map(|&b| [b, 0.0]).collect(),
                stiffness: 1000.0,
                relaxation_time_s: 3e-6,
            }],
            Some(vec![1e-12]),
            &gate(),
        )
        .unwrap()
        .with_radiation_load(
            &Model {
                ports: 2,
                poles: vec![Pole {
                    omega: 1e6,
                    zeta: 0.1,
                    coupling: vec![5000.0, -3000.0],
                }],
            },
            &[0, 3],
            1,
            &gate(),
        )
        .unwrap()
}
fn same_native(a: &CoupledModalSystem, b: &CoupledModalSystem) {
    for (x, y) in native_states(a).iter().zip(native_states(b)) {
        close(
            x.displacement_m_sqrt_kg,
            y.displacement_m_sqrt_kg,
            1e-26,
            3e-8,
        );
        close(
            x.velocity_m_sqrt_kg_per_s,
            y.velocity_m_sqrt_kg_per_s,
            1e-23,
            3e-8,
        );
    }
    for (x, y) in a.material_memory().iter().zip(b.material_memory()) {
        close(*x, *y, 1e-24, 3e-8);
    }
    for (x, y) in a.radiation_state().iter().zip(b.radiation_state()) {
        close(*x, *y, 1e-27, 3e-8);
    }
}

#[test]
fn one_original_junction_matches_the_single_native_consumer_on_all_dmt_phases() {
    let (n, s) = sources(0.01, 1e7, 6e6);
    for initial in [
        [0.2 * Z, 0.0, 0.0],
        [-0.6 * Z, 0.0, 0.0],
        [-3.0 * Z, 0.0, 0.0],
    ] {
        let velocity = [0.0, 1e-3, -0.5e-3];
        let mut single = make(&n, &s, velocity, initial, true, RATE);
        let mut several = MultipleAdhesiveDmtModalSystem::new(
            network(velocity, true, RATE),
            vec![contact(&n, &s, initial, columns())],
            controls(),
            &gate(),
        )
        .unwrap();
        assert!(std::ptr::eq(several.normal(0).unwrap(), &n));
        assert!(std::ptr::eq(several.shear(0).unwrap(), &s));
        assert_eq!(several.state(0).unwrap().regime(), single.state().regime());
        for _ in 0..5 {
            let a = *single.step(&[0.0; 6], [0.0; 3]).unwrap();
            let b = several.step(&[0.0; 6], &[[0.0; 3]]).unwrap().clone();
            let point = &b.contacts[0];
            for axis in 0..3 {
                close(
                    a.position_m[axis],
                    point.position_after_m[axis],
                    2e-21,
                    2e-8,
                );
                close(a.mean_force_n[axis], point.mean_force_n[axis], 2e-17, 2e-8);
            }
            assert_eq!(a.regime, point.regime);
            close(a.surface_energy_j, b.surface_energy_j, 1e-26, 2e-8);
            close(
                a.friction_dissipation_j,
                b.friction_dissipation_j,
                1e-26,
                2e-8,
            );
            close(
                single.total_energy_j().unwrap(),
                several.total_energy_j().unwrap(),
                1e-25,
                2e-8,
            );
            assert_eq!(b.sample, single.samples_rendered());
            assert!(b.contact_work_residual_j.abs() <= b.contact_work_tolerance_j);
            assert!(b.energy_residual_j.abs() <= b.energy_tolerance_j);
            assert!(b.evaluations <= controls().maximum_source_evaluations);
            for (x, y) in states(&single).iter().zip(native_states(several.network())) {
                close(
                    x.displacement_m_sqrt_kg,
                    y.displacement_m_sqrt_kg,
                    1e-26,
                    2e-8,
                );
                close(
                    x.velocity_m_sqrt_kg_per_s,
                    y.velocity_m_sqrt_kg_per_s,
                    1e-23,
                    2e-8,
                );
            }
        }
    }
    // A fresh solid junction can support a sub-static prescribed displacement
    // with zero relative slip. Native recoil cancels the support motion; the
    // remaining original coordinate-roundoff work stays inside the independent
    // authored contact-work tolerance and never becomes frictional heat.
    let initial = [0.2 * Z, 0.0, 0.0];
    let mut single = make(&n, &s, [0.0; 3], initial, false, RATE);
    let mut several = MultipleAdhesiveDmtModalSystem::new(
        network([0.0; 3], false, RATE),
        vec![contact(&n, &s, initial, columns())],
        controls(),
        &gate(),
    )
    .unwrap();
    let response = n.evaluate(initial[0]).unwrap();
    let reaction = 0.25 * s.law().static_shear_strength_pa() * response.contact_area_m2;
    let rate = [
        0.0,
        reaction * single.compliance_m_per_n()[1][1] / single.sample_period_s(),
        0.0,
    ];
    let external: Vec<_> = columns()[0]
        .iter()
        .map(|b| b * response.normal_force_n)
        .collect();
    let a = *single.step(&external, rate).unwrap();
    let b = several.step(&external, &[rate]).unwrap();
    assert_eq!(a.regime, DmtCoupledRegime::Sticking);
    assert_eq!(b.contacts[0].regime, DmtCoupledRegime::Sticking);
    assert_eq!(b.friction_dissipation_j, 0.0);
    close(b.contacts[0].mean_force_n[1], reaction, 1e-18, 1e-8);
    close(
        b.contacts[0].position_after_m[1],
        0.0,
        config().coordinate_absolute_tolerance_m,
        0.0,
    );
    assert!(b.contacts[0].work_residual_j.abs() <= b.contacts[0].work_tolerance_j);
    assert!(b.contact_work_residual_j.abs() <= b.contact_work_tolerance_j);
}

// Independent original potential work and analytic real-area integral. The
// manufactured paths below stay on one side of touch, so no numerical area
// approximation or copy of the coupled solver is needed.
fn path_force(
    n: &SourceDmtSphereContact,
    s: &SourceAdhesiveShear,
    before: [f64; 3],
    after: [f64; 3],
) -> [f64; 3] {
    let normal = n.reversible_work_between(before[0], after[0]).unwrap() / (after[0] - before[0]);
    let area = 0.5 * PI * n.effective_radius_m() * (before[0].max(0.0) + after[0].max(0.0));
    let slip = [after[1] - before[1], after[2] - before[2]];
    let length = slip[0].hypot(slip[1]);
    let shear = s.law().kinetic_shear_strength_pa() * area;
    [normal, shear * slip[0] / length, shear * slip[1] / length]
}

#[test]
fn shared_cross_contact_response_matches_manufactured_paths_and_all_native_memories() {
    let (n0, s0) = sources(0.01, 1e7, 6e6);
    let (n1, s1) = sources(0.018, 8e6, 4e6);
    for gap in [false, true] {
        let before = [
            [0.20 * Z, 0.0, 0.0],
            [if gap { -0.6 * Z } else { 0.14 * Z }, 0.0, 0.0],
        ];
        let target = [
            [0.24 * Z, 0.20 * Z, -0.13 * Z],
            [if gap { -0.7 * Z } else { 0.18 * Z }, -0.18 * Z, 0.09 * Z],
        ];
        let maps = [columns(), second_columns()];
        let make_owner = |reverse| {
            let mut inputs = vec![
                contact(&n0, &s0, before[0], maps[0].clone()),
                contact(&n1, &s1, before[1], maps[1].clone()),
            ];
            if reverse {
                inputs.reverse();
            }
            MultipleAdhesiveDmtModalSystem::new(rich_network(), inputs, controls(), &gate())
                .unwrap()
        };
        let mut owner = make_owner(false);
        let mut permuted = make_owner(true);
        let c = owner.compliance_m_per_n();
        assert_eq!(c.len(), 36);
        assert!(c[4] != 0.0 && c[9] != 0.0 && c[17] != 0.0);
        let force = [
            path_force(&n0, &s0, before[0], target[0]),
            path_force(&n1, &s1, before[1], target[1]),
        ];
        let external = [0.0; 6];
        let mut free_native = rich_network();
        let initial_native = native_states(&free_native);
        free_native.step(&external).unwrap();
        let delta: Vec<_> = native_states(&free_native)
            .iter()
            .zip(&initial_native)
            .map(|(a, b)| a.displacement_m_sqrt_kg - b.displacement_m_sqrt_kg)
            .collect();
        let mut rates = [[0.0; 3]; 2];
        let mut omitted_cross = 0.0_f64;
        for i in 0..2 {
            for axis in 0..3 {
                let response: f64 = (0..2)
                    .flat_map(|j| (0..3).map(move |a| (j, a)))
                    .map(|(j, a)| c[(3 * i + axis) * 6 + 3 * j + a] * force[j][a])
                    .sum();
                let cross: f64 = (0..3)
                    .map(|a| c[(3 * i + axis) * 6 + 3 * (1 - i) + a] * force[1 - i][a])
                    .sum();
                omitted_cross = omitted_cross.max(cross.abs());
                let free: f64 = maps[i][axis].iter().zip(&delta).map(|(b, q)| b * q).sum();
                rates[i][axis] =
                    (target[i][axis] - before[i][axis] + response - free) / owner.sample_period_s();
            }
        }
        assert!(omitted_cross > 1e-4 * Z);
        let mut reference = rich_network();
        let mut held = external;
        for mode in 0..6 {
            for i in 0..2 {
                for axis in 0..3 {
                    held[mode] -= maps[i][axis][mode] * force[i][axis];
                }
            }
        }
        let expected_native = reference.step(&held).unwrap().clone();
        let frame = owner.step(&external, &rates).unwrap().clone();
        let reversed = permuted
            .step(&external, &[rates[1], rates[0]])
            .unwrap()
            .clone();
        assert_eq!(owner.samples_rendered(), 1);
        assert!(frame.sweeps > 1);
        same_native(owner.network(), &reference);
        same_native(owner.network(), permuted.network());
        let mut heat = 0.0;
        let mut normal_work = 0.0;
        let mut prescribed_work = 0.0;
        for i in 0..2 {
            let point = &frame.contacts[i];
            let other = &reversed.contacts[1 - i];
            for axis in 0..3 {
                close(point.position_after_m[axis], target[i][axis], 2e-21, 2e-8);
                close(point.mean_force_n[axis], force[i][axis], 2e-17, 2e-8);
                close(
                    other.position_after_m[axis],
                    point.position_after_m[axis],
                    2e-21,
                    2e-8,
                );
                close(
                    other.mean_force_n[axis],
                    point.mean_force_n[axis],
                    2e-17,
                    2e-8,
                );
                prescribed_work += force[i][axis] * rates[i][axis] * owner.sample_period_s();
            }
            let displacement = [target[i][1] - before[i][1], target[i][2] - before[i][2]];
            let expected_heat = force[i][1] * displacement[0] + force[i][2] * displacement[1];
            close(point.friction_dissipation_j, expected_heat, 1e-27, 2e-8);
            assert!(point.work_residual_j.abs() <= point.work_tolerance_j);
            heat += expected_heat;
            normal_work += force[i][0] * (target[i][0] - before[i][0]);
        }
        assert!(frame.contacts[0].friction_dissipation_j > 0.0);
        if gap {
            assert!(frame.contacts[1].mean_force_n[0] < 0.0);
            assert_eq!(frame.contacts[1].path.mean_real_area_m2, 0.0);
            assert_eq!(frame.contacts[1].friction_dissipation_j, 0.0);
            assert_eq!(frame.contacts[1].regime, DmtCoupledRegime::Released);
        } else {
            assert!(frame.contacts[1].friction_dissipation_j > 0.0);
        }
        close(frame.friction_dissipation_j, heat, 1e-26, 2e-8);
        close(frame.prescribed_work_j, prescribed_work, 1e-26, 2e-8);
        close(
            frame.contact_work_on_network_j + normal_work + heat,
            prescribed_work,
            1e-26,
            2e-8,
        );
        close(
            frame.network_dissipation_j,
            expected_native.total_dissipation_j(),
            1e-28,
            3e-8,
        );
        let native = owner.network().last_frame().unwrap();
        assert!(native.material_dissipation_j > 0.0 && native.radiation_dissipation_j > 0.0);
        close(
            native.material_dissipation_j,
            expected_native.material_dissipation_j,
            1e-28,
            3e-8,
        );
        close(
            native.radiation_dissipation_j,
            expected_native.radiation_dissipation_j,
            1e-31,
            3e-8,
        );
        assert!(frame.contact_work_residual_j.abs() <= frame.contact_work_tolerance_j);
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
    }
}

#[test]
fn receiver_cancellation_and_joint_budgets_never_publish_partial_collisions() {
    let (n0, s0) = sources(0.01, 1e7, 6e6);
    let (n1, s1) = sources(0.018, 8e6, 4e6);
    let make_owner = |c| {
        MultipleAdhesiveDmtModalSystem::new(
            network([0.0, 1e-3, -0.5e-3], true, RATE),
            vec![
                contact(&n0, &s0, [0.2 * Z, 0.0, 0.0], columns()),
                contact(&n1, &s1, [-0.7 * Z, 0.0, 0.0], second_columns()),
            ],
            c,
            &gate(),
        )
        .unwrap()
    };
    let mut owner = make_owner(controls());
    let mut replay = make_owner(controls());
    owner.step(&[0.0; 6], &[[0.0; 3]; 2]).unwrap();
    replay.step(&[0.0; 6], &[[0.0; 3]; 2]).unwrap();
    let native = bits(owner.network());
    let frame = owner.last_frame().unwrap().clone();
    let original_native = owner.network().last_frame().unwrap().clone();
    let original_contacts = [owner.state(0), owner.state(1)];
    let refused: Result<(), MultipleAdhesiveDmtModalError> =
        owner.step_with_admission(&[0.0; 6], &[[0.0; 3]; 2], &gate(), |_, _| {
            Err(MultipleAdhesiveDmtModalError::EvaluationBudget)
        });
    assert!(refused.is_err());
    let cancelled = gate();
    let refused: Result<(), MultipleAdhesiveDmtModalError> = owner.step_with_admission(
        &[0.0; 6],
        &[[0.0; 3]; 2],
        &cancelled,
        |candidate, native| {
            assert_eq!(candidate.sample, 2);
            assert_eq!(candidate.contacts.len(), 2);
            assert_eq!(native.len(), 2);
            cancelled.request();
            Ok(())
        },
    );
    assert!(matches!(
        refused,
        Err(MultipleAdhesiveDmtModalError::Native(
            ModalCouplingError::Cancelled
        ))
    ));
    assert!(owner.step(&[0.0; 6], &[[0.0; 3]]).is_err());
    assert!(owner.step(&[0.0; 6], &[[1e9, 0.0, 0.0], [0.0; 3]]).is_err());
    assert_eq!(bits(owner.network()), native);
    assert_eq!(owner.last_frame().unwrap(), &frame);
    assert_eq!(owner.network().last_frame().unwrap(), &original_native);
    assert_eq!([owner.state(0), owner.state(1)], original_contacts);
    assert_eq!(
        owner.step(&[0.0; 6], &[[0.0; 3]; 2]).unwrap(),
        replay.step(&[0.0; 6], &[[0.0; 3]; 2]).unwrap()
    );
    assert_eq!(bits(owner.network()), bits(replay.network()));
    let mut budget = controls();
    budget.maximum_source_evaluations = 1;
    let mut bounded = make_owner(budget);
    let before = bits(bounded.network());
    assert!(matches!(
        bounded.step(&[0.0; 6], &[[0.0; 3]; 2]),
        Err(MultipleAdhesiveDmtModalError::EvaluationBudget)
    ));
    assert_eq!(bits(bounded.network()), before);
    assert!(bounded.last_frame().is_none());
    let mut budget = controls();
    budget.maximum_sweeps = 1;
    let mut bounded = make_owner(budget);
    let before = bits(bounded.network());
    assert!(matches!(
        bounded.step(&[0.0; 6], &[[0.0; 3]; 2]),
        Err(MultipleAdhesiveDmtModalError::NotConverged { sweeps: 1 })
    ));
    assert_eq!(bits(bounded.network()), before);
}

#[test]
fn shared_negative_curvature_is_screened_even_when_every_local_source_is_admissible() {
    let (n, s) = sources(0.01, 1e7, 6e6);
    let base = make(&n, &s, [0.0; 3], [0.2 * Z, 0.0, 0.0], false, RATE);
    let tail = -n.evaluate(0.0).unwrap().one_sided_tangents_n_per_m[0];
    let a = 0.5 * tail;
    let c = 0.25 * PI * n.effective_radius_m() * s.law().kinetic_shear_strength_pa();
    let eta = 0.5 * (a + a.hypot(2.0 * c));
    let scale = (0.55 / (eta * base.compliance_m_per_n()[0][0])).sqrt();
    let maps = || columns().map(|v| v.into_iter().map(|x| scale * x).collect());
    let local = base
        .compliance_m_per_n()
        .map(|v| v.map(|x| scale * scale * x));
    DmtCoupledLoadLine::new(&n, &s, local, config()).unwrap();
    let one = MultipleAdhesiveDmtModalSystem::new(
        network([0.0; 3], false, RATE),
        vec![contact(&n, &s, [0.2 * Z, 0.0, 0.0], maps())],
        controls(),
        &gate(),
    );
    assert!(one.is_ok());
    let two = MultipleAdhesiveDmtModalSystem::new(
        network([0.0; 3], false, RATE),
        vec![
            contact(&n, &s, [0.2 * Z, 0.0, 0.0], maps()),
            contact(&n, &s, [0.2 * Z, 0.0, 0.0], maps()),
        ],
        controls(),
        &gate(),
    );
    assert!(matches!(
        two,
        Err(MultipleAdhesiveDmtModalError::CoupledStability)
    ));
    let mut wrong = contact(&n, &s, [0.2 * Z, 0.0, 0.0], columns());
    wrong.surface_materials.swap(0, 1);
    assert!(
        MultipleAdhesiveDmtModalSystem::new(
            network([0.0; 3], false, RATE),
            vec![wrong],
            controls(),
            &gate()
        )
        .is_err()
    );
    let mut budget = controls();
    budget.maximum_setup_terms = 0;
    assert!(
        MultipleAdhesiveDmtModalSystem::new(
            network([0.0; 3], false, RATE),
            vec![contact(&n, &s, [0.2 * Z, 0.0, 0.0], columns())],
            budget,
            &gate()
        )
        .is_err()
    );
}
