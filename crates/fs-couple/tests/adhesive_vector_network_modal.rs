//! G1/G3/G4: sourced Coulomb disks act on all original mechanical coordinates.
#[path = "../../fs-contact/tests/adhesive_dmt_coupled/common.rs"]
mod common;
use common::{Z, close, sources};
use fs_contact::interface_binding::adhesive_shear::SourceAdhesiveShear;
use fs_contact::normal_patch::dmt::SourceDmtSphereContact;
use fs_couple::modal_acoustic_time::{ModalAcousticMode, ModalAcousticState,
    ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_couple::render::schedule::force::coupled::contact::adhesive_events::linear_network::vector::{
    VectorFrictionRegime, VectorLinearNetworkConfig};
use fs_couple::render::schedule::force::coupled::contact::adhesive_network_events::vector::*;
use fs_couple::render::schedule::force::coupled::{CoupledModalSystem, ModalAttachment,
    ModalConnection, ModalCouplingConfig, ModalCouplingError, PreparedModalBudget};
use fs_exec::CancelGate;
use fs_math::{c64::C64, det};
use fs_ports::{ModuleWindowAudit, PortedModule};

const MASS: [f64; 3] = [2e-9, 3e-9, 4e-9];
const DT: f64 = 1e-4;
fn controls(cells: usize) -> VectorAdhesiveNetworkModalConfig {
    VectorAdhesiveNetworkModalConfig {
        flow: VectorLinearNetworkConfig {
            maximum_modes: 9,
            maximum_connections: 2,
            maximum_contacts: 2,
            maximum_workspace_bytes: 1024 * 1024,
            maximum_cells: 4096,
            maximum_nonlinear_iterations: 100_000,
            maximum_backtracks: 30,
            maximum_step_s: DT,
            maximum_cell_s: DT / cells as f64,
            minimum_column_pivot: 1e-10,
            force_absolute_tolerance_n: 1e-16,
            force_relative_tolerance: 1e-10,
            velocity_absolute_tolerance_m_s: 1e-12,
            maximum_abs_modal_coordinate: 1.0,
            maximum_abs_modal_velocity: 1.0,
            energy_absolute_tolerance_j: 1e-23,
            energy_relative_tolerance: 1e-8,
            maximum_numerical_work_uncertainty_j: 1e-20,
            maximum_numerical_dissipation_j: 1e-10,
        },
        maximum_contact_force_n: 1.0,
        maximum_tangent_speed_m_s: 1000.0,
        maximum_accumulated_slip_m: 0.1,
    }
}
fn columns(elastic: bool) -> Vec<[Vec<f64>; 2]> {
    let width = if elastic { 3 } else { 2 };
    (0..2)
        .map(|site| {
            let mut b = [vec![0.0; 3 * width], vec![0.0; 3 * width]];
            for axis in 0..2 {
                b[axis][site * width + axis] = 1.0 / det::sqrt(MASS[site]);
                b[axis][(site + 1) * width + axis] = -1.0 / det::sqrt(MASS[site + 1]);
                if elastic {
                    b[axis][site * width + 2] = [0.2, -0.3][axis] / det::sqrt(MASS[site]);
                    b[axis][(site + 1) * width + 2] =
                        [-0.1, 0.25][axis] / det::sqrt(MASS[site + 1]);
                }
            }
            b
        })
        .collect()
}
fn network(
    elastic: bool,
    moving: bool,
    displacement_limit: f64,
    absolute: f64,
    relative: f64,
) -> CoupledModalSystem {
    let models = MASS
        .into_iter()
        .enumerate()
        .map(|(i, mass)| {
            let mut budget = ModalAcousticTimeBudget::audible_reference();
            budget.maximum_abs_displacement_m_sqrt_kg = displacement_limit;
            let mut modes = vec![
                ModalAcousticMode {
                    angular_frequency_rad_s: 0.0,
                    damping_ratio: 0.0,
                    pressure_per_modal_velocity: C64::ZERO
                };
                2
            ];
            if elastic {
                modes.push(ModalAcousticMode {
                    angular_frequency_rad_s: [1800.0, 2300.0, 2700.0][i],
                    damping_ratio: 0.03,
                    pressure_per_modal_velocity: C64 { re: 0.5, im: 0.1 },
                });
            }
            let mut model =
                ModalAcousticTimeModel::try_new_with_free_coordinates(10_000, modes, budget)
                    .unwrap();
            let mut states = (0..2)
                .map(|axis| ModalAcousticState {
                    displacement_m_sqrt_kg: if elastic {
                        det::sqrt(mass) * [[2e-8, -1e-8], [-1e-8, 1e-8], [1e-8, 2e-8]][i][axis]
                    } else {
                        0.0
                    },
                    velocity_m_sqrt_kg_per_s: if moving {
                        det::sqrt(mass)
                            * [[0.002, -0.001], [0.0, 0.001], [-0.001, -0.0005]][i][axis]
                    } else {
                        0.0
                    },
                })
                .collect::<Vec<_>>();
            if elastic {
                states.push(ModalAcousticState {
                    displacement_m_sqrt_kg: det::sqrt(mass) * [3e-8, -2e-8, 1e-8][i],
                    velocity_m_sqrt_kg_per_s: if moving {
                        det::sqrt(mass) * [0.0003, -0.0005, 0.0002][i]
                    } else {
                        0.0
                    },
                });
            }
            model.restore_states(&states).unwrap();
            model
        })
        .collect();
    let links = if elastic {
        (0..2)
            .map(|i| ModalConnection {
                left: ModalAttachment {
                    component: i,
                    shapes: [1.0, 0.3, 0.2].map(|x| x / det::sqrt(MASS[i])).to_vec(),
                },
                right: ModalAttachment {
                    component: i + 1,
                    shapes: [1.0, 0.3, -0.1]
                        .map(|x| x / det::sqrt(MASS[i + 1]))
                        .to_vec(),
                },
                stiffness_n_m: [0.04, 0.025][i],
                damping_n_s_m: [2e-7, 1e-7][i],
                rest_extension_m: [1e-8, -0.5e-8][i],
            })
            .collect()
    } else {
        vec![]
    };
    CoupledModalSystem::new_prepared(
        models,
        links,
        ModalCouplingConfig {
            max_modes: 9,
            max_connections: 2,
            max_setup_terms: 4096,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: 1e-5,
            maximum_abs_pressure_pa: 1.0,
            maximum_abs_connection_force_n: 1.0,
            solve_relative_tolerance: 1e-11,
            energy_absolute_tolerance_j: absolute,
            energy_relative_tolerance: relative,
        },
        PreparedModalBudget {
            max_modes: 9,
            max_dense_bytes: 1024 * 1024,
            max_multiply_terms: 100_000_000,
        },
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn junctions<'a>(
    normal: [&'a SourceDmtSphereContact; 2],
    shear: [&'a SourceAdhesiveShear; 2],
    elastic: bool,
) -> Vec<VectorAdhesiveNetworkJunction<'a>> {
    columns(elastic)
        .into_iter()
        .enumerate()
        .map(|(i, tangent_columns)| {
            let materials = normal[i]
                .normal_state()
                .elastic()
                .interface()
                .surface_materials();
            VectorAdhesiveNetworkJunction {
                normal: normal[i],
                shear: shear[i],
                surface_materials: [&materials[0], &materials[1]],
                indentation_m: 0.5 * Z,
                tangent_columns,
            }
        })
        .collect()
}
fn owner<'a>(
    normal: [&'a SourceDmtSphereContact; 2],
    shear: [&'a SourceAdhesiveShear; 2],
    elastic: bool,
    moving: bool,
    cells: usize,
) -> VectorAdhesiveNetworkModalSystem<'a> {
    VectorAdhesiveNetworkModalSystem::new(
        network(elastic, moving, 1.0, 1e-23, 2e-8),
        junctions(normal, shear, elastic),
        controls(cells),
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn states(system: &VectorAdhesiveNetworkModalSystem<'_>) -> Vec<ModalAcousticState> {
    system
        .components()
        .iter()
        .flat_map(|m| m.states().iter().copied())
        .collect()
}
fn momentum(system: &VectorAdhesiveNetworkModalSystem<'_>) -> [f64; 2] {
    core::array::from_fn(|a| {
        MASS.iter()
            .zip(system.components())
            .map(|(m, b)| det::sqrt(*m) * b.states()[a].velocity_m_sqrt_kg_per_s)
            .sum()
    })
}
fn same(a: f64, b: f64) {
    close(a, b, 2e-22, 5e-8);
}

#[test]
fn diagonal_load_outside_disk_slides_even_when_both_components_are_below_capacity() {
    let (n0, s0) = sources(1e-5, 1e5, 1e5);
    let (n1, s1) = sources(2e-5, 0.0, 0.0);
    let mut system = owner([&n0, &n1], [&s0, &s1], false, false, 8);
    let capacity = n0.evaluate(0.5 * Z).unwrap().contact_area_m2 * 1e5;
    let b = columns(false);
    let external: Vec<_> = (0..6)
        .map(|k| 0.8 * capacity * (b[0][0][k] + b[0][1][k]))
        .collect();
    let before = system.total_energy_j().unwrap();
    let frame = *system.step(&external).unwrap();
    assert_eq!(frame.path.regimes[0], VectorFrictionRegime::Sliding);
    let expected = (1.0 / MASS[0] + 1.0 / MASS[1]) * DT * capacity * (0.8 - 1.0 / det::sqrt(2.0));
    for a in 0..2 {
        same(frame.tangent_velocities_m_s[0][a], expected);
        same(
            frame.path.contact_impulse_n_s[0][a],
            capacity * DT / det::sqrt(2.0),
        );
        same(momentum(&system)[a], 0.0);
    }
    assert!(frame.friction_dissipation_j > 0.0 && frame.numerical_dissipation_j > 0.0);
    same(
        system.total_energy_j().unwrap()
            + frame.friction_dissipation_j
            + frame.numerical_dissipation_j
            - frame.external_work_j,
        before,
    );
    let native = system.network().last_frame().unwrap();
    assert_eq!(
        native.numerical_dissipation_j,
        frame.numerical_dissipation_j
    );
    assert_eq!(native.total_dissipation_j(), 0.0);
    let ModuleWindowAudit::Available(audit) = system.network().audit_last_window().unwrap() else {
        panic!("missing native audit")
    };
    assert!(audit.violations().is_empty());
    assert!(format!("{audit:?}").contains("numerical-dissipation"));
}

#[test]
fn joint_sticking_keeps_internal_elastic_motion_and_all_original_loss_channels() {
    let (n0, s0) = sources(1e-5, 2e8, 2e8);
    let (n1, s1) = sources(2e-5, 3e8, 3e8);
    let mut system = owner([&n0, &n1], [&s0, &s1], true, false, 8);
    let old = states(&system);
    let b = columns(true);
    let old_y: Vec<[f64; 2]> = b
        .iter()
        .map(|cols| {
            core::array::from_fn(|a| {
                cols[a]
                    .iter()
                    .zip(&old)
                    .map(|(b, s)| b * s.displacement_m_sqrt_kg)
                    .sum()
            })
        })
        .collect();
    for _ in 0..8 {
        let before = system.total_energy_j().unwrap();
        let frame = *system.step(&[0.0; 9]).unwrap();
        assert_eq!(
            &frame.path.regimes[..2],
            &[VectorFrictionRegime::Sticking; 2]
        );
        assert!(frame.modal_damping_dissipation_j > 0.0);
        assert!(frame.connection_damping_dissipation_j > 0.0);
        assert!(frame.numerical_dissipation_j > 0.0);
        assert!(frame.friction_dissipation_j < frame.energy_tolerance_j);
        for (i, y) in old_y.iter().enumerate() {
            for a in 0..2 {
                close(frame.tangent_positions_m[i][a], y[a], 1e-17, 1e-9);
            }
        }
        same(
            system.total_energy_j().unwrap()
                + frame.friction_dissipation_j
                + frame.modal_damping_dissipation_j
                + frame.connection_damping_dissipation_j
                + frame.numerical_dissipation_j,
            before,
        );
        let ModuleWindowAudit::Available(audit) = system.network().audit_last_window().unwrap()
        else {
            panic!("missing native audit")
        };
        assert!(audit.violations().is_empty());
    }
    assert_ne!(states(&system), old);
    for p in momentum(&system) {
        same(p, 0.0);
    }
}

#[test]
fn tangent_chart_rotation_and_contact_permutation_preserve_original_physics() {
    let (n0, s0) = sources(1e-5, 1e5, 1e5);
    let (n1, s1) = sources(2e-5, 2e5, 2e5);
    let mut a = owner([&n0, &n1], [&s0, &s1], true, true, 8);
    let mut sites = junctions([&n0, &n1], [&s0, &s1], true);
    let (c, s) = (0.6, 0.8);
    for site in &mut sites {
        let original = site.tangent_columns.clone();
        for k in 0..9 {
            site.tangent_columns[0][k] = c * original[0][k] + s * original[1][k];
            site.tangent_columns[1][k] = -s * original[0][k] + c * original[1][k];
        }
    }
    sites.swap(0, 1);
    let mut b = VectorAdhesiveNetworkModalSystem::new(
        network(true, true, 1.0, 1e-23, 2e-8),
        sites,
        controls(8),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    for _ in 0..8 {
        let fa = *a.step(&[0.0; 9]).unwrap();
        let fb = *b.step(&[0.0; 9]).unwrap();
        for (a, b) in states(&a).iter().zip(states(&b)) {
            same(a.displacement_m_sqrt_kg, b.displacement_m_sqrt_kg);
            same(a.velocity_m_sqrt_kg_per_s, b.velocity_m_sqrt_kg_per_s);
        }
        same(fa.friction_dissipation_j, fb.friction_dissipation_j);
        same(fa.numerical_dissipation_j, fb.numerical_dissipation_j);
        for i in 0..2 {
            let v = fa.tangent_velocities_m_s[i];
            same(fb.tangent_velocities_m_s[1 - i][0], c * v[0] + s * v[1]);
            same(fb.tangent_velocities_m_s[1 - i][1], -s * v[0] + c * v[1]);
        }
    }
}

#[test]
fn weaker_individual_source_changes_coupled_motion_without_replacing_normal_area() {
    let (n0, s0) = sources(1e-5, 2e5, 2e5);
    let (n1, s1) = sources(2e-5, 2e5, 2e5);
    let (nw, sw) = sources(2e-5, 0.0, 0.0);
    let mut a = owner([&n0, &n1], [&s0, &s1], false, true, 8);
    let mut b = owner([&n0, &nw], [&s0, &sw], false, true, 8);
    let p = momentum(&a);
    let fa = *a.step(&[0.0; 6]).unwrap();
    let fb = *b.step(&[0.0; 6]).unwrap();
    assert_ne!(states(&a), states(&b));
    assert!(fa.path.contact_friction_dissipation_j[1] > 0.0);
    assert_eq!(fb.path.contact_friction_dissipation_j[1], 0.0);
    same(
        a.normal_responses()[1].contact_area_m2,
        b.normal_responses()[1].contact_area_m2,
    );
    for axis in 0..2 {
        same(momentum(&a)[axis], p[axis]);
        same(momentum(&b)[axis], p[axis]);
    }
}

#[test]
fn zero_strength_refines_to_the_complete_original_prepared_network() {
    let (n0, s0) = sources(1e-5, 0.0, 0.0);
    let (n1, s1) = sources(2e-5, 0.0, 0.0);
    let external = [2e-5, -3e-6, 1e-6, -1e-5, 7e-6, -2e-6, 4e-6, 3e-6, 2e-6];
    let mut reference = network(true, true, 1.0, 1e-23, 2e-8);
    let expected = reference.step(&external).unwrap().clone();
    let exact: Vec<_> = reference
        .components()
        .iter()
        .flat_map(|b| b.states())
        .collect();
    let mut errors = vec![];
    let mut numerical = vec![];
    for cells in [4, 16, 64] {
        let mut system = owner([&n0, &n1], [&s0, &s1], true, true, cells);
        let frame = *system.step(&external).unwrap();
        assert_eq!(frame.friction_dissipation_j, 0.0);
        let error = states(&system)
            .iter()
            .zip(&exact)
            .map(|(a, b)| {
                (a.velocity_m_sqrt_kg_per_s - b.velocity_m_sqrt_kg_per_s).powi(2)
                    + (2000.0 * (a.displacement_m_sqrt_kg - b.displacement_m_sqrt_kg)).powi(2)
            })
            .sum::<f64>()
            .sqrt();
        errors.push(error);
        numerical.push(frame.numerical_dissipation_j);
        if cells == 64 {
            close(
                frame.observer_pressure_pa,
                expected.observer_pressure_pa,
                1e-12,
                0.01,
            );
        }
    }
    assert!(
        errors[1] < 0.3 * errors[0] && errors[2] < 0.3 * errors[1],
        "{errors:?}"
    );
    assert!(
        numerical[1] < 0.3 * numerical[0] && numerical[2] < 0.3 * numerical[1],
        "{numerical:?}"
    );
}

#[test]
fn strict_native_accuracy_and_numerical_loss_limits_refuse_without_publication() {
    let (n0, s0) = sources(1e-5, 1e5, 1e5);
    let (n1, s1) = sources(2e-5, 2e5, 2e5);
    let mut ordinary = owner([&n0, &n1], [&s0, &s1], true, true, 8);
    let accepted = *ordinary.step(&[0.0; 9]).unwrap();
    assert!(accepted.numerical_work_uncertainty_j > 0.0);
    for defect in 0..2 {
        let mut c = controls(8);
        if defect == 1 {
            c.flow.maximum_numerical_dissipation_j = 0.5 * accepted.numerical_dissipation_j;
        }
        let mut strict = VectorAdhesiveNetworkModalSystem::new(
            network(
                true,
                true,
                1.0,
                if defect == 0 {
                    0.25 * accepted.numerical_work_uncertainty_j
                } else {
                    1e-23
                },
                1e-30,
            ),
            junctions([&n0, &n1], [&s0, &s1], true),
            c,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
        let before = states(&strict);
        let failure = strict.step(&[0.0; 9]).unwrap_err();
        if defect == 0 {
            assert!(
                matches!(
                    failure,
                    VectorAdhesiveNetworkModalError::EnergyAccuracy { .. }
                ),
                "{failure:?}"
            );
        }
        assert_eq!(states(&strict), before);
        assert_eq!(strict.samples_rendered(), 0);
        assert!(strict.last_frame().is_none() && strict.network().last_frame().is_none());
    }
}

#[test]
fn cancellation_and_late_original_body_refusal_leave_identical_retry() {
    let (n0, s0) = sources(1e-5, 2e8, 2e8);
    let (n1, s1) = sources(2e-5, 2e8, 2e8);
    let make = || {
        VectorAdhesiveNetworkModalSystem::new(
            network(true, false, 2e-12, 1e-23, 2e-8),
            junctions([&n0, &n1], [&s0, &s1], true),
            controls(8),
            &CancelGate::new_clock_free(),
        )
        .unwrap()
    };
    let (mut a, mut b) = (make(), make());
    a.step(&[0.0; 9]).unwrap();
    b.step(&[0.0; 9]).unwrap();
    let before = (
        states(&a),
        *a.last_frame().unwrap(),
        a.network().last_frame().cloned(),
    );
    let gate = CancelGate::new_clock_free();
    gate.request();
    assert!(a.step_under_gate(&[0.0; 9], &gate).is_err());
    let force: Vec<_> = MASS
        .iter()
        .flat_map(|m| [100.0 * det::sqrt(*m), 0.0, 0.0])
        .collect();
    let failure = a.step(&force).unwrap_err();
    assert!(
        matches!(
            failure,
            VectorAdhesiveNetworkModalError::Network(ModalCouplingError::Component { .. })
        ),
        "{failure:?}"
    );
    assert_eq!(
        (
            states(&a),
            *a.last_frame().unwrap(),
            a.network().last_frame().cloned()
        ),
        before
    );
    assert_eq!(*a.step(&[0.0; 9]).unwrap(), *b.step(&[0.0; 9]).unwrap());
    assert_eq!(states(&a), states(&b));
}

#[test]
fn source_order_monotone_rung_and_complete_independent_tangent_planes_are_admitted() {
    let (n0, s0) = sources(1e-5, 1e5, 1e5);
    let (n1, s1) = sources(2e-5, 2e5, 2e5);
    let (nv, sv) = sources(2e-5, 2e5, 1e5);
    for defect in 0..5 {
        let mut sites = junctions([&n0, &n1], [&s0, &s1], true);
        match defect {
            0 => sites[1].shear = &s0,
            1 => sites[1].surface_materials.swap(0, 1),
            2 => sites[1].tangent_columns[1] = sites[1].tangent_columns[0].clone(),
            3 => {
                sites[1].tangent_columns[1].pop();
            }
            _ => {
                sites = junctions([&n0, &nv], [&s0, &sv], true);
            }
        }
        assert!(
            VectorAdhesiveNetworkModalSystem::new(
                network(true, false, 1.0, 1e-23, 2e-8),
                sites,
                controls(8),
                &CancelGate::new_clock_free()
            )
            .is_err()
        );
    }
}
