//! G1/G3/G4: two original Hertz patches share a finite body and caloric cell.
//! Independent algebra below uses actual finite masses/inertia and eliminates
//! the three linear caloric equations; it never consumes production matrices.

use super::{A, AREA, CAP, DT, K, LINK, MU, SLOPE, close, energy, law, routes, thermal};
use fs_contact::interface_binding::thermal_friction::SourceThermalCoulomb;
use fs_couple::friction_heat::FrictionHeatRoute;
use fs_couple::friction_heat::hertz::RigidHertzHeatContact;
use fs_couple::friction_heat::hertz::multiple::*;
use fs_exec::CancelGate;
use fs_ga::So3;
use fs_mbd::tangential_friction::RigidContactBinding;
use fs_mbd::tangential_friction::hertz::RigidHertzControls;
use fs_mbd::tangential_friction::hertz::multiple::RigidHertzIslandControls;
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_thermal::enthalpy::backward_euler::{BackwardEulerStep, StepCoefficients};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::sliding_patch::contact::{HertzWrenchConfig, HertzWrenchMode};
use std::cell::Cell;
use std::collections::BTreeMap;

const MASS: [f64; 3] = [2.0, 3.0, 4.0];
const INERTIA: [f64; 3] = [0.08, 0.12, 0.15];
const VELOCITY: [f64; 3] = [4.0, 0.0, 1.0];
const ARM: [f64; 2] = [-A, A];
const LOAD: [f64; 2] = [10.0, 12.0];
const PATCH: [&str; 2] = ["lower-original-Hertz", "upper-original-Hertz"];

#[derive(Clone)]
struct Island {
    bodies: [RigidBodyState; 3],
    properties: [MassProperties; 3],
}
impl Island {
    fn new(moving: bool) -> Self {
        Self {
            properties: std::array::from_fn(|i| {
                MassProperties::new(
                    MASS[i],
                    Vec3::ZERO,
                    Vec3::new(INERTIA[i], INERTIA[i], INERTIA[i]),
                )
                .unwrap()
            }),
            bodies: std::array::from_fn(|i| {
                RigidBodyState::new(
                    Pose::new(
                        Vec3::new(0.0, 0.0, [0.0, -2.0 * A, 2.0 * A][i]),
                        UnitQuaternion::IDENTITY,
                    )
                    .unwrap(),
                    Vec3::new(if moving { MASS[i] * VELOCITY[i] } else { 0.0 }, 0.0, 0.0),
                    Vec3::ZERO,
                )
                .unwrap()
            }),
        }
    }
    fn contacts<'a>(
        &self,
        sources: &'a [SourceThermalCoulomb; 2],
        routes: &'a [FrictionHeatRoute; 2],
    ) -> [RigidHertzIslandHeatContact<'a>; 2] {
        std::array::from_fn(|i| RigidHertzIslandHeatContact {
            patch_id: PATCH[i],
            contact: RigidHertzHeatContact {
                binding: RigidContactBinding {
                    body_a: 0,
                    body_b: i + 1,
                    arm_a_body: Vec3::new(0.0, 0.0, ARM[i]),
                    arm_b_body: Vec3::new(0.0, 0.0, -ARM[i]),
                },
                frame: TangentFrame::new([0.0, 0.0, [1.0, -1.0][i]], [1.0, 0.0, 0.0]).unwrap(),
                source: &sources[i],
                normal_force_n: LOAD[i],
                contact_radius_m: (AREA / core::f64::consts::PI).sqrt(),
                nominal_pressure_pa: None,
            },
            route: &routes[i],
            temperature_guess_k: 300.0,
        })
    }
}

fn controls() -> RigidHertzIslandHeatControls {
    let shared = super::controls();
    RigidHertzIslandHeatControls {
        mechanical: RigidHertzIslandControls {
            mechanical: RigidHertzControls {
                maximum_bodies: 3,
                geometry_tolerance_m: 1e-10,
                velocity_tolerance_m_s: 1e-10,
                angular_velocity_tolerance_rad_s: 1e-9,
                linear_momentum_tolerance_ns: 1e-10,
                angular_momentum_tolerance_nms: 1e-10,
                energy_tolerance_j: 1e-10,
                contact: HertzWrenchConfig {
                    maximum_static_iterations: 96,
                    maximum_iterations: 256,
                    maximum_evaluations: 512,
                    force_tolerance_n: 2e-9,
                    moment_tolerance_nm: 2e-10,
                    displacement_tolerance_m: 2e-13,
                    energy_tolerance_j: 2e-11,
                    relative_matrix_tolerance: 1e-13,
                    maximum_force_n: 100.0,
                    maximum_spin_moment_nm: 10.0,
                },
            },
            maximum_contacts: 2,
            maximum_matrix_terms: 36,
            maximum_sweeps: 96,
            maximum_evaluations: 8192,
        },
        maximum_cells: 3,
        maximum_unknowns: 5,
        maximum_trials: shared.maximum_trials,
        maximum_newton_iterations: shared.maximum_newton_iterations,
        difference_step_j: shared.difference_step_j,
        temperature_scale_j_per_k: shared.temperature_scale_j_per_k,
        heat: shared.heat,
        combined_energy_tolerance_j: shared.combined_energy_tolerance_j,
        newton: shared.newton,
        audit_entropy: true,
    }
}
fn source(
    _: &DepositionSimulation<'_>,
    ids: &[u64],
    temperatures: &[f64],
) -> Result<RigidHertzIslandHeatSource<Vec<f64>>, &'static str> {
    assert_eq!(ids, &[10, 20]);
    Ok(RigidHertzIslandHeatSource {
        interfaces: ids
            .iter()
            .map(|&link_id| RigidHertzIslandHeatConductivity {
                link_id,
                conductivity_w_m_k: [K[0], K[if link_id == 20 { 1 } else { 2 }]],
            })
            .collect(),
        transport: StepCoefficients {
            sources_w: BTreeMap::new(),
            conductances_w_k: BTreeMap::new(),
            payload: temperatures.to_vec(),
        },
    })
}
fn run<'m>(
    island: &Island,
    contacts: &[RigidHertzIslandHeatContact<'_>],
    thermal: &DepositionSimulation<'m>,
    dt: f64,
    controls: RigidHertzIslandHeatControls,
    resolve: impl Fn(
        &DepositionSimulation<'m>,
        &[u64],
        &[f64],
    ) -> Result<RigidHertzIslandHeatSource<Vec<f64>>, &'static str>,
    gate: &CancelGate,
) -> Result<HeatedRigidHertzIslandCandidate<'m, Vec<f64>>, RigidHertzIslandHeatError<&'static str>>
{
    propose_rigid_hertz_island_heat(
        &island.bodies,
        &island.properties,
        contacts,
        thermal,
        thermal.time_s(),
        dt,
        controls,
        &[0.0; 3],
        resolve,
        gate,
    )
}

fn mobility() -> [[f64; 2]; 2] {
    std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            1.0 / MASS[0]
                + ARM[i] * ARM[j] / INERTIA[0]
                + if i == j {
                    1.0 / MASS[i + 1] + ARM[i] * ARM[i] / INERTIA[i + 1]
                } else {
                    0.0
                }
        })
    })
}
// Solve the independent 2x2 face-temperature root after analytically eliminating
// the three original finite cells. Native forces are the exact pure-slip limit.
fn reference(
    initial: [f64; 3],
    outside_w: [f64; 3],
    dt: f64,
    moving: bool,
) -> ([f64; 2], [f64; 3], [f64; 2], [f64; 2]) {
    let d = K.map(|k| dt * k * AREA / A);
    let base: [f64; 3] = std::array::from_fn(|i| CAP[i] * (initial[i] - 300.0) + dt * outside_w[i]);
    let denom = [CAP[0] + 2.0 * d[0], CAP[1] + d[1], CAP[2] + d[2]];
    let cross = -d[0] * d[0] / denom[0];
    let matrix: [[f64; 2]; 2] = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            cross
                + if i == j {
                    d[0] + d[i + 1] - d[i + 1] * d[i + 1] / denom[i + 1]
                } else {
                    0.0
                }
        })
    });
    let rhs: [f64; 2] =
        std::array::from_fn(|i| d[0] * base[0] / denom[0] + d[i + 1] * base[i + 1] / denom[i + 1]);
    let w = mobility();
    let slope: [f64; 2] = std::array::from_fn(|i| LOAD[i] * SLOPE[i]);
    let mut theta = [0.0; 2];
    for _ in 0..20 {
        let effort: [f64; 2] = std::array::from_fn(|i| {
            if moving {
                LOAD[i] * (MU[i] - SLOPE[i] * theta[i])
            } else {
                0.0
            }
        });
        let slip: [f64; 2] = std::array::from_fn(|i| {
            if moving {
                dt * (VELOCITY[0] - VELOCITY[i + 1])
                    - 0.5 * dt * dt * (w[i][0] * effort[0] + w[i][1] * effort[1])
            } else {
                0.0
            }
        });
        let residual: [f64; 2] = std::array::from_fn(|i| {
            matrix[i][0] * theta[0] + matrix[i][1] * theta[1] - rhs[i] - effort[i] * slip[i]
        });
        let jac: [[f64; 2]; 2] = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                matrix[i][j]
                    - if moving {
                        -if i == j { slope[i] * slip[i] } else { 0.0 }
                            + 0.5 * dt * dt * effort[i] * w[i][j] * slope[j]
                    } else {
                        0.0
                    }
            })
        });
        let det = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];
        theta[0] -= (jac[1][1] * residual[0] - jac[0][1] * residual[1]) / det;
        theta[1] -= (jac[0][0] * residual[1] - jac[1][0] * residual[0]) / det;
    }
    let cells = [
        300.0 + (base[0] + d[0] * (theta[0] + theta[1])) / denom[0],
        300.0 + (base[1] + d[1] * theta[0]) / denom[1],
        300.0 + (base[2] + d[2] * theta[1]) / denom[2],
    ];
    let force: [f64; 2] = std::array::from_fn(|i| {
        if moving {
            -LOAD[i] * (MU[i] - SLOPE[i] * theta[i])
        } else {
            0.0
        }
    });
    let heat = std::array::from_fn(|i| {
        -force[i]
            * (dt * (VELOCITY[0] - VELOCITY[i + 1])
                + 0.5 * dt * dt * (w[i][0] * force[0] + w[i][1] * force[1]))
    });
    (theta.map(|t| t + 300.0), cells, heat, force)
}
fn kinetic(bodies: &[RigidBodyState]) -> f64 {
    bodies
        .iter()
        .enumerate()
        .map(|(i, b)| {
            let p = b.linear_momentum_world();
            let l = b.angular_momentum_body();
            0.5 * (p.dot(p) / MASS[i] + l.dot(l) / INERTIA[i])
        })
        .sum()
}
fn compare(
    candidate: &HeatedRigidHertzIslandCandidate<'_, Vec<f64>>,
    island: &Island,
    accepted: &DepositionSimulation<'_>,
    initial: [f64; 3],
    outside: [f64; 3],
    dt: f64,
    moving: bool,
) {
    let (face, cells, heat, force) = reference(initial, outside, dt, moving);
    assert_eq!(candidate.mechanical.body_indices, [0, 1, 2]);
    for i in 0..2 {
        let native = &candidate.mechanical.report.contacts[i];
        assert_eq!(native.patch_id, PATCH[i]);
        close(native.contact.force_n[0], force[i]);
        close(native.contact.force_n[1], 0.0);
        close(native.contact.spin_moment_nm, 0.0);
        close(native.dissipated_energy_j, heat[i]);
        let report = candidate
            .heat
            .interfaces
            .iter()
            .find(|r| r.link_id == LINK[i])
            .unwrap();
        close(report.interface_temperature_k, face[i]);
        close(report.dissipated_j, heat[i]);
        close(report.body_heat_j[0] + report.body_heat_j[1], heat[i]);
        let receipt = candidate
            .friction_sources
            .iter()
            .find(|r| r.patch_id == PATCH[i])
            .unwrap();
        close(
            receipt
                .friction_source
                .friction_source()
                .state_coordinate("T")
                .unwrap(),
            face[i],
        );
        close(
            receipt.friction_source.law().parameters().kinetic_mu,
            MU[i] - SLOPE[i] * (face[i] - 300.0),
        );
    }
    for i in 0..3 {
        close(
            candidate.thermal.temperature_k(CellId(i as u64)).unwrap(),
            cells[i],
        );
        let b = &candidate.mechanical.bodies_after[i];
        let p = island.bodies[i].linear_momentum_world();
        close(
            b.linear_momentum_world().x,
            p.x + dt
                * if i == 0 {
                    force[0] + force[1]
                } else {
                    -force[i - 1]
                },
        );
        close(b.linear_momentum_world().y, p.y);
        close(
            b.angular_momentum_body().y,
            dt * if i == 0 {
                ARM[0] * force[0] + ARM[1] * force[1]
            } else {
                ARM[i - 1] * force[i - 1]
            },
        );
    }
    let outside_j = dt * outside.iter().sum::<f64>();
    close(candidate.heat.external_heat_j, outside_j);
    close(
        energy(&candidate.thermal) - energy(accepted) - outside_j,
        heat.iter().sum(),
    );
    close(
        kinetic(&candidate.mechanical.bodies_after) - kinetic(&island.bodies)
            + heat.iter().sum::<f64>(),
        0.0,
    );
    close(
        kinetic(&candidate.mechanical.bodies_after) - kinetic(&island.bodies)
            + energy(&candidate.thermal)
            - energy(accepted)
            - outside_j,
        0.0,
    );
    close(candidate.combined_energy_residual_j, 0.0);
    close(candidate.summed_energy_residual_j, 0.0);
    assert_eq!(candidate.thermal.time_s(), accepted.time_s() + dt);
}

#[test]
fn original_shared_body_and_cell_match_independent_joint_friction_thermal_root() {
    let island = Island::new(true);
    let saved = island.bodies;
    let sources = [law(0), law(1)];
    let original = sources.clone();
    let initial = [302.0, 298.0, 304.0];
    let thermal = thermal(initial, So3::identity());
    let routes = routes(&thermal);
    let candidate = run(
        &island,
        &island.contacts(&sources, &routes),
        &thermal,
        DT,
        controls(),
        source,
        &CancelGate::new(),
    )
    .unwrap();
    compare(&candidate, &island, &thermal, initial, [0.0; 3], DT, true);
    assert!(candidate.nonlinear.converged);
    assert!(
        candidate
            .mechanical
            .report
            .contacts
            .iter()
            .all(|p| p.contact.mode == HertzWrenchMode::Sliding)
    );
    for (i, source) in sources.iter().enumerate() {
        let receipt = &candidate
            .friction_sources
            .iter()
            .find(|r| r.patch_id == PATCH[i])
            .unwrap()
            .friction_source;
        assert_eq!(receipt.owner_identity(), source.identity());
        for property in receipt.friction_source().properties() {
            let initial = source
                .initial_source()
                .property(property.requirement().name())
                .unwrap();
            assert_eq!(
                property.answer().receipt.selected,
                initial.answer().receipt.selected
            );
            source
                .interface()
                .card()
                .claims()
                .verify_receipt(&property.answer().receipt)
                .unwrap();
        }
    }
    let entropy = candidate.heat.entropy.unwrap();
    close(
        entropy.friction_generation_j_k,
        candidate
            .heat
            .interfaces
            .iter()
            .map(|r| r.dissipated_j / r.interface_temperature_k)
            .sum(),
    );
    assert!(entropy.interface_conduction_generation_j_k > 0.0);
    assert!(!entropy.has_unclassified_sources);
    assert_eq!(entropy.network_conduction_generation_j_k, 0.0); // Both links replaced once.
    assert_eq!(island.bodies, saved);
    assert_eq!(sources, original);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn original_patch_permutation_is_exact_and_signed_outside_heat_stays_separate() {
    let island = Island::new(true);
    let sources = [law(0), law(1)];
    let initial = [302.0, 298.0, 304.0];
    let thermal = thermal(initial, So3::identity());
    let routes = routes(&thermal);
    let contacts = island.contacts(&sources, &routes);
    let baseline = run(
        &island,
        &contacts,
        &thermal,
        DT,
        controls(),
        source,
        &CancelGate::new(),
    )
    .unwrap();
    let permuted = run(
        &island,
        &[contacts[1], contacts[0]],
        &thermal,
        DT,
        controls(),
        source,
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(baseline.mechanical, permuted.mechanical);
    assert_eq!(baseline.friction_sources, permuted.friction_sources);
    assert_eq!(baseline.source, permuted.source);
    for i in 0..3 {
        assert_eq!(
            baseline.thermal.temperature_k(CellId(i)).unwrap(),
            permuted.thermal.temperature_k(CellId(i)).unwrap()
        );
    }
    let outside = [1.0, -2.0, 0.4];
    let heated = run(
        &island,
        &contacts,
        &thermal,
        DT,
        controls(),
        |s, ids, t| {
            let mut result = source(s, ids, t)?;
            result.transport.sources_w = outside
                .iter()
                .enumerate()
                .map(|(i, &p)| (CellId(i as u64), p))
                .collect();
            Ok(result)
        },
        &CancelGate::new(),
    )
    .unwrap();
    compare(&heated, &island, &thermal, initial, outside, DT, true);
    assert!(heated.heat.entropy.unwrap().has_unclassified_sources);
}

#[test]
fn stationary_island_keeps_both_ordinary_conduction_paths_without_friction_heat() {
    let island = Island::new(false);
    let sources = [law(0), law(1)];
    let initial = [310.0, 290.0, 305.0];
    let thermal = thermal(initial, So3::identity());
    let routes = routes(&thermal);
    let candidate = run(
        &island,
        &island.contacts(&sources, &routes),
        &thermal,
        DT,
        controls(),
        source,
        &CancelGate::new(),
    )
    .unwrap();
    compare(&candidate, &island, &thermal, initial, [0.0; 3], DT, false);
    assert_eq!(candidate.mechanical.bodies_after, island.bodies);
    assert_eq!(candidate.heat.dissipated_j, 0.0);
    // The shared middle cell cools below the upper body at the implicit
    // endpoint. Actual endpoint temperatures, not their initial ranking,
    // therefore send heat into the lower body and out of the upper body.
    let lower = candidate
        .heat
        .interfaces
        .iter()
        .find(|r| r.link_id == LINK[0])
        .unwrap();
    let upper = candidate
        .heat
        .interfaces
        .iter()
        .find(|r| r.link_id == LINK[1])
        .unwrap();
    assert!(lower.body_heat_j[0] < 0.0 && lower.body_heat_j[1] > 0.0);
    assert!(upper.body_heat_j[0] > 0.0 && upper.body_heat_j[1] < 0.0);
    close(candidate.heat.entropy.unwrap().friction_generation_j_k, 0.0);
}

#[test]
fn original_routing_geometry_source_domains_and_boundaries_refuse_without_publication() {
    let island = Island::new(true);
    let saved = island.bodies;
    let sources = [law(0), law(1)];
    let thermal = thermal([300.0; 3], So3::identity());
    let routes = routes(&thermal);
    let contacts = island.contacts(&sources, &routes);
    let attempt = |c: &[RigidHertzIslandHeatContact<'_>], ctrl| {
        run(&island, c, &thermal, DT, ctrl, source, &CancelGate::new())
    };
    assert!(attempt(&[], controls()).is_err());
    let mut changed = contacts;
    changed[1].patch_id = PATCH[0];
    assert!(attempt(&changed, controls()).is_err());
    let mut changed = contacts;
    changed[1].route = &routes[0];
    assert!(attempt(&changed, controls()).is_err());
    let mut changed = contacts;
    changed[1].temperature_guess_k = 330.0;
    assert!(matches!(
        attempt(&changed, controls()),
        Err(RigidHertzIslandHeatError::FrictionSource(_))
    ));
    let mut changed = contacts;
    changed[1].contact.contact_radius_m *= 1.1;
    assert!(matches!(
        attempt(&changed, controls()),
        Err(RigidHertzIslandHeatError::Invalid(_))
    ));
    let mut changed = contacts;
    changed[1].contact.nominal_pressure_pa = Some(1e5);
    assert!(matches!(
        attempt(&changed, controls()),
        Err(RigidHertzIslandHeatError::FrictionSource(_))
    ));
    for which in 0..3 {
        let mut c = controls();
        match which {
            0 => c.maximum_unknowns = 4,
            1 => c.maximum_cells = 2,
            _ => c.mechanical.maximum_contacts = 1,
        }
        assert!(attempt(&contacts, c).is_err());
    }
    let bad_order = run(
        &island,
        &contacts,
        &thermal,
        DT,
        controls(),
        |s, ids, t| {
            let mut out = source(s, ids, t)?;
            out.interfaces.swap(0, 1);
            Ok(out)
        },
        &CancelGate::new(),
    );
    assert!(matches!(
        bad_order,
        Err(RigidHertzIslandHeatError::Invalid(_))
    ));
    let duplicate = run(
        &island,
        &contacts,
        &thermal,
        DT,
        controls(),
        |s, ids, t| {
            let mut out = source(s, ids, t)?;
            out.transport.conductances_w_k.insert(10, 1.0);
            Ok(out)
        },
        &CancelGate::new(),
    );
    assert!(matches!(
        duplicate,
        Err(RigidHertzIslandHeatError::Thermal(_))
    ));
    assert_eq!(island.bodies, saved);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(energy(&thermal), 0.0);
}

#[test]
fn loaded_late_source_failure_trial_budget_and_cancellation_retry_exactly() {
    let island = Island::new(true);
    let saved = island.bodies;
    let sources = [law(0), law(1)];
    let thermal = thermal([300.0; 3], So3::identity());
    let routes = routes(&thermal);
    let contacts = island.contacts(&sources, &routes);
    let baseline = run(
        &island,
        &contacts,
        &thermal,
        DT,
        controls(),
        source,
        &CancelGate::new(),
    )
    .unwrap();
    let calls = Cell::new(0);
    let late = run(
        &island,
        &contacts,
        &thermal,
        DT,
        controls(),
        |s, ids, t| {
            calls.set(calls.get() + 1);
            if calls.get() == baseline.trials {
                Err("last accepted-source evaluation")
            } else {
                source(s, ids, t)
            }
        },
        &CancelGate::new(),
    );
    assert!(matches!(late, Err(RigidHertzIslandHeatError::Source(_))));
    assert_eq!(calls.get(), baseline.trials);
    let mut c = controls();
    c.maximum_trials = baseline.trials - 1;
    assert!(matches!(
        run(
            &island,
            &contacts,
            &thermal,
            DT,
            c,
            source,
            &CancelGate::new()
        ),
        Err(RigidHertzIslandHeatError::TrialBudget)
    ));
    let gate = CancelGate::new();
    let calls = Cell::new(0);
    let cancelled = run(
        &island,
        &contacts,
        &thermal,
        DT,
        controls(),
        |s, ids, t| {
            calls.set(calls.get() + 1);
            if calls.get() == baseline.trials {
                gate.request();
            }
            source(s, ids, t)
        },
        &gate,
    );
    assert!(matches!(
        cancelled,
        Err(RigidHertzIslandHeatError::Thermal(
            DepositionError::Cancelled
        ))
    ));
    let retry = run(
        &island,
        &contacts,
        &thermal,
        DT,
        controls(),
        source,
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(retry.mechanical, baseline.mechanical);
    assert_eq!(retry.source, baseline.source);
    assert_eq!(retry.friction_sources, baseline.friction_sources);
    assert_eq!(island.bodies, saved);
    assert_eq!(energy(&thermal), 0.0);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn large_accepted_clock_uses_identical_representable_duration_for_all_patches_and_cells() {
    let island = Island::new(true);
    let sources = [law(0), law(1)];
    let original = thermal([300.0; 3], So3::identity());
    let clock = 2.0_f64.powi(40);
    let accepted = BackwardEulerStep::new(&original, clock, || Ok(()))
        .unwrap()
        .evaluate::<(), DepositionError>(
            &[0.0; 3],
            |_| {
                Ok(StepCoefficients {
                    sources_w: BTreeMap::new(),
                    conductances_w_k: BTreeMap::new(),
                    payload: (),
                })
            },
            || Ok(()),
        )
        .unwrap()
        .accept(0.0, 0.0, || Ok(()))
        .unwrap()
        .0;
    let routes = routes(&accepted);
    let contacts = island.contacts(&sources, &routes);
    let candidate = run(
        &island,
        &contacts,
        &accepted,
        DT,
        controls(),
        source,
        &CancelGate::new(),
    )
    .unwrap();
    let actual = (clock + DT) - clock;
    assert_ne!(actual.to_bits(), DT.to_bits());
    compare(
        &candidate, &island, &accepted, [300.0; 3], [0.0; 3], actual, true,
    );
    assert_eq!(candidate.target_time_s, clock + DT);
    assert_eq!(accepted.time_s(), clock);
}
