//! Source-resolved Hertz recoil and the original three-cell thermal fixture.

use super::{A, AREA, CAP, DT, K, MU, SLOPE, close, energy, law, routes, solve_dense, thermal};
use fs_contact::interface_binding::thermal_friction::SourceThermalCoulomb;
use fs_couple::friction_heat::FrictionHeatRoute;
use fs_couple::friction_heat::hertz::*;
use fs_exec::CancelGate;
use fs_ga::So3;
use fs_mbd::tangential_friction::RigidContactBinding;
use fs_mbd::tangential_friction::hertz::RigidHertzControls;
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_thermal::enthalpy::backward_euler::{BackwardEulerStep, StepCoefficients};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::sliding_patch::contact::{HertzWrenchConfig, HertzWrenchMode};
use std::cell::Cell;
use std::collections::BTreeMap;

const MASS: [f64; 2] = [2.0, 3.0];
const I: [[f64; 3]; 2] = [[0.08, 0.10, 0.12], [0.15, 0.20, 0.25]];
const ARMS: [f64; 2] = [-0.2, 0.15];
const SPEED: f64 = 4.0;
const SPIN: f64 = 100.0;

fn radius() -> f64 {
    (AREA / core::f64::consts::PI).sqrt()
}
fn mobility(spin: bool) -> f64 {
    if spin {
        1.0 / I[0][2] + 1.0 / I[1][2]
    } else {
        1.0 / MASS[0] + 1.0 / MASS[1] + ARMS[0] * ARMS[0] / I[0][1] + ARMS[1] * ARMS[1] / I[1][1]
    }
}
fn effort_scale(spin: bool) -> f64 {
    if spin {
        3.0 * core::f64::consts::PI * 10.0 * radius() / 16.0
    } else {
        10.0
    }
}

#[derive(Clone, Copy)]
struct Pair {
    bodies: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    binding: RigidContactBinding,
}
impl Pair {
    fn new(speed: f64, spin: f64) -> Self {
        let properties = std::array::from_fn(|i| {
            MassProperties::new(MASS[i], Vec3::ZERO, Vec3::new(I[i][0], I[i][1], I[i][2])).unwrap()
        });
        let bodies = std::array::from_fn(|i| {
            RigidBodyState::new(
                Pose::new(Vec3::new(0.0, 0.0, -A - ARMS[i]), UnitQuaternion::IDENTITY).unwrap(),
                Vec3::new(if i == 0 { MASS[i] * speed } else { 0.0 }, 0.0, 0.0),
                Vec3::new(0.0, 0.0, if i == 0 { I[i][2] * spin } else { 0.0 }),
            )
            .unwrap()
        });
        Self {
            bodies,
            properties,
            binding: RigidContactBinding {
                body_a: 0,
                body_b: 1,
                arm_a_body: Vec3::new(0.0, 0.0, ARMS[0]),
                arm_b_body: Vec3::new(0.0, 0.0, ARMS[1]),
            },
        }
    }
    fn contact<'a>(&self, source: &'a SourceThermalCoulomb) -> RigidHertzHeatContact<'a> {
        RigidHertzHeatContact {
            binding: self.binding,
            frame: TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap(),
            source,
            normal_force_n: 10.0,
            contact_radius_m: radius(),
            nominal_pressure_pa: None,
        }
    }
}

fn controls() -> RigidHertzHeatControls {
    let shared = super::controls();
    RigidHertzHeatControls {
        mechanical: RigidHertzControls {
            maximum_bodies: 2,
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
        maximum_cells: shared.maximum_cells,
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
    temperature: f64,
) -> Result<RigidHertzHeatSource<f64>, &'static str> {
    Ok(RigidHertzHeatSource {
        conductivity_w_m_k: [K[0], K[1]],
        transport: StepCoefficients {
            sources_w: BTreeMap::new(),
            conductances_w_k: BTreeMap::new(),
            payload: temperature,
        },
    })
}
fn run<'m>(
    pair: &Pair,
    owner: &SourceThermalCoulomb,
    thermal: &DepositionSimulation<'m>,
    route: &FrictionHeatRoute,
    controls: RigidHertzHeatControls,
    resolve: impl Fn(&DepositionSimulation<'m>, f64) -> Result<RigidHertzHeatSource<f64>, &'static str>,
    gate: &CancelGate,
) -> Result<HeatedRigidHertzCandidate<'m, f64>, RigidHertzHeatError<&'static str>> {
    propose_rigid_hertz_heat(
        &pair.bodies,
        &pair.properties,
        pair.contact(owner),
        thermal,
        route,
        0.0,
        DT,
        controls,
        &[0.0; 3],
        300.0,
        resolve,
        gate,
    )
}

// Independent scalar physical root. Cell 2 remains on the original direct
// series link 10, so eliminating it changes the capacity seen by cell 0.
fn reference(spin: bool, duration_s: f64) -> (f64, [f64; 3], f64, f64) {
    let d0 = duration_s * K[0] * AREA / A;
    let d1 = duration_s * K[1] * AREA / A;
    let d02 = duration_s * AREA / (A / K[0] + A / K[2]);
    let c0 = CAP[0] + d02 * CAP[2] / (CAP[2] + d02);
    let thermal_response = d0 * c0 / (c0 + d0) + d1 * CAP[1] / (CAP[1] + d1);
    let heat = |theta: f64| {
        let effort = effort_scale(spin) * (MU[0] - SLOPE[0] * theta);
        effort
            * (duration_s * if spin { SPIN } else { SPEED }
                - 0.5 * duration_s * duration_s * mobility(spin) * effort)
    };
    let mut lower = 0.0;
    let mut upper = 20.0;
    assert!(thermal_response * upper > heat(upper));
    for _ in 0..80 {
        let mid = 0.5 * (lower + upper);
        if thermal_response * mid < heat(mid) {
            lower = mid;
        } else {
            upper = mid;
        }
    }
    let theta = 0.5 * (lower + upper);
    let t0 = d0 * theta / (c0 + d0);
    (
        300.0 + theta,
        [
            300.0 + t0,
            300.0 + d1 * theta / (CAP[1] + d1),
            300.0 + d02 * t0 / (CAP[2] + d02),
        ],
        heat(theta),
        effort_scale(spin) * (MU[0] - SLOPE[0] * theta),
    )
}
fn kinetic(bodies: &[RigidBodyState; 2]) -> f64 {
    (0..2)
        .map(|i| {
            let p = bodies[i].linear_momentum_world();
            let l = bodies[i].angular_momentum_body();
            0.5 * (p.dot(p) / MASS[i]
                + l.x * l.x / I[i][0]
                + l.y * l.y / I[i][1]
                + l.z * l.z / I[i][2])
        })
        .sum()
}
fn audit(
    candidate: &HeatedRigidHertzCandidate<'_, f64>,
    pair: &Pair,
    old: &DepositionSimulation<'_>,
) {
    assert_eq!(candidate.heat.interfaces.len(), 1);
    assert_eq!(candidate.heat.interfaces[0].link_id, 20);
    close(candidate.combined_energy_residual_j, 0.0);
    close(candidate.mechanical.total_energy_residual_j, 0.0);
    close(candidate.heat.total_energy_residual_j, 0.0);
    close(
        kinetic(&candidate.bodies_after) - kinetic(&pair.bodies) + energy(&candidate.thermal)
            - energy(old)
            - candidate.heat.external_heat_j,
        0.0,
    );
    close(
        candidate.heat.dissipated_j,
        candidate.mechanical.dissipated_energy_j,
    );
    close(
        candidate.heat.interfaces[0].body_heat_j.iter().sum(),
        candidate.heat.dissipated_j,
    );
    assert_eq!(
        candidate.source,
        candidate.heat.interfaces[0].interface_temperature_k
    );
    assert_eq!(candidate.target_time_s, DT);
    assert_eq!(candidate.thermal.time_s(), DT);
    assert_eq!(old.time_s(), 0.0);
}
fn sliding_or_spin_reference(spin: bool) {
    let pair = Pair::new(
        if spin { 0.0 } else { SPEED },
        if spin { SPIN } else { 0.0 },
    );
    let saved = pair.bodies;
    let owner = law(0);
    let original = owner.clone();
    let thermal = thermal([300.0; 3], So3::identity());
    let routes = routes(&thermal);
    let candidate = run(
        &pair,
        &owner,
        &thermal,
        &routes[0],
        controls(),
        source,
        &CancelGate::new(),
    )
    .unwrap();
    let (interface, cells, heat, effort) = reference(spin, DT);
    close(
        candidate.heat.interfaces[0].interface_temperature_k,
        interface,
    );
    close(candidate.mechanical.dissipated_energy_j, heat);
    assert!(effort < effort_scale(spin) * MU[0]);
    assert!(candidate.trials > 2);
    assert!(candidate.nonlinear.converged);
    assert_eq!(candidate.mechanical.contact.mode, HertzWrenchMode::Sliding);
    for i in 0..3 {
        close(
            candidate.thermal.temperature_k(CellId(i as u64)).unwrap(),
            cells[i],
        );
    }
    assert!(cells[2] > 300.0); // Original ordinary link 10 still conducts.
    if spin {
        close(candidate.mechanical.contact.spin_moment_nm, -effort);
        assert_eq!(candidate.mechanical.contact.force_n, [0.0; 2]);
        close(
            candidate.mechanical.spin_after_rad_s,
            SPIN - DT * mobility(true) * effort,
        );
        close(
            candidate.bodies_after[0].angular_momentum_body().z,
            I[0][2] * SPIN - DT * effort,
        );
        close(
            candidate.bodies_after[1].angular_momentum_body().z,
            DT * effort,
        );
        close(
            candidate.mechanical.scaled_mobility_per_kg[2][2],
            radius() * radius() * mobility(true),
        );
    } else {
        close(candidate.mechanical.contact.force_n[0], -effort);
        assert_eq!(candidate.mechanical.contact.spin_moment_nm, 0.0);
        close(
            candidate.mechanical.slip_velocity_after_m_s[0],
            SPEED - DT * mobility(false) * effort,
        );
        close(
            candidate.bodies_after[0].linear_momentum_world().x,
            MASS[0] * SPEED - DT * effort,
        );
        close(
            candidate.bodies_after[1].linear_momentum_world().x,
            DT * effort,
        );
        close(
            candidate.bodies_after[0].angular_momentum_body().y,
            -ARMS[0] * DT * effort,
        );
        close(
            candidate.bodies_after[1].angular_momentum_body().y,
            ARMS[1] * DT * effort,
        );
        close(
            candidate.mechanical.scaled_mobility_per_kg[0][0],
            mobility(false),
        );
    }
    assert_eq!(candidate.friction_source.owner_identity(), owner.identity());
    assert_eq!(
        candidate
            .friction_source
            .friction_source()
            .state_coordinate("T"),
        Some(candidate.source)
    );
    close(
        candidate.friction_source.law().parameters().kinetic_mu,
        MU[0] - SLOPE[0] * (interface - 300.0),
    );
    for property in candidate.friction_source.friction_source().properties() {
        let initial = owner
            .initial_source()
            .property(property.requirement().name())
            .unwrap();
        assert_eq!(
            property.answer().receipt.selected,
            initial.answer().receipt.selected
        );
        owner
            .interface()
            .card()
            .claims()
            .verify_receipt(&property.answer().receipt)
            .unwrap();
    }
    let entropy = candidate.heat.entropy.unwrap();
    assert!(entropy.friction_generation_j_k > 0.0);
    assert!(entropy.network_conduction_generation_j_k > 0.0);
    assert!(entropy.interface_conduction_generation_j_k > 0.0);
    assert!(!entropy.has_unclassified_sources);
    audit(&candidate, &pair, &thermal);
    assert_eq!(pair.bodies, saved);
    assert_eq!(owner, original);
}

#[test]
fn original_temperature_claims_close_translation_recoil_and_finite_cell_heat() {
    sliding_or_spin_reference(false);
}
#[test]
fn original_temperature_claims_close_spin_recoil_and_actual_torsional_heat() {
    sliding_or_spin_reference(true);
}

#[test]
fn static_hertz_contact_preserves_ordinary_conduction_and_exact_caloric_entropy() {
    let pair = Pair::new(0.0, 0.0);
    let owner = law(0);
    let initial = [310.0, 290.0, 300.0];
    let thermal = thermal(initial, So3::identity());
    let routes = routes(&thermal);
    let candidate = run(
        &pair,
        &owner,
        &thermal,
        &routes[0],
        controls(),
        source,
        &CancelGate::new(),
    )
    .unwrap();
    let g = [1, 2].map(|i| DT * AREA / (A / K[0] + A / K[i]));
    let reference = solve_dense(
        vec![
            vec![CAP[0] + g[0] + g[1], -g[0], -g[1]],
            vec![-g[0], CAP[1] + g[0], 0.0],
            vec![-g[1], 0.0, CAP[2] + g[1]],
        ],
        (0..3).map(|i| CAP[i] * initial[i]).collect(),
    );
    for i in 0..3 {
        close(
            candidate.thermal.temperature_k(CellId(i as u64)).unwrap(),
            reference[i],
        );
    }
    let heat = &candidate.heat.interfaces[0];
    close(
        heat.interface_temperature_k,
        (K[0] * reference[0] + K[1] * reference[1]) / (K[0] + K[1]),
    );
    assert!(heat.body_heat_j[0] < 0.0);
    assert!(heat.body_heat_j[1] > 0.0);
    assert_eq!(candidate.mechanical.contact.mode, HertzWrenchMode::Sticking);
    assert_eq!(candidate.mechanical.dissipated_energy_j, 0.0);
    assert_eq!(candidate.bodies_after, pair.bodies);
    let entropy = candidate.heat.entropy.unwrap();
    assert_eq!(entropy.friction_generation_j_k, 0.0);
    assert!(entropy.network_conduction_generation_j_k > 0.0);
    assert!(entropy.interface_conduction_generation_j_k > 0.0);
    close(
        entropy.thermal_storage_change_j_k,
        (0..3)
            .map(|i| CAP[i] * (reference[i] / initial[i]).ln())
            .sum(),
    );
    audit(&candidate, &pair, &thermal);
}

#[test]
fn physical_face_radius_clock_source_domain_and_competing_conduction_refuse() {
    let pair = Pair::new(SPEED, 0.0);
    let saved = pair.bodies;
    let owner = law(0);
    let thermal = thermal([300.0; 3], So3::identity());
    let routes = routes(&thermal);
    let contact = pair.contact(&owner);
    let attempt = |bodies: &[RigidBodyState], contact: RigidHertzHeatContact<'_>, clock, guess| {
        propose_rigid_hertz_heat(
            bodies,
            &pair.properties,
            contact,
            &thermal,
            &routes[0],
            clock,
            DT,
            controls(),
            &[0.0; 3],
            guess,
            source,
            &CancelGate::new(),
        )
    };
    let wrong_radius = RigidHertzHeatContact {
        contact_radius_m: 1.01 * radius(),
        ..contact
    };
    assert!(matches!(
        attempt(&pair.bodies, wrong_radius, 0.0, 300.0),
        Err(RigidHertzHeatError::Invalid(_))
    ));
    let wrong_normal = RigidHertzHeatContact {
        frame: TangentFrame::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).unwrap(),
        ..contact
    };
    assert!(matches!(
        attempt(&pair.bodies, wrong_normal, 0.0, 300.0),
        Err(RigidHertzHeatError::Invalid(_))
    ));
    let displaced = pair.bodies.map(|b| {
        RigidBodyState::new(
            Pose::new(
                b.pose().position_world().add(Vec3::new(0.001, 0.0, 0.0)),
                b.pose().orientation(),
            )
            .unwrap(),
            b.linear_momentum_world(),
            b.angular_momentum_body(),
        )
        .unwrap()
    });
    assert!(matches!(
        attempt(&displaced, contact, 0.0, 300.0),
        Err(RigidHertzHeatError::Invalid(_))
    ));
    assert!(matches!(
        attempt(&pair.bodies, contact, DT, 300.0),
        Err(RigidHertzHeatError::Invalid(_))
    ));
    assert!(matches!(
        attempt(&pair.bodies, contact, 0.0, 330.0),
        Err(RigidHertzHeatError::FrictionSource(_))
    ));
    assert!(matches!(
        attempt(
            &pair.bodies,
            RigidHertzHeatContact {
                nominal_pressure_pa: Some(1e5),
                ..contact
            },
            0.0,
            300.0
        ),
        Err(RigidHertzHeatError::FrictionSource(_))
    ));
    let duplicate = run(
        &pair,
        &owner,
        &thermal,
        &routes[0],
        controls(),
        |t, k| {
            let mut s = source(t, k)?;
            s.transport.conductances_w_k.insert(20, 1.0);
            Ok(s)
        },
        &CancelGate::new(),
    );
    assert!(matches!(duplicate, Err(RigidHertzHeatError::Thermal(_))));
    assert_eq!(pair.bodies, saved);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(energy(&thermal), 0.0);

    // An exactly isothermal native interval establishes a representable nonzero
    // accepted clock. At this scale, target - accepted differs materially from
    // the requested DT, so independent recoil and work expose a split clock.
    let accepted_time_s = 2.0_f64.powi(40);
    let accepted = BackwardEulerStep::new(&thermal, accepted_time_s, || Ok(()))
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
    let target_time_s = accepted_time_s + DT;
    let duration_s = target_time_s - accepted_time_s;
    assert_ne!(duration_s.to_bits(), DT.to_bits());
    assert!((duration_s - DT).abs() > 1e-6);
    let candidate = propose_rigid_hertz_heat(
        &pair.bodies,
        &pair.properties,
        contact,
        &accepted,
        &routes[0],
        accepted_time_s,
        DT,
        controls(),
        &[0.0; 3],
        300.0,
        source,
        &CancelGate::new(),
    )
    .unwrap();
    let (interface, cells, heat, effort) = reference(false, duration_s);
    close(candidate.source, interface);
    close(candidate.mechanical.contact.force_n[0], -effort);
    close(candidate.mechanical.dissipated_energy_j, heat);
    close(
        candidate.mechanical.slip_velocity_after_m_s[0],
        SPEED - duration_s * mobility(false) * effort,
    );
    close(
        candidate.bodies_after[0].linear_momentum_world().x,
        MASS[0] * SPEED - duration_s * effort,
    );
    close(
        candidate.bodies_after[1].linear_momentum_world().x,
        duration_s * effort,
    );
    close(
        candidate.bodies_after[0].angular_momentum_body().y,
        -ARMS[0] * duration_s * effort,
    );
    close(
        candidate.bodies_after[1].angular_momentum_body().y,
        ARMS[1] * duration_s * effort,
    );
    for i in 0..3 {
        close(
            candidate.thermal.temperature_k(CellId(i as u64)).unwrap(),
            cells[i],
        );
    }
    close(
        kinetic(&candidate.bodies_after) - kinetic(&pair.bodies) + energy(&candidate.thermal)
            - energy(&accepted),
        0.0,
    );
    close(candidate.combined_energy_residual_j, 0.0);
    assert_eq!(candidate.target_time_s, target_time_s);
    assert_eq!(candidate.thermal.time_s(), target_time_s);
    assert_eq!(accepted.time_s(), accepted_time_s);
    assert_eq!(energy(&accepted), 0.0);
    assert_eq!(pair.bodies, saved);
}

#[test]
fn final_source_cancellation_and_budget_failures_leave_exact_joint_retry() {
    let pair = Pair::new(SPEED, 0.0);
    let saved = pair.bodies;
    let owner = law(0);
    let original = owner.clone();
    let thermal = thermal([300.0; 3], So3::identity());
    let routes = routes(&thermal);
    let calls = Cell::new(0usize);
    let reference = run(
        &pair,
        &owner,
        &thermal,
        &routes[0],
        controls(),
        |t, k| {
            calls.set(calls.get() + 1);
            source(t, k)
        },
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(calls.get(), reference.trials);
    assert!(reference.trials > 2);
    let calls = Cell::new(0usize);
    let failed = run(
        &pair,
        &owner,
        &thermal,
        &routes[0],
        controls(),
        |t, k| {
            calls.set(calls.get() + 1);
            if calls.get() == reference.trials {
                Err("late bulk source")
            } else {
                source(t, k)
            }
        },
        &CancelGate::new(),
    );
    assert!(matches!(
        failed,
        Err(RigidHertzHeatError::Source("late bulk source"))
    ));
    assert_eq!(calls.get(), reference.trials);
    let calls = Cell::new(0usize);
    let gate = CancelGate::new();
    let cancelled = run(
        &pair,
        &owner,
        &thermal,
        &routes[0],
        controls(),
        |t, k| {
            calls.set(calls.get() + 1);
            if calls.get() == reference.trials {
                gate.request();
            }
            source(t, k)
        },
        &gate,
    );
    assert!(cancelled.is_err());
    assert_eq!(calls.get(), reference.trials);
    assert!(matches!(
        run(
            &pair,
            &owner,
            &thermal,
            &routes[0],
            RigidHertzHeatControls {
                maximum_trials: reference.trials - 1,
                ..controls()
            },
            source,
            &CancelGate::new()
        ),
        Err(RigidHertzHeatError::TrialBudget)
    ));
    assert!(matches!(
        run(
            &pair,
            &owner,
            &thermal,
            &routes[0],
            RigidHertzHeatControls {
                maximum_cells: 2,
                ..controls()
            },
            source,
            &CancelGate::new()
        ),
        Err(RigidHertzHeatError::Invalid(_))
    ));
    let retry = run(
        &pair,
        &owner,
        &thermal,
        &routes[0],
        controls(),
        source,
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(retry.bodies_after, reference.bodies_after);
    assert_eq!(retry.mechanical, reference.mechanical);
    assert_eq!(retry.friction_source, reference.friction_source);
    assert_eq!(retry.source, reference.source);
    for id in thermal.cells().keys() {
        assert_eq!(
            retry.thermal.temperature_k(*id).unwrap(),
            reference.thermal.temperature_k(*id).unwrap()
        );
    }
    audit(&retry, &pair, &thermal);
    assert_eq!(pair.bodies, saved);
    assert_eq!(owner, original);
    assert_eq!(thermal.time_s(), 0.0);
}
