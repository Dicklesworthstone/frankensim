use super::*;
use crate::cavitation::{BubbleGas, BubbleLimits};
use core::f64::consts::PI;

const R: f64 = 50e-6;
const D: f64 = 1e-3;
const AMBIENT: f64 = 4000.0;

fn model(viscosity: f64) -> RayleighPlessetBubble {
    let liquid = CavitationState::try_new(300.0, 1000.0, viscosity, 0.072, 2000.0, 1500.0).unwrap();
    RayleighPlessetBubble::try_new(
        liquid,
        BubbleGas::equilibrium(liquid, R, AMBIENT, 1.4).unwrap(),
        BubbleLimits {
            minimum_radius_m: 0.5 * R,
            maximum_radius_m: 2.0 * R,
            maximum_wall_mach_number: 0.01,
            maximum_vapor_density_ratio: 0.01,
        },
    )
    .unwrap()
}
fn omega0() -> f64 {
    let bubble = model(0.0);
    ((3.0 * 1.4 * bubble.gas().reference_pressure_pa() - 2.0 * 0.072 / R) / (1000.0 * R * R)).sqrt()
}
fn limits() -> BubbleClusterLimits {
    BubbleClusterLimits {
        maximum_bubbles: 16,
        maximum_radius_to_separation: 0.1,
        maximum_inertia_row_sum: 0.5,
        maximum_liquid_pressure_ratio: 0.01,
        maximum_source_angular_frequency_rad_s: 1.5 * omega0(),
        maximum_retardation_phase_rad: 0.1,
    }
}
fn numerics(duration: f64) -> BubbleClusterNumerics {
    BubbleClusterNumerics {
        radial: BubbleNumerics {
            relative_tolerance: 1e-11,
            radius_absolute_tolerance_m: 1e-18,
            velocity_absolute_tolerance_m_s: 1e-12,
            initial_step_s: duration / 1000.0,
            minimum_step_s: duration * 1e-12,
            maximum_step_s: duration / 100.0,
            maximum_attempts: 20_000,
        },
        absolute_energy_tolerance_j: 1e-18,
    }
}
fn state(radius: f64, velocity: f64) -> BubbleState {
    BubbleState::try_new(0.0, radius, velocity).unwrap()
}
fn pair(states: [BubbleState; 2], viscosity: f64) -> BubbleCluster {
    let bubble = model(viscosity);
    BubbleCluster::try_new(
        vec![
            ClusterBubble {
                bubble,
                center_m: [0.0, 0.0, 0.0],
            },
            ClusterBubble {
                bubble,
                center_m: [D, 0.0, 0.0],
            },
        ],
        states.to_vec(),
        0.0,
        AMBIENT,
        limits(),
        AMBIENT,
    )
    .unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}, tolerance={tolerance:.4e}"
    );
}

#[test]
fn single_bubble_uses_the_exact_original_radial_trajectory_and_source() {
    let bubble = model(0.0);
    let initial = state(0.96 * R, -0.02);
    let mut isolated = initial;
    let mut cluster = BubbleCluster::try_new(
        vec![ClusterBubble {
            bubble,
            center_m: [3.0, 4.0, 5.0],
        }],
        vec![initial],
        0.0,
        AMBIENT,
        limits(),
        AMBIENT,
    )
    .unwrap();
    let duration = 2e-5;
    let controls = numerics(duration);
    let expected = bubble
        .advance(&mut isolated, duration, &controls.radial, &|_| AMBIENT)
        .unwrap();
    let actual = cluster.advance(duration, &controls, &|_| AMBIENT).unwrap();
    assert_eq!(cluster.states(), &[isolated]);
    assert_eq!(actual.observation.bubbles[0], expected.observation);
    assert_eq!(actual.accepted_steps, expected.accepted_steps);
    assert_eq!(actual.rejected_steps, expected.rejected_steps);
    assert_eq!(actual.rhs_evaluations, expected.rhs_evaluations);
    assert_eq!(actual.pressure_work_j, 0.0);
    assert_eq!(actual.viscous_loss_j, 0.0);
    assert_eq!(cluster.members()[0].bubble, bubble);
}

#[test]
fn pair_accelerations_match_the_independent_two_by_two_radial_equations() {
    let states = [state(0.94 * R, 0.08), state(1.04 * R, -0.05)];
    let cluster = pair(states, 0.002);
    let pressure = AMBIENT + 700.0;
    let actual = cluster.observe(pressure).unwrap();
    let isolated: Vec<_> = cluster
        .members()
        .iter()
        .zip(states)
        .map(|(member, s)| member.bubble.observe(&s, pressure).unwrap())
        .collect();
    let r1 = states[0].radius_m();
    let r2 = states[1].radius_m();
    let u1 = states[0].wall_velocity_m_s();
    let u2 = states[1].wall_velocity_m_s();
    let a = r1;
    let b = r2 * r2 / D;
    let c = r1 * r1 / D;
    let d = r2;
    let f1 = r1 * isolated[0].wall_acceleration_m_s2 - 2.0 * r2 * u2 * u2 / D;
    let f2 = r2 * isolated[1].wall_acceleration_m_s2 - 2.0 * r1 * u1 * u1 / D;
    let determinant = a * d - b * c;
    let expected = [
        (d * f1 - b * f2) / determinant,
        (a * f2 - c * f1) / determinant,
    ];
    for i in 0..2 {
        near(
            actual.bubbles[i].wall_acceleration_m_s2,
            expected[i],
            2e-14 * expected[i].abs(),
        );
        assert_ne!(
            actual.bubbles[i].wall_acceleration_m_s2,
            isolated[i].wall_acceleration_m_s2
        );
        let j = 1 - i;
        let neighbor_pressure = 1000.0
            * (states[j].radius_m().powi(2) * expected[j]
                + 2.0 * states[j].radius_m() * states[j].wall_velocity_m_s().powi(2))
            / D;
        near(
            actual.bubbles[i].imposed_absolute_pressure_pa,
            pressure + neighbor_pressure,
            2e-11,
        );
    }
    let kinetic = 2.0 * PI * 1000.0 * (r1.powi(3) * u1 * u1 + r2.powi(3) * u2 * u2)
        + 4.0 * PI * 1000.0 * r1 * r1 * u1 * r2 * r2 * u2 / D;
    near(actual.liquid_kinetic_energy_j, kinetic, 3e-15 * kinetic);
    assert!(
        actual.liquid_kinetic_energy_j
            < actual
                .bubbles
                .iter()
                .map(|o| o.liquid_kinetic_energy_j)
                .sum::<f64>()
    );
}

#[test]
fn collective_in_phase_and_antisymmetric_frequencies_match_independent_modes() {
    let amplitude = R * 1e-5;
    for sign in [1.0, -1.0] {
        let omega = omega0() / (1.0 + sign * R / D).sqrt();
        let phase = 4.3 * PI;
        let duration = phase / omega;
        let mut cluster = pair(
            [state(R + amplitude, 0.0), state(R + sign * amplitude, 0.0)],
            0.0,
        );
        let mut controls = numerics(duration);
        controls.radial.relative_tolerance = 1e-12;
        let step = cluster.advance(duration, &controls, &|_| AMBIENT).unwrap();
        for (i, direction) in [1.0, sign].into_iter().enumerate() {
            near(
                cluster.states()[i].radius_m() - R,
                direction * amplitude * phase.cos(),
                amplitude * 3e-4,
            );
            near(
                cluster.states()[i].wall_velocity_m_s(),
                -direction * amplitude * omega * phase.sin(),
                amplitude * omega * 3e-4,
            );
        }
        assert_eq!(step.pressure_work_j, 0.0);
        assert_eq!(step.viscous_loss_j, 0.0);
        assert!(step.energy_balance_residual_j.abs() < 1e-18);
    }
}

// Independent energy formula in radius coordinates, including the pair cross
// kinetic energy and the analytic kappa!=1 gas potential.
fn physical_energy(states: &[BubbleState; 2], bubble: RayleighPlessetBubble) -> f64 {
    let mut total = 0.0;
    let rho = bubble.liquid().liquid_density_kg_m3();
    let sigma = bubble.liquid().surface_tension_n_m();
    let vapor = bubble.liquid().saturation_pressure_pa();
    let gas = bubble.gas();
    let kappa = gas.polytropic_exponent();
    let v0 = 4.0 * PI / 3.0 * gas.reference_radius_m().powi(3);
    for state in states {
        let r = state.radius_m();
        let u = state.wall_velocity_m_s();
        let volume = 4.0 * PI / 3.0 * r.powi(3);
        total += 2.0 * PI * rho * r.powi(3) * u * u
            + 4.0 * PI * sigma * r * r
            + (AMBIENT - vapor) * volume
            + gas.reference_pressure_pa() * v0 * ((v0 / volume).powf(kappa - 1.0) - 1.0)
                / (kappa - 1.0);
    }
    total
        + 4.0
            * PI
            * rho
            * states[0].radius_m().powi(2)
            * states[0].wall_velocity_m_s()
            * states[1].radius_m().powi(2)
            * states[1].wall_velocity_m_s()
            / D
}

#[test]
fn coupled_acceleration_closes_the_independent_lagrange_energy_power_identity() {
    let states = [state(0.92 * R, 0.08), state(1.06 * R, -0.05)];
    let bubble = model(0.002);
    let cluster = pair(states, 0.002);
    let pressure = AMBIENT + 700.0;
    let observation = cluster.observe(pressure).unwrap();
    near(
        observation.total_mechanical_energy_j,
        physical_energy(&states, bubble),
        5e-23,
    );
    let h = 1e-10;
    let displaced = |direction: f64| -> [BubbleState; 2] {
        core::array::from_fn(|i| {
            state(
                states[i].radius_m() + direction * h * states[i].wall_velocity_m_s(),
                states[i].wall_velocity_m_s()
                    + direction * h * observation.bubbles[i].wall_acceleration_m_s2,
            )
        })
    };
    let energy_rate = (physical_energy(&displaced(1.0), bubble)
        - physical_energy(&displaced(-1.0), bubble))
        / (2.0 * h);
    let mut expected = 0.0;
    for s in states {
        let r = s.radius_m();
        let u = s.wall_velocity_m_s();
        expected += -(pressure - AMBIENT) * 4.0 * PI * r * r * u - 16.0 * PI * 0.002 * r * u * u;
    }
    near(energy_rate, expected, 2e-6 * expected.abs());
    near(
        observation.drive_power_w - observation.viscous_dissipation_w,
        expected,
        2e-20,
    );
}

#[test]
fn nonlinear_joint_trajectory_accounts_for_pressure_work_and_viscous_loss() {
    let initial = [state(0.96 * R, 0.01), state(1.02 * R, -0.03)];
    let mut cluster = pair(initial, 0.002);
    let before = physical_energy(&initial, model(0.002));
    let duration = 2e-5;
    let drive = |time: f64| AMBIENT + 200.0 * (0.2 * omega0() * time).sin();
    let step = cluster
        .advance(duration, &numerics(duration), &drive)
        .unwrap();
    let final_states = [cluster.states()[0], cluster.states()[1]];
    let after = physical_energy(&final_states, model(0.002));
    assert!(step.viscous_loss_j > 0.0);
    assert_ne!(step.pressure_work_j, 0.0);
    near(
        after - before,
        step.pressure_work_j - step.viscous_loss_j,
        1e-18,
    );
    assert!(cluster.states().iter().all(|s| s.time_s() == duration));
    assert_eq!(cluster.time_s(), duration);
}

#[test]
fn exchanging_member_order_preserves_the_same_physical_pair() {
    let initial = [state(0.97 * R, 0.04), state(1.02 * R, -0.03)];
    let mut forward = pair(initial, 0.001);
    let mut members = forward.members().to_vec();
    members.reverse();
    let mut reverse = BubbleCluster::try_new(
        members,
        vec![initial[1], initial[0]],
        0.0,
        AMBIENT,
        limits(),
        AMBIENT,
    )
    .unwrap();
    let duration = 1e-5;
    forward
        .advance(duration, &numerics(duration), &|_| AMBIENT)
        .unwrap();
    reverse
        .advance(duration, &numerics(duration), &|_| AMBIENT)
        .unwrap();
    for i in 0..2 {
        near(
            forward.states()[i].radius_m(),
            reverse.states()[1 - i].radius_m(),
            R * 5e-12,
        );
        near(
            forward.states()[i].wall_velocity_m_s(),
            reverse.states()[1 - i].wall_velocity_m_s(),
            5e-9,
        );
    }
}

#[test]
fn cancellation_budget_and_late_member_domain_failure_never_publish_partial_state() {
    let initial = [state(0.97 * R, 0.01), state(1.02 * R, -0.02)];
    let duration = 1e-5;
    let controls = numerics(duration);
    let mut reference = pair(initial, 0.001);
    let mut polls = 0;
    let expected = reference
        .advance_controlled(duration, &controls, &|_| AMBIENT, &mut || {
            polls += 1;
            false
        })
        .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut cluster = pair(initial, 0.001);
        let originals = cluster.members().to_vec();
        let mut count = 0;
        let error = cluster
            .advance_controlled(duration, &controls, &|_| AMBIENT, &mut || {
                count += 1;
                count == stop
            })
            .unwrap_err();
        assert_eq!(error, BubbleClusterError::Cancelled);
        assert_eq!(cluster.states(), &initial);
        assert_eq!(cluster.members(), originals);
        assert_eq!(cluster.time_s(), 0.0);
        let retry = cluster.advance(duration, &controls, &|_| AMBIENT).unwrap();
        assert_eq!(retry, expected);
        assert_eq!(cluster.states(), reference.states());
    }
    let mut cluster = pair(initial, 0.001);
    let mut limited = controls;
    limited.radial.maximum_attempts = 1;
    assert!(matches!(
        cluster.advance(duration, &limited, &|_| AMBIENT),
        Err(BubbleClusterError::Native(
            CavitationError::StepBudgetExceeded { attempts: 1, .. }
        ))
    ));
    assert_eq!(cluster.states(), &initial);

    let bubble = model(0.001);
    let tight = RayleighPlessetBubble::try_new(
        bubble.liquid(),
        bubble.gas(),
        BubbleLimits {
            minimum_radius_m: R * (1.0 - 1e-6),
            ..bubble.limits()
        },
    )
    .unwrap();
    let inward = [state(R, -0.1), state(R, -0.1)];
    let mut cluster = BubbleCluster::try_new(
        vec![
            ClusterBubble {
                bubble,
                center_m: [0.0, 0.0, 0.0],
            },
            ClusterBubble {
                bubble: tight,
                center_m: [D, 0.0, 0.0],
            },
        ],
        inward.to_vec(),
        0.0,
        AMBIENT,
        limits(),
        AMBIENT,
    )
    .unwrap();
    let error = cluster
        .advance(duration, &controls, &|_| AMBIENT)
        .unwrap_err();
    assert!(
        matches!(
            error,
            BubbleClusterError::Member {
                index: 1,
                cause: CavitationError::RadiusOutsideDomain { .. }
            }
        ),
        "{error:?}"
    );
    assert_eq!(cluster.states(), &inward);
    assert_eq!(cluster.time_s(), 0.0);
}

#[test]
fn separation_collective_inertia_bandwidth_material_and_pressure_domains_refuse() {
    let initial = vec![state(R, 0.0), state(R, 0.0)];
    let bubble = model(0.0);
    let members = vec![
        ClusterBubble {
            bubble,
            center_m: [0.0, 0.0, 0.0],
        },
        ClusterBubble {
            bubble,
            center_m: [D, 0.0, 0.0],
        },
    ];
    let mut close = members.clone();
    close[1].center_m[0] = 3.0 * R;
    assert!(matches!(
        BubbleCluster::try_new(close, initial.clone(), 0.0, AMBIENT, limits(), AMBIENT,),
        Err(BubbleClusterError::Geometry { .. })
    ));
    let mut weak = limits();
    weak.maximum_inertia_row_sum = 0.05;
    assert!(matches!(
        BubbleCluster::try_new(
            members.clone(),
            initial.clone(),
            0.0,
            AMBIENT,
            weak,
            AMBIENT,
        ),
        Err(BubbleClusterError::OutsideDomain(_))
    ));
    let mut fast = limits();
    fast.maximum_source_angular_frequency_rad_s *= 10.0;
    assert!(matches!(
        BubbleCluster::try_new(
            members.clone(),
            initial.clone(),
            0.0,
            AMBIENT,
            fast,
            AMBIENT,
        ),
        Err(BubbleClusterError::OutsideDomain(_))
    ));
    let other = CavitationState::try_new(300.0, 999.0, 0.0, 0.072, 2000.0, 1500.0).unwrap();
    let mut changed = members.clone();
    changed[1].bubble =
        RayleighPlessetBubble::try_new(other, bubble.gas(), bubble.limits()).unwrap();
    assert!(matches!(
        BubbleCluster::try_new(changed, initial.clone(), 0.0, AMBIENT, limits(), AMBIENT,),
        Err(BubbleClusterError::MaterialMismatch { index: 1 })
    ));
    assert!(matches!(
        BubbleCluster::try_new(
            members.clone(),
            initial.clone(),
            0.0,
            AMBIENT,
            limits(),
            AMBIENT + 3e7,
        ),
        Err(BubbleClusterError::Member {
            cause: CavitationError::PressureRatioExceeded { .. },
            ..
        })
    ));
    let mut budget = limits();
    budget.maximum_bubbles = 1;
    assert!(BubbleCluster::try_new(members, initial, 0.0, AMBIENT, budget, AMBIENT).is_err());
}
