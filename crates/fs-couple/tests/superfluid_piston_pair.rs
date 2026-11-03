//! G0/G1 momentum, work and co-acceleration; G4 joint-owner rollback.

use std::cell::Cell;

use fs_couple::superfluid::piston::{BalancedPistonBoundary, PistonControls, PistonError};
use fs_couple::superfluid::piston_pair::advance_balanced_piston_pair;
use fs_duct::superfluid::{
    BodyForce, CellSample, ChannelEnd, ComponentForceDensity, FaceSample, LinearizationLimits,
    MutualFriction, TwoFluidChannel,
};
use fs_material::superfluid::{TwoFluidProperties, TwoFluidState};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

const RHO: f64 = 8.0;
const AREA: f64 = 0.01;

fn channel(
    n: usize,
    alpha: f64,
    pressure: f64,
    velocity: f64,
    relative: f64,
    pressure_limit: f64,
) -> TwoFluidChannel {
    let material = TwoFluidState::try_new(TwoFluidProperties {
        temperature_k: 2.0,
        density_kg_m3: RHO,
        superfluid_fraction: 0.75,
        entropy_j_kg_k: 4.0,
        specific_heat_cp_j_kg_k: 16.0,
        adiabatic_sound_speed_m_s: 12.0,
        thermal_expansion_per_k: alpha,
        normal_viscosity_pa_s: 0.03,
    })
    .unwrap();
    let mut channel = TwoFluidChannel::try_new(
        material,
        1.0,
        AREA,
        n,
        LinearizationLimits {
            max_pressure_perturbation_pa: pressure_limit,
            max_temperature_perturbation_k: 0.1,
            max_component_speed_m_s: 0.2,
        },
    )
    .unwrap();
    let mut normal = vec![velocity; n + 1];
    let mut superfluid = normal.clone();
    for i in 1..n {
        let w = relative * (core::f64::consts::PI * i as f64 / n as f64).sin();
        normal[i] += 0.75 * w;
        superfluid[i] -= 0.25 * w;
    }
    channel
        .initialize(
            &vec![pressure; n],
            &vec![0.0; n],
            &normal,
            &superfluid,
            || false,
        )
        .unwrap();
    channel
}

fn mechanics(
    masses: [f64; 2],
    velocity: f64,
) -> (
    [RigidBodyState; 2],
    [MassProperties; 2],
    [BalancedPistonBoundary; 2],
) {
    let properties =
        masses.map(|m| MassProperties::new(m, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap());
    let positions = [Vec3::ZERO, Vec3::new(0.0, 1.0, 0.0)];
    let bodies = [0, 1].map(|i| {
        RigidBodyState::new(
            Pose::new(positions[i], UnitQuaternion::IDENTITY).unwrap(),
            Vec3::new(0.0, masses[i] * velocity, 0.0),
            Vec3::ZERO,
        )
        .unwrap()
    });
    let boundaries = [ChannelEnd::Left, ChannelEnd::Right].map(|end| BalancedPistonBoundary {
        end,
        positive_axis_world: Vec3::new(0.0, 2.0, 0.0),
        equilibrium_position_world: positions[usize::from(end == ChannelEnd::Right)],
        max_displacement_m: 0.05,
    });
    (bodies, properties, boundaries)
}

fn controls() -> PistonControls {
    PistonControls {
        momentum_tolerance_n_s: 1e-12,
        energy_tolerance_j: 1e-12,
        velocity_tolerance_m_s: 1e-11,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}, tol={tolerance:e}"
    );
}

fn total_momentum(
    channel: &TwoFluidChannel,
    bodies: [RigidBodyState; 2],
    properties: [MassProperties; 2],
) -> f64 {
    let interior = (1..channel.cells())
        .map(|i| channel.face(i).unwrap().barycentric_velocity_m_s)
        .sum::<f64>();
    let volume = AREA / channel.cells() as f64;
    interior * RHO * volume
        + (0..2)
            .map(|i| {
                bodies[i].linear_momentum_world().y
                    + 0.5
                        * RHO
                        * volume
                        * bodies[i]
                            .center_of_mass_velocity_world(properties[i])
                            .unwrap()
                            .y
            })
            .sum::<f64>()
}

fn total_energy(
    channel: &TwoFluidChannel,
    bodies: [RigidBodyState; 2],
    properties: [MassProperties; 2],
) -> f64 {
    let half_mass = 0.5 * RHO * AREA / channel.cells() as f64;
    channel.wave_energy_j().unwrap()
        + (0..2)
            .map(|i| {
                let v = bodies[i]
                    .center_of_mass_velocity_world(properties[i])
                    .unwrap()
                    .y;
                bodies[i].kinetic_energy(properties[i]).unwrap() + 0.5 * half_mass * v * v
            })
            .sum::<f64>()
}

#[test]
fn freely_translating_channel_and_unequal_bodies_keep_momentum_without_false_pressure() {
    // Tight pressure admission is intentional: an artificial blocked-end
    // intermediate would produce O(rho*c*v) pressure and refuse this solution.
    for (n, alpha) in [(2, 0.0), (17, 0.1)] {
        let v = 0.02;
        let mut fluid = channel(n, alpha, 0.0, v, 0.0, 1e-10);
        let (mut bodies, properties, boundaries) = mechanics([0.02, 0.035], v);
        let p0 = total_momentum(&fluid, bodies, properties);
        let e0 = total_energy(&fluid, bodies, properties);
        let dt = 0.75 * fluid.maximum_step_s();
        for _ in 0..20 {
            let report = advance_balanced_piston_pair(
                &mut fluid,
                &mut bodies,
                properties,
                boundaries,
                dt,
                BodyForce::None,
                [0.0; 2],
                MutualFriction::None,
                controls(),
                || false,
            )
            .unwrap();
            close(report.external_impulse_world_n_s.y, 0.0, 0.0);
            close(report.total_external_work_j, 0.0, 0.0);
            close(total_momentum(&fluid, bodies, properties), p0, 2e-16);
            close(total_energy(&fluid, bodies, properties), e0, 2e-17);
            for i in 0..2 {
                close(
                    bodies[i]
                        .center_of_mass_velocity_world(properties[i])
                        .unwrap()
                        .y,
                    v,
                    2e-16,
                );
                close(
                    bodies[i].pose().position_world().y
                        - boundaries[i].equilibrium_position_world.y,
                    v * fluid.time_s(),
                    2e-14,
                );
            }
            for i in 0..n {
                close(fluid.cell(i).unwrap().pressure_perturbation_pa, 0.0, 1e-12);
                close(
                    fluid.cell(i).unwrap().temperature_perturbation_k,
                    0.0,
                    1e-14,
                );
            }
        }
    }
}

#[test]
fn common_acceleration_moves_all_fluid_mass_and_both_native_bodies_without_spurious_waves() {
    for n in [2, 32] {
        let acceleration = 0.35;
        let v0 = 0.01;
        let masses = [0.02, 0.035];
        let mut fluid = channel(n, 0.1, 0.0, v0, 0.0, 1e-10);
        let (mut bodies, properties, boundaries) = mechanics(masses, v0);
        let p0 = total_momentum(&fluid, bodies, properties);
        let e0 = total_energy(&fluid, bodies, properties);
        let total_mass = RHO * AREA + masses[0] + masses[1];
        let dt = 0.6 * fluid.maximum_step_s();
        let mut work = 0.0;
        for _ in 0..16 {
            let report = advance_balanced_piston_pair(
                &mut fluid,
                &mut bodies,
                properties,
                boundaries,
                dt,
                BodyForce::Acceleration(acceleration),
                [-masses[0] * acceleration, masses[1] * acceleration],
                MutualFriction::None,
                controls(),
                || false,
            )
            .unwrap();
            work += report.total_external_work_j;
            let expected_velocity = v0 + acceleration * fluid.time_s();
            close(
                report.external_impulse_world_n_s.y,
                total_mass * acceleration * dt,
                3e-18,
            );
            close(
                total_momentum(&fluid, bodies, properties),
                p0 + total_mass * acceleration * fluid.time_s(),
                3e-16,
            );
            close(total_energy(&fluid, bodies, properties) - e0, work, 3e-17);
            for i in 0..2 {
                close(
                    bodies[i]
                        .center_of_mass_velocity_world(properties[i])
                        .unwrap()
                        .y,
                    expected_velocity,
                    5e-16,
                );
                close(report.ends[i].fluid_contact_impulse_outward_n_s, 0.0, 2e-16);
            }
            for i in 0..n {
                close(fluid.cell(i).unwrap().pressure_perturbation_pa, 0.0, 3e-12);
                close(
                    fluid.cell(i).unwrap().temperature_perturbation_k,
                    0.0,
                    1e-13,
                );
            }
            for i in 1..n {
                close(
                    fluid.face(i).unwrap().normal_velocity_m_s,
                    expected_velocity,
                    1e-15,
                );
                close(
                    fluid.face(i).unwrap().superfluid_velocity_m_s,
                    expected_velocity,
                    1e-15,
                );
            }
        }
    }
}

#[test]
fn pressure_recoil_and_finite_expansion_friction_preserve_isolated_momentum_and_energy() {
    // Two cells retain a substantial cross-end response; 24 cells exercise the
    // same native wave owner after pressure and second-sound propagation.
    for n in [2, 24] {
        let v0 = 0.007;
        let mut fluid = channel(n, 0.1, 0.4, v0, 0.02, 4.0);
        let (mut bodies, properties, boundaries) = mechanics([0.02, 0.035], v0);
        let p0 = total_momentum(&fluid, bodies, properties);
        let e0 = total_energy(&fluid, bodies, properties);
        let dt = 0.6 * fluid.maximum_step_s();
        let mut heat = 0.0;
        for _ in 0..16 {
            let report = advance_balanced_piston_pair(
                &mut fluid,
                &mut bodies,
                properties,
                boundaries,
                dt,
                BodyForce::None,
                [0.0; 2],
                MutualFriction::GorterMellink { a_m_s_kg: 2.0 },
                controls(),
                || false,
            )
            .unwrap();
            heat += report.channel.mutual_friction_heat_j;
            close(total_momentum(&fluid, bodies, properties), p0, 3e-16);
            close(total_energy(&fluid, bodies, properties) + heat, e0, 3e-17);
            close(report.combined_momentum_residual_world_n_s.y, 0.0, 2e-16);
            close(report.combined_energy_residual_j, 0.0, 2e-17);
            for end in report.ends {
                close(
                    end.body_impulse_outward_n_s + end.half_cell_impulse_outward_n_s,
                    end.port_impulse_outward_n_s,
                    2e-16,
                );
            }
        }
        assert!(
            bodies[0]
                .center_of_mass_velocity_world(properties[0])
                .unwrap()
                .y
                < v0
        );
        assert!(
            bodies[1]
                .center_of_mass_velocity_world(properties[1])
                .unwrap()
                .y
                > v0
        );
        assert!(heat > 0.0);
        assert!((0..n).any(|i| fluid.cell(i).unwrap().temperature_perturbation_k.abs() > 1e-8));
    }
}

#[test]
fn spatial_component_and_piston_forces_keep_full_volume_impulse_and_actual_work() {
    let n = 8;
    let mut fluid = channel(n, 0.1, 0.1, 0.004, 0.01, 4.0);
    let (mut bodies, properties, boundaries) = mechanics([1e-9, 1e3], 0.004);
    let field: Vec<_> = (0..=n)
        .map(|i| ComponentForceDensity {
            normal_n_m3: 0.2 * (1.0 + i as f64 / n as f64),
            superfluid_n_m3: -0.07 * (0.5 + i as f64 / n as f64),
        })
        .collect();
    let external = [0.006, -0.004];
    let fluid_force = AREA / n as f64
        * field
            .iter()
            .enumerate()
            .map(|(i, f)| {
                (f.normal_n_m3 + f.superfluid_n_m3) * if i == 0 || i == n { 0.5 } else { 1.0 }
            })
            .sum::<f64>();
    let total_force = fluid_force - external[0] + external[1];
    let p0 = total_momentum(&fluid, bodies, properties);
    let e0 = total_energy(&fluid, bodies, properties);
    let dt = 0.5 * fluid.maximum_step_s();
    let mut heat = 0.0;
    let mut work = 0.0;
    for _ in 0..20 {
        let report = advance_balanced_piston_pair(
            &mut fluid,
            &mut bodies,
            properties,
            boundaries,
            dt,
            BodyForce::Faces(&field),
            external,
            MutualFriction::Linear { gamma_kg_m3_s: 1.5 },
            controls(),
            || false,
        )
        .unwrap();
        heat += report.channel.mutual_friction_heat_j;
        work += report.total_external_work_j;
        close(report.external_impulse_world_n_s.y, total_force * dt, 3e-18);
        close(
            total_momentum(&fluid, bodies, properties),
            p0 + total_force * fluid.time_s(),
            2e-14,
        );
        close(
            total_energy(&fluid, bodies, properties) - e0 + heat,
            work,
            3e-16,
        );
        for i in 0..2 {
            let sign = if i == 0 { -1.0 } else { 1.0 };
            let f = field[i * n];
            close(
                report.ends[i].half_cell_source_impulse_outward_n_s,
                sign * 0.5 * AREA / n as f64 * (f.normal_n_m3 + f.superfluid_n_m3) * dt,
                1e-20,
            );
            close(
                report.ends[i].fluid_contact_impulse_outward_n_s,
                report.ends[i].body_impulse_outward_n_s - external[i] * dt,
                1e-20,
            );
        }
    }
    assert!(heat > 0.0);
}

#[derive(Debug, PartialEq)]
struct Snapshot {
    time: f64,
    energy: f64,
    heat: f64,
    cells: Vec<CellSample>,
    faces: Vec<FaceSample>,
    last_heat: Vec<f64>,
    bodies: [RigidBodyState; 2],
}

fn snapshot(fluid: &TwoFluidChannel, bodies: [RigidBodyState; 2]) -> Snapshot {
    Snapshot {
        time: fluid.time_s(),
        energy: fluid.wave_energy_j().unwrap(),
        heat: fluid.cumulative_mutual_friction_heat_j(),
        cells: (0..fluid.cells()).map(|i| fluid.cell(i).unwrap()).collect(),
        faces: (0..=fluid.cells())
            .map(|i| fluid.face(i).unwrap())
            .collect(),
        last_heat: fluid.last_mutual_friction_heat_j().to_vec(),
        bodies,
    }
}

#[test]
fn second_body_refusal_and_late_cancellation_preserve_all_owners_and_retry() {
    let make = || {
        (
            channel(8, 0.1, 0.4, 0.0, 0.02, 4.0),
            mechanics([0.02, 0.035], 0.0),
        )
    };
    let (mut baseline, (mut baseline_bodies, properties, boundaries)) = make();
    let dt = 0.6 * baseline.maximum_step_s();
    let calls = Cell::new(0usize);
    let expected_report = advance_balanced_piston_pair(
        &mut baseline,
        &mut baseline_bodies,
        properties,
        boundaries,
        dt,
        BodyForce::Acceleration(0.1),
        [0.001, -0.002],
        MutualFriction::Linear { gamma_kg_m3_s: 1.5 },
        controls(),
        || {
            calls.set(calls.get() + 1);
            false
        },
    )
    .unwrap();
    let expected = snapshot(&baseline, baseline_bodies);
    for stop in [1, calls.get()] {
        let (mut fluid, (mut bodies, _, _)) = make();
        let before = snapshot(&fluid, bodies);
        let count = Cell::new(0usize);
        let result = advance_balanced_piston_pair(
            &mut fluid,
            &mut bodies,
            properties,
            boundaries,
            dt,
            BodyForce::Acceleration(0.1),
            [0.001, -0.002],
            MutualFriction::Linear { gamma_kg_m3_s: 1.5 },
            controls(),
            || {
                count.set(count.get() + 1);
                count.get() >= stop
            },
        );
        assert!(
            matches!(result, Err(PistonError::Cancelled)),
            "stop={stop}, result={result:?}"
        );
        assert_eq!(snapshot(&fluid, bodies), before);
        let retried = advance_balanced_piston_pair(
            &mut fluid,
            &mut bodies,
            properties,
            boundaries,
            dt,
            BodyForce::Acceleration(0.1),
            [0.001, -0.002],
            MutualFriction::Linear { gamma_kg_m3_s: 1.5 },
            controls(),
            || false,
        )
        .unwrap();
        assert_eq!(retried, expected_report);
        assert_eq!(snapshot(&fluid, bodies), expected);
    }
    let (mut fluid, (mut bodies, _, mut limited)) = make();
    let before = snapshot(&fluid, bodies);
    limited[1].max_displacement_m = 1e-12;
    assert!(matches!(
        advance_balanced_piston_pair(
            &mut fluid,
            &mut bodies,
            properties,
            limited,
            dt,
            BodyForce::None,
            [0.0; 2],
            MutualFriction::None,
            controls(),
            || false
        ),
        Err(PistonError::Invalid(_))
    ));
    assert_eq!(snapshot(&fluid, bodies), before);
}
