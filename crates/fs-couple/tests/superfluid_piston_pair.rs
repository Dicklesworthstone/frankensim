//! G0/G1 momentum, work and co-acceleration; G4 joint-owner rollback.

use std::cell::Cell;

use fs_couple::superfluid::piston::{BalancedPistonBoundary, PistonControls, PistonError};
use fs_couple::superfluid::piston_pair::advance_balanced_piston_pair;
use fs_duct::superfluid::{
    BodyForce, CellSample, ChannelEnd, ComponentForceDensity, FaceSample, LinearizationLimits,
    MutualFriction, TwoFluidChannel,
};
use fs_material::superfluid::{
    NormalViscousSpectrum, NormalViscousWall, TwoFluidProperties, TwoFluidState,
};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

const RHO: f64 = 8.0;
const AREA: f64 = 0.01;

#[derive(Clone, Copy)]
enum NormalWallModel {
    None,
    Darcy(NormalViscousWall),
    CircularProfile,
    SlitProfile,
}

impl NormalWallModel {
    fn is_profile(self) -> bool {
        matches!(self, Self::CircularProfile | Self::SlitProfile)
    }
}

fn channel(
    n: usize,
    alpha: f64,
    pressure: f64,
    velocity: f64,
    relative: f64,
    pressure_limit: f64,
) -> TwoFluidChannel {
    channel_with_wall(n, alpha, pressure, velocity, relative, pressure_limit, None)
}

fn channel_with_wall(
    n: usize,
    alpha: f64,
    pressure: f64,
    velocity: f64,
    relative: f64,
    pressure_limit: f64,
    wall: Option<NormalViscousWall>,
) -> TwoFluidChannel {
    channel_with_model(
        n,
        alpha,
        pressure,
        velocity,
        relative,
        pressure_limit,
        wall.map_or(NormalWallModel::None, NormalWallModel::Darcy),
    )
}

fn channel_with_model(
    n: usize,
    alpha: f64,
    pressure: f64,
    velocity: f64,
    relative: f64,
    pressure_limit: f64,
    model: NormalWallModel,
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
    let limits = LinearizationLimits {
        max_pressure_perturbation_pa: pressure_limit,
        max_temperature_perturbation_k: 0.1,
        max_component_speed_m_s: 0.2,
    };
    let mut channel = match model {
        NormalWallModel::Darcy(wall) => {
            TwoFluidChannel::try_new_with_normal_wall(material, 1.0, AREA, n, limits, wall)
        }
        NormalWallModel::None => TwoFluidChannel::try_new(material, 1.0, AREA, n, limits),
        NormalWallModel::CircularProfile | NormalWallModel::SlitProfile => {
            // Both spectra have the scalar case's DC permeability, but their
            // independent profile memory and kinetic storage remain active.
            let spectrum = match model {
                NormalWallModel::CircularProfile => NormalViscousSpectrum::circular_capillary(
                    material,
                    (8.0_f64 * 0.0003).sqrt(),
                    8,
                ),
                _ => {
                    NormalViscousSpectrum::parallel_plates(material, (12.0_f64 * 0.0003).sqrt(), 8)
                }
            }
            .unwrap();
            TwoFluidChannel::try_new_with_normal_profile(material, 1.0, AREA, n, limits, spectrum)
        }
    }
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
    if model.is_profile() {
        channel.initialize_normal_profile_steady(|| false).unwrap();
        assert!(channel.normal_profile_excess_energy_j().unwrap() > 0.0);
    }
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
    let fluid_energy = if let Some(spectrum) = channel.normal_viscous_spectrum() {
        // Independently integrate physical thermodynamic storage and the FULL
        // weighted normal profile, instead of trusting a mean-only face norm.
        let volume = AREA / channel.cells() as f64;
        let material = channel.material();
        let mut energy = 0.0;
        for i in 0..channel.cells() {
            let cell = channel.cell(i).unwrap();
            let p = cell.pressure_perturbation_pa;
            let theta = cell.temperature_perturbation_k;
            energy += 0.5
                * volume
                * (material.isothermal_compressibility_per_pa() * p * p
                    - 2.0 * material.thermal_expansion_per_k() * p * theta
                    + material.density_kg_m3() * material.specific_heat_cp_j_kg_k()
                        / material.temperature_k()
                        * theta
                        * theta);
        }
        for i in 1..channel.cells() {
            let normal_squared: f64 = spectrum
                .modes()
                .iter()
                .enumerate()
                .map(|(mode, coefficients)| {
                    let velocity = channel.normal_profile_mode_velocity_m_s(i, mode).unwrap();
                    coefficients.weight * velocity * velocity
                })
                .sum();
            let superfluid = channel.face(i).unwrap().superfluid_velocity_m_s;
            energy += 0.5
                * volume
                * (material.normal_density_kg_m3() * normal_squared
                    + material.superfluid_density_kg_m3() * superfluid * superfluid);
        }
        close(
            channel.wave_energy_j().unwrap(),
            energy,
            5e-13 * energy.abs() + 1e-21,
        );
        energy
    } else {
        channel.wave_energy_j().unwrap()
    };
    fluid_energy
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

#[test]
fn stationary_wall_reaction_closes_native_pair_momentum_and_separate_heat() {
    let wall = NormalViscousWall::try_new(0.0003).unwrap();
    for (n, model) in [2, 24].into_iter().flat_map(|n| {
        [
            NormalWallModel::Darcy(wall),
            NormalWallModel::CircularProfile,
            NormalWallModel::SlitProfile,
        ]
        .map(|model| (n, model))
    }) {
        let v0 = 0.007;
        let mut fluid = channel_with_model(n, 0.1, 0.1, v0, 0.02, 4.0, model);
        let profile_before = fluid.normal_profile_excess_energy_j().unwrap();
        let (mut bodies, properties, boundaries) = mechanics([0.02, 0.035], v0);
        let p0 = total_momentum(&fluid, bodies, properties);
        let e0 = total_energy(&fluid, bodies, properties);
        let dt = 0.6 * fluid.maximum_step_s();
        let mut wall_impulse = 0.0;
        let mut mutual_heat = 0.0;
        let mut viscous_heat = 0.0;
        for _ in 0..12 {
            let report = advance_balanced_piston_pair(
                &mut fluid,
                &mut bodies,
                properties,
                boundaries,
                dt,
                BodyForce::None,
                [0.0; 2],
                MutualFriction::Linear { gamma_kg_m3_s: 1.5 },
                controls(),
                || false,
            )
            .unwrap();
            wall_impulse += report.channel.interior_wall_impulse_n_s;
            mutual_heat += report.channel.mutual_friction_heat_j;
            viscous_heat += report.channel.normal_viscous_heat_j;
            // The stationary support supplies momentum but no mechanical work.
            close(
                report.external_impulse_world_n_s.y,
                report.channel.interior_wall_impulse_n_s,
                1e-20,
            );
            close(report.total_external_work_j, 0.0, 0.0);
            close(
                total_momentum(&fluid, bodies, properties),
                p0 + wall_impulse,
                3e-16,
            );
            close(
                total_energy(&fluid, bodies, properties) + mutual_heat + viscous_heat,
                e0,
                3e-16,
            );
            close(report.combined_momentum_residual_world_n_s.y, 0.0, 2e-16);
            close(report.combined_energy_residual_j, 0.0, 2e-16);
            close(
                fluid.cumulative_mutual_friction_heat_j(),
                mutual_heat,
                1e-18,
            );
            close(
                fluid.cumulative_normal_viscous_heat_j(),
                viscous_heat,
                1e-18,
            );
            for end in report.ends {
                close(
                    end.body_impulse_outward_n_s + end.half_cell_impulse_outward_n_s,
                    end.port_impulse_outward_n_s,
                    2e-16,
                );
            }
        }
        assert!(wall_impulse < 0.0);
        assert!(mutual_heat > 0.0 && viscous_heat > 0.0);
        if model.is_profile() {
            let profile_after = fluid.normal_profile_excess_energy_j().unwrap();
            assert!((profile_after - profile_before).abs() > 1e-6 * profile_before);
        }
    }
}

#[derive(Debug, PartialEq)]
struct Snapshot {
    time: f64,
    energy: f64,
    heat: f64,
    viscous_heat: f64,
    cells: Vec<CellSample>,
    faces: Vec<FaceSample>,
    last_heat: Vec<f64>,
    last_viscous_heat: Vec<f64>,
    profile_modes: Vec<f64>,
    bodies: [RigidBodyState; 2],
}

fn snapshot(fluid: &TwoFluidChannel, bodies: [RigidBodyState; 2]) -> Snapshot {
    let profile_modes = fluid
        .normal_viscous_spectrum()
        .map_or_else(Vec::new, |spectrum| {
            (1..fluid.cells())
                .flat_map(|face| {
                    (0..spectrum.modes().len()).map(move |mode| {
                        fluid.normal_profile_mode_velocity_m_s(face, mode).unwrap()
                    })
                })
                .collect()
        });
    Snapshot {
        time: fluid.time_s(),
        energy: fluid.wave_energy_j().unwrap(),
        heat: fluid.cumulative_mutual_friction_heat_j(),
        viscous_heat: fluid.cumulative_normal_viscous_heat_j(),
        cells: (0..fluid.cells()).map(|i| fluid.cell(i).unwrap()).collect(),
        faces: (0..=fluid.cells())
            .map(|i| fluid.face(i).unwrap())
            .collect(),
        last_heat: fluid.last_mutual_friction_heat_j().to_vec(),
        last_viscous_heat: fluid.last_normal_viscous_heat_j().to_vec(),
        profile_modes,
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

mod thermal_tests {
    use super::*;
    use fs_couple::superfluid::piston_pair::{
        PistonPairHeatError, advance_balanced_piston_pair_with_heat,
    };
    use fs_couple::superfluid::{FrictionHeatControls, FrictionHeatError};
    use fs_thermal::enthalpy::{
        Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
        DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, ThermalStorageBasis,
        TimedEvent,
    };

    const N: usize = 8;
    const CELL_MASS: f64 = RHO * AREA / N as f64;

    struct Caloric {
        maximum_specific_energy_j_kg: f64,
    }

    impl ThermalLaw for Caloric {
        fn storage_basis(&self) -> ThermalStorageBasis {
            ThermalStorageBasis::InternalEnergy
        }
        fn temperature_k(&self, energy: f64) -> Result<f64, DepositionError> {
            if !energy.is_finite() || !(0.0..=self.maximum_specific_energy_j_kg).contains(&energy) {
                return Err(DepositionError::MaterialDomain);
            }
            // Authored calorics differ from the acoustic cp: the composition
            // must use the original receiving owner, not copy a wave-law cp.
            Ok(2.0 + energy / 8.0)
        }
        fn maximum_temperature_slope(&self) -> f64 {
            0.125
        }
    }

    fn heat_controls() -> FrictionHeatControls {
        FrictionHeatControls {
            thermal: StepConfig {
                max_dt_s: 0.001,
                energy_tolerance_j: 1e-13,
            },
            cell_heat_tolerance_j: 1e-15,
            heat_energy_tolerance_j: 1e-14,
            combined_energy_tolerance_j: 1e-12,
        }
    }

    fn destinations() -> Vec<CellId> {
        (0..N).map(|i| CellId(i as u64 + 10)).collect()
    }

    fn thermal(law: &dyn ThermalLaw, source_w: f64) -> DepositionSimulation<'_> {
        let schedule = destinations()
            .iter()
            .enumerate()
            .map(|(i, &id)| TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id,
                    segment_id: i as u64,
                    body_id: 1,
                    material_id: MaterialId(1),
                    mass_kg: CELL_MASS,
                    specific_h_j_kg: 0.0,
                    geometry: CellGeometry {
                        start_m: [i as f64 / N as f64, 0.0, 0.0],
                        end_m: [(i + 1) as f64 / N as f64, 0.0, 0.0],
                        volume_m3: AREA / N as f64,
                        surface_area_m2: AREA,
                    },
                    boundary: Boundary {
                        ambient_k: 2.0,
                        htc_w_m2_k: 0.0,
                        reservoir: None,
                        source_w,
                    },
                }),
            })
            .collect();
        let original = DepositionSimulation::new(
            schedule,
            &[(MaterialId(1), law)],
            DepositionBudget {
                max_cells: N,
                max_links: 0,
                max_events: N,
                max_substeps: 100,
            },
        )
        .unwrap();
        original
            .advance_to(0.0, heat_controls().thermal, || Ok(()))
            .unwrap()
            .0
    }

    fn thermal_snapshot(thermal: &DepositionSimulation<'_>) -> (f64, Vec<DepositedCell>) {
        (
            thermal.time_s(),
            thermal.cells().values().cloned().collect(),
        )
    }

    #[test]
    fn driven_pair_deposits_each_friction_increment_once_and_closes_full_storage() {
        check_driven_pair_heat(NormalWallModel::None);
    }

    #[test]
    fn viscous_driven_pair_routes_both_losses_and_external_wall_impulse() {
        check_driven_pair_heat(NormalWallModel::Darcy(
            NormalViscousWall::try_new(0.0003).unwrap(),
        ));
    }

    #[test]
    fn unsteady_profiles_preserve_native_pair_storage_and_both_caloric_heat_routes() {
        for model in [
            NormalWallModel::CircularProfile,
            NormalWallModel::SlitProfile,
        ] {
            check_driven_pair_heat(model);
        }
    }

    fn check_driven_pair_heat(model: NormalWallModel) {
        let law = Caloric {
            maximum_specific_energy_j_kg: 2.0,
        };
        let source = 1e-5;
        let mut thermal = thermal(&law, source);
        let mut fluid = channel_with_model(N, 0.1, 0.4, 0.003, 0.02, 4.0, model);
        let original_material = *fluid.material();
        let (mut bodies, properties, boundaries) = mechanics([0.02, 0.035], 0.003);
        let ids = destinations();
        let force = BodyForce::Uniform(ComponentForceDensity {
            normal_n_m3: 0.08,
            superfluid_n_m3: -0.08,
        });
        let external = [0.001, -0.002];
        let dt = 0.4 * fluid.maximum_step_s();
        let e0 = total_energy(&fluid, bodies, properties) + thermal.energy_j();
        let p0 = total_momentum(&fluid, bodies, properties);
        let mut source_energy = 0.0;
        let mut total_heat = 0.0;
        let mut wall_impulse = 0.0;
        for _ in 0..12 {
            let old_energy: Vec<_> = ids
                .iter()
                .map(|id| thermal.cells()[id].specific_h_j_kg)
                .collect();
            let report = advance_balanced_piston_pair_with_heat(
                &mut fluid,
                &mut bodies,
                properties,
                boundaries,
                &mut thermal,
                &ids,
                dt,
                force,
                external,
                MutualFriction::Linear { gamma_kg_m3_s: 1.5 },
                controls(),
                heat_controls(),
                || false,
            )
            .unwrap();
            total_heat += report.pair.channel.dissipated_heat_j();
            wall_impulse += report.pair.channel.interior_wall_impulse_n_s;
            source_energy +=
                report.pair.total_external_work_j + report.heat.thermal_march.external_heat_j;
            assert_eq!(thermal.time_s(), fluid.time_s());
            assert_eq!(fluid.material(), &original_material);
            assert!(thermal.uses_material_law(MaterialId(1), &law));
            close(
                report.heat.thermal_march.external_heat_j,
                N as f64 * source * dt,
                1e-18,
            );
            close(
                report.heat.heat_intake.external_heat_j,
                report.pair.channel.dissipated_heat_j(),
                1e-18,
            );
            close(report.combined_storage_residual_j, 0.0, 1e-16);
            for (i, id) in ids.iter().enumerate() {
                let energy = thermal.cells()[id].specific_h_j_kg;
                close(
                    CELL_MASS * (energy - old_energy[i]),
                    source * dt
                        + fluid.last_mutual_friction_heat_j()[i]
                        + fluid.last_normal_viscous_heat_j()[i],
                    1e-18,
                );
                close(
                    thermal.temperature_k(*id).unwrap(),
                    2.0 + energy / 8.0,
                    1e-15,
                );
                assert_eq!(thermal.cells()[id].boundary.source_w, source);
            }
            close(
                total_energy(&fluid, bodies, properties) + thermal.energy_j() - e0,
                source_energy,
                5e-16,
            );
            close(
                total_momentum(&fluid, bodies, properties),
                p0 + (-external[0] + external[1]) * fluid.time_s() + wall_impulse,
                5e-16,
            );
            close(
                thermal.energy_j(),
                N as f64 * source * fluid.time_s() + total_heat,
                1e-17,
            );
        }
        assert!(total_heat > 0.0);
        if !matches!(model, NormalWallModel::None) {
            assert!(fluid.cumulative_normal_viscous_heat_j() > 0.0);
            assert!(wall_impulse != 0.0);
        }
        assert!(
            bodies
                .iter()
                .zip(boundaries)
                .any(|(body, boundary)| body.pose().position_world().y
                    != boundary.equilibrium_position_world.y)
        );
    }

    #[test]
    fn caloric_and_final_cancellation_refusals_publish_none_of_the_four_owners() {
        check_pair_heat_refusal(NormalWallModel::None);
    }

    #[test]
    fn viscous_caloric_and_final_cancellation_refusals_preserve_all_owners() {
        check_pair_heat_refusal(NormalWallModel::Darcy(
            NormalViscousWall::try_new(0.0003).unwrap(),
        ));
    }

    #[test]
    fn unsteady_profile_caloric_refusal_and_final_cancellation_preserve_all_owners() {
        check_pair_heat_refusal(NormalWallModel::CircularProfile);
    }

    fn check_pair_heat_refusal(model: NormalWallModel) {
        let make = || {
            (
                channel_with_model(N, 0.1, 0.4, 0.003, 0.02, 4.0, model),
                mechanics([0.02, 0.035], 0.003),
            )
        };
        let good = Caloric {
            maximum_specific_energy_j_kg: 2.0,
        };
        let refused = Caloric {
            maximum_specific_energy_j_kg: 0.0,
        };
        let ids = destinations();
        let (mut fluid, (mut bodies, properties, boundaries)) = make();
        let dt = 0.4 * fluid.maximum_step_s();
        let before = snapshot(&fluid, bodies);
        let mut cold = thermal(&refused, 0.0);
        let cold_before = thermal_snapshot(&cold);
        let result = advance_balanced_piston_pair_with_heat(
            &mut fluid,
            &mut bodies,
            properties,
            boundaries,
            &mut cold,
            &ids,
            dt,
            BodyForce::Acceleration(0.1),
            [0.001, -0.002],
            MutualFriction::Linear { gamma_kg_m3_s: 1.5 },
            controls(),
            heat_controls(),
            || false,
        );
        assert!(matches!(
            result,
            Err(PistonPairHeatError::Heat(FrictionHeatError::Thermal(
                DepositionError::MaterialDomain
            )))
        ));
        assert_eq!(snapshot(&fluid, bodies), before);
        assert_eq!(thermal_snapshot(&cold), cold_before);

        let (mut baseline, (mut baseline_bodies, _, _)) = make();
        let mut baseline_thermal = thermal(&good, 0.0);
        let calls = Cell::new(0usize);
        let expected_report = advance_balanced_piston_pair_with_heat(
            &mut baseline,
            &mut baseline_bodies,
            properties,
            boundaries,
            &mut baseline_thermal,
            &ids,
            dt,
            BodyForce::Acceleration(0.1),
            [0.001, -0.002],
            MutualFriction::Linear { gamma_kg_m3_s: 1.5 },
            controls(),
            heat_controls(),
            || {
                calls.set(calls.get() + 1);
                false
            },
        )
        .unwrap();
        let expected = snapshot(&baseline, baseline_bodies);
        let expected_thermal = thermal_snapshot(&baseline_thermal);
        for stop in [1, calls.get()] {
            let (mut fluid, (mut bodies, _, _)) = make();
            let mut thermal = thermal(&good, 0.0);
            let before = snapshot(&fluid, bodies);
            let thermal_before = thermal_snapshot(&thermal);
            let count = Cell::new(0usize);
            let result = advance_balanced_piston_pair_with_heat(
                &mut fluid,
                &mut bodies,
                properties,
                boundaries,
                &mut thermal,
                &ids,
                dt,
                BodyForce::Acceleration(0.1),
                [0.001, -0.002],
                MutualFriction::Linear { gamma_kg_m3_s: 1.5 },
                controls(),
                heat_controls(),
                || {
                    count.set(count.get() + 1);
                    count.get() >= stop
                },
            );
            assert!(
                matches!(result, Err(PistonPairHeatError::Cancelled)),
                "stop={stop}, result={result:?}"
            );
            assert_eq!(snapshot(&fluid, bodies), before);
            assert_eq!(thermal_snapshot(&thermal), thermal_before);
            let retried = advance_balanced_piston_pair_with_heat(
                &mut fluid,
                &mut bodies,
                properties,
                boundaries,
                &mut thermal,
                &ids,
                dt,
                BodyForce::Acceleration(0.1),
                [0.001, -0.002],
                MutualFriction::Linear { gamma_kg_m3_s: 1.5 },
                controls(),
                heat_controls(),
                || false,
            )
            .unwrap();
            assert_eq!(retried, expected_report);
            assert_eq!(snapshot(&fluid, bodies), expected);
            assert_eq!(thermal_snapshot(&thermal), expected_thermal);
        }
    }
}
