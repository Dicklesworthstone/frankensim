//! Native pressure/piston reciprocity, physical endpoint inertia and refusal.

use std::cell::Cell;

use fs_couple::superfluid::piston::{
    BalancedPistonBoundary, PistonControls, PistonError, advance_balanced_piston,
};
use fs_duct::superfluid::{
    BoundaryVelocities, CellSample, ChannelEnd, FaceSample, LinearizationLimits, MutualFriction,
    TwoFluidChannel,
};
use fs_material::superfluid::{TwoFluidProperties, TwoFluidState};
use fs_mbd::{MassProperties, Pose, RigidBodyState, Vec3};

const RHO: f64 = 8.0;
const AREA: f64 = 0.01;
const C: f64 = 12.0;

fn material(alpha: f64) -> TwoFluidState {
    TwoFluidState::try_new(TwoFluidProperties {
        temperature_k: 2.0,
        density_kg_m3: RHO,
        superfluid_fraction: 0.75,
        entropy_j_kg_k: 4.0,
        specific_heat_cp_j_kg_k: 16.0,
        adiabatic_sound_speed_m_s: C,
        thermal_expansion_per_k: alpha,
        normal_viscosity_pa_s: 0.03,
    })
    .unwrap()
}

fn channel(n: usize, alpha: f64, pressure: f64, incident: bool) -> TwoFluidChannel {
    let mut channel = TwoFluidChannel::try_new(
        material(alpha),
        1.0,
        AREA,
        n,
        LinearizationLimits {
            max_pressure_perturbation_pa: 4.0,
            max_temperature_perturbation_k: 0.1,
            max_component_speed_m_s: 0.2,
        },
    )
    .unwrap();
    let mean = if incident { pressure / (RHO * C) } else { 0.0 };
    let mut normal = vec![mean; n + 1];
    let mut superfluid = normal.clone();
    if alpha != 0.0 {
        for i in 1..n {
            let w = 0.02 * (core::f64::consts::PI * i as f64 / n as f64).sin();
            normal[i] += 0.75 * w;
            superfluid[i] -= 0.25 * w;
        }
    }
    normal[0] = 0.0;
    normal[n] = 0.0;
    superfluid[0] = 0.0;
    superfluid[n] = 0.0;
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

fn body(mass: f64, velocity_y: f64) -> (RigidBodyState, MassProperties) {
    let properties = MassProperties::new(mass, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap();
    let state = RigidBodyState::new(
        Pose::identity(),
        Vec3::new(0.0, mass * velocity_y, 0.0),
        Vec3::ZERO,
    )
    .unwrap();
    (state, properties)
}

fn boundary(end: ChannelEnd) -> BalancedPistonBoundary {
    BalancedPistonBoundary {
        end,
        positive_axis_world: Vec3::new(0.0, 2.0, 0.0),
        equilibrium_position_world: Vec3::ZERO,
        max_displacement_m: 0.05,
    }
}

fn controls() -> PistonControls {
    PistonControls {
        momentum_tolerance_n_s: 1e-12,
        energy_tolerance_j: 1e-12,
        velocity_tolerance_m_s: 1e-11,
    }
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}, tol {tolerance:e}"
    );
}

fn fluid_momentum(channel: &TwoFluidChannel) -> f64 {
    // Independent quadrature of the actual accepted staggered velocity field.
    (1..channel.cells())
        .map(|i| channel.face(i).unwrap().barycentric_velocity_m_s)
        .sum::<f64>()
        * (RHO * AREA * channel.length_m() / channel.cells() as f64)
}

#[test]
fn finite_expansion_piston_transfers_native_momentum_position_and_energy_at_either_end() {
    for end in [ChannelEnd::Left, ChannelEnd::Right] {
        let mut channel = channel(24, 0.1, 0.4, false);
        let (mut body, properties) = body(0.02, 0.0);
        let initial_momentum = fluid_momentum(&channel);
        let initial_energy = channel.wave_energy_j().unwrap();
        let dt = 0.6 * channel.maximum_step_s();
        let sign = if end == ChannelEnd::Left { -1.0 } else { 1.0 };
        let mut support_impulse = 0.0;
        let mut heat = 0.0;
        for _ in 0..12 {
            let report = advance_balanced_piston(
                &mut channel,
                &mut body,
                properties,
                boundary(end),
                dt,
                MutualFriction::Linear { gamma_kg_m3_s: 1.5 },
                controls(),
                || false,
            )
            .unwrap();
            support_impulse += report.stationary_support_impulse_x_n_s;
            heat += report.channel.mutual_friction_heat_j;
            close(report.boundary_half_mass_kg, 0.5 * RHO * AREA / 24.0, 1e-18);
            close(
                report.body_impulse_outward_n_s + report.half_cell_impulse_outward_n_s,
                report.port_impulse_outward_n_s,
                2e-16,
            );
            let velocity = body.center_of_mass_velocity_world(properties).unwrap().y;
            close(
                fluid_momentum(&channel)
                    + body.linear_momentum_world().y
                    + report.boundary_half_mass_kg * velocity
                    + support_impulse,
                initial_momentum,
                3e-15,
            );
            close(
                channel.wave_energy_j().unwrap()
                    + body.kinetic_energy(properties).unwrap()
                    + 0.5 * report.boundary_half_mass_kg * velocity * velocity
                    + heat,
                initial_energy,
                3e-15,
            );
            close(report.combined_energy_residual_j, 0.0, 2e-16);
        }
        assert!(sign * body.linear_momentum_world().y > 0.0);
        assert!(sign * body.pose().position_world().y > 0.0);
        assert!(heat > 0.0);
        assert!(
            (0..channel.cells())
                .any(|i| { channel.cell(i).unwrap().temperature_perturbation_k.abs() > 1e-8 })
        );
    }
}

#[test]
fn heavy_piston_approaches_the_original_closed_end() {
    let mut fixed = channel(24, 0.1, 0.4, false);
    let mut moving = channel(24, 0.1, 0.4, false);
    let (mut body, properties) = body(1e9, 0.0);
    let dt = 0.5 * fixed.maximum_step_s();
    fixed
        .step(dt, BoundaryVelocities::CLOSED, MutualFriction::None, || {
            false
        })
        .unwrap();
    let report = advance_balanced_piston(
        &mut moving,
        &mut body,
        properties,
        boundary(ChannelEnd::Right),
        dt,
        MutualFriction::None,
        controls(),
        || false,
    )
    .unwrap();
    assert!(report.midpoint_velocity_outward_m_s.abs() < 1e-14);
    for i in 0..fixed.cells() {
        close(
            moving.cell(i).unwrap().pressure_perturbation_pa,
            fixed.cell(i).unwrap().pressure_perturbation_pa,
            2e-12,
        );
    }
}

#[test]
fn incident_pressure_step_refines_to_continuum_piston_motion() {
    // Before the opposite-end reflection can arrive, incoming characteristic
    // p + rho*c*v = 2*p_inc gives M*dU/dt=A*(2*p_inc-rho*c*U).
    let pressure = 0.2;
    let mass = 0.02;
    let target = 0.02; // < length/c: the fixed left end cannot yet affect the piston.
    let exact =
        2.0 * pressure / (RHO * C) * (1.0 - fs_math::det::exp(-RHO * C * AREA * target / mass));
    let mut errors = Vec::new();
    for n in [32, 64, 128] {
        let mut channel = channel(n, 0.0, pressure, true);
        let (mut body, properties) = body(mass, 0.0);
        let steps = (target / (0.6 * channel.maximum_step_s())).ceil() as usize;
        let dt = target / steps as f64;
        for _ in 0..steps {
            advance_balanced_piston(
                &mut channel,
                &mut body,
                properties,
                boundary(ChannelEnd::Right),
                dt,
                MutualFriction::None,
                controls(),
                || false,
            )
            .unwrap();
        }
        errors.push((body.center_of_mass_velocity_world(properties).unwrap().y - exact).abs());
    }
    assert!(errors[1] < 0.7 * errors[0], "refinement errors: {errors:?}");
    assert!(errors[2] < 0.7 * errors[1], "refinement errors: {errors:?}");
    assert!(
        errors[2] < 0.01 * exact,
        "refinement errors: {errors:?}, exact {exact:e}"
    );
}

#[derive(Debug, PartialEq)]
struct Snapshot {
    time: f64,
    energy: f64,
    heat: f64,
    cells: Vec<CellSample>,
    faces: Vec<FaceSample>,
    last_heat: Vec<f64>,
    body: RigidBodyState,
}
fn snapshot(channel: &TwoFluidChannel, body: RigidBodyState) -> Snapshot {
    Snapshot {
        time: channel.time_s(),
        energy: channel.wave_energy_j().unwrap(),
        heat: channel.cumulative_mutual_friction_heat_j(),
        cells: (0..channel.cells())
            .map(|i| channel.cell(i).unwrap())
            .collect(),
        faces: (0..=channel.cells())
            .map(|i| channel.face(i).unwrap())
            .collect(),
        last_heat: channel.last_mutual_friction_heat_j().to_vec(),
        body,
    }
}

#[test]
fn late_geometry_and_native_mechanics_refusals_preserve_both_owners() {
    let mut channel = channel(16, 0.1, 0.4, false);
    let (mut body, properties) = body(0.02, 0.0);
    let dt = 0.6 * channel.maximum_step_s();
    let before = snapshot(&channel, body);
    let mut tiny = boundary(ChannelEnd::Right);
    tiny.max_displacement_m = 1e-12;
    assert!(matches!(
        advance_balanced_piston(
            &mut channel,
            &mut body,
            properties,
            tiny,
            dt,
            MutualFriction::None,
            controls(),
            || false,
        ),
        Err(PistonError::Invalid(_))
    ));
    assert_eq!(snapshot(&channel, body), before);

    // The original native owner refuses an unrepresentable p^2 diagnostic.
    // The coupled path must preserve that refusal, rather than replace its step.
    let (mut enormous, enormous_properties) = self::body(1e200, 1e-5);
    let before = snapshot(&channel, enormous);
    assert!(matches!(
        advance_balanced_piston(
            &mut channel,
            &mut enormous,
            enormous_properties,
            boundary(ChannelEnd::Right),
            dt,
            MutualFriction::None,
            controls(),
            || false,
        ),
        Err(PistonError::Mechanical(_))
    ));
    assert_eq!(snapshot(&channel, enormous), before);

    // Changing port type must not discard momentum/energy already present in
    // the now-fixed end or in a counterflowing selected half-cell.
    for moving_opposite_end in [true, false] {
        let mut incompatible = self::channel(16, 0.1, 0.4, false);
        let mut normal = vec![0.0; 17];
        let mut superfluid = normal.clone();
        if moving_opposite_end {
            normal[0] = 0.001;
            superfluid[0] = 0.001;
        } else {
            normal[16] = 0.001;
            superfluid[16] = -0.001;
        }
        incompatible
            .initialize(&[0.4; 16], &[0.0; 16], &normal, &superfluid, || false)
            .unwrap();
        let (mut body, properties) = self::body(0.02, 0.0);
        let before = snapshot(&incompatible, body);
        assert!(matches!(
            advance_balanced_piston(
                &mut incompatible,
                &mut body,
                properties,
                boundary(ChannelEnd::Right),
                dt,
                MutualFriction::None,
                controls(),
                || false,
            ),
            Err(PistonError::Invalid(_))
        ));
        assert_eq!(snapshot(&incompatible, body), before);
    }
}

#[test]
fn cancellation_at_entry_and_final_acceptance_does_not_publish_and_retry_matches() {
    let mut baseline = channel(16, 0.1, 0.4, false);
    let (mut baseline_body, properties) = body(0.02, 0.0);
    let dt = 0.6 * baseline.maximum_step_s();
    let calls = Cell::new(0usize);
    let baseline_report = advance_balanced_piston(
        &mut baseline,
        &mut baseline_body,
        properties,
        boundary(ChannelEnd::Right),
        dt,
        MutualFriction::GorterMellink { a_m_s_kg: 2.0 },
        controls(),
        || {
            calls.set(calls.get() + 1);
            false
        },
    )
    .unwrap();
    let expected = snapshot(&baseline, baseline_body);
    for stop in [1, calls.get()] {
        let mut channel = channel(16, 0.1, 0.4, false);
        let (mut body, properties) = body(0.02, 0.0);
        let before = snapshot(&channel, body);
        let count = Cell::new(0usize);
        let result = advance_balanced_piston(
            &mut channel,
            &mut body,
            properties,
            boundary(ChannelEnd::Right),
            dt,
            MutualFriction::GorterMellink { a_m_s_kg: 2.0 },
            controls(),
            || {
                count.set(count.get() + 1);
                count.get() >= stop
            },
        );
        assert!(
            matches!(result, Err(PistonError::Cancelled)),
            "stop {stop}: {result:?}"
        );
        assert_eq!(snapshot(&channel, body), before);
        let retried = advance_balanced_piston(
            &mut channel,
            &mut body,
            properties,
            boundary(ChannelEnd::Right),
            dt,
            MutualFriction::GorterMellink { a_m_s_kg: 2.0 },
            controls(),
            || false,
        )
        .unwrap();
        assert_eq!(retried, baseline_report);
        assert_eq!(snapshot(&channel, body), expected);
    }
}
