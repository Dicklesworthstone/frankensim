use super::*;
use fs_aeroac::curle3d::refraction::full_wave::{
    FullWaveForceLimits, array::FullWaveCompactSourceLimits,
};
use fs_duct::network::{
    JunctionCorrections, NetworkEdge, NetworkNode, NodeBoundary, PortExcitation,
};
use fs_duct::{Duct, Segment};
use fs_material::gas::GasSpec;
use std::f64::consts::{PI, TAU};

const HZ: f64 = 700.0;
fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}
fn limits() -> CompactWallRadiationLimits {
    CompactWallRadiationLimits {
        maximum_outlets: 4,
        maximum_total_evaluations: 262_144,
        maximum_ka: 0.3,
        maximum_radius_over_distance: 0.1,
        monopole: FullWaveMonopoleLimits {
            absolute_pressure_tolerance_pa_per_m3_s: 1e-7,
            absolute_velocity_tolerance_per_m2: 1e-10,
            relative_tolerance: 1e-10,
            maximum_evaluations: 65_536,
            maximum_depth: 36,
        },
    }
}
fn wall(z: Option<C64>) -> PlanarImpedanceWall {
    let state = gas().acoustic_state().unwrap();
    PlanarImpedanceWall::new(
        state,
        [0.0; 3],
        [0.0, 0.0, 1.0],
        z.map(|v| v.scale(state.characteristic_impedance_pa_s_m())),
    )
    .unwrap()
}
fn outlets() -> [CompactWallOutlet; 2] {
    [
        CompactWallOutlet {
            label: "left".into(),
            position_m: [0.0, 0.0, -0.03],
            radius_m: 0.002,
        },
        CompactWallOutlet {
            label: "right".into(),
            position_m: [0.12, -0.02, -0.045],
            radius_m: 0.0015,
        },
    ]
}
fn close(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * b.abs().max(1e-12),
        "{a:?} != {b:?}"
    );
}
fn green(distance: f64) -> C64 {
    let gas = gas();
    let phase = TAU * HZ * distance / gas.sound_speed;
    let (sin, cos) = phase.sin_cos();
    C64::new(sin, -cos).scale(TAU * HZ * gas.density / (4.0 * PI * distance))
}
fn self_load(radius: f64) -> C64 {
    let gas = gas();
    let ka = TAU * HZ * radius / gas.sound_speed;
    C64::new(0.25 * ka * ka, -0.6133 * ka)
        .scale(gas.density * gas.sound_speed / (PI * radius * radius))
}

#[test]
fn ideal_wall_load_keeps_original_self_mutual_and_signed_reflected_impedance() {
    let ports = outlets();
    for (z, sign) in [(None, 1.0), (Some(C64::ZERO), -1.0)] {
        let load = PreparedCompactWallRadiation::new(
            wall(z),
            HZ,
            &ports,
            "original two-mouth fixture",
            limits(),
            &mut || false,
        )
        .unwrap();
        for i in 0..2 {
            for j in 0..2 {
                let image = [
                    ports[j].position_m[0],
                    ports[j].position_m[1],
                    -ports[j].position_m[2],
                ];
                let returned = green(distance(ports[i].position_m, image)).scale(sign);
                let direct = if i == j {
                    self_load(ports[i].radius_m)
                } else {
                    green(distance(ports[i].position_m, ports[j].position_m))
                };
                close(load.reflected_impedance()[i * 2 + j], returned, 3e-12);
                close(load.free_space_termination().z[i * 2 + j], direct, 3e-12);
                close(load.termination().z[i * 2 + j], direct + returned, 3e-12);
            }
        }
        assert_eq!(load.termination().z[1], load.termination().z[2]);
        assert_eq!(load.evaluations(), 0);
        assert_eq!(load.reflected_diagnostics().len(), 3);
        load.termination().validate().unwrap();
        if sign < 0.0 {
            assert!(load.reflected_impedance()[0].re < 0.0);
            assert!(load.termination().z[0].re > 0.0);
        }
    }
    // A passive rigid wall can also have NEGATIVE added resistance. The free
    // self term is essential: R_total/R_self = 1 + sin(2kh)/(2kh).
    let x = 1.5 * PI;
    let mut port = ports[0].clone();
    port.position_m[2] = -x * gas().sound_speed / (2.0 * TAU * HZ);
    let load = PreparedCompactWallRadiation::new(
        wall(None),
        HZ,
        &[port],
        "signed return",
        limits(),
        &mut || false,
    )
    .unwrap();
    assert!(load.reflected_impedance()[0].re < 0.0);
    close(
        C64::from_re(load.termination().z[0].re / self_load(ports[0].radius_m).re),
        C64::from_re(1.0 + x.sin() / x),
        3e-12,
    );
}

fn port(label: &str) -> NetworkNode {
    NetworkNode {
        id: label.into(),
        boundary: NodeBoundary::Port {
            label: label.into(),
        },
    }
}
fn edge(id: &str, from: usize, to: usize, radius: f64, length: f64) -> NetworkEdge {
    NetworkEdge {
        id: id.into(),
        from,
        to,
        duct: Duct {
            segments: vec![Segment::Cylinder { radius, length }],
        },
        corrections: JunctionCorrections::Disabled,
    }
}
fn assert_work(response: &NetworkSourceResponse) {
    let f = &response.fields;
    let scale = response.generator_power_w.abs().max(1e-20);
    assert!(
        (response.generator_power_w - response.source_dissipation_w - f.source_power_w).abs()
            < 2e-10 * scale
    );
    assert!((f.source_power_w - f.external_load_power_w - f.network_power_w).abs() < 2e-10 * scale);
    assert!(f.network_power_w.abs() < 2e-10 * scale);
    assert!(f.external_load_power_w > 0.0);
}

#[test]
fn pressure_driven_native_duct_reacts_to_wall_and_observes_actual_outward_flow() {
    let gas = gas();
    let outlet = outlets()[0].clone();
    let radius = outlet.radius_m;
    let length = 0.08;
    let graph = DuctNetwork {
        nodes: vec![port("driver"), port(&outlet.label)],
        edges: vec![edge("original bore", 0, 1, radius, length)],
    };
    let zc = gas.density * gas.sound_speed / (PI * radius * radius);
    let (s, c) = (TAU * HZ * length / gas.sound_speed).sin_cos();
    let free = self_load(radius);
    let zin =
        (free.scale(c) - C64::new(0.0, zc * s)) / (C64::from_re(c) - C64::new(0.0, s / zc) * free);
    // Original passive pressure generator with compliance cancelling the
    // free-space input reactance at this declared frequency, not a forced Q.
    let zs = C64::new(zin.re, -zin.im);
    let drive = C64::new(3.0, -1.2);
    let sources = [PortSource {
        label: "driver".into(),
        excitation: PortExcitation::Pressure {
            open_circuit_pressure_pa: drive,
            impedance_pa_s_m3: zs,
        },
    }];
    let load = PreparedCompactWallRadiation::new(
        wall(Some(C64::new(0.7, 0.9))),
        HZ,
        &[outlet.clone()],
        "original bore mouth",
        limits(),
        &mut || false,
    )
    .unwrap();
    let response = load
        .drive_network(&graph, &gas, LossModel::Lossless, &sources, &[])
        .unwrap();
    let z = load.termination().z[0];
    let pressure_per_flow = z.scale(c) - C64::new(0.0, zc * s);
    let inlet_per_outlet = C64::from_re(c) - C64::new(0.0, s / zc) * z;
    let q = drive / (pressure_per_flow + zs * inlet_per_outlet);
    close(response.fields.port_inflows[1], -q, 2e-10);
    close(response.fields.port_inflows[0], inlet_per_outlet * q, 2e-10);
    close(response.fields.port_pressures[1], z * q, 2e-10);
    assert_work(&response);
    let unreflected = network_drive_sources(
        &graph,
        &gas,
        TAU * HZ,
        LossModel::Lossless,
        &sources,
        &[],
        Some(load.free_space_termination()),
    )
    .unwrap();
    assert!(
        (response.fields.port_inflows[1] - unreflected.fields.port_inflows[1]).abs()
            > 0.01 * unreflected.fields.port_inflows[1].abs()
    );
    let receiver = [0.15, 0.04, -0.08];
    let observer = PreparedImpedanceWallCompactSourceObserver::new_total(
        load.wall(),
        HZ,
        receiver,
        &[outlet.position_m],
        FullWaveCompactSourceLimits {
            max_sources: 1,
            maximum_total_evaluations: 262_144,
            monopole: limits().monopole,
            force: FullWaveForceLimits {
                absolute_pressure_tolerance_pa_per_n: 1e-8,
                absolute_velocity_tolerance_m_s_per_n: 1e-11,
                relative_tolerance: 1e-9,
                maximum_evaluations: 65_536,
                maximum_depth: 36,
            },
        },
        &mut || false,
    )
    .unwrap();
    let actual = load.observe(&response.fields, &observer).unwrap();
    let unit = PreparedImpedanceWallMonopole::new(
        load.wall(),
        HZ,
        outlet.position_m,
        receiver,
        limits().monopole,
        &mut || false,
    )
    .unwrap()
    .evaluate_total(q)
    .unwrap();
    close(actual.pressure_pa, unit.pressure_pa, 2e-10);
    for i in 0..3 {
        close(
            actual.particle_velocity_m_s[i],
            unit.particle_velocity_m_s[i],
            2e-10,
        );
    }
}

#[test]
fn mutually_loaded_native_outlets_retain_coupled_flow_and_complete_power() {
    let gas = gas();
    let ports = outlets();
    let lengths = [0.071, 0.103];
    let graph = DuctNetwork {
        nodes: vec![port("driver"), port("left"), port("right")],
        edges: vec![
            edge("left bore", 0, 1, ports[0].radius_m, lengths[0]),
            edge("right bore", 0, 2, ports[1].radius_m, lengths[1]),
        ],
    };
    let zs = C64::new(10_000.0, 200_000.0);
    let drive = C64::new(7.0, 2.0);
    let sources = [PortSource {
        label: "driver".into(),
        excitation: PortExcitation::Pressure {
            open_circuit_pressure_pa: drive,
            impedance_pa_s_m3: zs,
        },
    }];
    let load = PreparedCompactWallRadiation::new(
        wall(None),
        HZ,
        &ports,
        "original bifurcating bore",
        limits(),
        &mut || false,
    )
    .unwrap();
    let response = load
        .drive_network(&graph, &gas, LossModel::Lossless, &sources, &[])
        .unwrap();
    // Independent branch equations: A Q = p_driver [1,1], followed by the
    // pressure generator's own impedance and the sum of inlet branch flows.
    let z = &load.termination().z;
    let sin = lengths.map(|l| (TAU * HZ * l / gas.sound_speed).sin());
    let cos = lengths.map(|l| (TAU * HZ * l / gas.sound_speed).cos());
    let zc = ports
        .each_ref()
        .map(|p| gas.density * gas.sound_speed / (PI * p.radius_m * p.radius_m));
    let a = [
        z[0].scale(cos[0]) - C64::new(0.0, zc[0] * sin[0]),
        z[1].scale(cos[0]),
        z[2].scale(cos[1]),
        z[3].scale(cos[1]) - C64::new(0.0, zc[1] * sin[1]),
    ];
    let determinant = a[0] * a[3] - a[1] * a[2];
    let v = [(a[3] - a[1]) / determinant, (a[0] - a[2]) / determinant];
    let inlet_per_pressure = (0..2).fold(C64::ZERO, |sum, i| {
        sum + v[i].scale(cos[i])
            - C64::new(0.0, sin[i] / zc[i]) * (z[i * 2] * v[0] + z[i * 2 + 1] * v[1])
    });
    let p = drive / (C64::ONE + zs * inlet_per_pressure);
    let q = v.map(|value| value * p);
    for i in 0..2 {
        close(response.fields.port_inflows[i + 1], -q[i], 2e-10);
    }
    close(response.fields.port_pressures[0], p, 2e-10);
    let load_power = (0..2)
        .map(|i| 0.5 * ((z[i * 2] * q[0] + z[i * 2 + 1] * q[1]) * q[i].conj()).re)
        .sum::<f64>();
    close(
        C64::from_re(response.fields.external_load_power_w),
        C64::from_re(load_power),
        2e-10,
    );
    assert_work(&response);
}

#[test]
fn finite_wall_loading_reuses_one_budget_and_refuses_unsupported_geometry() {
    let ports = outlets();
    let boundary = wall(Some(C64::new(0.4, 0.8)));
    let owner = PreparedCompactWallRadiation::new(
        boundary,
        HZ,
        &ports,
        "finite wall",
        limits(),
        &mut || false,
    )
    .unwrap();
    let work = owner
        .reflected_diagnostics()
        .iter()
        .map(|d| d.evaluations)
        .sum::<usize>();
    assert_eq!(owner.evaluations(), work);
    assert!(work > 0);
    assert!(
        PreparedCompactWallRadiation::new(
            boundary,
            HZ,
            &ports,
            "finite wall",
            CompactWallRadiationLimits {
                maximum_total_evaluations: work - 1,
                ..limits()
            },
            &mut || false
        )
        .is_err()
    );
    let mut bad = ports.clone();
    bad[1].position_m = bad[0].position_m;
    assert!(
        PreparedCompactWallRadiation::new(boundary, HZ, &bad, "overlap", limits(), &mut || false)
            .is_err()
    );
    let mut bad = ports.clone();
    bad[0].position_m[2] = -bad[0].radius_m;
    assert!(
        PreparedCompactWallRadiation::new(
            boundary,
            HZ,
            &bad,
            "wall clearance",
            limits(),
            &mut || false
        )
        .is_err()
    );
    assert!(
        PreparedCompactWallRadiation::new(
            boundary,
            HZ,
            &ports,
            "ka",
            CompactWallRadiationLimits {
                maximum_ka: 1e-4,
                ..limits()
            },
            &mut || false
        )
        .is_err()
    );
    assert!(matches!(
        PreparedCompactWallRadiation::new(boundary, HZ, &ports, "cancel", limits(), &mut || true),
        Err(CompactWallRadiationError::Cancelled)
    ));
    let active = CoupledTermination {
        omega: TAU * HZ,
        labels: vec!["a".into(), "b".into()],
        z: vec![C64::ONE, C64::from_re(2.0), C64::from_re(2.0), C64::ONE],
        source_id: "active mutual matrix".into(),
    };
    assert!(active.validate().is_err());
}
