//! Physical fields, interference and original constitutive admission.

use super::*;
use crate::network::field::{NetworkFieldBudget, NetworkFieldProbe, NetworkFieldSection};
use fs_fft::harmonic_window::HarmonicWindow;

fn pressure(label: &str, value: C64, impedance: C64) -> PortSource {
    PortSource {
        label: label.into(),
        excitation: PortExcitation::Pressure {
            open_circuit_pressure_pa: value,
            impedance_pa_s_m3: impedance,
        },
    }
}

fn request<'a>(
    omega: f64,
    sources: &'a [PortSource],
    loads: &'a [PortTermination],
    probes: &'a [NetworkFieldProbe],
) -> MolecularNetworkFieldRequest<'a> {
    MolecularNetworkFieldRequest {
        omega,
        source_clock: "original-molecular-source-clock",
        sources,
        source_bank: None,
        loads,
        coupled_load: None,
        probes,
        budget: NetworkFieldBudget::default(),
    }
}

fn probe(id: &str, x: f64) -> NetworkFieldProbe {
    NetworkFieldProbe {
        edge_id: id.into(),
        section: NetworkFieldSection::Cylinder {
            segment: 0,
            offset_m: x,
        },
    }
}

fn phase(k: C64, x: f64) -> C64 {
    C64::new((k.re * x).cos(), (k.re * x).sin()).scale((-k.im * x).exp())
}

// Direct solution of the homogeneous p/U equations. Neither a network nor
// scattering/transfer reconstruction appears in this reference.
fn analytic(k: C64, z: C64, length: f64, x: f64, a: C64, b: C64) -> [C64; 2] {
    let f = a * phase(k, x);
    let b = b * phase(k, length - x);
    [f + b, (f - b) / z]
}

#[test]
fn molecular_internal_fields_match_physical_waves_at_loaded_poles_and_opaque_lengths() {
    let radius = 0.008;
    let area = core::f64::consts::PI * radius * radius;
    let (a, b) = (C64::new(1.3, -0.4), C64::new(-0.2, 0.6));
    for (reservoirs, opaque) in [(0, false), (2, false), (2, true)] {
        let state = law(gas(293.15), reservoirs, 1.0);
        let omega = if reservoirs == 0 {
            core::f64::consts::PI * state.gas().sound_speed / 0.27
        } else {
            17_000.0
        };
        let k = state.relaxing_wavenumber(omega).unwrap();
        let z = C64::from_re(omega * state.gas().density / area) / k;
        let length = if opaque { 600.0 / k.im } else { 0.27 };
        let net = DuctNetwork {
            nodes: vec![port("left"), port("right")],
            edges: vec![edge("tube", 0, 1, radius, length)],
        };
        let original = model(net, &[state.clone()], None);
        // Matched physical generators select A at x=0 and B at x=L, including
        // the loaded half-wave pole; the reference impedance is still real.
        let sources = [
            pressure("right", b.scale(2.0), z),
            pressure("left", a.scale(2.0), z),
        ];
        let probes: Vec<_> = [0.0, 0.2, 0.5, 0.8, 1.0]
            .map(|f| probe("tube", f * length))
            .into();
        let owner = original
            .prepare_fields(&request(omega, &sources, &[], &probes), || false)
            .unwrap();
        assert!(core::ptr::eq(owner.network(), &original));
        assert_eq!(owner.network().edge("tube").unwrap().state(), &state);
        let fields = owner.fields();
        let values = fields
            .sample(fields.port_response().default_excitations(), || false)
            .unwrap();
        let ports = original
            .prepare_sources(omega, &sources, None, &[], None)
            .unwrap();
        assert_eq!(
            fields.port_response().pressure_transfer(),
            ports.pressure_transfer()
        );
        assert_eq!(
            fields.port_response().flow_transfer(),
            ports.flow_transfer()
        );
        for (value, location) in values.iter().zip(fields.locations()) {
            let x = location.axial_position_m;
            let [p, u] = analytic(k, z, length, x, a, b);
            near(value.pressure_pa, p, p.abs().max(1e-280));
            near(value.volume_flow_m3_s, u, u.abs().max(1e-280));
            near(
                value.particle_velocity_m_s,
                u.scale(1.0 / area),
                u.abs().max(1e-280) / area,
            );
            let expected = 0.5 * (p * u.conj()).re;
            assert!(
                (value.mean_power_w - expected).abs() < 2e-10 * (p.abs() * u.abs()).max(1e-280)
            );
        }
        assert!(
            fields
                .prepare_cylinder_absorption(&[[0, 4]], 1, || false)
                .is_err()
        );
    }
}

#[test]
fn molecular_branch_fields_retain_regional_materials_full_source_bank_and_coherence() {
    let states = [law(gas(293.15), 2, 1.0), law(gas(343.15), 2, 0.2)];
    let lengths = [0.11, 0.17];
    let radius = 0.008;
    let net = DuctNetwork {
        nodes: vec![port("left"), port("right"), port("receiver")],
        edges: vec![
            edge("left-arm", 0, 2, radius, lengths[0]),
            edge("right-arm", 1, 2, radius, lengths[1]),
        ],
    };
    let original = model(net.clone(), &states, None);
    let mut shuffled = net.clone();
    shuffled.edges.reverse();
    let shuffled = MolecularNetwork::new(
        shuffled,
        records(&net, &states),
        bindings(&net, &states, None),
        16,
    )
    .unwrap();
    let omega = 5107.0;
    let r = states[0].gas().characteristic_impedance / (core::f64::consts::PI * radius * radius);
    let source_z = [
        C64::new(0.7, -0.2),
        C64::new(-0.35, 0.6),
        C64::new(-0.35, 0.6),
        C64::new(0.9, 0.4),
    ]
    .map(|z| z.scale(r));
    let load_z = C64::new(0.8, -0.3).scale(r);
    let loads = [PortTermination {
        label: "receiver".into(),
        load: NetworkLoad::Impedance(load_z),
    }];
    let bank = CoupledSourceBank {
        impedance: CoupledTermination {
            omega,
            labels: vec!["left".into(), "right".into()],
            z: source_z.to_vec(),
            source_id: "original-synthetic-reciprocal-source-pair".into(),
        },
        open_circuit_pressures_pa: vec![C64::ZERO; 2],
    };
    let probes: Vec<_> = ["left-arm", "right-arm"]
        .iter()
        .enumerate()
        .flat_map(|(i, id)| [0.0, 0.4, 1.0].map(|f| probe(id, f * lengths[i])))
        .collect();
    let request = MolecularNetworkFieldRequest {
        source_bank: Some(&bank),
        ..request(omega, &[], &loads, &probes)
    };
    let owner = original.prepare_fields(&request, || false).unwrap();
    let reordered = shuffled.prepare_fields(&request, || false).unwrap();
    let fields = owner.fields();
    assert_eq!(
        fields.pressure_transfer(),
        reordered.fields().pressure_transfer()
    );
    assert_eq!(fields.flow_transfer(), reordered.fields().flow_transfer());
    let m = [
        line(&states[0], omega, radius, lengths[0]),
        line(&states[1], omega, radius, lengths[1]),
    ];
    let receiver_z = (load_z.recip() + m[0][2] / m[0][3] + m[1][2] / m[1][3]).recip();
    let mut matrix = source_z;
    for i in 0..2 {
        for j in 0..2 {
            matrix[2 * i + j] = matrix[2 * i + j] + receiver_z / (m[i][3] * m[j][3]);
        }
        matrix[2 * i + i] = matrix[2 * i + i] + m[i][1] / m[i][3];
    }
    let determinant = matrix[0] * matrix[3] - matrix[1] * matrix[2];
    let mut received = Vec::new();
    for angle in [0.0_f64, 0.8, core::f64::consts::PI] {
        let amplitude = [C64::new(1.2, -0.3), C64::new(angle.cos(), angle.sin())];
        let u = [
            (matrix[3] * amplitude[0] - matrix[1] * amplitude[1]) / determinant,
            (matrix[0] * amplitude[1] - matrix[2] * amplitude[0]) / determinant,
        ];
        let receiver = receiver_z * (u[0] / m[0][3] + u[1] / m[1][3]);
        let values = fields.sample(&amplitude, || false).unwrap();
        let moments = fields
            .observe_window(
                HarmonicWindow::new(0.0017, 0.00037).unwrap(),
                None,
                &[&amplitude],
                || false,
            )
            .unwrap();
        for i in 0..2 {
            let q_end = (u[i] - m[i][2] * receiver) / m[i][3];
            for j in 0..3 {
                let offset = lengths[i] * [0.0, 0.4, 1.0][j];
                let local = line(&states[i], omega, radius, lengths[i] - offset);
                let p = local[0] * receiver + local[1] * q_end;
                let q = local[2] * receiver + local[3] * q_end;
                let index = 3 * i + j;
                near(values[index].pressure_pa, p, 4.0);
                near(values[index].volume_flow_m3_s, q, 4.0 / r);
                assert_eq!(
                    fields.locations()[index].density_kg_m3,
                    states[i].gas().density
                );
                assert!((moments[index].mean_square_pressure_pa2 - 0.5 * p.norm_sq()).abs() < 2e-9);
                assert!((moments[index].mean_power_w - 0.5 * (p * q.conj()).re).abs() < 2e-9 / r);
                assert_eq!(moments[index].mean_pressure_pa, 0.0);
            }
        }
        near(values[2].pressure_pa, values[5].pressure_pa, 4.0);
        near(
            values[2].volume_flow_m3_s + values[5].volume_flow_m3_s,
            receiver / load_z,
            4.0 / r,
        );
        received.push(moments[2].mean_square_pressure_pa2);
    }
    assert!(
        (received[0] - received[2]).abs() > 1e-4,
        "coherent source phases must survive inside branches"
    );
}

#[test]
fn joint_molecular_wall_profile_fields_retain_original_diffusion_and_native_panel_cuts() {
    let state = law(gas(293.15), 2, 1.0);
    let gas = state.gas();
    let (radius, length, omega) = (0.001, 0.13, 5107.0);
    let resolution = MolecularWallResolution {
        poles_per_field: 8,
        band_rad_s: [0.1, 7000.0],
        maximum_relative_error: 0.02,
    };
    let wall = fs_phs::ZkDiffusive::new(
        radius,
        gas.density,
        gas.dynamic_viscosity,
        gas.gamma,
        gas.prandtl,
        8,
        resolution.band_rad_s,
        0.02,
    )
    .unwrap();
    let (cp, bulk_k) = crate::molecular::local::caloric(&state, omega).unwrap();
    let wave =
        crate::molecular::local::wave(&state, radius, 0.0, Some(&wall), omega, cp, bulk_k).unwrap();
    let (a, b) = (C64::new(0.7, -0.3), C64::new(-0.2, 0.4));
    let profile = crate::RadiusProfile::exponential(radius, 0.0, 0.04).unwrap();
    let mut tube = edge("profile", 0, 1, radius, length);
    // Three authored segments remain one physical uniform tube, including
    // all 16 requested cone slices and eight original Magnus panels.
    tube.duct.segments = vec![
        Segment::Cylinder {
            radius,
            length: 0.04,
        },
        Segment::Cone {
            inlet_radius: radius,
            outlet_radius: radius,
            length: 0.05,
        },
        Segment::Profile {
            profile: std::sync::Arc::new(profile),
            steps: crate::RiccatiSteps(8),
        },
    ];
    let net = DuctNetwork {
        nodes: vec![port("left"), port("right")],
        edges: vec![tube],
    };
    let original = model(net, &[state], Some(resolution));
    let sources = [
        pressure("left", a.scale(2.0), wave.characteristic_impedance),
        pressure("right", b.scale(2.0), wave.characteristic_impedance),
    ];
    let mut probes: Vec<_> = [
        NetworkFieldSection::Inlet { segment: 0 },
        NetworkFieldSection::Cylinder {
            segment: 0,
            offset_m: 0.017,
        },
        NetworkFieldSection::Outlet { segment: 0 },
        NetworkFieldSection::Inlet { segment: 1 },
        NetworkFieldSection::Cone {
            segment: 1,
            offset_m: 0.023,
        },
        NetworkFieldSection::Outlet { segment: 1 },
    ]
    .map(|section| NetworkFieldProbe {
        edge_id: "profile".into(),
        section,
    })
    .into();
    probes.extend([0, 1, 4, 7, 8].map(|boundary| NetworkFieldProbe {
        edge_id: "profile".into(),
        section: NetworkFieldSection::Profile {
            segment: 2,
            boundary,
        },
    }));
    let owner = original
        .prepare_fields(&request(omega, &sources, &[], &probes), || false)
        .unwrap();
    let fields = owner.fields();
    let values = fields
        .sample(fields.port_response().default_excitations(), || false)
        .unwrap();
    for (value, location) in values.iter().zip(fields.locations()) {
        let [p, q] = analytic(
            wave.wavenumber,
            wave.characteristic_impedance,
            length,
            location.axial_position_m,
            a,
            b,
        );
        near(value.pressure_pa, p, 2.0);
        near(
            value.volume_flow_m3_s,
            q,
            2.0 / wave.characteristic_impedance.abs(),
        );
    }
    // This joint response differs from bare bulk relaxation; the test would
    // fail if observations silently dropped the original wall memory.
    assert!((wave.wavenumber - bulk_k).abs() > 0.01);
    assert!(
        original
            .prepare_fields(&request(7500.0, &sources, &[], &probes), || false)
            .is_err()
    );
}

#[test]
fn molecular_cone_fields_match_independent_spherical_pressure_and_flow() {
    let state = law(gas(293.15), 2, 1.0);
    let (length, omega) = (0.37, 1700.0);
    let k = state.relaxing_wavenumber(omega).unwrap();
    let z = C64::from_re(omega * state.gas().density) / k;
    let (a, b) = (C64::new(1.1, 0.2), C64::new(-0.35, 0.4));
    for (inlet_radius, outlet_radius) in [(0.008, 0.015), (0.015, 0.008)] {
        let slope = (outlet_radius - inlet_radius) / length;
        // Independent spherical-wave solution p=(A exp(ikx)+B exp(-ikx))/r.
        // Momentum gives U=A_section*p'/(i*k*Z_specific), retaining the
        // spherical 1/r induction term as well as the traveling-wave term.
        let physical = |x: f64| {
            let r = inlet_radius + slope * x;
            let f = a.scale(inlet_radius / r) * phase(k, x);
            let b = b.scale(outlet_radius / r) * phase(k, length - x);
            let p = f + b;
            let q =
                (f - b + p * C64::new(0.0, slope / r) / k).scale(core::f64::consts::PI * r * r) / z;
            [p, q]
        };
        let mut tube = edge("cone", 0, 1, inlet_radius, length);
        tube.duct.segments = vec![Segment::Cone {
            inlet_radius,
            outlet_radius,
            length,
        }];
        let net = DuctNetwork {
            nodes: vec![port("left"), port("right")],
            edges: vec![tube],
        };
        let original = model(net, &[state.clone()], None);
        let sources = [
            pressure("right", physical(length)[0], C64::ZERO),
            pressure("left", physical(0.0)[0], C64::ZERO),
        ];
        let probes: Vec<_> = [0.0, 0.19, 0.53, 0.84, 1.0]
            .map(|f| NetworkFieldProbe {
                edge_id: "cone".into(),
                section: NetworkFieldSection::Cone {
                    segment: 0,
                    offset_m: length * f,
                },
            })
            .into();
        let owner = original
            .prepare_fields(&request(omega, &sources, &[], &probes), || false)
            .unwrap();
        let fields = owner.fields();
        let values = fields
            .sample(fields.port_response().default_excitations(), || false)
            .unwrap();
        for (value, location) in values.iter().zip(fields.locations()) {
            let [p, q] = physical(location.axial_position_m);
            near(value.pressure_pa, p, 4.0);
            near(value.volume_flow_m3_s, q, 4.0 * location.area_m2 / z.abs());
        }
    }
}

#[test]
fn molecular_fields_refuse_budget_geometry_cancellation_and_ignored_mechanical_wall() {
    let state = law(gas(293.15), 2, 1.0);
    let net = DuctNetwork {
        nodes: vec![port("drive"), port("load")],
        edges: vec![edge("tube", 0, 1, 0.008, 0.27)],
    };
    let original = model(net.clone(), &[state.clone()], None);
    let sources = [PortSource {
        label: "drive".into(),
        excitation: PortExcitation::VolumeFlow(C64::from_re(1e-7)),
    }];
    let loads = [PortTermination {
        label: "load".into(),
        load: NetworkLoad::Closed,
    }];
    let probes = [probe("tube", 0.07)];
    let mut request = request(1700.0, &sources, &loads, &probes);
    request.budget.maximum_transfer_entries = 5; // one source needs 2+4.
    assert!(original.prepare_fields(&request, || false).is_err());
    request.budget = NetworkFieldBudget::default();
    request.source_clock = "";
    assert!(original.prepare_fields(&request, || false).is_err());
    request.source_clock = "original-molecular-source-clock";
    let invalid = [probe("tube", 0.3)];
    request.probes = &invalid;
    assert!(original.prepare_fields(&request, || false).is_err());
    request.probes = &probes;
    let mut calls = 0;
    assert!(matches!(
        original.prepare_fields(&request, || {
            calls += 1;
            calls == 5
        }),
        Err(MolecularDuctError::Duct(DuctError::Cancelled { .. }))
    ));
    let first = original.prepare_fields(&request, || false).unwrap();
    let replay = original.prepare_fields(&request, || false).unwrap();
    assert_eq!(
        first.fields().pressure_transfer(),
        replay.fields().pressure_transfer()
    );
    let media = records(&net, &[state.clone()])
        .with_edge_walls(
            &net,
            vec![(
                "tube".into(),
                fs_phs::WallPin {
                    surface_density: 1.2,
                    stiffness_per_area: 1e9,
                    resistance: 150.0,
                },
            )],
        )
        .unwrap();
    assert!(matches!(
        MolecularNetwork::new(net.clone(), media, bindings(&net, &[state], None), 16),
        Err(MolecularDuctError::Duct(DuctError::BadParameter {
            what: "molecular network mechanical wall pins require their original coupled owner"
        }))
    ));
}
