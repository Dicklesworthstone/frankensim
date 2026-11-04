//! Physical wavelength, coherent boundary work and original material controls.

use super::*;
use crate::network::scattering::{PreparedNetworkScattering, WavePort};
use crate::network::{
    JunctionCorrections, NetworkEdge, NetworkLoad, NetworkNode, NodeBoundary, PortExcitation,
};
use crate::{Duct, DuctError, Segment, Termination};
use fs_la::eigen_complex::lu_complex;
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_math::c64::C64;
use std::cell::Cell;

fn liquid() -> FluidAcousticState {
    // Explicit analytic coefficients, not a measured-water claim.
    FluidAcousticState::try_new(1000.0, 1500.0).unwrap()
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
        corrections: JunctionCorrections::ChainLocal,
    }
}

fn line(radius: f64, length: f64) -> DuctNetwork {
    DuctNetwork {
        nodes: vec![port("input"), port("output")],
        edges: vec![edge("pipe", 0, 1, radius, length)],
    }
}

fn reference(state: &FluidAcousticState, radius: f64) -> f64 {
    state.characteristic_impedance_pa_s_m() / (core::f64::consts::PI * radius * radius)
}

fn pressure(label: &str, open: C64, impedance: C64) -> PortSource {
    PortSource {
        label: label.into(),
        excitation: PortExcitation::Pressure {
            open_circuit_pressure_pa: open,
            impedance_pa_s_m3: impedance,
        },
    }
}

fn near(actual: C64, expected: C64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 2e-10 * scale.max(expected.abs()).max(1e-30),
        "actual {actual:?}, expected {expected:?}, scale {scale:e}"
    );
}

fn balanced(response: &NetworkSourceResponse) {
    let fields = &response.fields;
    let scale = response
        .generator_powers_w
        .iter()
        .map(|p| p.abs())
        .sum::<f64>()
        + fields.port_powers_w.iter().map(|p| p.abs()).sum::<f64>();
    near(
        C64::from_re(response.generator_power_w),
        C64::from_re(
            response.source_dissipation_w + fields.external_load_power_w + fields.network_power_w,
        ),
        scale,
    );
    assert!(response.source_dissipation_w >= -2e-10 * scale);
}

#[test]
fn material_scattering_keeps_actual_wavelength_power_and_half_wave_poles() {
    let (radius, length) = (0.008, 0.2);
    let net = line(radius, length);
    let omega = 7000.0;
    for state in [
        liquid(),
        FluidAcousticState::try_new(850.0, 1200.0).unwrap(),
    ] {
        let r = reference(&state, radius);
        let ports = ["input", "output"].map(|label| WavePort {
            label: label.into(),
            reference_impedance_pa_s_m3: r,
        });
        for frequency in [
            omega,
            core::f64::consts::PI * state.sound_speed_m_s() / length,
        ] {
            let response = PreparedNetworkScattering::from_material(
                &net,
                &state,
                frequency,
                &ports,
                &[],
                None,
            )
            .unwrap();
            let phase = frequency * length / state.sound_speed_m_s();
            let delay = C64::new(phase.cos(), phase.sin());
            for (&actual, expected) in
                response
                    .pressure_matrix()
                    .iter()
                    .zip([C64::ZERO, delay, delay, C64::ZERO])
            {
                near(actual, expected, 1.0);
            }
            let incident = [C64::new(0.8, 0.2), C64::new(-0.3, 0.6)];
            let mut outgoing = [C64::ZERO; 2];
            response
                .scatter_power_normalized_into(&incident, &mut outgoing)
                .unwrap();
            near(
                C64::from_re(outgoing.iter().map(|z| z.norm_sq()).sum()),
                C64::from_re(incident.iter().map(|z| z.norm_sq()).sum()),
                1.0,
            );
        }
        let frequency = core::f64::consts::PI * state.sound_speed_m_s() / length;
        let open = C64::new(2.0, -0.4);
        let sources = [pressure("input", open, C64::from_re(r))];
        let loads = [PortTermination {
            label: "output".into(),
            load: NetworkLoad::Impedance(C64::from_re(r)),
        }];
        // The closed-circuit line's open-port impedance has a true pole here;
        // the matched pressure generator and receiver still give finite waves.
        let result = drive_sources(&net, &state, frequency, &sources, None, &loads, None).unwrap();
        near(result.fields.port_pressures[0], open.scale(0.5), open.abs());
        near(
            result.fields.port_pressures[1],
            open.scale(-0.5),
            open.abs(),
        );
        near(
            result.fields.port_inflows[0],
            open.scale(0.5 / r),
            open.abs() / r,
        );
        near(
            result.fields.port_inflows[1],
            open.scale(0.5 / r),
            open.abs() / r,
        );
        balanced(&result);
        let flow = [PortSource {
            label: "input".into(),
            excitation: PortExcitation::VolumeFlow(C64::ONE),
        }];
        let closed = [PortTermination {
            label: "output".into(),
            load: NetworkLoad::Closed,
        }];
        assert!(matches!(
            drive_sources(&net, &state, frequency, &flow, None, &closed, None),
            Err(FluidDuctError::Duct(DuctError::Singular))
        ));
    }
}

fn star() -> DuctNetwork {
    let mut nodes = vec![NetworkNode {
        id: "junction".into(),
        boundary: NodeBoundary::Junction,
    }];
    nodes.extend(["a", "b", "c", "d"].map(port));
    nodes.push(NetworkNode {
        id: "radiating branch".into(),
        boundary: NodeBoundary::Termination(NetworkLoad::Radiation {
            radius: 0.009,
            termination: Termination::UnflangedOpen,
        }),
    });
    let radii = [0.008, 0.011, 0.007, 0.013, 0.009];
    let lengths = [0.12, 0.17, 0.09, 0.21, 0.15];
    DuctNetwork {
        nodes,
        edges: (0..5)
            .map(|i| edge(&format!("arm-{i}"), 0, i + 1, radii[i], lengths[i]))
            .collect(),
    }
}

// Independent exact cylinder ODE expressed as nodal pressure equations.
// No production network, scattering, source or composition helper is used.
fn star_oracle(
    net: &DuctNetwork,
    state: &FluidAcousticState,
    omega: f64,
    bank: &CoupledSourceBank,
    coupled: &CoupledTermination,
    radiation: C64,
) -> (Vec<C64>, Vec<C64>) {
    let n = net.nodes.len();
    let mut y = vec![C64::ZERO; n * n];
    for edge in &net.edges {
        let Segment::Cylinder { radius, length } = edge.duct.segments[0] else {
            unreachable!()
        };
        let phase = omega * length / state.sound_speed_m_s();
        let r = reference(state, radius);
        let diagonal = C64::new(0.0, phase.cos() / (r * phase.sin()));
        let mutual = C64::new(0.0, -1.0 / (r * phase.sin()));
        for node in [edge.from, edge.to] {
            y[node * n + node] = y[node * n + node] + diagonal;
        }
        y[edge.from * n + edge.to] = y[edge.from * n + edge.to] + mutual;
        y[edge.to * n + edge.from] = y[edge.to * n + edge.from] + mutual;
    }
    let mut a = y.clone();
    let mut rhs = vec![C64::ZERO; n];
    for (first, impedance) in [(1, &bank.impedance.z), (3, &coupled.z)] {
        for row in 0..2 {
            for col in 0..n {
                a[(first + row) * n + col] = (0..2).fold(C64::ZERO, |z, other| {
                    z + impedance[row * 2 + other] * y[(first + other) * n + col]
                });
            }
            a[(first + row) * n + first + row] = a[(first + row) * n + first + row] + C64::ONE;
        }
    }
    rhs[1..3].copy_from_slice(&bank.open_circuit_pressures_pa);
    for col in 0..n {
        a[5 * n + col] = radiation * y[5 * n + col];
    }
    a[5 * n + 5] = a[5 * n + 5] + C64::ONE;
    for (row, value) in rhs.iter_mut().enumerate() {
        let scale = a[row * n..(row + 1) * n]
            .iter()
            .map(|z| z.abs())
            .fold(0.0_f64, f64::max);
        for z in &mut a[row * n..(row + 1) * n] {
            *z = z.scale(1.0 / scale);
        }
        *value = value.scale(1.0 / scale);
    }
    lu_complex(&a, n).unwrap().solve(&mut rhs);
    let currents = y
        .chunks_exact(n)
        .map(|row| {
            row.iter()
                .zip(&rhs)
                .fold(C64::ZERO, |sum, (&a, &p)| sum + a * p)
        })
        .collect();
    (rhs, currents)
}

#[test]
fn material_branches_keep_mutual_sources_loads_and_radiated_work() {
    let state = liquid();
    let net = star();
    let omega = 12.0 * state.sound_speed_m_s();
    let r = reference(&state, 0.008);
    let bank = CoupledSourceBank {
        impedance: CoupledTermination {
            omega,
            labels: vec!["a".into(), "b".into()],
            z: [
                C64::new(0.8, 0.2),
                C64::new(-0.18, 0.35),
                C64::new(-0.18, 0.35),
                C64::new(0.7, -0.1),
            ]
            .map(|z| z.scale(r))
            .to_vec(),
            source_id: "analytic mutual sources".into(),
        },
        open_circuit_pressures_pa: vec![C64::new(1.0, 0.3), C64::new(-0.2, 0.6)],
    };
    let coupled = CoupledTermination {
        omega,
        labels: vec!["c".into(), "d".into()],
        z: [
            C64::new(0.7, -0.3),
            C64::new(0.16, 0.22),
            C64::new(0.16, 0.22),
            C64::new(0.6, 0.1),
        ]
        .map(|z| z.scale(r))
        .to_vec(),
        source_id: "analytic reciprocal exterior".into(),
    };
    let radiation = crate::compact_radiation_impedance_from_material(
        Termination::UnflangedOpen,
        &state,
        0.009,
        omega,
    )
    .unwrap();
    let (p, u) = star_oracle(&net, &state, omega, &bank, &coupled, radiation);
    let result = drive_sources(&net, &state, omega, &[], Some(&bank), &[], Some(&coupled)).unwrap();
    for i in 0..4 {
        near(result.fields.port_pressures[i], p[i + 1], 2.0);
        near(result.fields.port_inflows[i], u[i + 1], 2.0 / r);
    }
    let radiated = -0.5 * (p[5] * u[5].conj()).re;
    assert!(radiated > 0.0);
    near(
        C64::from_re(result.fields.network_power_w),
        C64::from_re(radiated),
        1.0 / r,
    );
    balanced(&result);

    // Material and stable edge identities are retained under authored order.
    let mut reordered = net.clone();
    reordered.edges.reverse();
    let again = drive_sources(
        &reordered,
        &state,
        omega,
        &[],
        Some(&bank),
        &[],
        Some(&coupled),
    )
    .unwrap();
    assert_eq!(result.fields.port_pressures, again.fields.port_pressures);
    assert_eq!(result.fields.port_inflows, again.fields.port_inflows);
}

struct CountedFluid {
    state: ThermoviscousFluidState,
    acoustic: Cell<usize>,
    complete: Cell<usize>,
}
impl FluidAcoustics for CountedFluid {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        self.acoustic.set(self.acoustic.get() + 1);
        Ok(self.state.acoustic())
    }
}
impl FluidTransport for CountedFluid {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        panic!("network must retain one complete snapshot, not resample transport")
    }
}
impl ThermoviscousFluid for CountedFluid {
    fn thermoviscous_state(&self) -> Result<ThermoviscousFluidState, MaterialError> {
        self.complete.set(self.complete.get() + 1);
        Ok(self.state)
    }
}

#[test]
fn material_wall_loss_keeps_complex_travel_and_one_original_snapshot() {
    let acoustic = FluidAcousticState::try_new(950.0, 1300.0).unwrap();
    let transport = FluidTransportState::try_new(950.0, 4200.0, 0.6, 0.2).unwrap();
    let state = ThermoviscousFluidState::try_new(acoustic, transport, 4100.0).unwrap();
    let counted = CountedFluid {
        state,
        acoustic: Cell::new(0),
        complete: Cell::new(0),
    };
    let (radius, length, omega) = (0.0003, 0.8, 50_000.0);
    let net = line(radius, length);
    let wave =
        crate::fluid::segment_wave_with_loss(&state, radius, omega, LossModel::Bessel).unwrap();
    let open = C64::new(1.0, -0.3);
    let sources = [pressure("input", open, wave.characteristic_impedance)];
    let loads = [PortTermination {
        label: "output".into(),
        load: NetworkLoad::Impedance(wave.characteristic_impedance),
    }];
    let prepared = prepare_sources_with_loss(
        &net,
        &counted,
        omega,
        LossModel::Bessel,
        ConeSlices::Fixed(8),
        &sources,
        None,
        &loads,
        None,
    )
    .unwrap();
    assert_eq!((counted.acoustic.get(), counted.complete.get()), (0, 1));
    let attenuation = (-wave.wavenumber.im * length).exp();
    let phase = wave.wavenumber.re * length;
    let transmitted = open.scale(0.5) * C64::new(phase.cos(), phase.sin()).scale(attenuation);
    let mut result = prepared.response_buffer();
    for factor in [C64::ONE, C64::new(-0.2, 0.7)] {
        prepared.drive_into(&[open * factor], &mut result).unwrap();
        near(
            result.fields.port_pressures[0],
            open.scale(0.5) * factor,
            open.abs(),
        );
        near(
            result.fields.port_pressures[1],
            transmitted * factor,
            transmitted.abs(),
        );
        near(
            result.fields.port_inflows[1],
            -transmitted * factor / wave.characteristic_impedance,
            transmitted.abs() / wave.characteristic_impedance.abs(),
        );
        assert!(result.fields.network_power_w > 0.0);
        balanced(&result);
    }
    assert_eq!((counted.acoustic.get(), counted.complete.get()), (0, 1));
    let r = reference(&acoustic, radius);
    let ports = ["input", "output"].map(|label| WavePort {
        label: label.into(),
        reference_impedance_pa_s_m3: r,
    });
    let scattering = PreparedNetworkScattering::from_material_with_loss(
        &net,
        &counted,
        omega,
        LossModel::Bessel,
        ConeSlices::Fixed(8),
        &ports,
        &[],
        None,
    )
    .unwrap();
    let incident = [C64::ONE, C64::new(0.1, 0.3)];
    let mut outgoing = [C64::ZERO; 2];
    scattering
        .scatter_power_normalized_into(&incident, &mut outgoing)
        .unwrap();
    assert!(
        outgoing.iter().map(|z| z.norm_sq()).sum::<f64>()
            < incident.iter().map(|z| z.norm_sq()).sum::<f64>()
    );
    assert_eq!((counted.acoustic.get(), counted.complete.get()), (0, 2));
    assert!(prepare_sources(&net, &counted, 0.0, &sources, None, &loads, None).is_err());
    assert_eq!((counted.acoustic.get(), counted.complete.get()), (0, 2));
}

#[test]
fn porous_material_boundary_keeps_same_fluid_and_its_physical_load_work() {
    use fs_material::porous::{PorousCard, PorousFluidLaw};
    let acoustic = liquid();
    let transport = FluidTransportState::try_new(1000.0, 4200.0, 0.6, 0.001).unwrap();
    let state = ThermoviscousFluidState::try_new(acoustic, transport, 4100.0).unwrap();
    let counted = CountedFluid {
        state,
        acoustic: Cell::new(0),
        complete: Cell::new(0),
    };
    let (radius, length, omega) = (0.008, 0.2, core::f64::consts::TAU * 550.0);
    let net = line(radius, length);
    let law = PorousFluidLaw::try_new(PorousCard::melamine_uf_nguyen2024(), 1.81e-5).unwrap();
    let area = core::f64::consts::PI * radius * radius;
    let surface = law
        .surface_impedance_rigid_fluid(&state, omega, 0.05)
        .unwrap();
    let load = C64::new(surface.re / area, -surface.im / area);
    let r = reference(&acoustic, radius);
    let open = C64::new(1.2, -0.4);
    let sources = [pressure("input", open, C64::from_re(r))];
    let loads = [PortTermination {
        label: "output".into(),
        load: NetworkLoad::PorousRigid {
            material: law,
            thickness_m: 0.05,
            area_m2: area,
        },
    }];
    // An acoustic snapshot alone cannot invent porous drag or heat diffusion.
    assert!(prepare_sources(&net, &acoustic, omega, &sources, None, &loads, None).is_err());
    let response = drive_sources_with_loss(
        &net,
        &counted,
        omega,
        LossModel::Lossless,
        ConeSlices::Auto,
        &sources,
        None,
        &loads,
        None,
    )
    .unwrap();
    assert_eq!((counted.acoustic.get(), counted.complete.get()), (0, 1));
    let phase = omega * length / acoustic.sound_speed_m_s();
    let delay = C64::new(phase.cos(), phase.sin());
    let reflected = (load - C64::from_re(r)) / (load + C64::from_re(r));
    let pressure = open.scale(0.5) * delay * (C64::ONE + reflected);
    let outward = open.scale(0.5 / r) * delay * (C64::ONE - reflected);
    near(response.fields.port_pressures[1], pressure, open.abs());
    near(response.fields.port_inflows[1], -outward, open.abs() / r);
    let absorbed = 0.5 * load.re * outward.norm_sq();
    assert!(absorbed > 0.0);
    near(
        C64::from_re(response.fields.external_load_power_w),
        C64::from_re(absorbed),
        1.0 / r,
    );
    near(
        C64::from_re(response.fields.network_power_w),
        C64::ZERO,
        1.0 / r,
    );
    balanced(&response);
}

#[test]
fn original_gas_material_network_matches_native_loaded_sources_and_scattering() {
    use fs_material::gas::{GasSpec, GasState};
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let state = gas.acoustic_state().unwrap();
    let net = star();
    let r = reference(&state, 0.008);
    let sources = [
        pressure("a", C64::new(1.0, 0.2), C64::new(0.6, 0.1).scale(r)),
        PortSource {
            label: "b".into(),
            excitation: PortExcitation::VolumeFlow(C64::new(0.3, -0.2).scale(1.0 / r)),
        },
    ];
    let loads = [
        PortTermination {
            label: "c".into(),
            load: NetworkLoad::PressureRelease,
        },
        PortTermination {
            label: "d".into(),
            load: NetworkLoad::Impedance(C64::new(0.8, -0.2).scale(r)),
        },
    ];
    let ports = ["a", "b"].map(|label| WavePort {
        label: label.into(),
        reference_impedance_pa_s_m3: r,
    });
    for loss in [LossModel::Lossless, LossModel::WideTube, LossModel::Bessel] {
        let omega = 1300.0;
        let native =
            crate::network::network_drive_sources(&net, &gas, omega, loss, &sources, &loads, None)
                .unwrap();
        let material = drive_sources_with_loss(
            &net,
            &gas,
            omega,
            loss,
            ConeSlices::Auto,
            &sources,
            None,
            &loads,
            None,
        )
        .unwrap();
        for i in 0..4 {
            near(
                material.fields.port_pressures[i],
                native.fields.port_pressures[i],
                2.0,
            );
            near(
                material.fields.port_inflows[i],
                native.fields.port_inflows[i],
                2.0 / r,
            );
        }
        balanced(&material);
        let native =
            PreparedNetworkScattering::new(&net, &gas, omega, loss, &ports, &loads, None).unwrap();
        let material = PreparedNetworkScattering::from_material_with_loss(
            &net,
            &gas,
            omega,
            loss,
            ConeSlices::Auto,
            &ports,
            &loads,
            None,
        )
        .unwrap();
        for (&a, &b) in material
            .pressure_matrix()
            .iter()
            .zip(native.pressure_matrix())
        {
            near(a, b, 1.0);
        }
    }
}

#[test]
fn original_water_state_changes_network_phase_and_characteristic_reference() {
    use fs_qty::{Pressure, Temperature};
    use fs_thermochem::water::if97::{If97Region, If97State};
    let net = line(0.008, 0.2);
    let mut responses = Vec::new();
    for temperature in [280.0, 350.0] {
        let water = If97State::at_pressure_temperature(
            Pressure::new(3_000_000.0),
            Temperature::new(temperature),
            If97Region::Region1,
        )
        .unwrap();
        let state = water.acoustic_state().unwrap();
        let r = reference(&state, 0.008);
        let ports = ["input", "output"].map(|label| WavePort {
            label: label.into(),
            reference_impedance_pa_s_m3: r,
        });
        let response =
            PreparedNetworkScattering::from_material(&net, &water, 7000.0, &ports, &[], None)
                .unwrap();
        let phase = 7000.0 * 0.2 / state.sound_speed_m_s();
        near(
            response.pressure_matrix()[1],
            C64::new(phase.cos(), phase.sin()),
            1.0,
        );
        responses.push((r, response.pressure_matrix()[1]));
    }
    assert!((responses[0].0 / responses[1].0 - 1.0).abs() > 0.01);
    assert!((responses[0].1 - responses[1].1).abs() > 0.01);
}
