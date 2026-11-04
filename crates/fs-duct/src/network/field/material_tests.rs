//! G1/G3: real material coefficients in native fields and finite-window work.
//! Authored acoustic/transport tuples establish implementation behavior only.

use super::*;
use crate::fluid::FluidDuctError;
use crate::{ConeSlices, Duct, Segment};
use fs_material::MaterialError;
use fs_material::fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState};
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use std::cell::Cell;

fn line(radius: f64, length: f64) -> DuctNetwork {
    DuctNetwork {
        nodes: ["left", "right"]
            .map(|label| NetworkNode {
                id: label.into(),
                boundary: NodeBoundary::Port {
                    label: label.into(),
                },
            })
            .to_vec(),
        edges: vec![NetworkEdge {
            id: "liquid tube".into(),
            from: 0,
            to: 1,
            duct: Duct {
                segments: vec![Segment::Cylinder { radius, length }],
            },
            corrections: JunctionCorrections::ChainLocal,
        }],
    }
}

fn points(length: f64) -> [NetworkFieldProbe; 5] {
    core::array::from_fn(|i| NetworkFieldProbe {
        edge_id: "liquid tube".into(),
        section: match i {
            0 => NetworkFieldSection::Inlet { segment: 0 },
            4 => NetworkFieldSection::Outlet { segment: 0 },
            _ => NetworkFieldSection::Cylinder {
                segment: 0,
                offset_m: length * [0.0, 0.2, 0.5, 0.8, 1.0][i],
            },
        },
    })
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

fn request<'a>(
    omega: f64,
    loss: LossModel,
    sources: &'a [PortSource],
    loads: &'a [PortTermination],
    probes: &'a [NetworkFieldProbe],
) -> NetworkFieldRequest<'a> {
    NetworkFieldRequest {
        omega,
        loss,
        source_clock: "material-field-clock",
        sources,
        source_bank: None,
        loads,
        coupled_load: None,
        probes,
        budget: NetworkFieldBudget::default(),
    }
}

fn frequency<'a>(
    omega: f64,
    loss: LossModel,
    sources: &'a [PortSource],
    loads: &'a [PortTermination],
) -> NetworkFieldFrequency<'a> {
    NetworkFieldFrequency {
        omega,
        loss,
        sources,
        source_bank: None,
        loads,
        coupled_load: None,
    }
}

fn spectrum_request<'a>(
    frequencies: &'a [NetworkFieldFrequency<'a>],
    probes: &'a [NetworkFieldProbe],
) -> NetworkFieldSpectrumRequest<'a> {
    NetworkFieldSpectrumRequest {
        source_clock: "material-field-clock",
        probes,
        frequencies,
        field_budget: NetworkFieldBudget::default(),
        spectrum_budget: NetworkFieldSpectrumBudget {
            maximum_frequencies: 8,
            maximum_retained_transfer_entries: 4096,
        },
    }
}

fn near(actual: C64, expected: C64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 4e-10 * scale.max(expected.abs()).max(1e-30),
        "{actual:?} != {expected:?}; scale {scale:e}"
    );
}

fn travel(k: C64, x: f64) -> C64 {
    let (s, c) = (k.re * x).sin_cos();
    C64::new(c, s).scale((-k.im * x).exp())
}

// Direct homogeneous-cylinder ODE solution, with waves referenced to the two
// physical ends. This uses no field reconstruction or loaded-network solver.
fn analytic(k: C64, z: C64, length: f64, x: f64, a: C64, b: C64) -> [C64; 2] {
    let forward = a * travel(k, x);
    let backward = b * travel(k, length - x);
    [forward + backward, (forward - backward) / z]
}

#[test]
fn liquid_fields_keep_counterpropagating_pressure_velocity_and_loaded_halfwave_power() {
    let (radius, length) = (0.008, 0.2);
    let net = line(radius, length);
    let probes = points(length);
    let area = core::f64::consts::PI * radius * radius;
    let (a, b) = (C64::new(1.3, -0.4), C64::new(-0.2, 0.6));
    for (rho, c) in [(1000.0, 1500.0), (850.0, 1200.0)] {
        let material = FluidAcousticState::try_new(rho, c).unwrap();
        let z = rho * c / area;
        // Matched finite-impedance generators set incident amplitudes a and b,
        // even at the exact halfwave pole of an open-circuit impedance matrix.
        let sources = [
            pressure("right", b.scale(2.0), C64::from_re(z)),
            pressure("left", a.scale(2.0), C64::from_re(z)),
        ];
        for omega in [6000.0, core::f64::consts::PI * c / length] {
            let owner = PreparedNetworkFields::from_material(
                &net,
                &material,
                &request(omega, LossModel::Lossless, &sources, &[], &probes),
                || false,
            )
            .unwrap();
            let samples = owner
                .sample(owner.port_response().default_excitations(), || false)
                .unwrap();
            let expected_power = 0.5 * (a.norm_sq() - b.norm_sq()) / z;
            for (sample, location) in samples.iter().zip(owner.locations()) {
                assert_eq!(location.density_kg_m3, rho);
                assert_eq!(location.sound_speed_m_s, c);
                let [p, u] = analytic(
                    C64::from_re(omega / c),
                    C64::from_re(z),
                    length,
                    location.axial_position_m,
                    a,
                    b,
                );
                near(sample.pressure_pa, p, 2.0);
                near(sample.volume_flow_m3_s, u, 2.0 / z);
                near(
                    sample.particle_velocity_m_s,
                    u.scale(1.0 / area),
                    2.0 / (rho * c),
                );
                near(
                    C64::from_re(sample.mean_power_w),
                    C64::from_re(expected_power),
                    2.0 / z,
                );
            }
            let mut boundary = owner.port_response().response_buffer();
            owner
                .port_response()
                .drive_into(owner.port_response().default_excitations(), &mut boundary)
                .unwrap();
            near(
                samples[0].volume_flow_m3_s,
                boundary.fields.port_inflows[0],
                2.0 / z,
            );
            near(
                samples[4].volume_flow_m3_s,
                -boundary.fields.port_inflows[1],
                2.0 / z,
            );
        }
    }
}

struct CountedFluid {
    state: ThermoviscousFluidState,
    acoustic: Cell<usize>,
    complete: Cell<usize>,
}
impl CountedFluid {
    fn new() -> Self {
        let acoustic = FluidAcousticState::try_new(950.0, 1300.0).unwrap();
        let transport = FluidTransportState::try_new(950.0, 4200.0, 0.6, 0.2).unwrap();
        Self {
            state: ThermoviscousFluidState::try_new(acoustic, transport, 4100.0).unwrap(),
            acoustic: Cell::new(0),
            complete: Cell::new(0),
        }
    }
}
impl FluidAcoustics for CountedFluid {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        self.acoustic.set(self.acoustic.get() + 1);
        Ok(self.state.acoustic())
    }
}
impl FluidTransport for CountedFluid {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        panic!("fields must retain one complete state, not resample its transport")
    }
}
impl ThermoviscousFluid for CountedFluid {
    fn thermoviscous_state(&self) -> Result<ThermoviscousFluidState, MaterialError> {
        self.complete.set(self.complete.get() + 1);
        Ok(self.state)
    }
}

#[test]
fn full_material_spectrum_uses_one_snapshot_and_reuses_original_complex_waves() {
    let material = CountedFluid::new();
    let (radius, length) = (0.0003, 0.08);
    let net = line(radius, length);
    let probes = points(length);
    let area = core::f64::consts::PI * radius * radius;
    let omegas = [5000.0, 11_000.0];
    let (a, b) = (C64::new(1.0, -0.3), C64::new(0.2, 0.4));
    for loss in [LossModel::Lossless, LossModel::Bessel] {
        // Original constitutive owner supplies k and Z; the independent
        // propagation oracle is the decaying two-wave cylinder solution.
        let waves = omegas.map(|omega| {
            crate::fluid::segment_wave_with_loss(&material.state, radius, omega, loss).unwrap()
        });
        let sources: [[PortSource; 2]; 2] = core::array::from_fn(|i| {
            [
                pressure("left", a.scale(2.0), waves[i].characteristic_impedance),
                pressure("right", b.scale(2.0), waves[i].characteristic_impedance),
            ]
        });
        let frequencies: [_; 2] =
            core::array::from_fn(|i| frequency(omegas[i], loss, &sources[i], &[]));
        let query = spectrum_request(&frequencies, &probes);
        let owner = if loss == LossModel::Lossless {
            let original: &dyn FluidAcoustics = &material;
            PreparedNetworkFieldSpectrum::from_material(&net, original, &query, || false).unwrap()
        } else {
            let original: &dyn ThermoviscousFluid = &material;
            PreparedNetworkFieldSpectrum::from_material_with_loss(
                &net,
                original,
                ConeSlices::Auto,
                &query,
                || false,
            )
            .unwrap()
        };
        assert_eq!(owner.frequency_count(), 2);
        let expected_calls = if loss == LossModel::Lossless {
            (1, 0)
        } else {
            (1, 1)
        };
        assert_eq!(
            (material.acoustic.get(), material.complete.get()),
            expected_calls
        );
        for (i, wave) in waves.iter().enumerate() {
            let field = owner.frequency_fields(i).unwrap();
            for gain in [C64::ONE, C64::new(-0.2, 0.7)] {
                let amplitudes: Vec<_> = field
                    .port_response()
                    .default_excitations()
                    .iter()
                    .map(|&value| value * gain)
                    .collect();
                let samples = field.sample(&amplitudes, || false).unwrap();
                let z = wave.characteristic_impedance.abs();
                for (sample, location) in samples.iter().zip(field.locations()) {
                    let [p, u] = analytic(
                        wave.wavenumber,
                        wave.characteristic_impedance,
                        length,
                        location.axial_position_m,
                        a * gain,
                        b * gain,
                    );
                    near(sample.pressure_pa, p, 2.0);
                    near(sample.volume_flow_m3_s, u, 2.0 / z);
                    near(
                        sample.particle_velocity_m_s,
                        u.scale(1.0 / area),
                        2.0 / (z * area),
                    );
                    near(
                        C64::from_re(sample.mean_power_w),
                        C64::from_re(0.5 * (p * u.conj()).re),
                        2.0 / z,
                    );
                }
                if loss == LossModel::Bessel {
                    assert!(wave.wavenumber.im > 0.0);
                    assert!(samples[0].mean_power_w > samples[4].mean_power_w);
                }
            }
        }
        assert_eq!(
            (material.acoustic.get(), material.complete.get()),
            expected_calls
        );
    }

    // Eight native cone subdivisions must still expose ONE authored cone
    // outlet, followed by the actual inlet of the larger cylindrical segment.
    // The matched outlet makes that last segment a single forward wave, so
    // its original boundary solution gives an independent interface oracle.
    let cone_length = 0.05;
    let cylinder_length = length - cone_length;
    let cylinder_radius = 2.0 * radius;
    let mut tapered = line(radius, length);
    tapered.edges[0].duct.segments = vec![
        Segment::Cone {
            inlet_radius: radius,
            outlet_radius: 1.5 * radius,
            length: cone_length,
        },
        Segment::Cylinder {
            radius: cylinder_radius,
            length: cylinder_length,
        },
    ];
    let boundaries = [
        NetworkFieldSection::Inlet { segment: 0 },
        NetworkFieldSection::Outlet { segment: 0 },
        NetworkFieldSection::Inlet { segment: 1 },
        NetworkFieldSection::Outlet { segment: 1 },
    ]
    .map(|section| NetworkFieldProbe {
        edge_id: "liquid tube".into(),
        section,
    });
    let inlet_wave =
        crate::fluid::segment_wave_with_loss(&material.state, radius, omegas[0], LossModel::Bessel)
            .unwrap();
    let outlet_wave = crate::fluid::segment_wave_with_loss(
        &material.state,
        cylinder_radius,
        omegas[0],
        LossModel::Bessel,
    )
    .unwrap();
    let sources = [pressure(
        "left",
        a.scale(2.0),
        inlet_wave.characteristic_impedance,
    )];
    let loads = [PortTermination {
        label: "right".into(),
        load: NetworkLoad::Impedance(outlet_wave.characteristic_impedance),
    }];
    let owner = PreparedNetworkFields::from_material_with_loss(
        &tapered,
        &material.state,
        ConeSlices::Fixed(8),
        &request(omegas[0], LossModel::Bessel, &sources, &loads, &boundaries),
        || false,
    )
    .unwrap();
    let samples = owner
        .sample(owner.port_response().default_excitations(), || false)
        .unwrap();
    let mut ports = owner.port_response().response_buffer();
    owner
        .port_response()
        .drive_into(owner.port_response().default_excitations(), &mut ports)
        .unwrap();
    let z = outlet_wave.characteristic_impedance.abs();
    let interface_p =
        ports.fields.port_pressures[1] * travel(outlet_wave.wavenumber, -cylinder_length);
    let interface_u = interface_p / outlet_wave.characteristic_impedance;
    for i in [1, 2] {
        let location = &owner.locations()[i];
        near(samples[i].pressure_pa, interface_p, 2.0);
        near(samples[i].volume_flow_m3_s, interface_u, 2.0 / z);
        near(
            samples[i].particle_velocity_m_s,
            interface_u.scale(1.0 / location.area_m2),
            2.0 / (z * location.area_m2),
        );
        near(
            C64::from_re(location.axial_position_m),
            C64::from_re(cone_length),
            length,
        );
    }
    assert!(owner.locations()[1].area_m2 < owner.locations()[2].area_m2);
    near(samples[0].pressure_pa, ports.fields.port_pressures[0], 2.0);
    near(
        samples[3].volume_flow_m3_s,
        -ports.fields.port_inflows[1],
        2.0 / z,
    );
    assert!(samples[0].mean_power_w > samples[3].mean_power_w);
}

#[test]
fn material_two_frequency_exposure_retains_beats_in_pressure_velocity_and_power() {
    let material = FluidAcousticState::try_new(1000.0, 1500.0).unwrap();
    let (radius, length) = (0.008, 0.2);
    let area = core::f64::consts::PI * radius * radius;
    let z = material.characteristic_impedance_pa_s_m() / area;
    let net = line(radius, length);
    let probes = points(length);
    let sources = [pressure("left", C64::from_re(2.0), C64::from_re(z))];
    let loads = [PortTermination {
        label: "right".into(),
        load: NetworkLoad::Impedance(C64::from_re(z)),
    }];
    let frequencies =
        [900.0, 901.0].map(|omega| frequency(omega, LossModel::Lossless, &sources, &loads));
    let owner = PreparedNetworkFieldSpectrum::from_material(
        &net,
        &material,
        &spectrum_request(&frequencies, &probes),
        || false,
    )
    .unwrap();
    let positive = [C64::from_re(2.0)];
    let negative = [-positive[0]];
    let lines = [
        NetworkFieldSpectralLine {
            frequency_index: 0,
            excitation: &positive,
        },
        NetworkFieldSpectralLine {
            frequency_index: 1,
            excitation: &negative,
        },
    ];
    let limits = NetworkFieldExposureLimits {
        maximum_lines: 8,
        maximum_factors: 4,
        maximum_transfer_products: 128,
        maximum_probe_frequency_pairs: 128,
        maximum_workspace_bytes: 65_536,
    };
    let coherent = [NetworkFieldSpectralFactor {
        lines: &lines,
        statistics: NetworkFieldFactorStatistics::ProperComplex,
    }];
    let bank = owner
        .prepare_observation(&[], &coherent, limits, || false)
        .unwrap();
    let independent = [
        NetworkFieldSpectralFactor {
            lines: &lines[..1],
            statistics: NetworkFieldFactorStatistics::ProperComplex,
        },
        NetworkFieldSpectralFactor {
            lines: &lines[1..],
            statistics: NetworkFieldFactorStatistics::ProperComplex,
        },
    ];
    let uncorrelated = owner
        .prepare_observation(&[], &independent, limits, || false)
        .unwrap();
    for (start, duration) in [(-0.05, 0.1), (0.2, 0.0037)] {
        let window = HarmonicWindow::new(start, duration).unwrap();
        let exposure = bank.observe_window(window, || false).unwrap();
        let separate = uncorrelated.observe_window(window, || false).unwrap();
        // One common proper factor multiplies exp(i*900*x/c)-exp(i*901*x/c).
        // The direct real-field ensemble mean has this analytic beat variance.
        let half_width = 0.5 * duration;
        for ((point, other), location) in exposure
            .moments()
            .iter()
            .zip(separate.moments())
            .zip(owner.locations())
        {
            let center =
                start + half_width - location.axial_position_m / material.sound_speed_m_s();
            let variance = 1.0 - half_width.sin() / half_width * center.cos();
            near(
                C64::from_re(point.mean_square_pressure_pa2),
                C64::from_re(variance),
                variance,
            );
            near(
                C64::from_re(point.mean_square_volume_flow_m6_s2),
                C64::from_re(variance / z.powi(2)),
                variance / z.powi(2),
            );
            near(
                C64::from_re(point.mean_square_particle_velocity_m2_s2),
                C64::from_re(variance / (z * area).powi(2)),
                variance / (z * area).powi(2),
            );
            near(
                C64::from_re(point.mean_power_w),
                C64::from_re(variance / z),
                variance / z,
            );
            near(
                C64::from_re(point.mean_intensity_w_m2),
                C64::from_re(variance / (z * area)),
                variance / (z * area),
            );
            near(C64::from_re(other.mean_square_pressure_pa2), C64::ONE, 1.0);
            assert!(point.mean_square_pressure_pa2 < 0.03 * other.mean_square_pressure_pa2);
        }
    }
}

#[test]
fn wrong_loss_and_initial_cancellation_refuse_before_any_material_query() {
    let material = CountedFluid::new();
    let (radius, length) = (0.008, 0.2);
    let z = material.state.acoustic().characteristic_impedance_pa_s_m()
        / (core::f64::consts::PI * radius * radius);
    let net = line(radius, length);
    let probes = points(length);
    let sources = [pressure("left", C64::from_re(2.0), C64::from_re(z))];
    let loads = [PortTermination {
        label: "right".into(),
        load: NetworkLoad::Impedance(C64::from_re(z)),
    }];
    for loss in [LossModel::Bessel, LossModel::WideTube] {
        let query = request(1000.0, loss, &sources, &loads, &probes);
        assert!(matches!(
            PreparedNetworkFields::from_material(&net, &material, &query, || false),
            Err(FluidDuctError::Duct(DuctError::BadParameter { .. }))
        ));
        let frequencies = [
            frequency(900.0, LossModel::Lossless, &sources, &loads),
            frequency(1000.0, loss, &sources, &loads),
        ];
        assert!(matches!(
            PreparedNetworkFieldSpectrum::from_material(
                &net,
                &material,
                &spectrum_request(&frequencies, &probes),
                || false,
            ),
            Err(FluidDuctError::Duct(DuctError::BadParameter { .. }))
        ));
    }
    assert_eq!((material.acoustic.get(), material.complete.get()), (0, 0));
    let query = request(1000.0, LossModel::Lossless, &sources, &loads, &probes);
    let frequencies = [
        frequency(900.0, LossModel::Lossless, &sources, &loads),
        frequency(1000.0, LossModel::Lossless, &sources, &loads),
    ];
    let spectrum = spectrum_request(&frequencies, &probes);
    for refused in [
        spectrum_request(&[], &probes),
        NetworkFieldSpectrumRequest {
            source_clock: " ",
            ..spectrum_request(&frequencies, &probes)
        },
        NetworkFieldSpectrumRequest {
            spectrum_budget: NetworkFieldSpectrumBudget {
                maximum_frequencies: 1,
                maximum_retained_transfer_entries: 4096,
            },
            ..spectrum_request(&frequencies, &probes)
        },
    ] {
        assert!(matches!(
            PreparedNetworkFieldSpectrum::from_material(&net, &material, &refused, || false,),
            Err(FluidDuctError::Duct(DuctError::BadParameter { .. }))
        ));
        assert!(matches!(
            PreparedNetworkFieldSpectrum::from_material_with_loss(
                &net,
                &material,
                ConeSlices::Auto,
                &refused,
                || false,
            ),
            Err(FluidDuctError::Duct(DuctError::BadParameter { .. }))
        ));
    }
    assert_eq!((material.acoustic.get(), material.complete.get()), (0, 0));
    assert!(matches!(
        PreparedNetworkFields::from_material(&net, &material, &query, || true),
        Err(FluidDuctError::Duct(DuctError::Cancelled { .. }))
    ));
    assert!(matches!(
        PreparedNetworkFields::from_material_with_loss(
            &net,
            &material,
            ConeSlices::Auto,
            &query,
            || true
        ),
        Err(FluidDuctError::Duct(DuctError::Cancelled { .. }))
    ));
    assert!(matches!(
        PreparedNetworkFieldSpectrum::from_material(&net, &material, &spectrum, || true),
        Err(FluidDuctError::Duct(DuctError::Cancelled { .. }))
    ));
    assert!(matches!(
        PreparedNetworkFieldSpectrum::from_material_with_loss(
            &net,
            &material,
            ConeSlices::Auto,
            &spectrum,
            || true
        ),
        Err(FluidDuctError::Duct(DuctError::Cancelled { .. }))
    ));
    assert_eq!((material.acoustic.get(), material.complete.get()), (0, 0));
    let prepared = PreparedNetworkFields::from_material(&net, &material, &query, || false).unwrap();
    let amplitudes = prepared.port_response().default_excitations();
    let accepted = prepared.sample(amplitudes, || false).unwrap();
    assert!(matches!(
        prepared.sample(amplitudes, || true),
        Err(DuctError::Cancelled { .. })
    ));
    assert_eq!(prepared.sample(amplitudes, || false).unwrap(), accepted);
    assert_eq!((material.acoustic.get(), material.complete.get()), (1, 0));
}
