use super::*;
use crate::network::field::*;
use fs_fft::harmonic_window::HarmonicWindow;

mod spectrum;
mod welch;

fn pressure_source(label: &str, value: C64) -> PortSource {
    PortSource {
        label: label.into(),
        excitation: PortExcitation::Pressure {
            open_circuit_pressure_pa: value,
            impedance_pa_s_m3: C64::ZERO,
        },
    }
}

fn cylinder_net(length: f64) -> DuctNetwork {
    DuctNetwork {
        nodes: vec![port("left"), port("right")],
        edges: vec![edge("tube", 0, 1, 0.008, length)],
    }
}

fn probes(length: f64, fractions: &[f64]) -> Vec<NetworkFieldProbe> {
    fractions
        .iter()
        .map(|fraction| NetworkFieldProbe {
            edge_id: "tube".into(),
            section: NetworkFieldSection::Cylinder {
                segment: 0,
                offset_m: length * fraction,
            },
        })
        .collect()
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
        source_clock: "manufactured-common-clock",
        sources,
        source_bank: None,
        loads,
        coupled_load: None,
        probes,
        budget: NetworkFieldBudget::default(),
    }
}

fn assert_close(a: C64, b: C64, scale: f64) {
    assert!((a - b).abs() <= 4.0e-10 * scale, "{a:?} != {b:?}");
}

fn phase(k: C64, distance: f64) -> C64 {
    let arg = k.re * distance;
    C64::new(arg.cos(), arg.sin()).scale((-k.im * distance).exp())
}

fn analytic(k: C64, z: C64, length: f64, x: f64, forward: C64, backward: C64) -> [C64; 2] {
    let a = forward * phase(k, x);
    let b = backward * phase(k, length - x);
    [a + b, (a - b) / z]
}

fn cone_net(r0: f64, r1: f64, length: f64) -> DuctNetwork {
    let mut net = cylinder_net(length);
    net.edges[0].duct.segments[0] = Segment::Cone {
        inlet_radius: r0,
        outlet_radius: r1,
        length,
    };
    net
}

fn cone_points(length: f64, fractions: &[f64]) -> Vec<NetworkFieldProbe> {
    cone_points_at(fractions.iter().map(|fraction| length * fraction))
}

fn cone_points_at(offsets: impl IntoIterator<Item = f64>) -> Vec<NetworkFieldProbe> {
    offsets
        .into_iter()
        .map(|offset_m| NetworkFieldProbe {
            edge_id: "tube".into(),
            section: NetworkFieldSection::Cone {
                segment: 0,
                offset_m,
            },
        })
        .collect()
}

fn assert_conical_field(
    sample: &NetworkFieldSample,
    expected: [C64; 2],
    area: f64,
    specific_impedance: f64,
) {
    let [p, u] = expected;
    // Scale with BOTH physical fields, so a pressure or velocity null is
    // checked against the nonzero wave field rather than a large floor.
    let pressure_scale = p.abs().max(u.abs() * specific_impedance / area);
    let flow_scale = pressure_scale * area / specific_impedance;
    for (actual, expected, scale) in [
        (sample.pressure_pa, p, pressure_scale),
        (sample.volume_flow_m3_s, u, flow_scale),
        (
            sample.particle_velocity_m_s,
            u.scale(1.0 / area),
            flow_scale / area,
        ),
    ] {
        assert!(
            (actual - expected).abs() <= 4e-11 * scale,
            "{actual:?} != {expected:?}; scale {scale:e}"
        );
    }
    let power = 0.5 * (p * u.conj()).re;
    assert!((sample.mean_power_w - power).abs() <= 4e-11 * pressure_scale * flow_scale);
}

// Independent solution of (r p)'' + k² (r p) = 0 and
// U = A p'/(i k Z_specific). The curvature term is essential:
// pressure and flow nodes do not have the cylinder's phase relationship.
fn cone_analytic(
    k: C64,
    z: C64,
    radii: [f64; 2],
    length: f64,
    x: f64,
    a: C64,
    b: C64,
) -> [C64; 2] {
    let [r0, r1] = radii;
    let radius = r0 + (r1 - r0) * (x / length);
    let forward = a * phase(k, x).scale(r0 / radius);
    let backward = b * phase(k, length - x).scale(r1 / radius);
    let p = forward + backward;
    let curvature = C64::new(0.0, (r1 - r0) / length / radius) / k;
    let u = (forward - backward + p * curvature)
        .scale(core::f64::consts::PI * radius * radius)
        / z;
    [p, u]
}

#[test]
fn conical_interference_recovers_spherical_fields_and_pressure_nulls() {
    let gas = gas();
    let length = 0.37;
    let omega = 2971.0;
    let k = C64::from_re(omega / gas.sound_speed);
    let z = C64::from_re(gas.density * gas.sound_speed);
    let node = length * 0.41;
    let a = C64::new(1.4, -0.3);
    for radii in [[0.006, 0.015], [0.015, 0.006]] {
        let b = -a * phase(k, 2.0 * node - length).scale(radii[0] / radii[1]);
        let field = |x| cone_analytic(k, z, radii, length, x, a, b);
        let sources = [
            pressure_source("right", field(length)[0]),
            pressure_source("left", field(0.0)[0]),
        ];
        let net = cone_net(radii[0], radii[1], length);
        let points = cone_points(length, &[0.0, 0.13, 0.41, 0.63, 0.91, 1.0]);
        let owner = PreparedNetworkFields::new(
            &net,
            &gas,
            &request(omega, LossModel::Lossless, &sources, &[], &points),
            || false,
        )
        .unwrap();
        let actual = owner
            .sample(owner.port_response().default_excitations(), || false)
            .unwrap();
        for (sample, location) in actual.iter().zip(owner.locations()) {
            assert_conical_field(
                sample,
                field(location.axial_position_m),
                location.area_m2,
                z.abs(),
            );
            assert!(sample.mean_power_w.abs() < 2e-15);
        }
        assert!(actual[2].pressure_pa.abs() < 2e-13);
        assert!(actual[2].volume_flow_m3_s.abs() > 1e-7);
        // Endpoints retain exactly the original native cut and area arithmetic.
        let endpoints = [
            NetworkFieldProbe {
                edge_id: "tube".into(),
                section: NetworkFieldSection::Inlet { segment: 0 },
            },
            NetworkFieldProbe {
                edge_id: "tube".into(),
                section: NetworkFieldSection::Outlet { segment: 0 },
            },
        ];
        let endpoint_owner = PreparedNetworkFields::new(
            &net,
            &gas,
            &request(omega, LossModel::Lossless, &sources, &[], &endpoints),
            || false,
        )
        .unwrap();
        let endpoint_values = endpoint_owner
            .sample(endpoint_owner.port_response().default_excitations(), || false)
            .unwrap();
        assert_eq!([actual[0], actual[5]].as_slice(), endpoint_values.as_slice());
    }
}

// Back-propagate a finite state with the direct spherical ODE solution.
// Production reconstruction instead uses decaying scattering from both ends.
fn cone_backwards(
    wave: &crate::SegmentWave,
    radii: [f64; 2],
    length: f64,
    x: f64,
    end: [C64; 2],
) -> [C64; 2] {
    let [r0, r1] = radii;
    let area = core::f64::consts::PI * r1 * r1;
    let curvature = C64::new(0.0, (r1 - r0) / length / r1) / wave.wavenumber;
    let difference = end[1] * wave.specific_impedance.scale(1.0 / area) - end[0] * curvature;
    let forward = (end[0] + difference).scale(0.5);
    let backward = (end[0] - difference).scale(0.5);
    let a = forward * phase(wave.wavenumber, -length).scale(r1 / r0);
    cone_analytic(
        wave.wavenumber,
        wave.specific_impedance,
        radii,
        length,
        x,
        a,
        backward,
    )
}

#[test]
fn lossy_conical_interference_keeps_original_midpoint_slices() {
    let gas = gas();
    let (length, omega) = (0.18, 1700.0);
    for radii in [[0.001, 0.004], [0.004, 0.001]] {
        let count = crate::cone_loss_slices(radii[0], radii[1]).unwrap();
        let dx = length / count as f64;
        let radius = |i| radii[0] + (radii[1] - radii[0]) * (i as f64 / count as f64);
        let waves: Vec<_> = (0..count)
            .map(|i| {
                crate::segment_wave(
                    &gas,
                    f64::midpoint(radius(i), radius(i + 1)),
                    omega,
                    LossModel::Bessel,
                )
                .unwrap()
            })
            .collect();
        let mut states = vec![[C64::ZERO; 2]; count + 1];
        states[count] = [C64::new(0.8, -0.2), C64::new(4e-9, 3e-9)];
        for i in (0..count).rev() {
            states[i] = cone_backwards(
                &waves[i],
                [radius(i), radius(i + 1)],
                dx,
                0.0,
                states[i + 1],
            );
        }
        let mut offsets = Vec::new();
        let mut expected = Vec::new();
        for i in 0..count {
            let x0 = length * (i as f64 / count as f64);
            let x1 = length * ((i + 1) as f64 / count as f64);
            let mut local_points = vec![x0, x0 + 0.37 * (x1 - x0), x1.next_down()];
            if i > 0 {
                local_points.push(x0.next_up());
            }
            for x in local_points {
                offsets.push(x);
                expected.push(cone_backwards(
                    &waves[i],
                    [radius(i), radius(i + 1)],
                    dx,
                    dx * ((x - x0) / (x1 - x0)),
                    states[i + 1],
                ));
            }
        }
        offsets.push(length);
        expected.push(states[count]);
        let points = cone_points_at(offsets);
        let sources = [
            pressure_source("left", states[0][0]),
            pressure_source("right", states[count][0]),
        ];
        let net = cone_net(radii[0], radii[1], length);
        let owner = PreparedNetworkFields::new(
            &net,
            &gas,
            &request(omega, LossModel::Bessel, &sources, &[], &points),
            || false,
        )
        .unwrap();
        let actual = owner
            .sample(owner.port_response().default_excitations(), || false)
            .unwrap();
        for ((sample, location), field) in actual.iter().zip(owner.locations()).zip(expected) {
            assert_conical_field(
                sample,
                field,
                location.area_m2,
                gas.density * gas.sound_speed,
            );
        }
    }
}

#[test]
fn material_conical_fields_retain_explicit_loss_slice_resolution() {
    use fs_material::fluid_acoustics::FluidAcousticState;
    use fs_material::fluid_acoustics::thermoviscous::ThermoviscousFluidState;
    use fs_material::fluid_transport::FluidTransportState;

    let material = ThermoviscousFluidState::try_new(
        FluidAcousticState::try_new(950.0, 1300.0).unwrap(),
        FluidTransportState::try_new(950.0, 4200.0, 0.6, 0.2).unwrap(),
        4100.0,
    )
    .unwrap();
    let (length, omega) = (0.19, 1700.0);
    for radii in [[0.001, 0.004], [0.004, 0.001]] {
        for count in [1, 4, 11] {
            let dx = length / count as f64;
            let radius = |i| radii[0] + (radii[1] - radii[0]) * (i as f64 / count as f64);
            let waves: Vec<_> = (0..count)
                .map(|i| {
                    crate::fluid::segment_wave_with_loss(
                        &material,
                        f64::midpoint(radius(i), radius(i + 1)),
                        omega,
                        LossModel::Bessel,
                    )
                    .unwrap()
                })
                .collect();
            let mut states = vec![[C64::ZERO; 2]; count + 1];
            states[count] = [C64::new(0.8, -0.2), C64::new(4e-12, 3e-12)];
            for i in (0..count).rev() {
                states[i] = cone_backwards(
                    &waves[i],
                    [radius(i), radius(i + 1)],
                    dx,
                    0.0,
                    states[i + 1],
                );
            }
            let fractions: Vec<_> = (0..count)
                .flat_map(|i| [i as f64 / count as f64, (i as f64 + 0.37) / count as f64])
                .chain(core::iter::once(1.0))
                .collect();
            let points = cone_points(length, &fractions);
            let sources = [
                pressure_source("left", states[0][0]),
                pressure_source("right", states[count][0]),
            ];
            let net = cone_net(radii[0], radii[1], length);
            let owner = PreparedNetworkFields::from_material_with_loss(
                &net,
                &material,
                crate::ConeSlices::Fixed(count),
                &request(omega, LossModel::Bessel, &sources, &[], &points),
                || false,
            )
            .unwrap();
            let actual = owner
                .sample(owner.port_response().default_excitations(), || false)
                .unwrap();
            for i in 0..count {
                for (j, local) in [(2 * i, 0.0), (2 * i + 1, 0.37)] {
                    let field = cone_backwards(
                        &waves[i],
                        [radius(i), radius(i + 1)],
                        dx,
                        local * dx,
                        states[i + 1],
                    );
                    let location = &owner.locations()[j];
                    assert_eq!(location.density_kg_m3, 950.0);
                    assert_eq!(location.sound_speed_m_s, 1300.0);
                    assert_conical_field(&actual[j], field, location.area_m2, 950.0 * 1300.0);
                }
            }
            assert_conical_field(
                actual.last().unwrap(),
                states[count],
                owner.locations().last().unwrap().area_m2,
                950.0 * 1300.0,
            );
            assert!(actual[0].mean_power_w > actual.last().unwrap().mean_power_w);
        }
    }
}

#[test]
fn conical_uniform_limit_matches_original_cylinder_at_slice_neighbors() {
    let gas = gas();
    let length = 0.37_f64;
    let cut = length * 0.5;
    let offsets = [
        0.0,
        length * f64::EPSILON,
        0.13 * length,
        cut.next_down(),
        cut,
        cut.next_up(),
        0.81 * length,
        length.next_down(),
        length,
    ];
    let conical_points = cone_points_at(offsets);
    let cylindrical_points: Vec<_> = offsets
        .into_iter()
        .map(|offset_m| NetworkFieldProbe {
            edge_id: "tube".into(),
            section: NetworkFieldSection::Cylinder {
                segment: 0,
                offset_m,
            },
        })
        .collect();
    let sources = [
        pressure_source("left", C64::new(1.1, -0.2)),
        pressure_source("right", C64::new(-0.3, 0.7)),
    ];
    for loss in [LossModel::Lossless, LossModel::Bessel] {
        let conical = PreparedNetworkFields::new(
            &cone_net(0.008, 0.008, length),
            &gas,
            &request(1700.0, loss, &sources, &[], &conical_points),
            || false,
        )
        .unwrap();
        let cylindrical = PreparedNetworkFields::new(
            &cylinder_net(length),
            &gas,
            &request(1700.0, loss, &sources, &[], &cylindrical_points),
            || false,
        )
        .unwrap();
        // Each source column and their coherent sum must agree; one drive
        // alone can hide a source-dependent reconstruction error.
        for excitation in [
            [C64::ONE, C64::ZERO],
            [C64::ZERO, C64::ONE],
            [C64::new(1.1, -0.2), C64::new(-0.3, 0.7)],
        ] {
            let actual = conical.sample(&excitation, || false).unwrap();
            let expected = cylindrical.sample(&excitation, || false).unwrap();
            for ((sample, other), location) in actual.iter().zip(expected).zip(conical.locations()) {
                assert_conical_field(
                    sample,
                    [other.pressure_pa, other.volume_flow_m3_s],
                    location.area_m2,
                    gas.density * gas.sound_speed,
                );
            }
        }
    }
}

#[test]
fn conical_fields_refuse_offsets_and_mismatched_authored_geometry() {
    let gas = gas();
    let length = 0.37_f64;
    let net = cone_net(0.006, 0.015, length);
    let sources = [
        pressure_source("left", C64::ONE),
        pressure_source("right", C64::ZERO),
    ];
    for offset in [-f64::EPSILON, length.next_up(), f64::NAN, f64::INFINITY] {
        let points = cone_points_at([offset]);
        let query = request(1700.0, LossModel::Bessel, &sources, &[], &points);
        assert!(PreparedNetworkFields::new(&net, &gas, &query, || false).is_err());
    }
    let points = cone_points_at([length * 0.5]);
    let query = request(1700.0, LossModel::Lossless, &sources, &[], &points);
    assert!(PreparedNetworkFields::new(&cylinder_net(length), &gas, &query, || false).is_err());
    let points = [NetworkFieldProbe {
        edge_id: "tube".into(),
        section: NetworkFieldSection::Cone {
            segment: 1,
            offset_m: 0.0,
        },
    }];
    let query = request(1700.0, LossModel::Lossless, &sources, &[], &points);
    assert!(PreparedNetworkFields::new(&net, &gas, &query, || false).is_err());
    let points = probes(length, &[0.5]);
    let query = request(1700.0, LossModel::Lossless, &sources, &[], &points);
    assert!(PreparedNetworkFields::new(&net, &gas, &query, || false).is_err());
}

// G1: prescribed endpoint pressures are manufactured from independently
// evaluated counterpropagating waves, including lossy complex impedance.
#[test]
fn spatial_fields_recover_both_physical_waves_and_directed_flow() {
    let gas = gas();
    let length = 0.37;
    let net = cylinder_net(length);
    let omega = 2971.0;
    let locations = probes(length, &[0.0, 0.07, 0.3, 0.5, 0.9, 1.0]);
    let (a, b) = (C64::new(1.4, -0.3), C64::new(-0.4, 0.9));
    for loss in [LossModel::Lossless, LossModel::Bessel] {
        let wave = crate::segment_wave(&gas, 0.008, omega, loss).unwrap();
        let field = |x| {
            analytic(
                wave.wavenumber,
                wave.characteristic_impedance,
                length,
                x,
                a,
                b,
            )
        };
        let sources = [
            pressure_source("right", field(length)[0]),
            pressure_source("left", field(0.0)[0]),
        ];
        let query = request(omega, loss, &sources, &[], &locations);
        let prepared = PreparedNetworkFields::new(&net, &gas, &query, || false).unwrap();
        let output = prepared
            .sample(prepared.port_response().default_excitations(), || false)
            .unwrap();
        let z = wave.characteristic_impedance.abs();
        for (sample, location) in output.iter().zip(prepared.locations()) {
            let [p, q] = field(location.axial_position_m);
            assert_close(sample.pressure_pa, p, 2.0);
            assert_close(sample.volume_flow_m3_s, q, 2.0 / z);
            assert_close(
                sample.particle_velocity_m_s,
                q.scale(1.0 / location.area_m2),
                2.0 / z / location.area_m2,
            );
            let expected = 0.5 * (p * q.conj()).re;
            assert!((sample.mean_power_w - expected).abs() < 5.0e-10 / z);
        }
        // Exterior flow convention is into the graph, while every spatial
        // probe's U follows the original from->to edge orientation.
        let mut boundary = prepared.port_response().response_buffer();
        prepared
            .port_response()
            .drive_into(
                prepared.port_response().default_excitations(),
                &mut boundary,
            )
            .unwrap();
        assert_close(
            output[0].volume_flow_m3_s,
            boundary.fields.port_inflows[0],
            2.0 / z,
        );
        assert_close(
            output[5].volume_flow_m3_s,
            -boundary.fields.port_inflows[1],
            2.0 / z,
        );
    }
}

// G1: the loaded solution and all interior fields remain finite when the
// open-circuit impedance fails at the exact half-wave pole.
#[test]
fn spatial_half_wave_uses_the_actual_loaded_graph() {
    let gas = gas();
    let length = 0.2;
    let net = cylinder_net(length);
    let omega = core::f64::consts::PI * gas.sound_speed / length;
    assert!(matches!(
        network_port_impedance(&net, &gas, omega, LossModel::Lossless),
        Err(DuctError::Singular)
    ));
    let z = gas.density * gas.sound_speed / (core::f64::consts::PI * 0.008_f64.powi(2));
    let q = C64::new(1.0e-6, -0.3e-6);
    let sources = [PortSource {
        label: "left".into(),
        excitation: PortExcitation::VolumeFlow(q),
    }];
    let loads = [PortTermination {
        label: "right".into(),
        load: NetworkLoad::Impedance(C64::from_re(z)),
    }];
    let points = probes(length, &[0.0, 0.25, 0.5, 0.75, 1.0]);
    let query = request(omega, LossModel::Lossless, &sources, &loads, &points);
    let prepared = PreparedNetworkFields::new(&net, &gas, &query, || false).unwrap();
    let output = prepared.sample(&[q], || false).unwrap();
    for (sample, location) in output.iter().zip(prepared.locations()) {
        let want = q * phase(
            C64::from_re(omega / gas.sound_speed),
            location.axial_position_m,
        );
        assert_close(sample.volume_flow_m3_s, want, q.abs());
        assert_close(sample.pressure_pa, want.scale(z), z * q.abs());
        assert!((sample.mean_power_w - 0.5 * z * q.norm_sq()).abs() < 2.0e-14);
    }
}

// G1: exact pressure and flow nodes are complementary. Dividing by a local
// impedance/admittance or adding independently computed source powers loses it.
#[test]
fn spatial_standing_wave_retains_complementary_interference_nodes() {
    let gas = gas();
    let length = 0.31;
    let omega = 1741.0;
    let net = cylinder_net(length);
    let wave = crate::segment_wave(&gas, 0.008, omega, LossModel::Lossless).unwrap();
    let points = probes(length, &[0.5]);
    for sign in [1.0, -1.0] {
        let field = |x| {
            analytic(
                wave.wavenumber,
                wave.characteristic_impedance,
                length,
                x,
                C64::ONE,
                C64::from_re(sign),
            )
        };
        let sources = [
            pressure_source("left", field(0.0)[0]),
            pressure_source("right", field(length)[0]),
        ];
        let query = request(omega, LossModel::Lossless, &sources, &[], &points);
        let prepared = PreparedNetworkFields::new(&net, &gas, &query, || false).unwrap();
        let window = HarmonicWindow::new(-0.00017, 0.00031).unwrap();
        let mean = prepared.port_response().default_excitations();
        let moments = prepared
            .observe_window(window, Some(mean), &[], || false)
            .unwrap()[0];
        let z = wave.characteristic_impedance.abs();
        if sign > 0.0 {
            assert!(moments.mean_square_pressure_pa2 > 0.1);
            assert!(moments.mean_square_volume_flow_m6_s2 < 1.0e-20 / (z * z));
        } else {
            assert!(moments.mean_square_pressure_pa2 < 1.0e-20);
            assert!(moments.mean_square_volume_flow_m6_s2 > 0.1 / (z * z));
        }
        assert!(moments.mean_power_w.abs() < 1.0e-10 / z);
    }
}

// G1/G3: integrate the actual analytic real fields in time and average four
// unit proper latent states. Mean/factor cross terms vanish statistically,
// while every correlated source cross term within the factor survives.
#[test]
fn spatial_finite_window_ensemble_matches_independent_real_signal_integral() {
    let gas = gas();
    let length = 0.22;
    let omega = 731.0;
    let net = cylinder_net(length);
    let wave = crate::segment_wave(&gas, 0.008, omega, LossModel::Lossless).unwrap();
    let points = probes(length, &[0.17, 0.61]);
    let sources = [
        pressure_source("left", C64::ZERO),
        pressure_source("right", C64::ZERO),
    ];
    let query = request(omega, LossModel::Lossless, &sources, &[], &points);
    let prepared = PreparedNetworkFields::new(&net, &gas, &query, || false).unwrap();
    let mean_waves = [C64::new(0.8, -0.3), C64::new(-0.2, 0.1)];
    let factor_waves = [C64::new(-0.4, 0.7), C64::new(0.6, 0.2)];
    let endpoints = |waves: [C64; 2]| {
        [
            analytic(
                wave.wavenumber,
                wave.characteristic_impedance,
                length,
                0.0,
                waves[0],
                waves[1],
            )[0],
            analytic(
                wave.wavenumber,
                wave.characteristic_impedance,
                length,
                length,
                waves[0],
                waves[1],
            )[0],
        ]
    };
    let mean = endpoints(mean_waves);
    let factor = endpoints(factor_waves);
    let window = HarmonicWindow::new(-0.0013, 0.0047).unwrap();
    let actual = prepared
        .observe_window(window, Some(&mean), &[&factor], || false)
        .unwrap();
    let states = [C64::ONE, -C64::ONE, C64::new(0.0, 1.0), C64::new(0.0, -1.0)];
    let n = 20_000;
    for (index, location) in prepared.locations().iter().enumerate() {
        let m = analytic(
            wave.wavenumber,
            wave.characteristic_impedance,
            length,
            location.axial_position_m,
            mean_waves[0],
            mean_waves[1],
        );
        let f = analytic(
            wave.wavenumber,
            wave.characteristic_impedance,
            length,
            location.axial_position_m,
            factor_waves[0],
            factor_waves[1],
        );
        let mut sum = [0.0; 5];
        for j in 0..n {
            let time = window.start_s() + (j as f64 + 0.5) * window.duration_s() / n as f64;
            let rotation = C64::new((omega * time).cos(), -(omega * time).sin());
            for state in states {
                let p = ((m[0] + f[0] * state) * rotation).re;
                let q = ((m[1] + f[1] * state) * rotation).re;
                for (slot, value) in sum.iter_mut().zip([p, q, p * p, q * q, p * q]) {
                    *slot += value;
                }
            }
        }
        let expected = sum.map(|v| v / (4 * n) as f64);
        let got = actual[index];
        for (a, b, scale) in [
            (got.mean_pressure_pa, expected[0], 1.0),
            (
                got.mean_volume_flow_m3_s,
                expected[1],
                1.0 / wave.characteristic_impedance.abs(),
            ),
            (got.mean_square_pressure_pa2, expected[2], 2.0),
            (
                got.mean_square_volume_flow_m6_s2,
                expected[3],
                2.0 / wave.characteristic_impedance.norm_sq(),
            ),
            (
                got.mean_power_w,
                expected[4],
                2.0 / wave.characteristic_impedance.abs(),
            ),
        ] {
            assert!((a - b).abs() < 2.0e-8 * scale, "{a:e} != {b:e}");
        }
        let independent = prepared
            .observe_window(
                window,
                None,
                &[&[factor[0], C64::ZERO], &[C64::ZERO, factor[1]]],
                || false,
            )
            .unwrap();
        let coherent = prepared
            .observe_window(window, None, &[&factor], || false)
            .unwrap();
        assert!(
            (independent[index].mean_square_pressure_pa2
                - coherent[index].mean_square_pressure_pa2)
                .abs()
                > 0.01
        );
    }
}

// G1: exponentials always decay from their own incident end. An opaque line
// has two finite boundary layers and a dark interior, not an ABCD overflow.
#[test]
fn spatial_opaque_lossy_line_keeps_both_incident_boundary_layers() {
    let gas = gas();
    let length = 1.0e6;
    let omega = 3001.0;
    let net = cylinder_net(length);
    let wave = crate::segment_wave(&gas, 0.008, omega, LossModel::Bessel).unwrap();
    assert!(wave.wavenumber.im * length > 1000.0);
    let a = C64::new(1.0, 0.3);
    let b = C64::new(0.2, -0.8);
    let sources = [pressure_source("left", a), pressure_source("right", b)];
    let points = probes(length, &[0.0, 1.0e-6, 0.5, 1.0 - 1.0e-6, 1.0]);
    let query = request(omega, LossModel::Bessel, &sources, &[], &points);
    let prepared = PreparedNetworkFields::new(&net, &gas, &query, || false).unwrap();
    let actual = prepared
        .sample(prepared.port_response().default_excitations(), || false)
        .unwrap();
    for (sample, location) in actual.iter().zip(prepared.locations()) {
        let [p, q] = analytic(
            wave.wavenumber,
            wave.characteristic_impedance,
            length,
            location.axial_position_m,
            a,
            b,
        );
        assert_close(sample.pressure_pa, p, 2.0);
        assert_close(
            sample.volume_flow_m3_s,
            q,
            2.0 / wave.characteristic_impedance.abs(),
        );
    }
    assert_eq!(actual[2].pressure_pa, C64::ZERO);
    assert_eq!(actual[2].volume_flow_m3_s, C64::ZERO);
}

// G1/G3: zero-length side branches have separate physical states on each side,
// with original correction policy and unchanged position.
#[test]
fn spatial_compact_hole_keeps_both_sides_and_removed_flow() {
    let gas = gas();
    let omega = 1301.0;
    let mut net = cylinder_net(0.1);
    net.edges[0].duct.segments.extend([
        Segment::ToneHole {
            hole_radius: 0.002,
            chimney_height: 0.003,
            bore_radius: 0.008,
            state: HoleState::Open,
        },
        Segment::Cylinder {
            radius: 0.008,
            length: 0.12,
        },
    ]);
    let sources = [PortSource {
        label: "left".into(),
        excitation: PortExcitation::VolumeFlow(C64::from_re(1.0e-6)),
    }];
    let z = gas.density * gas.sound_speed / (core::f64::consts::PI * 0.008_f64.powi(2));
    let loads = [PortTermination {
        label: "right".into(),
        load: NetworkLoad::Impedance(C64::from_re(z)),
    }];
    let points = [
        NetworkFieldProbe {
            edge_id: "tube".into(),
            section: NetworkFieldSection::Inlet { segment: 1 },
        },
        NetworkFieldProbe {
            edge_id: "tube".into(),
            section: NetworkFieldSection::Outlet { segment: 1 },
        },
    ];
    let query = request(omega, LossModel::Lossless, &sources, &loads, &points);
    let prepared = PreparedNetworkFields::new(&net, &gas, &query, || false).unwrap();
    let actual = prepared.sample(&[C64::from_re(1.0e-6)], || false).unwrap();
    assert_eq!(
        prepared.locations()[0].axial_position_m,
        prepared.locations()[1].axial_position_m
    );
    assert_close(actual[0].pressure_pa, actual[1].pressure_pa, z * 1.0e-6);
    assert!((actual[0].volume_flow_m3_s - actual[1].volume_flow_m3_s).abs() > 1.0e-9);
    assert!(actual[0].mean_power_w > actual[1].mean_power_w);
    let mut ports = prepared.port_response().response_buffer();
    prepared
        .port_response()
        .drive_into(&[C64::from_re(1.0e-6)], &mut ports)
        .unwrap();
    assert!(
        (actual[0].mean_power_w - actual[1].mean_power_w - ports.fields.network_power_w).abs()
            < 1.0e-14
    );
}

// G3: node/edge/source list order is not the phase or summation order.
#[test]
fn spatial_reindexing_preserves_source_coordinates_and_fields() {
    let gas = gas();
    let mut a = DuctNetwork {
        nodes: vec![
            port("left"),
            port("right"),
            node("join", NodeBoundary::Junction),
        ],
        edges: vec![
            edge("tube", 0, 2, 0.008, 0.12),
            edge("other", 2, 1, 0.008, 0.19),
        ],
    };
    let points = probes(0.12, &[0.4]);
    let mut sources = [
        pressure_source("right", C64::new(-0.4, 0.2)),
        pressure_source("left", C64::new(0.7, -0.1)),
    ];
    let query = request(1723.0, LossModel::Bessel, &sources, &[], &points);
    let first = PreparedNetworkFields::new(&a, &gas, &query, || false).unwrap();
    a.nodes.swap(0, 1);
    for edge in &mut a.edges {
        if edge.from < 2 {
            edge.from = 1 - edge.from;
        }
        if edge.to < 2 {
            edge.to = 1 - edge.to;
        }
    }
    a.edges.reverse();
    sources.reverse();
    let query = request(1723.0, LossModel::Bessel, &sources, &[], &points);
    let second = PreparedNetworkFields::new(&a, &gas, &query, || false).unwrap();
    assert_eq!(
        first.port_response().source_coordinates(),
        second.port_response().source_coordinates()
    );
    assert_eq!(first.pressure_transfer(), second.pressure_transfer());
    assert_eq!(first.flow_transfer(), second.flow_transfer());
}

// G0/G4: a probe must describe real authored geometry; refusing a point,
// budget or cancellation publishes no partially usable prepared bank.
#[test]
fn spatial_fields_refuse_wrong_geometry_budget_and_cancelled_queries() {
    let gas = gas();
    let net = cylinder_net(0.2);
    let sources = [
        pressure_source("left", C64::ONE),
        pressure_source("right", C64::ZERO),
    ];
    let mut points = probes(0.2, &[0.5]);
    let query = request(1401.0, LossModel::Lossless, &sources, &[], &points);
    assert!(matches!(
        PreparedNetworkFields::new(&net, &gas, &query, || true),
        Err(DuctError::Cancelled { .. })
    ));
    let prepared = PreparedNetworkFields::new(&net, &gas, &query, || false).unwrap();
    assert!(prepared.sample(&[C64::ONE], || false).is_err());
    assert!(prepared.sample(&[C64::ONE, C64::ZERO], || true).is_err());
    let window = HarmonicWindow::new(0.0, 0.01).unwrap();
    let factor = [C64::ONE, C64::ZERO];
    let factors = vec![factor.as_slice(); MAX_FIELD_FACTORS + 1];
    assert!(
        prepared
            .observe_window(window, None, &factors, || false)
            .is_err()
    );
    let mut short = request(1401.0, LossModel::Lossless, &sources, &[], &points);
    short.budget.maximum_transfer_entries = 11; // (2*one probe + 4*one edge)*2 sources = 12
    assert!(PreparedNetworkFields::new(&net, &gas, &short, || false).is_err());
    points[0].section = NetworkFieldSection::Cylinder {
        segment: 0,
        offset_m: 0.3,
    };
    let query = request(1401.0, LossModel::Lossless, &sources, &[], &points);
    assert!(PreparedNetworkFields::new(&net, &gas, &query, || false).is_err());
}
