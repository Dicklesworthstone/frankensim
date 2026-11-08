//! Native smooth-profile cuts preserve the original field solve and interference.

use super::*;
use crate::{ConeSlices, RadiusProfile, RiccatiSteps};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_material::fluid_acoustics::thermoviscous::ThermoviscousFluidState;
use fs_material::fluid_transport::FluidTransportState;
use std::sync::Arc;

fn profile_net(profile: &RadiusProfile, steps: usize) -> DuctNetwork {
    let mut net = cylinder_net(profile.length());
    net.edges[0].duct.segments = vec![Segment::Profile {
        profile: Arc::new(profile.clone()),
        steps: RiccatiSteps(steps),
    }];
    net
}

fn panel_probes(segment: usize, boundaries: &[usize]) -> Vec<NetworkFieldProbe> {
    boundaries
        .iter()
        .map(|&boundary| NetworkFieldProbe {
            edge_id: "tube".into(),
            section: NetworkFieldSection::Profile { segment, boundary },
        })
        .collect()
}

fn material() -> ThermoviscousFluidState {
    ThermoviscousFluidState::try_new(
        FluidAcousticState::try_new(950.0, 1300.0).unwrap(),
        FluidTransportState::try_new(950.0, 4200.0, 0.6, 0.2).unwrap(),
        4100.0,
    )
    .unwrap()
}

// Independent exact solution of p'' + flare*p' + k²*p = 0,
// U = A*p'/(i*k*rho*c), with forward-left and backward-right amplitudes.
fn exponential_field(
    gas: &GasState,
    omega: f64,
    flare: f64,
    length: f64,
    x: f64,
    amplitudes: [C64; 2],
) -> [C64; 2] {
    let k = omega / gas.sound_speed;
    let q = (k * k - 0.25 * flare * flare).sqrt();
    let forward = amplitudes[0]
        * phase(C64::from_re(q), x).scale((-0.5 * flare * x).exp());
    let backward = amplitudes[1]
        * phase(C64::from_re(q), length - x).scale((0.5 * flare * (length - x)).exp());
    let p = forward + backward;
    let area = core::f64::consts::PI * 0.008_f64.powi(2) * (flare * x).exp();
    let u = ((forward - backward).scale(q / k) + p * C64::new(0.0, flare / (2.0 * k)))
        .scale(area / (gas.density * gas.sound_speed));
    [p, u]
}

// G1: the ORIGINAL Magnus resolution converges to a manufactured standing
// wave, including a pressure null with nonzero volume flow. Observation
// points are original panel cuts; adding them must not refine the physics.
#[test]
fn exponential_profile_interference_converges_at_original_panel_boundaries() {
    let gas = gas();
    let (length, omega) = (0.37, 1700.0);
    let node = 3.0 * length / 8.0;
    let a = C64::new(1.4, -0.3);
    for flare in [-2.4, 2.4] {
        let profile = RadiusProfile::exponential(0.008, flare, length).unwrap();
        let k = omega / gas.sound_speed;
        let q = (k * k - 0.25 * flare * flare).sqrt();
        let b = -a * phase(C64::from_re(q), 2.0 * node - length)
            .scale((-0.5 * flare * length).exp());
        let field = |x| exponential_field(&gas, omega, flare, length, x, [a, b]);
        let sources = [
            pressure_source("left", field(0.0)[0]),
            pressure_source("right", field(length)[0]),
        ];
        let mut errors = Vec::new();
        for steps in [8, 16, 32] {
            let net = profile_net(&profile, steps);
            let points = panel_probes(
                0,
                &[0, steps / 4, 3 * steps / 8, steps / 2, 3 * steps / 4, steps],
            );
            let owner = PreparedNetworkFields::new(
                &net,
                &gas,
                &request(omega, LossModel::Lossless, &sources, &[], &points),
                || false,
            )
            .unwrap();
            let drive = owner.port_response().default_excitations();
            let samples = owner.sample(drive, || false).unwrap();
            let mut error = 0.0_f64;
            for (sample, location) in samples.iter().zip(owner.locations()) {
                let [p, u] = field(location.axial_position_m);
                let z = gas.density * gas.sound_speed / location.area_m2;
                error = error.max((sample.pressure_pa - p).abs())
                    .max((sample.volume_flow_m3_s - u).abs() * z);
                assert!(sample.mean_power_w.abs() < 2e-16);
            }
            errors.push(error);
            if steps == 32 {
                assert!(samples[2].pressure_pa.abs() < 4e-8);
                assert!(samples[2].volume_flow_m3_s.abs() > 1e-7);
                // The finite-window p*U transport includes the reactive
                // standing wave; its cycle-mean power alone is insufficient.
                let window = HarmonicWindow::new(0.00017, 0.00083).unwrap();
                let moments = owner
                    .observe_window(window, Some(drive), &[], || false)
                    .unwrap();
                let center = window.start_s() + 0.5 * window.duration_s();
                let carrier = phase(C64::from_re(omega), -center);
                let t = omega * window.duration_s();
                let mean_gain = (0.5 * t).sin() / (0.5 * t);
                let pair_gain = t.sin() / t;
                for (point, location) in moments.iter().zip(owner.locations()) {
                    let [p, u] = field(location.axial_position_m);
                    let inverse_z = location.area_m2 / (gas.density * gas.sound_speed);
                    let pair = |a: C64, b: C64| {
                        0.5 * (a * b.conj()).re
                            + 0.5 * pair_gain * (a * b * carrier * carrier).re
                    };
                    for (actual, expected, scale) in [
                        (point.mean_pressure_pa, (p * carrier).re * mean_gain, 1.0),
                        (
                            point.mean_volume_flow_m3_s,
                            (u * carrier).re * mean_gain,
                            inverse_z,
                        ),
                        (point.mean_square_pressure_pa2, pair(p, p), 4.0),
                        (
                            point.mean_square_volume_flow_m6_s2,
                            pair(u, u),
                            4.0 * inverse_z * inverse_z,
                        ),
                        (point.mean_power_w, pair(p, u), 4.0 * inverse_z),
                    ] {
                        assert!((actual - expected).abs() <= 1e-6 * scale);
                    }
                }
                assert!(moments[1].mean_power_w.abs() > 1e-9);
            }
        }
        for adjacent in errors.windows(2) {
            assert!(adjacent[0] > 12.0 * adjacent[1], "{errors:?}");
            assert!(adjacent[0] < 20.0 * adjacent[1], "{errors:?}");
        }
        assert!(errors[2] < 3e-7, "{errors:?}");
    }
}

// G1/G3: the cylindrical profile limit with original liquid diffusion has
// an independent two-wave solution. Include a preceding segment and a
// nonzero authored spline origin to test actual axial coordinates.
#[test]
fn material_profile_cuts_keep_complex_interference_and_exact_endpoints() {
    let material = material();
    let (radius, prefix, omega) = (0.002, 0.07, 1700.0);
    let profile = RadiusProfile::spline(
        vec![(4.0, radius), (4.07, radius), (4.2, radius)],
        None,
    )
    .unwrap();
    let steps = 7;
    let length = prefix + profile.length();
    let mut net = profile_net(&profile, steps);
    net.edges[0].duct.segments.insert(
        0,
        Segment::Cylinder {
            radius,
            length: prefix,
        },
    );
    let indices = [0, 1, 3, 6, steps];
    let points = panel_probes(1, &indices);
    let wave = crate::fluid::segment_wave_with_loss(&material, radius, omega, LossModel::Bessel)
        .unwrap();
    let field = |x| {
        analytic(
            wave.wavenumber,
            wave.characteristic_impedance,
            length,
            x,
            C64::new(1.2, -0.3),
            C64::new(-0.4, 0.7),
        )
    };
    let sources = [
        pressure_source("left", field(0.0)[0]),
        pressure_source("right", field(length)[0]),
    ];
    let owner = PreparedNetworkFields::from_material_with_loss(
        &net, &material, ConeSlices::Auto,
        &request(omega, LossModel::Bessel, &sources, &[], &points), || false,
    )
    .unwrap();
    let samples = owner
        .sample(owner.port_response().default_excitations(), || false)
        .unwrap();
    for ((sample, location), boundary) in samples.iter().zip(owner.locations()).zip(indices) {
        let offset = if boundary == 0 {
            0.0
        } else if boundary == steps {
            profile.length()
        } else {
            profile.length() + (steps - boundary) as f64 * (-profile.length() / steps as f64)
        };
        assert_eq!(location.axial_position_m, prefix + offset);
        assert_eq!(location.density_kg_m3, 950.0);
        assert_eq!(location.sound_speed_m_s, 1300.0);
        assert_conical_field(
            sample,
            field(location.axial_position_m),
            location.area_m2,
            950.0 * 1300.0,
        );
    }
    assert!(samples[0].mean_power_w > samples[4].mean_power_w);
    let endpoints = [
        NetworkFieldProbe {
            edge_id: "tube".into(),
            section: NetworkFieldSection::Inlet { segment: 1 },
        },
        NetworkFieldProbe {
            edge_id: "tube".into(),
            section: NetworkFieldSection::Outlet { segment: 1 },
        },
    ];
    let endpoint_owner = PreparedNetworkFields::from_material_with_loss(
        &net, &material, ConeSlices::Auto,
        &request(omega, LossModel::Bessel, &sources, &[], &endpoints), || false,
    )
    .unwrap();
    let endpoint_samples = endpoint_owner
        .sample(endpoint_owner.port_response().default_excitations(), || false)
        .unwrap();
    assert_eq!([samples[0], samples[4]].as_slice(), endpoint_samples.as_slice());
}

// G3: native curved-profile loss quadrature is the model being observed.
// Direct small-panel transfer from the original owner is an independent
// reconstruction route, not an independent discretization-error oracle.
#[test]
fn curved_material_profile_fields_reuse_every_original_magnus_panel() {
    let material = material();
    let profile = RadiusProfile::spline(
        vec![(2.0, 0.003), (2.04, 0.004), (2.11, 0.0035), (2.2, 0.005)],
        None,
    )
    .unwrap();
    let (steps, omega) = (9, 2300.0);
    let reference = 950.0 * 1300.0
        / (core::f64::consts::PI * profile.inlet_radius().powi(2));
    let (blocks, _) = crate::profile::profile_blocks_with_waves_controlled(
        &profile, RiccatiSteps(steps), reference, false, &mut || false,
        &mut |radius, _| crate::fluid::thermoviscous_wave_from_state(
            &material, radius, omega, LossModel::Bessel,
        ),
    )
    .unwrap();
    assert_eq!(blocks.len(), steps);
    let points = panel_probes(0, &(0..=steps).collect::<Vec<_>>());
    let net = profile_net(&profile, steps);
    let sources = [
        pressure_source("left", C64::ZERO),
        pressure_source("right", C64::ZERO),
    ];
    let owner = PreparedNetworkFields::from_material_with_loss(
        &net, &material, ConeSlices::Fixed(1),
        &request(omega, LossModel::Bessel, &sources, &[], &points), || false,
    )
    .unwrap();
    for end in [
        [C64::new(0.8, -0.2), C64::new(4e-11, 3e-11)],
        [C64::new(-0.4, 0.5), C64::new(-2e-11, 1e-11)],
    ] {
        let mut expected = vec![[C64::ZERO; 2]; steps + 1];
        expected[steps] = end;
        for i in (0..steps).rev() {
            let [p, u] = expected[i + 1];
            let [a, b, c, d] = blocks[i];
            expected[i] = [a * p + b * u, c * p + d * u];
        }
        let samples = owner
            .sample(&[expected[0][0], expected[steps][0]], || false)
            .unwrap();
        for ((sample, location), field) in samples.iter().zip(owner.locations()).zip(expected) {
            assert_conical_field(sample, field, location.area_m2, 950.0 * 1300.0);
        }
    }
}

// G0/G4: a requested boundary belongs to the original authored profile and
// original resolution; oversized indices do not create new propagation.
#[test]
fn profile_field_cuts_refuse_wrong_geometry_indices_work_and_cancellation() {
    let gas = gas();
    let profile = RadiusProfile::power(0.008, 0.3, 0.4, 0.2).unwrap();
    let net = profile_net(&profile, 8);
    let sources = [
        pressure_source("left", C64::ONE),
        pressure_source("right", C64::ZERO),
    ];
    for (segment, boundary) in [(0, 9), (0, usize::MAX), (1, 0)] {
        let points = panel_probes(segment, &[boundary]);
        let query = request(1700.0, LossModel::Bessel, &sources, &[], &points);
        assert!(PreparedNetworkFields::new(&net, &gas, &query, || false).is_err());
    }
    let points = panel_probes(0, &[0, 3, 8]);
    let query = request(1700.0, LossModel::Bessel, &sources, &[], &points);
    assert!(PreparedNetworkFields::new(&cylinder_net(0.2), &gas, &query, || false).is_err());
    assert!(matches!(
        PreparedNetworkFields::new(&net, &gas, &query, || true),
        Err(DuctError::Cancelled { .. })
    ));
    let mut limited = request(1700.0, LossModel::Bessel, &sources, &[], &points);
    limited.budget.maximum_chain_blocks = 7;
    assert!(PreparedNetworkFields::new(&net, &gas, &limited, || false).is_err());
    let points = probes(0.2, &[0.5]);
    let query = request(1700.0, LossModel::Bessel, &sources, &[], &points);
    assert!(PreparedNetworkFields::new(&net, &gas, &query, || false).is_err());
}
