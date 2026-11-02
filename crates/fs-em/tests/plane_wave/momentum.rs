//! Original-film radiation pressure and independent exterior Maxwell stress.
use super::*;
use fs_em::plane_wave::{
    angular_spectrum::{
        AngularSpectrumCoherence, AngularSpectrumCoherenceLimits, AngularSpectrumError,
        AngularSpectrumGrid, AngularSpectrumLimits, AngularSpectrumSide as Side,
        PlanarAngularSpectrum,
    },
    layers::{DielectricLayer, PlanarDielectricLayers},
};

fn close(a: f64, b: f64) {
    near(a, b, 2e-11 * a.abs().max(b.abs()).max(1e-24));
}
fn limits() -> AngularSpectrumLimits {
    AngularSpectrumLimits {
        maximum_grid_points: 64,
        maximum_workspace_bytes: 1 << 20,
        maximum_fft_work: 1 << 20,
    }
}

#[test]
fn original_film_force_retains_fresnel_resonances_and_oblique_absorption() {
    let outside = material(E, M);
    let film = material(4.0 * E, M);
    let passive = finite_layers::debye(3.0 * E, 4.0 * E, 1.0, 1.0);
    let inverse_speed = (E * M).sqrt();
    let omega = 1e9;
    for response in [
        DielectricResponse::Lossless(&film),
        DielectricResponse::Debye(&passive),
    ] {
        for depth in [0.0, 0.3, core::f64::consts::PI / 2.0, 20.0] {
            let stack = PlanarDielectricLayers::try_new(
                &outside,
                DielectricResponse::Lossless(&outside),
                omega,
                &[DielectricLayer {
                    response,
                    thickness_m: depth / (omega * inverse_speed),
                }],
            )
            .unwrap();
            for cosine in [1.0, 0.6] {
                for pol in [Polarization::Te, Polarization::Tm] {
                    let source = C64::new(0.8, -0.2);
                    let r = stack.response(cosine, pol).unwrap();
                    let y = (E / M).sqrt()
                        * match pol {
                            Polarization::Te => cosine,
                            Polarization::Tm => 1.0 / cosine,
                        };
                    let incoming = 0.5 * y * source.norm_sq();
                    let pressure = stack.radiation_pressure(cosine, pol, source).unwrap();
                    close(
                        pressure[0],
                        incoming
                            * inverse_speed
                            * (1.0 - cosine * cosine).sqrt()
                            * r.absorbed_power_fraction,
                    );
                    near(
                        pressure[2],
                        incoming
                            * inverse_speed
                            * cosine
                            * (1.0 + r.reflected_power_fraction - r.transmitted_power_fraction),
                        256.0 * f64::EPSILON * incoming * inverse_speed,
                    );
                    assert_eq!(pressure[1], 0.0);
                    let dual = stack.prepare_illumination(cosine, pol).unwrap();
                    let left = dual.radiation_pressure([source, C64::ZERO]).unwrap();
                    let right = dual.radiation_pressure([C64::ZERO, source]).unwrap();
                    for i in 0..3 {
                        close(left[i], pressure[i]);
                    }
                    close(right[0], left[0]);
                    close(right[2], -left[2]);
                }
            }
        }
    }
}

#[test]
fn simultaneous_sources_retain_phase_dependent_signed_normal_force() {
    let outside = material(1.0, 1.0);
    let film = material(4.0, 1.0);
    let stack = PlanarDielectricLayers::try_new(
        &outside,
        DielectricResponse::Lossless(&outside),
        1.0,
        &[DielectricLayer {
            response: DielectricResponse::Lossless(&film),
            thickness_m: 0.4,
        }],
    )
    .unwrap();
    for pol in [Polarization::Te, Polarization::Tm] {
        let operator = stack.prepare_illumination(0.7, pol).unwrap();
        let mut observed = [0.0; 2];
        for (i, phase) in [C64::new(0.0, 1.0), C64::new(0.0, -1.0)]
            .into_iter()
            .enumerate()
        {
            let input = [C64::ONE, phase];
            let response = operator.evaluate(input, &mut [0.0]).unwrap();
            let p = response.power;
            let force = operator.radiation_pressure(input).unwrap();
            let expected = 0.7
                * (p.incoming_w_m2[0] + p.outgoing_w_m2[0]
                    - p.incoming_w_m2[1]
                    - p.outgoing_w_m2[1]);
            close(force[2], expected);
            assert_eq!(force[0], 0.0);
            observed[i] = force[2];
        }
        assert!(observed[0].abs() > 1e-4);
        close(observed[0], -observed[1]);
    }
}

// T_iz = 1/2 Re(epsilon E_i E_z* + mu H_i H_z*)
//        - delta_iz (epsilon |E|² + mu |H|²)/4.
fn stress_column(e: [C64; 3], h: [C64; 3]) -> [f64; 3] {
    let energy = 0.25
        * (E * e.iter().map(|z| z.norm_sq()).sum::<f64>()
            + M * h.iter().map(|z| z.norm_sq()).sum::<f64>());
    std::array::from_fn(|i| {
        0.5 * (E * (e[i] * e[2].conj()).re + M * (h[i] * h[2].conj()).re)
            - if i == 2 { energy } else { 0.0 }
    })
}

#[test]
fn vector_beam_force_matches_full_cell_maxwell_stress_and_coherence() {
    let outside = material(E, M);
    let passive = finite_layers::debye(3.0 * E, 4.0 * E, 1.0, 1.0);
    let omega = 1e9;
    let stack = PlanarDielectricLayers::try_new(
        &outside,
        DielectricResponse::Lossless(&outside),
        omega,
        &[DielectricLayer {
            response: DielectricResponse::Debye(&passive),
            thickness_m: 0.2,
        }],
    )
    .unwrap();
    let grid = AngularSpectrumGrid {
        width: 8,
        height: 4,
        spacing_m: [2.0, 3.0],
    };
    let mut input = Vec::new();
    for row in 0..grid.height {
        for column in 0..grid.width {
            let phase = -core::f64::consts::TAU * (column as f64 / 8.0 - row as f64 / 4.0);
            let z = C64::new(phase.cos(), phase.sin());
            input.push([
                C64::new(0.2, 0.4) + z.scale(0.7),
                C64::ONE + z * C64::new(-0.3, 0.2),
            ]);
        }
    }
    let beam = PlanarAngularSpectrum::prepare(&stack, grid, &input, limits(), || false).unwrap();
    let force = beam.radiation_force(|| false).unwrap();
    assert!(force[0] > 0.0 && force[1] < 0.0 && force[2] > 0.0);
    for distance in [0.0, 1.1] {
        let left = beam
            .field(Side::TotalIncidentSide, distance, || false)
            .unwrap();
        let right = beam.field(Side::Transmitted, distance, || false).unwrap();
        let mut integral = [0.0; 3];
        for i in 0..input.len() {
            let a = stress_column(left.electric_v_per_m()[i], left.magnetic_a_per_m()[i]);
            let b = stress_column(right.electric_v_per_m()[i], right.magnetic_a_per_m()[i]);
            for j in 0..3 {
                integral[j] += (b[j] - a[j]) * grid.spacing_m[0] * grid.spacing_m[1];
            }
        }
        for j in 0..3 {
            close(force[j], integral[j]);
        }
    }
    let scaled = input
        .iter()
        .map(|x| x.map(|z| z.scale(0.5)))
        .collect::<Vec<_>>();
    let factors = [
        beam,
        PlanarAngularSpectrum::prepare(&stack, grid, &scaled, limits(), || false).unwrap(),
    ];
    let budget = AngularSpectrumCoherenceLimits {
        maximum_modes: 2,
        maximum_point_mode_pairs: 64,
        maximum_workspace_bytes: 2 << 20,
        maximum_fft_work: 1 << 20,
    };
    let bank = AngularSpectrumCoherence::try_new(&factors[0], &factors, budget, || false).unwrap();
    let ensemble = bank.radiation_force(|| false).unwrap();
    for j in 0..3 {
        close(ensemble[j], 1.25 * force[j]);
    }
    let dark = AngularSpectrumCoherence::try_new(&factors[0], &[], budget, || false).unwrap();
    assert_eq!(dark.radiation_force(|| false).unwrap(), [0.0; 3]);
    assert!(matches!(
        bank.radiation_force(|| true),
        Err(AngularSpectrumError::Cancelled)
    ));
    assert!(matches!(
        factors[0].radiation_force(|| true),
        Err(AngularSpectrumError::Cancelled)
    ));
    assert_eq!(factors[0].radiation_force(|| false).unwrap(), force);
}

#[test]
fn force_rejects_a_different_exterior_and_unrepresentable_source() {
    let outside = material(1.0, 1.0);
    let other = material(2.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&outside, DielectricResponse::Lossless(&other), 1.0, &[])
            .unwrap();
    assert!(
        stack
            .radiation_pressure(1.0, Polarization::Te, C64::ONE)
            .is_err()
    );
    assert!(
        stack
            .prepare_illumination(1.0, Polarization::Te)
            .unwrap()
            .radiation_pressure([C64::ONE, C64::ZERO])
            .is_err()
    );
    let grid = AngularSpectrumGrid {
        width: 1,
        height: 1,
        spacing_m: [1.0; 2],
    };
    let beam =
        PlanarAngularSpectrum::prepare(&stack, grid, &[[C64::ONE, C64::ZERO]], limits(), || false)
            .unwrap();
    assert!(beam.radiation_force(|| false).is_err());
    let identity =
        PlanarDielectricLayers::try_new(&outside, DielectricResponse::Lossless(&outside), 1.0, &[])
            .unwrap();
    assert_eq!(
        identity
            .radiation_pressure(0.3, Polarization::Tm, C64::ONE)
            .unwrap(),
        [0.0; 3]
    );
    assert!(
        identity
            .radiation_pressure(1.0, Polarization::Te, C64::new(f64::NAN, 0.0))
            .is_err()
    );
}
