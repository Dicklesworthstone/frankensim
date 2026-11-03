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
use fs_math::det;

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
fn stress_column(e: [C64; 3], h: [C64; 3], epsilon: f64, mu: f64) -> [f64; 3] {
    let energy = 0.25
        * (epsilon * e.iter().map(|z| z.norm_sq()).sum::<f64>()
            + mu * h.iter().map(|z| z.norm_sq()).sum::<f64>());
    std::array::from_fn(|i| {
        0.5 * (epsilon * (e[i] * e[2].conj()).re + mu * (h[i] * h[2].conj()).re)
            - if i == 2 { energy } else { 0.0 }
    })
}

// Independent propagating-wave reconstruction: transversality fixes Ez and
// Faraday fixes every H component, including the sign of a backward wave.
fn fresnel_field(
    omega: f64,
    kx: f64,
    kz: f64,
    mu: f64,
    polarization: Polarization,
    amplitude: C64,
) -> ([C64; 3], [C64; 3]) {
    let e = match polarization {
        Polarization::Te => [C64::ZERO, amplitude, C64::ZERO],
        Polarization::Tm => [amplitude, C64::ZERO, amplitude.scale(-kx / kz)],
    };
    let h = [
        e[1].scale(-kz / (omega * mu)),
        (e[0].scale(kz) - e[2].scale(kx)).scale(1.0 / (omega * mu)),
        e[1].scale(kx / (omega * mu)),
    ];
    (e, h)
}

#[test]
fn reflectionless_refraction_transfers_momentum_from_either_face() {
    let left = material(1.0, 1.0);
    let right = material(2.0, 2.0);
    let stack =
        PlanarDielectricLayers::try_new(&left, DielectricResponse::Lossless(&right), 1.0, &[])
            .unwrap();
    // Matched impedance does not imply matched momentum per transmitted watt.
    // Both directions pull toward the same lower-index exterior.
    for polarization in [Polarization::Te, Polarization::Tm] {
        let response = stack.response(1.0, polarization).unwrap();
        assert_eq!(response.reflected_power_fraction, 0.0);
        assert_eq!(response.absorbed_power_fraction, 0.0);
        close(response.transmitted_power_fraction, 1.0);
        let force = stack
            .radiation_pressure(1.0, polarization, C64::ONE)
            .unwrap();
        close(force[2], -0.5);
        assert_eq!([force[0], force[1]], [0.0; 2]);
        let dual = stack.prepare_illumination(1.0, polarization).unwrap();
        for input in [[C64::ONE, C64::ZERO], [C64::ZERO, C64::ONE]] {
            let force = dual.radiation_pressure(input).unwrap();
            close(force[2], -0.5);
            assert_eq!([force[0], force[1]], [0.0; 2]);
        }
        for phase in [C64::ONE, C64::new(0.0, 1.0), -C64::ONE] {
            close(dual.radiation_pressure([C64::ONE, phase]).unwrap()[2], -1.0);
        }
    }

    // A bare nonmagnetic dielectric at TM Brewster incidence has R=A=0,
    // yet its refracted normal momentum produces an exactly known load.
    let dielectric = material(4.0, 1.0);
    let brewster =
        PlanarDielectricLayers::try_new(&left, DielectricResponse::Lossless(&dielectric), 1.0, &[])
            .unwrap();
    let cosine = 1.0 / 5.0_f64.sqrt();
    let response = brewster.response(cosine, Polarization::Tm).unwrap();
    assert!(response.reflected_power_fraction < 1e-28);
    close(response.transmitted_power_fraction, 1.0);
    let force = brewster
        .radiation_pressure(cosine, Polarization::Tm, C64::ONE)
        .unwrap();
    close(force[2], -1.5);
    assert_eq!([force[0], force[1]], [0.0; 2]);
}

#[test]
fn distinct_exterior_fresnel_forces_match_independent_six_component_stress() {
    let omega = 1.7;
    let amplitude = C64::new(0.8, -0.2);
    for (epsilon_l, mu_l, epsilon_r, mu_r) in [
        (1.0, 1.0, 4.0, 1.0),
        (1.0, 1.0, 1.0, 4.0),
        (1.0, 1.0, 4.0, 0.25), // Same index, different magnetic impedance.
        (2.0, 0.5, 3.0, 2.0),
    ] {
        let left = material(epsilon_l, mu_l);
        let right = material(epsilon_r, mu_r);
        let stack = PlanarDielectricLayers::try_new(
            &left,
            DielectricResponse::Lossless(&right),
            omega,
            &[],
        )
        .unwrap();
        for cosine in [1.0_f64, 0.8, 0.3] {
            let k_l = omega * (epsilon_l * mu_l).sqrt();
            let kx = k_l * (1.0 - cosine * cosine).sqrt();
            let kz_l = k_l * cosine;
            let kz_r = (omega * omega * epsilon_r * mu_r - kx * kx).sqrt();
            for polarization in [Polarization::Te, Polarization::Tm] {
                let (y_l, y_r) = match polarization {
                    Polarization::Te => (kz_l / (omega * mu_l), kz_r / (omega * mu_r)),
                    Polarization::Tm => (omega * epsilon_l / kz_l, omega * epsilon_r / kz_r),
                };
                let r = (y_l - y_r) / (y_l + y_r);
                let t_lr = 2.0 * y_l / (y_l + y_r);
                let t_rl = 2.0 * y_r / (y_l + y_r);
                let response = stack.response(cosine, polarization).unwrap();
                complex_near(response.reflected_tangential_e, C64::from_re(r), 3e-14);
                complex_near(response.transmitted_tangential_e, C64::from_re(t_lr), 3e-14);
                let dual = stack.prepare_illumination(cosine, polarization).unwrap();
                for input in [
                    [amplitude, C64::ZERO],
                    [C64::ZERO, amplitude],
                    [amplitude, amplitude],
                    [amplitude, amplitude * C64::new(0.0, 1.0)],
                ] {
                    let outgoing = [
                        input[0].scale(r) + input[1].scale(t_rl),
                        input[0].scale(t_lr) - input[1].scale(r),
                    ];
                    let mut stress = [[0.0; 3]; 2];
                    for (side, epsilon, mu, kz, forward, backward) in [
                        (0, epsilon_l, mu_l, kz_l, input[0], outgoing[0]),
                        (1, epsilon_r, mu_r, kz_r, outgoing[1], input[1]),
                    ] {
                        let (ef, hf) = fresnel_field(omega, kx, kz, mu, polarization, forward);
                        let (eb, hb) = fresnel_field(omega, kx, -kz, mu, polarization, backward);
                        stress[side] = stress_column(
                            std::array::from_fn(|i| ef[i] + eb[i]),
                            std::array::from_fn(|i| hf[i] + hb[i]),
                            epsilon,
                            mu,
                        );
                    }
                    let scale =
                        0.5 * (y_l * input[0].norm_sq() + y_r * input[1].norm_sq()) * (kz_l + kz_r)
                            / omega;
                    let force = dual.radiation_pressure(input).unwrap();
                    for c in 0..3 {
                        near(force[c], stress[1][c] - stress[0][c], 2e-12 * scale);
                    }
                    if input[1] == C64::ZERO {
                        let single = stack
                            .radiation_pressure(cosine, polarization, input[0])
                            .unwrap();
                        for c in 0..3 {
                            near(single[c], stress[1][c] - stress[0][c], 2e-12 * scale);
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn critical_and_evanescent_exit_force_matches_finite_maxwell_fields() {
    let grid = AngularSpectrumGrid {
        width: 4,
        height: 1,
        spacing_m: [core::f64::consts::PI / 0.6, 1.0],
    };
    let sine = core::f64::consts::PI / grid.spacing_m[0];
    let cosine = det::sqrt((1.0 - sine) * (1.0 + sine));
    // Match the owner's admitted cosine exactly, so the middle case is kz=0.
    let critical_epsilon = (1.0 - cosine) * (1.0 + cosine);
    let outside = material(1.0, 1.0);
    let film = material(1.7, 1.0);
    let amplitude = [C64::new(0.7, 0.2), C64::new(-0.3, 0.5)];
    let input: Vec<_> = (0..4)
        .map(|i| amplitude.map(|a| a.scale(if i % 2 == 0 { 1.0 } else { -1.0 })))
        .collect();
    for contrast_delta in [-1e-6, 0.0, 1e-6] {
        let epsilon_r = critical_epsilon + contrast_delta;
        let exit = material(epsilon_r, 1.0);
        for thickness in [0.0, 0.43] {
            let stack = PlanarDielectricLayers::try_new(
                &outside,
                DielectricResponse::Lossless(&exit),
                1.0,
                &[DielectricLayer {
                    response: DielectricResponse::Lossless(&film),
                    thickness_m: thickness,
                }],
            )
            .unwrap();
            for polarization in [Polarization::Te, Polarization::Tm] {
                let response = stack.response(cosine, polarization).unwrap();
                if contrast_delta <= 0.0 {
                    assert_eq!(response.transmitted_power_fraction, 0.0);
                    assert!(stack.prepare_illumination(cosine, polarization).is_err());
                    let y = match polarization {
                        Polarization::Te => cosine,
                        Polarization::Tm => 1.0 / cosine,
                    };
                    let force = stack
                        .radiation_pressure(cosine, polarization, C64::ONE)
                        .unwrap();
                    close(force[2], cosine * y);
                }
                if contrast_delta == 0.0 {
                    assert_eq!(response.exit_normal_wavenumber_rad_m, C64::ZERO);
                }
            }
            let beam =
                PlanarAngularSpectrum::prepare(&stack, grid, &input, limits(), || false).unwrap();
            let force = beam.radiation_force(|| false).unwrap();
            assert_eq!([force[0], force[1]], [0.0; 2]);
            if contrast_delta <= 0.0 {
                close(force[2], 2.0 * cosine * beam.power().incident_w);
            }
            for distance in [0.0, 0.7, 4.0] {
                let left = beam
                    .field(Side::TotalIncidentSide, distance, || false)
                    .unwrap();
                let right = beam.field(Side::Transmitted, distance, || false).unwrap();
                let mut integral = [0.0; 3];
                for i in 0..input.len() {
                    let a = stress_column(
                        left.electric_v_per_m()[i],
                        left.magnetic_a_per_m()[i],
                        1.0,
                        1.0,
                    );
                    let b = stress_column(
                        right.electric_v_per_m()[i],
                        right.magnetic_a_per_m()[i],
                        epsilon_r,
                        1.0,
                    );
                    for c in 0..3 {
                        integral[c] += (b[c] - a[c]) * grid.spacing_m[0] * grid.spacing_m[1];
                    }
                }
                for c in 0..3 {
                    near(force[c], integral[c], 2e-11 * beam.power().incident_w);
                }
            }
        }
    }
}

#[test]
fn vector_beam_force_matches_full_cell_maxwell_stress_and_coherence() {
    let outside = material(E, M);
    let refracting = material(2.25 * E, 0.8 * M);
    let passive = finite_layers::debye(3.0 * E, 4.0 * E, 1.0, 1.0);
    let omega = 1e9;
    for (exit, epsilon_r, mu_r) in [(&outside, E, M), (&refracting, 2.25 * E, 0.8 * M)] {
        let stack = PlanarDielectricLayers::try_new(
            &outside,
            DielectricResponse::Lossless(exit),
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
        let beam =
            PlanarAngularSpectrum::prepare(&stack, grid, &input, limits(), || false).unwrap();
        let force = beam.radiation_force(|| false).unwrap();
        assert!(force[0] > 0.0 && force[1] < 0.0);
        if epsilon_r == E && mu_r == M {
            assert!(force[2] > 0.0);
        }
        for distance in [0.0, 1.1] {
            let left = beam
                .field(Side::TotalIncidentSide, distance, || false)
                .unwrap();
            let right = beam.field(Side::Transmitted, distance, || false).unwrap();
            let mut integral = [0.0; 3];
            for i in 0..input.len() {
                let a = stress_column(left.electric_v_per_m()[i], left.magnetic_a_per_m()[i], E, M);
                let b = stress_column(
                    right.electric_v_per_m()[i],
                    right.magnetic_a_per_m()[i],
                    epsilon_r,
                    mu_r,
                );
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
        let bank =
            AngularSpectrumCoherence::try_new(&factors[0], &factors, budget, || false).unwrap();
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
}

#[test]
fn force_rejects_lossy_exterior_and_unrepresentable_source() {
    let outside = material(E, M);
    let lossy = finite_layers::debye(3.0 * E, 4.0 * E, 1.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&outside, DielectricResponse::Debye(&lossy), 1e9, &[])
            .unwrap();
    assert!(
        stack
            .radiation_pressure(1.0, Polarization::Te, C64::ONE)
            .is_err()
    );
    assert!(stack.prepare_illumination(1.0, Polarization::Te).is_err());
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
