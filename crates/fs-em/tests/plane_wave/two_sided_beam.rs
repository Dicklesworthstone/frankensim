//! Two-face beam fields, independent interface/characteristic-matrix oracles,
//! positive interior loss and coherent momentum. Synthetic constitutive data.
use super::*;
use fs_em::plane_wave::{
    angular_spectrum::{
        AngularSpectrumError, AngularSpectrumGrid, AngularSpectrumIlluminationLimits,
        AngularSpectrumIlluminationSide as DualSide, AngularSpectrumLayerLimits,
        AngularSpectrumLimits, AngularSpectrumSide as Side, PlanarAngularSpectrum,
    },
    layers::{DielectricLayer, PlanarDielectricLayers},
};
use fs_math::det;

fn limits() -> AngularSpectrumLimits {
    AngularSpectrumLimits {
        maximum_grid_points: 64,
        maximum_workspace_bytes: 1 << 20,
        maximum_fft_work: 1 << 20,
    }
}
fn dual_limits() -> AngularSpectrumIlluminationLimits {
    AngularSpectrumIlluminationLimits {
        maximum_layer_mode_pairs: 1024,
        maximum_workspace_bytes: 2 << 20,
        maximum_fft_work: 1 << 20,
    }
}
fn grid(spacing: f64) -> AngularSpectrumGrid {
    AngularSpectrumGrid {
        width: 4,
        height: 4,
        spacing_m: [spacing, spacing],
    }
}
fn phase(angle: f64) -> C64 {
    C64::new(det::cos(angle), det::sin(angle))
}
fn close(a: f64, b: f64) {
    near(a, b, 3e-11 * a.abs().max(b.abs()).max(1e-24));
}

#[test]
fn coherent_counterbeams_match_direct_fourier_maxwell_fields_at_both_faces() {
    let left = material(1.0, 1.0);
    let right = material(2.25, 0.8);
    let stack =
        PlanarDielectricLayers::try_new(&left, DielectricResponse::Lossless(&right), 1.0, &[])
            .unwrap();
    let g = grid(8.0);
    // Authored signed spatial modes, including polarization mixing and a
    // mode launched from only one face. Equal bins interfere coherently.
    let specs = [
        (
            [1.0, 0.0],
            [C64::new(0.3, 0.1), C64::new(-0.2, 0.7)],
            [C64::new(0.4, -0.3), C64::new(0.1, 0.2)],
        ),
        (
            [-1.0, 1.0],
            [C64::new(0.1, -0.4), C64::new(0.6, 0.2)],
            [C64::new(-0.5, 0.2), C64::new(0.2, 0.4)],
        ),
        (
            [0.0, 0.0],
            [C64::ZERO; 2],
            [C64::new(0.2, -0.1), C64::new(0.3, 0.2)],
        ),
    ];
    let mut input = [vec![[C64::ZERO; 2]; 16], vec![[C64::ZERO; 2]; 16]];
    for (bins, a, b) in specs {
        for y in 0..4 {
            for x in 0..4 {
                let p = phase(
                    -core::f64::consts::TAU * (bins[0] * x as f64 + bins[1] * y as f64) / 4.0,
                );
                for c in 0..2 {
                    input[0][y * 4 + x][c] = input[0][y * 4 + x][c] + a[c] * p;
                    input[1][y * 4 + x][c] = input[1][y * 4 + x][c] + b[c] * p;
                }
            }
        }
    }
    let beam = PlanarAngularSpectrum::prepare(&stack, g, &input[0], limits(), || false).unwrap();
    let dual = beam
        .prepare_illumination(&input[1], dual_limits(), || false)
        .unwrap();
    assert!(std::ptr::eq(dual.left_beam(), &beam));
    assert!(std::ptr::eq(dual.source(), &stack));
    assert_eq!(dual.power().film_absorbed_w, 0.0);
    assert_eq!(dual.layer_absorbed_power_w(0), None);
    let force = dual.radiation_force(|| false).unwrap();
    assert_eq!([force[0], force[1]], [0.0; 2]);
    for distance in [0.0, 0.7, 9.0] {
        let mut integrated_stress = [0.0; 3];
        for side in [
            DualSide::IncomingLeft,
            DualSide::OutgoingLeft,
            DualSide::TotalLeft,
            DualSide::IncomingRight,
            DualSide::OutgoingRight,
            DualSide::TotalRight,
        ] {
            let actual = dual.field(side, distance, || false).unwrap();
            let mut e = vec![[C64::ZERO; 3]; 16];
            let mut h = e.clone();
            for (bins, a, b) in specs {
                let kx = core::f64::consts::TAU * bins[0] / 32.0;
                let ky = core::f64::consts::TAU * bins[1] / 32.0;
                let kt2 = kx * kx + ky * ky;
                let kz = [(1.0 - kt2).sqrt(), (1.8 - kt2).sqrt()];
                let (px, py) = if kt2 == 0.0 {
                    (1.0, 0.0)
                } else {
                    (kx / kt2.sqrt(), ky / kt2.sqrt())
                };
                let split = |v: [C64; 2]| {
                    [
                        v[1].scale(px) - v[0].scale(py),
                        v[0].scale(px) + v[1].scale(py),
                    ]
                };
                let compose = |v: [C64; 2]| {
                    [
                        v[1].scale(px) - v[0].scale(py),
                        v[1].scale(py) + v[0].scale(px),
                    ]
                };
                let al = split(a);
                let ar = split(b);
                let yl = [kz[0], 1.0 / kz[0]];
                let yr = [kz[1] / 0.8, 2.25 / kz[1]];
                let out_l = std::array::from_fn(|p| {
                    (al[p].scale(yl[p] - yr[p]) + ar[p].scale(2.0 * yr[p]))
                        .scale(1.0 / (yl[p] + yr[p]))
                });
                let out_r = std::array::from_fn(|p| {
                    (al[p].scale(2.0 * yl[p]) + ar[p].scale(yr[p] - yl[p]))
                        .scale(1.0 / (yl[p] + yr[p]))
                });
                // Independent fixed-axis E transversality and H=k cross E/(omega mu).
                let terms = [
                    (
                        matches!(side, DualSide::IncomingLeft | DualSide::TotalLeft),
                        kz[0],
                        1.0,
                        a,
                        -distance,
                    ),
                    (
                        matches!(side, DualSide::OutgoingLeft | DualSide::TotalLeft),
                        -kz[0],
                        1.0,
                        compose(out_l),
                        -distance,
                    ),
                    (
                        matches!(side, DualSide::IncomingRight | DualSide::TotalRight),
                        -kz[1],
                        0.8,
                        b,
                        distance,
                    ),
                    (
                        matches!(side, DualSide::OutgoingRight | DualSide::TotalRight),
                        kz[1],
                        0.8,
                        compose(out_r),
                        distance,
                    ),
                ];
                for (active, q, mu, transverse, z) in terms {
                    if !active {
                        continue;
                    }
                    let em = [
                        transverse[0],
                        transverse[1],
                        -(transverse[0].scale(kx) + transverse[1].scale(ky)).scale(1.0 / q),
                    ];
                    let hm = [
                        (em[2].scale(ky) - em[1].scale(q)).scale(1.0 / mu),
                        (em[0].scale(q) - em[2].scale(kx)).scale(1.0 / mu),
                        (em[1].scale(kx) - em[0].scale(ky)).scale(1.0 / mu),
                    ];
                    for y in 0..4 {
                        for x in 0..4 {
                            let ph = phase(-(kx * 8.0 * x as f64 + ky * 8.0 * y as f64 + q * z));
                            for c in 0..3 {
                                e[y * 4 + x][c] = e[y * 4 + x][c] + em[c] * ph;
                                h[y * 4 + x][c] = h[y * 4 + x][c] + hm[c] * ph;
                            }
                        }
                    }
                }
            }
            for i in 0..16 {
                for c in 0..3 {
                    complex_near(actual.electric_v_per_m()[i][c], e[i][c], 3e-12);
                    complex_near(actual.magnetic_a_per_m()[i][c], h[i][c], 3e-12);
                }
            }
            // Integrate the independently reconstructed fields above, with
            // each face's own constitutive law. No scattering-power shortcut.
            let exterior = match side {
                DualSide::TotalLeft => Some((1.0, 1.0, -1.0)),
                DualSide::TotalRight => Some((2.25, 0.8, 1.0)),
                _ => None,
            };
            if let Some((epsilon, mu, normal)) = exterior {
                for i in 0..16 {
                    let energy = 0.25 * (epsilon * e[i].iter().map(|z| z.norm_sq()).sum::<f64>()
                        + mu * h[i].iter().map(|z| z.norm_sq()).sum::<f64>());
                    for c in 0..3 {
                        let traction = 0.5 * (epsilon * (e[i][c] * e[i][2].conj()).re
                            + mu * (h[i][c] * h[i][2].conj()).re)
                            - if c == 2 { energy } else { 0.0 };
                        integrated_stress[c] += normal * traction * g.spacing_m[0] * g.spacing_m[1];
                    }
                }
            }
            let p = dual.power();
            let expected = match side {
                DualSide::IncomingLeft => p.incoming_w[0],
                DualSide::OutgoingLeft => -p.outgoing_w[0],
                DualSide::TotalLeft => p.incoming_w[0] - p.outgoing_w[0],
                DualSide::IncomingRight => -p.incoming_w[1],
                DualSide::OutgoingRight => p.outgoing_w[1],
                DualSide::TotalRight => p.outgoing_w[1] - p.incoming_w[1],
            };
            near(
                actual.normal_power_w(),
                expected,
                1e-10 * (p.incoming_w[0] + p.incoming_w[1]),
            );
        }
        let scale =
            (dual.power().incoming_w[0] + dual.power().incoming_w[1]) * (1.0 + 1.8_f64.sqrt());
        for c in 0..3 {
            near(force[c], integrated_stress[c], 3e-11 * scale);
        }
    }
}

#[test]
fn counterbeam_phase_drives_original_interior_heat_and_exterior_maxwell_force() {
    let outside = material(E, M);
    let film = finite_layers::debye(3.0 * E, 4.0 * E, 1.0, 1.0);
    let omega = 1e9;
    let k = omega * (E * M).sqrt();
    let depth = 0.3 / k;
    let stack = PlanarDielectricLayers::try_new(
        &outside,
        DielectricResponse::Lossless(&outside),
        omega,
        &[DielectricLayer {
            response: DielectricResponse::Debye(&film),
            thickness_m: depth,
        }],
    )
    .unwrap();
    let g = grid(8.0 / k);
    let area = 16.0 * g.spacing_m[0] * g.spacing_m[1];
    let incoming = vec![[C64::ZERO, C64::new(0.8, 0.2)]; 16];
    let beam = PlanarAngularSpectrum::prepare(&stack, g, &incoming, limits(), || false).unwrap();
    let mut heats = Vec::new();
    for phi in [0.0, 0.8, core::f64::consts::PI] {
        let input = [incoming[0][1], phase(phi).scale(0.7)];
        let right = vec![[C64::ZERO, input[1]]; 16];
        let dual = beam
            .prepare_illumination(&right, dual_limits(), || false)
            .unwrap();
        let m = illumination::matrix(C64::new(5.0, -2.0), 1.0, 0.3, 1.0, Polarization::Te);
        let plus = illumination::act(m, [C64::ONE, C64::ONE]);
        let minus = illumination::act(m, [C64::ONE, -C64::ONE]);
        let out_r = (input[0].scale(2.0) - (minus[0] + minus[1]) * input[1]) / (plus[0] + plus[1]);
        let out_l = plus[0] * out_r + minus[0] * input[1] - input[0];
        let y = (E / M).sqrt();
        let incident = 0.5 * area * y * (input[0].norm_sq() + input[1].norm_sq());
        let outgoing = 0.5 * area * y * (out_l.norm_sq() + out_r.norm_sq());
        close(dual.power().film_absorbed_w, incident - outgoing);
        close(
            dual.layer_absorbed_power_w(0).unwrap(),
            dual.power().film_absorbed_w,
        );
        let boundary = [out_r + input[1], out_r - input[1]];
        let mut integral = 0.0;
        const STEPS: usize = 128;
        for i in 0..=STEPS {
            let z = depth * i as f64 / STEPS as f64;
            let expected = illumination::act(
                illumination::matrix(
                    C64::new(5.0, -2.0),
                    1.0,
                    k * (depth - z),
                    1.0,
                    Polarization::Te,
                ),
                boundary,
            );
            let field = dual.layer_field(0, z, || false).unwrap();
            for index in 0..16 {
                complex_near(
                    field.field().electric_v_per_m()[index][1],
                    expected[0],
                    2e-12,
                );
                complex_near(
                    field.field().magnetic_a_per_m()[index][0],
                    -expected[1].scale(y),
                    2e-14,
                );
                close(
                    field.absorbed_power_density_w_m3()[index],
                    omega * E * expected[0].norm_sq(),
                );
            }
            integral += field.absorbed_power_per_depth_w_m()
                * if i == 0 || i == STEPS {
                    1.0
                } else if i % 2 == 0 {
                    2.0
                } else {
                    4.0
                };
        }
        near(
            integral * depth / (3 * STEPS) as f64,
            dual.power().film_absorbed_w,
            2e-10 * incident,
        );
        let left = dual.field(DualSide::TotalLeft, 0.0, || false).unwrap();
        let right = dual.field(DualSide::TotalRight, 0.0, || false).unwrap();
        let stress = |f: &fs_em::plane_wave::angular_spectrum::AngularSpectrumField| {
            let e = f.electric_v_per_m()[0];
            let h = f.magnetic_a_per_m()[0];
            0.5 * (E * e[2].norm_sq() + M * h[2].norm_sq())
                - 0.25
                    * (E * e.iter().map(|v| v.norm_sq()).sum::<f64>()
                        + M * h.iter().map(|v| v.norm_sq()).sum::<f64>())
        };
        let force = dual.radiation_force(|| false).unwrap();
        assert_eq!(force[0], 0.0);
        assert_eq!(force[1], 0.0);
        near(
            force[2],
            area * (stress(&right) - stress(&left)),
            2e-12 * incident * (E * M).sqrt(),
        );
        heats.push(dual.power().film_absorbed_w);
    }
    assert!(
        heats[0] > 3.0 * heats[2],
        "source phase must change actual absorbed heat"
    );
}

#[test]
fn refracting_absorbing_counterbeams_match_each_face_maxwell_stress() {
    let left = material(E, M);
    let right = material(2.25 * E, 0.8 * M);
    let film = finite_layers::debye(3.0 * E, 4.0 * E, 1.0, 1.0);
    let omega = 1e9;
    let inverse_speed = (E * M).sqrt();
    let stack = PlanarDielectricLayers::try_new(
        &left,
        DielectricResponse::Lossless(&right),
        omega,
        &[DielectricLayer {
            response: DielectricResponse::Debye(&film),
            thickness_m: 0.3 / (omega * inverse_speed),
        }],
    )
    .unwrap();
    let g = grid(8.0 / (omega * inverse_speed));
    let area = 16.0 * g.spacing_m[0] * g.spacing_m[1];
    let incoming = vec![[C64::ZERO, C64::new(0.8, 0.2)]; 16];
    let beam = PlanarAngularSpectrum::prepare(&stack, g, &incoming, limits(), || false).unwrap();
    let mut observed = [0.0; 2];
    let mut scale = 0.0;
    for (phase_index, phi) in [0.0, core::f64::consts::PI].into_iter().enumerate() {
        let right_input = vec![[C64::ZERO, phase(phi).scale(0.7)]; 16];
        let dual = beam
            .prepare_illumination(&right_input, dual_limits(), || false)
            .unwrap();
        let power = dual.power();
        assert!(power.film_absorbed_w > 0.0);
        scale = (power.incoming_w[0] + power.incoming_w[1]) * inverse_speed;
        let force = dual.radiation_force(|| false).unwrap();
        assert_eq!([force[0], force[1]], [0.0; 2]);
        near(
            force[2],
            inverse_speed * (power.incoming_w[0] + power.outgoing_w[0])
                - (1.8_f64).sqrt() * inverse_speed * (power.incoming_w[1] + power.outgoing_w[1]),
            2e-12 * scale,
        );
        for distance in [0.0, 0.3] {
            let mut integral = [0.0; 3];
            for (side, epsilon, mu, normal) in [
                (DualSide::TotalLeft, E, M, -1.0),
                (DualSide::TotalRight, 2.25 * E, 0.8 * M, 1.0),
            ] {
                let field = dual.field(side, distance, || false).unwrap();
                for i in 0..16 {
                    let e = field.electric_v_per_m()[i];
                    let h = field.magnetic_a_per_m()[i];
                    let energy = 0.25
                        * (epsilon * e.iter().map(|z| z.norm_sq()).sum::<f64>()
                            + mu * h.iter().map(|z| z.norm_sq()).sum::<f64>());
                    for c in 0..3 {
                        let stress = 0.5
                            * (epsilon * (e[c] * e[2].conj()).re + mu * (h[c] * h[2].conj()).re)
                            - if c == 2 { energy } else { 0.0 };
                        integral[c] += normal * stress * area / 16.0;
                    }
                }
            }
            for c in 0..3 {
                near(force[c], integral[c], 2e-12 * scale);
            }
        }
        observed[phase_index] = force[2];
    }
    assert!((observed[0] - observed[1]).abs() > 1e-3 * scale);
}

#[test]
fn dark_counterbeam_recovers_native_fields_heat_and_force_with_exact_budget_boundaries() {
    let outside = material(E, M);
    let film = finite_layers::debye(3.0 * E, 4.0 * E, 1.0, 1.0);
    let stack = PlanarDielectricLayers::try_new(
        &outside,
        DielectricResponse::Lossless(&outside),
        1e9,
        &[DielectricLayer {
            response: DielectricResponse::Debye(&film),
            thickness_m: 0.13,
        }],
    )
    .unwrap();
    let g = grid(3.0);
    let input: Vec<_> = (0..16)
        .map(|i| [phase(i as f64 * 0.4), phase(i as f64 * -0.6).scale(0.3)])
        .collect();
    let dark = vec![[C64::ZERO; 2]; 16];
    let beam = PlanarAngularSpectrum::prepare(&stack, g, &input, limits(), || false).unwrap();
    let dual = beam
        .prepare_illumination(&dark, dual_limits(), || false)
        .unwrap();
    close(dual.power().incoming_w[0], beam.power().incident_w);
    assert_eq!(dual.power().incoming_w[1], 0.0);
    close(dual.power().outgoing_w[0], beam.power().reflected_w);
    close(dual.power().outgoing_w[1], beam.power().transmitted_w);
    close(dual.power().film_absorbed_w, beam.power().film_absorbed_w);
    for (s, d) in [
        (Side::Incident, DualSide::IncomingLeft),
        (Side::Reflected, DualSide::OutgoingLeft),
        (Side::TotalIncidentSide, DualSide::TotalLeft),
        (Side::Transmitted, DualSide::TotalRight),
    ] {
        let a = beam.field(s, 0.3, || false).unwrap();
        let b = dual.field(d, 0.3, || false).unwrap();
        for i in 0..16 {
            for c in 0..3 {
                complex_near(
                    a.electric_v_per_m()[i][c],
                    b.electric_v_per_m()[i][c],
                    3e-12,
                );
                complex_near(
                    a.magnetic_a_per_m()[i][c],
                    b.magnetic_a_per_m()[i][c],
                    3e-14,
                );
            }
        }
    }
    let original = beam
        .prepare_layer(
            0,
            AngularSpectrumLayerLimits {
                maximum_layer_mode_pairs: 1024,
                maximum_workspace_bytes: 2 << 20,
                maximum_fft_work: 1 << 20,
            },
            || false,
        )
        .unwrap();
    for z in [0.0, 0.047, 0.13] {
        let a = original.field(z, || false).unwrap();
        let b = dual.layer_field(0, z, || false).unwrap();
        for i in 0..16 {
            close(
                a.absorbed_power_density_w_m3()[i],
                b.absorbed_power_density_w_m3()[i],
            );
        }
    }
    let a = beam.radiation_force(|| false).unwrap();
    let b = dual.radiation_force(|| false).unwrap();
    for c in 0..3 {
        near(a[c], b[c], 1e-11 * beam.power().incident_w * (E * M).sqrt());
    }
    let exact = AngularSpectrumIlluminationLimits {
        maximum_layer_mode_pairs: dual.preparation_layer_mode_pairs(),
        maximum_workspace_bytes: dual.workspace_bytes(),
        maximum_fft_work: beam.field_fft_work(),
    };
    assert!(beam.prepare_illumination(&dark, exact, || false).is_ok());
    for bounded in [
        AngularSpectrumIlluminationLimits {
            maximum_layer_mode_pairs: exact.maximum_layer_mode_pairs - 1,
            ..exact
        },
        AngularSpectrumIlluminationLimits {
            maximum_workspace_bytes: exact.maximum_workspace_bytes - 1,
            ..exact
        },
        AngularSpectrumIlluminationLimits {
            maximum_fft_work: exact.maximum_fft_work - 1,
            ..exact
        },
    ] {
        assert!(beam.prepare_illumination(&dark, bounded, || false).is_err());
    }
    assert!(matches!(
        beam.prepare_illumination(&dark, exact, || true),
        Err(AngularSpectrumError::Cancelled)
    ));
    assert!(matches!(
        dual.field(DualSide::TotalRight, 0.0, || true),
        Err(AngularSpectrumError::Cancelled)
    ));
    assert!(matches!(
        dual.layer_field(0, 0.03, || true),
        Err(AngularSpectrumError::Cancelled)
    ));
    assert!(matches!(
        dual.radiation_force(|| true),
        Err(AngularSpectrumError::Cancelled)
    ));
    assert!(dual.layer_field(0, -0.01, || false).is_err());
    assert!(dual.layer_field(1, 0.0, || false).is_err());
    assert!(dual.field(DualSide::TotalLeft, f64::NAN, || false).is_err());
    assert!(
        beam.prepare_illumination(&dark[..15], exact, || false)
            .is_err()
    );
    let mut invalid = dark.clone();
    invalid[3][1] = C64::new(f64::NAN, 0.0);
    assert!(
        beam.prepare_illumination(&invalid, exact, || false)
            .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        beam.prepare_illumination(&dark, exact, || {
            polls += 1;
            polls > 15
        }),
        Err(AngularSpectrumError::Cancelled)
    ));
    assert_eq!(
        beam.prepare_illumination(&dark, exact, || false)
            .unwrap()
            .power(),
        dual.power()
    );
}

#[test]
fn counterbeam_refuses_lossy_or_nonpropagating_right_exterior_even_when_dark() {
    let dense = material(4.0, 1.0);
    let sparse = material(1.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&dense, DielectricResponse::Lossless(&sparse), 1.0, &[])
            .unwrap();
    let g = AngularSpectrumGrid {
        width: 4,
        height: 1,
        spacing_m: [2.5, 1.0],
    };
    let input = vec![[C64::ONE, C64::ZERO]; 4];
    let beam = PlanarAngularSpectrum::prepare(&stack, g, &input, limits(), || false).unwrap();
    assert!(beam.radiation_force(|| false).is_ok());
    assert!(
        beam.prepare_illumination(&[[C64::ZERO; 2]; 4], dual_limits(), || false)
            .is_err()
    );
    let outside = material(E, M);
    let lossy = finite_layers::debye(3.0 * E, 4.0 * E, 1.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&outside, DielectricResponse::Debye(&lossy), 1e9, &[])
            .unwrap();
    let beam = PlanarAngularSpectrum::prepare(
        &stack,
        grid(3.0),
        &vec![[C64::ZERO; 2]; 16],
        limits(),
        || false,
    )
    .unwrap();
    assert!(beam.radiation_force(|| false).is_err());
    assert!(
        beam.prepare_illumination(&vec![[C64::ZERO; 2]; 16], dual_limits(), || false)
            .is_err()
    );
}
