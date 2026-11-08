use super::*;
use fs_em::plane_wave::angular_spectrum::{
    AngularSpectrumLayerLimits, AngularSpectrumLayerObservation, AngularSpectrumObservation,
    AngularSpectrumProbeField, AngularSpectrumProbeLimits,
};

fn probe_limits() -> AngularSpectrumProbeLimits {
    AngularSpectrumProbeLimits {
        maximum_points: 1024,
        maximum_mode_evaluations: 1 << 16,
        maximum_workspace_bytes: 16 << 20,
    }
}

fn cross(e: [C64; 3], h: [C64; 3]) -> [C64; 3] {
    [(1, 2), (2, 0), (0, 1)].map(|(a, b)| (e[a] * h[b].conj() - e[b] * h[a].conj()).scale(0.5))
}

fn check_field(actual: AngularSpectrumProbeField, e: [C64; 3], h: [C64; 3]) {
    let scale = |v: [C64; 3]| {
        v.iter()
            .map(|z| z.re.abs() + z.im.abs())
            .fold(1e-20, f64::max)
    };
    let es = scale(e);
    let hs = scale(h);
    let flux = cross(e, h);
    for c in 0..3 {
        complex_near(actual.electric_v_per_m[c], e[c], 3e-12 * es);
        complex_near(actual.magnetic_a_per_m[c], h[c], 3e-12 * hs);
        complex_near(actual.complex_poynting_w_m2[c], flux[c], 6e-12 * es * hs);
    }
}

fn transverse(grid: AngularSpectrumGrid, carrier: [f64; 2], index: usize) -> [f64; 2] {
    let bin = |j: usize, count: usize| {
        if j < count.div_ceil(2) {
            j as f64
        } else {
            j as f64 - count as f64
        }
    };
    [
        carrier[0]
            + core::f64::consts::TAU * bin(index % grid.width, grid.width)
                / (grid.width as f64 * grid.spacing_m[0]),
        carrier[1]
            + core::f64::consts::TAU * bin(index / grid.width, grid.height)
                / (grid.height as f64 * grid.spacing_m[1]),
    ]
}

fn add_wave(
    sum: &mut ([C64; 3], [C64; 3]),
    k: [f64; 3],
    omega_mu: f64,
    tangential: [C64; 2],
    position: [f64; 3],
) {
    let e = [
        tangential[0],
        tangential[1],
        -(tangential[0].scale(k[0]) + tangential[1].scale(k[1])).scale(1.0 / k[2]),
    ];
    let h = [
        e[2].scale(k[1]) - e[1].scale(k[2]),
        e[0].scale(k[2]) - e[2].scale(k[0]),
        e[1].scale(k[0]) - e[0].scale(k[1]),
    ]
    .map(|v| v.scale(1.0 / omega_mu));
    let phasor = phase(-(k[0] * position[0] + k[1] * position[1] + k[2] * position[2]));
    for c in 0..3 {
        sum.0[c] = sum.0[c] + e[c] * phasor;
        sum.1[c] = sum.1[c] + h[c] * phasor;
    }
}

#[test]
fn g1_nonuniform_vector_probes_match_fresnel_maxwell_sums_and_original_fft_planes() {
    let incident = material(1.0, 1.4);
    let exit = material(2.25, 0.8);
    let omega = 1.3;
    let stack =
        PlanarDielectricLayers::try_new(&incident, DielectricResponse::Lossless(&exit), omega, &[])
            .unwrap();
    let grid = AngularSpectrumGrid {
        width: 8,
        height: 4,
        spacing_m: [7.0, 8.0],
    };
    let carrier = [0.17, -0.11];
    let authored = [
        (0, [C64::new(0.7, 0.2), C64::new(-0.3, 0.5)]),
        (1, [C64::new(-0.4, 0.1), C64::new(0.2, -0.6)]),
        (grid.width + 2, [C64::new(0.1, -0.3), C64::new(0.4, 0.2)]),
    ];
    let mut modes = vec![[C64::ZERO; 2]; grid.width * grid.height];
    for (index, amplitude) in authored {
        modes[index] = amplitude;
    }
    let beam = PlanarAngularSpectrum::prepare_modes_with_carrier(
        &stack,
        grid,
        carrier,
        &modes,
        limits(),
        || false,
    )
    .unwrap();
    // The transmitted receivers lie on an oblique plane; the same query also
    // observes the incident, reflected and standing fields at unrelated depths.
    let mut observations = Vec::new();
    for (x, y) in [(0.37, -1.4), (-4.1, 2.3), (3.2, 4.7), (11.3, -2.6)] {
        observations.push(AngularSpectrumObservation {
            side: Side::Transmitted,
            distance_from_face_m: 2.0 + 0.11 * x + 0.07 * y,
            position_m: [x, y],
        });
    }
    for side in [Side::Incident, Side::Reflected, Side::TotalIncidentSide] {
        observations.push(AngularSpectrumObservation {
            side,
            distance_from_face_m: 0.73,
            position_m: [-2.6, 1.81],
        });
    }
    let actual = beam
        .probe_fields(&observations, probe_limits(), || false)
        .unwrap();
    for (observation, actual) in observations.iter().zip(actual) {
        let mut sum = ([C64::ZERO; 3], [C64::ZERO; 3]);
        for (index, amplitude) in authored {
            let [kx, ky] = transverse(grid, carrier, index);
            let kt = det::sqrt(kx * kx + ky * ky);
            let p = [kx / kt, ky / kt];
            let a = [
                -amplitude[0].scale(p[1]) + amplitude[1].scale(p[0]),
                amplitude[0].scale(p[0]) + amplitude[1].scale(p[1]),
            ];
            let qi = det::sqrt(omega * omega * 1.4 - kt * kt);
            let qt = det::sqrt(omega * omega * 2.25 * 0.8 - kt * kt);
            let yi = [qi / (omega * 1.4), omega / qi];
            let yt = [qt / (omega * 0.8), omega * 2.25 / qt];
            let r = [0, 1].map(|c| (yi[c] - yt[c]) / (yi[c] + yt[c]));
            let t = [0, 1].map(|c| 2.0 * yi[c] / (yi[c] + yt[c]));
            let cartesian = |gain: [f64; 2]| {
                [
                    -a[0].scale(p[1] * gain[0]) + a[1].scale(p[0] * gain[1]),
                    a[0].scale(p[0] * gain[0]) + a[1].scale(p[1] * gain[1]),
                ]
            };
            let [x, y] = observation.position_m;
            let d = observation.distance_from_face_m;
            if matches!(observation.side, Side::Incident | Side::TotalIncidentSide) {
                add_wave(&mut sum, [kx, ky, qi], omega * 1.4, amplitude, [x, y, -d]);
            }
            if matches!(observation.side, Side::Reflected | Side::TotalIncidentSide) {
                add_wave(
                    &mut sum,
                    [kx, ky, -qi],
                    omega * 1.4,
                    cartesian(r),
                    [x, y, -d],
                );
            }
            if observation.side == Side::Transmitted {
                add_wave(&mut sum, [kx, ky, qt], omega * 0.8, cartesian(t), [x, y, d]);
            }
        }
        check_field(actual, sum.0, sum.1);
    }
    for side in [
        Side::Incident,
        Side::Reflected,
        Side::TotalIncidentSide,
        Side::Transmitted,
    ] {
        let plane = beam.field(side, 0.61, || false).unwrap();
        let observations: Vec<_> = (0..modes.len())
            .map(|index| AngularSpectrumObservation {
                side,
                distance_from_face_m: 0.61,
                position_m: [
                    (index % grid.width) as f64 * grid.spacing_m[0],
                    (index / grid.width) as f64 * grid.spacing_m[1],
                ],
            })
            .collect();
        let probes = beam
            .probe_fields(&observations, probe_limits(), || false)
            .unwrap();
        for (index, probe) in probes.into_iter().enumerate() {
            check_field(
                probe,
                plane.electric_v_per_m()[index],
                plane.magnetic_a_per_m()[index],
            );
        }
    }
}

#[test]
fn g1_evanescent_probes_keep_longitudinal_fields_reactive_flux_and_exact_decay() {
    let incident = material(4.0, 1.0);
    let exit = material(1.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&incident, DielectricResponse::Lossless(&exit), 1.0, &[])
            .unwrap();
    let grid = AngularSpectrumGrid {
        width: 1,
        height: 1,
        spacing_m: [1.0; 2],
    };
    let kx = 1.6;
    let input = [C64::new(0.7, -0.2), C64::ZERO];
    let beam = PlanarAngularSpectrum::prepare_modes_with_carrier(
        &stack,
        grid,
        [kx, 0.0],
        &[input],
        limits(),
        || false,
    )
    .unwrap();
    let kappa = det::sqrt(kx * kx - 1.0);
    let qt = C64::new(0.0, -kappa);
    let yi = 4.0 / det::sqrt(4.0 - kx * kx);
    let gain = C64::from_re(2.0 * yi) / (C64::from_re(yi) + C64::ONE / qt);
    let observations: Vec<_> = [0.0, 0.31, 1.2]
        .into_iter()
        .map(|d| AngularSpectrumObservation {
            side: Side::Transmitted,
            distance_from_face_m: d,
            position_m: [0.37 + 0.2 * d, -0.81],
        })
        .collect();
    let probes = beam
        .probe_fields(&observations, probe_limits(), || false)
        .unwrap();
    for (observation, probe) in observations.iter().zip(&probes) {
        let ex = input[0]
            * gain
            * phase(-kx * observation.position_m[0])
                .scale(det::exp(-kappa * observation.distance_from_face_m));
        check_field(
            *probe,
            [ex, C64::ZERO, -ex.scale(kx) / qt],
            [C64::ZERO, ex / qt, C64::ZERO],
        );
        assert!(probe.complex_poynting_w_m2[0].re > 0.0);
        assert!(probe.complex_poynting_w_m2[2].im < 0.0);
        near(probe.complex_poynting_w_m2[2].re, 0.0, 2e-14);
        relative(
            probe.complex_poynting_w_m2[2].im,
            probes[0].complex_poynting_w_m2[2].im
                * det::exp(-2.0 * kappa * observation.distance_from_face_m),
        );
    }
}

#[test]
fn g1_exact_critical_probes_preserve_finite_longitudinal_fields() {
    // This original Nyquist geometry admits cosine=3/4 exactly. Its exit
    // epsilon=7/16 gives kz=0, where reconstructing Ez by division would fail.
    let incident = material(1.0, 1.0);
    let exit = material(0.4375, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&incident, DielectricResponse::Lossless(&exit), 1.0, &[])
            .unwrap();
    assert_eq!(
        stack
            .response(0.75, Polarization::Tm)
            .unwrap()
            .exit_normal_wavenumber_rad_m,
        C64::ZERO
    );
    let grid = AngularSpectrumGrid {
        width: 4,
        height: 1,
        spacing_m: [core::f64::consts::PI / det::sqrt(0.4375), 1.0],
    };
    let amplitude = [C64::new(0.7, 0.2), C64::new(-0.3, 0.5)];
    let beam = PlanarAngularSpectrum::prepare_modes(
        &stack,
        grid,
        &[[C64::ZERO; 2], [C64::ZERO; 2], amplitude, [C64::ZERO; 2]],
        limits(),
        || false,
    )
    .unwrap();
    let observations = [0.0, 19.2].map(|d| AngularSpectrumObservation {
        side: Side::Transmitted,
        distance_from_face_m: d,
        position_m: [0.41, -0.7],
    });
    let probes = beam
        .probe_fields(&observations, probe_limits(), || false)
        .unwrap();
    let kx = -core::f64::consts::PI / grid.spacing_m[0];
    let weight = phase(-kx * observations[0].position_m[0]);
    let ey = amplitude[1].scale(2.0) * weight;
    let hy = amplitude[0].scale(2.0 / 0.75) * weight;
    for probe in probes {
        check_field(
            probe,
            [C64::ZERO, ey, hy.scale(-kx / 0.4375)],
            [C64::ZERO, hy, ey.scale(kx)],
        );
        assert!(probe.electric_v_per_m[2].norm_sq() > 0.0);
        near(probe.complex_poynting_w_m2[2].re, 0.0, 1e-14);
    }
}

fn film_matrix(q: C64, admittance: C64, depth: f64) -> [[C64; 2]; 2] {
    let z = q.scale(depth);
    let (s, c) = det::sin_cos(z.re);
    let ch = 0.5 * (det::exp(z.im) + det::exp(-z.im));
    let sh = 0.5 * (det::exp(z.im) - det::exp(-z.im));
    let cosine = C64::new(c * ch, -s * sh);
    let i_sine = C64::new(-c * sh, s * ch);
    [[cosine, i_sine / admittance], [i_sine * admittance, cosine]]
}

#[test]
fn g1_internal_probes_match_independent_film_matrices_interference_heat_and_traces() {
    let outside = material(E, M);
    let lossy = super::super::finite_layers::debye(3.0 * E, 4.0 * E, 0.03, 0.001);
    let omega = 1e9;
    let epsilon = C64::new(5.0 * E, -2.0 * E);
    let thickness = 0.27;
    let stack = PlanarDielectricLayers::try_new(
        &outside,
        DielectricResponse::Lossless(&outside),
        omega,
        &[DielectricLayer {
            response: DielectricResponse::Debye(&lossy),
            thickness_m: thickness,
        }],
    )
    .unwrap();
    let grid = AngularSpectrumGrid {
        width: 4,
        height: 2,
        spacing_m: [4.0, 6.0],
    };
    let carrier = [0.15, 0.12];
    let authored = [
        (0, [C64::new(0.7, 0.2), C64::new(-0.1, 0.6)]),
        (1, [C64::new(0.3, -0.4), C64::new(0.5, 0.2)]),
    ];
    let mut modes = vec![[C64::ZERO; 2]; grid.width * grid.height];
    for (index, value) in authored {
        modes[index] = value;
    }
    let beam = PlanarAngularSpectrum::prepare_modes_with_carrier(
        &stack,
        grid,
        carrier,
        &modes,
        limits(),
        || false,
    )
    .unwrap();
    let layer = beam
        .prepare_layer(
            0,
            AngularSpectrumLayerLimits {
                maximum_layer_mode_pairs: 1024,
                maximum_workspace_bytes: 8 << 20,
                maximum_fft_work: 1 << 20,
            },
            || false,
        )
        .unwrap();
    let observations = [
        AngularSpectrumLayerObservation {
            depth_m: 0.0,
            position_m: [0.37, -0.52],
        },
        AngularSpectrumLayerObservation {
            depth_m: 0.087,
            position_m: [-1.8, 2.4],
        },
        AngularSpectrumLayerObservation {
            depth_m: thickness,
            position_m: [3.1, -0.8],
        },
    ];
    let probes = layer
        .probe_fields(&observations, probe_limits(), || false)
        .unwrap();
    for (observation, probe) in observations.iter().zip(&probes) {
        let mut e = [C64::ZERO; 3];
        let mut h = [C64::ZERO; 3];
        let mut incoherent_heat = 0.0;
        for (index, input) in authored {
            let [kx, ky] = transverse(grid, carrier, index);
            let kt = det::sqrt(kx * kx + ky * ky);
            let p = [kx / kt, ky / kt];
            let amplitude = [
                -input[0].scale(p[1]) + input[1].scale(p[0]),
                input[0].scale(p[0]) + input[1].scale(p[1]),
            ];
            let qi = det::sqrt(omega * omega * E * M - kt * kt);
            let mut q = (epsilon.scale(omega * omega * M) - C64::from_re(kt * kt)).sqrt();
            if q.im > 0.0 {
                q = -q;
            }
            let y0 = [qi / (omega * M), omega * E / qi];
            let y = [q.scale(1.0 / (omega * M)), epsilon.scale(omega) / q];
            let mut et = [C64::ZERO; 2];
            let mut ht = [C64::ZERO; 2];
            for pol in 0..2 {
                let m = film_matrix(q, y[pol], thickness);
                let transmission = C64::from_re(2.0)
                    / (m[0][0] + m[0][1].scale(y0[pol]) + m[1][0].scale(1.0 / y0[pol]) + m[1][1]);
                let at = film_matrix(q, y[pol], thickness - observation.depth_m);
                et[pol] = (at[0][0] + at[0][1].scale(y0[pol])) * transmission * amplitude[pol];
                ht[pol] = (at[1][0] + at[1][1].scale(y0[pol])) * transmission * amplitude[pol];
            }
            let em = [
                -et[0].scale(p[1]) + et[1].scale(p[0]),
                et[0].scale(p[0]) + et[1].scale(p[1]),
                -ht[1].scale(kt) / epsilon.scale(omega),
            ];
            let hm = [
                -ht[0].scale(p[0]) - ht[1].scale(p[1]),
                -ht[0].scale(p[1]) + ht[1].scale(p[0]),
                et[0].scale(kt / (omega * M)),
            ];
            let weight = phase(-kx * observation.position_m[0] - ky * observation.position_m[1]);
            for c in 0..3 {
                e[c] = e[c] + em[c] * weight;
                h[c] = h[c] + hm[c] * weight;
            }
            incoherent_heat +=
                -0.5 * omega * epsilon.im * em.iter().map(|v| v.norm_sq()).sum::<f64>();
        }
        check_field(probe.field, e, h);
        let expected_heat = -0.5 * omega * epsilon.im * e.iter().map(|v| v.norm_sq()).sum::<f64>();
        relative(probe.absorbed_power_density_w_m3, expected_heat);
        assert!(probe.absorbed_power_density_w_m3 > 0.0);
        assert!((expected_heat - incoherent_heat).abs() > 0.02 * expected_heat);
    }
    for (index, side) in [(0, Side::TotalIncidentSide), (2, Side::Transmitted)] {
        let outside = beam
            .probe_fields(
                &[AngularSpectrumObservation {
                    side,
                    distance_from_face_m: 0.0,
                    position_m: observations[index].position_m,
                }],
                probe_limits(),
                || false,
            )
            .unwrap()[0];
        let inside = probes[index].field;
        for c in 0..2 {
            complex_near(
                inside.electric_v_per_m[c],
                outside.electric_v_per_m[c],
                3e-12,
            );
            complex_near(
                inside.magnetic_a_per_m[c],
                outside.magnetic_a_per_m[c],
                3e-14,
            );
        }
        complex_near(
            inside.electric_v_per_m[2] * epsilon,
            outside.electric_v_per_m[2].scale(E),
            3e-22,
        );
        complex_near(
            inside.magnetic_a_per_m[2],
            outside.magnetic_a_per_m[2],
            3e-14,
        );
    }
    let depth = 0.113;
    let plane = layer.field(depth, || false).unwrap();
    let observations: Vec<_> = (0..modes.len())
        .map(|index| AngularSpectrumLayerObservation {
            depth_m: depth,
            position_m: [
                (index % grid.width) as f64 * grid.spacing_m[0],
                (index / grid.width) as f64 * grid.spacing_m[1],
            ],
        })
        .collect();
    let probes = layer
        .probe_fields(&observations, probe_limits(), || false)
        .unwrap();
    for (index, probe) in probes.into_iter().enumerate() {
        check_field(
            probe.field,
            plane.field().electric_v_per_m()[index],
            plane.field().magnetic_a_per_m()[index],
        );
        relative(
            probe.absorbed_power_density_w_m3,
            plane.absorbed_power_density_w_m3()[index],
        );
    }
    for depth in [-0.1, thickness + 0.1, f64::NAN] {
        assert!(
            layer
                .probe_fields(
                    &[AngularSpectrumLayerObservation {
                        depth_m: depth,
                        position_m: [0.0; 2],
                    }],
                    probe_limits(),
                    || false
                )
                .is_err()
        );
    }
}

#[test]
fn g0_g4_sparse_probe_budgets_cancellation_and_zero_fields_are_atomic() {
    let material = material(1.0, 1.0);
    let stack = PlanarDielectricLayers::try_new(
        &material,
        DielectricResponse::Lossless(&material),
        1.0,
        &[],
    )
    .unwrap();
    let grid = AngularSpectrumGrid {
        width: 32,
        height: 1,
        spacing_m: [1.0; 2],
    };
    let mut input = vec![[C64::ZERO; 2]; 32];
    input[0][1] = C64::ONE;
    input[1][0] = C64::new(0.2, 0.3);
    let beam = PlanarAngularSpectrum::prepare_modes_with_carrier(
        &stack,
        grid,
        [0.17, 0.0],
        &input,
        limits(),
        || false,
    )
    .unwrap();
    let observations = [AngularSpectrumObservation {
        side: Side::Transmitted,
        distance_from_face_m: 0.7,
        position_m: [0.31, -0.2],
    }; 2];
    let exact = AngularSpectrumProbeLimits {
        maximum_points: 2,
        maximum_mode_evaluations: 4,
        ..probe_limits()
    };
    let before = beam.probe_fields(&observations, exact, || false).unwrap();
    for limits in [
        AngularSpectrumProbeLimits {
            maximum_points: 1,
            ..exact
        },
        AngularSpectrumProbeLimits {
            maximum_mode_evaluations: 3,
            ..exact
        },
        AngularSpectrumProbeLimits {
            maximum_workspace_bytes: 0,
            ..exact
        },
    ] {
        assert!(beam.probe_fields(&observations, limits, || false).is_err());
    }
    let mut polls = 0;
    assert!(matches!(
        beam.probe_fields(&observations, exact, || {
            polls += 1;
            polls == 8
        }),
        Err(AngularSpectrumError::Cancelled)
    ));
    for observation in [
        AngularSpectrumObservation {
            distance_from_face_m: -0.1,
            ..observations[0]
        },
        AngularSpectrumObservation {
            position_m: [f64::NAN, 0.0],
            ..observations[0]
        },
    ] {
        assert!(beam.probe_fields(&[observation], exact, || false).is_err());
    }
    assert_eq!(
        beam.probe_fields(&observations, exact, || false).unwrap(),
        before
    );
    let zero = PlanarAngularSpectrum::prepare_modes_with_carrier(
        &stack,
        grid,
        [0.17, 0.0],
        &[[C64::ZERO; 2]; 32],
        limits(),
        || false,
    )
    .unwrap();
    let fields = zero
        .probe_fields(
            &[AngularSpectrumObservation {
                side: Side::Transmitted,
                distance_from_face_m: f64::MAX,
                position_m: [f64::MAX; 2],
            }],
            AngularSpectrumProbeLimits {
                maximum_mode_evaluations: 0,
                ..exact
            },
            || false,
        )
        .unwrap();
    assert_eq!(
        fields[0],
        AngularSpectrumProbeField {
            electric_v_per_m: [C64::ZERO; 3],
            magnetic_a_per_m: [C64::ZERO; 3],
            complex_poynting_w_m2: [C64::ZERO; 3]
        }
    );
}
