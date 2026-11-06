use super::*;
use crate::resistive_film::shear::WallVelocityPorts;
use crate::resistive_film::{FilmCell, FilmChannel, FilmLimits, GapPort, ResistiveFilm};

const H: f64 = 0.0002;
const W: f64 = 0.01;
const L: f64 = 0.02;
const MU: f64 = 0.01;
const AREA: f64 = 0.0001;

fn near(a: f64, b: f64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * a.abs().max(b.abs()).max(1e-16),
        "{a:e} != {b:e}"
    );
}

fn base(alpha: f64, moving: bool) -> ShearedFilm {
    let gap = |reference_m| GapPort {
        reference_m,
        closure: vec![1.0, 0.0, 0.0],
    };
    let film = ResistiveFilm::new(
        vec![FilmCell {
            area_m2: AREA,
            gap: gap(H),
        }],
        (0..2)
            .map(|_| FilmChannel {
                from: 0,
                to: None,
                width_m: W,
                length_m: L,
                gap: gap(H),
            })
            .collect(),
        3,
        MU,
        FilmLimits {
            minimum_cell_gap_m: H / 10.0,
            maximum_gap_m: 2.0 * H,
            maximum_pressure_pa: 1e4,
        },
    )
    .unwrap()
    .with_barus_viscosity(alpha)
    .unwrap();
    let walls = [-1.0, 1.0].map(|sign| WallVelocityPorts {
        lower: vec![0.0, if moving { sign } else { 0.0 }, 0.0],
        upper: vec![0.0, 0.0, if moving { sign } else { 0.0 }],
    });
    ShearedFilm::new(film, walls.into()).unwrap()
}

fn response(
    model: &NavierSlipFilm,
    q: [f64; 3],
    velocity: [f64; 3],
) -> (PressureDrivenReport, [f64; 3], f64, [SlipChannelReport; 2]) {
    let mut force = [0.0; 3];
    let mut pressure = [0.0];
    let mut channels = [SlipChannelReport::default(); 2];
    let report = model
        .evaluate_into(&q, &velocity, &mut force, &mut pressure, &mut channels)
        .unwrap();
    (report, force, pressure[0], channels)
}

#[test]
fn zero_slip_recovers_original_pressure_wall_force_heat_and_tangent() {
    for alpha in [0.0, 3e-4] {
        let original = base(alpha, true);
        let native = super::super::PressureDrivenShearedFilm::new(
            original.clone(),
            &[Some(1500.0), Some(-500.0)],
        )
        .unwrap();
        let model = NavierSlipFilm::new(original, vec![WallSlipLengths::default(); 2])
            .unwrap()
            .with_reservoir_pressures(&[Some(1500.0), Some(-500.0)])
            .unwrap();
        let q = [H / 7.0, 0.0, 0.0];
        let velocity = [0.00001, 0.02, -0.006];
        let (report, actual, pressure, channels) = response(&model, q, velocity);
        let mut fast_force = [0.0; 3];
        let mut fast_pressure = [0.0];
        let fast = model
            .evaluate_force_into(&q, &velocity, &mut fast_force, &mut fast_pressure)
            .unwrap();
        assert_eq!(fast_force, actual);
        assert_eq!(fast_pressure, [pressure]);
        assert_eq!(fast.film.dissipated_power_w, report.film.dissipated_power_w);
        assert_eq!(model.sheared_film().pressure_operator().port_count(), 3);
        let mut expected = [0.0; 3];
        let mut expected_pressure = [0.0];
        let mut expected_channels = [ShearChannelReport::default(); 2];
        let reference = native
            .evaluate_into(
                &q,
                &velocity,
                &mut expected,
                &mut expected_pressure,
                &mut expected_channels,
            )
            .unwrap();
        for i in 0..3 {
            near(actual[i], expected[i], 5e-13);
        }
        near(pressure, expected_pressure[0], 5e-13);
        near(
            report.film.dissipated_power_w,
            reference.film.dissipated_power_w,
            5e-13,
        );
        for k in 0..2 {
            near(
                channels[k].film.volume_flow_m3_s,
                expected_channels[k].volume_flow_m3_s,
                5e-13,
            );
            near(
                channels[k].bulk_dissipated_power_w,
                expected_channels[k].dissipated_power_w(),
                5e-13,
            );
            assert_eq!(channels[k].lower_slip_dissipated_power_w, 0.0);
            assert_eq!(channels[k].upper_slip_dissipated_power_w, 0.0);
        }
        let dq = [0.1 * H, 0.0, 0.0];
        let dv = [-0.000003, 0.001, 0.002];
        model
            .tangent_into(&q, &velocity, &dq, &dv, &mut expected)
            .unwrap();
        let mut reference = [0.0; 3];
        native
            .tangent_into(&q, &velocity, &dq, &dv, &[0.0; 2], &mut reference)
            .unwrap();
        for i in 0..3 {
            near(expected[i], reference[i], 5e-12);
        }
    }
}

#[test]
fn unequal_slip_couette_drag_and_three_heat_shares_match_closed_form() {
    let slip = WallSlipLengths {
        lower_m: H / 4.0,
        upper_m: H,
    };
    let separation = H + slip.lower_m + slip.upper_m;
    for alpha in [0.0, 3e-4] {
        let model = NavierSlipFilm::new(base(alpha, true), vec![slip; 2])
            .unwrap()
            .with_reservoir_pressures(&[Some(2000.0); 2])
            .unwrap();
        let mu = MU * (alpha * 2000.0).exp();
        for velocity in [
            [0.0, 0.002, 0.006],
            [0.0, 0.006, 0.002],
            [0.0, 0.004, 0.004],
        ] {
            let (report, _, pressure, channels) = response(&model, [0.0; 3], velocity);
            near(pressure, 2000.0, 5e-14);
            let du = velocity[2] - velocity[1];
            let heat_per_length = W * L * mu * (du / separation).powi(2);
            for (channel, sign) in channels.into_iter().zip([-1.0, 1.0]) {
                near(
                    channel.film.upper_wall_traction_pa,
                    sign * mu * du / separation,
                    5e-13,
                );
                near(
                    channel.film.lower_wall_traction_pa,
                    -sign * mu * du / separation,
                    5e-13,
                );
                near(
                    channel.lower_slip_velocity_m_s,
                    sign * slip.lower_m * du / separation,
                    5e-13,
                );
                near(
                    channel.upper_slip_velocity_m_s,
                    -sign * slip.upper_m * du / separation,
                    5e-13,
                );
                near(channel.bulk_dissipated_power_w, heat_per_length * H, 5e-13);
                near(
                    channel.lower_slip_dissipated_power_w,
                    heat_per_length * slip.lower_m,
                    5e-13,
                );
                near(
                    channel.upper_slip_dissipated_power_w,
                    heat_per_length * slip.upper_m,
                    5e-13,
                );
            }
            near(
                report.film.dissipated_power_w,
                2.0 * heat_per_length * separation,
                5e-13,
            );
        }
    }
}

#[test]
fn pressure_mobility_has_exact_equal_and_single_wall_slip_limits() {
    for slip in [
        WallSlipLengths {
            lower_m: H,
            upper_m: H,
        },
        WallSlipLengths {
            lower_m: H,
            upper_m: 0.0,
        },
        WallSlipLengths {
            lower_m: 0.0,
            upper_m: 1e6 * H,
        },
    ] {
        let model = NavierSlipFilm::new(base(0.0, false), vec![slip; 2])
            .unwrap()
            .with_reservoir_pressures(&[Some(1500.0), Some(-500.0)])
            .unwrap();
        let (_, _, _, channels) = response(&model, [0.0; 3], [0.0; 3]);
        let factor = if slip.lower_m == slip.upper_m {
            1.0 + 6.0 * slip.lower_m / H
        } else {
            let b = slip.lower_m + slip.upper_m;
            (H + 4.0 * b) / (H + b)
        };
        near(
            channels[1].film.volume_flow_m3_s,
            W * H.powi(3) * 1000.0 / (12.0 * MU * L) * factor,
            5e-13,
        );
        assert!(channels[1].film.pressure_dissipated_power_w > 0.0);
        assert_eq!(channels[1].film.wall_power_w, 0.0);
    }
}

#[test]
fn mixed_pressure_and_shear_match_independent_local_velocity_and_heat_integrals() {
    let slip = WallSlipLengths {
        lower_m: 0.2 * H,
        upper_m: 0.8 * H,
    };
    for alpha in [0.0, 3e-4] {
        let model = NavierSlipFilm::new(base(alpha, true), vec![slip; 2])
            .unwrap()
            .with_reservoir_pressures(&[Some(1500.0), Some(-500.0)])
            .unwrap();
        let velocity = [0.0, 0.02, 0.06];
        let (_, force, pressure, channels) = response(&model, [0.0; 3], velocity);
        for (k, sign, boundary) in [(0, -1.0, 1500.0), (1, 1.0, -500.0)] {
            let potential = |p: f64| {
                if alpha == 0.0 {
                    p
                } else {
                    -(-alpha * p).exp_m1() / alpha
                }
            };
            let gradient = (potential(pressure) - potential(boundary)) / L;
            let c1 = (sign * (velocity[2] - velocity[1])
                + gradient * H * (H / 2.0 + slip.upper_m) / MU)
                / (H + slip.lower_m + slip.upper_m);
            let c0 = sign * velocity[1] + slip.lower_m * c1;
            let speed = |y: f64| c0 + c1 * y - gradient * y * y / (2.0 * MU);
            let flow = W * H / 6.0 * (speed(0.0) + 4.0 * speed(H / 2.0) + speed(H));
            near(channels[k].film.volume_flow_m3_s, flow, 5e-12);
            let mut mu_average = 0.0;
            for i in 0..=256 {
                let x = i as f64 / 256.0;
                let fluidity =
                    (1.0 - x) * (-alpha * pressure).exp() + x * (-alpha * boundary).exp();
                let weight = if i == 0 || i == 256 {
                    1.0
                } else if i % 2 == 0 {
                    2.0
                } else {
                    4.0
                };
                mu_average += weight * MU / fluidity / (3.0 * 256.0);
            }
            let shear = |y: f64| c1 - gradient * y / MU;
            let bulk = W * L * H * mu_average / 6.0
                * (shear(0.0).powi(2) + 4.0 * shear(H / 2.0).powi(2) + shear(H).powi(2));
            near(channels[k].bulk_dissipated_power_w, bulk, 3e-10);
            near(
                channels[k].lower_slip_dissipated_power_w,
                W * L * mu_average * slip.lower_m * shear(0.0).powi(2),
                3e-10,
            );
            near(
                channels[k].upper_slip_dissipated_power_w,
                W * L * mu_average * slip.upper_m * shear(H).powi(2),
                3e-10,
            );
            near(
                channels[k].film.lower_wall_traction_pa,
                -mu_average * shear(0.0),
                3e-10,
            );
            near(
                channels[k].film.upper_wall_traction_pa,
                mu_average * shear(H),
                3e-10,
            );
            near(
                channels[k].film.volume_flow_m3_s * (pressure - boundary)
                    + channels[k].film.wall_power_w,
                channels[k].film.dissipated_power_w(),
                5e-12,
            );
        }
        near(
            channels.iter().map(|c| c.film.wall_power_w).sum(),
            force[1] * velocity[1] + force[2] * velocity[2],
            5e-12,
        );
    }
}

#[test]
fn asymmetric_slip_tangent_differentiates_gap_squeeze_and_both_wall_rates() {
    let slips = vec![
        WallSlipLengths {
            lower_m: 0.3 * H,
            upper_m: 1.1 * H,
        },
        WallSlipLengths {
            lower_m: H,
            upper_m: 0.1 * H,
        },
    ];
    for alpha in [0.0, 3e-4] {
        let model = NavierSlipFilm::new(base(alpha, true), slips.clone())
            .unwrap()
            .with_reservoir_pressures(&[Some(1500.0), Some(-500.0)])
            .unwrap();
        let q = [H / 5.0, 0.0, 0.0];
        let velocity = [0.00002, 0.02, -0.03];
        let dq = [H / 8.0, 0.0, 0.0];
        let dv = [-0.00003, 0.004, 0.002];
        let mut tangent = [0.0; 3];
        model
            .tangent_into(&q, &velocity, &dq, &dv, &mut tangent)
            .unwrap();
        let eps = 1e-4;
        let plus = response(
            &model,
            std::array::from_fn(|i| q[i] + eps * dq[i]),
            std::array::from_fn(|i| velocity[i] + eps * dv[i]),
        )
        .1;
        let minus = response(
            &model,
            std::array::from_fn(|i| q[i] - eps * dq[i]),
            std::array::from_fn(|i| velocity[i] - eps * dv[i]),
        )
        .1;
        for i in 0..3 {
            near(tangent[i], (plus[i] - minus[i]) / (2.0 * eps), 2e-7);
        }
        let updated = model.with_channel_viscosities(&[2.0 * MU; 2]).unwrap();
        assert_eq!(updated.wall_slip_lengths(), &slips);
        assert!(response(&updated, q, velocity).0.film.dissipated_power_w > 0.0);
    }
}

#[test]
fn invalid_slip_and_collapsed_moving_channels_refuse_without_partial_outputs() {
    for value in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            NavierSlipFilm::new(
                base(0.0, true),
                vec![
                    WallSlipLengths {
                        lower_m: value,
                        upper_m: 0.0
                    };
                    2
                ]
            )
            .is_err()
        );
    }
    let model = NavierSlipFilm::new(
        base(0.0, true),
        vec![
            WallSlipLengths {
                lower_m: H,
                upper_m: H
            };
            2
        ],
    )
    .unwrap();
    let mut force = [42.0; 3];
    let mut pressure = [42.0];
    let sentinel = SlipChannelReport {
        bulk_dissipated_power_w: 42.0,
        ..SlipChannelReport::default()
    };
    let mut channels = [sentinel; 2];
    assert!(
        model
            .evaluate_into(
                &[H, 0.0, 0.0],
                &[0.0, 0.0, 1.0],
                &mut force,
                &mut pressure,
                &mut channels
            )
            .is_err()
    );
    assert_eq!(force, [42.0; 3]);
    assert_eq!(pressure, [42.0]);
    assert_eq!(channels, [sentinel; 2]);
}
