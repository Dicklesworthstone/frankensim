use super::*;
use crate::resistive_film::stratified::ViscosityLayer;
use crate::resistive_film::{FilmCell, FilmChannel, FilmChannelReport, FilmLimits, GapPort};

const H: f64 = 0.0002;
const W: f64 = 0.01;
const L: f64 = 0.02;
const AREA: f64 = 0.0001;
const MU: f64 = 0.01;

mod barus;

fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}

fn film(gaps: &[f64]) -> ResistiveFilm {
    let gap = |reference_m| GapPort {
        reference_m,
        closure: vec![1.0, 0.0, 0.0],
    };
    ResistiveFilm::new(
        vec![FilmCell {
            area_m2: AREA,
            gap: gap(H),
        }],
        gaps.iter()
            .map(|&h| FilmChannel {
                from: 0,
                to: None,
                width_m: W,
                length_m: L,
                gap: gap(h),
            })
            .collect(),
        3,
        MU,
        FilmLimits {
            minimum_cell_gap_m: H / 10.0,
            maximum_gap_m: 2.0 * H,
            maximum_pressure_pa: 10000.0,
        },
    )
    .unwrap()
}

fn walls(sign: f64) -> WallVelocityPorts {
    WallVelocityPorts {
        lower: vec![0.0, sign, 0.0],
        upper: vec![0.0, 0.0, sign],
    }
}

fn layers() -> GapViscosityProfile {
    GapViscosityProfile::new(vec![
        ViscosityLayer {
            upper_fraction: 0.5,
            dynamic_viscosity_pa_s: 0.01,
        },
        ViscosityLayer {
            upper_fraction: 1.0,
            dynamic_viscosity_pa_s: 0.04,
        },
    ])
    .unwrap()
}

#[test]
fn zero_wall_rows_recover_native_force_work_tangent_and_closed_passages() {
    // G3: zero wall rows are an exact native limit, including a shut passage.
    for gaps in [[H, H], [H, 0.0]] {
        for stratified in [false, true] {
            let original = film(&gaps);
            let fixed = vec![walls(0.0); 2];
            let (native, moving) = if stratified {
                let binding = StratifiedFilm::new(&original, vec![layers(), layers()]).unwrap();
                let native = binding.film().clone();
                (
                    native,
                    ShearedFilm::from_stratified(binding, fixed.clone()).unwrap(),
                )
            } else {
                (
                    original.clone(),
                    ShearedFilm::new(original, fixed.clone()).unwrap(),
                )
            };
            assert_eq!(moving.wall_ports(), fixed);
            assert_eq!(moving.stratified_profiles().is_some(), stratified);
            let q = [0.0; 3];
            let v = [0.001, 0.013, -0.02];
            let mut expected_force = [0.0; 3];
            let mut force = [0.0; 3];
            let mut expected_pressure = [0.0];
            let mut pressure = [0.0];
            let mut expected_channels = [FilmChannelReport::default(); 2];
            let mut channels = [ShearChannelReport::default(); 2];
            let expected = native
                .evaluate_with_channels_into(
                    &q,
                    &v,
                    &mut expected_force,
                    &mut expected_pressure,
                    &mut expected_channels,
                )
                .unwrap();
            let actual = moving
                .evaluate_with_channels_into(&q, &v, &mut force, &mut pressure, &mut channels)
                .unwrap();
            assert_eq!(force, expected_force);
            assert_eq!(pressure, expected_pressure);
            assert_eq!(actual.dissipated_power_w, expected.dissipated_power_w);
            for i in 0..2 {
                assert_eq!(
                    channels[i].volume_flow_m3_s,
                    expected_channels[i].volume_flow_m3_s
                );
                assert_eq!(
                    channels[i].dissipated_power_w(),
                    expected_channels[i].dissipated_power_w
                );
                assert_eq!(channels[i].entrained_flow_m3_s, 0.0);
                assert_eq!(channels[i].shear_dissipated_power_w, 0.0);
            }
            let dq = [0.00001, 0.0, 0.0];
            let dv = [0.0001, 0.004, -0.002];
            native
                .tangent_into(&q, &v, &dq, &dv, &mut expected_force)
                .unwrap();
            moving.tangent_into(&q, &v, &dq, &dv, &mut force).unwrap();
            assert_eq!(force, expected_force);
        }
    }
}

#[test]
fn open_couette_flow_has_exact_drag_entrainment_and_rigid_translation_limit() {
    // G2: two oriented halves connect the same ambient pressure through a cell.
    let moving = ShearedFilm::new(film(&[H, H]), vec![walls(1.0), walls(-1.0)]).unwrap();
    let drag = MU * W * L / H;
    for (lower, upper) in [(0.0005, 0.002), (0.002, 0.002)] {
        let mut force = [0.0; 3];
        let mut pressure = [0.0];
        let mut channels = [ShearChannelReport::default(); 2];
        let report = moving
            .evaluate_with_channels_into(
                &[0.0; 3],
                &[0.0, lower, upper],
                &mut force,
                &mut pressure,
                &mut channels,
            )
            .unwrap();
        let relative = upper - lower;
        let flow = W * H * 0.5 * (lower + upper);
        assert_eq!(pressure, [0.0]);
        assert_eq!(force[0], 0.0);
        near(force[1], -2.0 * drag * relative, 1e-18);
        near(force[2], 2.0 * drag * relative, 1e-18);
        near(channels[0].volume_flow_m3_s, flow, 1e-23);
        near(channels[1].volume_flow_m3_s, -flow, 1e-23);
        near(channels[0].entrained_flow_m3_s, flow, 1e-23);
        near(
            channels[0].dissipated_power_w(),
            drag * relative * relative,
            1e-20,
        );
        near(
            report.dissipated_power_w,
            2.0 * drag * relative * relative,
            1e-20,
        );
        assert_eq!(channels[0].pressure_dissipated_power_w, 0.0);
        if lower == upper {
            assert_eq!(report.dissipated_power_w, 0.0);
            assert_eq!(force, [0.0; 3]);
            assert!(channels[0].volume_flow_m3_s > 0.0);
        }
    }
}

#[test]
fn blocked_entrainment_generates_back_pressure_and_complete_wall_work() {
    // G2: zero net flow requires p=-6*mu*L*U/h^2 and upper drag=4*mu*w*L*U/h.
    // Q_total*delta_p is zero, although both shear and pressure dissipation are positive.
    let moving = ShearedFilm::new(film(&[H]), vec![walls(1.0)]).unwrap();
    let speed = 0.002;
    let drag = MU * W * L / H;
    let mut force = [0.0; 3];
    let mut pressure = [0.0];
    let mut channels = [ShearChannelReport::default()];
    let report = moving
        .evaluate_with_channels_into(
            &[0.0; 3],
            &[0.0, 0.0, speed],
            &mut force,
            &mut pressure,
            &mut channels,
        )
        .unwrap();
    let exact_pressure = -6.0 * MU * L * speed / (H * H);
    near(pressure[0], exact_pressure, 1e-12);
    near(force[0], AREA * exact_pressure, 1e-15);
    near(force[1], 2.0 * drag * speed, 1e-18);
    near(force[2], 4.0 * drag * speed, 1e-18);
    near(channels[0].volume_flow_m3_s, 0.0, 1e-22);
    near(channels[0].entrained_flow_m3_s, W * H * speed / 2.0, 1e-23);
    near(
        channels[0].pressure_dissipated_power_w,
        3.0 * drag * speed * speed,
        1e-20,
    );
    near(
        channels[0].shear_dissipated_power_w,
        drag * speed * speed,
        1e-20,
    );
    near(report.dissipated_power_w, 4.0 * drag * speed * speed, 1e-20);
    near(report.dissipated_power_w, force[2] * speed, 1e-20);
}

#[test]
fn stratified_couette_poiseuille_flow_matches_independent_velocity_and_shear_quadrature() {
    // G2: asymmetric layers distinguish hydraulic viscosity, harmonic shear
    // viscosity and wall entrainment. Integrate the physical velocity/stress.
    let original = layers();
    let binding = StratifiedFilm::new(&film(&[H]), vec![original.clone()]).unwrap();
    let moving = ShearedFilm::from_stratified(binding, vec![walls(1.0)]).unwrap();
    assert_eq!(moving.stratified_profiles().unwrap(), &[original.clone()]);
    let velocity = [0.0001, 0.001, 0.004];
    let mut force = [0.0; 3];
    let mut pressure = [0.0];
    let mut channels = [ShearChannelReport::default()];
    let report = moving
        .evaluate_with_channels_into(
            &[0.0; 3],
            &velocity,
            &mut force,
            &mut pressure,
            &mut channels,
        )
        .unwrap();
    let neutral = 7.0 / 20.0;
    let fluidity = 0.5 / 0.01 + 0.5 / 0.04;
    let hydraulic = 0.01 * 160.0 / 73.0;
    let shear_viscosity = 1.0 / fluidity;
    assert_ne!(hydraulic, shear_viscosity);
    let conductance = W * H.powi(3) / (12.0 * hydraulic * L);
    let entrained = W * H * (neutral * velocity[1] + (1.0 - neutral) * velocity[2]);
    let exact_pressure = (AREA * velocity[0] - entrained) / conductance;
    near(pressure[0], exact_pressure, 1e-11);
    near(channels[0].entrained_flow_m3_s, entrained, 1e-22);
    let relative = velocity[2] - velocity[1];
    let stress = |z: f64| -exact_pressure * H / L * (z - neutral) + relative / (H * fluidity);
    let mut lo = 0.0;
    let mut lower_velocity = velocity[1];
    let mut flux_integral = 0.0;
    let mut work_integral = 0.0;
    for layer in original.layers() {
        let hi = layer.upper_fraction;
        let half = 0.5 * (hi - lo);
        let points = [
            lo + half - half / 3.0_f64.sqrt(),
            lo + half + half / 3.0_f64.sqrt(),
        ];
        let u = |z: f64| {
            lower_velocity
                + (-exact_pressure * H * H / L * (0.5 * (z * z - lo * lo) - neutral * (z - lo))
                    + relative / fluidity * (z - lo))
                    / layer.dynamic_viscosity_pa_s
        };
        flux_integral += half * (u(points[0]) + u(points[1]));
        work_integral += half * (stress(points[0]).powi(2) + stress(points[1]).powi(2))
            / layer.dynamic_viscosity_pa_s;
        lower_velocity = u(hi);
        lo = hi;
    }
    near(lower_velocity, velocity[2], 1e-17);
    near(channels[0].volume_flow_m3_s, W * H * flux_integral, 1e-22);
    near(channels[0].volume_flow_m3_s, AREA * velocity[0], 1e-22);
    near(force[0], AREA * exact_pressure, 1e-15);
    near(force[1], -W * L * stress(0.0), 1e-18);
    near(force[2], W * L * stress(1.0), 1e-18);
    near(report.dissipated_power_w, W * L * H * work_integral, 1e-18);
    near(report.dissipated_power_w, dot(&velocity, &force), 1e-18);
    near(
        channels[0].shear_dissipated_power_w,
        W * L * shear_viscosity / H * relative * relative,
        1e-20,
    );
}

fn network() -> ShearedFilm {
    let gap = |reference_m, a, b| GapPort {
        reference_m,
        closure: vec![a, b, 0.0, 0.0],
    };
    let base = ResistiveFilm::new(
        vec![
            FilmCell {
                area_m2: AREA,
                gap: gap(H, 1.0, 0.0),
            },
            FilmCell {
                area_m2: 0.75 * AREA,
                gap: gap(1.1 * H, 0.0, 1.0),
            },
        ],
        vec![
            FilmChannel {
                from: 0,
                to: None,
                width_m: W,
                length_m: 0.015,
                gap: gap(H, 1.0, 0.0),
            },
            FilmChannel {
                from: 0,
                to: Some(1),
                width_m: 0.008,
                length_m: L,
                gap: gap(H, 0.4, 0.6),
            },
            FilmChannel {
                from: 1,
                to: None,
                width_m: 0.013,
                length_m: 0.012,
                gap: gap(1.1 * H, 0.0, 1.0),
            },
        ],
        4,
        MU,
        FilmLimits {
            minimum_cell_gap_m: H / 10.0,
            maximum_gap_m: 3.0 * H,
            maximum_pressure_pa: 10000.0,
        },
    )
    .unwrap();
    let profiles = vec![
        layers(),
        layers(),
        GapViscosityProfile::new(vec![ViscosityLayer {
            upper_fraction: 1.0,
            dynamic_viscosity_pa_s: 0.006,
        }])
        .unwrap(),
    ];
    let binding = StratifiedFilm::new(&base, profiles).unwrap();
    let ports = [-1.0, 1.0, 1.0]
        .into_iter()
        .map(|sign| WallVelocityPorts {
            lower: vec![0.0, 0.0, sign, 0.0],
            upper: vec![0.0, 0.0, 0.0, sign],
        })
        .collect();
    ShearedFilm::from_stratified(binding, ports).unwrap()
}

#[test]
fn complete_network_tangent_is_reciprocal_and_differentiates_gap_and_wall_motion() {
    // G0/G3: channel continuity, complete port work, Onsager symmetry and a
    // directional derivative across the shared multi-cell pressure factorization.
    let moving = network();
    let q = [0.00001, -0.00002, 0.0, 0.0];
    let velocity = [0.0003, -0.0002, 0.001, 0.006];
    let mut force = [0.0; 4];
    let mut pressure = [0.0; 2];
    let mut channels = [ShearChannelReport::default(); 3];
    let report = moving
        .evaluate_with_channels_into(&q, &velocity, &mut force, &mut pressure, &mut channels)
        .unwrap();
    near(
        channels[0].volume_flow_m3_s + channels[1].volume_flow_m3_s,
        AREA * velocity[0],
        1e-21,
    );
    near(
        channels[2].volume_flow_m3_s - channels[1].volume_flow_m3_s,
        0.75 * AREA * velocity[1],
        1e-21,
    );
    near(dot(&velocity, &force), report.dissipated_power_w, 1e-16);
    near(
        channels
            .iter()
            .map(ShearChannelReport::dissipated_power_w)
            .sum(),
        report.dissipated_power_w,
        1e-16,
    );
    assert!(report.dissipated_power_w > 0.0);
    let other = [-0.0002, 0.00015, -0.002, 0.005];
    let mut other_force = [0.0; 4];
    moving
        .evaluate_into(&q, &other, &mut other_force, &mut pressure)
        .unwrap();
    near(dot(&velocity, &other_force), dot(&other, &force), 1e-16);
    assert!(dot(&other, &other_force) > 0.0);
    let dq = [0.00002, -0.00001, 0.0, 0.0];
    let dv = [0.0001, -0.0002, 0.003, -0.002];
    let mut tangent = [0.0; 4];
    moving
        .tangent_into(&q, &velocity, &dq, &dv, &mut tangent)
        .unwrap();
    let epsilon = 1e-5;
    let mut forces = [[0.0; 4]; 2];
    for (i, sign) in [-1.0, 1.0].into_iter().enumerate() {
        let position = std::array::from_fn::<_, 4, _>(|j| q[j] + sign * epsilon * dq[j]);
        let rate = std::array::from_fn::<_, 4, _>(|j| velocity[j] + sign * epsilon * dv[j]);
        moving
            .evaluate_into(&position, &rate, &mut forces[i], &mut pressure)
            .unwrap();
    }
    for i in 0..4 {
        let finite_difference = (forces[1][i] - forces[0][i]) / (2.0 * epsilon);
        near(
            tangent[i],
            finite_difference,
            1e-8 * tangent[i].abs() + 1e-10,
        );
    }
}

#[test]
fn invalid_rows_collapsed_sliding_and_late_refusals_leave_outputs_untouched() {
    let base = film(&[H, H]);
    assert!(ShearedFilm::new(base.clone(), vec![walls(1.0)]).is_err());
    let mut invalid = vec![walls(1.0), walls(-1.0)];
    invalid[0].lower.pop();
    assert!(ShearedFilm::new(base.clone(), invalid).is_err());
    let mut invalid = vec![walls(1.0), walls(-1.0)];
    invalid[1].upper[2] = f64::NAN;
    assert!(ShearedFilm::new(base.clone(), invalid).is_err());
    let mut unrepresentable = vec![walls(1.0), walls(-1.0)];
    unrepresentable[0].lower[1] = -f64::MAX;
    unrepresentable[0].upper[1] = f64::MAX;
    assert!(ShearedFilm::new(base.clone(), unrepresentable).is_err());
    let moving = ShearedFilm::new(base, vec![walls(1.0); 2]).unwrap();
    let sentinel = ShearChannelReport {
        volume_flow_m3_s: 42.0,
        entrained_flow_m3_s: 43.0,
        pressure_dissipated_power_w: 44.0,
        shear_dissipated_power_w: 45.0,
    };
    for rate in [[0.0, 0.0, 1000.0], [0.0, f64::INFINITY, 0.0]] {
        let mut force = [42.0; 3];
        let mut pressure = [42.0];
        let mut channels = [sentinel; 2];
        assert!(
            moving
                .evaluate_with_channels_into(
                    &[0.0; 3],
                    &rate,
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
    let collapsed = ShearedFilm::new(film(&[0.0, H]), vec![walls(1.0); 2]).unwrap();
    let mut force = [42.0; 3];
    let mut pressure = [42.0];
    let mut channels = [sentinel; 2];
    assert!(
        collapsed
            .evaluate_with_channels_into(
                &[0.0; 3],
                &[0.0, 0.0, 0.002],
                &mut force,
                &mut pressure,
                &mut channels
            )
            .is_err()
    );
    assert_eq!(force, [42.0; 3]);
    assert_eq!(pressure, [42.0]);
    assert_eq!(channels, [sentinel; 2]);
    assert!(
        collapsed
            .tangent_into(&[0.0; 3], &[0.0; 3], &[0.0; 3], &[0.0; 3], &mut force)
            .is_err()
    );
    assert!(
        moving
            .tangent_into(&[0.0; 3], &[0.0; 3], &[f64::NAN; 3], &[0.0; 3], &mut force)
            .is_err()
    );
    assert_eq!(force, [42.0; 3]);
}

#[test]
fn common_translation_does_not_erase_a_sub_ulp_sliding_coordinate_or_its_tangent() {
    // G3: generalized coordinates retain a relative velocity even when adding
    // it to either physical wall's common velocity rounds away that increment.
    // Physical drag depends on the relative row, not subtraction of two rounded
    // endpoint velocities. The two open channel halves keep pressure exactly zero.
    let relative = 1e-17;
    assert_eq!(1.0 + relative, 1.0);
    let original_rows: Vec<_> = [1.0, -1.0]
        .into_iter()
        .map(|sign| WallVelocityPorts {
            lower: vec![0.0, sign, 0.0],
            upper: vec![0.0, sign, sign],
        })
        .collect();
    for stratified in [false, true] {
        let moving = if stratified {
            ShearedFilm::from_stratified(
                StratifiedFilm::new(&film(&[H, H]), vec![layers(), layers()]).unwrap(),
                original_rows.clone(),
            )
            .unwrap()
        } else {
            ShearedFilm::new(film(&[H, H]), original_rows.clone()).unwrap()
        };
        assert_eq!(moving.wall_ports(), original_rows);
        let shear_viscosity = if stratified {
            1.0 / (0.5 / 0.01 + 0.5 / 0.04)
        } else {
            MU
        };
        let drag = W * L * shear_viscosity / H;
        let expected_force = 2.0 * drag * relative;
        let expected_power = expected_force * relative;
        let mut reference = None;
        for translation in [0.0, 1.0, -1.0] {
            let velocity = [0.0, translation, relative];
            let mut force = [0.0; 3];
            let mut pressure = [0.0];
            let report = moving
                .evaluate_into(&[0.0; 3], &velocity, &mut force, &mut pressure)
                .unwrap();
            assert_eq!(pressure, [0.0]);
            assert_eq!(force[0], 0.0);
            assert_eq!(force[1], 0.0);
            near(force[2], expected_force, expected_force * 1e-14);
            near(
                report.dissipated_power_w,
                expected_power,
                expected_power * 1e-14,
            );
            assert!(report.dissipated_power_w > 0.0);
            if let Some((previous_force, previous_power)) = reference {
                assert_eq!(force, previous_force);
                assert_eq!(report.dissipated_power_w, previous_power);
            }
            reference = Some((force, report.dissipated_power_w));
            let mut tangent = [0.0; 3];
            moving
                .tangent_into(
                    &[0.0; 3],
                    &velocity,
                    &[0.0; 3],
                    &[0.0, 1.0, relative],
                    &mut tangent,
                )
                .unwrap();
            assert_eq!(tangent[0], 0.0);
            assert_eq!(tangent[1], 0.0);
            near(tangent[2], expected_force, expected_force * 1e-14);
        }
    }
}

#[test]
fn combined_layer_work_matches_independent_stress_quadrature() {
    // G2: integrate the actual affine shear stress at two Gauss points in each
    // layer. This oracle does not use normalized work fractions or moments.
    let profile = GapViscosityProfile::new(vec![
        ViscosityLayer {
            upper_fraction: 0.17,
            dynamic_viscosity_pa_s: 0.009,
        },
        ViscosityLayer {
            upper_fraction: 0.63,
            dynamic_viscosity_pa_s: 0.08,
        },
        ViscosityLayer {
            upper_fraction: 1.0,
            dynamic_viscosity_pa_s: 0.003,
        },
    ])
    .unwrap();
    let moving = ShearedFilm::from_stratified(
        StratifiedFilm::new(&film(&[H]), vec![profile.clone()]).unwrap(),
        vec![walls(1.0)],
    )
    .unwrap();
    assert_eq!(moving.layer_count(), 3);
    assert_eq!(moving.channel_layer_range(0), Some(0..3));
    assert_eq!(moving.channel_layer_range(1), None);
    let velocity = [0.0001, 0.001, 0.004];
    let mut force = [0.0; 3];
    let mut pressure = [0.0];
    let mut channels = [ShearChannelReport::default()];
    let mut actual = [0.0; 3];
    let report = moving
        .evaluate_with_layers_into(
            &[0.0; 3],
            &velocity,
            &mut force,
            &mut pressure,
            &mut channels,
            &mut actual,
        )
        .unwrap();
    let mut fluidity = 0.0;
    let mut first = 0.0;
    let mut lower = 0.0;
    for layer in profile.layers() {
        let upper = layer.upper_fraction;
        fluidity += (upper - lower) / layer.dynamic_viscosity_pa_s;
        first += 0.5 * (upper * upper - lower * lower) / layer.dynamic_viscosity_pa_s;
        lower = upper;
    }
    let neutral = first / fluidity;
    let couette_stress = (velocity[2] - velocity[1]) / (H * fluidity);
    let pressure_stress = pressure[0] * H / L;
    lower = 0.0;
    for (&power, layer) in actual.iter().zip(profile.layers()) {
        let thickness = layer.upper_fraction - lower;
        let midpoint = 0.5 * (lower + layer.upper_fraction);
        let offset = thickness / (2.0 * 3.0_f64.sqrt());
        let stress = |z: f64| couette_stress - pressure_stress * (z - neutral);
        let expected = W * L * H * thickness / (2.0 * layer.dynamic_viscosity_pa_s)
            * (stress(midpoint - offset).powi(2) + stress(midpoint + offset).powi(2));
        assert!(power > 0.0);
        near(power, expected, 2e-12 * report.dissipated_power_w);
        lower = layer.upper_fraction;
    }
    near(
        actual.iter().sum(),
        channels[0].dissipated_power_w(),
        2e-12 * report.dissipated_power_w,
    );
    near(
        dot(&force, &velocity),
        actual.iter().sum(),
        2e-12 * report.dissipated_power_w,
    );
}

#[test]
fn relative_reversal_exposes_layer_cross_terms_that_whole_gap_work_cancels() {
    // G2: equal half layers and tau_C = ±delta_p*h/(4L) give work fractions
    // (13/14, 1/14), exchanged by relative reversal at the SAME pressure drop.
    // Both pressure and Couette whole-gap powers stay unchanged.
    let profile = GapViscosityProfile::new(vec![
        ViscosityLayer {
            upper_fraction: 0.5,
            dynamic_viscosity_pa_s: MU,
        },
        ViscosityLayer {
            upper_fraction: 1.0,
            dynamic_viscosity_pa_s: MU,
        },
    ])
    .unwrap();
    let moving = ShearedFilm::from_stratified(
        StratifiedFilm::new(&film(&[H]), vec![profile]).unwrap(),
        vec![walls(1.0)],
    )
    .unwrap();
    let prescribed_pressure = 1000.0;
    let relative = prescribed_pressure * H * H / (4.0 * MU * L);
    let conductance = W * H.powi(3) / (12.0 * MU * L);
    let mut result = [[0.0; 2]; 2];
    let mut work = [ShearChannelReport::default(); 2];
    for (case, sign) in [1.0, -1.0].into_iter().enumerate() {
        let upper = sign * relative;
        let swept = (conductance * prescribed_pressure + 0.5 * W * H * upper) / AREA;
        let mut force = [0.0; 3];
        let mut pressure = [0.0];
        let mut channels = [ShearChannelReport::default()];
        moving
            .evaluate_with_layers_into(
                &[0.0; 3],
                &[swept, 0.0, upper],
                &mut force,
                &mut pressure,
                &mut channels,
                &mut result[case],
            )
            .unwrap();
        work[case] = channels[0];
        near(pressure[0], prescribed_pressure, 1e-10);
        let total = channels[0].dissipated_power_w();
        let expected = if sign > 0.0 {
            [13.0 / 14.0, 1.0 / 14.0]
        } else {
            [1.0 / 14.0, 13.0 / 14.0]
        };
        for (actual, fraction) in result[case].iter().zip(expected) {
            near(*actual, total * fraction, total * 2e-12);
        }
        assert!((result[case][0] - 0.5 * total).abs() > 0.4 * total);
    }
    near(result[0][0], result[1][1], 1e-15);
    near(result[0][1], result[1][0], 1e-15);
    near(
        work[0].pressure_dissipated_power_w,
        work[1].pressure_dissipated_power_w,
        1e-15,
    );
    near(
        work[0].shear_dissipated_power_w,
        work[1].shear_dissipated_power_w,
        1e-15,
    );
}

#[test]
fn layer_reporting_recovers_native_uniform_pressure_only_and_closed_limits() {
    // G3: implicit uniform channels are exactly the native channel account.
    let moving = ShearedFilm::new(film(&[H, 0.0]), vec![walls(1.0), walls(0.0)]).unwrap();
    assert_eq!(moving.layer_count(), 2);
    assert_eq!(moving.channel_layer_range(0), Some(0..1));
    assert_eq!(moving.channel_layer_range(1), Some(1..2));
    let velocity = [0.0001, 0.001, 0.004];
    let mut force = [0.0; 3];
    let mut pressure = [0.0];
    let mut channels = [ShearChannelReport::default(); 2];
    let mut layer_power = [0.0; 2];
    let report = moving
        .evaluate_with_layers_into(
            &[0.0; 3],
            &velocity,
            &mut force,
            &mut pressure,
            &mut channels,
            &mut layer_power,
        )
        .unwrap();
    for (layer, channel) in layer_power.iter().zip(channels) {
        assert_eq!(*layer, channel.dissipated_power_w());
    }
    assert_eq!(layer_power[1], 0.0);
    let expected_force = force;
    let expected_pressure = pressure;
    let expected_channels = channels;
    let legacy = moving
        .evaluate_with_channels_into(
            &[0.0; 3],
            &velocity,
            &mut force,
            &mut pressure,
            &mut channels,
        )
        .unwrap();
    assert_eq!(force, expected_force);
    assert_eq!(pressure, expected_pressure);
    assert_eq!(channels, expected_channels);
    assert_eq!(report.minimum_gap_m, legacy.minimum_gap_m);
    assert_eq!(report.peak_pressure_pa, legacy.peak_pressure_pa);
    assert_eq!(report.dissipated_power_w, legacy.dissipated_power_w);

    let profile = layers();
    let fixed = ShearedFilm::from_stratified(
        StratifiedFilm::new(&film(&[H, 0.0]), vec![profile.clone(); 2]).unwrap(),
        vec![walls(0.0); 2],
    )
    .unwrap();
    assert_eq!(fixed.layer_count(), 4);
    assert_eq!(fixed.channel_layer_range(0), Some(0..2));
    assert_eq!(fixed.channel_layer_range(1), Some(2..4));
    let mut layer_power = [0.0; 4];
    fixed
        .evaluate_with_layers_into(
            &[0.0; 3],
            &[0.001, 0.0, 0.0],
            &mut force,
            &mut pressure,
            &mut channels,
            &mut layer_power,
        )
        .unwrap();
    for (channel, report) in channels.iter().enumerate() {
        for (local, fraction) in profile.layer_pressure_work_fractions().iter().enumerate() {
            assert_eq!(
                layer_power[2 * channel + local],
                report.pressure_dissipated_power_w * fraction
            );
        }
    }
    assert_eq!(&layer_power[2..], &[0.0, 0.0]);
}

#[test]
fn layer_partition_and_velocity_rescaling_preserve_physical_work() {
    // G3: subdividing a physical layer changes only report grouping.
    let coarse = layers();
    let fine = GapViscosityProfile::new(vec![
        ViscosityLayer {
            upper_fraction: 0.2,
            dynamic_viscosity_pa_s: 0.01,
        },
        ViscosityLayer {
            upper_fraction: 0.5,
            dynamic_viscosity_pa_s: 0.01,
        },
        ViscosityLayer {
            upper_fraction: 0.8,
            dynamic_viscosity_pa_s: 0.04,
        },
        ViscosityLayer {
            upper_fraction: 1.0,
            dynamic_viscosity_pa_s: 0.04,
        },
    ])
    .unwrap();
    let build = |profile| {
        ShearedFilm::from_stratified(
            StratifiedFilm::new(&film(&[H]), vec![profile]).unwrap(),
            vec![walls(1.0)],
        )
        .unwrap()
    };
    let coarse = build(coarse);
    let fine = build(fine);
    let base_velocity = [0.0001, 0.001, 0.004];
    let mut reference = [0.0; 2];
    for scale in [1.0, -1.0, 1e-100] {
        let velocity = base_velocity.map(|v| v * scale);
        let mut force = [0.0; 3];
        let mut pressure = [0.0];
        let mut channels = [ShearChannelReport::default()];
        let mut grouped = [0.0; 2];
        let report = coarse
            .evaluate_with_layers_into(
                &[0.0; 3],
                &velocity,
                &mut force,
                &mut pressure,
                &mut channels,
                &mut grouped,
            )
            .unwrap();
        let mut split = [0.0; 4];
        fine.evaluate_with_layers_into(
            &[0.0; 3],
            &velocity,
            &mut force,
            &mut pressure,
            &mut channels,
            &mut split,
        )
        .unwrap();
        near(
            split[0] + split[1],
            grouped[0],
            report.dissipated_power_w * 3e-12,
        );
        near(
            split[2] + split[3],
            grouped[1],
            report.dissipated_power_w * 3e-12,
        );
        if scale == 1.0 {
            reference = grouped;
        }
        for (actual, base) in grouped.into_iter().zip(reference) {
            assert!(actual > 0.0);
            near(actual / (base * scale * scale), 1.0, 3e-12);
        }
    }
}

#[test]
fn layer_reporting_refusals_leave_all_outputs_unchanged() {
    // G4: a positive channel power may have an unrepresentably small layer
    // contribution. Refuse that report after the flow solve without publishing
    // the otherwise valid force, pressure or channel result.
    let profile = GapViscosityProfile::new(vec![
        ViscosityLayer {
            upper_fraction: 0.5,
            dynamic_viscosity_pa_s: 0.01,
        },
        ViscosityLayer {
            upper_fraction: 1.0,
            dynamic_viscosity_pa_s: 1e200,
        },
    ])
    .unwrap();
    let moving = ShearedFilm::from_stratified(
        StratifiedFilm::new(&film(&[H]), vec![profile]).unwrap(),
        vec![walls(0.0)],
    )
    .unwrap();
    let sentinel = ShearChannelReport {
        volume_flow_m3_s: 42.0,
        entrained_flow_m3_s: 43.0,
        pressure_dissipated_power_w: 44.0,
        shear_dissipated_power_w: 45.0,
    };
    let mut force = [42.0; 3];
    let mut pressure = [42.0];
    let mut channels = [sentinel];
    let mut wrong_shape = [42.0; 3];
    assert!(
        moving
            .evaluate_with_layers_into(
                &[0.0; 3],
                &[0.001, 0.0, 0.0],
                &mut force,
                &mut pressure,
                &mut channels,
                &mut wrong_shape,
            )
            .is_err()
    );
    assert_eq!(force, [42.0; 3]);
    assert_eq!(pressure, [42.0]);
    assert_eq!(channels, [sentinel]);
    assert_eq!(wrong_shape, [42.0; 3]);
    let mut layer_power = [42.0; 2];
    assert!(
        moving
            .evaluate_with_layers_into(
                &[0.0; 3],
                &[1e-100, 0.0, 0.0],
                &mut force,
                &mut pressure,
                &mut channels,
                &mut layer_power,
            )
            .is_err()
    );
    assert_eq!(force, [42.0; 3]);
    assert_eq!(pressure, [42.0]);
    assert_eq!(channels, [sentinel]);
    assert_eq!(layer_power, [42.0; 2]);
    let allowed = moving
        .evaluate_with_channels_into(
            &[0.0; 3],
            &[1e-100, 0.0, 0.0],
            &mut force,
            &mut pressure,
            &mut channels,
        )
        .unwrap();
    assert!(allowed.dissipated_power_w > 0.0);
}
