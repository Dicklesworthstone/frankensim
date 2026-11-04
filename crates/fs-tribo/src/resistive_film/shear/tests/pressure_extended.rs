//! Additional multicell, layer-integral and high-contrast physical oracles.
//! The pressure-fed owner and its published API remain the implementation.

use super::*;
use crate::resistive_film::pressure::PressureDrivenFilm;

fn close(a: f64, b: f64) {
    near(a, b, 3e-12 * a.abs().max(b.abs()).max(1e-14));
}

fn evaluate(
    model: &PressureDrivenShearedFilm,
    q: &[f64],
    v: &[f64],
) -> (
    PressureDrivenReport,
    Vec<f64>,
    Vec<f64>,
    Vec<ShearChannelReport>,
    Vec<f64>,
) {
    let mut force = vec![0.0; model.port_count()];
    let mut pressure = vec![0.0; model.cell_count()];
    let mut channels = vec![ShearChannelReport::default(); model.channel_count()];
    let mut layers = vec![0.0; model.layer_count()];
    let report = model
        .evaluate_with_layers_into(q, v, &mut force, &mut pressure, &mut channels, &mut layers)
        .unwrap();
    (report, force, pressure, channels, layers)
}

#[test]
fn supplied_layer_heat_matches_independent_local_stokes_stress_quadrature() {
    // G1: independently integrate tau^2/mu along and across each layer. This
    // tests the cross term, unequal wall speeds, nonzero supplies and Barus.
    let profile = layers();
    let a = 0.5 / 0.01 + 0.5 / 0.04;
    let z0 = (0.125 / 0.01 + 0.375 / 0.04) / a;
    let v = [0.0, -0.003, 0.009];
    for alpha in [0.0, 3e-4] {
        let base = film(&[H, H]).with_barus_viscosity(alpha).unwrap();
        let stratified = StratifiedFilm::new(&base, vec![profile.clone(); 2]).unwrap();
        let sliding =
            ShearedFilm::from_stratified(stratified, vec![walls(-1.0), walls(1.0)]).unwrap();
        let supplies = [1800.0, -300.0];
        let model = PressureDrivenShearedFilm::new(sliding, &supplies.map(Some)).unwrap();
        let (report, f, p, channels, powers) = evaluate(&model, &[0.0; 3], &v);
        let fluidity = |p: f64| (-alpha * p).exp();
        let center = if alpha == 0.0 {
            (supplies[0] + supplies[1]) / 2.0
        } else {
            -((fluidity(supplies[0]) + fluidity(supplies[1])) / 2.0).ln() / alpha
        };
        close(p[0], center);
        let mut expected_force = [AREA * center, 0.0, 0.0];
        for channel in 0..2 {
            let sign = if channel == 0 { -1.0 } else { 1.0 };
            let relative = sign * (v[2] - v[1]);
            let r0 = fluidity(center);
            let r1 = fluidity(supplies[channel]);
            let delta_psi = if alpha == 0.0 {
                center - supplies[channel]
            } else {
                (r1 - r0) / alpha
            };
            let base_stress = relative / (H * a);
            let gradient_stress = delta_psi * H / L;
            let mut expected = [0.0; 2];
            let mut wall_forces = [0.0; 2];
            // Composite Simpson in x; Simpson is exact for the quadratic
            // stress within each constant-mu layer in normalized gap z.
            let n = 256;
            for i in 0..=n {
                let x = i as f64 / n as f64;
                let r = r0 * (1.0 - x) + r1 * x;
                let weight = if i == 0 || i == n {
                    1.0
                } else if i % 2 == 0 {
                    2.0
                } else {
                    4.0
                };
                let dx_weight = weight * L / (3.0 * n as f64);
                let stress = |z: f64| (base_stress - gradient_stress * (z - z0)) / r;
                wall_forces[0] -= W * stress(0.0) * dx_weight;
                wall_forces[1] += W * stress(1.0) * dx_weight;
                for (layer, (low, high, mu)) in
                    [(0.0, 0.5, 0.01), (0.5, 1.0, 0.04)].into_iter().enumerate()
                {
                    let work = |z: f64| W * H * stress(z).powi(2) * r / mu;
                    expected[layer] += dx_weight * (high - low) / 6.0
                        * (work(low) + 4.0 * work((low + high) / 2.0) + work(high));
                }
            }
            for layer in 0..2 {
                near(
                    powers[2 * channel + layer],
                    expected[layer],
                    expected[layer] * 2e-11,
                );
            }
            near(
                channels[channel].dissipated_power_w(),
                expected.iter().sum(),
                expected.iter().sum::<f64>() * 2e-11,
            );
            expected_force[1] += sign * wall_forces[0];
            expected_force[2] += sign * wall_forces[1];
        }
        for j in 0..3 {
            near(f[j], expected_force[j], expected_force[j].abs() * 2e-11);
        }
        close(powers.iter().sum(), report.film.dissipated_power_w);
        close(
            report.mechanical_power_w + report.reservoir_power_w,
            report.film.dissipated_power_w,
        );
    }
}

fn network(alpha: f64) -> ShearedFilm {
    let gap = |h| GapPort {
        reference_m: h,
        closure: vec![1.0, 0.0, 0.0],
    };
    let base = ResistiveFilm::new(
        vec![
            FilmCell {
                area_m2: AREA,
                gap: gap(H)
            };
            2
        ],
        vec![
            FilmChannel {
                from: 0,
                to: None,
                width_m: W,
                length_m: L,
                gap: gap(H),
            },
            FilmChannel {
                from: 0,
                to: Some(1),
                width_m: W * 1.2,
                length_m: L * 0.8,
                gap: gap(H * 1.1),
            },
            FilmChannel {
                from: 1,
                to: None,
                width_m: W * 0.9,
                length_m: L * 1.1,
                gap: gap(H * 0.9),
            },
        ],
        3,
        MU,
        FilmLimits {
            minimum_cell_gap_m: H * 0.1,
            maximum_gap_m: H * 2.0,
            maximum_pressure_pa: 10000.0,
        },
    )
    .unwrap()
    .with_barus_viscosity(alpha)
    .unwrap();
    ShearedFilm::from_stratified(
        StratifiedFilm::new(&base, vec![layers(); 3]).unwrap(),
        vec![walls(-1.0), walls(1.0), walls(1.0)],
    )
    .unwrap()
}

#[test]
fn complete_multicell_tangent_includes_geometry_squeeze_walls_and_both_supplies() {
    // G1: differentiate the full driven force through independent central
    // differences, including the reservoir used as the numerical datum.
    for alpha in [0.0, 2e-4] {
        let supplies = [Some(1500.0), None, Some(-300.0)];
        let model = PressureDrivenShearedFilm::new(network(alpha), &supplies).unwrap();
        let q = [H * 0.1, 0.0, 0.0];
        let v = [0.00003, 0.013, -0.009];
        let dq = [H * 0.03, 0.0, 0.0];
        let dv = [-0.00001, 0.005, -0.002];
        let dp = [200.0, 0.0, -100.0];
        let mut tangent = [0.0; 3];
        model
            .tangent_into(&q, &v, &dq, &dv, &dp, &mut tangent)
            .unwrap();
        let eps = 1e-5;
        let shifted_force = |sign: f64| {
            let pressure =
                std::array::from_fn::<_, 3, _>(|k| supplies[k].map(|p| p + sign * eps * dp[k]));
            let shifted = model.with_reservoir_pressures(&pressure).unwrap();
            let qt = std::array::from_fn::<_, 3, _>(|j| q[j] + sign * eps * dq[j]);
            let vt = std::array::from_fn::<_, 3, _>(|j| v[j] + sign * eps * dv[j]);
            evaluate(&shifted, &qt, &vt).1
        };
        let plus = shifted_force(1.0);
        let minus = shifted_force(-1.0);
        for j in 0..3 {
            near(
                tangent[j],
                (plus[j] - minus[j]) / (2.0 * eps),
                tangent[j].abs() * 2e-7 + 1e-12,
            );
        }
        let (_, _, _, c, _) = evaluate(&model, &q, &v);
        close(c[0].volume_flow_m3_s + c[1].volume_flow_m3_s, AREA * v[0]);
        close(c[2].volume_flow_m3_s - c[1].volume_flow_m3_s, AREA * v[0]);
    }
}

#[test]
fn zero_supplies_preserve_native_bits_and_fixed_walls_recover_pressure_only() {
    // G3: neither the existing passive path nor stationary-wall mechanics is
    // weakened by composition. Layers and channel ordering remain original.
    for alpha in [0.0, 2e-4] {
        let native = network(alpha);
        let model =
            PressureDrivenShearedFilm::new(native.clone(), &[Some(0.0), None, Some(0.0)]).unwrap();
        let q = [0.0; 3];
        let v = [0.00003, 0.004, -0.001];
        let (r, f, p, c, l) = evaluate(&model, &q, &v);
        let mut nf = [0.0; 3];
        let mut np = [0.0; 2];
        let mut nc = [ShearChannelReport::default(); 3];
        let mut nl = [0.0; 6];
        let nr = native
            .evaluate_with_layers_into(&q, &v, &mut nf, &mut np, &mut nc, &mut nl)
            .unwrap();
        assert_eq!(f, nf);
        assert_eq!(p, np);
        assert_eq!(c, nc);
        assert_eq!(l, nl);
        assert_eq!(r.film.dissipated_power_w, nr.dissipated_power_w);
        assert_eq!(r.reservoir_power_w, 0.0);
        native
            .tangent_into(&q, &v, &[H * 0.05, 0.0, 0.0], &v, &mut nf)
            .unwrap();
        let mut df = [0.0; 3];
        model
            .tangent_into(&q, &v, &[H * 0.05, 0.0, 0.0], &v, &[0.0; 3], &mut df)
            .unwrap();
        assert_eq!(df, nf);

        let base = film(&[H, 0.0]).with_barus_viscosity(alpha).unwrap();
        let supplies = [Some(1600.0), Some(-200.0)];
        let pressure_only = PressureDrivenFilm::new(base.clone(), &supplies).unwrap();
        let fixed = ShearedFilm::new(base, vec![walls(0.0); 2]).unwrap();
        let driven = PressureDrivenShearedFilm::new(fixed, &supplies).unwrap();
        let (report, force, pressure, channels, _) = evaluate(&driven, &q, &v);
        let mut expected_f = [0.0; 3];
        let mut expected_p = [0.0];
        let mut expected_c = [FilmChannelReport::default(); 2];
        let expected = pressure_only
            .evaluate_into(&q, &v, &mut expected_f, &mut expected_p, &mut expected_c)
            .unwrap();
        assert_eq!(force, expected_f);
        assert_eq!(pressure, expected_p);
        close(
            report.film.dissipated_power_w,
            expected.film.dissipated_power_w,
        );
        close(report.reservoir_power_w, expected.reservoir_power_w);
        for k in 0..2 {
            assert_eq!(channels[k].volume_flow_m3_s, expected_c[k].volume_flow_m3_s);
            close(
                channels[k].dissipated_power_w(),
                expected_c[k].dissipated_power_w,
            );
        }
    }
}

#[test]
fn high_fluidity_contrast_retains_finite_pressure_drop_and_fast_force_outputs() {
    // A narrow outlet and strongly pressure-dependent viscosity leave a
    // perfectly finite high-pressure cell. In the former one-sided recovery,
    // alpha*delta_psi/r_out rounds to one and log1p falsely becomes singular.
    let base = film(&[H, H * 1e-6]).with_barus_viscosity(1.0).unwrap();
    let supplies = [Some(40.0), Some(0.0)];
    let pressure_only = PressureDrivenFilm::new(base.clone(), &supplies).unwrap();
    let sliding = ShearedFilm::new(base, vec![walls(-1.0), walls(1.0)]).unwrap();
    let model = PressureDrivenShearedFilm::new(sliding, &supplies).unwrap();
    let (report, force, pressure, channels, layers) = evaluate(&model, &[0.0; 3], &[0.0; 3]);
    let high_g = W * H.powi(3) / (12.0 * MU * L);
    let low_g = W * (H * 1e-6).powi(3) / (12.0 * MU * L);
    let expected_pressure = -((high_g * (-40.0_f64).exp() + low_g) / (high_g + low_g)).ln();
    near(pressure[0], expected_pressure, 2e-12);
    close(force[0], AREA * expected_pressure);
    assert!(report.film.dissipated_power_w > 0.0);
    close(report.reservoir_power_w, report.film.dissipated_power_w);
    close(layers.iter().sum(), report.film.dissipated_power_w);
    assert!(channels.iter().all(|c| c.dissipated_power_w().is_finite()));
    let mut fast_force = [0.0; 3];
    let mut fast_pressure = [0.0];
    let fast = model.evaluate_force_into(
        &[0.0; 3], &[0.0; 3], &mut fast_force, &mut fast_pressure,
    ).unwrap();
    assert_eq!(fast_force.as_slice(), force.as_slice());
    assert_eq!(fast_pressure.as_slice(), pressure.as_slice());
    assert_eq!(fast.film.dissipated_power_w, report.film.dissipated_power_w);
    assert_eq!(fast.reservoir_power_w, report.reservoir_power_w);
    let mut pressure_force = [0.0; 3];
    let mut pressure_cells = [0.0];
    let mut pressure_channels = [FilmChannelReport::default(); 2];
    let fixed = pressure_only.evaluate_into(
        &[0.0; 3], &[0.0; 3], &mut pressure_force, &mut pressure_cells, &mut pressure_channels,
    ).unwrap();
    near(pressure_cells[0], expected_pressure, 2e-12);
    close(fixed.film.dissipated_power_w, report.film.dissipated_power_w);
    for k in 0..2 {
        close(pressure_channels[k].dissipated_power_w, channels[k].dissipated_power_w());
    }

    // The force-only path retains the same atomic refusal as full reports.
    let mut bad_force = [17.0; 3];
    let mut bad_pressure = [19.0];
    assert!(model.evaluate_force_into(
        &[H; 3], &[0.0; 3], &mut bad_force, &mut bad_pressure,
    ).is_err());
    assert_eq!(bad_force, [17.0; 3]);
    assert_eq!(bad_pressure, [19.0]);
}

#[test]
fn finite_viscosity_at_subnormal_fluidity_preserves_pressure_only_and_fixed_walls() {
    // The reference viscosity is deliberately tiny: actual viscosity is finite
    // even though the standalone exp(alpha*p) multiplier would overflow.
    // This probes the previously admitted pressure-only arithmetic domain.
    let mu = 1e-300;
    let base = film(&[H, H]).with_viscosity(mu).unwrap()
        .with_barus_viscosity(1.0).unwrap();
    for low in [715.0, 714.0] {
        let supplies = [Some(715.0), Some(low)];
        let pressure_only = PressureDrivenFilm::new(base.clone(), &supplies).unwrap();
        let sliding = ShearedFilm::new(base.clone(), vec![walls(0.0); 2]).unwrap();
        let model = PressureDrivenShearedFilm::new(sliding, &supplies).unwrap();
        let (report, force, pressure, _, _) = evaluate(&model, &[0.0; 3], &[0.0; 3]);
        let expected_pressure = 715.0 - ((1.0 + (715.0-low).exp()) / 2.0).ln();
        near(pressure[0], expected_pressure, 1e-10);
        near(force[0], AREA * expected_pressure, 1e-14);
        let g = W * H.powi(3) / (12.0 * mu * L);
        let expected_heat = 0.5 * g * (-715.0_f64).exp()
            * (715.0-low).exp_m1() * (715.0-low);
        near(report.film.dissipated_power_w, expected_heat,
            1e-10 * expected_heat.max(f64::MIN_POSITIVE));
        let mut fixed_force = [0.0; 3];
        let mut fixed_pressure = [0.0];
        let mut channels = [FilmChannelReport::default(); 2];
        let fixed = pressure_only.evaluate_into(
            &[0.0; 3], &[0.0; 3], &mut fixed_force, &mut fixed_pressure, &mut channels,
        ).unwrap();
        near(fixed_pressure[0], expected_pressure, 1e-10);
        near(fixed.film.dissipated_power_w, expected_heat,
            1e-10 * expected_heat.max(f64::MIN_POSITIVE));
        if low == 715.0 {
            assert_eq!(fixed.film.dissipated_power_w, 0.0);
            assert_eq!(report.film.dissipated_power_w, 0.0);
        } else {
            assert!(fixed.film.dissipated_power_w > 0.0);
            assert!(report.film.dissipated_power_w > 0.0);
        }
    }
}
