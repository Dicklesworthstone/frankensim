use super::*;

#[test]
fn pressure_dependent_wall_force_and_layer_heat_match_local_stokes_work() {
    // G1: independently integrate mu(p,z)*(du/dy)^2 along and across a
    // passage. Pressure/Couette cross terms differ between retained layers.
    let bounds: [(f64, f64, f64); 2] = [(0.0, 0.5, 0.01), (0.5, 1.0, 0.04)];
    let fluidity = 0.5 / 0.01 + 0.5 / 0.04;
    let z0 = (0.125 / 0.01 + 0.375 / 0.04) / fluidity;
    let mut mobility = 0.0;
    for (lo, hi, mu) in bounds {
        mobility += ((hi - z0).powi(3) - (lo - z0).powi(3)) / (3.0 * mu);
    }
    let conductance = W * H.powi(3) * mobility / L;
    let lower = 0.001;
    let upper = 0.006;
    let relative = upper - lower;
    let entrainment = W * H * (z0 * lower + (1.0 - z0) * upper);
    for alpha in [0.0, 1e-12, 3e-4] {
        let native = film(&[H]).with_barus_viscosity(alpha).unwrap();
        let model = ShearedFilm::from_stratified(
            StratifiedFilm::new(&native, vec![layers()]).unwrap(),
            vec![walls(1.0)],
        )
        .unwrap();
        for normal in [-0.0002, 0.0002] {
            let velocity = [normal, lower, upper];
            let psi = (AREA * normal - entrainment) / conductance;
            let expected_pressure = if alpha == 0.0 {
                psi
            } else {
                -(-alpha * psi).ln_1p() / alpha
            };
            let mut force = [0.0; 3];
            let mut pressure = [0.0];
            let mut channels = [ShearChannelReport::default()];
            let mut heat = [0.0; 2];
            let report = model
                .evaluate_with_layers_into(
                    &[0.0; 3],
                    &velocity,
                    &mut force,
                    &mut pressure,
                    &mut channels,
                    &mut heat,
                )
                .unwrap();
            near(
                pressure[0],
                expected_pressure,
                expected_pressure.abs() * 2e-14,
            );
            near(channels[0].volume_flow_m3_s, AREA * normal, 1e-22);
            near(channels[0].entrained_flow_m3_s, entrainment, 1e-23);
            near(dot(&force, &velocity), report.dissipated_power_w, 1e-17);
            for (layer, (lo, hi, mu)) in bounds.into_iter().enumerate() {
                // Composite Simpson rule is exact in z (quadratic), and
                // independently resolves the rational streamwise viscosity.
                let mut integral = 0.0;
                let nx = 256;
                let nz = 16;
                for ix in 0..=nx {
                    let x = ix as f64 / nx as f64;
                    let streamwise = 1.0 / (1.0 - alpha * psi * (1.0 - x));
                    let wx = if ix == 0 || ix == nx {
                        1.0
                    } else if ix % 2 == 0 {
                        2.0
                    } else {
                        4.0
                    };
                    for iz in 0..=nz {
                        let z = lo + (hi - lo) * iz as f64 / nz as f64;
                        let wz = if iz == 0 || iz == nz {
                            1.0
                        } else if iz % 2 == 0 {
                            2.0
                        } else {
                            4.0
                        };
                        let gradient =
                            relative / (H * fluidity * mu) - psi * H * (z - z0) / (L * mu);
                        integral += wx * wz * mu * streamwise * gradient * gradient;
                    }
                }
                integral *= W * L * H * (hi - lo) / (9.0 * nx as f64 * nz as f64);
                near(heat[layer], integral, integral * 2e-10);
            }
            near(heat.iter().sum(), report.dissipated_power_w, 1e-17);
            // Both signs of pressure are covered with the same unequal wall speeds.
            assert_eq!(pressure[0].signum(), normal.signum());
        }
    }
}

#[test]
fn barus_network_conserves_flow_work_and_differentiates_full_wall_response() {
    // G0/G3: all existing ports participate, including gap and wall directions.
    for alpha in [1e-12, 3e-4] {
        let original = network();
        let model = ShearedFilm::compile(
            original.film.with_barus_viscosity(alpha).unwrap(),
            original.profiles,
            original.walls,
        )
        .unwrap();
        let q = [0.00001, -0.00002, 0.0, 0.0];
        let dq = [0.00002, -0.00001, 0.0, 0.0];
        let dv = [0.0001, -0.0002, 0.003, -0.002];
        for sign in [-1.0, 1.0] {
            let v = [0.0003 * sign, -0.0002 * sign, 0.001, 0.006];
            let mut f = [0.0; 4];
            let mut p = [0.0; 2];
            let mut channels = [ShearChannelReport::default(); 3];
            let report = model
                .evaluate_with_channels_into(&q, &v, &mut f, &mut p, &mut channels)
                .unwrap();
            near(
                channels[0].volume_flow_m3_s + channels[1].volume_flow_m3_s,
                AREA * v[0],
                1e-21,
            );
            near(
                channels[2].volume_flow_m3_s - channels[1].volume_flow_m3_s,
                0.75 * AREA * v[1],
                1e-21,
            );
            near(dot(&f, &v), report.dissipated_power_w, 1e-16);
            let mut tangent = [0.0; 4];
            model.tangent_into(&q, &v, &dq, &dv, &mut tangent).unwrap();
            let epsilon = 1e-5;
            let mut force = [[0.0; 4]; 2];
            for (i, direction) in [-1.0, 1.0].into_iter().enumerate() {
                let position =
                    std::array::from_fn::<_, 4, _>(|j| q[j] + direction * epsilon * dq[j]);
                let rate = std::array::from_fn::<_, 4, _>(|j| v[j] + direction * epsilon * dv[j]);
                model
                    .evaluate_into(&position, &rate, &mut force[i], &mut p)
                    .unwrap();
            }
            for i in 0..4 {
                let fd = (force[1][i] - force[0][i]) / (2.0 * epsilon);
                near(tangent[i], fd, 1e-8 * tangent[i].abs() + 1e-10);
            }
        }
    }
}

#[test]
fn equal_pressure_sliding_retains_ambient_pressure_viscosity_and_its_tangent() {
    // G1: the middle passage has zero pressure drop but nonzero absolute
    // gauge pressure. Its Couette drag must use mu(p), not reservoir mu_0.
    let gap = GapPort {
        reference_m: H,
        closure: vec![1.0, 0.0, 0.0],
    };
    let alpha = 3e-4;
    let native = ResistiveFilm::new(
        vec![
            FilmCell {
                area_m2: AREA,
                gap: gap.clone()
            };
            2
        ],
        [(0, None), (0, Some(1)), (1, None)]
            .into_iter()
            .map(|(from, to)| FilmChannel {
                from,
                to,
                width_m: W,
                length_m: L,
                gap: gap.clone(),
            })
            .collect(),
        3,
        MU,
        film(&[H]).limits,
    )
    .unwrap()
    .with_barus_viscosity(alpha)
    .unwrap();
    let model = ShearedFilm::new(native, vec![walls(0.0), walls(1.0), walls(0.0)]).unwrap();
    let v = [0.0002, -0.003, 0.003];
    let mut f = [0.0; 3];
    let mut p = [0.0; 2];
    let mut channels = [ShearChannelReport::default(); 3];
    model
        .evaluate_with_channels_into(&[0.0; 3], &v, &mut f, &mut p, &mut channels)
        .unwrap();
    near(p[0], p[1], 1e-11);
    let multiplier = (alpha * p[0]).exp();
    let expected = MU * W * L / H * multiplier * (v[2] - v[1]);
    near(f[2], expected, expected.abs() * 2e-14);
    near(f[1], -expected, expected.abs() * 2e-14);
    near(channels[1].pressure_dissipated_power_w, 0.0, 1e-25);
    let mut derivative = [0.0; 3];
    model
        .tangent_into(&[0.0; 3], &v, &[0.0; 3], &[1.0, 0.0, 0.0], &mut derivative)
        .unwrap();
    let conductance = W * H.powi(3) / (12.0 * MU * L);
    near(
        derivative[2],
        expected * alpha * multiplier * AREA / conductance,
        derivative[2].abs() * 2e-14,
    );
    // G3: the exact coincident limit also admits independently changing ends.
    let psi = 800.0;
    let (mean, tangent) = barus_mean(alpha, psi, psi, 3.0, -1.0).unwrap();
    let r = 1.0 - alpha * psi;
    near(mean, 1.0 / r, 1e-15);
    near(tangent, alpha / (r * r), 1e-18);
}

#[test]
fn fixed_walls_match_native_barus_and_pressure_barriers_remain_atomic() {
    // G3/G4: pressure-only compatibility and original finite-pressure refusal.
    let native = film(&[H, H]).with_barus_viscosity(3e-4).unwrap();
    let fixed = ShearedFilm::new(native.clone(), vec![walls(0.0); 2]).unwrap();
    let v = [0.0002, 0.0, 0.0];
    let mut expected = [0.0; 3];
    let mut actual = [0.0; 3];
    let mut expected_pressure = [0.0];
    let mut actual_pressure = [0.0];
    let a = native
        .evaluate_into(&[0.0; 3], &v, &mut expected, &mut expected_pressure)
        .unwrap();
    let b = fixed
        .evaluate_into(&[0.0; 3], &v, &mut actual, &mut actual_pressure)
        .unwrap();
    assert_eq!(expected, actual);
    assert_eq!(expected_pressure, actual_pressure);
    assert_eq!(a.dissipated_power_w, b.dissipated_power_w);
    let dq = [H * 0.05, 0.0, 0.0];
    native
        .tangent_into(&[0.0; 3], &v, &dq, &v, &mut expected)
        .unwrap();
    fixed
        .tangent_into(&[0.0; 3], &v, &dq, &v, &mut actual)
        .unwrap();
    assert_eq!(expected, actual);
    let moving = ShearedFilm::new(native, vec![walls(1.0); 2]).unwrap();
    actual = [42.0; 3];
    actual_pressure = [43.0];
    assert!(
        moving
            .evaluate_into(
                &[0.0; 3],
                &[1.0, 0.002, -0.001],
                &mut actual,
                &mut actual_pressure
            )
            .is_err()
    );
    assert_eq!(actual, [42.0; 3]);
    assert_eq!(actual_pressure, [43.0]);
    assert!(
        moving
            .tangent_into(&[0.0; 3], &[1.0, 0.002, -0.001], &[0.0; 3], &v, &mut actual)
            .is_err()
    );
    assert_eq!(actual, [42.0; 3]);
}
