//! Original regional media, complete rational memory, and independent physics.
use super::*;
use fs_couple::mm_line::continuous::rational::RationalModalError;

fn hot_air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 493.15, 101_325.0).unwrap()
}

fn regions(elements: usize) -> Vec<GasState> {
    assert_eq!(elements % 2, 0);
    (0..elements)
        .map(|element| {
            if element < elements / 2 {
                air()
            } else {
                hot_air()
            }
        })
        .collect()
}

fn runtime_budget() -> RationalModalLineBudget {
    RationalModalLineBudget {
        linear: LinearSplitBudget {
            max_modes: 124,
            max_omega_dt: 1.0,
        },
        exact: DenseExponentialBudget {
            max_dimension: 256,
            max_dense_bytes: 64 << 20,
            max_multiply_terms: 20_000_000_000,
        },
        relative_tolerance: 1e-8,
    }
}

fn runtime_source() -> MimoModel {
    let mut model = source(2, true);
    // Full-rank positive direct impedance admits nonzero original history.
    model.d[0] += 0.25 * reference();
    model.d[3] += 0.25 * reference();
    model
}

fn close_real(a: f64, b: f64, scale: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * scale.max(1e-100),
        "{a:e} != {b:e}; scale={scale:e}"
    );
}

#[test]
fn repeated_uniform_regions_preserve_every_state_history_clock_and_response() {
    let model = runtime_source();
    let gate = CancelGate::new();
    let gases = [air(); 4];
    let boundary = load(&model, 1e-5, true);
    let accepted = boundary.load().continuous_operator().real_coordinates();
    let mut uniform =
        RationalModalLine::new(&profile(), &air(), 4, boundary, runtime_budget(), &gate).unwrap();
    let mut regional = RationalModalLine::new_regional(
        &profile(),
        &gases,
        load(&model, 1e-5, true),
        runtime_budget(),
        &gate,
    )
    .unwrap();
    assert_eq!(regional.source().regional_gases(), Some(gases.as_slice()));
    assert!(uniform.source().regional_gases().is_none());
    assert_eq!(uniform.source().mass(), regional.source().mass());
    assert_eq!(uniform.source().stiffness(), regional.source().stiffness());
    assert_eq!(
        uniform.source().inverse_mass(),
        regional.source().inverse_mass()
    );
    assert_eq!(regional.load_state(), accepted);
    assert_eq!(regional.load().source_model(), &model);
    assert_eq!(regional.duration(), 1e-5);
    assert_eq!(regional.boundary().exterior(), &air());
    assert_eq!(regional.load().limits().band_rad_s, [1.0, 4e6]);
    for step in 0..16 {
        let wave = [20.0 * (0.17 * step as f64).sin(), -8.0];
        assert_eq!(
            uniform.push(&wave, &gate).unwrap(),
            regional.push(&wave, &gate).unwrap()
        );
        assert_eq!(uniform.state(), regional.state());
        assert_eq!(uniform.incoming(), regional.incoming());
        assert_eq!(uniform.terminal_accounts(), regional.terminal_accounts());
        assert_eq!(
            regional.load_state(),
            regional.load().continuous_operator().real_coordinates()
        );
    }
    let before = regional.state().to_vec();
    let uniform_analysis = uniform.prepare_analysis(budget(), &gate).unwrap();
    let regional_analysis = regional.prepare_analysis(budget(), &gate).unwrap();
    let direct = PreparedRationalModalAnalysis::new_regional(
        &profile(),
        &gases,
        regional.boundary(),
        budget(),
        &gate,
    )
    .unwrap();
    for omega in [6000.0, 19000.0] {
        let drive = [C64::new(1.0, 0.2), C64::new(-0.4, 0.3)];
        let a = uniform_analysis
            .factor(omega, &gate)
            .unwrap()
            .response(&drive, &gate)
            .unwrap();
        let b = regional_analysis
            .factor(omega, &gate)
            .unwrap()
            .response(&drive, &gate)
            .unwrap();
        let c = direct
            .factor(omega, &gate)
            .unwrap()
            .response(&drive, &gate)
            .unwrap();
        assert_eq!(a.state, b.state);
        assert_eq!(b.state, c.state);
        assert_eq!(a.reflected_pressure_pa, b.reflected_pressure_pa);
        assert_eq!(a.load_power_w, b.load_power_w);
        check_power(&b);
    }
    assert_eq!(regional.state(), before);
}

// Independent piecewise cylindrical pressure/flow transfer, including local
// propagating/evanescent cutoffs and the complete mutual mouth impedance.
fn layered_reflection(
    model: &MimoModel,
    layers: &[(GasState, f64)],
    gammas: &[f64],
    omega: f64,
    incident: &[C64],
) -> Vec<C64> {
    let p = model.ports;
    let area = PI * RADIUS * RADIUS;
    let inlet_z = layers[0].0.density * layers[0].0.sound_speed / area;
    let mouth = model.eval(C64::new(0.0, -omega));
    let mut pressure_map = vec![C64::ZERO; p * p];
    let mut lhs = vec![C64::ZERO; p * p];
    for k in 0..p {
        let (mut a, mut b, mut c, mut d) = (C64::ONE, C64::ZERO, C64::ZERO, C64::ONE);
        for &(gas, length) in layers.iter().rev() {
            let beta =
                C64::from_re((omega / gas.sound_speed).powi(2) - (gammas[k] / RADIUS).powi(2))
                    .sqrt();
            let z = C64::from_re(gas.density * omega / area) / beta;
            let (sin, cos) = sinusoid(beta.scale(length));
            let upper = C64::new(0.0, -1.0) * z * sin;
            let lower = C64::new(0.0, -1.0) * sin / z;
            (a, b, c, d) = (
                cos * a + upper * c,
                cos * b + upper * d,
                lower * a + cos * c,
                lower * b + cos * d,
            );
        }
        for j in 0..p {
            let identity = if k == j { C64::ONE } else { C64::ZERO };
            let pressure = a * mouth[k * p + j] + b * identity;
            let flow = c * mouth[k * p + j] + d * identity;
            pressure_map[k * p + j] = pressure;
            lhs[k * p + j] = pressure + flow.scale(inlet_z);
        }
    }
    let mut outward: Vec<_> = incident.iter().map(|a| a.scale(2.0)).collect();
    lu_complex(&lhs, p).unwrap().solve(&mut outward);
    (0..p)
        .map(|k| {
            (0..p).fold(C64::ZERO, |sum, j| {
                sum + pressure_map[k * p + j] * outward[j]
            }) - incident[k]
        })
        .collect()
}

// Independent second-order original weak equation in [q; U_mouth]. It never
// uses the analysis generator, exact stepper, or an inverse bare mouth Z.
fn weak_reflection(
    prepared: &PreparedRationalModalAnalysis<'_>,
    omega: f64,
    incident: &[C64],
) -> Vec<C64> {
    let source = prepared.source();
    let (n, p) = (source.coordinates(), incident.len());
    let dimension = n + p;
    let s = C64::new(0.0, -omega);
    let mouth = prepared.boundary().load().source_model().eval(s);
    let mut lhs = vec![C64::ZERO; dimension * dimension];
    let mut rhs = vec![C64::ZERO; dimension];
    for i in 0..n {
        for j in 0..n {
            lhs[i * dimension + j] = C64::from_re(
                source.stiffness()[i * n + j] - omega * omega * source.mass()[i * n + j],
            );
        }
    }
    for k in 0..p {
        lhs[k * dimension + k] =
            lhs[k * dimension + k] + s.scale(1.0 / source.inlet_reference_impedance());
        rhs[k] = incident[k].scale(2.0 / source.inlet_reference_impedance());
        let node = n - p + k;
        lhs[node * dimension + n + k] = C64::ONE;
        lhs[(n + k) * dimension + node] = s;
        for j in 0..p {
            lhs[(n + k) * dimension + n + j] = mouth[k * p + j].scale(-1.0);
        }
    }
    lu_complex(&lhs, dimension).unwrap().solve(&mut rhs);
    (0..p).map(|k| s * rhs[k] - incident[k]).collect()
}

#[test]
fn hot_cold_modes_cross_local_cutoff_and_converge_to_layered_transfer() {
    let gate = CancelGate::new();
    let incident = [C64::new(1.0, 0.2), C64::new(-0.4, 0.3)];
    let layers = [(air(), LENGTH / 2.0), (hot_air(), LENGTH / 2.0)];
    let average = GasState::try_new(&GasSpec::dry_air_ussa1976(), 393.15, 101_325.0).unwrap();
    for direct in [false, true] {
        // The direct case has a genuine rank-one D and full mutual memory.
        let model = source(2, direct);
        let boundary = load(&model, 1e-5, false);
        for omega in [6000.0, 19000.0] {
            let mut errors = Vec::new();
            for elements in [4, 8, 16] {
                let gases = regions(elements);
                let prepared = PreparedRationalModalAnalysis::new_regional(
                    &profile(),
                    &gases,
                    &boundary,
                    budget(),
                    &gate,
                )
                .unwrap();
                assert_eq!(prepared.source().regional_gases(), Some(gases.as_slice()));
                for (element, gas) in gases.iter().enumerate() {
                    assert_eq!(prepared.source().gas_at_element(element), Some(gas));
                }
                let gammas = prepared.source().radial_eigenvalues();
                assert!(gammas[1] * air().sound_speed / RADIUS < 19000.0);
                assert!(gammas[1] * hot_air().sound_speed / RADIUS > 19000.0);
                close_real(
                    prepared.source().inlet_reference_impedance(),
                    reference(),
                    reference(),
                    1e-14,
                );
                let mouth_z = hot_air().density * hot_air().sound_speed / (PI * RADIUS * RADIUS);
                close_real(
                    prepared.source().mouth_reference_impedance(),
                    mouth_z,
                    mouth_z,
                    1e-14,
                );
                let response = prepared
                    .factor(omega, &gate)
                    .unwrap()
                    .response(&incident, &gate)
                    .unwrap();
                let weak = weak_reflection(&prepared, omega, &incident);
                let exact = layered_reflection(&model, &layers, gammas, omega, &incident);
                let averaged =
                    layered_reflection(&model, &[(average, LENGTH)], gammas, omega, &incident);
                for k in 0..2 {
                    close_complex(response.reflected_pressure_pa[k], weak[k], 1e-8);
                }
                errors.push(
                    response
                        .reflected_pressure_pa
                        .iter()
                        .zip(&exact)
                        .map(|(a, b)| (*a - *b).abs())
                        .fold(0.0_f64, f64::max),
                );
                assert!(
                    exact
                        .iter()
                        .zip(&averaged)
                        .map(|(a, b)| (*a - *b).abs())
                        .fold(0.0_f64, f64::max)
                        > 0.05
                );
                check_power(&response);
            }
            assert!(
                errors[1] < 0.4 * errors[0] && errors[2] < 0.4 * errors[1],
                "{errors:?}"
            );
            assert!(errors[2] < 0.013, "{errors:?}");
        }
    }
}

fn backward_real(gas: GasState, length: f64, s: f64, pressure: f64, flow: f64) -> (f64, f64) {
    let z = gas.density * gas.sound_speed / (PI * RADIUS * RADIUS);
    let phase = s * length / gas.sound_speed;
    let (cosh, sinh) = (phase.cosh(), phase.sinh());
    (
        cosh * pressure + z * sinh * flow,
        sinh * pressure / z + cosh * flow,
    )
}

// Multiply away the original one-pole impedance denominator before the
// independent physical transfer. No bare-component pole is mistaken for a root.
fn layered_decay_pole(model: &MimoModel) -> f64 {
    let equation = |s: f64| {
        let pressure = model.d[0] * (s + 400.0) + model.residues[0][0].re;
        let flow = s + 400.0;
        let (pressure, flow) = backward_real(hot_air(), LENGTH / 2.0, s, pressure, flow);
        let (pressure, flow) = backward_real(air(), LENGTH / 2.0, s, pressure, flow);
        pressure + reference() * flow
    };
    let (mut left, mut right) = (-700.0, -410.0);
    let mut at_left = equation(left);
    assert!(at_left * equation(right) < 0.0);
    for _ in 0..80 {
        let middle = 0.5 * (left + right);
        let value = equation(middle);
        if at_left * value <= 0.0 {
            right = middle;
        } else {
            left = middle;
            at_left = value;
        }
    }
    0.5 * (left + right)
}

fn layered_pole_seed(prepared: &PreparedRationalModalAnalysis<'_>, s: f64) -> Vec<C64> {
    let source = prepared.source();
    let gases = source.regional_gases().unwrap();
    let coordinates = source.coordinates();
    let mut pressure = vec![0.0; coordinates];
    let mut p = prepared
        .boundary()
        .load()
        .source_model()
        .eval(C64::from_re(s))[0]
        .re;
    let mut flow = 1.0;
    pressure[coordinates - 1] = p;
    for element in (0..gases.len()).rev() {
        (p, flow) = backward_real(gases[element], LENGTH / gases.len() as f64, s, p, flow);
        pressure[element] = p;
    }
    let mut seed = vec![C64::ZERO; prepared.dimension()];
    for i in 0..coordinates {
        seed[i] = C64::from_re(pressure[i] / s);
        seed[coordinates + i] = C64::from_re(
            (0..coordinates)
                .map(|j| source.mass()[i * coordinates + j] * pressure[j])
                .sum(),
        );
    }
    let load = prepared
        .boundary()
        .load()
        .continuous_operator()
        .state_space(1 << 20, &mut || false)
        .unwrap();
    assert_eq!(load.n, 1);
    seed[prepared.duct_dimension()] = C64::from_re(load.b[0] / (s - load.a[0]));
    seed[prepared.dimension() - 1] = C64::ONE;
    seed
}

#[test]
fn regional_joint_decay_pole_converges_with_original_load_state_and_power() {
    let model = source(1, true);
    let boundary = load(&model, 1e-5, false);
    let gate = CancelGate::new();
    let expected = layered_decay_pole(&model);
    let uniform = -400.0 - model.residues[0][0].re / (model.d[0] + reference());
    assert!((expected - uniform).abs() > 1.0);
    let mut errors = Vec::new();
    for elements in [2, 4, 8] {
        let prepared = PreparedRationalModalAnalysis::new_regional(
            &profile(),
            &regions(elements),
            &boundary,
            budget(),
            &gate,
        )
        .unwrap();
        let seed_s = expected * 1.001;
        let seed = layered_pole_seed(&prepared, seed_s);
        let pole = prepared
            .resonance_pole(C64::from_re(seed_s), &seed, &poles(), &gate)
            .unwrap();
        assert!(pole.laplace_pole_per_s.im.abs() < 1e-8);
        assert!(pole.laplace_pole_per_s.re < 0.0);
        assert_eq!(pole.physical.state.len(), prepared.dimension());
        assert!(pole.physical.state[prepared.duct_dimension()].abs() > 0.0);
        errors.push((pole.laplace_pole_per_s.re - expected).abs());
        check_power(&pole.physical);
    }
    assert!(
        errors[1] < 0.4 * errors[0] && errors[2] < 0.4 * errors[1],
        "{errors:?}"
    );
    assert!(errors[2] < 0.0002, "{errors:?}");
}

#[test]
fn regional_playback_preserves_work_and_every_coordinate_under_clock_subdivision() {
    let model = runtime_source();
    let gate = CancelGate::new();
    let mut gases = regions(4);
    gases.reverse(); // The actual inlet reference now differs from the exterior.
    let dt = 1e-5;
    let mut coarse = RationalModalLine::new_regional(
        &profile(),
        &gases,
        load(&model, dt, false),
        runtime_budget(),
        &gate,
    )
    .unwrap();
    let mut half = RationalModalLine::new_regional(
        &profile(),
        &gases,
        load(&model, dt / 2.0, false),
        runtime_budget(),
        &gate,
    )
    .unwrap();
    let z = hot_air().density * hot_air().sound_speed / (PI * RADIUS * RADIUS);
    close_real(coarse.source().inlet_reference_impedance(), z, z, 1e-14);
    assert_eq!(coarse.source().regional_gases(), Some(gases.as_slice()));
    let (mut inlet, mut mouth) = (0.0, 0.0);
    for step in 0..64 {
        let wave = if step < 24 {
            [
                30.0 * (0.19 * step as f64).sin(),
                18.0 * (0.23 * step as f64).cos(),
            ]
        } else {
            [0.0; 2]
        };
        let record = coarse.push(&wave, &gate).unwrap();
        let first = half.push(&wave, &gate).unwrap();
        let second = half.push(&wave, &gate).unwrap();
        let scale = record.incident_energy_j
            + record.reflected_energy_j
            + coarse.stored_duct_energy_j()
            + record.terminal_work_j.abs();
        close_real(
            record.incident_energy_j,
            dt * (wave[0] * wave[0] + wave[1] * wave[1]) / z,
            scale,
            1e-12,
        );
        close_real(
            record.terminal_work_j,
            first.terminal_work_j + second.terminal_work_j,
            scale,
            3e-8,
        );
        close_real(
            record.reflected_energy_j,
            first.reflected_energy_j + second.reflected_energy_j,
            scale,
            3e-8,
        );
        let coordinates = coarse.source().coordinates();
        for range in [
            0..coordinates,
            coordinates..2 * coordinates,
            2 * coordinates..coarse.state().len(),
        ] {
            let state_scale = coarse.state()[range.clone()]
                .iter()
                .chain(&half.state()[range.clone()])
                .map(|x| x.abs())
                .fold(0.0_f64, f64::max);
            for index in range {
                close_real(
                    coarse.state()[index],
                    half.state()[index],
                    state_scale,
                    3e-8,
                );
            }
        }
        for account in coarse.terminal_accounts() {
            close_real(
                account.mean_duct_pressure_pa,
                account.mean_load_pressure_pa,
                60.0,
                2e-9,
            );
            close_real(account.duct_work_j, account.load_work_j, scale, 3e-8);
        }
        inlet += record.inlet_work_j;
        mouth += record.terminal_work_j;
        close_real(
            coarse.stored_duct_energy_j(),
            inlet - mouth,
            inlet.abs() + mouth.abs(),
            3e-8,
        );
        assert_eq!(
            coarse.load_state(),
            coarse.load().continuous_operator().real_coordinates()
        );
    }
    assert_eq!(coarse.duration(), dt);
    assert_eq!(half.duration(), dt / 2.0);
    assert_eq!(coarse.load().source_model(), &model);
}

#[test]
fn regional_pressure_capacity_band_and_cancellation_refusals_preserve_sources() {
    let model = runtime_source();
    let boundary = load(&model, 1e-5, true);
    let gate = CancelGate::new();
    let accepted = boundary.load().pole_states().to_vec();
    let mut different_pressure = regions(2);
    different_pressure[1] =
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 493.15, 90_000.0).unwrap();
    let mut incoherent = regions(2);
    incoherent[1].density *= 1.2;
    for gases in [vec![], different_pressure, incoherent] {
        assert!(
            PreparedRationalModalAnalysis::new_regional(
                &profile(),
                &gases,
                &boundary,
                budget(),
                &gate,
            )
            .is_err()
        );
        assert!(
            RationalModalLine::new_regional(
                &profile(),
                &gases,
                load(&model, 1e-5, true),
                runtime_budget(),
                &gate,
            )
            .is_err()
        );
    }
    let gases = regions(4);
    for limits in [
        RationalModalAnalysisBudget {
            maximum_dimension: 1,
            ..budget()
        },
        RationalModalAnalysisBudget {
            maximum_payload_bytes: 1,
            ..budget()
        },
        RationalModalAnalysisBudget {
            maximum_work: 1,
            ..budget()
        },
    ] {
        assert!(
            PreparedRationalModalAnalysis::new_regional(
                &profile(),
                &gases,
                &boundary,
                limits,
                &gate,
            )
            .is_err()
        );
    }
    let mut limits = runtime_budget();
    limits.exact.max_multiply_terms = 1;
    assert!(
        RationalModalLine::new_regional(
            &profile(),
            &gases,
            load(&model, 1e-5, true),
            limits,
            &gate,
        )
        .is_err()
    );
    let prepared =
        PreparedRationalModalAnalysis::new_regional(&profile(), &gases, &boundary, budget(), &gate)
            .unwrap();
    assert!(prepared.factor(0.5, &gate).is_err());
    assert!(prepared.factor(4e6 + 1.0, &gate).is_err());
    gate.request();
    assert!(matches!(
        PreparedRationalModalAnalysis::new_regional(&profile(), &[], &boundary, budget(), &gate,),
        Err(RationalModalAnalysisError::Cancelled)
    ));
    assert!(matches!(
        RationalModalLine::new_regional(
            &profile(),
            &[],
            load(&model, 1e-5, true),
            runtime_budget(),
            &gate,
        ),
        Err(RationalModalError::Cancelled)
    ));
    assert_eq!(boundary.load().pole_states(), accepted);
}
