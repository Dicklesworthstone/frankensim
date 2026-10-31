use super::*;
use crate::{LossModel, input_impedance_load};
use fs_evidence::ValidityDomain;
use fs_gas::dynamics::CaloricDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;
use fs_qty::{Length, MolarMass, Pressure, Temperature};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{
    Composition, CompositionBasis, ElementalReferenceIdV1, Nasa9RegionV1,
    Nasa9StandardStateModelV1, ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1,
    StandardStatePhaseV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K, required_nasa9_parameter_dims,
};
use fsci_runtime::RuntimeMode;
use fsci_special::SpecialTensor;

const STEP_S: f64 = 1.0 / 512_000.0;

fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}

// A manufactured constant-Cp NASA-9 card with exactly the material state's R
// and Cp, so the field and the transfer-matrix chain share density, sound
// speed and gamma. An analytic fixture, not measured gas data.
fn air_source(state: &GasState) -> FrozenGas {
    caloric_source(state, 0.0)
}

// Cp(T)/R = Cp(0)/R + slope*T, with an independently integrable Cv.
fn caloric_source(state: &GasState, slope: f64) -> FrozenGas {
    let mut a = [0.0; 9];
    a[2] = state.specific_heat_cp / state.specific_gas_constant;
    a[3] = slope;
    let mut parameters = std::collections::BTreeMap::new();
    for ((name, dims), value) in required_nasa9_parameter_dims().into_iter().take(9).zip(a) {
        parameters.insert(name.to_owned(), LawParameter { value, dims });
    }
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: 1.0,
            dims: Pressure::DIMS,
        },
    );
    let card = ConstitutiveModelCard {
        law: LawId("nasa9-standard-state".into()),
        law_version: 1,
        parameters,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 100.0, 1000.0),
        sources: Vec::new(),
        provenance: Provenance {
            source: "Manufactured caloric bore test".into(),
            license: "test fixture".into(),
            artifact: None,
        },
    };
    let convention = StandardStateConventionV1::new(
        StandardStatePhaseV1::Gas,
        ReferenceEquationOfStateV1::IdealGas,
        Pressure::new(1.0),
        ElementalReferenceIdV1::new("bore-fixture").unwrap(),
    )
    .unwrap();
    let species = Nasa9StandardStateModelV1::new(
        SpeciesId::new("bore-test-gas").unwrap(),
        MolarMass::new(UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / state.specific_gas_constant),
        convention,
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap();
    let source = FrozenIdealGasMixtureModelV1::new(
        vec![species],
        Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
    )
    .unwrap();
    FrozenGas::new(
        source,
        CaloricDomain::new(Temperature::new(100.0), Temperature::new(1000.0)),
        || false,
    )
    .unwrap()
}

fn wall_gas(state: &GasState) -> WallGas {
    WallGas {
        density_kg_m3: state.density,
        dynamic_viscosity_pa_s: state.dynamic_viscosity,
        gamma: state.gamma,
        prandtl: state.prandtl,
    }
}

// A single, exactly soluble thermal arm: gain=rate=1 and exp(-2*dt)=1/5.
// From zero history its pressure response is 3/5 of the initial perturbation.
fn caloric_walls(cells: usize, reference_pressure: f64, drag: f64, thermal: bool) -> Walls {
    let dt = 0.5 * det::ln(5.0);
    let thermal_arms: &[(f64, f64)] = if thermal { &[(1.0, 1.0)] } else { &[] };
    Walls {
        viscous_of: vec![0; cells],
        thermal_of: vec![1; cells],
        propagators: vec![
            Propagator::new(drag, &[], dt).unwrap(),
            Propagator::new(0.0, thermal_arms, dt).unwrap(),
        ],
        viscous_state: vec![0.0; 2 * cells],
        thermal_state: vec![0.0; 2 * cells],
        spare_viscous: vec![0.0; 2 * cells],
        spare_thermal: vec![0.0; 2 * cells],
        state_stride: 2,
        reference_pressure,
        receipt: WallReceipt {
            distinct_radii: 1,
            maximum_arms: usize::from(thermal),
            maximum_complex_error: 0.0,
            maximum_loss_error: 0.0,
        },
    }
}

fn caloric_field(gas: &FrozenGas, states: Vec<Conserved>) -> Euler1d {
    Euler1d::new(
        gas.clone(),
        Length::new(1.0),
        1.0,
        states,
        Boundaries::Ends {
            left: Boundary::Reflecting,
            right: Boundary::Reflecting,
        },
        || false,
    )
    .unwrap()
}

#[test]
fn g1_caloric_wall_heat_matches_the_exact_cv_integral_in_both_directions() {
    let state = air();
    let r = state.specific_gas_constant;
    let cv0 = state.specific_heat_cp - r;
    let rho = 1.7;
    for slope in [0.0, 0.003] {
        let gas = caloric_source(&state, slope);
        for initial_t in [180.0, 600.0] {
            let q = gas.conserved(rho, 0.0, rho * r * initial_t).unwrap();
            let mut solver = caloric_field(&gas, vec![q; 2]);
            let mut walls = caloric_walls(2, rho * r * 300.0, 0.0, true);
            walls.apply(&mut solver, &gas).unwrap();
            let expected_t = 300.0 + 0.6 * (initial_t - 300.0);
            let expected_heat = rho
                * (cv0 * (expected_t - initial_t)
                    + 0.5 * r * slope * (expected_t * expected_t - initial_t * initial_t));
            let next = solver.states()[0];
            let p = gas.primitive(next, || false).unwrap();
            assert_eq!(next.density, q.density);
            assert_eq!(next.momentum, q.momentum);
            assert!((p.temperature - expected_t).abs() < 1e-8);
            assert!((next.energy - q.energy - expected_heat).abs() < 1e-7);
        }
    }
}

#[test]
fn g1_caloric_wall_viscous_heat_precedes_thermal_exchange_without_frozen_gamma() {
    let state = air();
    let r = state.specific_gas_constant;
    let slope = 0.003;
    let gas = caloric_source(&state, slope);
    let (rho, initial_t, initial_u) = (1.7, 450.0, 600.0);
    let q = gas.conserved(rho, initial_u, rho * r * initial_t).unwrap();
    let next_u = initial_u / 5.0_f64.sqrt();
    let dissipated_per_mass = 0.5 * (initial_u * initial_u - next_u * next_u);
    let cv = state.specific_heat_cp - r + r * slope * initial_t;
    // Positive quadratic root of integral Cv(T) dT = lost kinetic energy.
    let delta_t =
        2.0 * dissipated_per_mass / (cv + (cv * cv + 2.0 * r * slope * dissipated_per_mass).sqrt());
    let heated_t = initial_t + delta_t;
    for thermal in [false, true] {
        let mut solver = caloric_field(&gas, vec![q; 2]);
        let mut walls = caloric_walls(2, rho * r * 300.0, 1.0, thermal);
        walls.apply(&mut solver, &gas).unwrap();
        let expected_t = if thermal {
            300.0 + 0.6 * (heated_t - 300.0)
        } else {
            heated_t
        };
        let next = solver.states()[0];
        let p = gas.primitive(next, || false).unwrap();
        let expected_heat = rho
            * ((state.specific_heat_cp - r) * (expected_t - heated_t)
                + 0.5 * r * slope * (expected_t * expected_t - heated_t * heated_t));
        assert_eq!(next.density, q.density);
        assert!((p.velocity - next_u).abs() < 1e-10);
        assert!((p.temperature - expected_t).abs() < 1e-8);
        assert!((next.energy - q.energy - expected_heat).abs() < 1e-7);
        if !thermal {
            assert_eq!(
                next.energy, q.energy,
                "viscous heat is an internal transfer"
            );
        }
    }
}

#[test]
fn g0_caloric_wall_domain_refusal_preserves_all_cells_and_memory() {
    let state = air();
    let gas = caloric_source(&state, 0.003);
    let cells: Vec<_> = [300.0, 160.0]
        .into_iter()
        .map(|t| {
            gas.conserved(1.0, 0.0, state.specific_gas_constant * t)
                .unwrap()
        })
        .collect();
    let mut solver = caloric_field(&gas, cells.clone());
    // The first cell cools to 180 K; the second would leave the admitted
    // [100, 1000] K interval at 96 K. Neither cell nor its arms may publish.
    let mut walls = caloric_walls(2, 0.0, 0.0, true);
    let viscous_before = walls.viscous_state.clone();
    let thermal_before = walls.thermal_state.clone();
    assert!(walls.apply(&mut solver, &gas).is_err());
    assert_eq!(solver.states(), cells);
    assert_eq!(solver.time().value(), 0.0);
    assert_eq!(walls.viscous_state, viscous_before);
    assert_eq!(walls.thermal_state, thermal_before);
}

fn j1(x: f64) -> f64 {
    match fsci_special::j1(&SpecialTensor::RealScalar(x), RuntimeMode::Strict).unwrap() {
        SpecialTensor::RealScalar(v) => v,
        other => panic!("non-scalar oracle: {other:?}"),
    }
}

// Rigid baffled piston, +i omega t convention: rho c / S [1 - J1(2ka)/ka + i H1(2ka)/ka].
fn piston(ka: f64) -> C64 {
    C64::new(
        1.0 - j1(2.0 * ka) / ka,
        fsci_special::struve(1.0, 2.0 * ka) / ka,
    )
}

#[test]
fn piston_weights_sum_to_one_and_reproduce_the_struve_impedance() {
    let (radius, c) = (0.061, 343.2);
    let weights = piston_weights(STEP_S, 2.0 * radius / c).unwrap();
    let total: f64 = weights.iter().sum();
    println!("piston: {} taps, sum-1 = {:e}", weights.len(), total - 1.0);
    assert!((total - 1.0).abs() < 1e-12);
    assert!(weights.iter().all(|&w| w >= 0.0));
    for f in [20.0, 100.0, 400.0, 1000.0, 2000.0, 4000.0, 6900.0] {
        let omega = core::f64::consts::TAU * f;
        let mut transfer = C64::new(0.0, 0.0);
        for (j, &w) in weights.iter().enumerate() {
            let phase = -omega * STEP_S * j as f64;
            transfer = transfer + C64::new(w * det::cos(phase), w * det::sin(phase));
        }
        let discrete = C64::new(1.0 - transfer.re, -transfer.im);
        let exact = piston(omega * radius / c);
        let error = (discrete - exact).abs();
        // Linear interpolation of a tone plus its images: (omega dt)^2 / 6.
        let bound = (omega * STEP_S).powi(2) / 6.0 + 1e-11;
        println!(
            "piston f={f:6.0} ka={:.3} exact={:.6}{:+.6}i error={error:.2e} bound={bound:.2e}",
            omega * radius / c,
            exact.re,
            exact.im
        );
        assert!(error <= bound);
    }
}

#[test]
fn step_mesh_keeps_exact_volume_and_moves_faces_off_the_step() {
    let state = air();
    let duct = Duct {
        segments: vec![
            Segment::Cone {
                inlet_radius: 1.825e-3,
                outlet_radius: 4.2e-3,
                length: 0.06597,
            },
            Segment::Cone {
                inlet_radius: 4.64e-3,
                outlet_radius: 5.0e-3,
                length: 0.055,
            },
            Segment::Cylinder {
                radius: 5.0e-3,
                length: 0.3,
            },
        ],
    };
    let bore = NonlinearBore::new(
        &duct,
        air_source(&state),
        293.15,
        101_325.0,
        None,
        BoreSpec {
            maximum_cell_m: 2.4e-3,
            inlet_volume_m3: 0.0,
            step_s: STEP_S,
            cfl: 0.45,
            wall: None,
        },
    )
    .unwrap();
    let mesh = bore.solver().mesh();
    let frustum =
        |a: f64, b: f64, l: f64| core::f64::consts::PI * l / 3.0 * (a * a + a * b + b * b);
    let exact = frustum(1.825e-3, 4.2e-3, 0.06597)
        + frustum(4.64e-3, 5.0e-3, 0.055)
        + core::f64::consts::PI * 5.0e-3 * 5.0e-3 * 0.3;
    let total: f64 = mesh.volumes().iter().sum();
    assert!((total - exact).abs() <= 1e-12 * exact, "{total} {exact}");
    let n = mesh.volumes().len();
    let h = 0.42097 / n as f64;
    let step = 0.06597;
    assert!(
        mesh.positions()
            .iter()
            .all(|&x| (x - step).abs() >= 0.25 * h * (1.0 - 1e-12))
    );
    let cell = mesh.positions().partition_point(|&x| x < step) - 1;
    let areas = mesh.face_areas();
    let a = |r: f64| core::f64::consts::PI * r * r;
    assert!(areas[cell] < a(4.2e-3) && areas[cell + 1] > a(4.64e-3));
    println!(
        "step cell {cell} of {n}, faces {:?}",
        &mesh.positions()[cell..cell + 2]
    );
}

fn multiply(a: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    let mut c = vec![0.0; n * n];
    for i in 0..n {
        for k in 0..n {
            for j in 0..n {
                c[i * n + j] += a[i * n + k] * b[k * n + j];
            }
        }
    }
    c
}

// Independent oracle: exp(G dt) by its Taylor series after scaling, then squaring.
fn taylor_exponential(generator: &[f64], n: usize, dt: f64) -> Vec<f64> {
    let norm = (0..n)
        .map(|r| {
            generator[r * n..(r + 1) * n]
                .iter()
                .map(|x| x.abs())
                .sum::<f64>()
        })
        .fold(0.0, f64::max);
    let mut squarings = 0;
    let mut step = dt;
    while norm * step > 0.5 {
        step *= 0.5;
        squarings += 1;
    }
    let g: Vec<f64> = generator.iter().map(|x| x * step).collect();
    let mut result = vec![0.0; n * n];
    let mut term = vec![0.0; n * n];
    for i in 0..n {
        result[i * n + i] = 1.0;
        term[i * n + i] = 1.0;
    }
    for order in 1..=24 {
        term = multiply(&term, &g, n);
        for x in &mut term {
            *x /= f64::from(order);
        }
        for (r, t) in result.iter_mut().zip(&term) {
            *r += t;
        }
    }
    for _ in 0..squarings {
        result = multiply(&result, &result, n);
    }
    result
}

#[test]
fn wall_propagators_are_contractive_semigroups_of_their_generator() {
    let arms = [
        (120.0, 30.0),
        (900.0, 4.0e3),
        (5.0e3, 2.0e5),
        (2.0e4, 3.0e6),
    ];
    let n = arms.len() + 1;
    let mut generator = vec![0.0; n * n];
    for r in 0..n {
        for c in 0..n {
            generator[r * n + c] = match (r, c) {
                (0, 0) => -(4.8 + arms.iter().map(|a| a.0).sum::<f64>()),
                (0, k) | (k, 0) => det::sqrt(arms[k - 1].0 * arms[k - 1].1),
                (k, l) if k == l => -arms[k - 1].1,
                _ => 0.0,
            };
        }
    }
    // The eigen propagator is the exponential of this test's own generator
    // (the Taylor oracle) at the bore's step, twice it, and a vanishing step.
    for dt in [STEP_S, 2.0 * STEP_S, 1e-12] {
        let propagator = Propagator::new(4.8, &arms, dt).unwrap();
        let oracle = taylor_exponential(&generator, n, dt);
        for (k, (a, b)) in propagator.matrix().iter().zip(&oracle).enumerate() {
            assert!((a - b).abs() <= 1e-13, "dt {dt:e} entry {k}: {a} {b}");
        }
    }
    let one = Propagator::new(4.8, &arms, STEP_S).unwrap();
    let two = Propagator::new(4.8, &arms, 2.0 * STEP_S).unwrap();
    let squared = multiply(&one.matrix(), &one.matrix(), n);
    for (a, b) in two.matrix().iter().zip(&squared) {
        assert!((a - b).abs() <= 1e-13, "{a} {b}");
    }
    // One advance with a replaced mean is the dense map of [mean; arms].
    let mut state = vec![1.0, -0.4, 0.3, 2.0, -1.0];
    let dense = |y: &[f64]| -> Vec<f64> {
        (0..n)
            .map(|i| (0..n).map(|k| one.basis[i * n + k] * y[k]).sum())
            .collect()
    };
    let mut x = dense(&state);
    x[0] = 0.7;
    let expected: Vec<f64> = (0..n)
        .map(|i| (0..n).map(|j| one.matrix()[i * n + j] * x[j]).sum())
        .collect();
    let (mean, before, after) = one.advance(&mut state, 0.7);
    let actual = dense(&state);
    for (a, b) in actual.iter().zip(&expected) {
        assert!((a - b).abs() <= 1e-13, "{a} {b}");
    }
    assert!((mean - expected[0]).abs() <= 1e-13);
    assert!((before - x.iter().map(|v| v * v).sum::<f64>()).abs() <= 1e-13);
    assert!((after - expected.iter().map(|v| v * v).sum::<f64>()).abs() <= 1e-13);
    // Held at its own mean, the state's energy never grows.
    let mut energy = after;
    for _ in 0..2000 {
        let current: f64 = one.head.iter().zip(&state).map(|(h, y)| h * y).sum();
        let (_, before, after) = one.advance(&mut state, current);
        assert!(before <= energy * (1.0 + 1e-14));
        assert!(after <= before * (1.0 + 1e-14));
        energy = after;
    }
}

#[test]
fn a_bore_refuses_a_step_its_field_cannot_take() {
    let state = air();
    let duct = Duct {
        segments: vec![Segment::Cylinder {
            radius: 5.0e-3,
            length: 0.2,
        }],
    };
    let mut bore = NonlinearBore::new(
        &duct,
        air_source(&state),
        293.15,
        101_325.0,
        None,
        BoreSpec {
            maximum_cell_m: 2.4e-3,
            inlet_volume_m3: 0.0,
            step_s: 1.0e-5,
            cfl: 0.45,
            wall: None,
        },
    )
    .unwrap();
    let error = bore.step(|_| Ok(None)).unwrap_err().to_string();
    assert!(error.contains("fixed step"), "{error}");
}

// A prescribed end volume velocity: the interior's own end wave supplies the pressure.
fn velocity_port(
    gas: &FrozenGas,
    stage: &PortStage,
    velocity: f64,
    guess: &mut f64,
) -> Result<(FaceFlux, f64), GasError> {
    let inside = stage.left.1;
    let mut pressure = *guess;
    for _ in 0..64 {
        let wave = gas.end_wave(inside, -1.0, pressure, || false)?;
        let next = pressure - (wave.primitive.velocity - velocity) / wave.velocity_per_pressure;
        let done = (next - pressure).abs() <= 1e-13 * pressure;
        pressure = next;
        if done {
            break;
        }
    }
    let wave = gas.end_wave(inside, -1.0, pressure, || false)?;
    *guess = pressure;
    let face = Primitive {
        velocity,
        ..wave.primitive
    };
    Ok((gas.face_flux(face, wave.signal_speed)?, pressure))
}

fn dft(signal: &[f64], omega_dt: f64) -> C64 {
    let (mut re, mut im) = (0.0, 0.0);
    let (s, c) = (det::sin(omega_dt), det::cos(omega_dt));
    let (mut pr, mut pi) = (1.0, 0.0);
    for &x in signal {
        re += x * pr;
        im -= x * pi;
        let next = pr * c - pi * s;
        pi = pr * s + pi * c;
        pr = next;
    }
    C64::new(re, im)
}

fn parabolic_peak(f: &[f64], m: &[f64], i: usize) -> (f64, f64) {
    let (a, b, c) = (m[i - 1], m[i], m[i + 1]);
    let offset = 0.5 * (a - c) / (a - 2.0 * b + c);
    let x = det::ln(f[i]) + offset * (det::ln(f[i + 1]) - det::ln(f[i]));
    (det::exp(x), b - 0.25 * (a - c) * offset)
}

#[test]
fn g3_small_amplitude_bore_matches_the_transfer_matrix_input_impedance() {
    let state = air();
    let duct = Duct {
        segments: vec![
            Segment::Cylinder {
                radius: 4.0e-3,
                length: 0.3,
            },
            Segment::Cone {
                inlet_radius: 4.0e-3,
                outlet_radius: 15.0e-3,
                length: 0.1,
            },
        ],
    };
    let mut bore = NonlinearBore::new(
        &duct,
        air_source(&state),
        293.15,
        101_325.0,
        Some(wall_gas(&state)),
        BoreSpec {
            maximum_cell_m: 2.4e-3,
            inlet_volume_m3: 0.0,
            step_s: STEP_S,
            cfl: 0.45,
            wall: Some(WallMemory {
                // 48 poles: the 15 mm bell end is admitted on 20-8000 Hz
                // (fs-phs measures 7.5e-4); 24 resolve only up to about 6 mm.
                poles: 48,
                band_hz: [20.0, 8000.0],
                maximum_relative_error: 1e-2,
            }),
        },
    )
    .unwrap();
    let room = bore.room();
    assert!((room.sound_speed / state.sound_speed - 1.0).abs() < 1e-9);
    assert!((room.density / state.density - 1.0).abs() < 1e-9);
    println!("g3 wall {:?}", bore.wall_receipt());
    let source = bore.gas().clone();
    let inlet_area = bore.inlet_area_m2();
    let samples = (0.25 / STEP_S) as usize;
    let width = 1.0e-4;
    let peak_flow = 1.0e-8;
    let flow = |t: f64| {
        if t < width {
            peak_flow * 0.5 * (1.0 - det::cos(core::f64::consts::TAU * t / width))
        } else {
            0.0
        }
    };
    let mut pressure = Vec::with_capacity(samples);
    let mut volume = Vec::with_capacity(samples);
    let mut guess = room.pressure;
    for n in 0..samples {
        let mut face = room.pressure;
        bore.step(|stage| {
            let (flux, p) =
                velocity_port(&source, stage, flow(stage.time) / inlet_area, &mut guess)?;
            if stage.stage == 0 {
                face = p;
            }
            Ok(Some(flux))
        })
        .unwrap();
        pressure.push(face - room.pressure);
        volume.push(flow(n as f64 * STEP_S));
    }
    let zc = state.density * state.sound_speed / (core::f64::consts::PI * 15.0e-3 * 15.0e-3);
    let reference = |f: f64| {
        let omega = core::f64::consts::TAU * f;
        // The transfer-matrix chain is e^{-i omega t}: the conjugate load.
        let load = piston(omega * 15.0e-3 / state.sound_speed);
        let z = input_impedance_load(
            &duct,
            &state,
            omega,
            LossModel::Bessel,
            C64::new(zc * load.re, -zc * load.im),
            None,
        )
        .unwrap();
        z.impedance.abs()
    };
    let grid: Vec<f64> = (0..3000).map(|i| 50.0 + i as f64).collect();
    let tmm: Vec<f64> = grid.iter().map(|&f| reference(f)).collect();
    let mut peaks = Vec::new();
    for i in 1..grid.len() - 1 {
        if tmm[i] > tmm[i - 1] && tmm[i] >= tmm[i + 1] && peaks.len() < 6 {
            let fine: Vec<f64> = (-20..=20).map(|k| grid[i] + 0.05 * f64::from(k)).collect();
            let values: Vec<f64> = fine
                .iter()
                .map(|&f| 20.0 * det::ln(reference(f)) / core::f64::consts::LN_10)
                .collect();
            let best = (1..fine.len() - 1)
                .max_by(|&a, &b| values[a].total_cmp(&values[b]))
                .unwrap();
            peaks.push(parabolic_peak(&fine, &values, best));
        }
    }
    assert_eq!(peaks.len(), 6, "transfer-matrix peaks {peaks:?}");
    let fv = |f: f64| {
        let omega_dt = core::f64::consts::TAU * f * STEP_S;
        let z = dft(&pressure, omega_dt) / dft(&volume, omega_dt);
        20.0 * det::ln(z.abs()) / core::f64::consts::LN_10
    };
    for (k, &(f_ref, db_ref)) in peaks.iter().enumerate() {
        let fine: Vec<f64> = (-200..=200)
            .map(|c| f_ref * det::exp(f64::from(c) * 0.5 / 1200.0 * core::f64::consts::LN_2))
            .collect();
        let values: Vec<f64> = fine.iter().map(|&f| fv(f)).collect();
        let best = (1..fine.len() - 1)
            .max_by(|&a, &b| values[a].total_cmp(&values[b]))
            .unwrap();
        let (f_fv, db_fv) = parabolic_peak(&fine, &values, best);
        let cents = 1200.0 * det::ln(f_fv / f_ref) / core::f64::consts::LN_2;
        println!(
            "g3 peak {}: tmm {f_ref:8.2} Hz {db_ref:6.2} dB | bore {f_fv:8.2} Hz {db_fv:6.2} dB | {cents:+6.2} cents {:+5.2} dB",
            k + 1,
            db_fv - db_ref
        );
        assert!(cents.abs() <= 10.0 && (db_fv - db_ref).abs() <= 1.5);
    }
}

#[test]
fn g2_sinusoid_steepens_to_the_lossless_shock_distance() {
    let state = air();
    let duct = Duct {
        segments: vec![Segment::Cylinder {
            radius: 10.0e-3,
            length: 3.0,
        }],
    };
    let step = 1.0e-6;
    let mut bore = NonlinearBore::new(
        &duct,
        air_source(&state),
        293.15,
        101_325.0,
        None,
        BoreSpec {
            maximum_cell_m: 1.2e-3,
            inlet_volume_m3: 0.0,
            step_s: step,
            cfl: 0.45,
            wall: None,
        },
    )
    .unwrap();
    let room = bore.room();
    let source = bore.gas().clone();
    let (frequency, amplitude) = (1000.0, 2000.0);
    let omega = core::f64::consts::TAU * frequency;
    let velocity = amplitude / (room.density * room.sound_speed);
    let centres = bore.centers_m().to_vec();
    let near = centres.partition_point(|&x| x < 0.05);
    let far = centres.partition_point(|&x| x < 1.3);
    let (x0, x1) = (centres[near], centres[far]);
    let mut series = [Vec::new(), Vec::new()];
    let mut guess = room.pressure;
    let end = 0.0105;
    let steps = (end / step) as usize;
    for _ in 0..steps {
        bore.step(|stage| {
            let u = velocity * det::sin(omega * stage.time);
            Ok(Some(velocity_port(&source, stage, u, &mut guess)?.0))
        })
        .unwrap();
        for (k, &cell) in [near, far].iter().enumerate() {
            let q = bore.solver().states()[cell];
            series[k].push(source.primitive(q, || false).unwrap().pressure - room.pressure);
        }
    }
    // The second full period to pass each point: away from the switch-on.
    let slope = |k: usize, x: f64| {
        let start = ((x / room.sound_speed + 1.0 / frequency) / step) as usize;
        let stop = start + (1.0 / frequency / step) as usize;
        let window = &series[k][start..stop];
        let peak = window
            .windows(3)
            .map(|w| (w[2] - w[0]) / (2.0 * step))
            .fold(0.0, f64::max);
        let swing = 0.5
            * (window.iter().copied().fold(f64::MIN, f64::max)
                - window.iter().copied().fold(f64::MAX, f64::min));
        (peak, swing)
    };
    let (s0, p0) = slope(0, x0);
    let (s1, _) = slope(1, x1);
    let ratio = s1 / s0;
    let gamma = state.gamma;
    let beta = 0.5 * (gamma + 1.0);
    let theory = room.density * room.sound_speed.powi(3) / (beta * omega * p0);
    let measured = (x0 - ratio * x1) / (1.0 - ratio);
    println!(
        "g2 amplitude {p0:.1} Pa slope ratio {ratio:.4}: shock distance {measured:.4} m vs {theory:.4} m ({:+.2} %)",
        100.0 * (measured / theory - 1.0)
    );
    assert!((measured / theory - 1.0).abs() <= 0.05);
}
