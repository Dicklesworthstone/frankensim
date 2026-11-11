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

fn reference_bore(
    gas: &FrozenGas,
    temperature: f64,
    pressure: f64,
    transport: WallGas,
) -> Result<NonlinearBore, NonlinearDuctError> {
    NonlinearBore::new(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: 1.0e-3,
                length: 0.02,
            }],
        },
        gas.clone(),
        temperature,
        pressure,
        Some(transport),
        BoreSpec {
            maximum_cell_m: 0.01,
            inlet_volume_m3: 0.0,
            step_s: 1.0e-7,
            cfl: 0.45,
            wall: Some(WallMemory {
                poles: 24,
                band_hz: [20.0, 8000.0],
                maximum_relative_error: 1e-2,
            }),
        },
    )
}

#[test]
fn g0_bore_refuses_stale_wall_gamma_and_density() {
    let state = air();
    let slope = 0.003;
    let gas = caloric_source(&state, slope);
    let (temperature, pressure) = (600.0, 150_000.0);
    let cp = state.specific_heat_cp + state.specific_gas_constant * slope * temperature;
    let transport = WallGas {
        density_kg_m3: pressure / (state.specific_gas_constant * temperature),
        dynamic_viscosity_pa_s: state.dynamic_viscosity,
        gamma: cp / (cp - state.specific_gas_constant),
        prandtl: state.dynamic_viscosity * cp / state.thermal_conductivity,
    };
    for (coefficient, stale) in [
        (
            "gamma",
            WallGas {
                gamma: state.gamma,
                ..transport
            },
        ),
        (
            "density",
            WallGas {
                density_kg_m3: 1.01 * transport.density_kg_m3,
                ..transport
            },
        ),
        (
            "gamma",
            WallGas {
                gamma: f64::NAN,
                ..transport
            },
        ),
        (
            "density",
            WallGas {
                density_kg_m3: f64::INFINITY,
                ..transport
            },
        ),
    ] {
        let error = reference_bore(&gas, temperature, pressure, stale)
            .err()
            .expect("inconsistent wall snapshot must refuse")
            .to_string();
        assert!(
            error.contains(coefficient) && error.contains("room state"),
            "{error}"
        );
    }
    // An independent analytic snapshot at the actual reference remains valid.
    reference_bore(&gas, temperature, pressure, transport).unwrap();
}

#[test]
fn g0_source_wall_gas_retains_the_field_domain_and_transport_admission() {
    let state = air();
    let source = caloric_source(&state, 0.003);
    let gas = FrozenGas::new(
        source.source().clone(),
        CaloricDomain::new(Temperature::new(200.0), Temperature::new(700.0)),
        || false,
    )
    .unwrap();
    for temperature in [150.0, 750.0, f64::NAN] {
        assert!(matches!(
            WallGas::from_source(&gas, temperature, 100_000.0, 2e-5, 0.03),
            Err(NonlinearDuctError::Gas(GasError::CaloricPath(_)))
        ));
    }
    for (viscosity, conductivity) in [(0.0, 0.03), (2e-5, 0.0)] {
        assert!(WallGas::from_source(&gas, 300.0, 100_000.0, viscosity, conductivity).is_err());
    }
}

#[test]
fn g1_source_wall_gas_uses_local_cp_in_gamma_prandtl_and_compliance() {
    let state = air();
    let r = state.specific_gas_constant;
    for slope in [0.0, 0.003] {
        let gas = caloric_source(&state, slope);
        for temperature in [250.0, 650.0] {
            let (pressure, viscosity, conductivity) = (175_000.0, 2e-5, 0.03);
            let cp = state.specific_heat_cp + r * slope * temperature;
            let transport =
                WallGas::from_source(&gas, temperature, pressure, viscosity, conductivity).unwrap();
            assert!((transport.density_kg_m3 / (pressure / (r * temperature)) - 1.0).abs() < 1e-12);
            assert!((transport.gamma / (cp / (cp - r)) - 1.0).abs() < 1e-12);
            assert!((transport.prandtl / (viscosity * cp / conductivity) - 1.0).abs() < 1e-12);

            let mut bore = reference_bore(&gas, temperature, pressure, transport).unwrap();
            let walls = bore.walls.as_ref().unwrap();
            assert_eq!(walls.wall_temperature, temperature);
            assert_wall_compliance(walls, transport, 1.0e-3, bore.step_s);
            // The source-bound wall and radiation also retain a quiescent
            // nonstandard-temperature room through actual accepted steps.
            let initial = bore.solver.states().to_vec();
            for _ in 0..3 {
                bore.step(|_| Ok(None)).unwrap();
            }
            for (before, after) in initial.iter().zip(bore.solver.states()) {
                assert!((after.density / before.density - 1.0).abs() < 1e-12);
                assert!((after.energy / before.energy - 1.0).abs() < 1e-12);
                assert!(
                    (gas.primitive(*after, || false).unwrap().temperature - temperature).abs()
                        < 1e-8
                );
            }
        }
    }
}

// A single, exactly soluble thermal arm: gain=rate=1 and exp(-2*dt)=1/5.
// From zero history its temperature response is 3/5 of the initial perturbation.
fn caloric_walls(cells: usize, wall_temperature: f64, drag: f64, thermal: bool) -> Walls {
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
        wall_temperature,
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
            let mut walls = caloric_walls(2, 300.0, 0.0, true);
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
        let mut walls = caloric_walls(2, 300.0, 1.0, thermal);
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
fn g1_caloric_wall_isothermal_density_variations_do_not_create_heat() {
    let state = air();
    let wall_temperature = 300.0;
    for slope in [0.0, 0.003] {
        let gas = caloric_source(&state, slope);
        let cells: Vec<_> = [0.2, 1.7, 7.0]
            .into_iter()
            .map(|rho| {
                gas.conserved(rho, 0.0, rho * gas.gas_constant() * wall_temperature)
                    .unwrap()
            })
            .collect();
        let mut solver = caloric_field(&gas, cells.clone());
        let mut walls = caloric_walls(cells.len(), wall_temperature, 0.0, true);
        for _ in 0..4 {
            walls.apply(&mut solver, &gas).unwrap();
        }
        for (&before, &after) in cells.iter().zip(solver.states()) {
            assert_eq!(after.density, before.density);
            assert_eq!(after.momentum, before.momentum);
            assert!((after.energy / before.energy - 1.0).abs() < 1e-12);
            let temperature = gas.primitive(after, || false).unwrap().temperature;
            assert!((temperature - wall_temperature).abs() < 1e-8);
        }
    }
}

#[test]
fn g1_caloric_wall_history_retains_temperature_meaning_across_density_changes() {
    let state = air();
    for slope in [0.0, 0.003] {
        let gas = caloric_source(&state, slope);
        for next_density in [0.425, 1.7, 6.8] {
            let mut walls = caloric_walls(2, 300.0, 0.0, true);
            let initial = gas
                .conserved(1.7, 0.0, 1.7 * gas.gas_constant() * 600.0)
                .unwrap();
            let mut solver = caloric_field(&gas, vec![initial; 2]);
            walls.apply(&mut solver, &gas).unwrap();
            let first_t = gas
                .primitive(solver.states()[0], || false)
                .unwrap()
                .temperature;
            assert!((first_t - 480.0).abs() < 1e-8);
            // An accepted transport step changes density while retaining the
            // cell's temperature. The stationary wall's history is unchanged.
            let transported = gas
                .conserved(
                    next_density,
                    0.0,
                    next_density * gas.gas_constant() * first_t,
                )
                .unwrap();
            let mut solver = caloric_field(&gas, vec![transported; 2]);
            walls.apply(&mut solver, &gas).unwrap();
            // exp(G dt)=[[.6,.4],[.4,.6]]: first history is 120 K,
            // so the second mean is .6*180 + .4*120 = 156 K above the wall.
            let expected_t = 456.0;
            let expected = gas
                .conserved(
                    next_density,
                    0.0,
                    next_density * gas.gas_constant() * expected_t,
                )
                .unwrap();
            let next = solver.states()[0];
            assert_eq!(next.density, next_density);
            assert_eq!(next.momentum, 0.0);
            assert!((gas.primitive(next, || false).unwrap().temperature - expected_t).abs() < 1e-8);
            assert!((next.energy / expected.energy - 1.0).abs() < 1e-10);
        }
    }
}

fn assert_wall_compliance(walls: &Walls, transport: WallGas, radius: f64, dt: f64) {
    let thermal = &walls.propagators[walls.thermal_of[0]];
    let gamma = transport.gamma;
    for frequency in [20.0, 100.0, 1000.0, 4000.0, 8000.0] {
        let omega = core::f64::consts::TAU * frequency;
        let s = C64::new(0.0, omega);
        // The mean response of the ACTUAL bound generator to a
        // heating source is e0^T(sI-G)^-1 e0 = 1/(s+H(s)).
        let response =
            thermal
                .head
                .iter()
                .zip(&thermal.decay)
                .fold(C64::ZERO, |total, (&head, &decay)| {
                    let rate = -det::ln(decay) / dt;
                    total + C64::from_re(head * head) / (s + C64::from_re(rate))
                });
        // T' is driven by (gamma-1) T0 rho'/rho0; the ideal-gas
        // relation then gives p'/(R T0 rho') below.
        let actual = C64::from_re(1.0) + (s * response).scale(gamma - 1.0);
        let shear = det::sqrt(
            omega * transport.density_kg_m3 * radius * radius * transport.prandtl
                / transport.dynamic_viscosity_pa_s,
        );
        let f = fs_phs::zwikker_kosten_f(shear).unwrap();
        let expected = C64::from_re(gamma) / (C64::from_re(1.0) + f.scale(gamma - 1.0));
        let error = (actual - expected).abs() / expected.abs();
        assert!(
            error < 1e-3,
            "wall gas {transport:?}, f {frequency}: {error}"
        );
    }
}

#[test]
fn g1_temperature_wall_preserves_the_exact_zk_compliance() {
    let state = air();
    let radius = 1.0e-3;
    let dt = 1.0e-7;
    for density_factor in [0.25, 1.0, 4.0] {
        for gamma in [state.gamma, 5.0 / 3.0] {
            let transport = WallGas {
                density_kg_m3: state.density * density_factor,
                gamma,
                ..wall_gas(&state)
            };
            let walls = Walls::new(
                &[radius],
                WallMemory {
                    poles: 24,
                    band_hz: [20.0, 8000.0],
                    maximum_relative_error: 1e-2,
                },
                transport,
                293.15,
                dt,
            )
            .unwrap();
            assert_wall_compliance(&walls, transport, radius, dt);
        }
    }
}

#[test]
fn g1_physical_wall_reaches_isothermal_equilibrium_at_each_density() {
    let state = air();
    let gas = air_source(&state);
    let reference = 293.15;
    let densities = [0.25 * state.density, state.density, 4.0 * state.density];
    let mut walls = Walls::new(
        &[1.0e-3; 3],
        WallMemory {
            poles: 24,
            band_hz: [20.0, 8000.0],
            maximum_relative_error: 1e-2,
        },
        wall_gas(&state),
        reference,
        0.01,
    )
    .unwrap();
    let cells = densities
        .iter()
        .map(|&rho| {
            gas.conserved(rho, 0.0, rho * gas.gas_constant() * 500.0)
                .unwrap()
        })
        .collect();
    let mut solver = caloric_field(&gas, cells);
    for _ in 0..200 {
        walls.apply(&mut solver, &gas).unwrap();
    }
    for (&rho, &next) in densities.iter().zip(solver.states()) {
        assert_eq!(next.density, rho);
        assert!((gas.primitive(next, || false).unwrap().temperature - reference).abs() < 1e-7);
        let isothermal_pressure = rho * gas.gas_constant() * reference;
        assert!(
            (gas.primitive(next, || false).unwrap().pressure / isothermal_pressure - 1.0).abs()
                < 1e-9
        );
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

// Independent constant-Cp shock/fan characteristic (Toro's pressure function).
// The production path evaluates the original NASA-9 source instead.
fn constant_cp_end_velocity(interior: Primitive, gamma: f64, pressure: f64) -> f64 {
    let jump = if pressure > interior.pressure {
        let beta = (gamma - 1.0) / (gamma + 1.0);
        (pressure - interior.pressure)
            * det::sqrt(
                2.0 / ((gamma + 1.0) * interior.density) / (pressure + beta * interior.pressure),
            )
    } else {
        2.0 * interior.sound_speed / (gamma - 1.0)
            * ((pressure / interior.pressure).powf((gamma - 1.0) / (2.0 * gamma)) - 1.0)
    };
    interior.velocity - jump
}

#[test]
fn radiation_pressure_solve_matches_shock_and_fan_roots_after_large_trace_changes() {
    let state = air();
    let gas = air_source(&state);
    let room_q = gas.conserved(state.density, 0.0, state.pressure).unwrap();
    let room = gas.primitive(room_q, || false).unwrap();
    // Pressure relative to room, interior Mach number, previous pressure seed,
    // and a steady past end velocity initializing the radiation memory.
    for (pressure_ratio, mach, seed_ratio, past_mach) in [
        (0.1, -0.5, 1.0, 0.0),
        (0.1, -0.5, 100.0, 0.0),
        (5.0, 0.5, 0.1, 0.0),
        (1.0, -0.1, 1.0, 0.05),
    ] {
        let mut radiation = PistonRadiation::new(0.006, room, STEP_S).unwrap();
        radiation.guess = seed_ratio * room.pressure;
        radiation
            .history
            .fill(radiation.area * past_mach * room.sound_speed);
        let pressure = pressure_ratio * room.pressure;
        let q = gas
            .conserved(
                pressure / (gas.gas_constant() * room.temperature),
                mach * room.sound_speed,
                pressure,
            )
            .unwrap();
        let interior = gas.primitive(q, || false).unwrap();
        let stage = PortStage {
            stage: 0,
            time: 0.0,
            dt: STEP_S,
            left: (q, interior),
            right: (q, interior),
        };
        let held = radiation.impedance * radiation.memory(0);
        let gain = radiation.rho_c * (1.0 - radiation.weights[0]);
        if pressure_ratio == 0.1 && seed_ratio == 1.0 {
            // This admitted initial source wave used to abort immediately:
            // its unguarded Newton iterate has negative absolute pressure.
            let seed = gas
                .end_wave(interior, 1.0, radiation.guess, || false)
                .unwrap();
            let residual = radiation.guess - room.pressure - gain * seed.primitive.velocity + held;
            assert!(radiation.guess - residual / (1.0 - gain * seed.velocity_per_pressure) < 0.0);
        } else if seed_ratio != 1.0 {
            // A previous trace's valid pressure can be outside this trace's
            // caloric wave path; it must not become a physical admission rule.
            assert!(
                gas.end_wave(interior, 1.0, radiation.guess, || false)
                    .is_err()
            );
        }
        let residual = |p: f64| {
            p - room.pressure - gain * constant_cp_end_velocity(interior, state.gamma, p) + held
        };
        let (mut lower, mut upper) = (room.pressure * 0.001, room.pressure * 10.0);
        assert!(residual(lower) < 0.0 && residual(upper) > 0.0);
        for _ in 0..96 {
            let midpoint = lower.midpoint(upper);
            if residual(midpoint) < 0.0 {
                lower = midpoint;
            } else {
                upper = midpoint;
            }
        }
        let reference_pressure = lower.midpoint(upper);
        let reference_velocity =
            constant_cp_end_velocity(interior, state.gamma, reference_pressure);
        let history = radiation.history.clone();
        let flux = radiation.flux(&gas, &stage).unwrap();
        let (flow, solved_pressure) = radiation.pending.unwrap();
        assert!((solved_pressure / reference_pressure - 1.0).abs() < 2e-10);
        assert!((flow / radiation.area - reference_velocity).abs() < 2e-8);
        assert!(residual(solved_pressure).abs() < 2e-7 * room.pressure);
        assert_eq!(
            radiation.history, history,
            "trial flow must remain uncommitted"
        );
        // The actual Euler flux retains the correct entropy donor on each
        // side of zero flow, including the pressure contribution to traction.
        let donor = if reference_velocity >= 0.0 {
            interior
        } else {
            room
        };
        let face = gas
            .end_wave(donor, 1.0, reference_pressure, || false)
            .unwrap()
            .primitive;
        let mass_flux = face.density * reference_velocity;
        let momentum_flux = mass_flux * reference_velocity + reference_pressure;
        assert!((flux.mass / mass_flux - 1.0).abs() < 2e-10);
        assert!((flux.momentum / momentum_flux - 1.0).abs() < 2e-10);
    }
}

#[test]
fn radiation_pressure_refuses_a_root_outside_the_original_caloric_path_without_committing() {
    let state = air();
    let gas = air_source(&state);
    let room_q = gas.conserved(state.density, 0.0, state.pressure).unwrap();
    let room = gas.primitive(room_q, || false).unwrap();
    let mut radiation = PistonRadiation::new(0.006, room, STEP_S).unwrap();
    // The outgoing characteristic requires compression above a trace that is
    // already at the source's maximum temperature. No supported root exists.
    let pressure = 0.5 * room.pressure;
    let q = gas
        .conserved(pressure / (gas.gas_constant() * 1000.0), 0.0, pressure)
        .unwrap();
    let interior = gas.primitive(q, || false).unwrap();
    let stage = PortStage {
        stage: 0,
        time: 0.0,
        dt: STEP_S,
        left: (q, interior),
        right: (q, interior),
    };
    let history = radiation.history.clone();
    let guess = radiation.guess;
    assert!(matches!(
        radiation.flux(&gas, &stage),
        Err(GasError::CaloricPath(_))
    ));
    assert_eq!(radiation.history, history);
    assert!(radiation.pending.is_none());
    assert_eq!(radiation.guess, guess);
}

#[test]
fn step_mesh_keeps_exact_volume_and_puts_a_junction_face_on_the_step() {
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
    // A junction face lies exactly on the step, so no cell holds it: the cell on its left sees
    // the first cone's outlet area, the cell on its right the second cone's inlet area, and the
    // cells of each stretch are equal and none is longer than the declared maximum.
    let step = 0.06597;
    let face = mesh.positions().iter().position(|&x| x == step).unwrap();
    let a = |r: f64| core::f64::consts::PI * r * r;
    assert!((mesh.face_areas()[face] / a(4.2e-3) - 1.0).abs() <= 1e-15);
    let steps = bore.solver().radius_steps();
    assert_eq!(steps.len(), 1);
    assert_eq!(steps[0].0, face);
    assert!((steps[0].1 / a(4.64e-3) - 1.0).abs() <= 1e-15);
    let widths: Vec<f64> = mesh.positions().windows(2).map(|w| w[1] - w[0]).collect();
    assert!(widths.iter().all(|&w| w <= 2.4e-3));
    assert!(widths[..face].windows(2).all(|w| (w[0] - w[1]).abs() <= 1e-15));
    assert!(widths[face..].windows(2).all(|w| (w[0] - w[1]).abs() <= 1e-15));
    // The cells beside the step carry exact piecewise volumes: the cone's own frustum.
    let left_cell = frustum(
        1.825e-3 + (4.2e-3 - 1.825e-3) * (mesh.positions()[face - 1] / 0.06597),
        4.2e-3,
        step - mesh.positions()[face - 1],
    );
    assert!((mesh.volumes()[face - 1] / left_cell - 1.0).abs() <= 1e-12);
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
fn a_128_pole_wall_memory_admits_a_100_mm_bell_radius_on_the_full_band() {
    // The old pole cap, 64, leaves the excess wall operators of a 0.1 m radius 28 % (complex) and 36 % (loss) off on
    // 20-8000 Hz (fs-phs's own admission; the excess is a tiny part of the whole operator there); 128 poles bring it under 1 %.
    let state = air();
    let duct = Duct {
        segments: vec![
            Segment::Cylinder {
                radius: 4.0e-3,
                length: 0.3,
            },
            Segment::Cone {
                inlet_radius: 4.0e-3,
                outlet_radius: 105.0e-3,
                length: 0.4,
            },
        ],
    };
    let build = |poles: usize, tolerance: f64| {
        NonlinearBore::new(
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
                    poles,
                    band_hz: [20.0, 8000.0],
                    maximum_relative_error: tolerance,
                }),
            },
        )
    };
    assert!(build(64, 1e-2).is_err(), "64 poles must not be admitted at 1 % on a 105 mm radius");
    let bore = build(128, 1e-2).unwrap();
    let receipt = bore.wall_receipt().unwrap();
    println!("128-pole wall {receipt:?}");
    assert!(receipt.maximum_complex_error < 1e-2 && receipt.maximum_loss_error < 1e-2);
    assert_eq!(receipt.maximum_arms, 128);
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

// ---- the first-order retry of a step that leaves the source domain (MOOD a posteriori fallback) ----

fn retry_bore() -> NonlinearBore {
    let state = air();
    let duct = Duct {
        segments: vec![Segment::Cylinder {
            radius: 5.0e-3,
            length: 0.2,
        }],
    };
    NonlinearBore::new(
        &duct,
        air_source(&state),
        293.15,
        101_325.0,
        None,
        BoreSpec {
            maximum_cell_m: 4.0e-3,
            inlet_volume_m3: 0.0,
            step_s: STEP_S,
            cfl: 0.45,
            wall: None,
        },
    )
    .unwrap()
}

// The step as it was before the retry existed: one second-order attempt, nothing else.
fn unretried_step(
    bore: &mut NonlinearBore,
    reconstruction: Reconstruction,
    mut inlet: impl FnMut(&PortStage) -> Result<Option<FaceFlux>, GasError>,
) -> Result<StepReport, GasError> {
    let radiation = &mut bore.radiation;
    let gas = &bore.gas;
    let controls = StepControl {
        cfl: bore.cfl,
        reconstruction,
        maximum_halvings: 0,
    };
    let report = bore.solver.step_with_ports(
        Time::new(bore.step_s),
        controls,
        |stage| Ok([inlet(stage)?, Some(radiation.flux(gas, stage)?)]),
        || false,
    )?;
    radiation.commit().unwrap();
    Ok(report)
}

fn bits(bore: &NonlinearBore) -> Vec<[u64; 3]> {
    bore.solver()
        .states()
        .iter()
        .map(|q| {
            [
                q.density.to_bits(),
                q.momentum.to_bits(),
                q.energy.to_bits(),
            ]
        })
        .collect()
}

// Where two fields first differ, and the largest relative difference of any component (0.0 when they are bit-identical).
fn difference(a: &[[u64; 3]], b: &[[u64; 3]]) -> (Option<(usize, usize, f64, f64)>, f64) {
    let (mut first, mut worst) = (None, 0.0_f64);
    for (cell, (x, y)) in a.iter().zip(b).enumerate() {
        for (k, (&p, &q)) in x.iter().zip(y).enumerate() {
            let (u, v) = (f64::from_bits(p), f64::from_bits(q));
            if p != q {
                first.get_or_insert((cell, k, u, v));
                worst = worst.max((u - v).abs() / u.abs().max(v.abs()));
            }
        }
    }
    (first, worst)
}

// A two-state field: the left half and the right half each at (temperature K, pressure over the room pressure, Mach number of its own sound speed).
fn planted_two_state(left: (f64, f64, f64), right: (f64, f64, f64)) -> Option<NonlinearBore> {
    let mut bore = retry_bore();
    let room_p = bore.room().pressure;
    let gas = bore.gas().clone();
    let make = |(t, ratio, mach): (f64, f64, f64)| {
        let p = room_p * ratio;
        let rho = p / (gas.gas_constant() * t);
        let rest = gas.conserved(rho, 0.0, p).ok()?;
        let c = gas.primitive(rest, || false).ok()?.sound_speed;
        gas.conserved(rho, mach * c, p).ok()
    };
    let (l, r) = (make(left)?, make(right)?);
    let cells = bore.solver().states().len();
    bore.solver
        .apply_cell_sources(|i, q, _| Ok(if i < cells / 2 { l } else { r }), || false)
        .ok()?;
    Some(bore)
}

#[test]
fn a_step_that_never_fails_admission_is_the_unretried_step_bit_for_bit() {
    let (mut with_retry, mut without) = (retry_bore(), retry_bore());
    let source = with_retry.gas().clone();
    let inlet_area = with_retry.inlet_area_m2();
    let width = 1.0e-4;
    let flow = |t: f64| {
        if t < width {
            1.0e-8 * 0.5 * (1.0 - det::cos(core::f64::consts::TAU * t / width))
        } else {
            0.0
        }
    };
    let (mut guess_a, mut guess_b) = (with_retry.room().pressure, without.room().pressure);
    for _ in 0..800 {
        with_retry
            .step(|stage| {
                let (flux, _) =
                    velocity_port(&source, stage, flow(stage.time) / inlet_area, &mut guess_a)?;
                Ok(Some(flux))
            })
            .unwrap();
        unretried_step(&mut without, Reconstruction::Muscl, |stage| {
            let (flux, _) =
                velocity_port(&source, stage, flow(stage.time) / inlet_area, &mut guess_b)?;
            Ok(Some(flux))
        })
        .unwrap();
    }
    assert_eq!(with_retry.first_order_fallbacks(), 0);
    assert_eq!(bits(&with_retry), bits(&without));
    assert_eq!(
        with_retry.solver().time().value().to_bits(),
        without.solver().time().value().to_bits()
    );
}

#[test]
fn a_step_that_fails_at_second_order_is_redone_at_first_order_and_counted() {
    // The left half 101 K gas at 100 x the room pressure, the right half the room: step 0 is admitted at second order, step 1 is not.
    let planted = || planted_two_state((101.0, 100.0, 0.0), (293.15, 1.0, 0.0)).unwrap();
    // `refused` takes the same two attempts by hand; `fresh` takes only the first-order attempt, with no refused attempt before it.
    let (mut retried, mut refused, mut fresh) = (planted(), planted(), planted());
    retried.step(|_| Ok(None)).unwrap();
    unretried_step(&mut refused, Reconstruction::Muscl, |_| Ok(None)).unwrap();
    unretried_step(&mut fresh, Reconstruction::Muscl, |_| Ok(None)).unwrap();
    assert_eq!(retried.first_order_fallbacks(), 0);
    assert_eq!(bits(&retried), bits(&refused));
    assert_eq!(bits(&retried), bits(&fresh));
    // The second-order step 1 leaves the source domain and leaves the field and the clock as they were ...
    let (field, clock) = (bits(&refused), refused.solver().time().value().to_bits());
    let guess = refused.radiation.guess;
    let refusal = unretried_step(&mut refused, Reconstruction::Muscl, |_| Ok(None)).unwrap_err();
    assert!(source_domain_refusal(&refusal), "{refusal:?}");
    assert_eq!(bits(&refused), field);
    assert_eq!(refused.solver().time().value().to_bits(), clock);
    // ... and the first-order step from the same field passes the same admission. The retried bore takes the same two attempts as
    // `refused` (the end's uncommitted flow and starting pressure put back between them), so its field and clock are bit-identical.
    refused.radiation.pending = None;
    refused.radiation.guess = guess;
    unretried_step(&mut refused, Reconstruction::FirstOrder, |_| Ok(None)).unwrap();
    unretried_step(&mut fresh, Reconstruction::FirstOrder, |_| Ok(None)).unwrap();
    let step = retried.step(|_| Ok(None)).unwrap();
    assert_eq!(retried.first_order_fallbacks(), 1);
    assert!(step.report.first_order_retry);
    let (first, worst) = difference(&bits(&retried), &bits(&refused));
    assert!(
        first.is_none(),
        "retried vs refused-then-first-order: first difference (cell, component, retried, refused) {first:?}, worst relative {worst:e}"
    );
    assert_eq!(
        retried.solver().time().value().to_bits(),
        refused.solver().time().value().to_bits()
    );
    // Against a first-order step that had no refused attempt before it: the shared EOS temperature solve starts from the
    // primitives the refused attempt left behind and stops inside its own energy tolerance, so only that spread is expected.
    let (first, worst) = difference(&bits(&retried), &bits(&fresh));
    eprintln!(
        "MOOD-RETRY-SPREAD vs a first-order step with no refused attempt: first difference (cell, component, retried, fresh) {first:?}, worst relative {worst:e}"
    );
    assert_eq!(
        retried.solver().time().value().to_bits(),
        fresh.solver().time().value().to_bits()
    );
}

#[test]
fn a_step_that_fails_at_both_orders_is_refused_naming_the_retry_and_changes_nothing() {
    // 101 K gas flowing apart at twice its sound speed leaves the source domain at the first RK stage, at either order.
    let mut bore = planted_two_state((101.0, 1.0, -2.0), (101.0, 1.0, 2.0)).unwrap();
    let (field, clock) = (bits(&bore), bore.solver().time().value().to_bits());
    let error = bore.step(|_| Ok(None)).unwrap_err();
    let text = error.to_string();
    assert!(
        matches!(error, NonlinearDuctError::GasAfterFirstOrderRetry(_)),
        "{text}"
    );
    assert!(text.starts_with("FS-DUCT-NONLINEAR-GAS: "), "{text}");
    assert!(text.contains("outside source domain"), "{text}");
    assert!(text.ends_with(" after first-order retry"), "{text}");
    assert_eq!(bore.first_order_fallbacks(), 0);
    assert_eq!(bits(&bore), field);
    assert_eq!(bore.solver().time().value().to_bits(), clock);
}

#[test]
fn a_refusal_that_is_not_a_source_domain_failure_is_not_retried() {
    // At Mach 4 the radiating end refuses the fixed step (the CFL bound is shorter than it): not a domain failure, so no retry and the plain text.
    let mut bore = planted_two_state((101.0, 1.0, -4.0), (101.0, 1.0, 4.0)).unwrap();
    let field = bits(&bore);
    let error = bore.step(|_| Ok(None)).unwrap_err();
    let text = error.to_string();
    assert!(matches!(error, NonlinearDuctError::Gas(_)), "{text}");
    assert!(
        text.contains("radiating end runs on its declared fixed step"),
        "{text}"
    );
    assert!(!text.contains("after first-order retry"), "{text}");
    assert_eq!(bore.first_order_fallbacks(), 0);
    assert_eq!(bits(&bore), field);
}

// ---- the radius-step junction face (MagentaBadger 7918, option A) ----

fn half_cosine(t: f64, width: f64, peak: f64) -> f64 {
    if t < width {
        peak * 0.5 * (1.0 - det::cos(core::f64::consts::TAU * t / width))
    } else {
        0.0
    }
}

fn junction_bore(duct: &Duct) -> NonlinearBore {
    let state = air();
    NonlinearBore::new(
        duct,
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
    .unwrap()
}

// Drive the inlet face with a prescribed volume flow [m^3/s] for `steps` fixed steps; the
// inlet-face pressure above the room pressure after each step.
fn drive_and_record(bore: &mut NonlinearBore, steps: usize, flow: impl Fn(f64) -> f64) -> Vec<f64> {
    let source = bore.gas().clone();
    let inlet_area = bore.inlet_area_m2();
    let room = bore.room().pressure;
    let mut guess = room;
    let mut record = Vec::with_capacity(steps);
    for _ in 0..steps {
        let mut face_pressure = room;
        bore.step(|stage| {
            let (flux, pressure) =
                velocity_port(&source, stage, flow(stage.time) / inlet_area, &mut guess)?;
            face_pressure = pressure;
            Ok(Some(flux))
        })
        .unwrap();
        record.push(face_pressure - room);
    }
    record
}

#[test]
fn a_radius_step_reflects_a_pulse_with_the_transfer_matrix_coefficient() {
    // Bound fixed before the first run (brasstap/predictions-step-junction.txt, T1): a 6 -> 7 mm step
    // at x = 0.5 m in a 2.0 m walls-off bore, 2.4 mm cells, a half-cosine inlet flow of 1e-6 m^3/s and
    // 1e-4 s. The time integral of the pressure the step sends back, over that of the pressure sent,
    // is (A1 - A2) / (A1 + A2) within 2e-3 absolute. The returning pulse meets the driven face with
    // its velocity held at zero after the pulse (a rigid end), where the pressure doubles.
    let duct = Duct {
        segments: vec![
            Segment::Cylinder {
                radius: 6.0e-3,
                length: 0.5,
            },
            Segment::Cylinder {
                radius: 7.0e-3,
                length: 1.5,
            },
        ],
    };
    let mut bore = junction_bore(&duct);
    let pressure = drive_and_record(&mut bore, 4300, |t| half_cosine(t, 1.0e-4, 1.0e-6));
    let integral = |from: f64, to: f64| {
        let first = (from / STEP_S).round() as usize;
        let last = (to / STEP_S).round() as usize;
        pressure[first..last].iter().sum::<f64>() * STEP_S
    };
    let sent = integral(0.0, 5.0e-4);
    let returned = 0.5 * integral(2.5e-3, 3.5e-3);
    let (a1, a2) = (6.0e-3_f64.powi(2), 7.0e-3_f64.powi(2));
    let expected = (a1 - a2) / (a1 + a2);
    let ratio = returned / sent;
    println!("T1: returned / sent = {ratio:.6}, (A1 - A2)/(A1 + A2) = {expected:.6}");
    assert!((ratio - expected).abs() <= 2.0e-3, "{ratio} against {expected}");
}

// The states and every recorded pressure of a driven bore, as one hash.
fn digest(bore: &NonlinearBore, pressure: &[f64]) -> String {
    let mut bytes = Vec::new();
    for cell in bits(bore) {
        for word in cell {
            bytes.extend_from_slice(&word.to_le_bytes());
        }
    }
    bytes.extend_from_slice(&bore.solver().time().value().to_bits().to_le_bytes());
    for p in pressure {
        bytes.extend_from_slice(&p.to_bits().to_le_bytes());
    }
    fs_blake3::hash_bytes(&bytes).to_hex()
}

fn driven_digest(duct: &Duct) -> String {
    let mut bore = junction_bore(duct);
    let pressure = drive_and_record(&mut bore, 2000, |t| half_cosine(t, 2.0e-4, 2.0e-5));
    digest(&bore, &pressure)
}

#[test]
fn a_bore_without_a_radius_step_is_the_scheme_before_the_junction_face_bit_for_bit() {
    // Stamped before the first run (T5): the digests were recorded from the scheme of d1de01bf6, which
    // has no junction face, by this same test; the new scheme must reproduce them bit for bit.
    let cylinder = Duct {
        segments: vec![Segment::Cylinder {
            radius: 5.0e-3,
            length: 0.5,
        }],
    };
    let taper = Duct {
        segments: vec![
            Segment::Cone {
                inlet_radius: 4.0e-3,
                outlet_radius: 9.0e-3,
                length: 0.4,
            },
            Segment::Cylinder {
                radius: 9.0e-3,
                length: 0.1,
            },
        ],
    };
    let (got_cylinder, got_taper) = (driven_digest(&cylinder), driven_digest(&taper));
    println!("T5 digests: cylinder {got_cylinder} taper {got_taper}");
    assert_eq!(
        got_cylinder,
        "543286989ddfaabe01f74d7ee6140136a21c1cecb1fdc590b4598509bd4e655f"
    );
    assert_eq!(
        got_taper,
        "dec58d42a8060a698ca4d645e22404667e385f66294b07172a4bf6802cd4f09b"
    );
}

// ---- the plane-mode equivalent of the flanged multimodal mouth as the radiating end (MB 14831 item 1; stamp
// brasstap/predictions-modal-mouth-kernel-v2.txt, bounds fixed before the first run) ----

fn mouth_duct() -> Duct {
    Duct {
        segments: vec![Segment::Cylinder {
            radius: 5.0e-3,
            length: 0.5,
        }],
    }
}

fn mouth_bore(mouth: Option<&MouthKernel>) -> Result<NonlinearBore, NonlinearDuctError> {
    let state = air();
    NonlinearBore::with_mouth_kernel(
        &mouth_duct(),
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
        mouth,
    )
}

// The piston's own weights wrapped as a kernel for the same mouth: the plumbing control.
fn piston_as_mouth_kernel(piston: &NonlinearBore) -> MouthKernel {
    let room = piston.room();
    let radius = 5.0e-3;
    MouthKernel {
        step_s: STEP_S,
        z_inf: room.density * room.sound_speed / (core::f64::consts::PI * radius * radius),
        mouth_radius_m: radius,
        weights: piston.radiation_weights().to_vec(),
        dropped_tail_energy: 0.0,
        dropped_sum: 0.0,
        bin_hz: 0.0,
        evaluated_hz: 0.0,
    }
}

fn worst_relative_difference(a: &[f64], b: &[f64]) -> (f64, f64) {
    let peak = a.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
    let worst = a.iter().zip(b).fold(0.0_f64, |m, (x, y)| m.max((x - y).abs()));
    (worst, peak)
}

#[test]
fn a_mouth_kernel_of_the_pistons_own_weights_is_the_piston() {
    // K1 (amended before any run: the last ulp of z_inf x area is not rho c, so the bound is 1e-12 of the peak, not the bits).
    let mut piston = mouth_bore(None).unwrap();
    let kernel = piston_as_mouth_kernel(&piston);
    let mut through_kernel = mouth_bore(Some(&kernel)).unwrap();
    assert_eq!(through_kernel.radiation_weights(), piston.radiation_weights());
    let drive = |t: f64| half_cosine(t, 2.0e-4, 2.0e-5);
    let a = drive_and_record(&mut piston, 2000, &drive);
    let b = drive_and_record(&mut through_kernel, 2000, &drive);
    let (worst, peak) = worst_relative_difference(&a, &b);
    println!("K1: worst |difference| {worst:.3e} Pa against a peak of {peak:.3e} Pa");
    assert!(peak > 0.0);
    assert!(worst <= 1.0e-12 * peak, "{worst} against {peak}");
}

#[test]
fn a_mouth_kernel_that_does_not_fit_the_bore_is_refused() {
    let piston = mouth_bore(None).unwrap();
    let kernel = piston_as_mouth_kernel(&piston);
    let mut wrong_step = kernel.clone();
    wrong_step.step_s *= 1.001;
    let mut wrong_gas = kernel.clone();
    wrong_gas.z_inf *= 1.01;
    let mut short = kernel.clone();
    short.weights.truncate(1);
    let mut not_finite = kernel.clone();
    not_finite.weights[1] = f64::NAN;
    for (name, bad) in [
        ("step", wrong_step),
        ("gas", wrong_gas),
        ("one tap", short),
        ("a NaN weight", not_finite),
    ] {
        assert!(mouth_bore(Some(&bad)).is_err(), "a kernel with the wrong {name} was admitted");
    }
    assert!(mouth_bore(Some(&kernel)).is_ok());
}

#[test]
fn a_kernel_naming_a_distinct_exterior_is_admitted_by_the_bore_that_holds_its_bore_gas() {
    // B2 (MB 24783): the bore end is matched to the gas the bore holds, so a kernel built with that gas as its bore gas and ANOTHER gas as the exterior fits it; the kernel
    // built the other way round (the room's gas as bore gas) names another characteristic impedance and is refused by that bore, and fits a bore that holds the room's gas.
    use crate::modal::{MouthKernelSpec, plane_mouth_kernel};
    let room = air();
    let bore_gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 301.0, 101_325.0).unwrap();
    let spec = MouthKernelSpec { step_s: STEP_S, bin_hz: 20.0, band_hz: 6000.0, tail_tolerance: 1.0e-10 };
    let build = |state: &GasState, kernel: &MouthKernel| {
        NonlinearBore::with_mouth_kernel(
            &mouth_duct(),
            air_source(state),
            state.temperature,
            state.pressure,
            None,
            BoreSpec { maximum_cell_m: 2.4e-3, inlet_volume_m3: 0.0, step_s: STEP_S, cfl: 0.45, wall: None },
            Some(kernel),
        )
    };
    let right = plane_mouth_kernel(&mouth_duct(), &bore_gas, &room, LossModel::Lossless, 1, 1, spec).unwrap();
    let wrong = plane_mouth_kernel(&mouth_duct(), &room, &bore_gas, LossModel::Lossless, 1, 1, spec).unwrap();
    assert!(build(&bore_gas, &wrong).is_err(), "a kernel referenced to the room's gas was admitted by a bore of another gas");
    assert!(build(&room, &wrong).is_ok(), "a kernel referenced to the room's gas must fit a bore that holds it");
    let mut bore = build(&bore_gas, &right).expect("a kernel referenced to the bore's own gas must fit");
    let drive = |t: f64| half_cosine(t, 2.0e-4, 2.0e-5);
    let recorded = drive_and_record(&mut bore, 2000, &drive);
    assert!(recorded.iter().all(|p| p.is_finite()) && recorded.iter().any(|p| *p != 0.0));
}

#[test]
fn the_one_mode_equivalent_kernel_drives_the_bore_like_the_piston() {
    // K4 (bore level, bound fixed before the first run): the one-mode plane-mode equivalent of the flanged mouth IS the baffled
    // piston's impedance (the equivalent module pins that to 1e-10); its realized time-domain kernel (bins 20 Hz, claim 6 kHz,
    // tail 1e-10, at the bore's own 1.953 us step) drives the bore's inlet pressure to within 5e-3 of the piston's peak.
    use crate::modal::{MouthKernelSpec, plane_mouth_kernel};
    let state = air();
    let spec = MouthKernelSpec {
        step_s: STEP_S,
        bin_hz: 20.0,
        band_hz: 6000.0,
        tail_tolerance: 1.0e-10,
    };
    let kernel = plane_mouth_kernel(&mouth_duct(), &state, &state, LossModel::Lossless, 1, 1, spec).unwrap();
    println!(
        "K4: {} taps, mouth radius {:.6e} m, dropped tail energy {:.2e}, dropped sum {:.2e}, evaluated to {:.0} Hz",
        kernel.weights.len(),
        kernel.mouth_radius_m,
        kernel.dropped_tail_energy,
        kernel.dropped_sum,
        kernel.evaluated_hz
    );
    assert!(kernel.weights.len() >= 2);
    let mut piston = mouth_bore(None).unwrap();
    let mut realized = mouth_bore(Some(&kernel)).unwrap();
    let drive = |t: f64| half_cosine(t, 2.0e-4, 2.0e-5);
    let a = drive_and_record(&mut piston, 2000, &drive);
    let b = drive_and_record(&mut realized, 2000, &drive);
    let (worst, peak) = worst_relative_difference(&a, &b);
    println!("K4: worst |difference| {worst:.3e} Pa against a peak of {peak:.3e} Pa ({:.3e})", worst / peak);
    assert!(peak > 0.0);
    assert!(worst <= 5.0e-3 * peak, "{worst} against {peak}");
}

// The breathing-wall shell (RubyPrairie, MB 23820): its tests share this module's fixtures.
#[path = "bore_shell_tests.rs"]
mod shell_tests;
