//! Finite wall power changes the continuing aperture through its original material.
//! Synthetic source laws exercise conservation and retained history, not measurements.
use super::*;
use fs_couple::acoustic_realize::AcousticRealizeError;
use fs_couple::bernoulli_aperture::dg::thermal_wall::FiniteWallThermalSolve;
use fs_couple::bernoulli_aperture::performance::thermal_wall::{
    FiniteWallHeatProgram, WallHeatControl,
};
use fs_couple::bernoulli_aperture::wall::material::{
    AxialConstraint, LocalMaterialWall, RegionalMaterialWallDg, WallMaterialLoss,
};
use fs_couple::cavity_phs::piston::thermal::HeatExchangeSolve;
use fs_couple::dg_port::DgPortWallSource;
use fs_duct::dg::PortInput;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::{
    EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState, HeatCapacitySourceQuery,
    ResolvedHeatCapacityCurve, SolidLiquidPhase, UniformEnthalpyStepInput,
};
use fs_material::state_point::{
    MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_phs::linear_exact::DenseExponentialBudget;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, DynViscosity, Pressure, QuantitySpec, Temperature};

const DENSITY: f64 = 4000.0;
const CP: f64 = 500.0;
const THICKNESS: f64 = 2e-5;
const MODULUS: f64 = 1e8;
const SLOPE: f64 = -0.001;

fn point(t: f64, lot: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            t,
        )
        .unwrap()
        .with("lot", lot)
        .unwrap()
}

fn material(density_slope: f64, elastic_upper: f64) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let mut claims = ClaimSet::new();
    for (name, dims, values, upper) in [
        (
            "specific_heat_capacity",
            SPECIFIC_HEAT_CAPACITY_DIMS,
            [CP; 3],
            800.0,
        ),
        (
            "density",
            Density::DIMS,
            [
                DENSITY - 100.0 * density_slope,
                DENSITY,
                DENSITY + 500.0 * density_slope,
            ],
            800.0,
        ),
        (
            "young_modulus",
            Pressure::DIMS,
            [1.1 * MODULUS, MODULUS, 0.5 * MODULUS],
            elastic_upper,
        ),
        ("poisson_ratio", Dims::NONE, [0.3; 3], 800.0),
        (
            "young_viscosity",
            DynViscosity::DIMS,
            [70.0, 100.0, 250.0],
            800.0,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: [200.0, 300.0, 800.0].into_iter().zip(values).collect(),
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        QuantitySpec::semantic(SemanticType::new(
                            QuantityKind::AbsoluteTemperature,
                            ValueForm::Static,
                        )),
                        200.0,
                        upper,
                    )
                    .with("lot", 7.0, 8.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic finite wall performance fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    // Heating enters an overlapping claim. The original selected source must
    // survive; it must not be replaced by this three-times-stiffer alternative.
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("young_modulus", Pressure::DIMS),
            value: PropertyValue::Scalar {
                value: 3.0 * MODULUS,
                dims: Pressure::DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "T",
                    QuantitySpec::semantic(SemanticType::new(
                        QuantityKind::AbsoluteTemperature,
                        ValueForm::Static,
                    )),
                    310.0,
                    800.0,
                )
                .with("lot", 7.0, 8.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "competing hot-wall source, never selected".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "one synthetic finite wall inventory".into(),
            phase: "solid".into(),
            process: "constant heat capacity and temperature-dependent elasticity".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &point(200.0, 7.0),
            upper: &point(800.0, 7.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-5,
    )
    .unwrap();
    (card, curve)
}

fn wall(card: &MaterialCard, lot: f64, thickness: f64, kelvin: bool) -> LocalMaterialWall {
    let mut properties = vec![
        ("density", Density::DIMS),
        ("young_modulus", Pressure::DIMS),
        ("poisson_ratio", Dims::NONE),
    ];
    if kelvin {
        properties.push(("young_viscosity", DynViscosity::DIMS));
    }
    let requirements = properties
        .into_iter()
        .map(|(name, dims)| {
            ScalarPropertyRequirement::try_with_quantity(
                name,
                QuantitySpec::dimensional(dims),
                ScalarAdmissibility::Finite,
            )
            .unwrap()
        })
        .collect::<Vec<_>>();
    let source = resolve_material_state_point(
        card,
        &point(300.0, lot),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    LocalMaterialWall::new(
        source,
        thickness,
        AxialConstraint::StressFree,
        if kelvin {
            WallMaterialLoss::KelvinVoigt
        } else {
            WallMaterialLoss::Elastic
        },
        [20.0, 100_000.0],
    )
    .unwrap()
}

fn wall_load(card: &MaterialCard, regional: bool, kelvin: bool) -> DgPortLoad {
    let gas = field(300.0, 1.0, 101_325.0, LENGTH);
    let duct = DgDuct::new_with_gas_field(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: RADIUS,
                length: LENGTH,
            }],
        },
        &gas,
        DgMesh {
            elements: if regional { 2 } else { 1 },
            order: 2,
        },
        DgBoundary::Port(PortId(0)),
        DgBoundary::Port(PortId(1)),
    )
    .unwrap();
    // Use the same causal rational terminal as the finite-gas suite.
    let terminal = super::super::load(&gas, false).terminal().unwrap().clone();
    if regional {
        DgPortLoad::new_with_regional_material_wall_and_terminal(
            RegionalMaterialWallDg::new(
                vec![
                    wall(card, 7.0, THICKNESS, kelvin),
                    wall(card, 8.0, 2.0 * THICKNESS, kelvin),
                ],
                duct,
            )
            .unwrap(),
            DT,
            1e-11,
            MAX_DG_PROPAGATOR_DOFS,
            terminal,
        )
        .unwrap()
    } else {
        DgPortLoad::new_with_material_wall_and_terminal(
            wall(card, 7.0, THICKNESS, kelvin).bind_dg(duct).unwrap(),
            DT,
            1e-11,
            MAX_DG_PROPAGATOR_DOFS,
            terminal,
        )
        .unwrap()
    }
}

fn wall_source(card: &MaterialCard, regional: bool, kelvin: bool) -> ApertureDg {
    let load = wall_load(card, regional, kelvin);
    ApertureDg::new(aperture(&load, load.impedance_pa_s_m3()), load).unwrap()
}

fn heat_budget() -> DenseExponentialBudget {
    DenseExponentialBudget {
        max_dimension: 256,
        max_dense_bytes: 64 << 20,
        max_multiply_terms: 40_000_000_000,
    }
}

fn solve(heat: bool) -> FiniteWallThermalSolve {
    FiniteWallThermalSolve {
        exchange: HeatExchangeSolve {
            heat_bracket_j: [-2e-4, 2e-4],
            source_tolerance_j: 1e-13,
            energy_tolerance_j: 1e-11,
            maximum_iterations: 80,
        },
        maximum_absolute_parameter_work_j: 1e-4,
        dissipation_heat: heat.then(heat_budget),
    }
}

fn body(
    curve: &EquilibriumEnthalpyPhaseCurve,
    state: EquilibriumPhaseState,
    regional: bool,
    heat: f64,
) -> UniformEnthalpyStepInput<'_> {
    let area = core::f64::consts::TAU * RADIUS * LENGTH;
    let volume = area * THICKNESS * if regional { 1.5 } else { 1.0 };
    UniformEnthalpyStepInput {
        curve,
        initial: state,
        mass_kg: DENSITY * volume,
        volume_m3: volume,
        surface_area_m2: area,
        internal_heat_j: heat,
        duration_s: DT,
    }
}

fn initial(curve: &EquilibriumEnthalpyPhaseCurve) -> EquilibriumPhaseState {
    curve.state_at_specific_enthalpy(100.0 * CP).unwrap()
}

fn equal_load(a: &DgPortLoad, b: &DgPortLoad) {
    assert_eq!(a.state(), b.state());
    assert_eq!(a.duct().operator_identity(), b.duct().operator_identity());
    assert_eq!(a.duct().gas_field(), b.duct().gas_field());
    assert_eq!(a.material_wall_source(), b.material_wall_source());
    assert_eq!(a.propagator().dt(), b.propagator().dt());
    assert_eq!(a.propagator().tolerance(), b.propagator().tolerance());
    assert_eq!(a.impedance_pa_s_m3(), b.impedance_pa_s_m3());
    assert_eq!(
        a.terminal().unwrap().pole_states(),
        b.terminal().unwrap().pole_states()
    );
    assert_eq!(a.last_terminal_frames(), b.last_terminal_frames());
    assert_eq!(a.stored_energy_j(), b.stored_energy_j());
}

fn equal_source(a: &ApertureDg, b: &ApertureDg) {
    equal_load(a.load(), b.load());
    assert_eq!(a.aperture().state(), b.aperture().state());
    assert_eq!(a.aperture().spec(), b.aperture().spec());
    assert_eq!(a.aperture().accepted_steps(), b.aperture().accepted_steps());
}

fn sources(load: &DgPortLoad) -> Vec<&LocalMaterialWall> {
    match load.material_wall_source().unwrap() {
        DgPortWallSource::Uniform(wall) => vec![wall],
        DgPortWallSource::Regional(walls) => walls.iter().collect(),
    }
}

// Independent Gauss integration of degree-two displacement: k=E*t/r^2,
// surface measure 2*pi*r*dx. Each element retains its own physical thickness.
fn spring_energy(displacement: &[f64], regional: bool) -> f64 {
    let elements = if regional { 2 } else { 1 };
    assert_eq!(displacement.len(), elements * 3);
    let a = (3.0_f64 / 5.0).sqrt();
    displacement
        .chunks_exact(3)
        .enumerate()
        .map(|(element, u)| {
            let integral = [(-a, 5.0 / 9.0), (0.0, 8.0 / 9.0), (a, 5.0 / 9.0)]
                .into_iter()
                .map(|(x, weight)| {
                    let value = 0.5 * x * (x - 1.0) * u[0]
                        + (1.0 - x * x) * u[1]
                        + 0.5 * x * (x + 1.0) * u[2];
                    weight * value * value
                })
                .sum::<f64>()
                * LENGTH
                / (2.0 * elements as f64);
            let thickness = THICKNESS * if regional && element == 1 { 2.0 } else { 1.0 };
            0.5 * MODULUS * thickness / RADIUS.powi(2) * core::f64::consts::TAU * RADIUS * integral
        })
        .sum()
}

#[test]
fn warm_uniform_and_regional_loads_retain_original_sources_motion_and_terminal_memory() {
    let (card, _) = material(0.0, 800.0);
    for regional in [false, true] {
        let mut source = wall_source(&card, regional, false);
        for _ in 0..24 {
            source
                .step(TubeDrive {
                    upstream_pressure_pa: 300.0,
                    body_flow_m3_s: 0.0,
                })
                .unwrap();
        }
        let mut load = source.load().clone();
        let before = load.clone();
        let motion = before.duct().compliant_wall_state(before.state()).unwrap();
        let spring = spring_energy(&motion.displacement_m, regional);
        assert!(spring > 0.0);
        let update = load
            .rebind_wall_temperature_controlled(&card, "T", 300.0, 320.0, 1.0, &mut || false)
            .unwrap();
        let expected = spring * SLOPE * 20.0;
        assert!((update.parameter_work_j - expected).abs() < 2e-8 * spring + 1e-25);
        assert!(
            (load.stored_energy_j() - before.stored_energy_j() - expected).abs()
                < 2e-8 * before.stored_energy_j() + 1e-25
        );
        assert_eq!(
            load.terminal().unwrap().pole_states(),
            before.terminal().unwrap().pole_states()
        );
        assert_eq!(load.last_terminal_frames(), before.last_terminal_frames());
        assert_eq!(load.state().step_count(), before.state().step_count());
        let after = load.duct().compliant_wall_state(load.state()).unwrap();
        for (a, b) in after.displacement_m.iter().zip(&motion.displacement_m) {
            assert!((a - b).abs() < 1e-18);
        }
        for (a, b) in after.velocity_m_s.iter().zip(&motion.velocity_m_s) {
            assert!((a - b).abs() < 1e-12);
        }
        for (a, b) in sources(&load).iter().zip(sources(&before)) {
            assert_eq!(a.material().card_identity(), card.content_hash());
            assert_eq!(a.material().query().axes()["T"], 320.0);
            assert_eq!(
                a.material().query().axes()["lot"],
                b.material().query().axes()["lot"]
            );
            for property in ["density", "young_modulus", "poisson_ratio"] {
                assert_eq!(
                    a.material()
                        .property(property)
                        .unwrap()
                        .answer()
                        .receipt
                        .selected,
                    b.material()
                        .property(property)
                        .unwrap()
                        .answer()
                        .receipt
                        .selected
                );
            }
            assert_eq!(
                a.material()
                    .property("young_modulus")
                    .unwrap()
                    .answer()
                    .receipt
                    .in_domain
                    .len(),
                2
            );
        }
        // An unrelated gas intervention must retain the selected wall source.
        let wall = load.material_wall_source().cloned();
        load.rebind_gas_field_controlled(
            &field(310.0, 1.0, 101_325.0, LENGTH),
            DgGasRefresh::Lossless,
            1.0,
            &mut || false,
        )
        .unwrap();
        assert_eq!(load.material_wall_source(), wall.as_ref());
        assert!(load.port_response().is_ok());
    }
    for (density_slope, upper) in [(1.0, 800.0), (0.0, 315.0)] {
        let (card, _) = material(density_slope, upper);
        let mut source = wall_source(&card, true, false);
        for _ in 0..8 {
            source
                .step(TubeDrive {
                    upstream_pressure_pa: 300.0,
                    body_flow_m3_s: 0.0,
                })
                .unwrap();
        }
        let mut load = source.load().clone();
        let before = load.clone();
        assert!(
            load.rebind_wall_temperature_controlled(&card, "T", 300.0, 320.0, 1.0, &mut || false)
                .is_err()
        );
        equal_load(&load, &before);
    }
}

fn program() -> FiniteWallHeatProgram {
    FiniteWallHeatProgram {
        controls: [(0, 0.0), (7, 600.0), (19, -300.0)]
            .into_iter()
            .map(|(sample, power_w)| WallHeatControl { sample, power_w })
            .collect(),
        coupling: solve(true),
    }
}

fn player(card: &MaterialCard, block: usize) -> AperturePerformance {
    let phrase = GestureSchedule::try_new(
        1,
        vec![GestureTrack {
            id: "mouth".into(),
            target: GestureTarget::BlowingPressure,
            initial: GestureValue::PressurePa(30.0),
            events: vec![],
        }],
    )
    .unwrap();
    AperturePerformance::new(
        CoupledAperture::Dg(wall_source(card, true, true)),
        ApertureObservation::Inlet,
        phrase,
        AperturePerformanceConfig {
            samples: 24,
            ..config(block)
        },
    )
    .unwrap()
}

fn attached(
    card: &Arc<MaterialCard>,
    curve: &EquilibriumEnthalpyPhaseCurve,
    block: usize,
    program: FiniteWallHeatProgram,
    gate: Arc<CancelGate>,
) -> AperturePerformance {
    player(card, block)
        .with_finite_wall_heat(
            body(curve, initial(curve), true, 0.0),
            card.clone(),
            "T".into(),
            program,
            gate,
        )
        .unwrap()
}

#[test]
fn scheduled_wall_heat_closes_independent_body_work_and_uses_the_solved_loaded_outlet() {
    let (card, caloric) = material(0.0, 800.0);
    let card = Arc::new(card);
    let curve = caloric.curve();
    let mut direct = wall_source(&card, true, true);
    let mut state = initial(curve);
    let program = program();
    let gate = Arc::new(CancelGate::new());
    let mut rendered = attached(&card, curve, 1, program.clone(), gate.clone());
    assert!(Arc::ptr_eq(
        rendered.finite_wall_cancel_gate().unwrap(),
        &gate
    ));
    assert_eq!(
        rendered.finite_wall_material().unwrap().content_hash(),
        card.content_hash()
    );
    assert_eq!(rendered.finite_wall_heat_program(), Some(&program));
    let drive = TubeDrive {
        upstream_pressure_pa: 30.0,
        body_flow_m3_s: 0.0,
    };
    let mut total_heat = 0.0;
    let mut outlet_matters = false;
    let mut heating_seen = false;
    let mut cooling_seen = false;
    for sample in 0..24 {
        let before = direct.clone();
        let old_state = state;
        let power = program
            .controls
            .iter()
            .rev()
            .find(|c| c.sample <= sample)
            .unwrap()
            .power_w;
        let input = body(curve, state, true, power * DT);
        let mut frozen = before.clone();
        let ordinary = frozen.step(drive).unwrap();
        let mut ports = vec![PortInput {
            id: PortId(0),
            effort: ordinary.aperture.bore_flow_m3_s,
        }];
        ports.extend(ordinary.terminals.iter().map(|terminal| PortInput {
            id: terminal.port,
            effort: -terminal.outward_flow_m3_s,
        }));
        let source = before
            .load()
            .duct()
            .prepare_compliant_wall_heat(before.load().propagator(), heat_budget(), &mut || false)
            .unwrap();
        let native_heat = source
            .integrated_heat_j(before.duct_state(), &ports, &mut || false)
            .unwrap();
        let inlet_only = source
            .integrated_heat_j(before.duct_state(), &ports[..1], &mut || false)
            .unwrap();
        if (native_heat - inlet_only).abs() > 1e-8 * native_heat.max(1e-30) {
            outlet_matters = true;
        }
        let motion = frozen
            .load()
            .duct()
            .compliant_wall_state(frozen.duct_state())
            .unwrap();
        let expected = direct
            .step_with_finite_wall_heat_controlled(
                &mut state,
                input,
                &card,
                "T",
                drive,
                0.0,
                solve(true),
                &mut || false,
            )
            .unwrap();
        let mut output = [f64::NAN];
        rendered.block(&mut output).unwrap();
        let actual = rendered.last_finite_wall_step().unwrap();
        assert_eq!(
            output[0].to_bits(),
            expected.coupled.aperture.bore_pressure_pa.to_bits()
        );
        assert_eq!(actual.coupled, expected.coupled);
        assert_eq!(actual.wall, expected.wall);
        assert_eq!(actual.body_state, state);
        assert_eq!(rendered.finite_wall_state(), Some(&state));
        equal_source(retained(&rendered), &direct);
        assert_eq!(
            direct.load().terminal().unwrap().pole_states(),
            frozen.load().terminal().unwrap().pole_states()
        );
        assert_eq!(actual.coupled.aperture, ordinary.aperture);
        assert_eq!(actual.wall_heat_j, native_heat);
        total_heat += native_heat;
        let change_t = state.temperature_k() - old_state.temperature_k();
        heating_seen |= power > 0.0 && change_t > 0.0;
        cooling_seen |= power < 0.0 && change_t < 0.0;
        let body_change = input.mass_kg * CP * change_t;
        let work = spring_energy(&motion.displacement_m, true) * SLOPE * change_t;
        assert!((actual.wall.parameter_work_j - work).abs() < 1e-14);
        assert!((body_change + work - input.internal_heat_j - native_heat).abs() < 1e-11);
        let terminal_work: f64 = actual.coupled.terminals.iter().map(|t| t.port_work_j).sum();
        let balance = direct.stored_energy_j() - before.stored_energy_j()
            + body_change
            + actual.exported_dissipation_j
            + terminal_work
            - input.internal_heat_j
            - actual.coupled.upstream_work_j
            - actual.coupled.body_work_j
            - actual.coupled.aperture.mechanical_work_j;
        assert!(balance.abs() <= solve(true).exchange.energy_tolerance_j);
        assert!(
            (actual.exported_dissipation_j + native_heat
                - actual.coupled.inferred_total_dissipation_j)
                .abs()
                < 1e-15
        );
        assert!(actual.body_energy_residual_j.abs() <= solve(true).exchange.energy_tolerance_j);
        assert!(actual.combined_energy_residual_j.abs() <= solve(true).exchange.energy_tolerance_j);
        assert!(actual.source_residual_j.abs() <= solve(true).exchange.source_tolerance_j);
        assert_eq!(
            direct.aperture().spec().impedance_pa_s_m3,
            direct.load().impedance_pa_s_m3()
        );
        let applied = program
            .controls
            .iter()
            .take_while(|event| event.sample <= sample)
            .count();
        assert_eq!(
            rendered.applied_wall_heat_controls(),
            &program.controls[..applied]
        );
    }
    assert!(total_heat > 0.0);
    assert!(
        outlet_matters,
        "the native solid heat must include the solved loaded outlet"
    );
    assert!(heating_seen && cooling_seen);
}

#[test]
fn wall_power_is_block_invariant_and_changes_the_continuing_pressure() {
    let (card, caloric) = material(0.0, 800.0);
    let card = Arc::new(card);
    let mut reference: Option<([f64; 24], AperturePerformance)> = None;
    for block in [1, 7, 24] {
        let mut renderer = attached(
            &card,
            caloric.curve(),
            block,
            program(),
            Arc::new(CancelGate::new()),
        );
        let mut pressure = [0.0; 24];
        for values in pressure.chunks_mut(block) {
            renderer.block(values).unwrap();
        }
        if let Some((expected, previous)) = &reference {
            assert_eq!(&pressure, expected);
            assert_eq!(renderer.finite_wall_state(), previous.finite_wall_state());
            equal_source(retained(&renderer), retained(previous));
        }
        reference = Some((pressure, renderer));
    }
    let mut unheated = player(&card, 24);
    let mut pressure = [0.0; 24];
    unheated.block(&mut pressure).unwrap();
    assert_ne!(pressure, reference.unwrap().0);
}

#[test]
fn refused_wall_power_and_cancelled_samples_preserve_accepted_physics_and_controls() {
    let (card, caloric) = material(0.0, 800.0);
    let card = Arc::new(card);
    let mut reference = attached(
        &card,
        caloric.curve(),
        7,
        program(),
        Arc::new(CancelGate::new()),
    );
    reference.block(&mut [0.0; 7]).unwrap();
    for cancel in [false, true] {
        let gate = Arc::new(CancelGate::new());
        let mut heat = program();
        if !cancel {
            heat.controls[1].power_w = 1e9;
        }
        let mut attempted = attached(&card, caloric.curve(), 7, heat, gate.clone());
        attempted.block(&mut [0.0; 7]).unwrap();
        if cancel {
            gate.request();
        }
        let mut output = [123.0];
        let error = attempted.block(&mut output).unwrap_err();
        if cancel {
            assert!(matches!(
                error,
                RenderError::Voice(AcousticRealizeError::Cancelled)
            ));
        }
        assert_eq!(output, [123.0]);
        assert_eq!(attempted.samples_rendered(), 7);
        assert_eq!(attempted.finite_wall_state(), reference.finite_wall_state());
        assert_eq!(
            attempted.applied_wall_heat_controls(),
            reference.applied_wall_heat_controls()
        );
        assert_eq!(attempted.applied_controls(), reference.applied_controls());
        assert_eq!(
            attempted.last_finite_wall_step().unwrap().coupled,
            reference.last_finite_wall_step().unwrap().coupled
        );
        equal_source(retained(&attempted), retained(&reference));
        assert!(matches!(
            attempted.block(&mut output),
            Err(RenderError::Poisoned)
        ));
    }
    let before = retained(&reference).clone();
    let initial = *reference.finite_wall_state().unwrap();
    let input = body(caloric.curve(), initial, true, 600.0 * DT);
    let drive = TubeDrive {
        upstream_pressure_pa: 30.0,
        body_flow_m3_s: 0.0,
    };
    let mut clean = before.clone();
    let mut clean_state = initial;
    let mut polls = 0;
    let expected = clean
        .step_with_finite_wall_heat_controlled(
            &mut clean_state,
            input,
            &card,
            "T",
            drive,
            0.0,
            solve(true),
            &mut || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut attempted = before.clone();
    let mut state = initial;
    let mut visited = 0;
    assert!(matches!(
        attempted.step_with_finite_wall_heat_controlled(
            &mut state,
            input,
            &card,
            "T",
            drive,
            0.0,
            solve(true),
            &mut || {
                visited += 1;
                visited == polls
            }
        ),
        Err(AcousticRealizeError::Cancelled)
    ));
    assert_eq!(visited, polls);
    equal_source(&attempted, &before);
    assert_eq!(state, initial);
    let retry = attempted
        .step_with_finite_wall_heat_controlled(
            &mut state,
            input,
            &card,
            "T",
            drive,
            0.0,
            solve(true),
            &mut || false,
        )
        .unwrap();
    assert_eq!(retry.coupled, expected.coupled);
    assert_eq!(retry.wall, expected.wall);
    assert_eq!(state, clean_state);
    equal_source(&attempted, &clean);
}
