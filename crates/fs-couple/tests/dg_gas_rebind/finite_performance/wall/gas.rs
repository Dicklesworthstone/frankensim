//! G1/G3/G4: reciprocal finite gas/wall heat on retained acoustic histories.
//! Synthetic source laws provide independent energy and midpoint oracles.
use super::*;
use fs_couple::bernoulli_aperture::dg::thermal_wall::gas::FiniteGasWallThermalSolve;
use fs_couple::bernoulli_aperture::performance::gas_wall::{
    FiniteGasWallHeatProgram, GasWallConduction, GasWallHeatControl,
};
use fs_couple::cavity_phs::piston::thermal::layers::MaterialConductionLayers;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;

fn joint_controls() -> FiniteGasWallThermalSolve {
    FiniteGasWallThermalSolve {
        gas: controls(),
        caloric: energy_solve(),
        wall: solve(true),
        exchange: HeatExchangeSolve {
            heat_bracket_j: [-0.001, 0.001],
            source_tolerance_j: 1e-13,
            energy_tolerance_j: 1e-11,
            maximum_iterations: 80,
        },
        maximum_material_trials: 20_000,
    }
}

fn parcel(t: f64) -> ClosedGasParcel {
    let cold = field(t, 1.0, 101_325.0, LENGTH);
    let point = QueryPoint::new()
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            t,
        )
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            101_325.0,
        )
        .unwrap();
    // This fixture has a rational terminal and NO finite terminal passage.
    ClosedGasParcel::new(
        cold.model().with_composition(composition(1.0)).unwrap(),
        point,
        core::f64::consts::PI * RADIUS * RADIUS * LENGTH,
    )
    .unwrap()
}

fn joint_source(card: &MaterialCard, t: f64) -> ApertureDg {
    let mut source = wall_source(card, true, true);
    source
        .rebind_gas_field_controlled(
            &field(t, 1.0, 101_325.0, LENGTH),
            DgGasRefresh::Lossless,
            1.0,
            &mut || false,
        )
        .unwrap();
    source
}

fn conduction() -> GasWallConduction {
    GasWallConduction {
        context: point(300.0, 7.0),
        conductivity: PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
        selection: MaterialPropertySelection::SingleClaimOnly,
        liquid_fraction_axis: None,
        area_m2: core::f64::consts::TAU * RADIUS * LENGTH,
        gas_distance_m: 2e-5,
        body_distance_m: THICKNESS,
    }
}

fn layers<'a>(
    card: &'a MaterialCard,
    input: &'a GasWallConduction,
) -> MaterialConductionLayers<'a> {
    MaterialConductionLayers {
        material: card,
        context: &input.context,
        temperature_axis: "T",
        liquid_fraction_axis: input.liquid_fraction_axis.as_deref(),
        conductivity: &input.conductivity,
        selection: input.selection.clone(),
        area_m2: input.area_m2,
        gas_distance_m: input.gas_distance_m,
        body_distance_m: input.body_distance_m,
    }
}

fn input_drive() -> TubeDrive {
    TubeDrive {
        upstream_pressure_pa: 30.0,
        body_flow_m3_s: 0.0,
    }
}

// Integral of the synthetic NASA cv/R = 2.5 + 0.001 T. Neither a production
// caloric inverse nor reported energy residual is used for this oracle.
fn gas_energy_change(gas: &ClosedGasParcel, t0: f64, t1: f64) -> f64 {
    gas.mass_kg()
        * gas.state().state().specific_gas_constant
        * (2.5 * (t1 - t0) + 0.0005 * (t1 * t1 - t0 * t0))
}

#[test]
fn g1_joint_exchange_uses_energy_midpoints_and_pays_both_native_material_works() {
    let (card, caloric) = material(0.0, 800.0);
    let declaration = conduction();
    let interface = layers(&card, &declaration);
    for t0 in [280.0, 360.0] {
        let mut live = joint_source(&card, t0);
        for _ in 0..20 {
            live.step(input_drive()).unwrap();
        }
        let before = live.clone();
        let mut frozen = before.clone();
        let ordinary = frozen.step_with_force(input_drive(), 2e-5).unwrap();
        let mut gas = parcel(t0);
        let initial_gas = gas.clone();
        let mut state = initial(caloric.curve());
        let original_body = state;
        let input = body(caloric.curve(), state, true, 0.001);
        let external_gas_heat = 0.0002;
        let actual = live
            .step_with_finite_gas_wall_heat_controlled(
                &mut gas,
                &mut state,
                input,
                &interface,
                input_drive(),
                2e-5,
                external_gas_heat,
                joint_controls(),
                &mut || false,
            )
            .unwrap();

        let tg = gas.state().state().temperature;
        let tb = state.temperature_k();
        let du = gas_energy_change(&initial_gas, t0, tg);
        let dh = input.mass_kg * CP * (tb - original_body.temperature_k());
        let q = actual.exchanged_heat_j;
        assert_eq!(q.is_sign_positive(), t0 < 300.0);
        assert!(q.abs() > 1e-6);
        assert!((du + actual.gas.refresh.parameter_work_j - external_gas_heat - q).abs() < 1e-11);
        assert!(
            (dh + actual.wall.parameter_work_j + q - input.internal_heat_j - actual.wall_heat_j)
                .abs()
                < 1e-11
        );

        // Gas midpoint halves INTERNAL ENERGY; body midpoint halves ENTHALPY.
        let tgm = actual.gas_midpoint.state().temperature;
        let tbm = actual.body_midpoint.temperature_k();
        assert!((gas_energy_change(&initial_gas, t0, tgm) - 0.5 * du).abs() < 1e-11);
        assert!((input.mass_kg * CP * (tbm - 300.0) - 0.5 * dh).abs() < 1e-11);
        assert!(
            (tgm - 0.5 * (t0 + tg)).abs() > 1e-8,
            "temperature averaging must not replace the nonlinear gas energy midpoint"
        );
        let kg = actual.gas_midpoint.state().thermal_conductivity;
        let kb = 3.0 + 0.01 * tbm;
        let g = declaration.area_m2
            / (declaration.gas_distance_m / kg + declaration.body_distance_m / kb);
        assert!((actual.interface.conductance_w_k - g).abs() < 1e-13);
        assert!((q - DT * g * (tbm - tgm)).abs() <= joint_controls().exchange.source_tolerance_j);
        assert_eq!(actual.interface.source.card_identity(), card.content_hash());
        assert_eq!(actual.interface.source.query().axes()["lot"], 7.0);
        assert_eq!(actual.interface.source.query().axes()["T"], tbm);

        let motion = frozen
            .load()
            .duct()
            .compliant_wall_state(frozen.duct_state())
            .unwrap();
        let expected_wall_work = spring_energy(&motion.displacement_m, true) * SLOPE * (tb - 300.0);
        assert!((actual.wall.parameter_work_j - expected_wall_work).abs() < 1e-14);
        let native_work = live.stored_energy_j() - frozen.stored_energy_j();
        assert!(
            (native_work - actual.gas.refresh.parameter_work_j - actual.wall.parameter_work_j)
                .abs()
                < 1e-14
        );
        assert_eq!(
            actual.coupled.parameter_work_j,
            actual.gas.refresh.parameter_work_j + actual.wall.parameter_work_j
        );
        let terminal_work: f64 = actual.coupled.terminals.iter().map(|t| t.port_work_j).sum();
        let balance = du + dh + live.stored_energy_j() - before.stored_energy_j()
            + actual.exported_dissipation_j
            + terminal_work
            - external_gas_heat
            - input.internal_heat_j
            - actual.coupled.upstream_work_j
            - actual.coupled.body_work_j
            - actual.coupled.aperture.mechanical_work_j;
        assert!(balance.abs() < 1e-11);
        assert!(
            (actual.exported_dissipation_j + actual.wall_heat_j
                - ordinary.inferred_total_dissipation_j)
                .abs()
                < 1e-15
        );
        assert!(actual.wall_heat_j > 0.0);
        assert!(actual.combined_energy_residual_j.abs() < 1e-11);
        assert!(actual.exchange_residual_j.abs() <= joint_controls().exchange.source_tolerance_j);
        assert!(actual.material_trials > 0);
        assert!(actual.material_trials <= joint_controls().maximum_material_trials);
        assert_eq!(actual.coupled.aperture, ordinary.aperture);
        assert_eq!(live.aperture().state(), frozen.aperture().state());
        assert_eq!(live.aperture().accepted_steps(), 21);
        assert_eq!(
            live.load().terminal().unwrap().pole_states(),
            frozen.load().terminal().unwrap().pole_states()
        );
        assert_eq!(gas.mass_kg(), initial_gas.mass_kg());
        assert_eq!(gas.volume_m3(), initial_gas.volume_m3());
        assert_eq!(
            live.load().duct().gas_field().unwrap().uniform_state(),
            Some(gas.state())
        );
        assert_eq!(
            live.aperture().spec().density_kg_m3,
            gas.state().state().density
        );
        for wall in sources(live.load()) {
            assert_eq!(wall.material().card_identity(), card.content_hash());
            assert_eq!(wall.material().query().axes()["T"], tb);
        }
    }
}

fn joint_program() -> FiniteGasWallHeatProgram {
    FiniteGasWallHeatProgram {
        controls: [(0, 0.0, 0.0), (7, 6.0, 600.0), (19, -3.0, -300.0)]
            .into_iter()
            .map(|(sample, gas_power_w, wall_power_w)| GasWallHeatControl {
                sample,
                gas_power_w,
                wall_power_w,
            })
            .collect(),
        coupling: joint_controls(),
    }
}

fn joint_player(card: &MaterialCard, block: usize) -> AperturePerformance {
    joint_player_config(
        card,
        AperturePerformanceConfig {
            samples: 24,
            ..config(block)
        },
    )
}

fn joint_player_config(
    card: &MaterialCard,
    settings: AperturePerformanceConfig,
) -> AperturePerformance {
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
        CoupledAperture::Dg(joint_source(card, 280.0)),
        ApertureObservation::Inlet,
        phrase,
        settings,
    )
    .unwrap()
}

fn attach(
    card: &Arc<MaterialCard>,
    curve: &EquilibriumEnthalpyPhaseCurve,
    block: usize,
    program: FiniteGasWallHeatProgram,
    gate: Arc<CancelGate>,
) -> AperturePerformance {
    joint_player(card, block)
        .with_finite_gas_wall_heat(
            parcel(280.0),
            body(curve, initial(curve), true, 0.0),
            card.clone(),
            "T".into(),
            conduction(),
            program,
            gate,
        )
        .unwrap()
}

#[test]
fn g3_scheduled_joint_powers_match_direct_steps_and_all_audio_block_sizes() {
    let (card, caloric) = material(0.0, 800.0);
    let card = Arc::new(card);
    let declaration = conduction();
    let interface = layers(&card, &declaration);
    let program = joint_program();
    let gate = Arc::new(CancelGate::new());
    let mut rendered = attach(&card, caloric.curve(), 1, program.clone(), gate.clone());
    assert!(Arc::ptr_eq(
        rendered.finite_gas_wall_cancel_gate().unwrap(),
        &gate
    ));
    assert_eq!(rendered.finite_gas_wall_conduction(), Some(&declaration));
    assert_eq!(rendered.finite_gas_wall_heat_program(), Some(&program));
    let mut direct = joint_source(&card, 280.0);
    let mut gas = parcel(280.0);
    let mut state = initial(caloric.curve());
    let mut expected_pressure = [0.0; 24];
    for (sample, value) in expected_pressure.iter_mut().enumerate() {
        let event = program
            .controls
            .iter()
            .rev()
            .find(|e| e.sample <= sample as u64)
            .unwrap();
        let input = body(caloric.curve(), state, true, event.wall_power_w * DT);
        let expected = direct
            .step_with_finite_gas_wall_heat_controlled(
                &mut gas,
                &mut state,
                input,
                &interface,
                input_drive(),
                0.0,
                event.gas_power_w * DT,
                program.coupling,
                &mut || false,
            )
            .unwrap();
        rendered.block(core::slice::from_mut(value)).unwrap();
        let actual = rendered.last_finite_gas_wall_step().unwrap();
        assert_eq!(*value, expected.coupled.aperture.bore_pressure_pa);
        assert_eq!(actual.coupled, expected.coupled);
        assert_eq!(actual.gas.refresh, expected.gas.refresh);
        assert_eq!(actual.wall, expected.wall);
        assert_eq!(actual.external_gas_heat_j, event.gas_power_w * DT);
        assert_eq!(actual.external_body_heat_j, event.wall_power_w * DT);
        assert_eq!(actual.exchanged_heat_j, expected.exchanged_heat_j);
        assert_eq!(rendered.finite_gas_wall_gas(), Some(&gas));
        assert_eq!(rendered.finite_gas_wall_state(), Some(&state));
        equal_source(retained(&rendered), &direct);
        let count = program
            .controls
            .iter()
            .take_while(|e| e.sample <= sample as u64)
            .count();
        assert_eq!(
            rendered.applied_gas_wall_heat_controls(),
            &program.controls[..count]
        );
    }
    for block in [7, 24] {
        let mut other = attach(
            &card,
            caloric.curve(),
            block,
            program.clone(),
            Arc::new(CancelGate::new()),
        );
        let mut values = [0.0; 24];
        for chunk in values.chunks_mut(block) {
            other.block(chunk).unwrap();
        }
        assert_eq!(values, expected_pressure);
        assert_eq!(other.finite_gas_wall_gas(), Some(&gas));
        assert_eq!(other.finite_gas_wall_state(), Some(&state));
        equal_source(retained(&other), &direct);
    }
}

#[test]
fn g4_joint_refusals_and_final_cancellation_publish_neither_inventory_nor_history() {
    let (card, caloric) = material(0.0, 800.0);
    let declaration = conduction();
    let interface = layers(&card, &declaration);
    let mut before = joint_source(&card, 280.0);
    for _ in 0..12 {
        before.step(input_drive()).unwrap();
    }
    let original_gas = parcel(280.0);
    let original_body = initial(caloric.curve());
    let input = body(caloric.curve(), original_body, true, 0.001);
    let mut clean = before.clone();
    let mut clean_gas = original_gas.clone();
    let mut clean_body = original_body;
    let mut polls = 0;
    let expected = clean
        .step_with_finite_gas_wall_heat_controlled(
            &mut clean_gas,
            &mut clean_body,
            input,
            &interface,
            input_drive(),
            0.0,
            0.0002,
            joint_controls(),
            &mut || {
                polls += 1;
                false
            },
        )
        .unwrap();
    for failure in 0..4 {
        let mut attempted = before.clone();
        let mut gas = original_gas.clone();
        let mut state = original_body;
        let mut limits = joint_controls();
        let mut source = conduction();
        match failure {
            0 => limits.maximum_material_trials = 1,
            1 => limits.exchange.heat_bracket_j = [0.0, 0.0],
            2 => source.context = point(300.0, 9.0),
            _ => {}
        }
        let interface = layers(&card, &source);
        let mut count = 0;
        let result = attempted.step_with_finite_gas_wall_heat_controlled(
            &mut gas,
            &mut state,
            input,
            &interface,
            input_drive(),
            0.0,
            0.0002,
            limits,
            &mut || {
                count += 1;
                failure == 3 && count == polls
            },
        );
        assert!(result.is_err());
        if failure == 3 {
            assert!(matches!(result, Err(AcousticRealizeError::Cancelled)));
            assert_eq!(count, polls);
        }
        assert_eq!(gas, original_gas);
        assert_eq!(state, original_body);
        equal_source(&attempted, &before);
        let retry = attempted
            .step_with_finite_gas_wall_heat_controlled(
                &mut gas,
                &mut state,
                input,
                &layers(&card, &declaration),
                input_drive(),
                0.0,
                0.0002,
                joint_controls(),
                &mut || false,
            )
            .unwrap();
        assert_eq!(retry.coupled, expected.coupled);
        assert_eq!(retry.gas.refresh, expected.gas.refresh);
        assert_eq!(retry.wall, expected.wall);
        assert_eq!(gas, clean_gas);
        assert_eq!(state, clean_body);
        equal_source(&attempted, &clean);
    }
}

#[test]
fn g4_playback_keeps_paired_control_cursors_at_the_last_accepted_sample() {
    let (card, caloric) = material(0.0, 800.0);
    let card = Arc::new(card);
    let bind = |settings, program| {
        joint_player_config(&card, settings).with_finite_gas_wall_heat(
            parcel(280.0),
            body(caloric.curve(), initial(caloric.curve()), true, 0.0),
            card.clone(),
            "T".into(),
            conduction(),
            program,
            Arc::new(CancelGate::new()),
        )
    };
    // A paired power event consumes two assignments alongside the mouth control.
    let settings = AperturePerformanceConfig {
        samples: 24,
        max_controls: 7,
        ..config(24)
    };
    assert!(bind(settings, joint_program()).is_ok());
    assert!(
        bind(
            AperturePerformanceConfig {
                max_controls: 6,
                ..settings
            },
            joint_program()
        )
        .is_err()
    );
    for powers in [(f64::from_bits(1), 0.0), (0.0, f64::NAN)] {
        let mut invalid = joint_program();
        invalid.controls[0].gas_power_w = powers.0;
        invalid.controls[0].wall_power_w = powers.1;
        assert!(bind(settings, invalid).is_err());
    }
    let mut reference = attach(
        &card,
        caloric.curve(),
        7,
        joint_program(),
        Arc::new(CancelGate::new()),
    );
    reference.block(&mut [0.0; 7]).unwrap();
    for cancel in [false, true] {
        let gate = Arc::new(CancelGate::new());
        let mut program = joint_program();
        if !cancel {
            program.controls[1].gas_power_w = 1e9;
        }
        let mut attempted = attach(&card, caloric.curve(), 7, program, gate.clone());
        attempted.block(&mut [0.0; 7]).unwrap();
        if cancel {
            gate.request();
        }
        let mut output = [123.0];
        assert!(attempted.block(&mut output).is_err());
        assert_eq!(output, [123.0]);
        assert_eq!(attempted.samples_rendered(), 7);
        assert_eq!(
            attempted.finite_gas_wall_gas(),
            reference.finite_gas_wall_gas()
        );
        assert_eq!(
            attempted.finite_gas_wall_state(),
            reference.finite_gas_wall_state()
        );
        assert_eq!(
            attempted.applied_gas_wall_heat_controls(),
            reference.applied_gas_wall_heat_controls()
        );
        assert_eq!(attempted.applied_controls(), reference.applied_controls());
        assert_eq!(
            attempted.last_finite_gas_wall_step().unwrap().coupled,
            reference.last_finite_gas_wall_step().unwrap().coupled
        );
        equal_source(retained(&attempted), retained(&reference));
        assert!(matches!(
            attempted.block(&mut output),
            Err(RenderError::Poisoned)
        ));
    }
}
