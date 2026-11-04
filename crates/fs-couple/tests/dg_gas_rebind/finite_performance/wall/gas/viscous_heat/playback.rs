//! The original renderer consumes the same joint dissipation/thermal step.
use super::*;

fn program() -> FiniteGasWallHeatProgram {
    FiniteGasWallHeatProgram {
        controls: [(0, 0.0, 0.0), (2, 6.0, 600.0), (4, -3.0, -300.0)]
            .into_iter()
            .map(|(sample, gas_power_w, wall_power_w)| GasWallHeatControl {
                sample,
                gas_power_w,
                wall_power_w,
            })
            .collect(),
        coupling: molecular_controls(),
    }
}

fn player(
    card: &Arc<MaterialCard>,
    curve: &EquilibriumEnthalpyPhaseCurve,
    block: usize,
    gate: Arc<CancelGate>,
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
        CoupledAperture::Dg(molecular_source(card)),
        ApertureObservation::Inlet,
        phrase,
        AperturePerformanceConfig {
            samples: 6,
            ..config(block)
        },
    )
    .unwrap()
    .with_finite_gas_wall_heat(
        parcel(280.0),
        body(curve, initial(curve), true, 0.0),
        card.clone(),
        "T".into(),
        conduction(),
        program(),
        gate,
    )
    .unwrap()
}

#[test]
fn scheduled_joint_viscous_heat_matches_direct_physics_and_audio_block_partitions() {
    let (card, caloric) = material(0.0, 800.0);
    let card = Arc::new(card);
    let declaration = conduction();
    let interface = layers(&card, &declaration);
    let gate = Arc::new(CancelGate::new());
    let original = player(&card, caloric.curve(), 1, gate.clone());
    assert_eq!(
        original.finite_gas_wall_viscous_molecular_heat_budget(),
        None
    );
    let before = retained(&original).clone();
    let original_gas = original.finite_gas_wall_gas().unwrap().clone();
    let original_body = *original.finite_gas_wall_state().unwrap();
    let mut rendered = original
        .with_finite_gas_wall_viscous_molecular_heat(budgets())
        .unwrap();
    assert_eq!(
        rendered.finite_gas_wall_viscous_molecular_heat_budget(),
        Some(budgets())
    );
    assert!(Arc::ptr_eq(
        rendered.finite_gas_wall_cancel_gate().unwrap(),
        &gate
    ));
    assert_eq!(rendered.finite_gas_wall_heat_program(), Some(&program()));
    assert_eq!(rendered.finite_gas_wall_gas(), Some(&original_gas));
    assert_eq!(rendered.finite_gas_wall_state(), Some(&original_body));
    assert!(rendered.applied_gas_wall_heat_controls().is_empty());
    assert_eq!(rendered.samples_rendered(), 0);
    equal_source(retained(&rendered), &before);
    let mut direct = before;
    let mut gas = original_gas;
    let mut state = original_body;
    let schedule = program();
    let mut pressure = [0.0; 6];
    for (sample, output) in pressure.iter_mut().enumerate() {
        let event = schedule
            .controls
            .iter()
            .rev()
            .find(|event| event.sample <= sample as u64)
            .unwrap();
        let input = body(caloric.curve(), state, true, event.wall_power_w * DT);
        let expected = direct
            .step_with_finite_gas_wall_viscous_molecular_heat_controlled(
                &mut gas,
                &mut state,
                input,
                &interface,
                input_drive(),
                0.0,
                event.gas_power_w * DT,
                schedule.coupling,
                budgets(),
                &mut || false,
            )
            .unwrap();
        rendered.block(core::slice::from_mut(output)).unwrap();
        let actual = rendered.last_finite_gas_wall_step().unwrap();
        assert_eq!(*output, expected.coupled.aperture.bore_pressure_pa);
        assert_eq!(actual.coupled, expected.coupled);
        assert_eq!(actual.gas.refresh, expected.gas.refresh);
        assert_eq!(actual.wall, expected.wall);
        assert_eq!(actual.gas.viscous_heat_j, expected.gas.viscous_heat_j);
        assert_eq!(actual.gas.molecular_heat_j, expected.gas.molecular_heat_j);
        assert!(actual.gas.viscous_heat_j > 0.0 && actual.gas.molecular_heat_j > 0.0);
        assert!(actual.wall_heat_j > 0.0);
        assert_eq!(actual.exchanged_heat_j, expected.exchanged_heat_j);
        assert_eq!(
            actual.exported_dissipation_j,
            expected.exported_dissipation_j
        );
        assert_eq!(actual.external_gas_heat_j, event.gas_power_w * DT);
        assert_eq!(actual.external_body_heat_j, event.wall_power_w * DT);
        assert_eq!(rendered.finite_gas_wall_gas(), Some(&gas));
        assert_eq!(rendered.finite_gas_wall_state(), Some(&state));
        equal_source(retained(&rendered), &direct);
    }
    assert_eq!(
        rendered.applied_gas_wall_heat_controls(),
        schedule.controls.as_slice()
    );
    for block in [2, 6] {
        let mut other = player(&card, caloric.curve(), block, Arc::new(CancelGate::new()))
            .with_finite_gas_wall_viscous_molecular_heat(budgets())
            .unwrap();
        let mut values = [0.0; 6];
        for chunk in values.chunks_mut(block) {
            other.block(chunk).unwrap();
        }
        assert_eq!(values, pressure);
        assert_eq!(other.finite_gas_wall_gas(), Some(&gas));
        assert_eq!(other.finite_gas_wall_state(), Some(&state));
        assert_eq!(
            other.applied_gas_wall_heat_controls(),
            schedule.controls.as_slice()
        );
        equal_source(retained(&other), &direct);
    }
}

#[test]
fn joint_heat_playback_admits_source_budgets_before_samples_and_keeps_cancelled_cursors() {
    let (card, caloric) = material(0.0, 800.0);
    let card = Arc::new(card);
    let denied = FiniteGasDgHeatBudget {
        molecular: heat_budget(),
        viscous: no_budget(),
    };
    assert!(
        player(&card, caloric.curve(), 1, Arc::new(CancelGate::new()))
            .with_finite_gas_wall_viscous_molecular_heat(denied)
            .is_err()
    );
    let gate = Arc::new(CancelGate::new());
    let cold = player(&card, caloric.curve(), 1, gate.clone());
    gate.request();
    assert!(
        cold.with_finite_gas_wall_viscous_molecular_heat(budgets())
            .is_err()
    );

    let gate = Arc::new(CancelGate::new());
    let mut live = player(&card, caloric.curve(), 2, gate.clone())
        .with_finite_gas_wall_viscous_molecular_heat(budgets())
        .unwrap();
    live.block(&mut [0.0; 2]).unwrap();
    let before = retained(&live).clone();
    let gas = live.finite_gas_wall_gas().unwrap().clone();
    let state = *live.finite_gas_wall_state().unwrap();
    let frame = live.last_finite_gas_wall_step().unwrap().coupled.clone();
    let controls = live.applied_gas_wall_heat_controls().to_vec();
    gate.request();
    let mut output = [123.0];
    assert!(live.block(&mut output).is_err());
    assert_eq!(output, [123.0]);
    assert_eq!(live.samples_rendered(), 2);
    assert_eq!(live.finite_gas_wall_gas(), Some(&gas));
    assert_eq!(live.finite_gas_wall_state(), Some(&state));
    assert_eq!(live.last_finite_gas_wall_step().unwrap().coupled, frame);
    assert_eq!(live.applied_gas_wall_heat_controls(), controls.as_slice());
    equal_source(retained(&live), &before);
    assert!(matches!(
        live.block(&mut output),
        Err(RenderError::Poisoned)
    ));
}
