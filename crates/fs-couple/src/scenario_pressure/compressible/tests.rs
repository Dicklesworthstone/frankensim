use super::*;
use fs_evidence::ValidityDomain;
use fs_gas::dynamics::{Boundaries, Boundary, CaloricDomain, Reconstruction};
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{Length, MolarMass, QtyAny};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{
    Composition, CompositionBasis, ElementalReferenceIdV1, Nasa9RegionV1,
    Nasa9StandardStateModelV1, ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1,
    StandardStatePhaseV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K, required_nasa9_parameter_dims,
};

#[test]
fn chirped_pressure_sweep_joins_are_breakpoints() {
    let signal = TimeSignal::Chirp(fs_scenario::signal::LinearChirp::new(
        1.0, 2.0, 2.0, 7.0, 0.0,
        QtyAny::new(100_000.0, Pressure::DIMS), QtyAny::new(100.0, Pressure::DIMS),
    ).unwrap());
    assert_eq!(next_breakpoint(&signal, 0.0), Some(1.0));
    assert_eq!(next_breakpoint(&signal, 1.0), Some(2.0));
    assert_eq!(next_breakpoint(&signal, 2.0), None);
}

// Manufactured SI source: R=1 J/(kg K), Cp(T)=3.5+slope*T and
// h(T)=3.5*T+slope*T^2/2. These are analytic fixtures, not material measurements.
fn gas(slope: f64) -> FrozenGas {
    let mut a = [0.0; 9];
    a[2] = 3.5;
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
        validity: ValidityDomain::unconstrained().with("T", 0.05, 50.0),
        sources: Vec::new(),
        provenance: Provenance {
            source: "Manufactured pressure-history caloric source".into(),
            license: "test fixture".into(),
            artifact: None,
        },
    };
    let convention = StandardStateConventionV1::new(
        StandardStatePhaseV1::Gas,
        ReferenceEquationOfStateV1::IdealGas,
        Pressure::new(1.0),
        ElementalReferenceIdV1::new("pressure-history-fixture").unwrap(),
    )
    .unwrap();
    let species = Nasa9StandardStateModelV1::new(
        SpeciesId::new("pressure-history-gas").unwrap(),
        MolarMass::new(UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K),
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
        CaloricDomain::new(Temperature::new(0.05), Temperature::new(50.0)),
        || false,
    )
    .unwrap()
}

fn closed() -> Boundaries {
    Boundaries::Ends {
        left: Boundary::Reflecting,
        right: Boundary::Reflecting,
    }
}

fn field(source: &FrozenGas, boundaries: Boundaries) -> Euler1d {
    let rest = source.conserved(1.0, 0.0, 1.0).unwrap();
    Euler1d::new(
        source.clone(),
        Length::new(0.5),
        0.2,
        vec![rest; 6],
        boundaries,
        || false,
    )
    .unwrap()
}

fn controls() -> StepControl {
    StepControl {
        cfl: 0.45,
        reconstruction: Reconstruction::FirstOrder,
        maximum_halvings: 12,
    }
}

fn constant(value: f64, dims: Dims) -> TimeSignal {
    TimeSignal::Constant(QtyAny::new(value, dims))
}

fn ramp(end: f64, from: f64, to: f64, dims: Dims) -> TimeSignal {
    TimeSignal::Ramp {
        t_start: 0.0,
        t_end: end,
        from: QtyAny::new(from, dims),
        to: QtyAny::new(to, dims),
    }
}

fn history<'a>(pressure: &'a TimeSignal, temperature: &'a TimeSignal) -> ReservoirHistory<'a> {
    ReservoirHistory {
        pressure,
        pressure_reference: ReservoirPressureReference::Absolute,
        temperature,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-10 * (1.0 + expected.abs()),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn constant_history_reproduces_the_native_reservoir_and_retains_other_boundaries() {
    let source = gas(0.2);
    let pressure = constant(2.0, Pressure::DIMS);
    let temperature = constant(1.2, Temperature::DIMS);
    let reservoir = source.conserved(2.0 / 1.2, 0.0, 2.0).unwrap();
    let mut native = field(
        &source,
        Boundaries::Ends {
            left: Boundary::Reservoir(reservoir),
            right: Boundary::Reflecting,
        },
    );
    let mut driven = field(&source, closed());
    let fixed = native.step(Time::new(0.002), controls(), || false).unwrap();
    let report = step_euler(
        &mut driven,
        &source,
        [Some(history(&pressure, &temperature)), None],
        Time::new(0.002),
        controls(),
        || false,
    )
    .unwrap();
    for (actual, expected) in driven.states().iter().zip(native.states()) {
        close(actual.density, expected.density);
        close(actual.momentum, expected.momentum);
        close(actual.energy, expected.energy);
    }
    close(report.gas.boundary_change.mass, fixed.boundary_change.mass);
    close(
        report.gas.boundary_change.energy,
        fixed.boundary_change.energy,
    );
    assert_eq!(driven.boundaries(), closed());
    assert!(report.stages.iter().all(|stage| stage[1].is_none()));
    assert_eq!(report.gas.right_boundary_change.mass, 0.0);
    assert_eq!(report.gas.right_boundary_change.energy, 0.0);
}

#[test]
fn varying_pressure_and_temperature_deliver_the_original_stage_stagnation_enthalpy() {
    let slope = 0.2;
    let source = gas(slope);
    let mut driven = field(&source, closed());
    let pressure = ramp(0.02, 2.0, 4.0, Pressure::DIMS);
    let temperature = ramp(0.02, 1.0, 1.4, Temperature::DIMS);
    let report = step_euler(
        &mut driven,
        &source,
        [Some(history(&pressure, &temperature)), None],
        Time::new(0.005),
        controls(),
        || false,
    )
    .unwrap();
    let [first, second] = report.stages.map(|stage| stage[0].unwrap());
    for stage in [first, second] {
        let t = 1.0 + 0.4 * stage.time_s / 0.02;
        close(stage.absolute_pressure_pa, 2.0 + 2.0 * stage.time_s / 0.02);
        close(stage.prescribed_temperature_k, t);
        assert!(stage.flux.mass > 0.0);
        let physical_enthalpy = stage.flux.energy / stage.flux.mass + source.energy_reference();
        close(physical_enthalpy, 3.5 * t + 0.5 * slope * t * t);
    }
    assert_ne!(first.flux.energy, second.flux.energy);
    let factor = 0.5 * report.gas.elapsed.value() * driven.mesh().face_areas()[0];
    close(
        report.gas.left_boundary_change.mass,
        factor * (first.flux.mass + second.flux.mass),
    );
    close(
        report.gas.left_boundary_change.energy,
        factor * (first.flux.energy + second.flux.energy),
    );
    close(
        report.gas.after.mass - report.gas.before.mass,
        report.gas.boundary_change.mass,
    );
    close(
        report.gas.after.energy - report.gas.before.energy,
        report.gas.boundary_change.energy,
    );
    close(report.gas.conservation_residual.energy, 0.0);
    assert!(
        driven
            .gas()
            .primitive(driven.states()[0], || false)
            .unwrap()
            .pressure
            > 1.0
    );
}

#[test]
fn g1_oscillating_reservoirs_resolve_pressure_and_temperature_stage_clocks() {
    let source = gas(0.0);
    let mut driven = field(&source, closed());
    let tone = |frequency, offset, amplitude, dims| {
        TimeSignal::Chirp(fs_scenario::signal::LinearChirp::new(
            0.0, 0.01, frequency, frequency, 0.0,
            QtyAny::new(offset, dims), QtyAny::new(amplitude, dims),
        ).unwrap())
    };
    let pressure = tone(1_000.0, 2.0, 0.1, Pressure::DIMS);
    let temperature = tone(2_000.0, 1.0, 0.1, Temperature::DIMS);
    let fixed_pressure = constant(2.0, Pressure::DIMS);
    let fixed_temperature = constant(1.0, Temperature::DIMS);
    let report = step_euler(
        &mut driven,
        &source,
        [Some(history(&pressure, &fixed_temperature)),
            Some(history(&fixed_pressure, &temperature))],
        Time::new(0.002),
        controls(),
        || false,
    ).unwrap();

    // The requested interval contains whole periods of both sources. Native
    // spatial CFL alone would accept it, aliasing both histories to their peaks.
    // The faster temperature source on the OTHER end must also limit the step.
    let dt = report.gas.elapsed.value();
    assert!(dt > 0.0 && dt * 2_000.0 <= DEFAULT_MAXIMUM_CYCLES_PER_STEP);
    close(dt, DEFAULT_MAXIMUM_CYCLES_PER_STEP / 2_000.0);
    assert_eq!(report.gas.halvings, 0);
    assert_eq!(driven.time().value(), dt);
    for end in 0..2 {
        assert_eq!(report.stages[0][end].unwrap().time_s, 0.0);
        assert_eq!(report.stages[1][end].unwrap().time_s, dt);
    }
    let left = report.stages[1][0].unwrap();
    let right = report.stages[1][1].unwrap();
    close(left.absolute_pressure_pa,
        2.0 + 0.1 * (core::f64::consts::TAU * 1_000.0 * dt).cos());
    close(right.prescribed_temperature_k,
        1.0 + 0.1 * (core::f64::consts::TAU * 2_000.0 * dt).cos());
    assert!(left.absolute_pressure_pa < report.stages[0][0].unwrap().absolute_pressure_pa);
    assert!(right.prescribed_temperature_k < report.stages[0][1].unwrap().prescribed_temperature_k);
    close(report.gas.after.mass - report.gas.before.mass, report.gas.boundary_change.mass);
    close(report.gas.after.energy - report.gas.before.energy, report.gas.boundary_change.energy);
    close(report.gas.conservation_residual.energy, 0.0);
    assert_eq!(driven.boundaries(), closed());
}

#[test]
fn pressure_reversal_changes_the_entropy_donor_at_either_end() {
    let source = gas(0.0);
    let pressure = ramp(0.001, 0.7, 1.5, Pressure::DIMS);
    let temperature = constant(1.0, Temperature::DIMS);
    for end in 0..2 {
        let mut driven = field(&source, closed());
        let mut histories = [None; 2];
        histories[end] = Some(history(&pressure, &temperature));
        let report = step_euler(
            &mut driven,
            &source,
            histories,
            Time::new(0.001),
            controls(),
            || false,
        )
        .unwrap();
        let first = report.stages[0][end].unwrap().flux;
        let second = report.stages[1][end].unwrap().flux;
        let inward = if end == 0 { 1.0 } else { -1.0 };
        assert!(inward * first.mass < 0.0, "initial outflow at end {end}");
        assert!(inward * second.mass > 0.0, "final inflow at end {end}");
        close(second.energy / second.mass + source.energy_reference(), 3.5);
        close(report.gas.conservation_residual.mass, 0.0);
        close(report.gas.conservation_residual.energy, 0.0);
    }
}

#[test]
fn retries_resample_the_actual_shorter_stage_clock_and_discard_rejected_samples() {
    let source = gas(0.0);
    let mut driven = field(&source, closed());
    // The first stage is exactly at rest. The hot endpoint raises the face
    // acoustic speed enough to violate the initial cold-field CFL estimate.
    let pressure = ramp(0.12, 1.0, 10.0, Pressure::DIMS);
    let temperature = ramp(0.12, 1.0, 20.0, Temperature::DIMS);
    let report = step_euler(
        &mut driven,
        &source,
        [Some(history(&pressure, &temperature)), None],
        Time::new(0.12),
        controls(),
        || false,
    )
    .unwrap();
    assert!(report.gas.halvings > 0);
    let dt = report.gas.elapsed.value();
    assert!(dt < 0.12);
    let first = report.stages[0][0].unwrap();
    let second = report.stages[1][0].unwrap();
    assert_eq!(first.time_s, 0.0);
    assert_eq!(second.time_s, dt);
    close(first.absolute_pressure_pa, 1.0);
    close(second.absolute_pressure_pa, 1.0 + 9.0 * dt / 0.12);
    close(second.prescribed_temperature_k, 1.0 + 19.0 * dt / 0.12);
    assert_eq!(driven.time().value(), dt);
}

#[test]
fn a_short_held_pressure_pulse_starts_at_its_event_without_a_premature_impulse() {
    let source = gas(0.0);
    let mut driven = field(&source, closed());
    let before = driven.states().to_vec();
    let pressure = TimeSignal::Table {
        times: vec![0.0, 0.01, 0.02],
        values: vec![1.0, 2.0, 1.0],
        dims: Pressure::DIMS,
        interp: Interp::Hold,
    };
    let temperature = constant(1.0, Temperature::DIMS);
    let histories = [Some(history(&pressure, &temperature)), None];
    let quiet = step_euler(
        &mut driven,
        &source,
        histories,
        Time::new(0.1),
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(driven.time().value(), 0.01);
    assert_eq!(driven.states(), before);
    assert!(
        quiet
            .stages
            .iter()
            .all(|stage| stage[0].unwrap().absolute_pressure_pa == 1.0)
    );
    assert_eq!(quiet.gas.boundary_change.mass, 0.0);
    assert_eq!(quiet.gas.boundary_change.energy, 0.0);
    let pulse = step_euler(
        &mut driven,
        &source,
        histories,
        Time::new(0.1),
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(driven.time().value(), 0.02);
    assert!(
        pulse
            .stages
            .iter()
            .all(|stage| stage[0].unwrap().absolute_pressure_pa == 2.0)
    );
    assert!(pulse.gas.boundary_change.mass > 0.0);
    assert!(pulse.gas.boundary_change.energy > 0.0);
}

#[test]
fn gauge_and_absolute_pressure_histories_produce_the_same_physical_exchange() {
    let source = gas(0.0);
    let absolute = constant(1.2, Pressure::DIMS);
    let gauge = constant(0.2, Pressure::DIMS);
    let temperature = constant(1.0, Temperature::DIMS);
    let mut absolute_field = field(&source, closed());
    let mut gauge_field = field(&source, closed());
    let a = step_euler(
        &mut absolute_field,
        &source,
        [None, Some(history(&absolute, &temperature))],
        Time::new(0.001),
        controls(),
        || false,
    )
    .unwrap();
    let g = step_euler(
        &mut gauge_field,
        &source,
        [
            None,
            Some(ReservoirHistory {
                pressure_reference: ReservoirPressureReference::Gauge {
                    reference_pressure_pa: 1.0,
                },
                ..history(&gauge, &temperature)
            }),
        ],
        Time::new(0.001),
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(a, g);
    assert_eq!(absolute_field.states(), gauge_field.states());
}

#[test]
fn invalid_late_pressure_units_source_temperature_and_cancellation_preserve_the_field() {
    let source = gas(0.0);
    let mut driven = field(&source, closed());
    let before = driven.states().to_vec();
    let temperature = constant(1.0, Temperature::DIMS);
    let invalid_later = ramp(0.001, 1.0, -1.0, Pressure::DIMS);
    assert!(matches!(
        step_euler(
            &mut driven,
            &source,
            [Some(history(&invalid_later, &temperature)), None],
            Time::new(0.001),
            controls(),
            || false,
        ),
        Err(ReservoirHistoryError::Invalid(_))
    ));
    let pressure = constant(2.0, Pressure::DIMS);
    for bad_temperature in [
        constant(1.0, Pressure::DIMS),
        constant(60.0, Temperature::DIMS),
    ] {
        assert!(
            step_euler(
                &mut driven,
                &source,
                [Some(history(&pressure, &bad_temperature)), None],
                Time::new(0.001),
                controls(),
                || false,
            )
            .is_err()
        );
    }
    let wrong_source = gas(0.2);
    assert!(matches!(
        step_euler(
            &mut driven,
            &wrong_source,
            [Some(history(&pressure, &temperature)), None],
            Time::new(0.001),
            controls(),
            || false,
        ),
        Err(ReservoirHistoryError::Invalid(_))
    ));
    let mut calls = 0;
    assert!(matches!(
        step_euler(
            &mut driven,
            &source,
            [Some(history(&pressure, &temperature)), None],
            Time::new(0.001),
            controls(),
            || {
                calls += 1;
                calls == 8
            },
        ),
        Err(ReservoirHistoryError::Cancelled)
    ));
    assert_eq!(driven.states(), before);
    assert_eq!(driven.time().value(), 0.0);
    assert_eq!(driven.boundaries(), closed());
}

#[cfg(feature = "nonlinear-brass")]
#[test]
fn a_pressure_supply_drives_the_original_nonlinear_duct_and_geometry_observations() {
    use fs_duct::nonlinear::{InitialGasState, NonlinearDuct};
    use fs_duct::{Duct, Segment};
    let source = gas(0.0);
    let geometry = Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.01,
            length: 1.0,
        }],
    };
    let mut duct = NonlinearDuct::new(&geometry, source.clone(), &[8], closed(), |_| {
        Ok(InitialGasState {
            density_kg_m3: 1.0,
            velocity_m_s: 0.0,
            pressure_pa: 1.0,
        })
    })
    .unwrap();
    let before = duct.sample(0).unwrap();
    let pressure = constant(2.0, Pressure::DIMS);
    let temperature = constant(1.0, Temperature::DIMS);
    let report = step_duct(
        &mut duct,
        &source,
        [Some(history(&pressure, &temperature)), None],
        Time::new(0.001),
        controls(),
        || false,
    )
    .unwrap();
    let after = duct.sample(0).unwrap();
    assert!(after.pressure_pa > before.pressure_pa);
    assert_eq!(after.area_m2, before.area_m2);
    assert_eq!(after.axial_position_m, before.axial_position_m);
    assert!(report.gas.boundary_change.mass > 0.0);
    close(report.gas.conservation_residual.energy, 0.0);
    assert_eq!(duct.time().value(), report.gas.elapsed.value());
}
