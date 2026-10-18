use super::*;
use crate::thermochemical_gas::tests::{model, point};

fn primitive(temperature_k: f64, fraction: f64) -> PrimitiveGasState {
    PrimitiveGasState {
        temperature_k,
        composition: Composition::new(
            CompositionBasis::MoleFraction,
            vec![fraction, 1.0 - fraction],
        )
        .unwrap(),
    }
}
fn field(states: Vec<PrimitiveGasState>) -> AxialGasField {
    AxialGasField::new(
        vec![0.0, 1.0],
        101325.0,
        states,
        model(vec![0.75, 0.25], 0.001),
        QueryPoint::new(),
        ThermochemicalPhase::FrozenGas,
    )
    .unwrap()
}

#[test]
fn g0_uniform_field_is_bitwise_single_state() {
    let state = primitive(300.0, 0.75);
    let field = field(vec![state.clone(), state]);
    let direct = field.model().evaluate(&point(300.0, 101325.0)).unwrap();
    for x in [0.0, 0.125, 0.5, 0.8, 1.0] {
        assert_eq!(field.state_at(x).unwrap(), direct);
    }
    assert_eq!(
        field.sample_states(&[0.0, 0.25, 1.0]).unwrap(),
        vec![*direct.state(); 3]
    );
    assert_eq!(field.uniform_state(), Some(&direct));
}

#[test]
fn g0_interpolates_primitive_state_only() {
    let field = field(vec![primitive(300.0, 0.75), primitive(600.0, 0.25)]);
    let middle = field.state_at(0.5).unwrap();
    let changed = field
        .model()
        .with_composition(primitive(450.0, 0.5).composition)
        .unwrap();
    assert_eq!(
        changed.transport_materials(),
        field.model().transport_materials()
    );
    assert_eq!(
        changed.thermodynamics().components(),
        field.model().thermodynamics().components()
    );
    let direct = changed.evaluate(&point(450.0, 101325.0)).unwrap();
    assert_eq!(middle, direct);
    assert_eq!(middle.state().temperature, 450.0);
    assert_eq!(
        middle.thermodynamics().properties().molar_mass().value(),
        0.016
    );
    middle.state().validate_coherence().unwrap();
    let left = field.state_at(0.0).unwrap();
    let right = field.state_at(1.0).unwrap();
    let rho_average = 0.5 * (left.state().density + right.state().density);
    let sound_average = 0.5 * (left.state().sound_speed + right.state().sound_speed);
    assert!((middle.state().density / rho_average - 1.0).abs() > 0.01);
    assert!((middle.state().sound_speed / sound_average - 1.0).abs() > 0.01);
}

#[test]
fn g0_species_can_enter_from_zero_without_invented_trace_gas() {
    let field = field(vec![primitive(300.0, 1.0), primitive(400.0, 0.0)]);
    let left = field.state_at(0.0).unwrap();
    let middle = field.state_at(0.5).unwrap();
    let right = field.state_at(1.0).unwrap();
    assert_eq!(left.components().len(), 1);
    assert_eq!(right.components().len(), 1);
    assert_eq!(middle.components().len(), 2);
    assert_eq!(
        left.thermodynamics().properties().molar_mass().value(),
        0.028
    );
    assert_eq!(
        right.thermodynamics().properties().molar_mass().value(),
        0.004
    );
    assert_eq!(
        middle.thermodynamics().properties().molar_mass().value(),
        0.016
    );
    assert_eq!(field.model().transport_materials().len(), 2);
    let direct = field
        .model()
        .with_composition(primitive(350.0, 0.5).composition)
        .unwrap()
        .evaluate(&point(350.0, 101325.0))
        .unwrap();
    assert_eq!(middle, direct);
}

#[test]
fn g0_station_and_sampling_refusals() {
    let state = primitive(300.0, 0.75);
    let make = |positions: Vec<f64>, states: Vec<PrimitiveGasState>| {
        AxialGasField::new(
            positions,
            101325.0,
            states,
            model(vec![0.75, 0.25], 0.001),
            QueryPoint::new(),
            ThermochemicalPhase::FrozenGas,
        )
    };
    for positions in [
        vec![],
        vec![0.0],
        vec![0.1, 1.0],
        vec![0.0, 0.0],
        vec![0.0, 1.0, 0.5],
        vec![0.0, f64::NAN],
    ] {
        assert!(make(positions.clone(), vec![state.clone(); positions.len()]).is_err());
    }
    assert!(make(vec![0.0, 1.0], vec![state.clone()]).is_err());
    assert!(make(vec![0.0, 1.0], vec![state.clone(), primitive(199.0, 0.75)]).is_err());
    let mut bad = state.clone();
    bad.composition = Composition::new(CompositionBasis::MassFraction, vec![0.75, 0.25]).unwrap();
    assert!(make(vec![0.0, 1.0], vec![state.clone(), bad]).is_err());
    let mut bad = state.clone();
    bad.composition = Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap();
    assert!(make(vec![0.0, 1.0], vec![state.clone(), bad]).is_err());
    let positions = (0..MAX_FIELD_STATIONS).map(|i| i as f64).collect();
    let at_cap = make(positions, vec![state.clone(); MAX_FIELD_STATIONS]).unwrap();
    assert!(
        make(
            (0..=MAX_FIELD_STATIONS).map(|i| i as f64).collect(),
            vec![state.clone(); MAX_FIELD_STATIONS + 1]
        )
        .is_err()
    );
    assert!(at_cap.sample_states(&vec![0.0; MAX_FIELD_STATIONS]).is_ok());
    assert!(
        at_cap
            .sample_states(&vec![0.0; MAX_FIELD_STATIONS + 1])
            .is_err()
    );
    for position in [
        -f64::EPSILON,
        f64::NAN,
        f64::INFINITY,
        at_cap.length_m() + 1e-9,
    ] {
        assert!(at_cap.state_at(position).is_err());
    }
}

#[test]
fn g0_source_conditions_and_phase_policy_survive_sampling() {
    let state = primitive(300.0, 0.75);
    let make = |pressure, conditions| {
        AxialGasField::new(
            vec![0.0, 1.0],
            pressure,
            vec![state.clone(), state.clone()],
            model(vec![0.75, 0.25], 0.001),
            conditions,
            ThermochemicalPhase::FrozenGas,
        )
    };
    assert!(make(0.0, QueryPoint::new()).is_err());
    assert!(make(1e7 + 1.0, QueryPoint::new()).is_err());
    assert!(
        make(
            101325.0,
            QueryPoint::new().with("pressure", 101325.0).unwrap()
        )
        .is_err()
    );
    assert!(
        make(
            101325.0,
            QueryPoint::new().with("temperature", 300.0).unwrap()
        )
        .is_err()
    );
    let original = model(vec![0.99, 0.01], 0.0);
    let water = fs_thermochem::SpeciesId::new("B").unwrap();
    let wet = ThermochemicalGasModel::new(
        original.thermodynamics().clone(),
        original.transport_materials().to_vec(),
        Some(&water),
    )
    .unwrap();
    let states = vec![primitive(300.0, 0.99), primitive(300.0, 0.9)];
    assert!(
        AxialGasField::new(
            vec![0.0, 1.0],
            101325.0,
            states.clone(),
            wet.clone(),
            QueryPoint::new(),
            ThermochemicalPhase::FrozenGas
        )
        .is_ok()
    );
    let error = AxialGasField::new(
        vec![0.0, 1.0],
        101325.0,
        states,
        wet,
        QueryPoint::new(),
        ThermochemicalPhase::LiquidWaterLimit,
    )
    .unwrap_err();
    assert!(error.to_string().contains("supersaturated"));
}
