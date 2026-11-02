use super::*;
use crate::thermochemical_gas::tests::{model, point, species, transport};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;

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

fn ternary() -> ThermochemicalGasModel {
    ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![species("A", 0.028, 3.5, 0.0), species("B", 0.004, 2.5, 0.0), species("C", 0.044, 4.0, 0.0)],
            Composition::new(CompositionBasis::MoleFraction, vec![0.5, 0.25, 0.25]).unwrap(),
        )
        .unwrap(),
        vec![transport("A", 1.7e-5, 0.025), transport("B", 2e-5, 0.15), transport("C", 1.5e-5, 0.02)],
        None,
    )
    .unwrap()
}

#[test]
fn g0_interpolated_fractions_that_one_residual_add_cannot_close_are_closed_not_refused() {
    // Two authored three-species stations (each sums to exactly 1.0) and the interpolation weight at which the old
    // residual loop (the residual added to the largest fraction up to four times, no final check) left sum bits
    // 0x3ff0000000000001, so the exact-unit-sum owner refused a sample between two exact stations.
    let a = vec![0.33462649384269494, 0.6298488836822311, 0.03552462247507393];
    let b = vec![0.36395405838550077, 0.3959793728157154, 0.24006656879878385];
    let w = 0.14331662586371663;
    assert_eq!(a.iter().sum::<f64>(), 1.0);
    assert_eq!(b.iter().sum::<f64>(), 1.0);
    let station = |temperature_k: f64, fractions: &[f64]| PrimitiveGasState {
        temperature_k,
        composition: Composition::new(CompositionBasis::MoleFraction, fractions.to_vec()).unwrap(),
    };
    let field = AxialGasField::new(
        vec![0.0, 1.0],
        101325.0,
        vec![station(300.0, &a), station(400.0, &b)],
        ternary(),
        QueryPoint::new(),
        ThermochemicalPhase::FrozenGas,
    )
    .unwrap();
    // Stations 0 and 1 make the interpolation weight exactly the sample position.
    let state = field.state_at(w).expect("an interpolation between exact stations must not be refused");
    let closed: Vec<f64> = state.thermodynamics().receipt().declared_fraction_bits().iter()
        .map(|bits| f64::from_bits(*bits)).collect();
    assert_eq!(closed.len(), 3);
    assert_eq!(closed.iter().sum::<f64>(), 1.0);
    let raw: Vec<f64> = a.iter().zip(&b).map(|(&a, &b)| a + w * (b - a)).collect();
    let moved = closed.iter().zip(&raw).filter(|(closed, raw)| closed != raw).count();
    assert!((1..=2).contains(&moved), "closure moves one or two entries: {closed:?} against {raw:?}");
    // The closure removes a sum error of about one ulp of unity, so it moves entries by that much in ABSOLUTE terms (here
    // 0.5 and 0.62 eps in two entries); relative to a small fraction that is a larger multiple of eps, which is why the
    // bound is absolute.
    for (closed, raw) in closed.iter().zip(&raw) {
        assert!((closed - raw).abs() <= 4.0 * f64::EPSILON, "{closed:e} moved from {raw:e}");
    }
    // Authored stations are still returned exactly, unrepaired.
    assert_eq!(field.state_at(0.0).unwrap().thermodynamics().receipt().declared_fraction_bits()
        .iter().map(|bits| f64::from_bits(*bits)).collect::<Vec<_>>(), a);
}
