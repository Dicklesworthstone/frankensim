//! Real material-source integration for locally nonisothermal molecular flux.
//! All diffusivity and Cp values are authored verification data, not measurements.

use super::*;
use fs_flux::miscibility::{
    BinaryMiscibilityState,
    nonisothermal::{
        BinaryNonisothermalError, BinaryNonisothermalFace, BinaryNonisothermalField,
        BinaryNonisothermalSolve,
    },
};
use fs_thermochem::{
    MolarThermalQuantityV1,
    solution::caloric::{BinaryLiquidCaloricModelV1, PureLiquidHeatCapacityV1},
};

fn calorics() -> BinaryLiquidCaloricModelV1 {
    // Equal constant pure Cp removes partial-enthalpy contrast, leaving the
    // explicitly authored Q as the only thermal drive for uniform composition.
    let capacities = ["A", "B"].map(|name| {
        PureLiquidHeatCapacityV1::new(
            SpeciesId::new(name).unwrap(),
            Temperature::new(300.0),
            [300.0, 600.0]
                .into_iter()
                .map(|t| (Temperature::new(t), MolarThermalQuantityV1::new(100.0)))
                .collect(),
        )
        .unwrap()
    });
    BinaryLiquidCaloricModelV1::new(pair(0.0), capacities).unwrap()
}

fn bind(
    net: &BinaryMiscibilityNetwork,
    state: &BinaryMiscibilityState,
    temperatures: [f64; 2],
) -> Result<BinaryNonisothermalField, BinaryNonisothermalError> {
    net.bind_nonisothermal_field(
        state,
        calorics(),
        temperatures.map(Temperature::new).to_vec(),
        vec![BinaryNonisothermalFace {
            // One chemical substage only: its D is frozen at the actual face
            // state, so a constant-D network provides an independent oracle.
            thermal_conductivity_w_m_k: 0.0,
            heat_of_transport_j_per_mol: 1500.0,
        }],
        [Temperature::new(300.0), Temperature::new(600.0)],
        1 << 20,
        || false,
    )
}

#[test]
fn local_temperature_source_matches_its_actual_face_diffusivity_and_retains_claim_pins() {
    let model = gradient(300.0, 10_000.0, 0.0, 0.0);
    let bare = network(model.clone(), 0.0);
    let fractions = vec![[0.4, 0.6]; 2];
    let seed = bare
        .state_from_mole_fractions(fractions.clone(), || false)
        .unwrap();
    let source = kinetics("T", vec![(300.0, 1e-5), (600.0, 4e-5)], 1.0);
    let (net, state) = bare
        .bind_diffusivity_source(&seed, source.clone(), || false)
        .unwrap()
        .into_parts();
    let temperatures = [350.0, 550.0];
    let resolved = source
        .resolve_nonisothermal_face(
            temperatures.map(Temperature::new),
            [fractions[0], fractions[1]],
            || false,
        )
        .unwrap();
    near(resolved.diffusivity_m2_s, 2.5e-5);
    near(net.faces().next().unwrap().diffusivity_m2_s, 1e-5);
    for (response, temperature) in [
        (&resolved.left, 350.0),
        (&resolved.face, 450.0),
        (&resolved.right, 550.0),
    ] {
        assert_eq!(response.state_coordinate("T"), Some(temperature));
        assert_eq!(response.state_coordinate("pressure"), Some(100_000.0));
        assert_eq!(
            response
                .property(source.property().name())
                .unwrap()
                .answer()
                .receipt
                .selected,
            source.selected_claim()
        );
    }

    let field = bind(&net, &state, temperatures).unwrap();
    let before = field.clone();
    let oracle_net = network(model, resolved.diffusivity_m2_s);
    let oracle_state = oracle_net
        .state_from_mole_fractions(fractions, || false)
        .unwrap();
    let oracle = bind(&oracle_net, &oracle_state, temperatures).unwrap();
    let controls = BinaryNonisothermalSolve::default();
    let expected = oracle.step_candidate(0.1, controls, || false).unwrap();
    let actual = field.step_candidate(0.1, controls, || false).unwrap();

    // Positive Q drives A toward the colder left cell despite uniform initial
    // composition. This would vanish if local temperatures were replaced by
    // the network's 300 K reference, and its magnitude detects stale D(300 K).
    assert!(actual.step().transferred_a_mol[0] < 0.0);
    assert_eq!(
        actual.field().component_moles(),
        expected.field().component_moles()
    );
    assert_eq!(
        actual.field().temperatures(),
        expected.field().temperatures()
    );
    assert_eq!(actual.step(), expected.step());
    assert_eq!(actual.field().network().diffusivity_source(), Some(&source));
    assert_eq!(field, before);
}

#[test]
fn an_unsupported_local_temperature_refuses_initialization_and_the_complete_candidate() {
    let model = gradient(300.0, 10_000.0, 0.0, 0.0);
    let bare = network(model.clone(), 2e-5);
    let fractions = vec![[0.8, 0.2], [0.2, 0.8]];
    let seed = bare
        .state_from_mole_fractions(fractions.clone(), || false)
        .unwrap();
    // The enclosing validity box extends to 600 K, but the original source
    // curve supports only 300..500 K. A supported face mean cannot authorize
    // an unsupported hot endpoint.
    let source = kinetics("T", vec![(300.0, 2e-5), (500.0, 2e-5)], 1.0);
    let (net, state) = bare
        .bind_diffusivity_source(&seed, source.clone(), || false)
        .unwrap()
        .into_parts();
    let state_before = state.clone();
    assert!(matches!(
        bind(&net, &state, [350.0, 550.0]),
        Err(BinaryNonisothermalError::Spatial(Error::DiffusivitySource(
            _
        )))
    ));
    assert_eq!(state, state_before);

    let field = bind(&net, &state, [350.0, 500.0]).unwrap();
    let before = field.clone();
    let oracle_net = network(model, 2e-5);
    let oracle_state = oracle_net
        .state_from_mole_fractions(fractions, || false)
        .unwrap();
    let oracle = bind(&oracle_net, &oracle_state, [350.0, 500.0]).unwrap();
    let controls = BinaryNonisothermalSolve::default();
    let valid_numeric_step = oracle.step_candidate(0.1, controls, || false).unwrap();
    // Strong composition drive transports A and its positive Q to the hotter
    // right cell. Calorics admits this endpoint, but the pinned D source does not.
    assert!(valid_numeric_step.field().temperatures()[1].value() > 500.0);
    assert!(matches!(
        field.step_candidate(0.1, controls, || false),
        Err(BinaryNonisothermalError::Spatial(Error::DiffusivitySource(
            _
        )))
    ));
    assert_eq!(field, before);
    assert_eq!(field.component_moles(), state.component_moles());
    assert_eq!(field.network().diffusivity_source(), Some(&source));
    assert_eq!(field.elapsed_s(), state.elapsed_s());
}
