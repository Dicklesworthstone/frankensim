//! Mutual drag must conserve momentum before it supplies current and heat.
use super::*;
use fs_qty::{Amount, MolarMass, Temperature, Volume};
use fs_thermochem::equilibrium::composition::{FrozenIonizedGasStateV1, IdealGasEquilibriumV1};
use fs_thermochem::{
    ChargeVector, ElementId, ElementalMatrix, Nasa9RegionV1, Nasa9StandardStateModelV1,
    required_nasa9_parameter_dims,
};

fn anchored(state: &FrozenIonizedGasStateV1, ion_rate: f64) -> IonizedGasConductivity {
    let th = state.heavy_temperature().value();
    let te = state.electron_temperature().value();
    let p = state.pressure().value();
    IonizedGasConductivity::from_state(
        state.clone(),
        vec![
            source_two_temperature("e-", th, te, p, 1e10),
            source_two_temperature("H+", th, te, p, ion_rate),
        ],
        limits(),
        || false,
    )
    .unwrap()
}

#[test]
fn g1_mutual_coulomb_drag_matches_independent_two_fluid_balance_and_heat() {
    let state = frozen_state(0.001, 6000.0, 12000.0, 1.0);
    let background = anchored(&state, 1e6);
    let gas = background
        .clone()
        .with_coulomb_drag(coulomb_model(&state), || false)
        .unwrap();
    let electron = gas.carriers().iter().find(|c| c.is_electron()).unwrap();
    let ion = gas.carriers().iter().find(|c| !c.is_electron()).unwrap();
    let pair = &gas.coulomb_pairs()[0];
    let ae =
        electron.number_density_per_m3() * electron.mass_kg() * electron.relaxation().rate_per_s();
    let ai = ion.number_density_per_m3() * ion.mass_kg() * ion.relaxation().rate_per_s();
    let k = pair.friction_coefficient_kg_m3_s();
    let b = electron.number_density_per_m3() * E;
    // Exact two-fluid inverse for b_i=-b_e=b, independently of Cholesky.
    let determinant = ae * ai + k * (ae + ai);
    let ue = -b * ai / determinant;
    let ui = b * ae / determinant;
    let sigma = b * b * (ae + ai) / determinant;
    relative(electron.drift_mobility_m2_v_s(), ue, 2e-12);
    relative(ion.drift_mobility_m2_v_s(), ui, 2e-12);
    relative(gas.conductivity_s_m(), sigma, 2e-12);
    assert!(gas.conductivity_s_m() < background.conductivity_s_m());
    let incorrect_independent_drag = b * b / (ae + k) + b * b / (ai + k);
    assert!((sigma / incorrect_independent_drag - 1.0).abs() > 0.05);

    let field = [0.01, -0.02];
    let e2 = field[0] * field[0] + field[1] * field[1];
    let w2 = (ue - ui).powi(2) * e2;
    let total_mass = electron.mass_kg() + ion.mass_kg();
    let correction = pair.thermal_conductance_w_m3_k()
        * (state.electron_temperature().value() - state.heavy_temperature().value())
        * w2
        / (2.0 * pair.relative_thermal_variance_m2_s2());
    let friction_e = ae * ue * ue * e2 + k * w2 * ion.mass_kg() / total_mass;
    let friction_h = ai * ui * ui * e2 + k * w2 * electron.mass_kg() / total_mass;
    let heat = gas.joule_heating_w_m3(field).unwrap();
    relative(
        heat.electron_friction_power_density_w_m3(),
        friction_e,
        2e-12,
    );
    relative(heat.heavy_friction_power_density_w_m3(), friction_h, 2e-12);
    relative(
        heat.thermal_exchange_correction_to_electrons_w_m3(),
        correction,
        2e-12,
    );
    relative(
        heat.electron_power_density_w_m3(),
        friction_e + correction,
        2e-12,
    );
    relative(
        heat.heavy_power_density_w_m3(),
        friction_h - correction,
        2e-12,
    );
    relative(heat.total_power_density_w_m3(), sigma * e2, 2e-12);
    assert!(
        (heat.electron_power_density_w_m3() - electron.conductivity_s_m() * e2).abs()
            > 0.05 * heat.total_power_density_w_m3()
    );
    assert!(gas.coulomb_force_balance_relative_residual().unwrap() < 1e-12);
    assert!(gas.coulomb_power_balance_relative_residual().unwrap() < 1e-12);
    assert!(
        gas.joule_heating_w_m3([2.0 * gas.maximum_electric_field_v_m(), 0.0])
            .is_err()
    );
}

fn multiple_ions() -> FrozenIonizedGasStateV1 {
    let original = chemistry();
    let mut card = original.models()[1].regions()[0].card().clone();
    card.provenance.source = "manufactured monatomic X+ carrier for mutual momentum transport; no empirical ionization or level-statistics claim".into();
    let keys = required_nasa9_parameter_dims();
    card.parameters.get_mut(keys[7].0).unwrap().value = 0.0;
    card.parameters.get_mut(keys[8].0).unwrap().value += 1.5 * (4.0 * MH / (MH - ME)).ln();
    let extra = Nasa9StandardStateModelV1::new(
        id("X+"),
        MolarMass::new(4.0 * MH * NA),
        original.convention().clone(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap();
    let models = vec![
        original.models()[0].clone(),
        original.models()[1].clone(),
        extra,
        original.models()[2].clone(),
    ];
    let ids = models
        .iter()
        .map(|m| m.species().clone())
        .collect::<Vec<_>>();
    let bank = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("H").unwrap(), ElementId::new("X").unwrap()],
            ids.clone(),
            vec![vec![1, 1, 0, 0], vec![0, 0, 1, 0]],
        )
        .unwrap(),
        ChargeVector::new(ids, vec![0, 1, 1, -1]).unwrap(),
        models,
    )
    .unwrap();
    FrozenIonizedGasStateV1::new(
        &bank,
        &[
            Amount::new(0.9998),
            Amount::new(0.0001),
            Amount::new(0.0001),
            Amount::new(0.0002),
        ],
        Volume::new(1.0),
        Temperature::new(6000.0),
        Temperature::new(12000.0),
        1e-12,
        || false,
    )
    .unwrap()
}

#[test]
fn g1_full_pair_transport_includes_ion_ion_entrainment() {
    let state = multiple_ions();
    let th = state.heavy_temperature().value();
    let te = state.electron_temperature().value();
    let p = state.pressure().value();
    let rates = vec![
        source_two_temperature("H+", th, te, p, 1e6),
        source_two_temperature("X+", th, te, p, 1e8),
        source_two_temperature("e-", th, te, p, 1e10),
    ];
    let gas = IonizedGasConductivity::from_state(state.clone(), rates, limits(), || false)
        .unwrap()
        .with_coulomb_drag(coulomb_model(&state), || false)
        .unwrap();
    assert_eq!(gas.coulomb_pairs().len(), 3);
    let carriers = gas.carriers(); // H+, X+, e-.
    let anchor = carriers
        .iter()
        .map(|c| c.number_density_per_m3() * c.mass_kg() * c.relaxation().rate_per_s())
        .collect::<Vec<_>>();
    let force = carriers
        .iter()
        .map(|c| c.number_density_per_m3() * c.charge_number() as f64 * E)
        .collect::<Vec<_>>();
    let coefficient = |first, second| {
        gas.coulomb_pairs()
            .iter()
            .find(|p| p.first_index() == first && p.second_index() == second)
            .unwrap()
            .friction_coefficient_kg_m3_s()
    };
    let ion_electron = [coefficient(1, 3), coefficient(2, 3)];
    assert!(coefficient(1, 2) > 0.0);
    // An independently eliminated electron-ion star intentionally omits
    // ion-ion drag. The complete source must change that relative velocity.
    let mut numerator = force[2];
    let mut denominator = anchor[2];
    for i in 0..2 {
        let fraction = ion_electron[i] / (anchor[i] + ion_electron[i]);
        numerator += fraction * force[i];
        denominator += fraction * anchor[i];
    }
    let ue = numerator / denominator;
    let without_ion_ion = [
        (force[0] + ion_electron[0] * ue) / (anchor[0] + ion_electron[0]),
        (force[1] + ion_electron[1] * ue) / (anchor[1] + ion_electron[1]),
        ue,
    ];
    let relative_before = (without_ion_ion[0] - without_ion_ion[1]).abs();
    let relative_after =
        (carriers[0].drift_mobility_m2_v_s() - carriers[1].drift_mobility_m2_v_s()).abs();
    assert!(relative_after < 0.95 * relative_before);
    let star_conductivity = force
        .iter()
        .zip(without_ion_ion)
        .map(|(b, u)| b * u)
        .sum::<f64>();
    assert!(gas.conductivity_s_m() < star_conductivity);
    relative(
        gas.joule_heating_w_m3([0.01, 0.0])
            .unwrap()
            .total_power_density_w_m3(),
        1e-4 * gas.conductivity_s_m(),
        2e-12,
    );
}

#[test]
fn g0_source_binding_and_signed_drift_exchange_remain_physical() {
    let state = frozen_state(0.0001, 2000.0, 20000.0, 1.0);
    let model = coulomb_model(&state);
    let background = anchored(&state, 10.0);
    assert_eq!(
        background.clone().with_coulomb_drag(model.clone(), || true),
        Err(IonizedGasError::Cancelled)
    );
    let changed = frozen_state(0.0002, 2000.0, 20000.0, 1.0);
    assert!(
        anchored(&changed, 10.0)
            .with_coulomb_drag(model.clone(), || false)
            .is_err()
    );
    let gas = background
        .with_coulomb_drag(model.clone(), || false)
        .unwrap();
    let heat = gas.joule_heating_w_m3([0.01, 0.0]).unwrap();
    assert!(heat.heavy_friction_power_density_w_m3() > 0.0);
    assert!(heat.heavy_power_density_w_m3() < 0.0);
    assert!(heat.thermal_exchange_correction_to_electrons_w_m3() > 0.0);
    relative(
        heat.total_power_density_w_m3(),
        1e-4 * gas.conductivity_s_m(),
        2e-12,
    );
    // Temperatures may evolve while the exact population/model owner remains.
    let refreshed = frozen_state(0.0001, 6000.0, 12000.0, 1.0);
    let next = anchored(&refreshed, 10.0)
        .with_coulomb_drag(model, || false)
        .unwrap();
    assert_eq!(next.state(), &refreshed);
    assert_eq!(next.coulomb_model().unwrap().initial_state(), &state);
    assert_ne!(
        next.coulomb_pairs()[0].friction_coefficient_kg_m3_s(),
        gas.coulomb_pairs()[0].friction_coefficient_kg_m3_s()
    );
}
