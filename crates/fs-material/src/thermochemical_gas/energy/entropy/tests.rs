//! Independent synthetic caloric identities, not measured gas qualification.
use super::*;
use crate::thermochemical_gas::tests::{close, model, point};

fn controls() -> GasEntropySolve {
    GasEntropySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 1500.0,
        entropy_tolerance_j_kg_k: 1e-9, maximum_iterations: 128 }
}

#[test]
fn density_entropy_roundtrip_preserves_complete_material_evaluation() {
    let model = model(vec![0.75, 0.25], 0.001);
    let context = point(600.0, 200000.0).with("material-lot", 12.0).unwrap();
    let original = model.evaluate(&context).unwrap();
    let result = model.reconstruct_density_entropy(&context, original.state().density,
        original.specific_entropy_j_kg_k().unwrap(), controls(), || false).unwrap();
    close(result.state().state().temperature, 600.0);
    close(result.state().state().pressure, 200000.0);
    close(result.state().state().dynamic_viscosity, original.state().dynamic_viscosity);
    close(result.state().state().thermal_conductivity, original.state().thermal_conductivity);
    assert!(result.residual_j_kg_k().abs() <= controls().entropy_tolerance_j_kg_k);
    assert_eq!(result.state().components().len(), 2);
    let energy = model.reconstruct_density_energy(&context, result.state().state().density,
        result.state().specific_internal_energy_j_kg(), super::super::GasEnergySolve {
            minimum_temperature_k: 200.0, maximum_temperature_k: 1500.0,
            energy_tolerance_j_kg: 1e-7, maximum_iterations: 128,
        }, || false).unwrap();
    close(energy.state().specific_entropy_j_kg_k().unwrap(), original.specific_entropy_j_kg_k().unwrap());
}

#[test]
fn variable_heat_capacity_uses_integrated_entropy_not_fixed_gamma() {
    let model = model(vec![0.75, 0.25], 0.001);
    let context = point(300.0, 100000.0);
    let initial = model.evaluate(&context).unwrap();
    // cp/R = a + b*T. Independently integrating ds at fixed composition gives
    // ln(rho1/rho0) = (a-1)*ln(T1/T0) + b*(T1-T0).
    let (a, b) = (3.25, 0.00075);
    let density_ratio = ((a - 1.0) * (600.0_f64/300.0).ln() + b * 300.0).exp();
    let result = model.reconstruct_density_entropy(&context, initial.state().density*density_ratio,
        initial.specific_entropy_j_kg_k().unwrap(), controls(), || false).unwrap();
    close(result.state().state().temperature, 600.0);
    close(result.state().state().pressure, 100000.0*density_ratio*2.0);
    assert!((300.0*density_ratio.powf(initial.state().gamma - 1.0) - 600.0).abs() > 10.0);
    assert!(result.state().state().gamma < initial.state().gamma);
}

#[test]
fn isentropic_pressure_derivative_matches_shared_acoustic_sound_speed() {
    let model = model(vec![0.25, 0.75], 0.001);
    let context = point(500.0, 100000.0);
    let initial = model.evaluate(&context).unwrap();
    let rho = initial.state().density;
    let epsilon = 1e-4;
    let solve = |density| model.reconstruct_density_entropy(&context, density,
        initial.specific_entropy_j_kg_k().unwrap(), controls(), || false).unwrap();
    let lo = solve(rho*(1.0-epsilon));
    let hi = solve(rho*(1.0+epsilon));
    let dp_drho = (hi.state().state().pressure-lo.state().state().pressure)/(2.0*rho*epsilon);
    let c2 = initial.state().sound_speed.powi(2);
    assert!((dp_drho/c2 - 1.0).abs() < 1e-7);
}

#[test]
fn bad_states_brackets_work_limits_and_cancellation_refuse() {
    let model = model(vec![0.75, 0.25], 0.001);
    let context = point(300.0, 100000.0);
    let initial = model.evaluate(&context).unwrap();
    let rho = initial.state().density;
    let s = initial.specific_entropy_j_kg_k().unwrap();
    assert!(matches!(model.reconstruct_density_entropy(&context, rho, s, controls(), || true),
        Err(GasEntropyError::Cancelled)));
    assert!(matches!(model.reconstruct_density_entropy(&context, rho, s+1e9, controls(), || false),
        Err(GasEntropyError::OutsideBracket { .. })));
    assert!(matches!(model.reconstruct_density_entropy(&context, rho, s,
        GasEntropySolve { maximum_iterations: 1, ..controls() }, || false),
        Err(GasEntropyError::NotConverged { .. })));
    for density in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(model.reconstruct_density_entropy(&context, density, s, controls(), || false).is_err());
    }
    for bad in [GasEntropySolve { maximum_iterations: 0, ..controls() },
        GasEntropySolve { entropy_tolerance_j_kg_k: 0.0, ..controls() },
        GasEntropySolve { minimum_temperature_k: 1600.0, ..controls() }] {
        assert!(model.reconstruct_density_entropy(&context, rho, s, bad, || false).is_err());
    }
    let mut polls = 0;
    assert!(matches!(model.reconstruct_density_entropy(&context, rho, s, controls(), || {
        polls += 1; polls == 4
    }), Err(GasEntropyError::Cancelled)));
}

fn with_entropy_reference(offset: f64) -> ThermochemicalGasModel {
    use fs_thermochem::{Nasa9RegionV1, Nasa9StandardStateModelV1};
    use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
    let original = model(vec![0.75, 0.25], 0.001);
    let components = original.thermodynamics().components().iter().map(|source| {
        let regions = source.regions().iter().map(|region| {
            let mut card = region.card().clone();
            card.parameters.get_mut("a8").unwrap().value += offset;
            Nasa9RegionV1::from_card(card).unwrap()
        }).collect();
        Nasa9StandardStateModelV1::new(source.species().clone(), source.molar_mass(),
            source.convention().clone(), regions).unwrap()
    }).collect();
    let mixture = FrozenIdealGasMixtureModelV1::new(components,
        original.thermodynamics().declared_composition().clone()).unwrap();
    ThermochemicalGasModel::new(mixture, original.transport.clone(), None).unwrap()
}

#[test]
fn g1_isentropic_reconstruction_is_independent_of_absolute_entropy_reference() {
    let plain = with_entropy_reference(0.0);
    let shifted = with_entropy_reference(1e24);
    let context = point(500.0, 100000.0);
    let density = plain.evaluate(&context).unwrap().state().density * 1.2;
    let solve = |source: &ThermochemicalGasModel| source.reconstruct_density_entropy_change(
        &context, density, 0.0, controls(), || false).unwrap();
    let a = solve(&plain);
    let b = solve(&shifted);
    assert_eq!(a.state().state(), b.state().state());
    assert_eq!(a.residual_j_kg_k(), b.residual_j_kg_k());
    assert!(b.state().state().temperature > 500.0);
    // Absolute source entropy deliberately cannot resolve the actual change.
    assert_eq!(shifted.evaluate(&point(600.0,100000.0)).unwrap().specific_entropy_j_kg_k().unwrap(),
        shifted.evaluate(&context).unwrap().specific_entropy_j_kg_k().unwrap());
    let mut parcel = super::super::parcel::ClosedGasParcel::new(shifted, context, 0.01).unwrap();
    let before = parcel.clone();
    parcel.set_volume_isentropic(0.005, controls(), || false).unwrap();
    let returned = parcel.set_volume_isentropic(0.01, controls(), || false).unwrap();
    close(parcel.state().state().temperature, before.state().state().temperature);
    assert!(returned.entropy_change_j_k.abs() <= parcel.mass_kg()*controls().entropy_tolerance_j_kg_k);
}

#[test]
fn g1_nozzle_choking_and_subsonic_flow_ignore_entropy_reference_offsets() {
    use super::super::parcel::ClosedGasParcel;
    use super::super::parcel::transfer::nozzle::GasNozzle;
    let parcel = |offset| ClosedGasParcel::new(with_entropy_reference(offset),
        point(500.0,200000.0),0.01).unwrap();
    let plain = parcel(0.0);
    let shifted = parcel(1e24);
    let nozzle = GasNozzle { area_m2: 1e-5, discharge_coefficient: 0.8 };
    let solve = super::super::GasEnergySolve { minimum_temperature_k:200.0,
        maximum_temperature_k:1500.0,energy_tolerance_j_kg:1e-7,maximum_iterations:128 };
    for back_pressure in [180000.0,50000.0] {
        let a = plain.nozzle_flow_to(back_pressure,nozzle,solve,||false).unwrap();
        let b = shifted.nozzle_flow_to(back_pressure,nozzle,solve,||false).unwrap();
        assert_eq!(a.mass_flow_kg_s,b.mass_flow_kg_s);
        assert_eq!(a.velocity_m_s,b.velocity_m_s);
        assert_eq!(a.choked,b.choked);
        assert_eq!(a.throat().state(),b.throat().state());
        assert!(b.entropy_energy_residual_j_kg.abs() <= 2.0*solve.energy_tolerance_j_kg);
    }
}

#[test]
fn g0_relative_entropy_admits_changes_but_retains_domain_and_cancellation_gates() {
    let source = with_entropy_reference(1e24);
    let context = point(500.0,100000.0);
    let initial = source.evaluate(&context).unwrap();
    let desired = source.specific_entropy_change_j_kg_k(&context,&point(600.0,120000.0)).unwrap();
    let solved = source.reconstruct_density_entropy_change(&context,initial.state().density,
        desired,controls(),||false).unwrap();
    close(solved.state().state().temperature,600.0);
    assert!(matches!(source.reconstruct_density_entropy_change(&context,initial.state().density,
        0.0,controls(),||true),Err(GasEntropyError::Cancelled)));
    assert!(source.reconstruct_density_entropy_change(&point(190.0,100000.0),
        initial.state().density,0.0,controls(),||false).is_err());
    assert!(source.reconstruct_density_entropy_change(&context,initial.state().density,
        f64::NAN,controls(),||false).is_err());
}

#[test]
fn g1_midpoint_heating_retains_small_entropy_gain_above_large_reference() {
    use super::super::parcel::ClosedGasParcel;
    let parcel = |offset| ClosedGasParcel::new(with_entropy_reference(offset),
        point(500.0,100000.0),0.01).unwrap();
    let solve = super::super::GasEnergySolve { minimum_temperature_k:200.0,
        maximum_temperature_k:1500.0,energy_tolerance_j_kg:1e-7,maximum_iterations:128 };
    let heat = 0.25;
    let a = parcel(0.0).midpoint_heat_volume(0.01,heat,solve,||false).unwrap();
    let b = parcel(1e24).midpoint_heat_volume(0.01,heat,solve,||false).unwrap();
    assert_eq!(a.endpoint().state().state(),b.endpoint().state().state());
    assert_eq!(a.step().entropy_change_j_k,b.step().entropy_change_j_k);
    let entropy = b.step().entropy_change_j_k;
    let final_temperature = b.endpoint().state().state().temperature;
    assert!(entropy > 0.0);
    assert!(entropy >= heat/final_temperature - 1e-12);
    assert!(entropy <= heat/500.0 + 1e-12);
}
