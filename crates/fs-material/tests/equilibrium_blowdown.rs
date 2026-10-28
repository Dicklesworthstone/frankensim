//! Independent finite-equilibrium-reservoir and shared density-query regressions.

use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::rocket::equilibrium::blowdown::{
    EquilibriumBlowdownSolve, EquilibriumGasReservoir,
};
use fs_thermochem::equilibrium::composition::{EquilibriumEntropySolveV1, EquilibriumSolveV1, IdealGasEquilibriumV1};
use fs_thermochem::equilibrium::composition::isentrope::{EquilibriumIsentropeV1, density::EquilibriumDensitySolveV1};
use fs_thermochem::{ChargeVector, ElementalMatrix, ElementId, ElementalReferenceIdV1,
    Nasa9StandardStateModelV1, Nasa9RegionV1, ReferenceEquationOfStateV1,
    StandardStateConventionV1, StandardStatePhaseV1, SpeciesId, required_nasa9_parameter_dims};
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{Amount, MolarMass, Temperature, Pressure};
use std::{collections::BTreeMap, sync::Arc};
const R: f64 = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
fn id(s: &str) -> SpeciesId { SpeciesId::new(s).unwrap() }
fn model(name: &str, mass: f64, cp: f64, h: f64, s: f64) -> Nasa9StandardStateModelV1 {
    let mut parameters = BTreeMap::new();
    for ((name, dims), value) in required_nasa9_parameter_dims().into_iter().take(9)
        .zip([0.0, 0.0, cp, 0.0, 0.0, 0.0, 0.0, h, s]) {
        parameters.insert(name.into(), LawParameter { value, dims });
    }
    parameters.insert("reference_pressure".into(), LawParameter { value: 100000.0, dims: Pressure::DIMS });
    let card = ConstitutiveModelCard { law: LawId(fs_thermochem::NASA9_LAW_ID_V1.into()),
        law_version: fs_thermochem::NASA9_LAW_VERSION_V1, parameters,
        state_schema_version: fs_thermochem::NASA9_STATE_SCHEMA_VERSION_V1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 200.0, 6000.0), sources: Vec::new(),
        provenance: Provenance { source: "manufactured reacting reservoir; not hydrogen calibration".into(),
            license: "test fixture".into(), artifact: None } };
    Nasa9StandardStateModelV1::new(id(name), MolarMass::new(mass),
        StandardStateConventionV1::new(StandardStatePhaseV1::Gas, ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(100000.0), ElementalReferenceIdV1::new("reservoir-test").unwrap()).unwrap(),
        vec![Nasa9RegionV1::from_card(card).unwrap()]).unwrap()
}
fn chemistry(reacting: bool, shift: f64) -> Arc<IdealGasEquilibriumV1> {
    let models = if reacting { vec![model("H", 0.001, 2.5, 12000.0+shift, 1.0),
        model("H2", 0.002, 3.5, 2.0*shift, 0.0)] }
        else { vec![model("H2", 0.002, 3.5, 2.0*shift, 0.0)] };
    let ids: Vec<_> = models.iter().map(|m| m.species().clone()).collect();
    let counts = if reacting { vec![1, 2] } else { vec![2] };
    Arc::new(IdealGasEquilibriumV1::new(
        ElementalMatrix::new(vec![ElementId::new("H").unwrap()], ids.clone(), vec![counts]).unwrap(),
        ChargeVector::new(ids.clone(), vec![0; ids.len()]).unwrap(), models).unwrap())
}
fn thermo(scale: f64) -> EquilibriumEntropySolveV1 {
    EquilibriumEntropySolveV1 { equilibrium: EquilibriumSolveV1 {
        relative_balance_tolerance: 1e-12, dimensionless_potential_tolerance: 1e-11,
        ..EquilibriumSolveV1::default() }, minimum_temperature_k: 200.0,
        maximum_temperature_k: 6000.0, entropy_tolerance_j_k: 1e-9*scale,
        maximum_temperature_iterations: 100 }
}
fn controls() -> EquilibriumBlowdownSolve {
    EquilibriumBlowdownSolve { minimum_pressure_pa: 1000.0, relative_density_tolerance: 2e-9,
        mass_tolerance_kg: 1e-11, energy_tolerance_j: 1e-5, maximum_density_iterations: 80 }
}
fn reservoir(reacting: bool, shift: f64, scale: f64) -> EquilibriumGasReservoir {
    EquilibriumGasReservoir::new(chemistry(reacting, shift), vec![(id("H2"), Amount::new(scale))],
        0.02, Temperature::new(2000.0), Pressure::new(100000.0), thermo(scale), || false).unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance*b.abs().max(1e-30), "{a:.17e} != {b:.17e}");
}

#[test]
fn g1_finite_withdrawal_matches_perfect_gas_mass_temperature_pressure_and_energy() {
    let initial = reservoir(false, 0.0, 1.0);
    for fraction in [0.95_f64, 0.8, 0.6, 0.4] {
        let candidate = initial.exhaust_mass_adiabatic((1.0-fraction)*initial.mass_kg(), controls(), || false).unwrap();
        let next = candidate.endpoint();
        near(next.mass_kg(), fraction*initial.mass_kg(), 1e-14);
        near(next.thermodynamics().equilibrium().temperature().value(), 2000.0*fraction.powf(0.4), 2e-8);
        near(next.thermodynamics().equilibrium().pressure().value(), 100000.0*fraction.powf(1.4), 2e-8);
        let energy = initial.mass_kg() * 2.5*R/0.002 * 2000.0 * (1.0-fraction.powf(1.4));
        near(candidate.step().exhausted_enthalpy_j, energy, 3e-8);
        near(candidate.step().exhausted_enthalpy_j, initial.internal_energy_j()-next.internal_energy_j(), 3e-8);
        assert!(candidate.step().density_mass_residual_kg.abs() <= controls().mass_tolerance_kg);
        assert!(candidate.step().specific_entropy_residual_j_kg_k.abs() <= thermo(1.0).entropy_tolerance_j_k/0.002);
    }
}

#[test]
fn g1_reacting_depletion_matches_independent_dissociation_density_and_entropy_roots() {
    let initial = reservoir(true, 0.0, 1.0);
    let candidate = initial.exhaust_mass_adiabatic(0.4*initial.mass_kg(), controls(), || false).unwrap();
    let next = candidate.endpoint();
    // Independently solve K=exp(1/2-24000/T+3/2 ln T),
    // alpha=sqrt(K/(K+4p/p0)), then fixed S and rho/rho0=0.6.
    near(initial.mass_kg(), 0.00016824287253346126, 2e-8);
    near(next.thermodynamics().equilibrium().temperature().value(), 1880.4422619335887, 3e-8);
    near(next.thermodynamics().equilibrium().pressure().value(), 54720.104592822645, 3e-8);
    near(candidate.step().exhausted_enthalpy_j, 4989.9439267931675, 5e-8);
    let alpha = next.thermodynamics().equilibrium().species().iter()
        .find(|s| s.species == id("H")).unwrap().amount.value()/2.0;
    near(alpha, 0.38683501030922224, 3e-8);
    let atoms = |r: &EquilibriumGasReservoir| r.thermodynamics().equilibrium().species().iter()
        .map(|s| s.amount.value()*if s.species == id("H") { 1.0 } else { 2.0 }).sum::<f64>()
        *r.physical_amount_scale();
    near(atoms(next)/atoms(&initial), 0.6, 3e-8);
    assert!(next.thermodynamics().equilibrium().species().iter().all(|s| s.source.is_some()));
}

#[test]
fn g3_partitioned_withdrawals_keep_the_birth_isentrope_and_total_exported_energy() {
    let initial = reservoir(true, 0.0, 1.0);
    let amount = 0.4*initial.mass_kg();
    let direct = initial.exhaust_mass_adiabatic(amount, controls(), || false).unwrap();
    let mut current = initial.clone(); let mut energy = 0.0;
    for _ in 0..8 {
        let (next, step) = current.exhaust_mass_adiabatic(amount/8.0, controls(), || false).unwrap().into_parts();
        energy += step.exhausted_enthalpy_j;
        assert!(step.specific_entropy_residual_j_kg_k.abs() <= thermo(1.0).entropy_tolerance_j_k/0.002);
        current = next;
    }
    near(current.thermodynamics().equilibrium().temperature().value(),
        direct.endpoint().thermodynamics().equilibrium().temperature().value(), 3e-8);
    near(energy, direct.step().exhausted_enthalpy_j, 3e-8);
    near(current.mass_kg(), direct.endpoint().mass_kg(), 1e-14);
}

#[test]
fn g3_feed_scaling_and_element_energy_datum_do_not_change_physical_depletion() {
    let base = reservoir(true, 0.0, 1.0);
    let expected = base.exhaust_mass_adiabatic(0.2*base.mass_kg(), controls(), || false).unwrap();
    for scale in [1e-6, 1e6] {
        let tank = reservoir(true, 0.0, scale);
        let result = tank.exhaust_mass_adiabatic(0.2*tank.mass_kg(), controls(), || false).unwrap();
        near(result.endpoint().mass_kg(), expected.endpoint().mass_kg(), 3e-8);
        near(result.step().exhausted_enthalpy_j, expected.step().exhausted_enthalpy_j, 3e-8);
    }
    let shifted = reservoir(true, 100000.0, 1.0);
    let result = shifted.exhaust_mass_adiabatic(0.2*shifted.mass_kg(), controls(), || false).unwrap();
    near(result.endpoint().thermodynamics().equilibrium().temperature().value(),
        expected.endpoint().thermodynamics().equilibrium().temperature().value(), 3e-8);
    near(result.step().exhausted_enthalpy_j-expected.step().exhausted_enthalpy_j,
        result.step().expelled_mass_kg*2.0*R*100000.0/0.002, 1e-7);
}

#[test]
fn g4_refused_or_cancelled_withdrawals_preserve_all_physical_inventory() {
    let initial = reservoir(true, 0.0, 1.0); let before = initial.clone();
    for dm in [-1.0, f64::NAN, initial.mass_kg(), initial.mass_kg()*2.0, 1e-300] {
        assert!(initial.exhaust_mass_adiabatic(dm, controls(), || false).is_err());
    }
    let mut restricted = controls(); restricted.minimum_pressure_pa = 90000.0;
    assert!(initial.exhaust_mass_adiabatic(0.4*initial.mass_kg(), restricted, || false).is_err());
    for limit in [0, 20, 200] {
        let mut polls = 0;
        assert!(initial.exhaust_mass_adiabatic(0.4*initial.mass_kg(), controls(), || {
            polls += 1; polls > limit
        }).is_err());
    }
    assert_eq!(initial, before);
    assert_eq!(initial.exhaust_mass_adiabatic(0.0, controls(), || false).unwrap().endpoint(), &initial);
}

#[test]
fn g1_g4_density_queries_reuse_the_same_isentrope_and_invalidate_failed_publication() {
    let source = chemistry(false, 0.0);
    let mut path = EquilibriumIsentropeV1::new(&source, Temperature::new(2000.0), Pressure::new(100000.0),
        &[(id("H2"), Amount::new(1.0))], thermo(1.0), || false).unwrap();
    let rho0 = path.reference_state().response().density_kg_m3;
    let c = EquilibriumDensitySolveV1 { minimum_pressure_pa: 1000.0, maximum_pressure_pa: 200000.0,
        relative_density_tolerance: 2e-9, maximum_pressure_iterations: 80 };
    for fraction in [0.95_f64, 0.4, 1.2, 0.8] {
        let result = path.sample_density(rho0*fraction, c, || false).unwrap();
        near(result.sample.pressure.value(), 100000.0*fraction.powf(1.4), 2e-8);
        assert!(path.publish_current(|| false).is_ok());
    }
    assert!(path.sample_density(0.0, c, || false).is_err());
    assert!(path.current_sample().is_none());
    assert!(path.publish_current(|| false).is_err());
    assert!(path.sample_density(rho0*0.8, c, || false).is_ok());
    assert!(path.sample_density(rho0*0.6, c, || true).is_err());
    assert!(path.publish_current(|| false).is_err());
}
