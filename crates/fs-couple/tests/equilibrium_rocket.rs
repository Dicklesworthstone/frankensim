//! Reacting finite propellant driving native point translation, with independent oracles.
use fs_couple::rocket_blowdown::{GasRocketEnvironment, equilibrium::EquilibriumGasRocketState};
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::rocket::{
    RocketNozzle, equilibrium::{EquilibriumRocketNozzleSolve, blowdown::{
        EquilibriumBlowdownSolve, EquilibriumGasReservoir, time::EquilibriumRocketBlowdownSolve}},
};
use fs_mbd::Vec3;
use fs_thermochem::{ChargeVector, ElementalMatrix, ElementId, ElementalReferenceIdV1,
    Nasa9StandardStateModelV1, Nasa9RegionV1, ReferenceEquationOfStateV1,
    StandardStateConventionV1, StandardStatePhaseV1, SpeciesId, required_nasa9_parameter_dims,
    equilibrium::composition::{EquilibriumEntropySolveV1, EquilibriumSolveV1, IdealGasEquilibriumV1}};
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{Amount, MolarMass, Temperature, Pressure};
use std::{collections::BTreeMap, sync::Arc};
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
        provenance: Provenance { source: "manufactured equilibrium flight, not real hydrogen calibration".into(),
            license: "test fixture".into(), artifact: None } };
    Nasa9StandardStateModelV1::new(id(name), MolarMass::new(mass),
        StandardStateConventionV1::new(StandardStatePhaseV1::Gas, ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(100000.0), ElementalReferenceIdV1::new("flight-test").unwrap()).unwrap(),
        vec![Nasa9RegionV1::from_card(card).unwrap()]).unwrap()
}
fn controls() -> EquilibriumRocketBlowdownSolve {
    EquilibriumRocketBlowdownSolve { reservoir: EquilibriumBlowdownSolve {
        minimum_pressure_pa: 1000.0, relative_density_tolerance: 2e-11,
        mass_tolerance_kg: 1e-13, energy_tolerance_j: 1e-5, maximum_density_iterations: 80 },
        nozzle: EquilibriumRocketNozzleSolve { thermochemistry: EquilibriumEntropySolveV1 {
            equilibrium: EquilibriumSolveV1 { relative_balance_tolerance: 1e-12,
                dimensionless_potential_tolerance: 1e-11, ..EquilibriumSolveV1::default() },
            minimum_temperature_k: 200.0, maximum_temperature_k: 6000.0,
            entropy_tolerance_j_k: 1e-9, maximum_temperature_iterations: 100 },
            minimum_pressure_pa: 1000.0, sonic_tolerance: 1e-9, relative_mass_flow_tolerance: 1e-9,
            energy_tolerance_j_kg: 1e-6, maximum_pressure_iterations: 100 },
        maximum_outgoing_fraction: 0.5, maximum_mass_iterations: 32 }
}
fn state(reacting: bool, velocity: Vec3) -> EquilibriumGasRocketState {
    let models = if reacting { vec![model("H", 0.001, 2.5, 12000.0, 1.0), model("H2", 0.002, 3.5, 0.0, 0.0)] }
        else { vec![model("H2", 0.002, 3.5, 0.0, 0.0)] };
    let ids: Vec<_> = models.iter().map(|m| m.species().clone()).collect();
    let counts = if reacting { vec![1, 2] } else { vec![2] };
    let chemistry = Arc::new(IdealGasEquilibriumV1::new(
        ElementalMatrix::new(vec![ElementId::new("H").unwrap()], ids.clone(), vec![counts]).unwrap(),
        ChargeVector::new(ids.clone(), vec![0; ids.len()]).unwrap(), models).unwrap());
    let reservoir = EquilibriumGasReservoir::new(chemistry, vec![(id("H2"), Amount::new(1.0))],
        0.02, Temperature::new(2000.0), Pressure::new(100000.0), controls().nozzle.thermochemistry, || false).unwrap();
    EquilibriumGasRocketState::new(reservoir, 0.001, Vec3 {x:-1.0,y:0.2,z:5.0}, velocity, 7).unwrap()
}
fn velocity() -> Vec3 { Vec3 {x:2.0,y:-1.0,z:3.0} }
fn nozzle() -> RocketNozzle { RocketNozzle { throat_area_m2:1e-5, exit_area_m2:5e-5 } }
fn environment() -> GasRocketEnvironment {
    GasRocketEnvironment { inertial_frame:7, thrust_direction_world:Vec3 {x:0.6,y:0.0,z:0.8},
        ambient_pressure_pa:0.0, gravity_m_s2:Vec3 {x:0.0,y:0.0,z:-9.81},
        external_force_n:Vec3 {x:0.0001,y:-0.0002,z:0.0003} }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance*b.abs().max(1.0), "{a:.17e} != {b:.17e}");
}
fn vector_near(a: Vec3, b: Vec3, tolerance: f64) {
    near(a.x,b.x,tolerance); near(a.y,b.y,tolerance); near(a.z,b.z,tolerance);
}

#[test]
fn g1_native_translation_refines_against_continuous_perfect_gas_flight() {
    // Independently integrated x'=v, v'=(direction*F0*xmass^1.4+Fext)/(dry+m)+g,
    // m=m0*(1+0.2*(mdot0/m0)*t)^-5, with DOP853 at rtol 2e-13.
    let exact_position = Vec3 {x:21.272745776833631,y:-0.055114574581530566,z:34.478027347929448};
    let exact_velocity = Vec3 {x:168.03318716779714,y:-1.0412073551328935,z:221.95942235300726};
    let initial = state(false, velocity());
    let mut errors = Vec::new();
    for n in [4,8,16] {
        let mut current = initial.clone(); let mut exported = 0.0;
        for _ in 0..n {
            let (next, step) = current.advance(environment(), nozzle(), 0.25/f64::from(n), controls(), || false).unwrap().into_parts();
            let propulsive = step.propulsion.as_ref().unwrap();
            assert!(step.remaining_mass_residual_kg.abs() <= controls().reservoir.mass_tolerance_kg);
            assert!(step.consumed_mass_residual_kg.abs() <= controls().reservoir.mass_tolerance_kg);
            vector_near(step.mechanics.thrust_impulse_n_s,
                environment().thrust_direction_world.scale(propulsive.axial_impulse_n_s), 2e-14);
            assert_eq!(next.vehicle().dry_mass_kg, initial.vehicle().dry_mass_kg);
            exported += propulsive.gas.exhausted_enthalpy_j; current = next;
        }
        near(exported, initial.reservoir().internal_energy_j()-current.reservoir().internal_energy_j(), 2e-8);
        errors.push((current.vehicle().position_m.sub(exact_position).norm_squared().sqrt(),
            current.vehicle().velocity_m_s.sub(exact_velocity).norm_squared().sqrt()));
    }
    for pair in errors.windows(2) {
        assert!((3.4..4.6).contains(&(pair[0].0/pair[1].0)), "position {errors:?}");
        assert!((3.4..4.6).contains(&(pair[0].1/pair[1].1)), "velocity {errors:?}");
    }
    assert!(errors[2].0 < 0.007 && errors[2].1 < 0.007, "{errors:?}");
}

#[test]
fn g1_g3_reacting_flight_counts_pressure_thrust_once_and_is_galilean_invariant() {
    let initial = state(true, velocity()); let dt = 0.0625;
    let base = initial.advance(environment(), nozzle(), dt, controls(), || false).unwrap();
    let mut pressure_environment = environment(); pressure_environment.ambient_pressure_pa = 1000.0;
    let pressurized = initial.advance(pressure_environment, nozzle(), dt, controls(), || false).unwrap();
    assert_eq!(base.endpoint().reservoir(), pressurized.endpoint().reservoir());
    let base_step = base.step().propulsion.as_ref().unwrap();
    let pressure_step = pressurized.step().propulsion.as_ref().unwrap();
    let difference = dt*1000.0*nozzle().exit_area_m2;
    near(base_step.axial_impulse_n_s-pressure_step.axial_impulse_n_s, difference, 1e-10);
    let m0 = initial.vehicle().dry_mass_kg+initial.vehicle().propellant_mass_kg;
    let m1 = base.endpoint().vehicle().dry_mass_kg+base.endpoint().vehicle().propellant_mass_kg;
    let expected = environment().thrust_direction_world.scale(difference/base_step.gas.expelled_mass_kg*(m0/m1).ln());
    vector_near(base.endpoint().vehicle().velocity_m_s.sub(pressurized.endpoint().vehicle().velocity_m_s), expected, 1e-8);
    let boost = Vec3 {x:123.0,y:-67.0,z:41.0};
    let boosted = state(true, velocity().add(boost)).advance(environment(), nozzle(), dt, controls(), || false).unwrap();
    assert_eq!(base.endpoint().reservoir(), boosted.endpoint().reservoir());
    vector_near(boosted.endpoint().vehicle().velocity_m_s.sub(base.endpoint().vehicle().velocity_m_s), boost, 1e-12);
    vector_near(boosted.endpoint().vehicle().position_m.sub(base.endpoint().vehicle().position_m), boost.scale(dt), 1e-12);
    let checkpoint = base.endpoint().clone();
    assert_eq!(base.endpoint().advance(environment(), nozzle(), dt, controls(), || false).unwrap(),
        checkpoint.advance(environment(), nozzle(), dt, controls(), || false).unwrap());
}

#[test]
fn g1_g4_closed_valve_and_late_failures_preserve_the_single_inventory() {
    let initial = state(true, velocity()); let before = initial.clone();
    let mut env = environment(); env.ambient_pressure_pa = 1000000.0;
    let dt = 0.3;
    let coast = initial.coast(env, dt, || false).unwrap();
    assert_eq!(coast.endpoint().reservoir(), initial.reservoir());
    assert!(coast.step().propulsion.is_none());
    assert_eq!(coast.step().mechanics.burned_propellant_mass_kg, 0.0);
    assert_eq!(coast.step().mechanics.thrust_impulse_n_s, Vec3::ZERO);
    let wet_mass = initial.vehicle().dry_mass_kg+initial.vehicle().propellant_mass_kg;
    let acceleration = env.gravity_m_s2.add(env.external_force_n.scale(1.0/wet_mass));
    vector_near(coast.endpoint().vehicle().velocity_m_s, velocity().add(acceleration.scale(dt)), 1e-13);
    vector_near(coast.endpoint().vehicle().position_m,
        initial.vehicle().position_m.add(velocity().scale(dt)).add(acceleration.scale(0.5*dt*dt)), 1e-13);
    let reopened = coast.endpoint().advance(environment(), nozzle(), 0.0625, controls(), || false).unwrap();
    let ordinary = initial.advance(environment(), nozzle(), 0.0625, controls(), || false).unwrap();
    assert_eq!(reopened.endpoint().reservoir(), ordinary.endpoint().reservoir());
    env = environment(); env.inertial_frame += 1;
    assert!(initial.advance(env, nozzle(), 0.0625, controls(), || false).is_err());
    env = environment(); env.thrust_direction_world = Vec3::ZERO;
    assert!(initial.advance(env, nozzle(), 0.0625, controls(), || false).is_err());
    env = environment(); env.external_force_n.x = 1e308;
    assert!(initial.advance(env, nozzle(), 0.0625, controls(), || false).is_err());
    let mut complete = 0;
    initial.advance(environment(), nozzle(), 0.0625, controls(), || {complete += 1; false}).unwrap();
    for limit in [0,complete/2,complete-1] {
        let mut polls = 0;
        assert!(initial.advance(environment(), nozzle(), 0.0625, controls(), || {polls += 1; polls > limit}).is_err());
    }
    assert!(initial.coast(environment(), 0.3, || true).is_err());
    assert_eq!(initial, before);
}
