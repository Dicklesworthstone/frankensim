//! Synthetic source data exercise the real shared gas, enthalpy and layer owners.
use super::*;
use crate::cavity_phs::piston::{GasPiston, tests::{gas, controls as piston_controls}};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use std::convert::Infallible;

fn input() -> DoublePistonInput {
    DoublePistonInput { duration_s: 0.0005, heat_j: [0.0; 2], external_force_n: 0.0 }
}
fn controls() -> DoublePistonSolve {
    DoublePistonSolve { piston: piston_controls(), opposing_gas: piston_controls().gas,
        volume_tolerance_m3: 1e-18 }
}
fn exchange() -> DoubleHeatExchangeSolve {
    DoubleHeatExchangeSolve { faces: [HeatExchangeSolve { heat_bracket_j: [-1.0, 1.0],
        source_tolerance_j: 1e-9, energy_tolerance_j: 1e-7, maximum_iterations: 128 }; 2],
        maximum_joint_trials: 4096 }
}
fn pair() -> DoubleActingGasPiston {
    DoubleActingGasPiston::new(GasPiston::new(gas(0.75, 0.001), 1.0, 0.01, 0.0).unwrap(),
        gas(0.25, 0.001), 0.01).unwrap()
}
fn curve(cp: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(fs_blake3::hash_bytes(format!("synthetic cp={cp}").as_bytes()),
        SolidLiquidPhase::Solid, [200.0, 800.0].into_iter().map(|t| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: cp*(t-200.0), temperature_k: t,
            liquid_mass_fraction: 0.0, bulk_density_kg_m3: 4000.0,
        }).collect()).unwrap()
}
fn body(curve: &EquilibriumEnthalpyPhaseCurve, h: f64) -> UniformEnthalpyStepInput<'_> {
    UniformEnthalpyStepInput { curve, initial: curve.state_at_specific_enthalpy(h).unwrap(),
        mass_kg: 0.001, volume_m3: 0.001/4000.0, surface_area_m2: 0.02,
        internal_heat_j: 0.0, duration_s: input().duration_s }
}
fn links(g: [f64; 2]) -> impl FnMut(&[ResolvedThermochemicalGasState; 2], EquilibriumPhaseState)
    -> Result<[ConductanceSample<EquilibriumPhaseState>; 2], Infallible>
{
    move |_, state| Ok(g.map(|conductance_w_k| ConductanceSample { conductance_w_k, source: state }))
}
fn check<S>(candidate: &DoubleThermalPistonCandidate<S>) {
    let r = candidate.step();
    for i in 0..2 { assert!(r.source_residual_j[i].abs() <= exchange().faces[i].source_tolerance_j); }
    assert!(r.total_energy_residual_j.abs() <= 1e-7);
    assert!(r.thermal_energy_residual_j.abs() <= 1e-7);
    assert_eq!(candidate.thermal().external_heat_j, -r.exchanged_heat_j.iter().sum::<f64>());
    assert!(r.joint_trials <= exchange().maximum_joint_trials);
    assert!(r.material_trials <= 2*(controls().piston.maximum_iterations+3)*r.joint_trials);
}

#[test]
fn equilibrium_and_two_insulated_faces_keep_original_dynamics() {
    let old = pair(); let chart = curve(500.0);
    let equilibrium = body(&chart, 100000.0);
    let c = old.thermal_candidate(input(), equilibrium, controls(), exchange(), links([1.0, 0.3]), || false).unwrap();
    assert_eq!(c.piston(), &old);
    assert_eq!(c.thermal().state, equilibrium.initial);
    assert_eq!(c.step().exchanged_heat_j, [0.0; 2]);
    assert_eq!(c.step().joint_trials, 1);
    check(&c);
    let drive = DoublePistonInput { heat_j: [0.02, 0.01], ..input() };
    let hot = body(&chart, 150000.0);
    let isolated = old.thermal_candidate(drive, hot, controls(), exchange(), links([0.0; 2]), || false).unwrap();
    let mut independent = old.clone();
    let report = independent.step(drive, controls(), || false).unwrap();
    assert_eq!(isolated.piston(), &independent);
    assert_eq!(isolated.step().piston, report);
    assert_eq!(isolated.thermal().state, hot.initial);
    check(&isolated);
}

#[test]
fn both_faces_debit_one_body_and_use_the_same_enthalpy_midpoint() {
    let old = pair(); let chart = curve(500.0); let thermal = body(&chart, 150000.0);
    let c = old.thermal_candidate(input(), thermal, controls(), exchange(), links([1.0, 0.3]), || false).unwrap();
    let r = c.step();
    assert!(r.exchanged_heat_j.iter().all(|q| *q > 0.0));
    assert!(c.piston().primary().velocity_m_s() > 0.0);
    assert!(c.piston().opposing().state().state().temperature > 400.0);
    assert_eq!(c.thermal().report, [r.thermal_midpoint; 2]);
    assert!((thermal.mass_kg*(c.thermal().state.specific_enthalpy_j_kg()
        -thermal.initial.specific_enthalpy_j_kg())+r.exchanged_heat_j.iter().sum::<f64>()).abs() < 1e-9);
    check(&c);
}

#[test]
fn phase_storage_can_take_heat_from_both_gases_without_temperature_rise() {
    let knot = |h,t,x| EnthalpyPhaseKnot { specific_enthalpy_j_kg: h, temperature_k: t,
        liquid_mass_fraction: x, bulk_density_kg_m3: 4000.0 };
    let chart = EquilibriumEnthalpyPhaseCurve::try_new(fs_blake3::hash_bytes(b"synthetic shared latent store"),
        vec![knot(0.0,200.0,0.0), knot(100000.0,350.0,0.0),
            knot(160000.0,350.0,1.0), knot(400000.0,800.0,1.0)]).unwrap();
    let thermal = body(&chart, 130000.0);
    let c = pair().thermal_candidate(input(), thermal, controls(), exchange(), links([1.0,0.3]), || false).unwrap();
    assert!(c.step().exchanged_heat_j.iter().all(|q| *q < 0.0));
    assert_eq!(c.thermal().state.temperature_k(),350.0);
    assert_eq!(c.step().thermal_midpoint.temperature_k(),350.0);
    assert!(c.thermal().state.liquid_mass_fraction() > thermal.initial.liquid_mass_fraction());
    check(&c);
}

#[test]
fn repeated_shared_heating_closes_all_storage_and_prescribed_work() {
    let initial = pair(); let mut piston = initial.clone(); let chart = curve(500.0);
    let mut thermal = body(&chart,150000.0); let h0 = thermal.initial.specific_enthalpy_j_kg();
    let mut supplied = 0.0;
    for _ in 0..20 {
        let drive = DoublePistonInput { heat_j: [0.001,0.002], external_force_n: 0.5, ..input() };
        thermal.internal_heat_j = 0.003;
        let c = piston.thermal_candidate(drive,thermal,controls(),exchange(),links([1.0,0.3]),||false).unwrap();
        check(&c);
        supplied += 0.006+c.step().piston.external_force_work_j;
        let (next, accepted, _) = c.into_parts(); piston = next; thermal.initial = accepted.state;
    }
    let gas_storage: f64 = [(initial.primary().gas(),piston.primary().gas()), (initial.opposing(),piston.opposing())]
        .into_iter().map(|(a,b)| a.mass_kg()*(b.state().specific_internal_energy_j_kg()
            -a.state().specific_internal_energy_j_kg())).sum();
    let kinetic = 0.5*piston.primary().mass_kg()*piston.primary().velocity_m_s().powi(2);
    let stored = thermal.mass_kg*(thermal.initial.specific_enthalpy_j_kg()-h0);
    assert!((gas_storage+kinetic+stored-supplied).abs() < 2e-6);
}

#[test]
fn trial_budget_source_failure_and_final_cancellation_return_no_partial_system() {
    let old = pair(); let saved = old.clone(); let chart = curve(500.0); let thermal = body(&chart,150000.0);
    assert!(old.thermal_candidate(input(),thermal,controls(),DoubleHeatExchangeSolve {
        maximum_joint_trials:1,..exchange() },links([1.0,0.3]),||false).is_err());
    let failure = |_: &[ResolvedThermochemicalGasState;2],_:EquilibriumPhaseState|
        -> Result<[ConductanceSample<()>;2],&str> { Err("opposing interface source missing") };
    assert!(matches!(old.thermal_candidate(input(),thermal,controls(),exchange(),failure,||false),
        Err(ThermalPistonError::Interface("opposing interface source missing"))));
    for g in [[1.0,-1.0],[f64::NAN,0.3]] {
        assert!(old.thermal_candidate(input(),thermal,controls(),exchange(),links(g),||false).is_err());
    }
    assert!(old.thermal_candidate(input(),UniformEnthalpyStepInput { duration_s:0.1,..thermal },
        controls(),exchange(),links([1.0,0.3]),||false).is_err());
    let mut polls = 0;
    old.thermal_candidate(input(),thermal,controls(),exchange(),links([1.0,0.3]),||{polls+=1;false}).unwrap();
    for at in [1,polls/2,polls] {
        let mut n = 0;
        assert!(matches!(old.thermal_candidate(input(),thermal,controls(),exchange(),links([1.0,0.3]),
            ||{n+=1;n==at}),Err(ThermalPistonError::Cancelled)));
    }
    assert_eq!(old,saved);
}

mod source_layers {
    use super::*;
    use crate::cavity_phs::piston::thermal::layers::MaterialConductionLayers;
    use fs_evidence::ValidityDomain;
    use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
        PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel};
    use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve};
    use fs_material::state_point::{MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_DIMS};
    use fs_qty::{Density, Temperature, QuantitySpec};
    use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

    fn point(t: f64) -> QueryPoint {
        QueryPoint::new().with_quantity("T",QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,ValueForm::Static)),t).unwrap().with("lot",7.0).unwrap()
    }
    fn material(k: f64) -> (MaterialCard,PropertyKey,ResolvedHeatCapacityCurve) {
        let conductivity = PropertyKey::new("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS);
        let cp = PropertyKey::new("specific_heat_capacity",SPECIFIC_HEAT_CAPACITY_DIMS);
        let density = PropertyKey::new("density",Density::DIMS);
        let mut claims = ClaimSet::new();
        for (key,a,b) in [(cp.clone(),500.0,500.0),(density.clone(),4000.0,4000.0),(conductivity.clone(),k,2.0*k)] {
            let dims = key.quantity().dims();
            claims.insert_claim(PropertyClaim { key, value:PropertyValue::Curve {
                abscissa:"T".into(),abscissa_dims:Temperature::DIMS,knots:vec![(200.0,a),(800.0,b)],dims },
                validity:ValidityDomain::unconstrained().with_quantity("T",QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::AbsoluteTemperature,ValueForm::Static)),200.0,800.0).with("lot",7.0,7.0),
                uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::LinearInside,
                observations:vec![],provenance:Provenance { source:"synthetic coupled body data".into(),
                    license:"CC0-1.0".into(),artifact:None } }).unwrap();
        }
        let card = MaterialCard::assemble(MaterialStateId { chemistry:"synthetic".into(),phase:"solid".into(),
            process:"authored".into(),revision:0 },claims,vec![]).unwrap();
        let chart = ResolvedHeatCapacityCurve::try_new(&card,HeatCapacitySourceQuery { heat_capacity:&cp,density:&density,
            temperature_axis:"T",lower:&point(200.0),upper:&point(800.0),selection:SelectionPolicy::SingleClaimOnly },
            SolidLiquidPhase::Solid,0.0,1e-5).unwrap();
        (card,conductivity,chart)
    }
    fn layer<'a>(card: &'a MaterialCard, key: &'a PropertyKey, context: &'a QueryPoint, area: f64) -> MaterialConductionLayers<'a> {
        MaterialConductionLayers { material:card,conductivity:key,context,temperature_axis:"T",liquid_fraction_axis:None,
            selection:MaterialPropertySelection::SingleClaimOnly,area_m2:area,gas_distance_m:0.001,body_distance_m:0.02 }
    }
    #[test]
    fn conductivity_only_source_swap_changes_both_gases_and_common_motion() {
        let run = |k| {
            let (card,key,chart) = material(k); let context = point(300.0);
            let thermal = body(chart.curve(),150000.0);
            let c = pair().thermal_layers_candidate(input(),thermal,controls(),exchange(),
                [layer(&card,&key,&context,0.02),layer(&card,&key,&context,0.01)],||false).unwrap();
            for i in 0..2 {
                let original = layer(&card,&key,&context,[0.02,0.01][i])
                    .resolve(&c.step().piston.midpoint[i],c.step().thermal_midpoint).unwrap();
                assert_eq!(c.thermal().report[i],original.source);
                assert_eq!(c.step().conductance_w_k[i],original.conductance_w_k);
            }
            assert_eq!(c.thermal().state.material_card_identity(),card.content_hash());
            check(&c);
            c
        };
        let slow = run(0.2); let fast = run(2.0);
        for i in 0..2 { assert!(fast.step().exchanged_heat_j[i] > slow.step().exchanged_heat_j[i]); }
        assert_ne!(fast.piston().primary().velocity_m_s(),slow.piston().primary().velocity_m_s());
        assert_ne!(fast.piston().opposing().state().state().sound_speed,slow.piston().opposing().state().state().sound_speed);
    }
    #[test]
    fn mismatched_opposing_body_source_is_rejected_even_when_insulated() {
        let (card,key,chart) = material(0.2); let (wrong,_,_) = material(2.0); let context = point(300.0);
        let old = pair(); let saved = old.clone();
        let result = old.thermal_layers_candidate(input(),body(chart.curve(),150000.0),controls(),exchange(),
            [layer(&card,&key,&context,0.02),layer(&wrong,&key,&context,0.0)],||false);
        assert!(matches!(result,Err(ThermalPistonError::Interface(_))));
        assert_eq!(old,saved);
    }
}
