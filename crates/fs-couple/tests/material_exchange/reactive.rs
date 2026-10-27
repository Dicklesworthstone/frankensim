//! Reaction extent is data, not an extra heat source or a Newton-iteration clock.
use super::*;
use fs_couple::cavity_phs::piston::{GasPiston, reactive::{GasReaction, GasReactionSolve}};
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::material::reactive::*;
use fs_material::thermochemical_gas::{ThermochemicalGasModel, energy::parcel::ClosedGasParcel};
use fs_qty::QuantitySpec;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_thermochem::*;
use fs_thermochem::equilibrium::*;
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;

pub(super) fn reaction(rate:f64) -> GasReaction {
    let (a,ta)=gas_fixture::component("A",0.028,0.0);
    let (b,tb)=gas_fixture::component("B",0.028,-50000.0);
    let species=vec![a.species().clone(),b.species().clone()];
    let matrix=StoichiometricMatrix::new(species.clone(),vec![ReactionId::new("isomerization").unwrap()],vec![vec![-1],vec![1]]).unwrap();
    let elements=ElementalMatrix::new(vec![ElementId::new("X").unwrap()],species.clone(),vec![vec![1,1]]).unwrap();
    let charge=ChargeVector::new(species,vec![0,0]).unwrap();
    let certificate=verify_conservation(&elements,&matrix,&charge).unwrap();
    let equilibrium=IdealGasReactionEquilibriumV1::new(matrix,certificate,ReactionId::new("isomerization").unwrap(),vec![a.clone(),b.clone()]).unwrap();
    let law=DeclaredStoichiometricMassActionV1::new(ThermodynamicReverseRateClosureV1::new(equilibrium,
        ReactionProgressRateScaleV1::new(ReactionProgressRateQuantityV1::new(rate)).unwrap()),
        MassActionKineticOrderConventionV1::CallerDeclaredStoichiometricCoefficients);
    let catalog=ThermochemicalGasModel::new(FrozenIdealGasMixtureModelV1::new(vec![a,b],
        Composition::new(CompositionBasis::MoleFraction,vec![0.5,0.5]).unwrap()).unwrap(),vec![ta,tb],None).unwrap();
    GasReaction::new(catalog,law).unwrap()
}
pub(super) fn piston(law:&GasReaction) -> GasPiston {
    let model=law.catalog().with_composition(Composition::new(CompositionBasis::MoleFraction,vec![1.0,0.0]).unwrap()).unwrap();
    let q=|k|QuantitySpec::semantic(SemanticType::new(k,ValueForm::Static));
    let point=QueryPoint::new().with_quantity("temperature",q(QuantityKind::AbsoluteTemperature),400.0).unwrap()
        .with_quantity("pressure",q(QuantityKind::Pressure),100000.0).unwrap();
    GasPiston::new(ClosedGasParcel::new(model,point,0.001).unwrap(),1.0,0.01,0.0).unwrap()
}
pub(super) fn chemical() -> GasReactionSolve {
    GasReactionSolve {gas:gas_fixture::chamber(1,&gas_fixture::pistons()[0],0.00005).piston_solve.gas,
        maximum_consumed_fraction:0.2,species_tolerance_mol:1e-12,entropy_tolerance_j_k:1e-10}
}

#[test]
fn isolated_reacting_chamber_reduces_to_existing_reactive_piston_and_advances_once() {
    fixture([0.2,0.2],150.0,|assembly,thermal,circuit|with_cx(|cx,_| {
        let prepared=assembly.prepare(cx,&thermal).unwrap();
        let law=reaction(50.0);let gas=piston(&law);let ch=gas_fixture::chamber(10,&gas,0.00005);
        let native=gas.reactive_candidate(&law,ch.input,chemical(),ch.piston_solve,||false).unwrap();
        let result=prepared.reactive_step_candidate(cx,&thermal,0.00005,&circuit,&[],&[ch],
            &[ChamberReaction {chamber:10,reaction:&law,controls:chemical()}],&[],controls(),&[0.0;2]).unwrap();
        assert_eq!(&result.chambers()[0].piston,native.piston());
        assert_eq!(result.reactions()[0].reaction,native.step().reaction);
        near(result.reactions()[0].gas_energy_change_j,native.step().gas_energy_change_j);
        assert!(result.reactions()[0].reaction.extent_mol>0.0);
        assert_eq!(gas.gas().model().thermodynamics().components().len(),1);
        assert_eq!(result.chambers()[0].piston.gas().model().thermodynamics().components().len(),2);
        assert!(result.step().total_energy_residual_j.abs()<=controls().total_energy_tolerance_j);
        let ordinary=prepared.step_candidate(cx,&thermal,0.00005,&circuit,&[],&[ch],&[],controls(),&[0.0;2]).unwrap();
        let inert=prepared.reactive_step_candidate(cx,&thermal,0.00005,&circuit,&[],&[ch],&[],&[],controls(),&[0.0;2]).unwrap();
        assert_eq!(inert.thermal().cells(),ordinary.thermal().cells());
        assert_eq!(inert.circuit(),ordinary.circuit());
        assert_eq!(&inert.chambers()[0].piston,ordinary.chambers()[0].boundary.piston());
        assert!(inert.reactions().is_empty());
    }));
}

#[test]
fn reaction_rate_data_changes_joint_material_wall_heat_motion_and_acoustic_state() {
    let run=|rate|fixture([0.8,0.8],150.0,|assembly,mut thermal,mut circuit|with_cx(|cx,_| {
        let prepared=assembly.prepare(cx,&thermal).unwrap();let law=reaction(rate);
        let mut gases=[piston(&law),gas_fixture::pistons()[1].clone()];let initial=gases.clone();
        let energy0=thermal.energy_j()+50.0;let mut work=0.0;
        for _ in 0..20 {
            let target=thermal.time_s()+0.00005;let dt=target-thermal.time_s();
            let chambers=[gas_fixture::chamber(20,&gases[1],dt),gas_fixture::chamber(10,&gases[0],dt)];
            let expected=law.parcel_candidate(gases[0].gas(),dt,chemical(),||false).unwrap();
            let result=prepared.reactive_step_candidate(cx,&thermal,target,&circuit,&[],&chambers,
                &[ChamberReaction {chamber:10,reaction:&law,controls:chemical()}],&gas_fixture::patches(),controls(),&[0.0;2]).unwrap();
            assert_eq!(result.reactions().len(),1);
            // Exact equality with ONE independent kinetic stage, regardless of
            // the number of electrical, radiation and gas Newton evaluations.
            assert_eq!(&result.reactions()[0].reaction,expected.step());
            near(result.step().external_network_heat_j,0.0);
            assert!(result.report().electrical.contacts[0].heat_j>0.0);
            let (t,c,r,_)=result.into_parts();
            for (i,chamber) in r.source.chambers.into_iter().enumerate() {
                work+=chamber.boundary.piston.exterior_pressure_work_j;
                assert_eq!(chamber.id,if i==0 {10} else {20});
                gases[i]=chamber.piston;
            }
            let gas_energy:f64=gases.iter().zip(&initial).map(|(g,o)|g.gas().mass_kg()
                *(g.gas().state().specific_internal_energy_j_kg()-o.gas().state().specific_internal_energy_j_kg())
                +0.5*g.mass_kg()*g.velocity_m_s().powi(2)).sum();
            assert!((t.energy_j()+r.electrical.stored_after_j+gas_energy+work-energy0).abs()<4e-5);
            thermal=t;circuit=c;
        }
        (thermal.temperature_k(CellId(30)).unwrap(),thermal.temperature_k(CellId(10)).unwrap(),circuit.charges_c()[0],gases[0].clone())
    }));
    let low=run(50.0);let high=run(200.0);
    assert!((high.0-low.0).abs()>1e-6);
    assert!((high.1-low.1).abs()>1e-7);
    assert!((high.2-low.2).abs()>1e-10);
    assert!(high.3.position_m()>low.3.position_m());
    assert!(high.3.gas().state().state().sound_speed>low.3.gas().state().state().sound_speed);
}

#[test]
fn reaction_mapping_consumption_and_late_exchange_failure_preserve_all_originals() {
    fixture([0.2,0.2],150.0,|assembly,thermal,circuit|with_cx(|cx,gate| {
        let prepared=assembly.prepare(cx,&thermal).unwrap();let law=reaction(50.0);let gas=piston(&law);
        let before=gas.clone();let cells=thermal.cells().clone();let electrical=circuit.clone();
        let ch=gas_fixture::chamber(10,&gas,0.00005);
        let binding=ChamberReaction {chamber:10,reaction:&law,controls:chemical()};
        let attempt=|chambers:&[_],reactions:&[_],settings|prepared.reactive_step_candidate(cx,&thermal,0.00005,
            &circuit,&[],chambers,reactions,&[],settings,&[0.0;2]);
        assert!(attempt(&[ch],&[binding,binding],controls()).is_err());
        assert!(attempt(&[ch],&[ChamberReaction {chamber:99,..binding}],controls()).is_err());
        assert!(attempt(&[ch,gas_fixture::chamber(20,&gas,0.00005)],&[binding],controls()).is_err());
        let rapid=reaction(1e9);
        assert!(matches!(attempt(&[ch],&[ChamberReaction {reaction:&rapid,..binding}],controls()),
            Err(ReactiveExchangeError::Reaction {error:fs_couple::cavity_phs::piston::reactive::GasReactionError::ConsumptionLimit {..},..})));
        let mut limited=controls();limited.network.maximum_circuit_trials=1;
        assert!(attempt(&[ch],&[binding],limited).is_err());
        // A discarded proposal cannot change the consumed amount on retry.
        let expected=law.parcel_candidate(gas.gas(),0.00005,chemical(),||false).unwrap();
        let result=attempt(&[ch],&[binding],controls()).unwrap();
        assert_eq!(&result.reactions()[0].reaction,expected.step());
        gate.request();assert!(matches!(attempt(&[ch],&[binding],controls()),Err(ReactiveExchangeError::Cancelled)));
        assert_eq!(gas,before);assert_eq!(thermal.cells(),&cells);assert_eq!(circuit,electrical);assert_eq!(thermal.time_s(),0.0);
    }));
}

#[path = "state.rs"]
mod state_tests;

#[path = "gas_law.rs"]
mod gas_law_tests;

#[path = "subcycling.rs"]
mod subcycling_tests;
