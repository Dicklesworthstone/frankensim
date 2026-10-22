//! Real shared owners with synthetic material data; no surrogate time integrator.
use super::*;
use fs_airflow::graph::thermal::medium::convection::finite_gas::network::exchange::{
    self, ExchangeContact, GasHeatLaw, GasHeatSample,
};
use fs_couple::cavity_phs::piston::thermal::boundary::network::material::GasConductionGeometry;
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::material::{
    exchange::MaterialExchangeNetwork, reactive::{ChamberReaction, ReactiveExchangeError},
};
use fs_couple::circuit_enthalpy::contacts::material::ElectricalContactNetwork;
use fs_couple::cavity_phs::piston::reactive::{GasReaction,GasReactionSolve};
use fs_exec::{Cx,CancelGate};
use fs_thermochem::*;
use fs_thermochem::equilibrium::*;
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_material::thermochemical_gas::ThermochemicalGasModel;
#[path = "../../../fs-couple/tests/material_exchange/gas.rs"]
#[allow(dead_code)]
mod chemical_source;

fn scoped<R>(run:impl FnOnce(&Cx<'_>,&CancelGate)->R)->R {
    let gate=CancelGate::new();let pool=fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx=Cx::new(&gate,arena,fs_exec::StreamKey {seed:71,kernel_id:3,tile:0,iteration:0},
            fs_exec::Budget::INFINITE,fs_exec::ExecMode::Deterministic);
        run(&cx,&gate)
    })
}
fn mixed(speed:f64)->[ExchangeContact;3] {
    let base=contacts(speed);
    [ExchangeContact {contact:base[0].contact,law:GasHeatLaw::ForcedConvection(base[0].surface)},
        ExchangeContact {contact:base[1].contact,law:GasHeatLaw::Conduction(GasConductionGeometry {
            area_m2:base[1].surface.area_m2,gas_distance_m:0.01,cell_distance_m:0.02})},
        ExchangeContact {contact:base[2].contact,law:GasHeatLaw::ForcedConvection(base[2].surface)}]
}
fn chemical(rate:f64)->GasReaction {
    let (a,ta)=chemical_source::component("A",0.028,0.0);
    let (b,tb)=chemical_source::component("B",0.028,-5000.0);
    let species=vec![a.species().clone(),b.species().clone()];
    let id=ReactionId::new("synthetic-convective-isomers").unwrap();
    let matrix=StoichiometricMatrix::new(species.clone(),vec![id.clone()],vec![vec![-1],vec![1]]).unwrap();
    let elements=ElementalMatrix::new(vec![ElementId::new("X").unwrap()],species.clone(),vec![vec![1,1]]).unwrap();
    let conserved=verify_conservation(&elements,&matrix,&ChargeVector::new(species,vec![0,0]).unwrap()).unwrap();
    let equilibrium=IdealGasReactionEquilibriumV1::new(matrix,conserved,id,vec![a.clone(),b.clone()]).unwrap();
    let kinetics=DeclaredStoichiometricMassActionV1::new(ThermodynamicReverseRateClosureV1::new(equilibrium,
        ReactionProgressRateScaleV1::new(ReactionProgressRateQuantityV1::new(rate)).unwrap()),
        MassActionKineticOrderConventionV1::CallerDeclaredStoichiometricCoefficients);
    let catalog=ThermochemicalGasModel::new(FrozenIdealGasMixtureModelV1::new(vec![a,b],
        Composition::new(CompositionBasis::MoleFraction,vec![0.5,0.5]).unwrap()).unwrap(),vec![ta,tb],None).unwrap();
    GasReaction::new(catalog,kinetics).unwrap()
}
fn reacting_gases(law:&GasReaction)->BTreeMap<u64,GasPiston> {
    let model=law.catalog().with_composition(Composition::new(CompositionBasis::MoleFraction,vec![1.0,0.0]).unwrap()).unwrap();
    let mut gases=gases(0.75);
    gases.insert(10,GasPiston::new(ClosedGasParcel::new(model,gas_source::point(400.0,100000.0),0.001).unwrap(),1.0,0.01,0.0).unwrap());
    gases
}
fn chemistry(gases:&BTreeMap<u64,GasPiston>)->GasReactionSolve {
    GasReactionSolve {gas:chambers(gases,DT)[0].piston_solve.gas,
        maximum_consumed_fraction:0.2,species_tolerance_mol:1e-12,entropy_tolerance_j_k:1e-10}
}

#[test]
fn complete_assembly_reduces_exactly_to_original_convection_without_extra_physics() {
    fixture(|binding,thermal,circuit|scoped(|cx,_| {
        let electrical=ElectricalContactNetwork {bulk:binding,contacts:&[],electrical_bases:&BASES};
        let assembly=MaterialExchangeNetwork {electrical:&electrical,thermal_contacts:&[],radiation:None};
        let prepared=assembly.prepare(cx,&thermal).unwrap();let gas=gases(0.75);let drives=chambers(&gas,DT);
        let patches=contacts(1.0);let mixed=patches.map(|p|ExchangeContact {contact:p.contact,law:GasHeatLaw::ForcedConvection(p.surface)});
        let original=network::step_candidate(binding,&thermal,DT,&circuit,&[],&BASES,&drives,&patches,controls(),&[0.0;2],||false).unwrap();
        let updated=exchange::step_candidate(&prepared,cx,&thermal,DT,&circuit,&[],&drives,&[],&mixed,controls(),&[0.0;2]).unwrap();
        assert_eq!(original.thermal().cells(),updated.thermal().cells());
        assert_eq!(original.circuit(),updated.circuit());assert_eq!(original.step(),updated.step());
        for (a,b) in original.chambers().iter().zip(updated.chambers()) {
            assert_eq!(a.boundary.piston(),&b.piston);assert_eq!(a.boundary.step(),&b.boundary);
            for (old,new) in a.boundary.source().iter().zip(&b.surfaces) {
                let GasHeatSample::ForcedConvection(sample)=&new.source.1 else {panic!("lost convection source")};
                assert_eq!(sample.gas,old.source.1.gas);assert_eq!(sample.body,old.source.1.body);
                assert_eq!(new.heat_j,old.heat_j);
            }
        }
    }));
}

#[test]
fn mixed_laws_use_reacted_midpoints_and_material_changes_reach_charge_heat_and_motion() {
    let run=|rate|fixture(|binding,mut thermal,mut circuit|scoped(|cx,_| {
        let electrical=ElectricalContactNetwork {bulk:binding,contacts:&[],electrical_bases:&BASES};
        let assembly=MaterialExchangeNetwork {electrical:&electrical,thermal_contacts:&[],radiation:None};
        let prepared=assembly.prepare(cx,&thermal).unwrap();let law=chemical(rate);let mut gas=reacting_gases(&law);
        let initial=storage(&thermal,&circuit,&gas);let mut outside=0.0;let patches=mixed(1.0);
        for _ in 0..12 {
            let target=thermal.time_s()+DT;let dt=target-thermal.time_s();let solve=chemistry(&gas);
            let expected=law.parcel_candidate(gas[&10].gas(),dt,solve,||false).unwrap();
            let result=exchange::step_candidate(&prepared,cx,&thermal,target,&circuit,&[],&chambers(&gas,dt),
                &[ChamberReaction {chamber:10,reaction:&law,controls:solve}],&patches,controls(),&[0.0;2]).unwrap();
            assert_eq!(&result.reactions()[0].reaction,expected.step());
            assert!(result.step().total_energy_residual_j.abs()<=controls().total_energy_tolerance_j);
            for entry in result.chambers() {for patch in &entry.surfaces {
                let cell=result.report().source.material.material.bulk.cells.iter().find(|c|c.cell==patch.source.0.cell).unwrap();
                match &patch.source.1 {
                    GasHeatSample::ForcedConvection(sample)=> {
                        assert_eq!(sample.gas,entry.boundary.piston.midpoint);assert_eq!(sample.body,cell.phase);
                        if entry.id==10 {assert_eq!(sample.gas.components().len(),2);}
                    }
                    GasHeatSample::Conduction(material)=>assert_eq!(Some(material),cell.conductivity.as_ref()),
                }
            }}
            let (t,c,report,_)=result.into_parts();
            gas=report.source.chambers.into_iter().map(|entry| {
                outside+=entry.boundary.piston.exterior_pressure_work_j;(entry.id,entry.piston)
            }).collect();thermal=t;circuit=c;
            assert!((storage(&thermal,&circuit,&gas)+outside-initial).abs()<5e-5);
        }
        (thermal.temperature_k(CellId(30)).unwrap(),circuit.charges_c()[0],gas[&10].clone())
    }));
    let low=run(50.0);let high=run(200.0);
    assert_ne!(low.0,high.0);assert_ne!(low.1,high.1);
    assert_ne!(low.2.velocity_m_s(),high.2.velocity_m_s());
    assert!(high.2.gas().state().state().sound_speed>low.2.gas().state().state().sound_speed);
}

#[test]
fn conduction_and_convection_share_area_budget_order_and_failure_atomicity() {
    fixture(|binding,thermal,circuit|scoped(|cx,gate| {
        let electrical=ElectricalContactNetwork {bulk:binding,contacts:&[],electrical_bases:&BASES};
        let assembly=MaterialExchangeNetwork {electrical:&electrical,thermal_contacts:&[],radiation:None};
        let prepared=assembly.prepare(cx,&thermal).unwrap();let gas=gases(0.75);let saved=gas.clone();
        let cells=thermal.cells().clone();let electrical_before=circuit.clone();let patches=mixed(1.0);
        let a=exchange::step_candidate(&prepared,cx,&thermal,DT,&circuit,&[],&chambers(&gas,DT),&[],&patches,controls(),&[0.0;2]).unwrap();
        let reversed=[patches[2],patches[1],patches[0]];
        let b=exchange::step_candidate(&prepared,cx,&thermal,DT,&circuit,&[],&chambers(&gas,DT),&[],&reversed,controls(),&[0.0;2]).unwrap();
        assert_eq!(a.thermal().cells(),b.thermal().cells());assert_eq!(a.circuit(),b.circuit());
        for (a,b) in a.chambers().iter().zip(b.chambers()) {assert_eq!(a.piston,b.piston);}
        let too_large=[ExchangeContact {contact:patches[0].contact,law:GasHeatLaw::ForcedConvection(
            ForcedGasSurface {area_m2:3.0,..contacts(1.0)[0].surface})},
            ExchangeContact {contact:patches[2].contact,law:GasHeatLaw::Conduction(GasConductionGeometry {
                area_m2:3.0,gas_distance_m:0.01,cell_distance_m:0.02})}];
        assert!(exchange::step_candidate(&prepared,cx,&thermal,DT,&circuit,&[],&chambers(&gas,DT),&[],&too_large,controls(),&[0.0;2]).is_err());
        let mut limited=controls();limited.network.maximum_circuit_trials=1;
        assert!(exchange::step_candidate(&prepared,cx,&thermal,DT,&circuit,&[],&chambers(&gas,DT),&[],&patches,limited,&[0.0;2]).is_err());
        gate.request();assert!(matches!(exchange::step_candidate(&prepared,cx,&thermal,DT,&circuit,&[],
            &chambers(&gas,DT),&[],&patches,controls(),&[0.0;2]),Err(ReactiveExchangeError::Cancelled)));
        assert_eq!(gas,saved);assert_eq!(thermal.cells(),&cells);assert_eq!(circuit,electrical_before);
    }));
}
