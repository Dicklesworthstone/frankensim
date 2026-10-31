//! Real shared owners with synthetic material data; no surrogate time integrator.
use super::*;
use fs_airflow::graph::thermal::medium::convection::finite_gas::natural::NaturalGasSurface;
use fs_airflow::graph::thermal::medium::convection::finite_gas::network::exchange::{
    self, ExchangeContact, GasHeatLaw, GasHeatSample,
};
use fs_couple::cavity_phs::piston::thermal::boundary::network::material::GasConductionGeometry;
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::material::state::{ChamberDrive, MaterialExchangeState};
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
use fs_material::thermochemical_gas::energy::GasEnergySolve;
use fs_material::thermochemical_gas::energy::parcel::transfer::GasMassTransferSolve;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::{GasNozzle, GasNozzleStepSolve,
    network::{GasNetworkLink, GasNetworkSolve}};
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
    [ExchangeContact {contact:base[0].contact,law:GasHeatLaw::NaturalConvection(
            NaturalGasSurface {area_m2:base[0].surface.area_m2,
                geometry:fs_convection::NaturalConvectionGeometry::VerticalPlate,
                characteristic_length:base[0].surface.characteristic_length,gravity_m_s2:9.81})},
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
                    GasHeatSample::ForcedConvection(sample) | GasHeatSample::NaturalConvection(sample)=> {
                        assert_eq!(sample.gas,entry.boundary.piston.midpoint);assert_eq!(sample.body,cell.phase);
                        let coefficient=sample.evaluation.as_ref().unwrap().coefficient.value.value();
                        let expected=dt*coefficient*sample.area_m2
                            *(sample.body.temperature_k()-sample.gas.state().temperature);
                        assert!((patch.heat_j-expected).abs()<1e-9);
                        assert!(patch.heat_j<0.0); // Actual hotter gas heats the finite cells.
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


fn live_drives(gas:&BTreeMap<u64,GasPiston>,dt:f64)->Vec<ChamberDrive> {
    chambers(gas,dt).into_iter().map(|c|ChamberDrive {
        chamber:c.id,heat_j:c.input.heat_j,exterior_pressure_pa:c.input.exterior_pressure_pa,
        external_force_n:c.input.external_force_n,piston_solve:c.piston_solve,exchange:c.exchange,
    }).collect()
}
fn transported_gases()->BTreeMap<u64,GasPiston> {
    let mut gas=gases(0.75);
    gas.insert(20,GasPiston::new(ClosedGasParcel::new(
        gas_source::model(0.25),gas_source::point(400.0,180000.0),0.001).unwrap(),
        1.0,0.01,0.0).unwrap());
    gas
}
fn flow_link(gas:&BTreeMap<u64,GasPiston>)->([GasNetworkLink;1],GasNetworkSolve) {
    let energy=chambers(gas,DT)[0].piston_solve.gas;
    ([GasNetworkLink {id:77,parcels:[20,10],
        nozzle:GasNozzle {area_m2:2e-6,discharge_coefficient:0.8},
        solve:GasNozzleStepSolve {nozzle:energy,transfer:GasMassTransferSolve {
            donor:energy,receiver:energy,mass_tolerance_kg:1e-14,energy_tolerance_j:1e-8},
            maximum_donor_fraction:0.05}}],
     GasNetworkSolve {maximum_parcels:2,maximum_links:1,maximum_withdrawn_fraction:0.05,
        mass_tolerance_kg:1e-13,energy_tolerance_j:1e-7})
}

#[test]
fn live_transport_changes_the_mixture_seen_by_mixed_heat_before_atomic_publish() {
    fixture(|binding,thermal,circuit|scoped(|cx,_| {
        let electrical=ElectricalContactNetwork {bulk:binding,contacts:&[],electrical_bases:&BASES};
        let assembly=MaterialExchangeNetwork {electrical:&electrical,thermal_contacts:&[],radiation:None};
        let prepared=assembly.prepare(cx,&thermal).unwrap();
        let mut state=MaterialExchangeState::from_parts(thermal,circuit,transported_gases());
        let drives=live_drives(state.chambers(),DT);
        let patches=mixed(1.0);
        let (links,flow_controls)=flow_link(state.chambers());
        let initial_fraction=state.chambers()[&10].gas().model()
            .thermodynamics().mole_composition().fractions()[0];

        let (preview,report)=exchange::flow_step_candidate(&state,&prepared,cx,DT,&[],
            &drives,&[],&patches,&links,flow_controls,controls(),&[0.0;2]).unwrap();
        assert!(report.flow.links[0].transfer.transferred_mass_kg>0.0);
        let preview_fraction=preview.chambers()[&10].gas().model()
            .thermodynamics().mole_composition().fractions()[0];
        assert_ne!(preview_fraction,initial_fraction);
        assert_eq!(state.time_s(),0.0);
        for chamber in &report.exchange.report.source.chambers {
            for surface in &chamber.surfaces {
                if chamber.id==10 {
                    match &surface.source.1 {
                        GasHeatSample::ForcedConvection(sample) | GasHeatSample::NaturalConvection(sample)=>
                            assert_eq!(&sample.gas,&chamber.boundary.piston.midpoint),
                        GasHeatSample::Conduction(_) => {}
                    }
                }
            }
        }

        let accepted=exchange::advance_with_flow(&mut state,&prepared,cx,DT,&[],
            &drives,&[],&patches,&links,flow_controls,controls(),&[0.0;2]).unwrap();
        assert_eq!(state.time_s(),DT);
        assert_eq!(state.chambers(),preview.chambers());
        assert_eq!(state.thermal().cells(),preview.thermal().cells());
        assert_eq!(state.circuit(),preview.circuit());
        assert_eq!(accepted.total_mass_residual_kg,report.total_mass_residual_kg);
        assert_eq!(accepted.total_energy_residual_j,report.total_energy_residual_j);
    }));
}
