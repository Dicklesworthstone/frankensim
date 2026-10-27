//! Joint ownership preserves the original solvers and their complete histories.
use super::*;
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::material::state::*;
use std::collections::BTreeMap;
use fs_ports::PortedModule;

fn loads(id:u64, gas:&GasPiston) -> ChamberDrive {
    let c=gas_fixture::chamber(id,gas,0.00005);
    ChamberDrive {chamber:id,heat_j:c.input.heat_j,exterior_pressure_pa:c.input.exterior_pressure_pa,
        external_force_n:c.input.external_force_n,piston_solve:c.piston_solve,exchange:c.exchange}
}
fn gas_map(law:&GasReaction) -> BTreeMap<u64,GasPiston> {
    [(20,gas_fixture::pistons()[1].clone()),(10,piston(law))].into_iter().collect()
}

#[test]
fn joint_advance_moves_the_exact_original_candidate_endpoints_and_reports() {
    fixture([0.8,0.8],150.0,|assembly,thermal,circuit|with_cx(|cx,_| {
        let prepared=assembly.prepare(cx,&thermal).unwrap();let law=reaction(50.0);
        let mut state=MaterialExchangeState::from_parts(thermal,circuit,gas_map(&law));
        let chambers=[gas_fixture::chamber(20,&state.chambers()[&20],0.00005),
            gas_fixture::chamber(10,&state.chambers()[&10],0.00005)];
        let reactions=[ChamberReaction {chamber:10,reaction:&law,controls:chemical()}];
        let drives=[loads(20,&state.chambers()[&20]),loads(10,&state.chambers()[&10])];
        let direct=prepared.reactive_step_candidate(cx,state.thermal(),0.00005,state.circuit(),&[],
            &chambers,&reactions,&gas_fixture::patches(),controls(),&[0.0;2]).unwrap();
        let result=state.advance(&prepared,cx,0.00005,&[],&drives,&reactions,
            &gas_fixture::patches(),controls(),&[0.0;2]).unwrap();
        assert_eq!(state.time_s(),0.00005);
        assert_eq!(state.thermal().cells(),direct.thermal().cells());
        assert_eq!(state.circuit(),direct.circuit());
        assert_eq!(result.step,*direct.step());
        assert_eq!(result.report.source.reactions.as_slice(),direct.reactions());
        for (account,chamber) in result.report.source.chambers.iter().zip(direct.chambers()) {
            assert_eq!(account.id,chamber.id);
            assert_eq!(&state.chambers()[&chamber.id],&chamber.piston);
            assert_eq!(state.chambers()[&chamber.id].last_window_duration_s(),chamber.piston.last_window_duration_s());
            assert_eq!(state.chambers()[&chamber.id].last_window_transfers(),chamber.piston.last_window_transfers());
            assert_eq!(account.surfaces,chamber.surfaces);
            assert_eq!(account.boundary,chamber.boundary);
        }
        assert_eq!(result.report.electrical.contacts,direct.report().electrical.contacts);
        assert!(result.report.electrical.contacts[0].heat_j>0.0);
    }));
}

#[test]
fn preview_and_resuming_the_owned_state_do_not_double_react_or_repeat_a_tick() {
    fixture([0.8,0.8],150.0,|assembly,thermal,circuit|with_cx(|cx,_| {
        let prepared=assembly.prepare(cx,&thermal).unwrap();let law=reaction(50.0);
        let mut state=MaterialExchangeState::from_parts(thermal,circuit,gas_map(&law));
        let drives=[loads(10,&state.chambers()[&10]),loads(20,&state.chambers()[&20])];
        let reactions=[ChamberReaction {chamber:10,reaction:&law,controls:chemical()}];
        for _ in 0..8 {
            let target=state.time_s()+0.00005;let dt=target-state.time_s();
            let cells=state.thermal().cells().clone();let electrical=state.circuit().clone();
            let gases=state.chambers().clone();let time=state.time_s();
            let chemistry=law.parcel_candidate(state.chambers()[&10].gas(),dt,chemical(),||false).unwrap();
            let (preview,preview_report)=state.step_candidate(&prepared,cx,target,&[],&drives,&reactions,
                &gas_fixture::patches(),controls(),&[0.0;2]).unwrap();
            assert_eq!(state.thermal().cells(),&cells);assert_eq!(state.circuit(),&electrical);
            assert_eq!(state.chambers(),&gases);assert_eq!(state.time_s(),time);
            let report=state.advance(&prepared,cx,target,&[],&drives,&reactions,
                &gas_fixture::patches(),controls(),&[0.0;2]).unwrap();
            assert_eq!(&report.report.source.reactions[0].reaction,chemistry.step());
            assert_eq!(state.thermal().cells(),preview.thermal().cells());
            assert_eq!(state.circuit(),preview.circuit());assert_eq!(state.chambers(),preview.chambers());
            assert_eq!(report.step,preview_report.step);
            assert_eq!(state.time_s(),target);
            assert!(state.circuit().charges_c()[0]<electrical.charges_c()[0]);
            // Move, do not reconstruct any gas or reset electrical/port history.
            let (t,c,g)=state.into_parts();state=MaterialExchangeState::from_parts(t,c,g);
            assert_eq!(state.chambers(),preview.chambers());
            for (&id,gas) in state.chambers() {
                assert_eq!(gas.last_window_duration_s(),Some(dt));
                assert_eq!(gas.last_window_transfers(),preview.chambers()[&id].last_window_transfers());
            }
        }
    }));
}

#[test]
fn bad_complete_drive_maps_failed_solves_and_cancellation_preserve_joint_state() {
    fixture([0.2,0.2],150.0,|assembly,thermal,circuit|with_cx(|cx,gate| {
        let prepared=assembly.prepare(cx,&thermal).unwrap();let law=reaction(50.0);
        let mut state=MaterialExchangeState::from_parts(thermal,circuit,gas_map(&law));
        let drives=[loads(10,&state.chambers()[&10]),loads(20,&state.chambers()[&20])];
        let reactions=[ChamberReaction {chamber:10,reaction:&law,controls:chemical()}];
        let cells=state.thermal().cells().clone();let circuit=state.circuit().clone();let gases=state.chambers().clone();
        for bad in [vec![drives[0]],vec![drives[0];2],
            vec![drives[0],ChamberDrive {chamber:99,..drives[1]}]] {
            assert!(state.advance(&prepared,cx,0.00005,&[],&bad,&reactions,&gas_fixture::patches(),controls(),&[0.0;2]).is_err());
        }
        for target in [0.0,-1.0,f64::NAN,f64::INFINITY] {
            assert!(state.advance(&prepared,cx,target,&[],&drives,&reactions,&gas_fixture::patches(),controls(),&[0.0;2]).is_err());
        }
        let mut limited=controls();limited.network.maximum_circuit_trials=1;
        assert!(state.advance(&prepared,cx,0.00005,&[],&drives,&reactions,&gas_fixture::patches(),limited,&[0.0;2]).is_err());
        assert_eq!(state.thermal().cells(),&cells);assert_eq!(state.circuit(),&circuit);
        assert_eq!(state.chambers(),&gases);assert_eq!(state.time_s(),0.0);
        let report=state.advance(&prepared,cx,0.00005,&[],&drives,&reactions,&gas_fixture::patches(),controls(),&[0.0;2]).unwrap();
        assert!(report.report.source.reactions[0].reaction.extent_mol>0.0);
        let cells=state.thermal().cells().clone();let circuit=state.circuit().clone();let gases=state.chambers().clone();
        gate.request();
        assert!(matches!(state.advance(&prepared,cx,0.0001,&[],&drives,&reactions,&gas_fixture::patches(),controls(),&[0.0;2]),
            Err(ReactiveExchangeError::Cancelled)));
        assert_eq!(state.thermal().cells(),&cells);assert_eq!(state.circuit(),&circuit);
        assert_eq!(state.chambers(),&gases);assert_eq!(state.time_s(),0.00005);
        for (&id,gas) in state.chambers() {
            assert_eq!(gas.last_window_transfers(),gases[&id].last_window_transfers());
            assert_eq!(gas.last_window_duration_s(),gases[&id].last_window_duration_s());
        }
    }));
}

#[test]
fn chamber_free_state_keeps_the_same_finite_electrical_thermal_evolution() {
    fixture([0.2,0.2],150.0,|assembly,thermal,circuit|with_cx(|cx,_| {
        let prepared=assembly.prepare(cx,&thermal).unwrap();
        let direct=prepared.step_candidate(cx,&thermal,0.00005,&circuit,&[],&[],&[],controls(),&[0.0;2]).unwrap();
        let mut state=MaterialExchangeState::from_parts(thermal,circuit,BTreeMap::new());
        let report=state.advance(&prepared,cx,0.00005,&[],&[],&[],&[],controls(),&[0.0;2]).unwrap();
        assert_eq!(state.thermal().cells(),direct.thermal().cells());assert_eq!(state.circuit(),direct.circuit());
        assert!(state.chambers().is_empty());assert!(report.report.source.reactions.is_empty());
        assert_eq!(report.step,*direct.step());assert_eq!(state.time_s(),0.00005);
    }));
}


#[path = "flow.rs"]
mod flow;

#[path = "state/coupled.rs"]
mod coupled;
