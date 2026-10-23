//! Transported species reach the same live typed gas-heat assembly before publication.
use super::*;
use fs_couple::cavity_phs::piston::thermal::ConductanceSample;
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::material::exchange::{
    MaterialExchangeSample, gas_law::GasSurface,
};
use fs_material::state_point::ResolvedMaterialStatePoint;

fn surfaces() -> Vec<GasSurface> {
    gas_fixture::patches().iter().map(|p| GasSurface {
        contact: p.contact, area_m2: p.geometry.area_m2,
    }).collect()
}

#[test]
fn transported_composition_reaches_heat_law_and_live_advance_once() {
    fixture([0.8,0.8],150.0,|assembly,thermal,circuit|with_cx(|cx,_| {
        let prepared=assembly.prepare(cx,&thermal).unwrap();
        let mut state=MaterialExchangeState::from_parts(thermal,circuit,gases(0.004));
        let drives=drives(state.chambers());
        let patches=gas_fixture::patches();
        let geometry=surfaces();
        let saw_mixed=std::cell::Cell::new(false);

        let (preview,report)=state.flow_step_with_gas_law_candidate(
            &prepared,cx,0.00005,&[],&drives,&[],&geometry,&links(),
            flow_controls(),controls(),&[0.0;2],
            |contact,state,source,gas| {
                if contact.chamber==10 && gas.components().len()==2 {
                    saw_mixed.set(true);
                }
                let patch=*patches.iter().find(|p|p.contact==contact).unwrap();
                prepared.resolve_conductive_gas_contact(patch,state,source,gas)
            }).unwrap();
        assert!(saw_mixed.get());
        assert!(report.flow.links[0].transfer.transferred_mass_kg>0.0);
        assert_eq!(preview.chambers()[&10].gas().model().thermodynamics().components().len(),2);
        assert_eq!(state.time_s(),0.0);

        let accepted=state.advance_with_flow_and_gas_law(
            &prepared,cx,0.00005,&[],&drives,&[],&geometry,&links(),
            flow_controls(),controls(),&[0.0;2],
            |contact,state,source,gas| {
                let patch=*patches.iter().find(|p|p.contact==contact).unwrap();
                prepared.resolve_conductive_gas_contact(patch,state,source,gas)
            }).unwrap();
        assert_eq!(state.time_s(),0.00005);
        assert_eq!(state.chambers(),preview.chambers());
        assert_eq!(state.thermal().cells(),preview.thermal().cells());
        assert_eq!(state.circuit(),preview.circuit());
        assert_eq!(accepted.total_mass_residual_kg,report.total_mass_residual_kg);
        assert_eq!(accepted.total_energy_residual_j,report.total_energy_residual_j);
        for gas in state.chambers().values() {
            assert_eq!(gas.last_window_duration_s(),Some(0.00005));
            assert!(matches!(gas.audit_last_window().unwrap(),ModuleWindowAudit::Unavailable {..}));
        }
    }));
}

#[test]
fn typed_heat_source_failure_after_flow_preserves_every_live_owner() {
    fixture([0.2,0.2],150.0,|assembly,thermal,circuit|with_cx(|cx,_| {
        let prepared=assembly.prepare(cx,&thermal).unwrap();
        let mut state=MaterialExchangeState::from_parts(thermal,circuit,gases(0.004));
        let drives=drives(state.chambers());
        let geometry=surfaces();
        let before_cells=state.thermal().cells().clone();
        let before_circuit=state.circuit().clone();
        let before_gases=state.chambers().clone();

        let result=state.advance_with_flow_and_gas_law::<ResolvedMaterialStatePoint,_>(
            &prepared,cx,0.00005,&[],&drives,&[],&geometry,&links(),
            flow_controls(),controls(),&[0.0;2],
            |_,_,_,_| Err(MaterialNetworkError::Invalid("deliberate gas heat refusal")));
        assert!(result.is_err());
        assert_eq!(state.time_s(),0.0);
        assert_eq!(state.thermal().cells(),&before_cells);
        assert_eq!(state.circuit(),&before_circuit);
        assert_eq!(state.chambers(),&before_gases);
    }));
}
