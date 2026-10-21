//! Existing complete-assembly fixtures; no replacement thermal/electrical solver.
use super::*;
use fs_material::thermochemical_gas::energy::parcel::transfer::{
    nozzle::network::snapshot::SnapshotNozzleLink,
    species::{network::{SpeciesNetworkSolve, transport::{self as shared, GasTransportSolve}},
        diffusion::{BinaryDiffusionMaterial, GasDiffusionGeometry, BINARY_DIFFUSIVITY_DIMS,
            BINARY_DIFFUSIVITY_PROPERTY, network::{DiffusionFaceSolve, DiffusionNetworkLink,
                DiffusionNetworkSolve}}},
};

const DT: f64 = 0.00005;
fn admission() -> GasTransportSolve {
    GasTransportSolve { exchange: SpeciesNetworkSolve { maximum_species_fraction: 0.05,
        mass_tolerance_kg: 1e-13, energy_tolerance_j: 1e-7 }, entropy_tolerance_j_k: 1e-10 }
}
fn coupled() -> NetworkSolve {
    NetworkSolve { maximum_chambers: 3, maximum_boundary_solves: 1200, ..controls() }
}
fn pair_card(value: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(BINARY_DIFFUSIVITY_PROPERTY, BINARY_DIFFUSIVITY_DIMS),
        value: PropertyValue::Scalar { value, dims: BINARY_DIFFUSIVITY_DIMS },
        validity: ValidityDomain::unconstrained()
            .with_quantity("temperature", QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature, ValueForm::Static)), 200.0, 1800.0)
            .with_quantity("pressure", QuantitySpec::semantic(SemanticType::new(
                QuantityKind::Pressure, ValueForm::Static)), 1000.0, 1e7),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![], provenance: Provenance { source: "synthetic joint transport pair".into(),
            license: "CC0-1.0".into(), artifact: None },
    }).unwrap();
    MaterialCard::assemble(MaterialStateId { chemistry: "synthetic A/B".into(), phase: "gas".into(),
        process: "test".into(), revision: 0 }, claims, vec![]).unwrap()
}
fn binding<'a>(card: &'a MaterialCard, point: &'a QueryPoint) -> [BinaryDiffusionMaterial<'a>; 1] {
    [BinaryDiffusionMaterial { species: [SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()],
        material: card, context: point, selection: MaterialPropertySelection::SingleClaimOnly }]
}
fn face<'a>(pairs: &'a [BinaryDiffusionMaterial<'a>]) -> DiffusionNetworkLink<'a> {
    DiffusionNetworkLink { id: 27, parcels: [10,20], pairs,
        geometry: GasDiffusionGeometry { area_m2: 0.1, distance_m: 0.0001 },
        solve: DiffusionFaceSolve { maximum_temperature_difference_k: 5.0,
            maximum_pressure_difference_pa: 2000.0, gradient_tolerance_per_m: 1e-8,
            molar_flux_tolerance_mol_m2_s: 1e-10 } }
}
fn restriction() -> SnapshotNozzleLink {
    SnapshotNozzleLink { id: 901, parcels: [10,30], nozzle: GasNozzle {
        area_m2: 2e-6, discharge_coefficient: 0.8 }, nozzle_solve: chemical().gas }
}
fn chambers() -> BTreeMap<u64, GasPiston> {
    BTreeMap::from([(10,source_gas("A",0.028,100000.0)),
        (20,source_gas("B",0.004,100000.0)),(30,source_gas("A",0.028,60000.0))])
}
fn storage(state: &MaterialExchangeState<'_>) -> f64 {
    state.thermal().energy_j() + state.circuit().charges_c()[0].powi(2)/(2.0*0.01)
        + state.chambers().values().map(|g| g.gas().mass_kg()*g.gas().state().specific_internal_energy_j_kg()
            + 0.5*g.mass_kg()*g.velocity_m_s().powi(2)).sum::<f64>()
}

#[test]
fn original_material_stage_enters_one_joint_assembly_solve() {
    fixture([0.2,0.2],150.0,|assembly,thermal,circuit|with_cx(|cx,_| {
        let d=pair_card(2e-5);let q=QueryPoint::new();let pairs=binding(&d,&q);let faces=[face(&pairs)];
        let prepared=assembly.prepare(cx,&thermal).unwrap();
        let state=MaterialExchangeState::from_parts(thermal,circuit,chambers());let drives=drives(state.chambers());
        let parcels=state.chambers().iter().map(|(&id,g)|(id,g.gas().clone())).collect();
        let endpoints=drives.iter().map(|d|(d.chamber,d.piston_solve.gas)).collect();
        let (transferred,flux)=shared::step_candidate(&parcels,&[restriction()],&faces,DT,&endpoints,
            admission(),||false).unwrap().into_parts();
        let staged=MaterialExchangeState::from_parts(state.thermal().clone(),state.circuit().clone(),
            transferred.into_iter().map(|(id,g)|(id,GasPiston::new(g,1.0,0.01,0.0).unwrap())).collect());
        let (expected,_)=staged.step_candidate(&prepared,cx,DT,&[],&drives,&[],
            &gas_fixture::patches(),coupled(),&[0.0;2]).unwrap();
        let (actual,report)=state.transport_step_candidate(&prepared,cx,DT,&[],&drives,&[],
            &gas_fixture::patches(),&[restriction()],&faces,admission(),coupled(),&[0.0;2]).unwrap();
        assert_eq!(report.flow,flux);
        // Independent staging has a different port-account reference. Compare
        // physical endpoints, then independently check the actual full interval.
        for (&id,gas) in actual.chambers() {
            assert_eq!(gas.gas(),expected.chambers()[&id].gas());
            assert_eq!(gas.velocity_m_s(),expected.chambers()[&id].velocity_m_s());
            assert_eq!(gas.position_m(),expected.chambers()[&id].position_m());
            assert_eq!(gas.last_window_duration_s(),Some(DT));
            assert!(matches!(gas.audit_last_window().unwrap(),ModuleWindowAudit::Unavailable {..}));
        }
        assert_eq!(actual.circuit(),expected.circuit());assert_eq!(actual.thermal().cells(),expected.thermal().cells());
        let work:f64=report.exchange.report.source.chambers.iter().map(|c|c.boundary.piston.exterior_pressure_work_j).sum();
        assert!((storage(&actual)+work-storage(&state)).abs()<4e-6);
        let (diff,_)=state.diffusion_step_candidate(&prepared,cx,DT,&[],&drives,&[],&gas_fixture::patches(),
            &faces,DiffusionNetworkSolve {exchange:admission().exchange,
                entropy_tolerance_j_k:admission().entropy_tolerance_j_k},coupled(),&[0.0;2]).unwrap();
        let (only,_)=state.transport_step_candidate(&prepared,cx,DT,&[],&drives,&[],&gas_fixture::patches(),
            &[],&faces,admission(),coupled(),&[0.0;2]).unwrap();
        assert_eq!(only.chambers(),diff.chambers());assert_eq!(only.circuit(),diff.circuit());
        assert_eq!(only.thermal().cells(),diff.thermal().cells());
        let (plain,_)=state.step_candidate(&prepared,cx,DT,&[],&drives,&[],&gas_fixture::patches(),coupled(),&[0.0;2]).unwrap();
        let (idle,_)=state.transport_step_candidate(&prepared,cx,DT,&[],&drives,&[],&gas_fixture::patches(),
            &[],&[],admission(),coupled(),&[0.0;2]).unwrap();
        assert_eq!(idle.chambers(),plain.chambers());assert_eq!(idle.circuit(),plain.circuit());
        assert_eq!(idle.thermal().cells(),plain.thermal().cells());assert_eq!(state.time_s(),0.0);
    }));
}

#[test]
fn pair_data_changes_downstream_flow_heat_charge_motion_and_sound() {
    let run=|coefficient|fixture([0.8,0.8],150.0,|assembly,thermal,circuit|with_cx(|cx,_| {
        let d=pair_card(coefficient);let q=QueryPoint::new();let pairs=binding(&d,&q);let faces=[face(&pairs)];
        let prepared=assembly.prepare(cx,&thermal).unwrap();
        let mut state=MaterialExchangeState::from_parts(thermal,circuit,chambers());
        let drives=drives(state.chambers());let start=storage(&state);let mut work=0.0;
        for _ in 0..16 {
            let target=state.time_s()+DT;
            let report=state.advance_with_transport(&prepared,cx,target,&[],&drives,&[],&gas_fixture::patches(),
                &[restriction()],&faces,admission(),coupled(),&[0.0;2]).unwrap();
            work+=report.exchange.report.source.chambers.iter().map(|c|c.boundary.piston.exterior_pressure_work_j).sum::<f64>();
            assert!(report.total_energy_residual_j.abs()<=coupled().total_energy_tolerance_j);
            assert!((storage(&state)+work-start).abs()<4e-5);
        }
        let sink=state.chambers()[&30].gas();assert_eq!(sink.model().thermodynamics().components().len(),2);
        (sink.mass_kg(),state.thermal().temperature_k(CellId(30)).unwrap(),state.circuit().charges_c()[0],
            state.chambers()[&10].position_m(),state.chambers()[&10].gas().state().state().sound_speed)
    }));
    let slow=run(2e-5);let fast=run(4e-5);
    assert_ne!(slow.0,fast.0);assert_ne!(slow.1,fast.1);assert_ne!(slow.2,fast.2);
    assert_ne!(slow.3,fast.3);assert!((slow.4-fast.4).abs()>0.01);
}

#[test]
fn chemistry_consumes_combined_arrivals_once_not_per_nonlinear_trial() {
    fixture([0.2,0.2],150.0,|assembly,thermal,circuit|with_cx(|cx,_| {
        let law=reaction(50.0);let d=pair_card(2e-5);let q=QueryPoint::new();let pairs=binding(&d,&q);let faces=[face(&pairs)];
        let gases=[(10,0.75,100000.0),(20,0.25,100000.0),(30,0.75,60000.0)].into_iter().map(|(id,x,p)| {
            let model=law.catalog().with_composition(Composition::new(CompositionBasis::MoleFraction,
                vec![x,1.0-x]).unwrap()).unwrap();(id,parcel_piston(model,p))
        }).collect();
        let prepared=assembly.prepare(cx,&thermal).unwrap();
        let state=MaterialExchangeState::from_parts(thermal,circuit,gases);let drives=drives(state.chambers());
        let parcels=state.chambers().iter().map(|(&id,g)|(id,g.gas().clone())).collect();
        let endpoints=drives.iter().map(|d|(d.chamber,d.piston_solve.gas)).collect();
        let moved=shared::step_candidate(&parcels,&[restriction()],&faces,DT,&endpoints,admission(),||false).unwrap();
        let expected=law.parcel_candidate(&moved.parcels()[&10],DT,chemical(),||false).unwrap();
        let (next,report)=state.transport_step_candidate(&prepared,cx,DT,&[],&drives,
            &[ChamberReaction {chamber:10,reaction:&law,controls:chemical()}],&gas_fixture::patches(),
            &[restriction()],&faces,admission(),coupled(),&[0.0;2]).unwrap();
        assert_eq!(report.exchange.report.source.reactions[0].reaction,*expected.step());
        assert!(report.total_energy_residual_j.abs()<=coupled().total_energy_tolerance_j);
        assert_eq!(next.time_s(),DT);assert_eq!(state.time_s(),0.0);
    }));
}

#[test]
fn late_feedback_failure_and_cancellation_publish_no_partial_assembly() {
    fixture([0.2,0.2],150.0,|assembly,thermal,circuit|with_cx(|cx,gate| {
        let d=pair_card(2e-5);let q=QueryPoint::new();let pairs=binding(&d,&q);let faces=[face(&pairs)];
        let prepared=assembly.prepare(cx,&thermal).unwrap();
        let mut state=MaterialExchangeState::from_parts(thermal,circuit,chambers());let drives=drives(state.chambers());
        let gases=state.chambers().clone();let cells=state.thermal().cells().clone();let electrical=state.circuit().clone();
        let mut limited=coupled();limited.network.maximum_circuit_trials=1;
        assert!(state.advance_with_transport(&prepared,cx,DT,&[],&drives,&[],&gas_fixture::patches(),
            &[restriction()],&faces,admission(),limited,&[0.0;2]).is_err());
        assert_eq!(state.chambers(),&gases);assert_eq!(state.thermal().cells(),&cells);
        assert_eq!(state.circuit(),&electrical);assert_eq!(state.time_s(),0.0);
        let (preview,_)=state.transport_step_candidate(&prepared,cx,DT,&[],&drives,&[],&gas_fixture::patches(),
            &[restriction()],&faces,admission(),coupled(),&[0.0;2]).unwrap();
        state.advance_with_transport(&prepared,cx,DT,&[],&drives,&[],&gas_fixture::patches(),
            &[restriction()],&faces,admission(),coupled(),&[0.0;2]).unwrap();
        assert_eq!(state.chambers(),preview.chambers());assert_eq!(state.circuit(),preview.circuit());
        assert_eq!(state.thermal().cells(),preview.thermal().cells());gate.request();
        assert!(matches!(state.advance_with_transport(&prepared,cx,2.0*DT,&[],&drives,&[],&gas_fixture::patches(),
            &[restriction()],&faces,admission(),coupled(),&[0.0;2]),Err(MaterialFlowError::Cancelled)));
        assert_eq!(state.time_s(),DT);assert_eq!(state.chambers(),preview.chambers());
        assert_eq!(state.circuit(),preview.circuit());assert_eq!(state.thermal().cells(),preview.thermal().cells());
    }));
}
