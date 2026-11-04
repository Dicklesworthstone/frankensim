//! G1/G3/G4: source diffusion changes the original complete physical assembly.
//! Manufactured properties; no measured engine or diffusion validation claim.
use super::*;
use fs_material::thermochemical_gas::energy::parcel::transfer::species::{
    network::SpeciesNetworkSolve,
    diffusion::{BinaryDiffusionMaterial, GasDiffusionDriving, GasDiffusionGeometry, BINARY_DIFFUSIVITY_DIMS,
        BINARY_DIFFUSIVITY_PROPERTY, network::{self as diffusion, DiffusionFaceSolve,
            DiffusionNetworkLink, DiffusionNetworkSolve}},
};
use diffusion::adaptive::{DiffusionNetworkAdvanceControls, DiffusionNetworkAdvanceError};

const DT: f64 = 0.00005;

fn pair_card(diffusivity: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(BINARY_DIFFUSIVITY_PROPERTY, BINARY_DIFFUSIVITY_DIMS),
        value: PropertyValue::Scalar { value: diffusivity, dims: BINARY_DIFFUSIVITY_DIMS },
        validity: ValidityDomain::unconstrained()
            .with_quantity("temperature", QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature, ValueForm::Static)), 200.0, 1800.0)
            .with_quantity("pressure", QuantitySpec::semantic(SemanticType::new(
                QuantityKind::Pressure, ValueForm::Static)), 1000.0, 1e7),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![], provenance: Provenance { source: "synthetic diffusion pair".into(),
            license: "CC0-1.0".into(), artifact: None },
    }).unwrap();
    MaterialCard::assemble(MaterialStateId { chemistry: "synthetic A/B".into(), phase: "gas".into(),
        process: "diffusion fixture".into(), revision: 0 }, claims, vec![]).unwrap()
}
fn binding<'a>(card: &'a MaterialCard, point: &'a QueryPoint) -> BinaryDiffusionMaterial<'a> {
    BinaryDiffusionMaterial { species: [SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()],
        material: card, context: point, selection: MaterialPropertySelection::SingleClaimOnly }
}
fn face<'a>(pairs: &'a [BinaryDiffusionMaterial<'a>]) -> DiffusionNetworkLink<'a> {
    DiffusionNetworkLink { id: 27, parcels: [10,20], pairs,
        geometry: GasDiffusionGeometry { area_m2: 0.1, distance_m: 0.0001 },
        solve: DiffusionFaceSolve { maximum_temperature_difference_k: 5.0,
            maximum_pressure_difference_pa: 2000.0, gradient_tolerance_per_m: 1e-8,
            molar_flux_tolerance_mol_m2_s: 1e-10 } }
}
fn diffusion_controls() -> DiffusionNetworkSolve {
    DiffusionNetworkSolve { exchange: SpeciesNetworkSolve { maximum_species_fraction: 0.05,
        mass_tolerance_kg: 1e-13, energy_tolerance_j: 1e-7 }, entropy_tolerance_j_k: 1e-10 }
}
fn mixtures() -> BTreeMap<u64, GasPiston> {
    let (a,ta) = gas_fixture::component("A",0.028,0.0);
    let (b,tb) = gas_fixture::component("B",0.004,0.0);
    [0.75,0.25].into_iter().enumerate().map(|(i,x)| {
        let model = ThermochemicalGasModel::new(FrozenIdealGasMixtureModelV1::new(vec![a.clone(),b.clone()],
            Composition::new(CompositionBasis::MoleFraction,vec![x,1.0-x]).unwrap()).unwrap(),
            vec![ta.clone(),tb.clone()],None).unwrap();
        (10+10*i as u64,parcel_piston(model,100000.0))
    }).collect()
}
fn storage(state: &MaterialExchangeState<'_>) -> f64 {
    state.thermal().energy_j() + state.circuit().charges_c()[0].powi(2)/(2.0*0.01)
        + state.chambers().values().map(|g| g.gas().mass_kg()*g.gas().state().specific_internal_energy_j_kg()
            + 0.5*g.mass_kg()*g.velocity_m_s().powi(2)).sum::<f64>()
}

#[test]
fn equal_pressure_diffusion_matches_original_owners_and_closed_faces_preserve_them() {
    fixture([0.2,0.2],150.0,|assembly,thermal,circuit|with_cx(|cx,_| {
        let card = pair_card(2e-5); let context = QueryPoint::new();
        let pairs = [binding(&card,&context)]; let links = [face(&pairs)];
        let prepared = assembly.prepare(cx,&thermal).unwrap();
        let state = MaterialExchangeState::from_parts(thermal,circuit,mixtures());
        let drives = drives(state.chambers());
        let parcels = state.chambers().iter().map(|(&id,g)|(id,g.gas().clone())).collect();
        let endpoints = drives.iter().map(|d|(d.chamber,d.piston_solve.gas)).collect();
        let (transferred, flux) = diffusion::step_candidate(&parcels,&links,DT,&endpoints,
            diffusion_controls(),||false).unwrap().into_parts();
        // Independent binary/Fick limit: c D (x_left-x_right)/distance.
        let expected = 100000.0/(UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K*400.0)*2e-5*0.5/0.0001;
        assert!((flux.links[0].flux.species[0].molar_flux_mol_m2_s-expected).abs()<1e-12);
        assert!((flux.links[0].flux.species[1].molar_flux_mol_m2_s+expected).abs()<1e-12);
        let staged = MaterialExchangeState::from_parts(state.thermal().clone(),state.circuit().clone(),
            transferred.into_iter().map(|(id,g)|(id,GasPiston::new(g,1.0,0.01,0.0).unwrap())).collect());
        let (expected,_) = staged.step_candidate(&prepared,cx,DT,&[],&drives,&[],
            &gas_fixture::patches(),controls(),&[0.0;2]).unwrap();
        let (actual,report) = state.diffusion_step_candidate(&prepared,cx,DT,&[],&drives,&[],
            &gas_fixture::patches(),&links,diffusion_controls(),controls(),&[0.0;2]).unwrap();
        assert_eq!(report.flow,flux);
        assert_eq!(actual.chambers(),expected.chambers());
        assert_eq!(actual.circuit(),expected.circuit());
        assert_eq!(actual.thermal().cells(),expected.thermal().cells());
        let work: f64 = report.exchange.report.source.chambers.iter()
            .map(|c|c.boundary.piston.exterior_pressure_work_j).sum();
        assert!((storage(&actual)+work-storage(&state)).abs()<4e-6);
        for gas in actual.chambers().values() {
            assert_eq!(gas.last_window_duration_s(),Some(DT));
            assert!(matches!(gas.audit_last_window().unwrap(),ModuleWindowAudit::Unavailable {..}));
        }
        let (plain,_) = state.step_candidate(&prepared,cx,DT,&[],&drives,&[],
            &gas_fixture::patches(),controls(),&[0.0;2]).unwrap();
        let closed = DiffusionNetworkLink { pairs:&[], geometry:GasDiffusionGeometry {
            area_m2:0.0,..links[0].geometry },..links[0].clone() };
        for links in [&[][..],std::slice::from_ref(&closed)] {
            let (idle,report) = state.diffusion_step_candidate(&prepared,cx,DT,&[],&drives,&[],
                &gas_fixture::patches(),links,diffusion_controls(),controls(),&[0.0;2]).unwrap();
            assert_eq!(idle.chambers(),plain.chambers());assert_eq!(idle.circuit(),plain.circuit());
            assert_eq!(idle.thermal().cells(),plain.thermal().cells());
            assert!(report.flow.exchange.fluxes.is_empty());
        }
        assert_eq!(state.time_s(),0.0);
    }));
}

#[test]
fn changing_only_pair_diffusivity_changes_heat_charge_motion_and_sound_speed() {
    let run = |d|fixture([0.8,0.8],150.0,|assembly,thermal,circuit|with_cx(|cx,_| {
        let card=pair_card(d);let context=QueryPoint::new();let pairs=[binding(&card,&context)];
        let links=[face(&pairs)];let prepared=assembly.prepare(cx,&thermal).unwrap();
        let mut state=MaterialExchangeState::from_parts(thermal,circuit,mixtures());
        let drives=drives(state.chambers());let before=storage(&state);let mut work=0.0;
        let mass:f64=state.chambers().values().map(|g|g.gas().mass_kg()).sum();
        for _ in 0..32 {
            let report=state.advance_with_diffusion(&prepared,cx,state.time_s()+DT,&[],&drives,&[],
                &gas_fixture::patches(),&links,diffusion_controls(),controls(),&[0.0;2]).unwrap();
            work+=report.exchange.report.source.chambers.iter()
                .map(|c|c.boundary.piston.exterior_pressure_work_j).sum::<f64>();
            assert!((storage(&state)+work-before).abs()<4e-5);
            assert!(report.total_energy_residual_j.abs()<=controls().total_energy_tolerance_j);
            assert!(report.flow.exchange.maximum_species_mass_residual_kg<1e-12);
            assert!(report.flow.entropy_change_j_k>=-diffusion_controls().entropy_tolerance_j_k);
        }
        assert!((state.chambers().values().map(|g|g.gas().mass_kg()).sum::<f64>()-mass).abs()<1e-12);
        (state.thermal().temperature_k(CellId(30)).unwrap(),state.circuit().charges_c()[0],
            state.chambers()[&10].position_m(),state.chambers()[&10].gas().state().state().sound_speed)
    }));
    let slow=run(2e-5);let fast=run(4e-5);
    assert!((slow.0-fast.0).abs()>1e-9);assert!((slow.1-fast.1).abs()>1e-14);
    assert!((slow.2-fast.2).abs()>1e-13);assert!((slow.3-fast.3).abs()>0.01);
}

#[test]
fn reactions_use_diffused_inventory_once_even_with_zero_net_mass_transfer() {
    fixture([0.2,0.2],150.0,|assembly,thermal,circuit|with_cx(|cx,_| {
        let law=reaction(50.0);let card=pair_card(2e-5);let context=QueryPoint::new();
        let pairs=[binding(&card,&context)];let links=[face(&pairs)];
        let gases=[0.75,0.25].into_iter().enumerate().map(|(i,x)| {
            let model=law.catalog().with_composition(Composition::new(CompositionBasis::MoleFraction,
                vec![x,1.0-x]).unwrap()).unwrap();
            (10+10*i as u64,parcel_piston(model,100000.0))
        }).collect();
        let prepared=assembly.prepare(cx,&thermal).unwrap();
        let state=MaterialExchangeState::from_parts(thermal,circuit,gases);let drives=drives(state.chambers());
        let parcels=state.chambers().iter().map(|(&id,g)|(id,g.gas().clone())).collect();
        let endpoints=drives.iter().map(|d|(d.chamber,d.piston_solve.gas)).collect();
        let diffused=diffusion::step_candidate(&parcels,&links,DT,&endpoints,diffusion_controls(),||false).unwrap();
        let expected=law.parcel_candidate(&diffused.parcels()[&10],DT,chemical(),||false).unwrap();
        let original=law.parcel_candidate(state.chambers()[&10].gas(),DT,chemical(),||false).unwrap();
        assert_ne!(expected.step().extent_mol,original.step().extent_mol);
        let (next,report)=state.diffusion_step_candidate(&prepared,cx,DT,&[],&drives,
            &[ChamberReaction {chamber:10,reaction:&law,controls:chemical()}],&gas_fixture::patches(),
            &links,diffusion_controls(),controls(),&[0.0;2]).unwrap();
        assert_eq!(report.exchange.report.source.reactions[0].reaction,*expected.step());
        assert!(report.total_energy_residual_j.abs()<=controls().total_energy_tolerance_j);
        for parcel in &report.flow.exchange.parcels {
            assert!(parcel.incoming_mass_kg>0.0);
            assert!((parcel.incoming_mass_kg-parcel.outgoing_mass_kg).abs()<1e-15);
            assert!(matches!(next.chambers()[&parcel.id].audit_last_window().unwrap(),
                ModuleWindowAudit::Unavailable {..}));
        }
    }));
}

#[test]
fn missing_sources_and_late_feedback_refusals_preserve_every_owner_and_retry_once() {
    fixture([0.2,0.2],150.0,|assembly,thermal,circuit|with_cx(|cx,gate| {
        let card=pair_card(2e-5);let context=QueryPoint::new();let pairs=[binding(&card,&context)];
        let links=[face(&pairs)];let prepared=assembly.prepare(cx,&thermal).unwrap();
        let mut state=MaterialExchangeState::from_parts(thermal,circuit,mixtures());
        let drives=drives(state.chambers());let old_gas=state.chambers().clone();
        let old_cells=state.thermal().cells().clone();let old_circuit=state.circuit().clone();
        let missing=DiffusionNetworkLink {pairs:&[],..links[0].clone()};
        assert!(matches!(state.advance_with_diffusion(&prepared,cx,DT,&[],&drives,&[],
            &gas_fixture::patches(),&[missing],diffusion_controls(),controls(),&[0.0;2]),
            Err(MaterialFlowError::Diffusion(_))));
        let mut limited=controls();limited.network.maximum_circuit_trials=1;
        assert!(state.advance_with_diffusion(&prepared,cx,DT,&[],&drives,&[],
            &gas_fixture::patches(),&links,diffusion_controls(),limited,&[0.0;2]).is_err());
        assert_eq!(state.chambers(),&old_gas);assert_eq!(state.circuit(),&old_circuit);
        assert_eq!(state.thermal().cells(),&old_cells);assert_eq!(state.time_s(),0.0);
        let (preview,_)=state.diffusion_step_candidate(&prepared,cx,DT,&[],&drives,&[],
            &gas_fixture::patches(),&links,diffusion_controls(),controls(),&[0.0;2]).unwrap();
        state.advance_with_diffusion(&prepared,cx,DT,&[],&drives,&[],&gas_fixture::patches(),
            &links,diffusion_controls(),controls(),&[0.0;2]).unwrap();
        assert_eq!(state.chambers(),preview.chambers());assert_eq!(state.circuit(),preview.circuit());
        assert_eq!(state.thermal().cells(),preview.thermal().cells());
        gate.request();
        assert!(matches!(state.advance_with_diffusion(&prepared,cx,2.0*DT,&[],&drives,&[],
            &gas_fixture::patches(),&links,diffusion_controls(),controls(),&[0.0;2]),
            Err(MaterialFlowError::Cancelled)));
        assert_eq!(state.time_s(),DT);assert_eq!(state.chambers(),preview.chambers());
        assert_eq!(state.circuit(),preview.circuit());assert_eq!(state.thermal().cells(),preview.thermal().cells());
    }));
}

fn adaptive_controls() -> DiffusionNetworkAdvanceControls {
    DiffusionNetworkAdvanceControls {
        maximum_step_s: DT,
        minimum_step_s: 1e-12,
        maximum_evaluations: 10000,
        absolute_species_tolerance_kg: 1e-13,
        relative_species_tolerance: 1e-4,
        temperature_tolerance_k: 1e-6,
    }
}

#[test]
fn adaptive_diffusion_completes_a_stiff_window_before_whole_assembly_feedback() {
    fixture([0.2,0.2],150.0,|assembly,thermal,circuit|with_cx(|cx,_| {
        let card = pair_card(0.002); let context = QueryPoint::new();
        let pairs = [binding(&card,&context)]; let links = [face(&pairs)];
        let prepared = assembly.prepare(cx,&thermal).unwrap();
        let mut state = MaterialExchangeState::from_parts(thermal,circuit,mixtures());
        let drives = drives(state.chambers());
        // A one-shot transfer exceeds the actual donor allowance. The adaptive
        // assembly must complete this SAME interval without changing that gate.
        assert!(matches!(state.diffusion_step_candidate(&prepared,cx,DT,&[],&drives,&[],
            &gas_fixture::patches(),&links,diffusion_controls(),controls(),&[0.0;2]),
            Err(MaterialFlowError::Diffusion(_))));
        let parcels = state.chambers().iter().map(|(&id,g)|(id,g.gas().clone())).collect();
        let endpoints = drives.iter().map(|d|(d.chamber,d.piston_solve.gas)).collect();
        let (diffused, original) = diffusion::adaptive::advance_candidate(&parcels,&links,DT,
            &endpoints,diffusion_controls(),adaptive_controls(),||false).unwrap().into_parts();
        let x = |id| diffused[&id].model().thermodynamics().mole_composition().fractions()[0];
        let exact_gap = 0.5 * (-2.0_f64 * 0.002 * 0.1 / (0.0001 * 0.001) * DT).exp();
        assert!(((x(10)-x(20))-exact_gap).abs()<5e-4);
        let staged = MaterialExchangeState::from_parts(state.thermal().clone(),state.circuit().clone(),
            diffused.into_iter().map(|(id,g)|(id,GasPiston::new(g,1.0,0.01,0.0).unwrap())).collect());
        let (expected,_) = staged.step_candidate(&prepared,cx,DT,&[],&drives,&[],
            &gas_fixture::patches(),controls(),&[0.0;2]).unwrap();
        let report = state.advance_with_adaptive_diffusion(&prepared,cx,DT,&[],&drives,&[],
            &gas_fixture::patches(),&links,diffusion_controls(),adaptive_controls(),
            GasDiffusionDriving::IsobaricMolar,controls(),&[0.0;2]).unwrap();
        assert_eq!(report.flow,original);
        assert_eq!(state.chambers(),expected.chambers());
        assert_eq!(state.circuit(),expected.circuit());
        assert_eq!(state.thermal().cells(),expected.thermal().cells());
        assert_eq!(state.time_s(),DT);
        assert!(report.flow.rejected_intervals>0 && report.flow.accepted_intervals>1);
        assert_eq!(report.flow.accepted_substeps,2*report.flow.accepted_intervals);
        assert!(report.flow.maximum_species_error_ratio<=1.0);
        assert!(report.total_mass_residual_kg.abs()<=diffusion_controls().exchange.mass_tolerance_kg);
        assert!(report.total_energy_residual_j.abs()<=controls().total_energy_tolerance_j);
        for gas in state.chambers().values() {
            assert_eq!(gas.last_window_duration_s(),Some(DT));
            assert!(matches!(gas.audit_last_window().unwrap(),ModuleWindowAudit::Unavailable {..}));
        }
    }));
}

#[test]
fn adaptive_pressure_diffusion_changes_uniform_composition_in_the_mass_frame() {
    fixture([0.2,0.2],150.0,|assembly,thermal,circuit|with_cx(|cx,_| {
        let card = pair_card(0.002); let context = QueryPoint::new();
        let pairs = [binding(&card,&context)]; let links = [face(&pairs)];
        let prepared = assembly.prepare(cx,&thermal).unwrap();
        let model = mixtures()[&10].gas().model().clone();
        let chambers = BTreeMap::from([(10,parcel_piston(model.clone(),100500.0)),
            (20,parcel_piston(model,99500.0))]);
        let state = MaterialExchangeState::from_parts(thermal,circuit,chambers);
        let drives = drives(state.chambers());
        let (_,molar) = state.adaptive_diffusion_step_candidate(&prepared,cx,DT,&[],&drives,&[],
            &gas_fixture::patches(),&links,diffusion_controls(),adaptive_controls(),
            GasDiffusionDriving::IsobaricMolar,controls(),&[0.0;2]).unwrap();
        assert!(molar.flow.transfers.is_empty());
        let driving = GasDiffusionDriving::CompositionAndPressureMass {
            mass_flux_tolerance_kg_m2_s: 1e-12 };
        let (next,mass) = state.adaptive_diffusion_step_candidate(&prepared,cx,DT,&[],&drives,&[],
            &gas_fixture::patches(),&links,diffusion_controls(),adaptive_controls(),driving,
            controls(),&[0.0;2]).unwrap();
        assert!(!mass.flow.transfers.is_empty());
        for face in &mass.flow.last_substep.links {
            assert_eq!(face.flux.driving,driving);
            assert!(face.flux.mass_frame_residual_kg_m2_s.abs()<=1e-12);
        }
        assert!(mass.flow.entropy_change_j_k>0.0);
        for (&id,before) in state.chambers() {
            assert!((next.chambers()[&id].gas().mass_kg()-before.gas().mass_kg()).abs()<1e-13);
            assert_ne!(next.chambers()[&id].gas().model().thermodynamics().mole_composition(),
                before.gas().model().thermodynamics().mole_composition());
            assert!(matches!(next.chambers()[&id].audit_last_window().unwrap(),
                ModuleWindowAudit::Unavailable {..}));
        }
        assert!(mass.total_energy_residual_j.abs()<=controls().total_energy_tolerance_j);
        assert_eq!(state.time_s(),0.0);
    }));
}

#[test]
fn adaptive_arrivals_reach_a_third_live_chamber_during_the_same_window() {
    fixture([0.2,0.2],150.0,|assembly,thermal,circuit|with_cx(|cx,_| {
        let card = pair_card(0.002); let context = QueryPoint::new();
        let pairs = [binding(&card,&context)];
        let links = [face(&pairs),DiffusionNetworkLink { id:28,parcels:[20,30],..face(&pairs) }];
        let prepared = assembly.prepare(cx,&thermal).unwrap();
        let mut chambers = mixtures();
        chambers.insert(30,chambers[&20].clone());
        let state = MaterialExchangeState::from_parts(thermal,circuit,chambers);
        let drives = drives(state.chambers());
        let mut settings = controls(); settings.maximum_chambers = 3;
        let parcels = state.chambers().iter().map(|(&id,g)|(id,g.gas().clone())).collect();
        let endpoints = drives.iter().map(|d|(d.chamber,d.piston_solve.gas)).collect();
        let original = diffusion::step_candidate(&parcels,&links,DT/128.0,&endpoints,
            diffusion_controls(),||false).unwrap();
        assert!(original.step().links[1].flux.species.iter()
            .all(|row|row.mass_rate_left_to_right_kg_s==0.0));
        let (next,report) = state.adaptive_diffusion_step_candidate(&prepared,cx,DT,&[],&drives,&[],
            &gas_fixture::patches(),&links,diffusion_controls(),adaptive_controls(),
            GasDiffusionDriving::IsobaricMolar,settings,&[0.0;2]).unwrap();
        assert!(report.flow.transfers.iter().any(|row|row.donor==20 && row.receiver==30
            && row.species.as_str()=="A" && row.mass_kg>0.0));
        // Independent three-cell line Laplacian eigenmodes, rates k and 3k.
        let q = 0.002 * 0.1 / (0.0001 * 0.001) * DT;
        let exact = 5.0/12.0 - 0.25*(-q).exp() + (-3.0*q).exp()/12.0;
        let actual = next.chambers()[&30].gas().model().thermodynamics()
            .mole_composition().fractions()[0];
        assert!(actual>0.25 && (actual-exact).abs()<3e-4);
        assert_eq!(next.time_s(),DT);
        assert!(report.flow.maximum_species_mass_residual_kg<=diffusion_controls().exchange.mass_tolerance_kg);
        assert!(report.flow.endpoint_energy_change_j.abs()<=diffusion_controls().exchange.energy_tolerance_j);
        assert!(report.total_energy_residual_j.abs()<=settings.total_energy_tolerance_j);
        assert!(matches!(next.chambers()[&30].audit_last_window().unwrap(),
            ModuleWindowAudit::Unavailable {..}));
    }));
}

#[test]
fn adaptive_work_exhaustion_and_late_feedback_failure_publish_nothing() {
    fixture([0.2,0.2],150.0,|assembly,thermal,circuit|with_cx(|cx,gate| {
        let card = pair_card(0.002); let context = QueryPoint::new();
        let pairs = [binding(&card,&context)]; let links = [face(&pairs)];
        let prepared = assembly.prepare(cx,&thermal).unwrap();
        let mut state = MaterialExchangeState::from_parts(thermal,circuit,mixtures());
        let drives = drives(state.chambers());
        let saved_gas = state.chambers().clone();
        let saved_cells = state.thermal().cells().clone();
        let saved_circuit = state.circuit().clone();
        let limited_time = DiffusionNetworkAdvanceControls { maximum_evaluations: 3,
            ..adaptive_controls() };
        assert!(matches!(state.advance_with_adaptive_diffusion(&prepared,cx,DT,&[],&drives,&[],
            &gas_fixture::patches(),&links,diffusion_controls(),limited_time,
            GasDiffusionDriving::IsobaricMolar,controls(),&[0.0;2]),
            Err(MaterialFlowError::AdaptiveDiffusion(DiffusionNetworkAdvanceError::WorkLimit {..}))));
        let mut late = controls(); late.network.maximum_circuit_trials = 1;
        assert!(matches!(state.advance_with_adaptive_diffusion(&prepared,cx,DT,&[],&drives,&[],
            &gas_fixture::patches(),&links,diffusion_controls(),adaptive_controls(),
            GasDiffusionDriving::IsobaricMolar,late,&[0.0;2]),Err(MaterialFlowError::Exchange(_))));
        gate.request();
        assert!(matches!(state.advance_with_adaptive_diffusion(&prepared,cx,DT,&[],&drives,&[],
            &gas_fixture::patches(),&links,diffusion_controls(),adaptive_controls(),
            GasDiffusionDriving::IsobaricMolar,controls(),&[0.0;2]),Err(MaterialFlowError::Cancelled)));
        assert_eq!(state.chambers(),&saved_gas);
        assert_eq!(state.circuit(),&saved_circuit);
        assert_eq!(state.thermal().cells(),&saved_cells);
        assert_eq!(state.time_s(),0.0);
    }));
}
