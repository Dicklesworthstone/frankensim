//! Full split accuracy, physical feedback and accepted-only source ledgers.
use super::*;
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::material::state::flow::diffusion::coupled::*;

fn loose() -> CoupledDiffusionControls {
    CoupledDiffusionControls {
        absolute: CoupledDiffusionScales { species_mass_kg: 1.0, gas_temperature_k: 1000.0,
            gas_pressure_pa: 1e6, gas_volume_m3: 1.0, piston_velocity_m_s: 100.0,
            cell_energy_j: 1e6, charge_c: 1.0, flux_wb: 1.0 },
        relative_tolerance: 0.0, maximum_step_s: DT, minimum_step_s: DT / 1024.0,
        maximum_evaluations: 1000,
    }
}
fn species_mass(gas: &GasPiston, species: &str) -> f64 {
    let n = gas.gas().mass_kg() / gas.gas().state().thermodynamics().properties().molar_mass().value();
    gas.gas().model().thermodynamics().components().iter()
        .zip(gas.gas().model().thermodynamics().mole_composition().fractions())
        .find(|(s, _)| s.species().as_str() == species)
        .map_or(0.0, |(s, x)| (n * x) * s.molar_mass().value())
}
fn half_drives(drives: &[ChamberDrive]) -> Vec<ChamberDrive> {
    drives.iter().map(|d| ChamberDrive { heat_j: 0.5 * d.heat_j, ..*d }).collect()
}

#[test]
fn coupled_two_half_steps_match_native_owners_and_complete_reaction_ledgers() {
    fixture([0.2,0.2],150.0,|assembly,thermal,circuit|with_cx(|cx,_| {
        let law = reaction(50.0); let card = pair_card(2e-5); let point = QueryPoint::new();
        let pairs = [binding(&card,&point)]; let links = [face(&pairs)];
        let gases = [0.75,0.25].into_iter().enumerate().map(|(i,x)| {
            let model = law.catalog().with_composition(Composition::new(
                CompositionBasis::MoleFraction,vec![x,1.0-x]).unwrap()).unwrap();
            (10+10*i as u64,parcel_piston(model,100000.0))
        }).collect();
        let prepared = assembly.prepare(cx,&thermal).unwrap();
        let state = MaterialExchangeState::from_parts(thermal,circuit,gases);
        let mut drives = drives(state.chambers());
        drives[0].heat_j = 0.01; drives[1].heat_j = 0.02;
        drives[0].external_force_n = 2.0;
        let reactions = [ChamberReaction { chamber:10,reaction:&law,controls:chemical() }];
        let half = half_drives(&drives);
        let (mid, first) = state.diffusion_step_candidate(&prepared,cx,DT/2.0,&[],&half,&reactions,
            &gas_fixture::patches(),&links,diffusion_controls(),controls(),&[0.0;2]).unwrap();
        let (expected, second) = mid.diffusion_step_candidate(&prepared,cx,DT,&[],&half,&reactions,
            &gas_fixture::patches(),&links,diffusion_controls(),controls(),&[0.0;2]).unwrap();
        let (actual, report) = state.coupled_diffusion_candidate(&prepared,cx,DT,&[],&drives,&reactions,
            &gas_fixture::patches(),&links,diffusion_controls(),GasDiffusionDriving::IsobaricMolar,
            controls(),&[0.0;2],CoupledDiffusionControls { maximum_evaluations:3,..loose() }).unwrap();
        assert_eq!(actual.chambers(),expected.chambers());
        assert_eq!(actual.circuit(),expected.circuit());
        assert_eq!(actual.thermal().cells(),expected.thermal().cells());
        assert_eq!(report.steps[0].advance.flow,first.flow);
        assert_eq!(report.steps[1].advance.flow,second.flow);
        assert_eq!(report.steps[0].advance.exchange.report.source.reactions,first.exchange.report.source.reactions);
        assert_eq!(report.steps[1].advance.exchange.report.source.reactions,second.exchange.report.source.reactions);
        assert_eq!(report.evaluations,3); assert_eq!(report.accepted_intervals,1);
        assert_eq!(report.steps.len(),2); assert_eq!(report.rejected_intervals,0);
        assert!(report.estimated_error.maximum()<=1.0);
        assert!((report.balance.prescribed_gas_heat_j-0.03).abs()<1e-16);
        let extent = first.exchange.report.source.reactions[0].reaction.extent_mol
            + second.exchange.report.source.reactions[0].reaction.extent_mol;
        assert!((report.balance.reaction_extents[0].extent_mol-extent).abs()<1e-16);
        assert!(report.balance.species.iter().any(|s|s.reaction_change_kg!=0.0));
        assert!(report.balance.species.iter().any(|s|s.diffusive_change_kg!=0.0));
        assert!(report.balance.maximum_species_mass_residual_kg<=diffusion_controls().exchange.mass_tolerance_kg);
        assert!(report.balance.total_energy_residual_j.abs()<=controls().total_energy_tolerance_j);
        assert!(report.balance.maximum_cell_residual_bound_j<=controls().network.cell_energy_tolerance_j);
        let cell_bound = first.exchange.report.maximum_cell_residual_j + second.exchange.report.maximum_cell_residual_j;
        assert!((report.balance.maximum_cell_residual_bound_j-cell_bound).abs()<1e-20);
        for (&id, gas) in actual.chambers() {
            assert_eq!(gas.last_window_duration_s(),Some(DT));
            assert!(matches!(gas.audit_last_window().unwrap(),ModuleWindowAudit::Unavailable {..}));
            let account = report.balance.chambers.iter().find(|r|r.chamber==id).unwrap();
            assert!(account.momentum_residual_n_s.abs()<=drives[0].piston_solve.momentum_tolerance_n_s);
        }
        assert_eq!(state.time_s(),0.0);
    }));
}

fn differences(a: &MaterialExchangeState<'_>, b: &MaterialExchangeState<'_>) -> [f64;8] {
    let mut error = [0.0_f64;8];
    for (&id,left) in a.chambers() {
        let right = &b.chambers()[&id];
        for species in ["A","B"] {
            error[0] = error[0].max((species_mass(left,species)-species_mass(right,species)).abs());
        }
        let x = left.gas().state().state(); let y = right.gas().state().state();
        error[1] = error[1].max((x.temperature-y.temperature).abs());
        error[2] = error[2].max((x.pressure-y.pressure).abs());
        error[3] = error[3].max((left.gas().volume_m3()-right.gas().volume_m3()).abs());
        error[4] = error[4].max((left.velocity_m_s()-right.velocity_m_s()).abs());
    }
    for (id,left) in a.thermal().cells() {
        error[5] = error[5].max(left.mass_kg*(left.specific_h_j_kg-b.thermal().cells()[id].specific_h_j_kg).abs());
    }
    for (left,right) in a.circuit().charges_c().iter().zip(b.circuit().charges_c()) {
        error[6] = error[6].max((left-right).abs());
    }
    for (left,right) in a.circuit().fluxes_wb().iter().zip(b.circuit().fluxes_wb()) {
        error[7] = error[7].max((left-right).abs());
    }
    error
}
fn select_scale(scales: &mut CoupledDiffusionScales, i: usize, value: f64) {
    match i {
        0 => scales.species_mass_kg=value, 1 => scales.gas_temperature_k=value,
        2 => scales.gas_pressure_pa=value, 3 => scales.gas_volume_m3=value,
        4 => scales.piston_velocity_m_s=value, 5 => scales.cell_energy_j=value,
        6 => scales.charge_c=value, 7 => scales.flux_wb=value, _ => unreachable!(),
    }
}

#[test]
fn coupled_estimate_cannot_ignore_thermal_electrical_or_piston_coordinates() {
    fixture([0.8,0.8],150.0,|assembly,thermal,circuit|with_cx(|cx,_| {
        let card = pair_card(2e-5); let point = QueryPoint::new();
        let pairs = [binding(&card,&point)]; let links = [face(&pairs)];
        // Keep original resistor/contact branches; add a real native inductive
        // state so both electrical storage coordinate families are exercised.
        let mut graph = circuit.graph().clone();
        graph.branches.push((2,0,Branch::Inductor { henries:0.01 }));
        let circuit = MaterialCircuit::new_with_contacts(graph,circuit.resistors().to_vec(),
            circuit.contacts().to_vec(),vec![0.001],vec![1.0],&[]).unwrap();
        let prepared = assembly.prepare(cx,&thermal).unwrap();
        let state = MaterialExchangeState::from_parts(thermal,circuit,mixtures());
        let drives = drives(state.chambers());
        let (coarse,_) = state.diffusion_step_candidate(&prepared,cx,DT,&[],&drives,&[],
            &gas_fixture::patches(),&links,diffusion_controls(),controls(),&[0.0;2]).unwrap();
        let (mid,_) = state.diffusion_step_candidate(&prepared,cx,DT/2.0,&[],&drives,&[],
            &gas_fixture::patches(),&links,diffusion_controls(),controls(),&[0.0;2]).unwrap();
        let (fine,_) = mid.diffusion_step_candidate(&prepared,cx,DT,&[],&drives,&[],
            &gas_fixture::patches(),&links,diffusion_controls(),controls(),&[0.0;2]).unwrap();
        let delta = differences(&coarse,&fine);
        for (i,difference) in delta.into_iter().enumerate() {
            assert!(difference>0.0,"fixture must excite coordinate {i}: {delta:?}");
            let mut policy = CoupledDiffusionControls { minimum_step_s:DT,maximum_evaluations:3,..loose() };
            select_scale(&mut policy.absolute,i,difference/4.0);
            let result = state.coupled_diffusion_candidate(&prepared,cx,DT,&[],&drives,&[],
                &gas_fixture::patches(),&links,diffusion_controls(),GasDiffusionDriving::IsobaricMolar,
                controls(),&[0.0;2],policy);
            let Err(CoupledDiffusionError::Accuracy { estimate,.. }) = result else {
                panic!("coordinate {i} did not reject an inaccurate complete assembly");
            };
            assert!(estimate.maximum()>3.9,"coordinate {i}: {estimate:?}");
        }
        let mut policy = loose();
        policy.absolute.charge_c = delta[6]/8.0;
        let (accepted,report) = state.coupled_diffusion_candidate(&prepared,cx,DT,&[],&drives,&[],
            &gas_fixture::patches(),&links,diffusion_controls(),GasDiffusionDriving::IsobaricMolar,
            controls(),&[0.0;2],policy).unwrap();
        assert!(report.rejected_intervals>0 && report.accepted_intervals>1);
        assert!(report.estimated_error.maximum()<=1.0);
        assert_eq!(accepted.time_s(),DT);
        assert!(report.balance.total_energy_residual_j.abs()<=controls().total_energy_tolerance_j);
    }));
}

#[test]
fn coupled_refinement_retries_withdrawals_and_converges_to_binary_relaxation() {
    fixture([0.2,0.2],150.0,|assembly,thermal,circuit|with_cx(|cx,_| {
        let card = pair_card(0.002); let point = QueryPoint::new();
        let pairs = [binding(&card,&point)]; let links = [face(&pairs)];
        let prepared = assembly.prepare(cx,&thermal).unwrap();
        let chambers = mixtures().into_iter().map(|(id,gas)|(id,gas.with_fixed_support().unwrap())).collect();
        let state = MaterialExchangeState::from_parts(thermal,circuit,chambers);
        let drives = drives(state.chambers());
        assert!(matches!(state.diffusion_step_candidate(&prepared,cx,DT,&[],&drives,&[],&[],
            &links,diffusion_controls(),controls(),&[0.0;2]),Err(MaterialFlowError::Diffusion(_))));
        // Constant-volume binary/Fick eigenmode; the finite circuit and solids
        // still advance, with no gas-wall heat to change this independent oracle.
        let exact = 0.5*(-2.0_f64*0.002*0.1/(0.0001*0.001)*DT).exp();
        let mut errors = Vec::new(); let mut counts = Vec::new();
        for relative in [0.02,0.002] {
            let mut policy = loose(); policy.relative_tolerance=relative;
            policy.absolute.species_mass_kg=1e-12;
            let (next,report) = state.coupled_diffusion_candidate(&prepared,cx,DT,&[],&drives,&[],&[],
                &links,diffusion_controls(),GasDiffusionDriving::IsobaricMolar,controls(),&[0.0;2],policy).unwrap();
            let x = |id|next.chambers()[&id].gas().model().thermodynamics().mole_composition().fractions()[0];
            errors.push(((x(10)-x(20))-exact).abs()); counts.push(report.steps.len());
            assert!(report.rejected_intervals>0 && report.accepted_intervals>1);
            assert!(report.estimated_error.maximum()<=1.0);
            assert_eq!(report.steps.len(),2*report.accepted_intervals);
            assert!(report.evaluations>report.steps.len());
            assert!(report.balance.maximum_species_mass_residual_kg<=diffusion_controls().exchange.mass_tolerance_kg);
            assert!(report.balance.total_energy_residual_j.abs()<=controls().total_energy_tolerance_j);
            assert_eq!(next.time_s(),DT);
            let gross: f64 = report.balance.transfers.iter().map(|t|t.mass_kg).sum();
            let accepted: f64 = report.steps.iter().flat_map(|s|&s.advance.flow.exchange.fluxes)
                .map(|f|f.transfer.mass_kg).sum();
            assert!((gross-accepted).abs()<1e-16);
        }
        assert!(errors[1]<0.3*errors[0],"coupled refinement errors {errors:?}");
        assert!(counts[1]>counts[0]); assert_eq!(state.time_s(),0.0);
    }));
}

#[test]
fn coupled_pressure_diffusion_samples_updated_live_feedback_on_the_second_half() {
    fixture([0.2,0.2],150.0,|assembly,thermal,circuit|with_cx(|cx,_| {
        let card = pair_card(0.002); let point = QueryPoint::new();
        let pairs = [binding(&card,&point)]; let links = [face(&pairs)];
        let model = mixtures()[&10].gas().model().clone();
        let gases = BTreeMap::from([(10,parcel_piston(model.clone(),100500.0)),(20,parcel_piston(model,99500.0))]);
        let prepared = assembly.prepare(cx,&thermal).unwrap();
        let state = MaterialExchangeState::from_parts(thermal,circuit,gases);
        let mut drives = drives(state.chambers()); drives[0].external_force_n=10.0;
        let driving = GasDiffusionDriving::CompositionAndPressureMass { mass_flux_tolerance_kg_m2_s:1e-12 };
        let (mid,_) = state.diffusion_step_with_driving_candidate(&prepared,cx,DT/2.0,&[],&drives,&[],
            &gas_fixture::patches(),&links,diffusion_controls(),driving,controls(),&[0.0;2]).unwrap();
        let (expected,last) = mid.diffusion_step_with_driving_candidate(&prepared,cx,DT,&[],&drives,&[],
            &gas_fixture::patches(),&links,diffusion_controls(),driving,controls(),&[0.0;2]).unwrap();
        let (next,report) = state.coupled_diffusion_candidate(&prepared,cx,DT,&[],&drives,&[],
            &gas_fixture::patches(),&links,diffusion_controls(),driving,controls(),&[0.0;2],loose()).unwrap();
        assert_eq!(next.chambers(),expected.chambers());
        assert_eq!(report.steps[1].advance.flow,last.flow);
        let first = &report.steps[0].advance.flow.links[0].flux;
        let second = &report.steps[1].advance.flow.links[0].flux;
        assert_ne!(first.face_pressure_pa,second.face_pressure_pa);
        assert_ne!(state.chambers()[&10].gas().volume_m3(),mid.chambers()[&10].gas().volume_m3());
        assert_eq!(second.driving,driving);
        assert!(second.mass_frame_residual_kg_m2_s.abs()<=1e-12);
        assert!(!report.balance.transfers.is_empty());
        for (&id,gas) in state.chambers() {
            assert!((next.chambers()[&id].gas().mass_kg()-gas.gas().mass_kg()).abs()<1e-13);
        }
    }));
}

#[test]
fn coupled_failed_coarse_and_fine_work_sources_and_cancellation_are_atomic() {
    fixture([0.2,0.2],150.0,|assembly,thermal,circuit|with_cx(|cx,gate| {
        let card = pair_card(0.002); let point = QueryPoint::new();
        let pairs = [binding(&card,&point)]; let links = [face(&pairs)];
        let prepared = assembly.prepare(cx,&thermal).unwrap();
        let mut state = MaterialExchangeState::from_parts(thermal,circuit,mixtures());
        let drives = drives(state.chambers());
        let gas = state.chambers().clone(); let cells = state.thermal().cells().clone(); let circuit = state.circuit().clone();
        let limited = CoupledDiffusionControls { maximum_evaluations:3,..loose() };
        assert!(matches!(state.advance_with_coupled_diffusion(&prepared,cx,DT,&[],&drives,&[],
            &gas_fixture::patches(),&links,diffusion_controls(),GasDiffusionDriving::IsobaricMolar,
            controls(),&[0.0;2],limited),Err(CoupledDiffusionError::WorkLimit { evaluations:3,accepted_substeps:0,.. })));
        let missing = DiffusionNetworkLink { pairs:&[],..links[0].clone() };
        assert!(matches!(state.advance_with_coupled_diffusion(&prepared,cx,DT,&[],&drives,&[],
            &gas_fixture::patches(),&[missing],diffusion_controls(),GasDiffusionDriving::IsobaricMolar,
            controls(),&[0.0;2],loose()),Err(CoupledDiffusionError::Step { error:MaterialFlowError::Diffusion(_),.. })));
        let mut late = controls(); late.network.maximum_circuit_trials=1;
        assert!(matches!(state.advance_with_coupled_diffusion(&prepared,cx,DT/8.0,&[],&drives,&[],
            &gas_fixture::patches(),&links,diffusion_controls(),GasDiffusionDriving::IsobaricMolar,
            late,&[0.0;2],loose()),Err(CoupledDiffusionError::Step { error:MaterialFlowError::Exchange(_),.. })));
        gate.request();
        assert!(matches!(state.advance_with_coupled_diffusion(&prepared,cx,DT,&[],&drives,&[],
            &gas_fixture::patches(),&links,diffusion_controls(),GasDiffusionDriving::IsobaricMolar,
            controls(),&[0.0;2],loose()),Err(CoupledDiffusionError::Cancelled)));
        assert_eq!(state.time_s(),0.0); assert_eq!(state.chambers(),&gas);
        assert_eq!(state.thermal().cells(),&cells); assert_eq!(state.circuit(),&circuit);
    }));
}
