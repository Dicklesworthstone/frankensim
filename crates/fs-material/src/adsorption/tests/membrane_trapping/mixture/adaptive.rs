//! Actual source/gas histories, independent temporal limits and atomic refusal.
use super::*;
use crate::thermochemical_gas::energy::parcel::transfer::species::network::transport::membrane::transient::{
    mixture::{adaptive::MembraneMixtureAdvanceControls, competitive::adaptive::*},
    slab::adaptive::MembraneDiffusionAdvanceError,
};

const BYTES: usize = 4 * 1024 * 1024;
fn controls(step: f64, pressure: f64, amount: f64) -> MembraneMixtureAdvanceControls {
    MembraneMixtureAdvanceControls { initial_step_s: step, minimum_step_s: 1e-6, maximum_step_s: step,
        maximum_evaluations: 20000, maximum_total_cells: 100,
        absolute_pressure_tolerance_pa: pressure, absolute_amount_tolerance_mol: amount,
        relative_tolerance: 0.0, species_balance_tolerance_mol: 1e-13,
        mass_balance_tolerance_kg: 1e-14, gas_energy_tolerance_j: 1e-9 }
}
fn audit_window(initial: &IsothermalCompetitiveGasMembrane<'_>, result: &IsothermalCompetitiveGasAdvanceCandidate<'_>) {
    let report = result.report();
    let next = result.state();
    assert_eq!(next.controls(), initial.controls());
    assert_eq!(next.gas_controls(), initial.gas_controls());
    assert!(report.total_mass_residual_kg.abs() <= 1e-14);
    assert!(report.direct_total_mass_residual_kg.abs() <= 1e-14);
    assert_eq!(report.source_relative_gas_energy_residual_j, [0.0; 2]);
    for (s,row) in report.species.iter().enumerate() {
        assert!(row.closed_species_residual_mol.abs() <= 1e-13);
        assert!(row.film_boundary_residual_mol.abs() <= 1e-13);
        for side in 0..2 {
            assert!(row.gas_boundary_residuals_mol[side].abs() <= 1e-13);
            assert!(row.gross_amount_into_gas_mol[side] >= row.amount_into_gas_mol[side].max(0.0));
            assert!(row.gross_amount_out_of_gas_mol[side] >= (-row.amount_into_gas_mol[side]).max(0.0));
            assert_eq!(next.finite_state().reservoir_partial_pressures_pa()[s][side],
                next.gases()[side].species_partial_pressure_pa(&row.species).unwrap());
        }
        assert_eq!(next.source_support(&row.species).unwrap().layers[0][0].card_identity(),
            initial.source_support(&row.species).unwrap().layers[0][0].card_identity());
    }
    for side in 0..2 {
        assert_eq!(next.gases()[side].state().state().temperature, T0);
        assert_eq!(&next.gases()[side].model().evaluate(next.gases()[side].point()).unwrap(), next.gases()[side].state());
        assert_eq!(next.gases()[side].model().thermodynamics().components(), initial.gases()[side].model().thermodynamics().components());
    }
}

fn reference_catalog(reference: f64) -> ThermochemicalGasModel {
    ThermochemicalGasModel::new(FrozenIdealGasMixtureModelV1::new(vec![
        species_with_energy_reference("A",0.028,3.5,0.0,reference),
        species_with_energy_reference("B",0.004,2.5,0.0,reference)],
        Composition::new(CompositionBasis::MoleFraction,vec![0.5,0.5]).unwrap()).unwrap(),
        vec![transport("A",1.7e-5,0.025),transport("B",2e-5,0.15)],None).unwrap()
}

#[test]
fn g1_actual_mixture_refines_to_finite_bath_eigenmodes_without_formation_reference_drift() {
    let card = membrane_card();
    let context = point(T0,PRESSURE);
    let duration = 50.0_f64;
    // Equal original C=V/(R*T) gas capacities and S*L*A film capacity.
    let capacity = 2e-5*1e-4*0.01;
    let volume = capacity * (R*T0);
    let rate = 2.0*1e-10/(1e-4*1e-4);
    let slow = (-rate*duration).exp();
    let fast = (-3.0*rate*duration).exp();
    let mean = 200_000.0/3.0;
    let exact = [mean-25_000.0*slow+(25_000.0/3.0)*fast,
        mean-(50_000.0/3.0)*fast,mean+25_000.0*slow+(25_000.0/3.0)*fast];
    let mut physical = None;
    for reference in [0.0,1e24] {
        let catalog = reference_catalog(reference);
        let gases = [PRESSURE,2.0*PRESSURE].map(|p| ClosedGasParcel::new(catalog.clone(),point(T0,p),volume).unwrap());
        let sources: Vec<_> = ["A","B"].into_iter().map(|name| MembraneMixtureChannelSource {
            layers: vec![MembraneDiffusionLayerSource { thickness_m: 1e-4, cells: 1,
                material: MembranePermeabilityMaterial { species: id(name),material:&card,context:&context,
                    selection:MaterialPropertySelection::SingleClaimOnly }}],
            interfaces:vec![],initial_partial_pressures_pa:vec![50_000.0],
        }).collect();
        let state = IsothermalGasMembraneMixture::new(T0,0.01,&sources,gases,gas_solve(),||false).unwrap()
            .with_competitive_sites(&[None],&[(id("A"),vec![0.0]),(id("B"),vec![0.0])],solve(),||false).unwrap();
        let run = |pressure,amount| state.advance_candidate(duration,controls(25.0,pressure,amount),
            &mut state.adaptive_workspace(BYTES).unwrap(),||false).unwrap();
        let fine = run(10.0,2e-10);
        audit_window(&state,&fine);
        let values = |c: &IsothermalCompetitiveGasAdvanceCandidate<'_>| {
            [c.state().gases()[0].species_partial_pressure_pa(&id("A")).unwrap(),
                c.state().membrane().mobile(0).unwrap().partial_pressures_pa()[0],
                c.state().gases()[1].species_partial_pressure_pa(&id("A")).unwrap()]
        };
        let error = |c:&IsothermalCompetitiveGasAdvanceCandidate<'_>| values(c).into_iter().zip(exact)
            .map(|(a,b)| (a-b).abs()).fold(0.0_f64,f64::max);
        if reference == 0.0 {
            let coarse = run(1000.0,2e-8);
            assert!(error(&fine) < 0.2*error(&coarse),"{} -> {}",error(&coarse),error(&fine));
        }
        assert_eq!(fine.state().elapsed_s(),duration);
        assert!(fine.report().diffusion.species[0].1.transport.rejected_intervals > 0);
        let signature=(values(&fine),fine.report().species.iter().map(|r|r.amount_into_gas_mol).collect::<Vec<_>>(),
            fine.report().thermostat_heat_into_gas_j);
        if let Some(expected)=&physical {assert_eq!(&signature,expected);}else{physical=Some(signature);}
        if reference > 0.0 {assert!(fine.report().enthalpy_into_gas_j[0].abs()>1e10);}
    }
}

#[test]
fn g3_competitive_displacement_retains_actual_gases_sources_and_accepted_only_accounts() {
    let catalog=model(vec![0.5,0.5],0.0);
    let card=membrane_card();let context=point(T0,PRESSURE);
    let source=sites(&catalog,[0.0,1e7]);
    let state=mixture(&card,&context,gases(&catalog,PRESSURE),[50_000.0;2],false)
        .with_competitive_sites(&[Some(source)],&[(id("A"),vec![2e-7;8]),(id("B"),vec![0.0;8])],solve(),||false).unwrap();
    let saved=state.finite_state().clone();let saved_gases=state.gases().clone();
    let c=controls(0.5,200.0,2e-9);
    let mut work=state.adaptive_workspace(BYTES).unwrap();
    let mut polls=0;
    let next=state.advance_candidate(1.0,c,&mut work,||{polls+=1;false}).unwrap();
    audit_window(&state,&next);
    assert!(next.report().diffusion.species[0].1.trapped_change_mol<0.0);
    assert!(next.report().diffusion.species[1].1.trapped_change_mol>0.0);
    assert!(next.report().species[0].actual_gas_amount_changes_mol.iter().sum::<f64>()>0.0);
    assert!(next.report().species[1].actual_gas_amount_changes_mol.iter().sum::<f64>()<0.0);
    for name in ["A","B"] {
        let a=state.trap_source_support(&id(name)).unwrap()[0].as_ref().unwrap();
        let b=next.state().trap_source_support(&id(name)).unwrap()[0].as_ref().unwrap();
        assert_eq!(a[0].source.card_identity(),b[0].source.card_identity());
    }
    // Exhaust work shortly before the known successful trajectory's endpoint.
    let mut limited = c;
    limited.maximum_evaluations = next.report().diffusion.species[0].1.transport.evaluations - 3;
    let limit = limited.maximum_evaluations;
    let refusal = state.advance_candidate(1.0, limited, &mut work, || false).unwrap_err();
    assert!(matches!(refusal,
        IsothermalCompetitiveGasAdvanceError::Diffusion(MembraneDiffusionAdvanceError::WorkLimit {time_s,evaluations}) if time_s>0.0 && evaluations==limit), "{refusal:?}");
    let mut seen=0;
    assert!(matches!(state.advance_candidate(1.0,c,&mut work,||{seen+=1;seen>=polls}),Err(IsothermalCompetitiveGasAdvanceError::Cancelled)));
    assert_eq!(state.finite_state(),&saved);assert_eq!(state.gases(),&saved_gases);
    let replay=state.advance_candidate(1.0,c,&mut work,||false).unwrap();
    assert_eq!(next.state().finite_state(),replay.state().finite_state());
    assert_eq!(next.state().gases(),replay.state().gases());assert_eq!(next.report(),replay.report());
}

#[test]
fn g3_every_adaptive_probe_keeps_original_source_domain_and_native_work_refusals() {
    let catalog=model(vec![0.5,0.5],0.0);
    let card=membrane_card();let context=point(T0,PRESSURE);
    let state=mixture(&card,&context,gases(&catalog,PRESSURE),[50_000.0;2],false)
        .with_competitive_sites(&[Some(sites(&catalog,[0.0,PRESSURE]))],
            &[(id("A"),vec![2.5e-7;8]),(id("B"),vec![0.0;8])],solve(),||false).unwrap();
    let saved=state.finite_state().clone();let saved_gases=state.gases().clone();
    assert!(matches!(state.advance_candidate(10.0,controls(10.0,1e6,1.0),
        &mut state.adaptive_workspace(BYTES).unwrap(),||false),
        Err(IsothermalCompetitiveGasAdvanceError::Native(IsothermalCompetitiveGasError::Source(_)))));
    assert_eq!(state.finite_state(),&saved);assert_eq!(state.gases(),&saved_gases);
    let limited=mixture(&card,&context,gases(&catalog,PRESSURE),[50_000.0;2],false)
        .with_competitive_sites(&[Some(sites(&catalog,[0.0,1e7]))],
            &[(id("A"),vec![2e-7;8]),(id("B"),vec![0.0;8])],
            CompetitiveMembraneSolve{maximum_linear_sweeps:1,..solve()},||false).unwrap();
    assert!(matches!(limited.advance_candidate(1.0,controls(1.0,1e6,1.0),
        &mut limited.adaptive_workspace(BYTES).unwrap(),||false),
        Err(IsothermalCompetitiveGasAdvanceError::Native(IsothermalCompetitiveGasError::Competitive(_)))));
    assert_eq!(limited.elapsed_s(),0.0);
}
