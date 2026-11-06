//! Finite-step use of the source-conjugate interior flux and native entropy gate.
use super::*;

fn fixture(formation: f64, absent: bool, one_sided_zero: bool) -> MixtureEuler1d {
    let mut models = vec![
        source("A",R,3.5,0.0,formation),
        source("B",R,6.5,0.0,-2.0*formation),
    ];
    if absent { models.push(source("C",R,4.0,0.0,3.0*formation)); }
    let g = MixtureGas::new(
        bank(models),
        CaloricDomain::new(Temperature::new(0.001),Temperature::new(50.0)),
        ||false,
    ).unwrap();
    field(g,2,Boundaries::Ends {
        left:Boundary::Transmissive,right:Boundary::Transmissive,
    },|x| {
        let a = if x<0.5 { if one_sided_zero {1.0} else {0.1} } else {0.9};
        let mut rho = vec![a,1.0-a];
        if absent { rho.push(0.0); }
        (rho,0.0,if x<0.5 {1.0} else {10.0})
    })
}

fn independent_entropy(flow: &MixtureEuler1d) -> f64 {
    (0..flow.cell_count()).map(|i| {
        let q = flow.cell(i).unwrap();
        let t = flow.gas().primitive(q,||false).unwrap().temperature;
        flow.mesh().volumes()[i]*q.partial_densities.iter().zip([2.5,5.5,3.0])
            .filter(|(rho,_)|**rho>0.0)
            .map(|(&rho,cv)|rho*(cv*t.ln()-rho.ln())).sum::<f64>()
    }).sum()
}

#[test]
fn actual_ssp_step_repairs_the_negative_spatial_fixture_without_energy_injection() {
    let mut reference = None;
    for formation in [0.0,1e18] {
        let mut flow = fixture(formation,false,false);
        let before = independent_entropy(&flow);
        let initial = flow.states.clone();
        let mut op = operator(&flow,0.0,0.0,0.0,0.01);
        let (report,entropy) = op.advance_entropy_consistent_for(
            &mut flow,Time::new(1e-4),
            DissipationControl { maximum_step_s:1e-4,..Default::default() },
            MolecularEntropyControl { absolute_tolerance_j_k:1e-14,..Default::default() },
            ||false,
        ).unwrap();
        assert_eq!(report.accepted_steps,1);
        assert_eq!(report.rejected_steps,0);
        assert_ne!(flow.states,initial);
        assert_eq!(flow.time().value(),0.0);
        assert_eq!(report.boundary_change.energy,0.0);
        assert!(entropy.entropy_change_j_k>0.0);
        let measured = independent_entropy(&flow)-before;
        assert!((entropy.source_entropy_change_j_k-measured).abs()<2e-13);
        // Original constant-Cp sources give the independently derived divided
        // h: Cp*TL*TR*ln(TR/TL)/(TR-TL)-Cv*Tmin. The resulting face has
        // j_A=-0.021972245773362195 and referenced E=0.16857730185383066.
        // Its initial total entropy rate is +0.06758914180550454 W/K.
        assert!((measured/1e-4-0.06758914180550454).abs()<3e-5);
        // Independently applying SSP-RK2 to those analytic face equations.
        assert!((measured-6.758763986258742e-6).abs()<2e-12);
        assert!(flow.cell(0).unwrap().partial_densities[0]>0.1);
        assert!(flow.cell(1).unwrap().partial_densities[0]<0.9);
        if let Some((states,account)) = &reference {
            assert_eq!(&flow.states,states);
            assert_eq!(&entropy,account);
        } else {
            assert!(report.absolute_energy_residual_j.abs()<2e-13);
            reference=Some((flow.states.clone(),entropy));
        }
    }
}

#[test]
fn common_absence_is_retained_and_one_sided_zero_refuses_atomically() {
    let mut flow = fixture(0.0,true,false);
    let mut op = operator(&flow,0.0,0.0,0.0,0.01);
    let (report,entropy) = op.advance_entropy_consistent_for(
        &mut flow,Time::new(1e-4),DissipationControl::default(),
        MolecularEntropyControl::default(),||false,
    ).unwrap();
    assert!(entropy.entropy_change_j_k>0.0);
    assert_eq!(report.after.species_masses[2],0.0);
    for i in 0..2 { assert_eq!(flow.cell(i).unwrap().partial_densities[2],0.0); }

    let mut refused = fixture(0.0,false,true);
    let before = refused.states.clone();
    let mut op = operator(&refused,0.0,0.0,0.0,0.01);
    assert!(matches!(op.advance_entropy_consistent_for(
        &mut refused,Time::new(1e-4),DissipationControl::default(),
        MolecularEntropyControl::default(),||false,
    ),Err(DissipationError::Invalid(_))));
    assert_eq!(refused.states,before);
    assert_eq!(refused.time().value(),0.0);
}
