//! Native open molecular histories with actual heat and chemical bath entropy.
use super::*;

fn ends() -> Boundaries {
    Boundaries::Ends { left: Boundary::Transmissive, right: Boundary::Transmissive }
}

fn physical_entropy(flow: &MixtureEuler1d) -> f64 {
    // Original test source: M=R, cp/R=3.5, s0/R=3.5 ln T, p0=1.
    (0..flow.cell_count()).map(|i| {
        let c = flow.cell(i).unwrap();
        let t = flow.gas().primitive(c, || false).unwrap().temperature;
        flow.mesh().volumes()[i] * c.partial_densities.iter().filter(|&&rho| rho > 0.0)
            .map(|&rho| rho * (2.5 * t.ln() - rho.ln())).sum::<f64>()
    }).sum()
}

#[test]
fn stationary_fourier_profile_uses_each_baths_own_temperature() {
    let g = gas(1, 0.0, false);
    let mesh = Mesh1d::new(
        vec![0.0, 0.2, 0.5, 1.0], vec![0.7;4], vec![0.14,0.21,0.35],
    ).unwrap();
    let states = [2.1,2.35,2.75].iter().map(|&t| g.conserved(&[1.0],0.0,t).unwrap()).collect();
    let left = g.conserved(&[1.0],0.0,2.0).unwrap();
    let right = g.conserved(&[1.0],0.0,3.0).unwrap();
    let mut flow = MixtureEuler1d::new_on_mesh(g.clone(),mesh.clone(),states,ends(),||false).unwrap();
    let original = flow.states.clone();
    let original_states = [2.1,2.35,2.75].iter().map(|&t| g.conserved(&[1.0],0.0,t).unwrap()).collect();
    let mut native = MixtureEuler1d::new_on_mesh(g,mesh,original_states,ends(),||false).unwrap();
    let mut plain = operator(&native,0.0,0.0,0.05,0.0);
    let mut op = operator(&flow,0.0,0.0,0.05,0.0);
    let baths = [Some(left.as_cell()),Some(right.as_cell())];
    let expected = plain.advance_with_reservoirs_for(
        &mut native,Time::new(0.04),DissipationControl::default(),baths,||false,
    ).unwrap();
    let (thermal,entropy) = op.advance_with_reservoirs_and_entropy_for(
        &mut flow,Time::new(0.04),DissipationControl::default(),
        MolecularEntropyControl::default(),baths,||false,
    ).unwrap();
    assert_eq!(thermal,expected);
    assert_eq!(flow.states,native.states);
    assert_eq!(flow.time().value(),0.0);
    for (&a,&b) in flow.states.iter().zip(&original) { assert!((a-b).abs()<2e-12); }
    let heat = 0.05*0.7*0.04;
    assert!((entropy.boundary_entropy_input_j_k[0]+heat/2.0).abs()<2e-13);
    assert!((entropy.boundary_entropy_input_j_k[1]-heat/3.0).abs()<2e-13);
    assert!(entropy.source_entropy_change_j_k.abs()<2e-12);
    assert!((entropy.entropy_generation_j_k-heat*(1.0/2.0-1.0/3.0)).abs()<2e-12);
}

#[test]
fn chemical_baths_can_reduce_body_entropy_while_total_generation_is_positive() {
    let mut baseline = None;
    for formation in [0.0,1e18] {
        let g = gas(2,formation,false);
        let left = g.conserved(&[0.8,0.2],0.0,2.0).unwrap();
        let right = g.conserved(&[0.2,0.8],0.0,2.0).unwrap();
        let make = || field(g.clone(),8,ends(),|_|(vec![0.5,0.5],0.0,2.0));
        let mut native = make();
        let mut flow = make();
        let before = physical_entropy(&flow);
        let mut plain = operator(&native,0.0,0.0,0.0,0.03);
        let mut op = operator(&flow,0.0,0.0,0.0,0.03);
        let controls = DissipationControl { maximum_step_s:0.01,..Default::default() };
        let baths = [Some(left.as_cell()),Some(right.as_cell())];
        let expected = plain.advance_with_reservoirs_for(
            &mut native,Time::new(0.03),controls,baths,||false,
        ).unwrap();
        let (report,entropy) = op.advance_with_reservoirs_and_entropy_for(
            &mut flow,Time::new(0.03),controls,MolecularEntropyControl::default(),baths,||false,
        ).unwrap();
        assert_eq!(flow.states,native.states);
        assert_eq!(report,expected);
        let state_change = physical_entropy(&flow)-before;
        assert!(state_change<0.0);
        assert!(entropy.entropy_generation_j_k>0.0);
        assert!((entropy.source_entropy_change_j_k-state_change).abs()<2e-12);
        let incoming = [report.left_boundary_change,report.right_boundary_change];
        let densities: [[f64;2];2] = [[0.8,0.2],[0.2,0.8]];
        let mut input = 0.0;
        for end in 0..2 {
            // q_reduced=0; the independently known bath partial entropies carry
            // species exchange. Equal partial enthalpies make net energy zero.
            let expected: f64 = (0..2).map(|k| incoming[end].species_masses[k]
                *(2.5*2.0_f64.ln()-densities[end][k].ln())).sum();
            assert!((entropy.boundary_entropy_input_j_k[end]-expected).abs()<2e-12);
            input += expected;
        }
        assert!(input<0.0);
        assert!((entropy.source_entropy_generation_j_k-(state_change-input)).abs()<3e-12);
        assert!((entropy.source_entropy_generation_j_k-entropy.entropy_generation_j_k
            -entropy.reference_entropy_residual_j_k).abs()<2e-14);
        assert!(entropy.maximum_relative_species_residual<1e-12);
        if let Some((states,account)) = &baseline {
            assert_eq!(&flow.states,states);
            assert_eq!(&entropy,account);
        } else { baseline = Some((flow.states.clone(),entropy)); }
    }
}

#[test]
fn admitted_bath_can_introduce_an_initially_absent_species() {
    let g = gas(2,0.0,false);
    let bath = g.conserved(&[0.5,0.5],0.0,2.0).unwrap();
    let mut flow = field(g,4,ends(),|_|(vec![1.0,0.0],0.0,2.0));
    let mut op = operator(&flow,0.0,0.0,0.0,0.03);
    let (report,entropy) = op.advance_with_reservoirs_and_entropy_for(
        &mut flow,Time::new(0.01),DissipationControl::default(),
        MolecularEntropyControl::default(),[Some(bath.as_cell()),None],||false,
    ).unwrap();
    assert_eq!(report.before.species_masses[1],0.0);
    assert!(report.after.species_masses[1]>0.0);
    assert!(report.left_boundary_change.species_masses[1]>0.0);
    assert!(entropy.entropy_generation_j_k>0.0);
    assert!(entropy.maximum_relative_species_residual<1e-12);
}

#[test]
fn zero_activity_species_and_late_cancellation_refuse_without_partial_publication() {
    let g = gas(2,0.0,false);
    let pure = g.conserved(&[1.0,0.0],0.0,3.0).unwrap();
    let make = || field(g.clone(),4,ends(),|_|(vec![0.5,0.5],0.0,2.0));
    let mut flow = make();
    let before = flow.states.clone();
    let mut op = operator(&flow,0.0,0.0,0.01,0.03);
    assert!(op.advance_with_reservoirs_and_entropy_for(
        &mut flow,Time::new(0.01),
        DissipationControl { maximum_halvings:1,..Default::default() },
        MolecularEntropyControl::default(),[Some(pure.as_cell()),None],||false,
    ).is_err());
    assert_eq!(flow.states,before);
    // A pure bath is admissible when the absent species has no molecular flux.
    let mut heat_only = operator(&flow,0.0,0.0,0.01,0.0);
    let (_,entropy) = heat_only.advance_with_reservoirs_and_entropy_for(
        &mut flow,Time::new(0.01),DissipationControl::default(),
        MolecularEntropyControl::default(),[Some(pure.as_cell()),None],||false,
    ).unwrap();
    assert!(entropy.entropy_generation_j_k>0.0);

    let bath = g.conserved(&[0.8,0.2],0.0,2.0).unwrap();
    let controls = DissipationControl { maximum_step_s:0.005,..Default::default() };
    let mut reference = make();
    let mut reference_op = operator(&reference,0.0,0.0,0.0,0.03);
    let mut polls = 0;
    reference_op.advance_with_reservoirs_and_entropy_for(
        &mut reference,Time::new(0.02),controls,MolecularEntropyControl::default(),
        [Some(bath.as_cell()),None],||{polls+=1;false},
    ).unwrap();
    let mut flow = make();
    let before = flow.states.clone();
    let mut op = operator(&flow,0.0,0.0,0.0,0.03);
    let mut visits = 0;
    let result = op.advance_with_reservoirs_and_entropy_for(
        &mut flow,Time::new(0.02),controls,MolecularEntropyControl::default(),
        [Some(bath.as_cell()),None],||{visits+=1;visits>=polls},
    );
    assert_eq!(result,Err(DissipationError::Cancelled));
    assert_eq!(flow.states,before);
    assert_eq!(flow.time().value(),0.0);
    op.advance_with_reservoirs_and_entropy_for(
        &mut flow,Time::new(0.02),controls,MolecularEntropyControl::default(),
        [Some(bath.as_cell()),None],||false,
    ).unwrap();
    for (actual,expected) in flow.states.iter().zip(&reference.states) {
        assert!((actual-expected).abs()<1e-11);
    }
}
