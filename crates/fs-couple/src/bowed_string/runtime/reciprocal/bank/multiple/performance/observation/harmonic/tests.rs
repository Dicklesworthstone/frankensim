use super::*;

#[test]
fn g1_receiver_projection_retains_complex_phase_body_cancellation_and_physical_mouth_shares() {
    // Two body modes followed by one neck. Their indices are original flattened
    // coordinates, including unrelated strings and subsequent internal modes.
    let map = flatten_coordinates(&[2,3,2,4],&[
        (1,0,1.0),(1,2,1.0),(2,1,0.125),(2,1,0.375),
    ]).unwrap();
    assert_eq!(map,vec![(2,1.0),(4,1.0),(6,0.125),(6,0.375)]);
    let mut velocity=vec![C64::ZERO;11];
    velocity[2]=C64::new(3.0,-2.0);
    velocity[4]=-velocity[2];
    velocity[6]=C64::new(2.0,4.0);
    let body=[C64::new(4.0,-1.0),C64::new(-2.0,0.5),C64::from_re(7.0),C64::new(1.0,2.0)];
    let mouth=[C64::new(2.0,3.0),C64::new(-1.0,0.5),C64::new(0.0,-2.0),C64::from_re(4.0)];
    let rows:Vec<_>=body.into_iter().chain(body).chain(mouth).chain(mouth).collect();
    let observed=project_field(&rows,&map,&velocity).unwrap();
    // Structural fields cancel. The two apertures carry shares of the SAME
    // volume rate (1+2i), with no second neck and no complex conjugation.
    assert_eq!(observed,[C64::new(-4.0,7.0),C64::new(-2.0,-1.5),
        C64::new(4.0,-2.0),C64::new(4.0,8.0)]);
    let scale=C64::new(-0.5,0.25);
    let changed:Vec<_>=velocity.iter().map(|&v|v*scale).collect();
    let rotated=project_field(&rows,&map,&changed).unwrap();
    for (actual,old) in rotated.into_iter().zip(observed) {
        assert_eq!(actual,old*scale,"force phase must survive source-to-receiver composition");
    }
    // Permuting complete source rows with their physical addresses leaves the
    // field unchanged; projecting magnitudes would fail the cancellation above.
    let permutation=[2,0,3,1];
    let shuffled:Vec<_>=permutation.iter().flat_map(|&i|rows[4*i..4*i+4].iter().copied()).collect();
    let reordered:Vec<_>=permutation.iter().map(|&i|map[i]).collect();
    assert_eq!(project_field(&shuffled,&reordered,&velocity).unwrap(),observed);
}

#[test]
fn g0_receiver_refuses_missing_physical_inputs_and_reserves_observation_before_mechanics() {
    assert!(flatten_coordinates(&[2,3],&[(1,3,1.0)]).is_err());
    assert!(flatten_coordinates(&[2,3],&[(2,0,1.0)]).is_err());
    assert!(flatten_coordinates(&[2,3],&[(0,0,f64::NAN)]).is_err());
    assert!(flatten_coordinates(&[usize::MAX,1],&[(0,0,1.0)]).is_err());
    let map=[(1,0.5)];
    let rows=[C64::from_re(1.0);4];
    assert!(project_field(&rows[..3],&map,&[C64::ZERO;2]).is_err());
    assert!(project_field(&rows,&map,&[C64::ZERO]).is_err());
    assert!(project_field(&rows,&map,&[C64::ZERO,C64::new(f64::NAN,0.0)]).is_err());
    assert_eq!(remaining_mechanics_work(1500,1000,3,4).unwrap(),20);
    assert!(remaining_mechanics_work(1480,1000,3,4).is_err());
    assert!(remaining_mechanics_work(1000,1000,3,4).is_err());
    assert!(remaining_mechanics_work(usize::MAX,u64::MAX,1,1).is_err());
    assert!(remaining_mechanics_work(usize::MAX,1,usize::MAX,2).is_err());
    assert!(remaining_mechanics_work(1000,1,0,1).is_err());
    assert!(remaining_mechanics_work(1000,1,1,0).is_err());
}

#[test]
fn g0_skin_observation_budget_counts_every_actual_boundary_factorization() {
    let capacity=helmholtz::MAX_RADIATION_FIELDS_PER_BATCH;
    let one=one_way::frequency_terms(80,capacity).unwrap();
    let tail=one_way::frequency_terms(80,7).unwrap();
    assert_eq!(one_way::frequency_terms(80,2*capacity+7).unwrap(),2*one+tail);
    let old_single_factor=80*80*(80+2*capacity+7+1) as u64;
    assert!(2*one+tail>old_single_factor);
    assert!(one_way::frequency_terms(usize::MAX,capacity).is_none());
}
