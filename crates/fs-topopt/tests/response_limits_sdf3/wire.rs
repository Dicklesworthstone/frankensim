use super::*;

#[test]
fn bounded_bytes_rebuild_the_actual_search_and_preserve_exact_floating_bits() {
    let (mut study, mut f, _) = fixture(false); f.fill(0.0);
    let rows = lower(0.002);
    let cases = [ResponseLimitCase3 { force: &f, prescribed: Some(&motion), limits: &rows }];
    let rho = vec![0.7; study.cells()]; let options = MinimumMaterialOptions3::default();
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let mut session = MinimumMaterialStudy3::new(&mut study, &cases, &rho, options, &mut control).unwrap();
    session.run(2).unwrap();
    let checkpoint = session.checkpoint(); let expected = session.accepted().clone(); drop(session);
    let bytes = checkpoint.to_bytes(Default::default()).unwrap();
    let mut decoded = MinimumMaterialCheckpoint3::from_bytes(&bytes, Default::default()).unwrap();
    assert_eq!(bytes, decoded.to_bytes(Default::default()).unwrap());
    assert_eq!(checkpoint.optimizer, decoded.optimizer);
    assert_eq!(checkpoint.accepted, decoded.accepted);
    let restored = MinimumMaterialStudy3::restore(&mut study, &cases, &mut decoded, options, &mut control).unwrap();
    same(restored.accepted(), &expected);
    // A valid signed zero is preserved as bits rather than formatted/parsed text.
    let mut signed = checkpoint; signed.optimizer.multiplier = -0.0;
    let bytes = signed.to_bytes(Default::default()).unwrap();
    assert_eq!(MinimumMaterialCheckpoint3::from_bytes(&bytes, Default::default()).unwrap()
        .optimizer.multiplier.to_bits(), (-0.0_f64).to_bits());
}

#[test]
fn malformed_truncated_and_over_budget_payloads_never_return_a_checkpoint() {
    let (mut study, mut f, _) = fixture(false); f.fill(0.0); let rows = lower(0.002);
    let cases = [ResponseLimitCase3 { force: &f, prescribed: Some(&motion), limits: &rows }];
    let rho = vec![0.7; study.cells()];
    let mut callback = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(), &mut callback);
    let session = MinimumMaterialStudy3::new(&mut study, &cases, &rho, Default::default(), &mut c).unwrap();
    let saved = session.checkpoint(); let bytes = saved.to_bytes(Default::default()).unwrap();
    for end in [0, 7, 8, 16, bytes.len()-1] {
        assert!(MinimumMaterialCheckpoint3::from_bytes(&bytes[..end], Default::default()).is_err());
    }
    let mut bad = bytes.clone(); bad[7] = 255;
    assert!(MinimumMaterialCheckpoint3::from_bytes(&bad, Default::default()).is_err());
    let mut bad = bytes.clone(); bad.extend([0]);
    assert!(MinimumMaterialCheckpoint3::from_bytes(&bad, Default::default()).is_err());
    // First point-vector length: reject before allocating anything from this length.
    let mut bad = bytes.clone(); bad[8..16].copy_from_slice(&u64::MAX.to_le_bytes());
    assert!(MinimumMaterialCheckpoint3::from_bytes(&bad, Default::default()).is_err());
    let mut bad = bytes.clone(); bad[16..24].copy_from_slice(&f64::NAN.to_bits().to_le_bytes());
    assert!(MinimumMaterialCheckpoint3::from_bytes(&bad, Default::default()).is_err());
    for budget in [MinimumMaterialCodecBudget3 { max_bytes: bytes.len()-1, ..Default::default() },
        MinimumMaterialCodecBudget3 { max_items: 0, ..Default::default() }] {
        assert!(saved.to_bytes(budget).is_err());
        assert!(MinimumMaterialCheckpoint3::from_bytes(&bytes, budget).is_err());
    }
}
