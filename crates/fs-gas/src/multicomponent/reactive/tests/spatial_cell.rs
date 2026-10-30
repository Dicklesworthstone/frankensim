use super::*;

fn source_and_cell<const D: usize>(
    mechanism: &IdealGasKineticsV1,
    velocity: [f64; D],
) -> (MixtureGas, Vec<f64>) {
    let gas = MixtureGas::new(
        mechanism.thermodynamics().clone(),
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
        || false,
    )
    .unwrap();
    let base = gas.conserved(&[1.0, 0.0], velocity[0], 1.0).unwrap();
    let mut row = vec![0.0; 2 + D + 1];
    row[..2].copy_from_slice(&base.partial_densities);
    row[2..2 + D].copy_from_slice(&velocity);
    row[2 + D] = base.energy;
    for &u in &velocity[1..] {
        row[2 + D] += 0.5 * u * u;
    }
    (gas, row)
}

fn run<const D: usize>(velocity: [f64; D]) -> (Vec<f64>, ChemistryReport) {
    let mechanism = mechanism(20.0, 50.0, 0.0);
    let (gas, mut q) = source_and_cell(&mechanism, velocity);
    let original = q.clone();
    let references = gas.energy_references().to_vec();
    let mut reactor = CellChemistry::new(gas, mechanism, chemistry(), || false).unwrap();
    let mut report = ChemistryReport::default();
    let volume = Volume::new(0.025);
    let duration = Time::new(0.02);
    match D {
        1 => reactor.advance_axial(&mut q, volume, duration, 7, &mut report, || false),
        2 => reactor.advance_planar(&mut q, volume, duration, 7, &mut report, || false),
        3 => reactor.advance_spatial(&mut q, volume, duration, 7, &mut report, || false),
        _ => unreachable!(),
    }
    .unwrap();
    assert_eq!(&q[2..2 + D], &original[2..2 + D]);
    let absolute_residual = q[2 + D] - original[2 + D]
        + (0..2)
            .map(|k| (q[k] - original[k]) * references[k])
            .sum::<f64>();
    assert!(absolute_residual.abs() < 2e-13);
    assert!((q[0] + q[1] - 1.0).abs() < 2e-14);
    assert!(q[1] > 0.5);
    assert!(report.entropy_change_j_k > 0.0);
    assert_eq!(report.cell_solves, 1);
    (q, report)
}

#[test]
fn g0_g3_original_cell_chemistry_retains_all_momenta_and_full_chemical_energy() {
    let (one, _) = run([0.7]);
    let (two, _) = run([0.7, -1.1]);
    let (three, _) = run([0.7, -1.1, 1.9]);
    for other in [&two, &three] {
        assert!((other[0] - one[0]).abs() < 2e-12);
        assert!((other[1] - one[1]).abs() < 2e-12);
    }
    let temperature = |q: &[f64], velocity: &[f64]| {
        let thermal = q[q.len() - 1] - velocity.iter().map(|u| 0.5 * u * u).sum::<f64>();
        // Independent constant-Cv inversion, R_mix=1, Cv=2.5, Tmin=.001.
        0.001 + thermal / 2.5
    };
    let t1 = temperature(&one, &[0.7]);
    close_temperature(t1, temperature(&two, &[0.7, -1.1]));
    close_temperature(t1, temperature(&three, &[0.7, -1.1, 1.9]));
    let (rotated, _) = run([1.9, 0.7, -1.1]);
    close_temperature(t1, temperature(&rotated, &[1.9, 0.7, -1.1]));
}

fn close_temperature(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-12, "{a} != {b}");
}

#[test]
fn g4_prepared_spatial_cell_publishes_cell_report_and_final_cancel_atomically() {
    let mechanism = mechanism(20.0, 50.0, 0.0);
    let (gas, initial) = source_and_cell(&mechanism, [0.7, -1.1, 1.9]);
    let mut source = CellChemistry::new(gas, mechanism, chemistry(), || false).unwrap();
    let mut expected = initial.clone();
    let mut report = ChemistryReport::default();
    let mut polls = 0;
    let volume = Volume::new(0.025);
    let duration = Time::new(0.02);
    source
        .advance_spatial(&mut expected, volume, duration, 9, &mut report, || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls > 10);
    let mut actual = initial.clone();
    let mut actual_report = ChemistryReport::default();
    let mut count = 0;
    let error = source
        .advance_spatial(&mut actual, volume, duration, 9, &mut actual_report, || {
            count += 1;
            count >= polls
        })
        .unwrap_err();
    assert_eq!(error, ReactiveError::Cancelled);
    assert_eq!(actual, initial);
    assert_eq!(actual_report, ChemistryReport::default());
    let bytes = source.workspace_bytes();
    source
        .advance_spatial(&mut actual, volume, duration, 9, &mut actual_report, || {
            false
        })
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(actual_report, report);
    assert_eq!(source.workspace_bytes(), bytes);
    assert!(
        source
            .advance_spatial(
                &mut actual[..5],
                volume,
                duration,
                9,
                &mut actual_report,
                || false
            )
            .is_err()
    );
    assert_eq!(actual, expected);
    assert_eq!(actual_report, report);
}
