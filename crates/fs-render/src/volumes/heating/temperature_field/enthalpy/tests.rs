use super::*;

fn spectrum(value: f64) -> MassAbsorptionSpectrum {
    MassAbsorptionSpectrum::source_table(&[(400.0, value), (700.0, value)], ContentHash([3; 32]))
        .unwrap()
}

#[test]
fn g0_enthalpy_interpolation_keeps_negative_datum_and_distinct_domain() {
    let table = EnthalpyAbsorptionTable::try_new(
        ContentHash([1; 32]),
        ContentHash([2; 32]),
        (500.0, 500.0),
        &[
            (-10.0, spectrum(0.2)),
            (0.0, spectrum(0.4)),
            (20.0, spectrum(0.8)),
        ],
        3,
    )
    .unwrap();
    let mut hash = DomainHasher::new("test.enthalpy-optical-table");
    let compiled = table
        .0
        .compile(&[550.0], &mut hash, &mut || Ok(()))
        .unwrap();
    for (h, expected) in [
        (-10.0, 0.2),
        (-5.0, 0.3),
        (0.0, 0.4),
        (10.0, 0.6),
        (20.0, 0.8),
    ] {
        let value = compiled.at(CellId(7), 500.0, h).unwrap()[0];
        assert!(
            (value - expected).abs() <= 2e-16,
            "h={h}: {value} != {expected}"
        );
    }
    assert!(matches!(
        compiled.at(CellId(7), 500.0, 20.1),
        Err(OpticalFeedbackError::EnthalpyDomain {
            thermal_cell: CellId(7),
            ..
        })
    ));
    assert_eq!(table.specific_enthalpy_bounds_j_kg(), (-10.0, 20.0));
    assert_eq!(table.temperature_bounds_k(), (500.0, 500.0));
}

#[test]
fn g0_enthalpy_table_refuses_missing_identity_unrepresentable_spans_and_bad_support() {
    let admit = |source, caloric, bounds, h: [f64; 2], max| {
        EnthalpyAbsorptionTable::try_new(
            source,
            caloric,
            bounds,
            &[(h[0], spectrum(0.2)), (h[1], spectrum(0.4))],
            max,
        )
    };
    let s = ContentHash([1; 32]);
    let c = ContentHash([2; 32]);
    for result in [
        admit(ContentHash([0; 32]), c, (500.0, 500.0), [0.0, 1.0], 2),
        admit(s, ContentHash([0; 32]), (500.0, 500.0), [0.0, 1.0], 2),
        admit(s, c, (0.0, 500.0), [0.0, 1.0], 2),
        admit(s, c, (501.0, 500.0), [0.0, 1.0], 2),
        admit(s, c, (500.0, f64::NAN), [0.0, 1.0], 2),
        admit(s, c, (500.0, 500.0), [0.0, 0.0], 2),
        admit(s, c, (500.0, 500.0), [1.0, 0.0], 2),
        admit(s, c, (500.0, 500.0), [0.0, f64::INFINITY], 2),
        admit(s, c, (500.0, 500.0), [-f64::MAX, f64::MAX], 2),
        admit(s, c, (500.0, 500.0), [0.0, 1.0], 1),
    ] {
        assert!(matches!(result, Err(OpticalFeedbackError::Invalid(_))));
    }
    let other = MassAbsorptionSpectrum::source_table(&[(450.0, 0.4), (650.0, 0.4)], s).unwrap();
    assert!(
        EnthalpyAbsorptionTable::try_new(
            s,
            c,
            (500.0, 500.0),
            &[(0.0, spectrum(0.2)), (1.0, other)],
            2
        )
        .is_err()
    );
}

#[test]
fn g4_compilation_cancels_before_sampling_next_enthalpy_row() {
    let table = EnthalpyScatteringTable::try_new(
        ContentHash([8; 32]),
        ContentHash([2; 32]),
        (500.0, 500.0),
        &[
            (0.0, MassScatteringSpectrum::constant(0.0).unwrap()),
            (10.0, MassScatteringSpectrum::constant(1.0).unwrap()),
        ],
        2,
    )
    .unwrap();
    let mut polls = 0;
    let mut hash = DomainHasher::new("test.enthalpy-optical-cancellation");
    let result = table.0.compile(&[550.0], &mut hash, &mut || {
        polls += 1;
        if polls == 2 {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    });
    assert!(matches!(
        result,
        Err(OpticalFeedbackError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(polls, 2);
}
