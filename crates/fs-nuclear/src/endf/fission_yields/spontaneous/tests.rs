//! Manufactured format records exercise the source boundary, not isotope data.

use super::*;

fn selection() -> EndfFissionYieldSelection {
    EndfFissionYieldSelection {
        material_number: 777,
        nuclide: NuclideId::new(10, 20, 1).unwrap(),
    }
}

fn real(value: f64) -> String {
    for precision in (0..=6).rev() {
        let number = format!("{value:.precision$e}");
        let (mantissa, exponent) = number.split_once('e').unwrap();
        let exponent: i32 = exponent.parse().unwrap();
        let field = format!("{mantissa}{exponent:+}");
        if field.len() <= 11 {
            return format!("{field:>11}");
        }
    }
    panic!("manufactured value does not fit the source field");
}

fn record(fields: [String; 6], mf: i32, mt: i32) -> String {
    assert!(fields.iter().all(|field| field.len() == 11));
    format!("{}{:4}{mf:2}{mt:3}{:5}", fields.concat(), 777, 1)
}

fn control(c1: f64, c2: f64, integers: [i32; 4], mf: i32, mt: i32) -> String {
    record(
        [
            real(c1),
            real(c2),
            format!("{:11}", integers[0]),
            format!("{:11}", integers[1]),
            format!("{:11}", integers[2]),
            format!("{:11}", integers[3]),
        ],
        mf,
        mt,
    )
}

fn fixture(sublibrary: i32, energies: usize, mt: i32) -> String {
    let mut lines = vec![
        control(10020.0, 20.0, [-1, 0, 0, 3], 1, 451),
        control(0.0, 0.0, [1, 1, 0, 6], 1, 451),
        control(0.0, 0.0, [0, 0, sublibrary, 8], 1, 451),
        control(0.0, 0.0, [0, 0, 0, 0], 1, 451),
        control(0.0, 0.0, [0, 0, 0, 0], 1, 0),
        control(10020.0, 20.0, [energies as i32, 0, 0, 0], 8, mt),
    ];
    for row in 0..energies {
        let law = if row == 0 { energies as i32 - 1 } else { 2 };
        lines.push(control(row as f64, 0.0, [law, 0, 8, 2], 8, mt));
        lines.push(record(
            [4008.0, 0.0, 0.5, 0.01, 4008.0, 1.0].map(real),
            8,
            mt,
        ));
        lines.push(record([1.75, 0.02, 0.0, 0.0, 0.0, 0.0].map(real), 8, mt));
    }
    lines.push(control(0.0, 0.0, [0, 0, 0, 0], 8, 0));
    lines.join("\n")
}

fn parse(text: &str) -> Result<SpontaneousFissionYieldEvaluation, EndfError> {
    SpontaneousFissionYieldEvaluation::parse(
        "manufactured spontaneous products".into(),
        selection(),
        text,
        EndfFissionYieldReadLimits::default(),
        || Ok(()),
    )
}

#[test]
fn independent_spontaneous_products_retain_identity_yields_and_uncertainties() {
    let yields = parse(&fixture(5, 1, 454)).unwrap();
    assert_eq!(yields.metadata().nuclide, selection().nuclide);
    assert_eq!(yields.metadata().material_number, 777);
    assert_eq!(yields.metadata().library_version, 8);
    assert_eq!(yields.data_id(), "manufactured spontaneous products");
    assert_eq!(
        yields.products(),
        &[
            NuclideId::new(4, 8, 0).unwrap(),
            NuclideId::new(4, 8, 1).unwrap(),
        ]
    );
    assert_eq!(yields.values()[0].yield_per_fission, 0.5);
    assert_eq!(yields.values()[1].yield_per_fission, 1.75);
    assert_eq!(yields.values()[0].uncertainty, 0.01);
    assert_eq!(yields.values()[1].uncertainty, 0.02);
    assert_eq!(yields.source_energy_j(), 0.0);
    // Preserve supplied multiplicities; do not normalize them to one or two.
    assert_eq!(
        yields
            .values()
            .iter()
            .map(|x| x.yield_per_fission)
            .sum::<f64>(),
        2.25
    );
    assert_eq!(yields, yields.clone());
}

#[test]
fn spontaneous_and_induced_sources_cannot_substitute_for_each_other() {
    for sublibrary in [1, 4, 6, 11] {
        assert!(matches!(
            parse(&fixture(sublibrary, 1, 454)),
            Err(EndfError::Unsupported(_))
        ));
    }
    assert!(matches!(
        parse(&fixture(5, 2, 454)),
        Err(EndfError::Unsupported(_))
    ));
    assert_eq!(parse(&fixture(5, 1, 459)), Err(EndfError::MissingSelection));
    let induced = |text: &str| {
        FissionYieldEvaluation::parse(
            "manufactured induced products".into(),
            selection(),
            text,
            EndfFissionYieldReadLimits::default(),
            || Ok(()),
        )
    };
    assert!(induced(&fixture(11, 1, 454)).is_ok());
    assert!(matches!(
        induced(&fixture(5, 1, 454)),
        Err(EndfError::Unsupported(_))
    ));
    let text = fixture(5, 1, 454);
    let duplicate = format!("{text}\n{text}");
    assert_eq!(parse(&duplicate), Err(EndfError::AmbiguousSelection));
    assert!(
        SpontaneousFissionYieldEvaluation::parse(
            "changed parent".into(),
            EndfFissionYieldSelection {
                nuclide: NuclideId::new(10, 20, 0).unwrap(),
                ..selection()
            },
            &text,
            EndfFissionYieldReadLimits::default(),
            || Ok(()),
        )
        .is_err()
    );
}

#[test]
fn shared_source_limits_and_final_cancellation_refuse_before_publication() {
    let text = fixture(5, 1, 454);
    let limits = EndfFissionYieldReadLimits::default();
    for tightened in [
        EndfFissionYieldReadLimits {
            max_bytes: text.len() - 1,
            ..limits
        },
        EndfFissionYieldReadLimits {
            max_records: 9,
            ..limits
        },
        EndfFissionYieldReadLimits {
            max_energy_points: 0,
            ..limits
        },
        EndfFissionYieldReadLimits {
            max_products_per_energy: 1,
            ..limits
        },
        EndfFissionYieldReadLimits {
            max_product_values: 1,
            ..limits
        },
    ] {
        assert!(matches!(
            SpontaneousFissionYieldEvaluation::parse(
                "bounded source".into(),
                selection(),
                &text,
                tightened,
                || Ok(()),
            ),
            Err(EndfError::Physics(FissionError::Budget(_)))
        ));
    }
    let mut polls = 0;
    SpontaneousFissionYieldEvaluation::parse(
        "poll count".into(),
        selection(),
        &text,
        limits,
        || {
            polls += 1;
            Ok(())
        },
    )
    .unwrap();
    let mut seen = 0;
    assert_eq!(
        SpontaneousFissionYieldEvaluation::parse(
            "cancelled source".into(),
            selection(),
            &text,
            limits,
            || {
                seen += 1;
                if seen == polls {
                    Err(FissionError::Cancelled)
                } else {
                    Ok(())
                }
            },
        ),
        Err(EndfError::Physics(FissionError::Cancelled))
    );
    assert!(parse(&text).is_ok());
}
