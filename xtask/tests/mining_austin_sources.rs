//! Actual Table 2 experimental fits -> pack compiler -> original batch PBM.
#![deny(unsafe_code)]
use fs_matdb::{NormalizedPack, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_material::comminution::{
    austin::resolve_austin_grinding,
    kinetics::{
        RATE_PREFIX, UniformizationBudget,
        target::{GrindTarget, GrindTargetError, GrindTargetSettings},
    },
};
use fs_material::separation::{ParticleClass, ParticleStream};
use fs_qty::{Dimensionless, Frequency, Length, Mass, MassFlowRate, Power, QuantitySpec, Time};
use std::{collections::BTreeMap, path::Path, process::Command};

fn conditions(filling: f64) -> QueryPoint {
    [
        ("ball-charge-volume-fraction", Dimensionless::DIMS, filling),
        ("powder-charge-volume-fraction", Dimensionless::DIMS, 0.04),
        ("mill-inside-diameter", Length::DIMS, 0.204),
        ("mill-inside-length", Length::DIMS, 0.166),
        ("grinding-ball-diameter", Length::DIMS, 0.0254),
        ("rotation-frequency", Frequency::DIMS, 1.1),
        ("dry-charge-mass", Mass::DIMS, 0.345),
        ("dry-grinding", Dimensionless::DIMS, 1.0),
        ("source-test-temperature-known", Dimensionless::DIMS, 0.0),
        ("source-test-pressure-known", Dimensionless::DIMS, 0.0),
        ("fitted-reference-model", Dimensionless::DIMS, 1.0),
        ("operating-ore-qualified", Dimensionless::DIMS, 0.0),
    ]
    .into_iter()
    .fold(QueryPoint::new(), |point, (axis, dims, value)| {
        point
            .with_quantity(axis, QuantitySpec::dimensional(dims), value)
            .unwrap()
    })
}
fn compile() -> NormalizedPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let output = std::env::temp_dir().join(format!(
        "quartz-austin-{}-{nonce}.fsmatpk",
        std::process::id()
    ));
    let result = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(root.join("data/matdb/seed-v1/quartz-assiros-petrakis-2017-austin/manifest.tsv"))
        .arg("--out")
        .arg(&output)
        .env("CARGO_WORKSPACE_DIR", root)
        .current_dir(root)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let bytes = std::fs::read(output).unwrap();
    let raw = NormalizedPack::from_bytes(&bytes).unwrap();
    NormalizedPack::from_bytes_verified(raw.content_hash(), &bytes).unwrap()
}
fn budget() -> UniformizationBudget {
    UniformizationBudget {
        maximum_l1_tail: Dimensionless::new(1e-13),
        max_poisson_terms: 256,
        max_coefficient_products: 1_000_000,
    }
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-12 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn experimental_fits_produce_source_bound_rates_and_conservative_psd_evolution() {
    let pack = compile();
    assert_eq!(pack.claims().claim_count(), 24);
    assert_eq!(pack.claims().observation_ids().count(), 3);
    // Independently transcribed Table 2, DOI10.14356/kona.2017007:
    // J, a_T/min, alpha, mu/mm, Lambda, Phi, gamma, beta, reported S(.85mm)/min.
    // These are the original fits, not the incomplete Table 3 back-fit.
    for (j, a, alpha, mu, lambda, phi, gamma, beta, reported_s) in [
        (0.2, 0.71, 1.13, 2.11, 2.90, 0.76, 1.02, 5.80, 0.55),
        (0.1, 0.35, 1.15, 2.34, 3.15, 0.77, 1.00, 5.80, 0.28),
        (0.067, 0.23, 1.17, 2.38, 3.33, 0.77, 1.11, 5.80, 0.18),
    ] {
        let query = conditions(j);
        let model = resolve_austin_grinding(
            &pack,
            &query,
            Time::new(240.0),
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
        let source = model.austin_calibration_source().unwrap();
        assert_eq!(source.pack_hash(), pack.content_hash());
        assert_eq!(source.conditions(), &query);
        assert_eq!(source.claims().len(), 8);
        assert_eq!(source.answers().len(), 13);
        assert!(source.observation().caveats.contains("PSD mismatch"));
        assert!(source.observation().caveats.contains("Beta fixed"));
        assert_eq!(
            source.claims()[0].observations,
            source.claims()[7].observations
        );
        for envelope in source.answers() {
            for answer in [&envelope.lower, &envelope.upper] {
                pack.claims().verify_receipt(&answer.receipt).unwrap();
                assert_eq!(
                    answer.evidence.value.uncertainty,
                    UncertaintyModel::Unstated
                );
            }
        }
        let rate = |id: &str| model.source().parameters[&format!("{RATE_PREFIX}{id}")].value;
        // Published rates are rounded to .01/min; compare at that printed precision.
        assert!((rate("size-000") * 60.0 - reported_s).abs() < 0.005);
        assert_eq!(rate("fines-sink"), 0.0);
        let s0 = a * 0.85_f64.powf(alpha) / (1.0 + (0.85_f64 / mu).powf(lambda));
        let s1 = a * 0.60_f64.powf(alpha) / (1.0 + (0.60_f64 / mu).powf(lambda));
        let ratio = 0.60_f64 / 0.85;
        let b10 = 1.0 - phi * ratio.powf(gamma) - (1.0 - phi) * ratio.powf(beta);
        close(rate("size-000") * 60.0, s0);
        close(rate("size-001") * 60.0, s1);
        close(
            model
                .daughters()
                .fraction("size-001", "size-000")
                .unwrap()
                .value(),
            b10,
        );
        let mut last_sink = 0.0;
        for seconds in [30.0, 60.0, 120.0, 240.0] {
            let pass = model
                .prepare(Time::new(seconds), &query, budget(), || false)
                .unwrap();
            let transition = pass.transition();
            let minutes = seconds / 60.0;
            // Independent analytical first two populations of the triangular
            // birth/death system; the production owner uses uniformization.
            close(
                transition.fraction("size-000", "size-000").unwrap().value(),
                (-s0 * minutes).exp(),
            );
            close(
                transition.fraction("size-001", "size-000").unwrap().value(),
                b10 * s0 * ((-s0 * minutes).exp() - (-s1 * minutes).exp()) / (s1 - s0),
            );
            for column in transition.columns() {
                close(
                    column
                        .transfers
                        .iter()
                        .map(|x| x.mass_fraction.value())
                        .sum(),
                    1.0,
                );
                assert!(
                    column
                        .transfers
                        .iter()
                        .all(|x| x.mass_fraction.value() >= 0.0)
                );
            }
            let sink = transition
                .fraction("fines-sink", "size-000")
                .unwrap()
                .value();
            assert!(sink > last_sink);
            last_sink = sink;
            assert_eq!(pass.model().austin_calibration_source(), Some(source));
        }
        let feed = ParticleStream::new(
            model
                .daughters()
                .classes()
                .iter()
                .map(|class| {
                    ParticleClass::new(
                        class.id.clone(),
                        class.diameter,
                        MassFlowRate::new(if class.id == "size-000" { 1.0 } else { 0.0 }),
                        BTreeMap::new(),
                    )
                    .unwrap()
                })
                .collect(),
            MassFlowRate::new(0.0),
        )
        .unwrap();
        let settings = GrindTargetSettings {
            maximum_contact_time: Time::new(240.0),
            time_tolerance: Time::new(1e-6),
            fraction_roundoff_allowance: Dimensionless::new(1e-12),
            maximum_forward_evaluations: 64,
            forward_budget: budget(),
        };
        let target = GrindTarget {
            coarsest_class: "size-001".into(),
            minimum_fraction: Dimensionless::new(0.5),
        };
        let solved = model
            .time_to_target(&feed, target.clone(), &query, settings, || false)
            .unwrap();
        let bracket = solved.bracket.unwrap();
        // Independent analytic crossing: top-class survival is exp(-S0*t).
        let exact_seconds = 60.0 * 2.0_f64.ln() / s0;
        assert!(bracket.below.contact_time.value() < exact_seconds);
        assert!(bracket.meeting.contact_time.value() > exact_seconds);
        assert!(bracket.meeting.contact_time.value() - bracket.below.contact_time.value() <= 1e-6);
        assert!(bracket.below.upper_fraction() < 0.5);
        assert!(bracket.meeting.lower_fraction() >= 0.5);
        assert_eq!(
            solved.pass.model().austin_calibration_source(),
            Some(source)
        );
        close(solved.product.dry_solids().unwrap().value(), 1.0);
        // Explicit test power is arithmetic input, not source-calibrated mill power.
        let energy = solved
            .energy_at_constant_net_power(Power::new(100.0), Mass::new(0.345))
            .unwrap();
        close(
            energy.meeting_energy.value(),
            100.0 * solved.pass.elapsed().value(),
        );
        assert!(
            solved
                .energy_at_constant_net_power(Power::new(100.0), Mass::new(1.0))
                .is_err()
        );
        assert!(matches!(
            model.time_to_target(
                &feed,
                GrindTarget {
                    minimum_fraction: Dimensionless::new(0.99),
                    ..target
                },
                &query,
                settings,
                || false
            ),
            Err(GrindTargetError::Unreachable { .. })
        ));
    }
    for query in [
        conditions(0.15),
        conditions(0.2).with("rotation-frequency", 1.2).unwrap(),
        conditions(0.2).with("dry-grinding", 0.0).unwrap(),
        conditions(0.2)
            .with("operating-ore-qualified", 1.0)
            .unwrap(),
    ] {
        assert!(
            resolve_austin_grinding(
                &pack,
                &query,
                Time::new(240.0),
                SelectionPolicy::SingleClaimOnly
            )
            .is_err()
        );
    }
    assert!(
        resolve_austin_grinding(
            &pack,
            &conditions(0.2),
            Time::new(241.0),
            SelectionPolicy::SingleClaimOnly
        )
        .is_err()
    );
}
