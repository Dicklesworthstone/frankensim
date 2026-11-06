//! Synthetic format and equation checks, not evaluated-isotope benchmarks.

use super::*;
use crate::diffusion::{
    DelayedGroup, DiffusionBudget, DiffusionData, DiffusionLimits, DiffusionModel,
    DiffusionReaction, DiffusionReactionKind, DiffusionState,
};
use crate::endf::spectrum::{FluxSpectrumData, FluxSpectrumLimits};
use crate::endf::{EV_J, EndfInterpolation, EndfInterpolationRegion};

fn real(value: f64) -> String {
    let field = format!("{value:11.4e}");
    assert_eq!(field.len(), 11);
    field
}

fn record(fields: [String; 6], mt: i32) -> String {
    format!("{}{:4}{:2}{mt:3}{:5}", fields.concat(), 101, 1, 1)
}

fn control(c1: f64, c2: f64, ints: [i32; 4], mt: i32) -> String {
    record(
        [
            real(c1),
            real(c2),
            format!("{:11}", ints[0]),
            format!("{:11}", ints[1]),
            format!("{:11}", ints[2]),
            format!("{:11}", ints[3]),
        ],
        mt,
    )
}

fn selection(kind: FissionNeutronKind) -> EndfNeutronYieldSelection {
    EndfNeutronYieldSelection {
        material_number: 101,
        nuclide: NuclideId::new(10, 20, 1).unwrap(),
        temperature_k: 0.0,
        kind,
    }
}

fn fixture(kind: FissionNeutronKind, points: &[(f64, f64)], lambdas: &[f64]) -> String {
    let mt = kind.mt();
    let mut lines = vec![
        control(10020.0, 19.8, [1, 1, 0, 3], 451),
        control(1.0, 0.0, [1, 1, 0, 6], 451),
        control(1.0, 2e7, [2, 0, 10, 8], 451),
        control(0.0, 0.0, [0, 0, 0, 0], 451),
        control(0.0, 0.0, [0; 4], 0),
        control(10020.0, 19.8, [0, 2, 0, 0], mt),
    ];
    if kind == FissionNeutronKind::Delayed {
        lines.push(control(0.0, 0.0, [0, 0, lambdas.len() as i32, 0], mt));
        for chunk in lambdas.chunks(6) {
            let mut fields = std::array::from_fn(|_| real(0.0));
            for (i, &lambda) in chunk.iter().enumerate() {
                fields[i] = real(lambda);
            }
            lines.push(record(fields, mt));
        }
    }
    lines.push(control(0.0, 0.0, [0, 0, 1, points.len() as i32], mt));
    lines.push(record(
        [
            format!("{:11}", points.len()),
            format!("{:11}", 2),
            format!("{:11}", 0),
            format!("{:11}", 0),
            format!("{:11}", 0),
            format!("{:11}", 0),
        ],
        mt,
    ));
    for chunk in points.chunks(3) {
        let mut fields = std::array::from_fn(|_| real(0.0));
        for (i, &(energy, value)) in chunk.iter().enumerate() {
            fields[2 * i] = real(energy);
            fields[2 * i + 1] = real(value);
        }
        lines.push(record(fields, mt));
    }
    lines.push(control(0.0, 0.0, [0; 4], 0));
    lines.join("\n") + "\n"
}

fn parse(text: &str, kind: FissionNeutronKind) -> Result<FissionNeutronYield, EndfError> {
    FissionNeutronYield::parse(
        "manufactured-neutron-yield".into(),
        selection(kind),
        text,
        EndfNeutronYieldReadLimits::default(),
        || Ok(()),
    )
}

fn region(points: usize) -> Arc<[EndfInterpolationRegion]> {
    vec![EndfInterpolationRegion {
        end_point: points,
        interpolation: EndfInterpolation::Linear,
    }]
    .into()
}

fn cross(yield_: &FissionNeutronYield, values: &[f64]) -> PointwiseCrossSection {
    let mut metadata = yield_.metadata().clone();
    metadata.temperature_k = 300.0;
    metadata.resonance_flag = 2;
    metadata.derivation_flag = 1;
    PointwiseCrossSection {
        data_id: "manufactured-processed-fission".into(),
        metadata,
        reaction: NeutronReaction::Fission,
        energies_j: vec![0.0, 2.0 * EV_J].into(),
        cross_sections_m2: values.to_vec().into(),
        interpolation_regions: region(2),
        mass_difference_q_j: 0.0,
        reaction_q_j: 0.0,
    }
}

fn grid() -> Arc<NeutronGroups> {
    Arc::new(NeutronGroups::new("synthetic-group".into(), vec![0.0, 2.0 * EV_J]).unwrap())
}

fn spectrum() -> FluxSpectrum {
    FluxSpectrum::new(
        FluxSpectrumData {
            data_id: "synthetic-absolute-flux".into(),
            energies_j: vec![0.0, 2.0 * EV_J],
            flux_density_m2_s_j: vec![0.0, 2.0 / EV_J],
            interpolation_regions: region(2).to_vec(),
        },
        FluxSpectrumLimits::default(),
    )
    .unwrap()
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-10 * expected.abs().max(1e-290),
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn imports_distinct_multiplicities_and_ordered_delayed_lifetimes_from_raw_evaluation() {
    for kind in [
        FissionNeutronKind::Total,
        FissionNeutronKind::Prompt,
        FissionNeutronKind::Delayed,
    ] {
        let source = fixture(kind, &[(0.0, 2.0), (2.0, 4.0)], &[0.5, 0.1]);
        let law = parse(&source, kind).unwrap();
        assert_eq!(law.kind(), kind);
        assert_eq!(law.data_id(), "manufactured-neutron-yield");
        assert_eq!(law.metadata().nuclide, selection(kind).nuclide);
        assert_eq!(law.metadata().resonance_flag, 1);
        assert_eq!(law.metadata().temperature_k, 0.0);
        assert_eq!(law.curve().values(), &[2.0, 4.0]);
        near(law.curve().value_at_energy(EV_J).unwrap(), 3.0);
        assert_eq!(
            law.decay_constants_s_inv(),
            if kind == FissionNeutronKind::Delayed {
                &[0.5, 0.1][..]
            } else {
                &[]
            }
        );
    }
}

#[test]
fn imported_prompt_law_preserves_production_and_drives_original_diffusion() {
    let law = parse(
        &fixture(FissionNeutronKind::Prompt, &[(0.0, 2.0), (2.0, 4.0)], &[]),
        FissionNeutronKind::Prompt,
    )
    .unwrap();
    let xs = cross(&law, &[0.0, 2e-28]);
    let condensed = law
        .collapse(
            &xs,
            &spectrum(),
            grid(),
            ReactionWeightedLimits::default(),
            || Ok(()),
        )
        .unwrap();
    // sigma=E barns, phi=E per eV, nu=2+E. Integral(sigma*phi*nu)=28/3 barns/s.
    let row = &condensed.collapse;
    let nu = row.yields[0].values[0].unwrap();
    near(nu, 3.5);
    near(
        row.cross_sections_m2[0] * row.flux.values_m2_s()[0] * nu,
        28.0 / 3.0 * 1e-28,
    );
    assert_eq!(condensed.kind, FissionNeutronKind::Prompt);
    assert_eq!(condensed.evaluation_metadata.temperature_k, 0.0);
    assert_eq!(row.metadata.temperature_k, 300.0);
    let limits = DiffusionLimits::default();
    let model = Arc::new(
        DiffusionModel::new(
            DiffusionData {
                data_id: "synthetic-imported-prompt".into(),
                groups: grid(),
                speeds_m_s: vec![1.0],
                cell_volumes_m3: vec![1.0],
                delayed_groups: vec![],
                links: vec![],
                leakage: vec![],
                reactions: vec![DiffusionReaction {
                    data_id: law.data_id().into(),
                    cell: 0,
                    incident_group: 0,
                    macroscopic_cross_section_m_inv: 7.5e27 * row.cross_sections_m2[0],
                    kind: DiffusionReactionKind::Fission {
                        prompt_yields: vec![nu],
                        delayed_yields: vec![],
                    },
                }],
            },
            limits,
        )
        .unwrap(),
    );
    let state = DiffusionState::new(model, vec![10.0], vec![]).unwrap();
    let (next, report) = state
        .trial(
            0.2,
            DiffusionBudget {
                limits,
                population_tolerance: 1e-8,
            },
            || Ok(()),
        )
        .unwrap();
    near(next.neutron_counts()[0], 10.0 * 0.5_f64.exp());
    near(report.reaction_counts[0], 10.0 * 0.5_f64.exp_m1() / 2.5);
}

#[test]
fn imported_delayed_lifetime_remains_a_neutron_history_not_instant_prompt_birth() {
    let law = parse(
        &fixture(
            FissionNeutronKind::Delayed,
            &[(0.0, 1.0), (2.0, 1.0)],
            &[0.5],
        ),
        FissionNeutronKind::Delayed,
    )
    .unwrap();
    let condensed = law
        .collapse(
            &cross(&law, &[1e-28, 1e-28]),
            &spectrum(),
            grid(),
            ReactionWeightedLimits::default(),
            || Ok(()),
        )
        .unwrap();
    // The manufactured single family explicitly has abundance 1 and a unit
    // outgoing-group spectrum; neither is inferred by the importer.
    let limits = DiffusionLimits::default();
    let model = Arc::new(
        DiffusionModel::new(
            DiffusionData {
                data_id: "synthetic-imported-delayed".into(),
                groups: grid(),
                speeds_m_s: vec![1.0],
                cell_volumes_m3: vec![1.0],
                links: vec![],
                leakage: vec![],
                delayed_groups: vec![DelayedGroup {
                    data_id: law.data_id().into(),
                    decay_constant_s_inv: condensed.decay_constants_s_inv[0],
                    spectrum: vec![1.0],
                }],
                reactions: vec![DiffusionReaction {
                    data_id: law.data_id().into(),
                    cell: 0,
                    incident_group: 0,
                    macroscopic_cross_section_m_inv: 1e28 * condensed.collapse.cross_sections_m2[0],
                    kind: DiffusionReactionKind::Fission {
                        prompt_yields: vec![0.0],
                        delayed_yields: vec![condensed.collapse.yields[0].values[0].unwrap()],
                    },
                }],
            },
            limits,
        )
        .unwrap(),
    );
    let state = DiffusionState::new(model, vec![10.0], vec![0.0]).unwrap();
    let (next, _) = state
        .trial(
            0.7,
            DiffusionBudget {
                limits,
                population_tolerance: 1e-8,
            },
            || Ok(()),
        )
        .unwrap();
    near(
        next.neutron_counts()[0],
        10.0 / 3.0 + 20.0 / 3.0 * (-1.05_f64).exp(),
    );
    near(
        next.precursor_counts()[0],
        20.0 / 3.0 * -(-1.05_f64).exp_m1(),
    );
}

#[test]
fn collapse_refuses_mixed_evaluations_and_preserves_undefined_zero_reaction_yields() {
    let law = parse(
        &fixture(FissionNeutronKind::Total, &[(0.0, 2.0), (2.0, 3.0)], &[]),
        FissionNeutronKind::Total,
    )
    .unwrap();
    let xs = cross(&law, &[1e-28, 1e-28]);
    for field in 0..7 {
        let mut wrong = xs.clone();
        match field {
            0 => wrong.reaction = NeutronReaction::Capture,
            1 => wrong.metadata.nuclide = NuclideId::new(10, 20, 0).unwrap(),
            2 => wrong.metadata.atomic_weight_ratio += 1.0,
            3 => wrong.metadata.library_id += 1,
            4 => wrong.metadata.library_version += 1,
            5 => wrong.metadata.library_release += 1,
            _ => wrong.metadata.material_revision += 1,
        }
        assert!(
            law.collapse(
                &wrong,
                &spectrum(),
                grid(),
                ReactionWeightedLimits::default(),
                || Ok(())
            )
            .is_err()
        );
    }
    let zero = law
        .collapse(
            &cross(&law, &[0.0, 0.0]),
            &spectrum(),
            grid(),
            ReactionWeightedLimits::default(),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(zero.collapse.yields[0].values, vec![None]);
}

#[test]
fn strict_selection_format_and_unsupported_representations_refuse() {
    let kind = FissionNeutronKind::Prompt;
    let source = fixture(kind, &[(0.0, 2.0), (2.0, 3.0)], &[]);
    // NSUB=10 incident-neutron data cannot declare LRP=-1, even when the
    // multiplicity reader does not require reconstructed MF3 cross sections.
    assert!(
        parse(
            &source.replace(
                &control(10020.0, 19.8, [1, 1, 0, 3], 451),
                &control(10020.0, 19.8, [-1, 1, 0, 3], 451),
            ),
            kind,
        )
        .is_err()
    );
    let head = control(10020.0, 19.8, [0, 2, 0, 0], kind.mt());
    for changed in [
        control(10021.0, 19.8, [0, 2, 0, 0], kind.mt()),
        control(10020.0, 19.9, [0, 2, 0, 0], kind.mt()),
        control(10020.0, 19.8, [0, 1, 0, 0], kind.mt()),
        control(10020.0, 19.8, [0, 2, 1, 0], kind.mt()),
    ] {
        assert!(parse(&source.replace(&head, &changed), kind).is_err());
    }
    assert!(matches!(
        parse(&(source.clone() + &source), kind),
        Err(EndfError::AmbiguousSelection)
    ));
    assert!(matches!(
        parse(&source, FissionNeutronKind::Total),
        Err(EndfError::MissingSelection)
    ));
    let mut wrong = selection(kind);
    wrong.nuclide = NuclideId::new(10, 20, 0).unwrap();
    assert!(
        FissionNeutronYield::parse(
            "test".into(),
            wrong,
            &source,
            EndfNeutronYieldReadLimits::default(),
            || Ok(())
        )
        .is_err()
    );
    let rows: Vec<_> = source.lines().collect();
    for length in 6..rows.len() {
        assert!(parse(&(rows[..length].join("\n") + "\n"), kind).is_err());
    }
    assert!(parse(&fixture(kind, &[(0.0, -1.0), (2.0, 3.0)], &[]), kind).is_err());
    let delayed = fixture(
        FissionNeutronKind::Delayed,
        &[(0.0, 1.0), (2.0, 1.0)],
        &[0.0],
    );
    assert!(parse(&delayed, FissionNeutronKind::Delayed).is_err());
    let delayed = fixture(
        FissionNeutronKind::Delayed,
        &[(0.0, 1.0), (2.0, 1.0)],
        &[0.5],
    );
    assert!(matches!(
        parse(
            &delayed.replace(
                &control(10020.0, 19.8, [0, 2, 0, 0], 455),
                &control(10020.0, 19.8, [1, 2, 0, 0], 455)
            ),
            FissionNeutronKind::Delayed
        ),
        Err(EndfError::Unsupported(_))
    ));
}

#[test]
fn bounded_import_and_late_cancellation_never_publish_partial_data() {
    let kind = FissionNeutronKind::Delayed;
    let source = fixture(kind, &[(0.0, 0.1), (2.0, 0.2)], &[0.1, 0.5]);
    for field in 0..5 {
        let mut limits = EndfNeutronYieldReadLimits::default();
        match field {
            0 => limits.records.max_bytes = source.len() - 1,
            1 => limits.records.max_records = source.lines().count() - 1,
            2 => limits.records.max_points = 1,
            3 => limits.records.max_interpolation_regions = 0,
            _ => limits.max_delayed_families = 1,
        }
        assert!(
            FissionNeutronYield::parse("test".into(), selection(kind), &source, limits, || Ok(()))
                .is_err()
        );
    }
    let mut polls = 0;
    FissionNeutronYield::parse(
        "test".into(),
        selection(kind),
        &source,
        EndfNeutronYieldReadLimits::default(),
        || {
            polls += 1;
            Ok(())
        },
    )
    .unwrap();
    let mut count = 0;
    let refused = FissionNeutronYield::parse(
        "test".into(),
        selection(kind),
        &source,
        EndfNeutronYieldReadLimits::default(),
        || {
            count += 1;
            if count == polls {
                Err(FissionError::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert!(matches!(
        refused,
        Err(EndfError::Physics(FissionError::Cancelled))
    ));
}
