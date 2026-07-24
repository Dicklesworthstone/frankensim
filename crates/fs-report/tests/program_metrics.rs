//! G0/G5 battery for the deterministic program metrics dashboard.
//!
//! The load-bearing distinction is between a measured zero and an unmeasured
//! `NO-DATA`: one is a known-bad outcome that must stay visible, the other is
//! an absent measurement that must never be rendered as a number. The trend
//! tests pin the arithmetic across two recorded generations in both
//! directions, including the seeded regression that must visibly flip a cell.
//! Determinism (G5) is byte-identity of both renders across independent
//! builds from independently constructed inputs.

use std::collections::BTreeMap;
use std::num::NonZeroU64;

use fs_blake3::ContentHash;
use fs_report::program_metrics::{
    DashboardError, HistoryGeneration, MAX_METRIC_ROWS, MetricCell, MetricDirection, MetricFamily,
    MetricHistory, MetricObservation, MetricRow, Trend, build_dashboard,
};

fn nonzero(value: u64) -> NonZeroU64 {
    NonZeroU64::new(value).expect("fixture denominators are nonzero")
}

fn ratio(numerator: u64, denominator: u64) -> MetricCell {
    MetricCell::Measured(MetricObservation::ratio(numerator, nonzero(denominator)))
}

fn count(value: u64) -> MetricCell {
    MetricCell::Measured(MetricObservation::count(value))
}

fn no_data(reason: &str, unblocked_by: Option<&str>) -> MetricCell {
    MetricCell::NoData {
        reason: reason.to_string(),
        unblocked_by: unblocked_by.map(str::to_string),
    }
}

fn row(id: &str, direction: MetricDirection, cell: MetricCell) -> MetricRow {
    let sources: &[&str] = if cell.is_measured() {
        &["vv-scorecard.json"]
    } else {
        &[]
    };
    MetricRow::try_new(
        id,
        &format!("title for {id}"),
        MetricFamily::Outcome,
        direction,
        cell,
        sources,
        "fixture row; captures nothing about the real program",
    )
    .expect("fixture rows admit cleanly")
}

fn history_of(pairs: &[(&str, MetricCell)], source_identity: ContentHash) -> MetricHistory {
    let metrics: BTreeMap<String, MetricCell> = pairs
        .iter()
        .map(|(id, cell)| ((*id).to_string(), cell.clone()))
        .collect();
    let generation = HistoryGeneration::try_new(1, "fixture generation", source_identity, metrics)
        .expect("fixture generation admits");
    MetricHistory::try_new(vec![generation]).expect("single-generation history admits")
}

/// A measured zero is a KNOWN outcome and must stay visible as `0`; only an
/// absent measurement is `NO-DATA`. Conflating them in either direction is the
/// central failure this artifact exists to prevent.
#[test]
fn measured_zero_is_not_no_data() {
    let rows = vec![
        row(
            "blind-datasets",
            MetricDirection::HigherIsBetter,
            ratio(0, 28),
        ),
        row(
            "decision-turnaround",
            MetricDirection::LowerIsBetter,
            no_data("no acceptance lane records timings", Some("f85xj.6.11")),
        ),
    ];
    let dashboard = build_dashboard(&rows, &MetricHistory::empty()).expect("dashboard builds");
    assert_eq!(dashboard.measured_count(), 1);
    assert_eq!(dashboard.no_data_count(), 1);

    let markdown = dashboard.render_markdown();
    let json = dashboard.render_json();

    // The measured zero is reported as a real zero, with its denominator.
    assert!(markdown.contains("0 of 28 (0.00%)"), "{markdown}");
    assert!(
        json.contains("\"status\":\"ratio\",\"numerator\":0,\"denominator\":28"),
        "{json}"
    );
    // The unmeasured row is loud and names the work that would make it live.
    assert!(
        markdown.contains("NO-DATA (needs f85xj.6.11)"),
        "{markdown}"
    );
    assert!(json.contains("\"status\":\"no-data\""), "{json}");

    // Negative assertions carry the doctrine: a measured zero must never be
    // laundered into NO-DATA, and an unmeasured row must never acquire a number.
    assert!(
        !markdown.contains("| NO-DATA | no prior generation | higher-is-better |"),
        "a measured zero was rendered as NO-DATA: {markdown}"
    );
    assert!(
        !json.contains("\"id\":\"decision-turnaround\",\"title\":\"title for decision-turnaround\",\"family\":\"outcome\",\"direction\":\"lower-is-better\",\"value\":{\"status\":\"count\""),
        "an unmeasured row acquired a count: {json}"
    );
}

/// An empty population cannot masquerade as `0%`: the denominator is
/// `NonZeroU64`, so the zero-denominator ratio is unrepresentable rather than
/// merely unchecked. Percentages are exact integer basis points.
#[test]
fn ratios_are_exact_and_denominators_are_nonzero() {
    assert_eq!(NonZeroU64::new(0), None);

    let third = MetricObservation::ratio(1, nonzero(3));
    assert_eq!(third.basis_points(), Some(3_333));
    assert_eq!(third.render(), "1 of 3 (33.33%)");
    assert_eq!(MetricObservation::count(7).basis_points(), None);

    // Cross-multiplied comparison: 13/25 > 1/3 exactly, with no float.
    assert_eq!(
        MetricObservation::ratio(13, nonzero(25)).compare(third),
        Some(core::cmp::Ordering::Greater)
    );
    // Different shapes claim no movement rather than inventing one.
    assert_eq!(third.compare(MetricObservation::count(1)), None);
}

/// Trend arithmetic across two generations, in both declared directions.
#[test]
fn trend_arithmetic_across_two_generations() {
    let prior = [
        ("higher", ratio(10, 25)),
        ("lower", ratio(15, 25)),
        ("neutral", count(4)),
        ("appearing", no_data("not wired", None)),
        ("vanishing", count(3)),
        ("shape-shift", ratio(1, 2)),
    ];
    let history = history_of(&prior, ContentHash([9; 32]));

    let rows = vec![
        // higher-is-better moving up is an improvement.
        row("higher", MetricDirection::HigherIsBetter, ratio(13, 25)),
        // lower-is-better moving up is a regression.
        row("lower", MetricDirection::LowerIsBetter, ratio(20, 25)),
        // a neutral metric reports movement without a judgement.
        row("neutral", MetricDirection::Neutral, count(9)),
        row("appearing", MetricDirection::HigherIsBetter, count(2)),
        row(
            "vanishing",
            MetricDirection::HigherIsBetter,
            no_data("source retired", None),
        ),
        row("shape-shift", MetricDirection::HigherIsBetter, count(1)),
    ];
    let dashboard = build_dashboard(&rows, &history).expect("dashboard builds");
    let trends: BTreeMap<&str, Trend> = dashboard
        .rows()
        .iter()
        .map(|entry| (entry.metric().id(), entry.trend()))
        .collect();

    assert_eq!(trends["higher"], Trend::Improved);
    assert_eq!(trends["lower"], Trend::Worsened);
    assert_eq!(trends["neutral"], Trend::Changed);
    assert_eq!(trends["appearing"], Trend::NewlyMeasured);
    assert_eq!(trends["vanishing"], Trend::LostData);
    assert_eq!(trends["shape-shift"], Trend::NotComparable);

    // An unrecorded metric has no prior, which is distinct from "unchanged".
    let fresh = vec![row(
        "never-recorded",
        MetricDirection::HigherIsBetter,
        count(1),
    )];
    let baseline = build_dashboard(&fresh, &history).expect("dashboard builds");
    assert_eq!(baseline.rows()[0].trend(), Trend::NoPrior);
}

/// The bead's required regression drill: a worsened metric must visibly flip
/// the trend cell, in both the machine-readable and human-readable renders.
#[test]
fn seeded_regression_flips_the_trend_cell() {
    let healthy = vec![row(
        "external-anchoring",
        MetricDirection::HigherIsBetter,
        ratio(13, 25),
    )];
    let first = build_dashboard(&healthy, &MetricHistory::empty()).expect("dashboard builds");
    assert_eq!(first.rows()[0].trend(), Trend::NoPrior);
    assert!(!first.render_markdown().contains("WORSENED"));

    // Record generation 1, then seed a regression against it.
    let recorded = history_of(
        &[("external-anchoring", ratio(13, 25))],
        first.source_identity(),
    );
    let regressed = vec![row(
        "external-anchoring",
        MetricDirection::HigherIsBetter,
        ratio(9, 25),
    )];
    let second = build_dashboard(&regressed, &recorded).expect("dashboard builds");

    assert_eq!(second.rows()[0].trend(), Trend::Worsened);
    assert!(
        second.render_markdown().contains("WORSENED"),
        "regression must be loud"
    );
    assert!(second.render_json().contains("\"trend\":\"worsened\""));
    assert_eq!(second.prior_generation(), Some(1));
    // The regression moved the values, so it cannot claim identical sources.
    assert!(!second.sources_identical_to_prior());
}

/// "Unchanged" is provable, not coincidental: identical source identity is
/// what distinguishes a genuinely static program from two renders that merely
/// happen to agree.
#[test]
fn unchanged_sources_are_provably_unchanged() {
    let rows = vec![row(
        "external-anchoring",
        MetricDirection::HigherIsBetter,
        ratio(13, 25),
    )];
    let first = build_dashboard(&rows, &MetricHistory::empty()).expect("dashboard builds");

    let recorded = history_of(
        &[("external-anchoring", ratio(13, 25))],
        first.source_identity(),
    );
    let second = build_dashboard(&rows, &recorded).expect("dashboard builds");

    assert_eq!(second.rows()[0].trend(), Trend::Unchanged);
    assert!(second.sources_identical_to_prior());
    assert!(
        second
            .render_markdown()
            .contains("provably unchanged, not coincidentally equal")
    );

    // A history whose recorded identity disagrees cannot claim identical sources.
    let mismatched = history_of(
        &[("external-anchoring", ratio(13, 25))],
        ContentHash([0xab; 32]),
    );
    let third = build_dashboard(&rows, &mismatched).expect("dashboard builds");
    assert_eq!(third.rows()[0].trend(), Trend::Unchanged);
    assert!(!third.sources_identical_to_prior());
}

/// The recorded history line carries values and source identity only — never
/// trends, which are derived — so recording cannot smuggle a computed verdict
/// into the next generation's basis.
#[test]
fn history_line_records_values_not_trends() {
    let rows = vec![
        row("measured", MetricDirection::HigherIsBetter, ratio(13, 25)),
        row(
            "absent",
            MetricDirection::LowerIsBetter,
            no_data("not wired", Some("f85xj.8.7")),
        ),
    ];
    let dashboard = build_dashboard(&rows, &MetricHistory::empty()).expect("dashboard builds");
    let line = dashboard.render_history_line(1, "2026-07-24 baseline");

    assert!(line.ends_with('\n'), "history is JSONL");
    assert_eq!(line.lines().count(), 1);
    assert!(line.contains("\"generation\":1"));
    assert!(line.contains("\"label\":\"2026-07-24 baseline\""));
    assert!(line.contains("\"status\":\"ratio\",\"numerator\":13,\"denominator\":25"));
    assert!(line.contains("\"status\":\"no-data\""));
    assert!(
        !line.contains("trend"),
        "trends are derived, never recorded: {line}"
    );
}

/// G5: both renders and the identity are byte-identical across independent
/// builds from independently constructed inputs.
#[test]
fn dashboard_regenerates_byte_identically() {
    let build = || {
        let rows = vec![
            row("alpha", MetricDirection::HigherIsBetter, ratio(13, 25)),
            row("beta", MetricDirection::LowerIsBetter, count(15)),
            row(
                "gamma",
                MetricDirection::Neutral,
                no_data("no machinery", Some("f85xj.14.2")),
            ),
        ];
        let history = history_of(&[("alpha", ratio(10, 25))], ContentHash([3; 32]));
        let dashboard = build_dashboard(&rows, &history).expect("dashboard builds");
        (
            dashboard.render_markdown(),
            dashboard.render_json(),
            dashboard.identity(),
        )
    };
    let (markdown_a, json_a, identity_a) = build();
    let (markdown_b, json_b, identity_b) = build();
    assert_eq!(markdown_a, markdown_b);
    assert_eq!(json_a, json_b);
    assert_eq!(identity_a, identity_b);
}

/// Caller order cannot move the artifact: rows sort by `(family, id)`.
#[test]
fn row_order_is_independent_of_caller_order() {
    let governance = MetricRow::try_new(
        "zzz-governance",
        "governance row",
        MetricFamily::Governance,
        MetricDirection::HigherIsBetter,
        ratio(11, 15),
        &["capability-maturity.json"],
        "registry levels are declarations, not proofs",
    )
    .expect("row admits");
    let outcome = row("aaa-outcome", MetricDirection::HigherIsBetter, count(1));

    let forward = build_dashboard(
        &[outcome.clone(), governance.clone()],
        &MetricHistory::empty(),
    )
    .expect("dashboard builds");
    let reversed =
        build_dashboard(&[governance, outcome], &MetricHistory::empty()).expect("dashboard builds");

    assert_eq!(forward.render_json(), reversed.render_json());
    assert_eq!(forward.identity(), reversed.identity());
    // Outcome precedes governance regardless of how the caller supplied them.
    assert_eq!(forward.rows()[0].metric().id(), "aaa-outcome");
    assert_eq!(forward.rows()[1].metric().id(), "zzz-governance");
}

/// Fail-closed admission, compared against exact error values.
#[test]
fn rows_and_history_fail_closed() {
    let valid_sources: &[&str] = &["vv-scorecard.json"];

    // An unqualified number is the failure mode this artifact prevents, so a
    // blank honesty statement is refused outright.
    assert_eq!(
        MetricRow::try_new(
            "id",
            "title",
            MetricFamily::Outcome,
            MetricDirection::Neutral,
            count(1),
            valid_sources,
            "   ",
        ),
        Err(DashboardError::InvalidMetricField {
            field: "honesty",
            reason: "must not be blank",
        })
    );

    assert_eq!(
        MetricRow::try_new(
            "",
            "title",
            MetricFamily::Outcome,
            MetricDirection::Neutral,
            count(1),
            valid_sources,
            "honest",
        ),
        Err(DashboardError::InvalidMetricField {
            field: "id",
            reason: "must not be blank",
        })
    );

    assert_eq!(
        MetricRow::try_new(
            "id",
            "title\nsplit",
            MetricFamily::Outcome,
            MetricDirection::Neutral,
            count(1),
            valid_sources,
            "honest",
        ),
        Err(DashboardError::InvalidMetricField {
            field: "title",
            reason: "must not contain control characters",
        })
    );

    // A measurement nobody can trace is not evidence.
    assert_eq!(
        MetricRow::try_new(
            "untraceable",
            "title",
            MetricFamily::Outcome,
            MetricDirection::Neutral,
            count(1),
            &[],
            "honest",
        ),
        Err(DashboardError::MissingSourceCitation {
            metric: "untraceable".to_string(),
        })
    );

    // A NO-DATA row needs no source, because it reports no measurement.
    assert!(
        MetricRow::try_new(
            "absent",
            "title",
            MetricFamily::Outcome,
            MetricDirection::Neutral,
            no_data("not wired", None),
            &[],
            "honest",
        )
        .is_ok()
    );

    let duplicate = row("same", MetricDirection::Neutral, count(1));
    assert_eq!(
        build_dashboard(&[duplicate.clone(), duplicate], &MetricHistory::empty()),
        Err(DashboardError::DuplicateMetric {
            metric: "same".to_string(),
        })
    );

    let overflow = vec![row("bulk", MetricDirection::Neutral, count(1)); MAX_METRIC_ROWS + 1];
    assert_eq!(
        build_dashboard(&overflow, &MetricHistory::empty()),
        Err(DashboardError::ResourceLimit {
            limit: MAX_METRIC_ROWS,
            observed: MAX_METRIC_ROWS + 1,
        })
    );

    // A history that can reorder is not a trend basis.
    let generation = |number: u64| {
        HistoryGeneration::try_new(number, "label", ContentHash([1; 32]), BTreeMap::new())
            .expect("generation admits")
    };
    assert_eq!(
        MetricHistory::try_new(vec![generation(2), generation(2)]),
        Err(DashboardError::NonMonotonicHistory { generation: 2 })
    );
    assert_eq!(
        MetricHistory::try_new(vec![generation(3), generation(1)]),
        Err(DashboardError::NonMonotonicHistory { generation: 1 })
    );
    assert!(MetricHistory::try_new(vec![generation(1), generation(4)]).is_ok());

    let empty = MetricHistory::empty();
    assert!(empty.is_empty());
    assert_eq!(empty.next_generation(), 1);
    assert_eq!(
        MetricHistory::try_new(vec![generation(7)])
            .expect("history admits")
            .next_generation(),
        8
    );
}

/// Every row must carry an honesty statement, and the render must surface all
/// of them plus the deliberate exclusions — an absent metric is explained, not
/// hidden.
#[test]
fn honesty_and_exclusions_are_always_rendered() {
    let rows = vec![
        row("measured", MetricDirection::HigherIsBetter, ratio(13, 25)),
        row(
            "absent",
            MetricDirection::LowerIsBetter,
            no_data("no acceptance lane exists", Some("f85xj.6.11")),
        ),
    ];
    let dashboard = build_dashboard(&rows, &MetricHistory::empty()).expect("dashboard builds");
    let markdown = dashboard.render_markdown();

    assert!(markdown.contains("## What each metric does not capture"));
    for entry in dashboard.rows() {
        assert!(
            markdown.contains(entry.metric().honesty()),
            "row {} lost its honesty statement",
            entry.metric().id()
        );
    }

    assert!(markdown.contains("## Why metrics are missing"));
    assert!(markdown.contains("no acceptance lane exists"));
    assert!(markdown.contains("(tracked: f85xj.6.11)"));

    // Diagnostics are excluded on the record, with the reason.
    assert!(markdown.contains("## Deliberately excluded"));
    assert!(markdown.contains("kernel throughput"));
    assert!(markdown.contains("crate count"));
    assert!(dashboard.render_json().contains("\"excluded\":["));

    // The trend basis is itself NO-DATA before anything is recorded.
    assert!(markdown.contains("trend_basis: NO-DATA"));
    assert_eq!(dashboard.prior_generation(), None);
}

/// The real seeded registries project the real metric set. This pins the
/// current honest program state: several rows are genuinely zero, and that is
/// the point — they must be visible, not hidden behind `NO-DATA`.
#[test]
fn real_registries_project_the_program_metric_set() {
    use fs_report::program_metrics::{ProgramSources, frankensim_rows};
    use fs_vvreg::adversarial::adversarial_registry;
    use fs_vvreg::corpus::corpus;
    use fs_vvreg::scorecard::build_scorecard;

    let scorecard = build_scorecard(corpus(), adversarial_registry(), &[], &[])
        .expect("the seeded scorecard builds");
    let levels: BTreeMap<String, String> = [("a", "L1"), ("b", "L2"), ("c", "L3")]
        .iter()
        .map(|(id, level)| ((*id).to_string(), (*level).to_string()))
        .collect();
    let sources =
        ProgramSources::from_registries(corpus(), adversarial_registry(), &scorecard, &levels)
            .expect("the real registries project cleanly");

    assert!(sources.corpus_seeded, "the public corpus must be seeded");
    assert_eq!(sources.corpus_datasets, scorecard_dataset_total());
    assert_eq!(
        sources.claim_cells,
        sources.externally_anchored_cells + sources.unanchored_cells,
        "every claim cell is either anchored or unanchored"
    );
    assert!(
        sources.external_datasets <= sources.corpus_datasets,
        "external datasets are a subset"
    );
    assert_eq!(sources.capabilities, 3);
    assert_eq!(sources.capabilities_verified, 2, "L2 and L3 are both >= L2");
    assert_eq!(sources.capabilities_integrated, 1, "only L3 is >= L3");

    let rows = frankensim_rows(sources).expect("the metric set projects");
    let dashboard = build_dashboard(&rows, &MetricHistory::empty()).expect("dashboard builds");
    assert!(dashboard.measured_count() >= 5);
    assert!(dashboard.no_data_count() > 0);

    // A caller-built corpus may never masquerade as the public program state.
    let mut synthetic = sources;
    synthetic.corpus_seeded = false;
    assert_eq!(
        frankensim_rows(synthetic),
        Err(DashboardError::UnseededCorpus)
    );

    // An empty capability registry is a READ FAILURE, never a program with
    // nothing in it — the trap that would otherwise render "0 capabilities".
    assert_eq!(
        ProgramSources::count_levels(&BTreeMap::new()),
        Err(DashboardError::EmptySourceRegistry {
            source: "capability-maturity",
        })
    );
    let bogus: BTreeMap<String, String> = [("a".to_string(), "L9".to_string())].into();
    assert_eq!(
        ProgramSources::count_levels(&bogus),
        Err(DashboardError::UninterpretableSource {
            source: "capability-maturity",
            value: "L9".to_string(),
        })
    );
}

fn scorecard_dataset_total() -> usize {
    fs_vvreg::corpus::corpus().datasets().len()
}
