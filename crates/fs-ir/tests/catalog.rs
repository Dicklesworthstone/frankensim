//! G0 battery for the operator catalog (bead gp3.6): generation
//! determinism, registry/dispatch drift in both directions, executable
//! examples (documentation rots loudly), deliberate-mismatch refusal,
//! structured query fixtures, signature-compatibility search,
//! error-model reference resolution, diffability, and the
//! serve-latency budget (strict concurrent-load variant gated for
//! quiet perf lanes).

use fs_ir::catalog::{
    ARITH_SAME_DIMS, COMPARE_FORMS, Catalog, CatalogQuery, OperatorEntry, OperatorKind,
    SUGAR_VERBS, SigType,
};
use fs_ir::lower::lower;
use fs_ir::sexpr;

fn head_symbol(node: &fs_ir::Node) -> Option<&str> {
    match &node.items()?.first()?.kind {
        fs_ir::NodeKind::Symbol(name) => Some(name.as_str()),
        _ => None,
    }
}

fn parse(program: &str) -> fs_ir::Node {
    sexpr::parse(program).expect("catalog example must parse")
}

#[test]
fn cat_001_generation_is_deterministic_and_validates() {
    let first = Catalog::builtin();
    let second = Catalog::builtin();
    assert_eq!(
        first.to_canonical_jsonl(),
        second.to_canonical_jsonl(),
        "two generations must export byte-identical canonical JSONL"
    );
    first.validate().expect("the built-in catalog validates");
    assert!(
        first.entries().len() >= SUGAR_VERBS.len() + ARITH_SAME_DIMS.len() + COMPARE_FORMS.len(),
        "catalog covers every registered surface"
    );
}

#[test]
fn cat_002_every_sugar_verb_lowers_and_dispatch_knows_no_others() {
    let catalog = Catalog::builtin();
    for entry in catalog.entries() {
        let OperatorKind::SugarVerb {
            lowers_to,
            injected_defaults,
        } = &entry.kind
        else {
            continue;
        };
        // The entry's FIRST example is the minimal invocation: lowering
        // it must succeed, expand to the declared target, and inject
        // exactly the declared defaults — executed drift proof.
        let example = entry.examples.first().expect("sugar entry has an example");
        let lowered = lower(&parse(example))
            .unwrap_or_else(|error| panic!("{} example must lower: {error:?}", entry.name));
        let trace = lowered.trace;
        let head = head_symbol(&lowered.node).unwrap_or_default();
        assert_eq!(
            head, *lowers_to,
            "{}: catalog lowers_to and the real expansion disagree",
            entry.name
        );
        let step = trace
            .iter()
            .find(|step| step.verb == entry.name)
            .unwrap_or_else(|| panic!("{}: lowering recorded no trace step", entry.name));
        assert_eq!(
            step.injected, *injected_defaults,
            "{}: catalog injected_defaults and the real trace disagree",
            entry.name
        );
    }
    // The reverse direction: a head the registry does not know must
    // pass through lower untouched (no hidden sugar dispatch).
    let untouched =
        lower(&parse("(not-a-registered-verb 1 2)")).expect("unknown heads pass through");
    assert_eq!(
        head_symbol(&untouched.node),
        Some("not-a-registered-verb"),
        "unknown heads must not be expanded"
    );
    assert!(
        untouched.trace.is_empty(),
        "no trace step for unregistered heads"
    );
}

#[test]
fn cat_003_every_example_parses() {
    let catalog = Catalog::builtin();
    for entry in catalog.entries() {
        for example in entry.examples {
            let parsed = parse(example);
            assert!(
                parsed.items().is_some(),
                "{}: example {example:?} must be a form",
                entry.name
            );
        }
    }
}

#[test]
fn cat_004_deliberate_mismatch_fails_validation() {
    // A sugar entry whose target is not cataloged: the drift lint must
    // refuse — this is the "deliberately mismatched hand-annotation
    // fails the build" acceptance criterion exercised at the unit
    // level (the build-level enforcement is cat_001 + cat_002 running
    // in CI).
    let catalog = Catalog::builtin();
    let mut entries: Vec<OperatorEntry> = catalog.entries().to_vec();
    for entry in &mut entries {
        if let OperatorKind::SugarVerb { lowers_to, .. } = &mut entry.kind {
            *lowers_to = "engine.that-does-not-exist";
            break;
        }
    }
    let tampered = Catalog::from_entries_for_test(entries);
    let drifts = tampered
        .validate()
        .expect_err("a mismatched lowers_to must fail validation");
    assert!(
        drifts
            .iter()
            .any(|d| d.detail.contains("engine.that-does-not-exist")),
        "the drift names the mismatched target: {drifts:?}"
    );
}

#[test]
fn cat_005_structured_queries_return_correct_sets() {
    let catalog = Catalog::builtin();

    let sugar = catalog.query(&CatalogQuery {
        namespace: Some("sugar".to_string()),
        ..CatalogQuery::default()
    });
    assert_eq!(
        sugar.iter().map(|e| e.name).collect::<Vec<_>>(),
        SUGAR_VERBS.to_vec(),
        "namespace=sugar returns exactly the registered sugar verbs"
    );

    let ascent = catalog.query(&CatalogQuery {
        name: Some("ascent.*".to_string()),
        ..CatalogQuery::default()
    });
    assert_eq!(ascent.len(), 1);
    assert_eq!(ascent[0].name, "ascent.optimize");

    let granted = catalog.query(&CatalogQuery {
        granted_capabilities: Some(vec!["flux.*".to_string()]),
        ..CatalogQuery::default()
    });
    assert!(
        granted.iter().any(|e| e.name == "flux.free-surface-lbm"),
        "flux.* grants admit the flux operator"
    );
    assert!(
        granted.iter().all(|e| e.name != "ascent.optimize"),
        "flux.* grants do not admit ascent.optimize"
    );

    let probabilistic = catalog.query(&CatalogQuery {
        namespace: Some("query".to_string()),
        ..CatalogQuery::default()
    });
    let exceedance = probabilistic
        .iter()
        .find(|e| e.name == "exceedance")
        .expect("exceedance cataloged");
    assert!(
        matches!(
            exceedance.kind,
            OperatorKind::QoiOperator {
                probabilistic: true,
                linear: false,
                ..
            }
        ),
        "exceedance metadata derives from the live Qoi menu"
    );

    let text = catalog.query(&CatalogQuery {
        text: Some("free-surface".to_string()),
        ..CatalogQuery::default()
    });
    assert!(text.iter().any(|e| e.name == "flux.free-surface-lbm"));
}

#[test]
fn cat_006_diff_reports_added_removed_changed() {
    let catalog = Catalog::builtin();
    let current = catalog.to_canonical_jsonl();

    let unchanged = catalog.diff(&current);
    assert!(unchanged.added.is_empty());
    assert!(unchanged.removed.is_empty());
    assert!(unchanged.changed.is_empty());

    // Remove one line -> reported as added (present now, absent then);
    // tamper one line -> reported as changed.
    let without_integral: String = current
        .lines()
        .filter(|line| !line.contains("\"name\":\"integral\""))
        .map(|line| format!("{line}\n"))
        .collect();
    let diff = catalog.diff(&without_integral);
    assert_eq!(diff.added, vec!["integral".to_string()]);
    assert!(diff.removed.is_empty());

    let tampered: String = current
        .lines()
        .map(|line| {
            if line.contains("\"name\":\"integral\"") {
                format!("{}\n", line.replace("\"linear\":true", "\"linear\":false"))
            } else {
                format!("{line}\n")
            }
        })
        .collect();
    let diff = catalog.diff(&tampered);
    assert_eq!(diff.changed, vec!["integral".to_string()]);
}

fn names(entries: &[&OperatorEntry]) -> Vec<&'static str> {
    entries.iter().map(|e| e.name).collect()
}

#[test]
fn cat_009_signature_compatibility_search() {
    let catalog = Catalog::builtin();

    // "What consumes a field?" — exactly the QoI menu.
    let field_consumers = catalog.query(&CatalogQuery {
        consumes: Some(vec!["field".to_string()]),
        ..CatalogQuery::default()
    });
    assert_eq!(
        names(&field_consumers),
        vec!["exceedance", "integral", "max-over-region"],
        "field consumers are the QoI operators, in canonical order"
    );

    // "What consumes a field and produces a probability?" — the bead's
    // canonical fixture: exactly the probabilistic QoI.
    let probabilistic = catalog.query(&CatalogQuery {
        consumes: Some(vec!["field".to_string()]),
        produces: Some("probability".to_string()),
        ..CatalogQuery::default()
    });
    assert_eq!(names(&probabilistic), vec!["exceedance"]);

    // "What gives me a campaign?" — the two core engines.
    let campaigns = catalog.query(&CatalogQuery {
        produces: Some("campaign".to_string()),
        ..CatalogQuery::default()
    });
    assert_eq!(
        names(&campaigns),
        vec!["ascent.optimize", "flux.free-surface-lbm"]
    );

    // "Quantity in, verdict out" — exactly the comparison forms.
    let verdicts = catalog.query(&CatalogQuery {
        consumes: Some(vec!["quantity".to_string()]),
        produces: Some("verdict".to_string()),
        ..CatalogQuery::default()
    });
    assert_eq!(names(&verdicts), vec!["<", "<=", "=", ">", ">="]);

    // Multi-token consumes is conjunctive over the entry's signature.
    let pours = catalog.query(&CatalogQuery {
        consumes: Some(vec!["vessel".to_string(), "schedule".to_string()]),
        ..CatalogQuery::default()
    });
    assert_eq!(
        names(&pours),
        vec!["flux.free-surface-lbm", "simulate-pour"]
    );

    // Signature filters compose conjunctively with the other filters.
    let granted_campaigns = catalog.query(&CatalogQuery {
        produces: Some("campaign".to_string()),
        granted_capabilities: Some(vec!["flux.*".to_string()]),
        ..CatalogQuery::default()
    });
    assert_eq!(names(&granted_campaigns), vec!["flux.free-surface-lbm"]);

    // A token outside the vocabulary matches nothing (fail-safe), in
    // either position.
    for query in [
        CatalogQuery {
            consumes: Some(vec!["warp-core".to_string()]),
            ..CatalogQuery::default()
        },
        CatalogQuery {
            produces: Some("warp-core".to_string()),
            ..CatalogQuery::default()
        },
    ] {
        assert!(
            catalog.query(&query).is_empty(),
            "out-of-vocabulary tokens must return the empty set"
        );
    }
}

#[test]
fn cat_010_signature_drift_fails_validation() {
    let catalog = Catalog::builtin();

    // A sugar verb claiming to produce anything but an IR form is a
    // structural lie about what expansion does.
    let mut entries: Vec<OperatorEntry> = catalog.entries().to_vec();
    for entry in &mut entries {
        if matches!(entry.kind, OperatorKind::SugarVerb { .. }) {
            entry.produces = SigType::Probability;
            break;
        }
    }
    let drifts = Catalog::from_entries_for_test(entries)
        .validate()
        .expect_err("a sugar verb not producing ir-form must fail validation");
    assert!(
        drifts.iter().any(|d| d.detail.contains("ir-form")),
        "the drift names the invariant: {drifts:?}"
    );

    // A probabilistic QoI must produce a probability — the signature
    // column and the planner metadata may not disagree.
    let mut entries: Vec<OperatorEntry> = catalog.entries().to_vec();
    for entry in &mut entries {
        if matches!(
            entry.kind,
            OperatorKind::QoiOperator {
                probabilistic: true,
                ..
            }
        ) {
            entry.produces = SigType::QoiValue;
            break;
        }
    }
    let drifts = Catalog::from_entries_for_test(entries)
        .validate()
        .expect_err("probabilistic/produces disagreement must fail validation");
    assert!(
        drifts.iter().any(|d| d.detail.contains("probabilistic")),
        "the drift names the disagreement: {drifts:?}"
    );
}

#[test]
fn cat_007_query_latency_smoke() {
    // Generous default-suite bound (shared hosts are noisy); the
    // strict ≤100ms certification lives in cat_008 for quiet lanes.
    let catalog = Catalog::builtin();
    let queries = latency_query_set();
    let start = std::time::Instant::now();
    let mut total = 0usize;
    for query in &queries {
        total += catalog.query(query).len();
    }
    let elapsed = start.elapsed();
    assert!(total > 0, "the latency sweep exercises real matches");
    assert!(
        elapsed.as_millis() < 1000,
        "catalog query sweep took {elapsed:?} (default-suite smoke bound)"
    );
}

#[test]
#[ignore = "quiet-lane certification of the plan §11.5 ≤100ms serve budget; run with --ignored on an unloaded host"]
fn cat_008_query_latency_certification() {
    // The acceptance budget holds UNDER CONCURRENT SWEEP LOAD: eight
    // threads hammer the same catalog and every individual sweep on
    // every thread must serve within 100ms — the worst observed sweep
    // is the certified number, not a per-thread average.
    let catalog = Catalog::builtin();
    let queries = latency_query_set();
    std::thread::scope(|scope| {
        for _ in 0..8 {
            scope.spawn(|| {
                let mut worst = std::time::Duration::ZERO;
                for _ in 0..100 {
                    let start = std::time::Instant::now();
                    for query in &queries {
                        std::hint::black_box(catalog.query(query));
                    }
                    worst = worst.max(start.elapsed());
                }
                assert!(
                    worst.as_millis() <= 100,
                    "worst full-catalog sweep {worst:?} under 8-thread contention \
                     exceeds the 100ms serve budget"
                );
            });
        }
    });
}

#[test]
fn cat_011_error_model_references_resolve() {
    // Every entry's error-model reference must name a real crate whose
    // CONTRACT.md carries an `## Error model` section — the mechanical
    // no-dangling-reference check (references are data, so they rot
    // loudly here instead of silently in prose).
    let catalog = Catalog::builtin();
    let crates_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("fs-ir lives under crates/");
    for entry in catalog.entries() {
        let contract = crates_dir.join(entry.error_model).join("CONTRACT.md");
        let text = std::fs::read_to_string(&contract).unwrap_or_else(|error| {
            panic!(
                "{}: error_model {:?} has no readable CONTRACT.md at {}: {error}",
                entry.name,
                entry.error_model,
                contract.display()
            )
        });
        assert!(
            text.contains("## Error model"),
            "{}: {}/CONTRACT.md lacks an `## Error model` section",
            entry.name,
            entry.error_model
        );
    }
}

fn latency_query_set() -> Vec<CatalogQuery> {
    vec![
        CatalogQuery::default(),
        CatalogQuery {
            namespace: Some("sugar".to_string()),
            ..CatalogQuery::default()
        },
        CatalogQuery {
            namespace: Some("query".to_string()),
            ..CatalogQuery::default()
        },
        CatalogQuery {
            name: Some("flux.*".to_string()),
            ..CatalogQuery::default()
        },
        CatalogQuery {
            granted_capabilities: Some(vec!["ascent.*".to_string(), "flux.*".to_string()]),
            ..CatalogQuery::default()
        },
        CatalogQuery {
            text: Some("optimize".to_string()),
            ..CatalogQuery::default()
        },
        CatalogQuery {
            consumes: Some(vec!["field".to_string()]),
            produces: Some("probability".to_string()),
            ..CatalogQuery::default()
        },
    ]
}
