//! Three-color schema conformance (bead qmao.1; CONTRACT.md addendum):
//! composition totality (G0), regime-exit auto-demotion, the adversarial
//! LAUNDERING gauntlet (G3, a Certifying-the-Certifiers gate), the
//! waiver-in-provenance path, the fs-evidence bridge, and determinism.
//! JSON-line verdicts; seeded cases carry seeds.

use fs_evidence::{
    Color, ColorError, ColorRank, IntervalOp, ModelEvidence, NumericalCertificate, ValidityDomain,
    check_regime, color_of, compose, verified_from,
};
use fs_ledger::{ColorGraph, ColorWriteError, Waiver};
use std::collections::BTreeMap;

fn verdict(case: &str, pass: bool, detail: &str) {
    println!(
        "{{\"suite\":\"fs-ledger/colors\",\"case\":\"{case}\",\"verdict\":\"{}\",\
         \"detail\":\"{detail}\"}}",
        if pass { "pass" } else { "fail" }
    );
    assert!(pass, "case {case}: {detail}");
}

struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0
    }

    fn unit(&mut self) -> f64 {
        ((self.next() >> 11) as f64) / (1u64 << 53) as f64
    }

    fn below(&mut self, n: u64) -> u64 {
        (self.next() >> 32) % n
    }
}

fn rand_color(rng: &mut Lcg) -> Color {
    match rng.below(3) {
        0 => {
            let lo = rng.unit() * 2.0 - 1.0;
            Color::Verified {
                lo,
                hi: lo + rng.unit(),
            }
        }
        1 => Color::Validated {
            regime: ValidityDomain::unconstrained().with(
                "reynolds",
                1e3 * (1.0 + rng.unit()),
                1e5 * (1.0 + rng.unit()),
            ),
            dataset: format!("ds-{}", rng.below(4)),
        },
        _ => Color::Estimated {
            estimator: format!("est-{}", rng.below(4)),
            dispersion: rng.unit(),
        },
    }
}

/// col-001 — G0 totality: every color pair composes to a defined,
/// CONSERVATIVE result (rank = min of operand ranks, never higher);
/// verified bounds compose per op; validated regimes intersect;
/// estimated dispersions add.
#[test]
fn col_001_composition_totality() {
    let mut rng = Lcg(0x1001_2026_0707_0041);
    let mut total_ok = true;
    let mut conservative_ok = true;
    for _ in 0..600 {
        let (a, b) = (rand_color(&mut rng), rand_color(&mut rng));
        let op = match rng.below(3) {
            0 => IntervalOp::Add,
            1 => IntervalOp::Mul,
            _ => IntervalOp::Hull,
        };
        let out = compose(&a, &b, op);
        conservative_ok &= out.rank() <= a.rank().min(b.rank());
        total_ok &= matches!(
            out,
            Color::Verified { .. } | Color::Validated { .. } | Color::Estimated { .. }
        );
    }
    // Verified interval arithmetic spot checks.
    let v = |lo: f64, hi: f64| Color::Verified { lo, hi };
    let add = compose(&v(1.0, 2.0), &v(10.0, 20.0), IntervalOp::Add);
    let mul = compose(&v(-1.0, 2.0), &v(3.0, 4.0), IntervalOp::Mul);
    let verified_math = add == v(11.0, 22.0) && mul == v(-4.0, 8.0);
    // Regime intersection: both anchors must hold.
    let val = |lo: f64, hi: f64| Color::Validated {
        regime: ValidityDomain::unconstrained().with("re", lo, hi),
        dataset: "wind-tunnel-a".to_string(),
    };
    let both = compose(&val(1e3, 1e5), &val(5e3, 5e5), IntervalOp::Add);
    let intersected = matches!(&both, Color::Validated { regime, .. }
        if regime.bounds()["re"] == (5e3, 1e5));
    // Estimated absorbs everything, dispersion conservatively.
    let est = compose(
        &Color::Estimated {
            estimator: "koopman".to_string(),
            dispersion: 0.1,
        },
        &v(0.0, 1.0),
        IntervalOp::Add,
    );
    let absorbed = matches!(est, Color::Estimated { dispersion, .. } if dispersion >= 0.1);
    verdict(
        "col-001",
        total_ok && conservative_ok && verified_math && intersected && absorbed,
        "all 600 random pairs compose totally with rank = min (never higher), \
         verified interval add/mul are exact, validated regimes INTERSECT, and \
         estimated absorbs everything with conservative dispersion; \
         seed 0x1001_2026_0707_0041",
    );
}

/// col-002 — regime-exit AUTO-DEMOTION: validated survives inside its
/// regime, demotes to estimated (with the flag naming dataset, axis,
/// value) the moment any axis exits or goes unreported.
#[test]
fn col_002_regime_demotion() {
    let validated = Color::Validated {
        regime: ValidityDomain::unconstrained()
            .with("reynolds", 1e3, 1e5)
            .with("mach", 0.0, 0.3),
        dataset: "wind-tunnel-a".to_string(),
    };
    let inside: BTreeMap<String, f64> =
        [("reynolds".to_string(), 5e4), ("mach".to_string(), 0.2)].into();
    let (c1, d1) = check_regime(&validated, &inside);
    let stays = c1 == validated && d1.is_none();

    let outside: BTreeMap<String, f64> =
        [("reynolds".to_string(), 2e5), ("mach".to_string(), 0.2)].into();
    let (c2, d2) = check_regime(&validated, &outside);
    let demotes = matches!(&c2, Color::Estimated { estimator, dispersion }
        if estimator.contains("regime-exit") && estimator.contains("wind-tunnel-a")
            && dispersion.is_infinite())
        && d2
            .as_ref()
            .is_some_and(|d| d.axis == "reynolds" && (d.value - 2e5).abs() < 1.0);

    let unreported: BTreeMap<String, f64> = [("reynolds".to_string(), 5e4)].into();
    let (c3, d3) = check_regime(&validated, &unreported);
    let unreported_demotes =
        matches!(c3, Color::Estimated { .. }) && d3.is_some_and(|d| d.axis == "mach");

    // Verified and estimated pass through untouched.
    let v = Color::Verified { lo: 0.0, hi: 1.0 };
    let (cv, dv) = check_regime(&v, &outside);
    let passthrough = cv == v && dv.is_none();

    verdict(
        "col-002",
        stays && demotes && unreported_demotes && passthrough,
        "validated survives inside its regime, AUTO-DEMOTES to estimated (infinite \
         dispersion, flag naming wind-tunnel-a/reynolds/2e5) on exit, demotes on an \
         UNREPORTED axis, and verified/estimated pass through untouched",
    );
}

/// col-003 — the LAUNDERING gauntlet (G3, security-critical): every
/// adversarial attempt to upgrade a color fails the type check with
/// the capping parents named.
#[test]
#[allow(clippy::too_many_lines)] // the gauntlet's five doors are one story
fn col_003_laundering_gauntlet() {
    // (a) The constructor door: estimate/no-claim certs refuse verified.
    let est_cert = NumericalCertificate::estimate(0.0, 1.0);
    let door = verified_from(&est_cert);
    let door_refuses = matches!(&door, Err(ColorError::LaunderingRefused { actual }) if *actual == "estimate")
        && door.unwrap_err().to_string().contains("waiver");
    let no_claim = verified_from(&NumericalCertificate::no_claim());
    let no_claim_refuses = no_claim.is_err();
    // Enclosures pass the door.
    let ok_door = verified_from(&NumericalCertificate::enclosure(0.0, 0.5)).is_ok();

    // (b) The write gate: estimated parent caps every claim.
    let state = BTreeMap::new();
    let mut g = ColorGraph::new();
    let clean = g.source("fem-bound", Color::Verified { lo: 0.9, hi: 1.1 });
    let dirty = g.source(
        "surrogate-drag",
        Color::Estimated {
            estimator: "deeponet-v2".to_string(),
            dispersion: 0.05,
        },
    );
    let attempt = g.derive(
        "total-drag",
        &[clean, dirty],
        IntervalOp::Add,
        Some(Color::Verified { lo: 0.0, hi: 2.0 }),
        &state,
        None,
    );
    let gate_refuses = matches!(&attempt,
        Err(ColorWriteError::LaunderingRefused { claimed: ColorRank::Verified,
            derived: ColorRank::Estimated, offending_parents })
            if offending_parents.contains(&dirty));

    // (c) Claiming validated over estimated parents refuses too.
    let attempt2 = g.derive(
        "calibrated-drag",
        &[dirty],
        IntervalOp::Hull,
        Some(Color::Validated {
            regime: ValidityDomain::unconstrained().with("re", 1e3, 1e5),
            dataset: "wishful".to_string(),
        }),
        &state,
        None,
    );
    let gate_refuses2 = attempt2.is_err();

    // (d) Post-demotion upgrade attempts refuse: a validated parent
    // whose regime the state has exited caps at estimated.
    let val = g.source(
        "turbulence-closure",
        Color::Validated {
            regime: ValidityDomain::unconstrained().with("reynolds", 1e3, 1e5),
            dataset: "wind-tunnel-a".to_string(),
        },
    );
    let exited: BTreeMap<String, f64> = [("reynolds".to_string(), 9e5)].into();
    let attempt3 = g.derive(
        "lift-coefficient",
        &[val],
        IntervalOp::Hull,
        Some(Color::Validated {
            regime: ValidityDomain::unconstrained().with("reynolds", 1e3, 1e6),
            dataset: "wind-tunnel-a".to_string(),
        }),
        &exited,
        None,
    );
    let post_demotion_refuses = attempt3.is_err();

    // (e) Seeded adversarial pipelines: random DAGs, random upgrade
    // attempts — every single one must refuse.
    let mut rng = Lcg(0x1001_2026_0707_0043);
    let mut attempts = 0u32;
    let mut refusals = 0u32;
    for _ in 0..60 {
        let mut gg = ColorGraph::new();
        let mut ids = Vec::new();
        for k in 0..5 {
            ids.push(gg.source(&format!("s{k}"), rand_color(&mut rng)));
        }
        for k in 0..6 {
            let a = ids[rng.below(ids.len() as u64) as usize];
            let b = ids[rng.below(ids.len() as u64) as usize];
            let derived_rank = {
                let (ca, _) = check_regime(&gg.node(a).color.clone(), &state);
                let (cb, _) = check_regime(&gg.node(b).color.clone(), &state);
                compose(&ca, &cb, IntervalOp::Hull).rank()
            };
            // Claim strictly ABOVE what the parents support.
            let claim = match derived_rank {
                ColorRank::Verified => None, // nothing to launder
                ColorRank::Validated | ColorRank::Estimated => {
                    Some(Color::Verified { lo: 0.0, hi: 1.0 })
                }
            };
            if let Some(c) = claim {
                attempts += 1;
                if gg
                    .derive(
                        &format!("d{k}"),
                        &[a, b],
                        IntervalOp::Hull,
                        Some(c),
                        &state,
                        None,
                    )
                    .is_err()
                {
                    refusals += 1;
                }
            } else if let Ok(id) = gg.derive(
                &format!("d{k}"),
                &[a, b],
                IntervalOp::Hull,
                None,
                &state,
                None,
            ) {
                ids.push(id);
            }
        }
    }
    verdict(
        "col-003",
        door_refuses
            && no_claim_refuses
            && ok_door
            && gate_refuses
            && gate_refuses2
            && post_demotion_refuses
            && attempts > 100
            && refusals == attempts,
        &format!(
            "every laundering path refuses: the constructor door (estimate/no-claim \
             -> verified), the write gate (estimated parent caps all claims, capping \
             parents NAMED), validated-over-estimated, post-demotion re-claims, and \
             {refusals}/{attempts} seeded adversarial upgrade attempts; \
             seed 0x1001_2026_0707_0043"
        ),
    );
}

/// col-004 — the waiver path: a signed waiver authorizes the upgrade,
/// appears in the ledger row, AND participates in the provenance hash
/// (dropping it changes the hash — it cannot vanish quietly).
#[test]
fn col_004_waiver_in_provenance() {
    let state = BTreeMap::new();
    let mut g = ColorGraph::new();
    let dirty = g.source(
        "surrogate",
        Color::Estimated {
            estimator: "pod-deim".to_string(),
            dispersion: 0.02,
        },
    );
    let waiver = Waiver {
        id: "WVR-2026-041".to_string(),
        signer: "chief-engineer".to_string(),
        reason: "surrogate validated offline against holdout campaign 7".to_string(),
    };
    let waived = g.derive(
        "release-metric",
        &[dirty],
        IntervalOp::Hull,
        Some(Color::Verified { lo: 0.0, hi: 1.0 }),
        &state,
        Some(waiver.clone()),
    );
    let succeeded = waived.is_ok();
    let id = waived.expect("waived write");
    let node = g.node(id);
    let in_row = g
        .rows()
        .iter()
        .any(|r| r.contains("WVR-2026-041") && r.contains("chief-engineer"));
    // Provenance: an identical write WITHOUT the waiver (in a fresh
    // graph, forced through by claiming only what parents support)
    // hashes DIFFERENTLY.
    let mut g2 = ColorGraph::new();
    let dirty2 = g2.source(
        "surrogate",
        Color::Estimated {
            estimator: "pod-deim".to_string(),
            dispersion: 0.02,
        },
    );
    let _ = dirty2;
    let plain = g2
        .derive(
            "release-metric",
            &[dirty2],
            IntervalOp::Hull,
            None,
            &state,
            None,
        )
        .expect("plain");
    let hash_differs = g2.node(plain).hash.to_hex() != node.hash.to_hex();
    verdict(
        "col-004",
        succeeded && node.waiver.is_some() && in_row && hash_differs,
        "the signed waiver authorizes the upgrade, appears verbatim in the ledger \
         row, and participates in the provenance hash (the same write without it \
         hashes differently — waivers cannot vanish quietly)",
    );
}

/// col-005 — the fs-evidence bridge: existing receipts color honestly
/// (enclosure -> verified; carded model with bounded validity ->
/// validated with THAT regime; estimates -> estimated).
#[test]
fn col_005_evidence_bridge() {
    let verified = color_of(
        &NumericalCertificate::enclosure(0.9, 1.1),
        &ModelEvidence::none(),
    );
    // Bounds pass through by BITS (no arithmetic on this path).
    let v_ok = matches!(verified, Color::Verified { lo, hi }
        if lo.to_bits() == 0.9f64.to_bits() && hi.to_bits() == 1.1f64.to_bits());
    let validated = color_of(
        &NumericalCertificate::estimate(0.0, 1.0),
        &ModelEvidence {
            cards: vec!["k-epsilon".to_string()],
            assumptions: vec![],
            validity: ValidityDomain::unconstrained().with("reynolds", 1e3, 1e5),
            discrepancy_rel: 0.03,
            in_domain: true,
        },
    );
    let val_ok = matches!(&validated, Color::Validated { regime, dataset }
        if dataset == "k-epsilon" && regime.bounds().contains_key("reynolds"));
    let estimated = color_of(
        &NumericalCertificate::estimate(0.0, 1.0),
        &ModelEvidence::none(),
    );
    let est_ok = matches!(estimated, Color::Estimated { .. });
    verdict(
        "col-005",
        v_ok && val_ok && est_ok,
        "existing fs-evidence receipts color honestly: enclosures become verified \
         with their bounds, carded models with bounded validity become validated \
         with THAT regime, uncarded estimates stay estimated",
    );
}

/// col-006 — determinism: identical write sequences give bitwise
/// identical rows and hashes.
#[test]
fn col_006_determinism() {
    let build = || -> Vec<String> {
        let state: BTreeMap<String, f64> = [("re".to_string(), 5e4)].into();
        let mut g = ColorGraph::new();
        let a = g.source("a", Color::Verified { lo: 0.0, hi: 1.0 });
        let b = g.source(
            "b",
            Color::Validated {
                regime: ValidityDomain::unconstrained().with("re", 1e3, 1e5),
                dataset: "ds".to_string(),
            },
        );
        let c = g
            .derive("c", &[a, b], IntervalOp::Add, None, &state, None)
            .expect("c");
        let _ = c;
        g.rows().to_vec()
    };
    let (r1, r2) = (build(), build());
    verdict(
        "col-006",
        r1 == r2 && !r1.is_empty(),
        "identical write sequences produce bitwise-identical rows and provenance \
         hashes",
    );
}

/// v3 migration regression (bead lmp4.3): the speculation extension
/// table exists, round-trips the four solve-node fields, and every
/// pre-existing table still answers queries (nothing broke).
#[test]
fn speculation_schema_migration() {
    let ledger = fs_ledger::Ledger::open(":memory:").expect("open");
    assert_eq!(ledger.schema_version().expect("version"), 3);
    let body = "{\"proposer_id\":\"neighbor-extrapolation\",\"accepted\":true,\
                \"bound\":3.2e-4,\"iterations_saved\":4}";
    ledger
        .put_extension(fs_ledger::ExtensionTable::Speculation, "solve-op-17", body)
        .expect("put");
    let back = ledger
        .get_extension(fs_ledger::ExtensionTable::Speculation, "solve-op-17")
        .expect("get")
        .expect("present");
    assert!(back.contains("iterations_saved"), "{back}");
    // Existing tables unbroken.
    for table in fs_ledger::ALL_TABLES {
        let _ = ledger.table_count(table).expect("old queries still work");
    }
}
