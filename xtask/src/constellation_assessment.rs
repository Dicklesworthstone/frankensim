//! Deterministic trust-cone assessment for FrankenSim's seven pinned siblings.
//!
//! Governance judgments are explicit below, while the usage map is derived
//! from workspace manifests and tokenized Rust sources. This prevents the
//! rendered assessment from laundering the plan's intended integrations into
//! claims about code that actually exists.

use super::{LockRow, Violation};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;
use std::path::Path;

pub(crate) const CHECK: &str = "constellation-trust-cone";
const DATA_PATH: &str = "constellation-trust-assessment.json";
const DOC_PATH: &str = "docs/CONSTELLATION_TRUST_CONE.md";
const SCHEMA: &str = "frankensim-constellation-trust-assessment-v1";
const IDENTITY_DOMAIN: &str = "org.frankensim.xtask.constellation-trust-assessment.v1";
const BEAD_ID: &str = "frankensim-extreal-program-f85xj.13.1";
const MAX_SAMPLES_PER_SIBLING: usize = 5;

#[derive(Clone, Copy)]
struct Policy {
    lib: &'static str,
    role: &'static str,
    layers: &'static [&'static str],
    claims: &'static [&'static str],
    correctness_risk: &'static str,
    availability_risk: &'static str,
    security_surface: &'static str,
    review_status: &'static str,
    verification: &'static [&'static str],
    gaps: &'static [&'static str],
    priority: u8,
    review_targets: &'static [&'static str],
}

const POLICIES: &[Policy] = &[
    Policy {
        lib: "asupersync",
        role: "Structured-concurrency, cancellation, scope, clock, and latency-lane substrate.",
        layers: &["L0 SUBSTRATE", "all cancellable kernels through fs-exec"],
        claims: &[
            "bounded cancellation and request-drain-finalize semantics",
            "task-scoped capability and budget propagation",
            "deterministic pause, drain, and replay boundaries",
        ],
        correctness_risk: "critical",
        availability_risk: "high",
        security_surface: "Capability, deadline, cancellation, executor, and fault-propagation inputs; no FrankenSim network parser is delegated to it.",
        review_status: "Pinned and exercised; independent review of load-bearing protocols is pending.",
        verification: &[
            "constellation.lock pin and clean-tree verification",
            "fs-exec unit, conformance, cancellation, race, and constellation-smoke batteries",
            "consumer tests construct real asupersync Cx and Budget values",
        ],
        gaps: &[
            "asupersync model checks do not prove FrankenSim TilePool, arena-lease, parked-crew, or injected-fault protocols",
            "no independently adjudicated cancellation-latency and drain proof spans both repositories",
        ],
        priority: 1,
        review_targets: &[
            "scope cancellation",
            "deadline propagation",
            "runtime shutdown",
            "fs-exec adapter boundary",
        ],
    },
    Policy {
        lib: "franken_networkx",
        role: "Graph representation and compatibility substrate for topology, voxel, sparse, and task-graph surfaces.",
        layers: &["L1 BEDROCK", "L2 MORPH", "L3/L4 graph consumers"],
        claims: &[
            "ground-structure and lattice graph semantics",
            "graph interchange compatibility",
            "planned tropical critical-path analytics",
        ],
        correctness_risk: "high",
        availability_risk: "medium",
        security_surface: "Graph keys, attributes, compatibility modes, and serialized graph-shaped inputs.",
        review_status: "Pinned with active production and optional interop call sites; cross-repo compatibility coverage is partial.",
        verification: &[
            "constellation.lock pin and clean-tree verification",
            "fs-truss and fs-rep-voxel graph construction tests",
            "fs-sparse optional fnx interoperability tests",
        ],
        gaps: &[
            "no release-train gate proves every active graph consumer against a proposed sibling bump",
            "planned critical-path analytics are not evidence for current runtime use",
        ],
        priority: 3,
        review_targets: &[
            "Graph identity",
            "attribute conversion",
            "compatibility-mode behavior",
        ],
    },
    Policy {
        lib: "franken_numpy",
        role: "Array dtype and ufunc interoperability membrane; intended zero-copy results exchange.",
        layers: &["L1 BEDROCK interop", "planned HELM and notebook membranes"],
        claims: &[
            "optional sparse-array interchange",
            "planned zero-copy cochain and lattice views",
        ],
        correctness_risk: "high",
        availability_risk: "medium",
        security_surface: "Shape, stride, dtype, ownership, and buffer-boundary metadata.",
        review_status: "Pinned; measured use is narrower than the plan-level zero-copy contract.",
        verification: &[
            "constellation.lock pin and clean-tree verification",
            "fs-sparse optional fnp dtype and ufunc conversion tests",
        ],
        gaps: &[
            "no repository-wide zero-copy cochain or lattice-view conformance suite exists",
            "aliasing, lifetime, stride, and device-boundary compatibility are not independently reviewed",
        ],
        priority: 2,
        review_targets: &[
            "dtype mapping",
            "shape/stride preservation",
            "ownership and zero-copy claims",
        ],
    },
    Policy {
        lib: "frankenpandas",
        role: "Planned dataframe results plane for ledger tables, studies, and generated reports.",
        layers: &["planned L6 HELM reporting"],
        claims: &["planned dataframe materialization and report analytics"],
        correctness_risk: "medium",
        availability_risk: "medium",
        security_surface: "Planned table schema, query, display, and export inputs; no measured FrankenSim API surface today.",
        review_status: "Pinned but unused by measured FrankenSim manifests and Rust API references.",
        verification: &["constellation.lock pin and clean-tree verification only"],
        gaps: &[
            "the planned results-plane contract has no FrankenSim consumer or compatibility test",
            "no runtime claim may cite FrankenPandas until an actual routed surface and tests land",
        ],
        priority: 4,
        review_targets: &[
            "first-consumer admission",
            "table identity",
            "deterministic rendering",
        ],
    },
    Policy {
        lib: "frankenscipy",
        role: "Independent numerical oracle surface for special functions, FFTs, integration, sparse kernels, and optimizers.",
        layers: &["L1 BEDROCK oracle lanes", "L4 ASCENT oracle lanes"],
        claims: &[
            "cross-checks against independently implemented numerical algorithms",
            "Gauntlet comparison evidence for overlapping semantics",
        ],
        correctness_risk: "high",
        availability_risk: "medium",
        security_surface: "Numerical options, tolerances, shapes, iteration budgets, and oracle result parsing.",
        review_status: "Pinned and used by retained oracle casebooks; not blanket authority for all production kernels.",
        verification: &[
            "constellation.lock pin and clean-tree verification",
            "retained FFT, integration, linalg, sparse, special-function, ODE, and optimizer oracle casebooks",
        ],
        gaps: &[
            "oracle agreement covers named cases and tolerances rather than all inputs",
            "shared assumptions or copied algorithms can defeat independence and require review per casebook",
        ],
        priority: 2,
        review_targets: &[
            "oracle independence",
            "tolerance semantics",
            "failure and non-finite behavior",
        ],
    },
    Policy {
        lib: "frankensqlite",
        role: "Durable storage engine under the Design Ledger, session receipts, checkpoints, and replay metadata.",
        layers: &["L6 HELM", "ledger-backed cross-layer evidence"],
        claims: &[
            "artifact and lineage durability",
            "session and idempotency receipt persistence",
            "schema migration, crash recovery, and time-travel state",
        ],
        correctness_risk: "critical",
        availability_risk: "high",
        security_surface: "SQL, schema migration, database-file, blob, transaction, and recovery inputs.",
        review_status: "Pinned and heavily exercised through fs-ledger; independent durability review remains pending.",
        verification: &[
            "constellation.lock pin and clean-tree verification",
            "fs-ledger conformance, migration, time-travel, identity-guard, and artifact batteries",
            "fs-vskeleton ledger integration and reopen tests",
        ],
        gaps: &[
            "no independent cross-repo WAL, power-loss, and checkpoint adjudication backs FrankenSim durability claims",
            "multi-GiB artifact and long-running concurrent-reader stress remains outside routine focused tests",
        ],
        priority: 1,
        review_targets: &[
            "WAL durability",
            "transaction boundaries",
            "migration refusal",
            "blob and checkpoint behavior",
        ],
    },
    Policy {
        lib: "frankentorch",
        role: "Optional reverse-mode tape and learned-component substrate for surrogate and neural representation work.",
        layers: &[
            "L1 AD bridge",
            "planned L2 neural charts",
            "planned L4 surrogates",
        ],
        claims: &[
            "optional fs-ad tape bridge",
            "planned learned-component differentiation and training",
        ],
        correctness_risk: "high",
        availability_risk: "medium",
        security_surface: "Tape graphs, tensor shapes, execution modes, checkpoints, and learned model artifacts.",
        review_status: "Pinned with an optional source-level bridge; default production authority is deliberately limited.",
        verification: &[
            "constellation.lock pin and clean-tree verification",
            "feature-gated fs-ad bridge compilation and local bridge tests",
        ],
        gaps: &[
            "no default-path Gauntlet evidence promotes learned gradients or neural charts to certified authority",
            "device, checkpoint, mixed-precision, and custom-adjoint compatibility are not cross-repo gated",
        ],
        priority: 3,
        review_targets: &[
            "tape identity",
            "gradient semantics",
            "checkpoint compatibility",
            "custom adjoint boundary",
        ],
    },
];

#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct Usage {
    runtime_consumers: BTreeSet<String>,
    dev_consumers: BTreeSet<String>,
    production_references: usize,
    test_references: usize,
    samples: Vec<Sample>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct Sample {
    path: String,
    line: usize,
    symbol: String,
}

#[derive(Clone)]
struct Sibling {
    policy: Policy,
    lock: LockRow,
    usage: Usage,
}

struct Assessment {
    lock_hash: String,
    siblings: Vec<Sibling>,
}

fn repo_for_dependency(name: &str) -> Option<&'static str> {
    if name == "asupersync" {
        Some("asupersync")
    } else if name == "fsqlite" || name.starts_with("fsqlite-") {
        Some("frankensqlite")
    } else if name == "fnx" || name.starts_with("fnx-") {
        Some("franken_networkx")
    } else if name == "fnp" || name.starts_with("fnp-") {
        Some("franken_numpy")
    } else if name == "ft" || name.starts_with("ft-") {
        Some("frankentorch")
    } else if name == "fsci" || name.starts_with("fsci-") {
        Some("frankenscipy")
    } else if name == "fp" || name.starts_with("fp-") {
        Some("frankenpandas")
    } else {
        None
    }
}

fn repo_for_source_identifier(name: &str) -> Option<&'static str> {
    if name == "asupersync" {
        Some("asupersync")
    } else if name == "fsqlite" || name.starts_with("fsqlite_") {
        Some("frankensqlite")
    } else if name == "fnx" || name.starts_with("fnx_") {
        Some("franken_networkx")
    } else if name == "fnp" || name.starts_with("fnp_") {
        Some("franken_numpy")
    } else if name == "ft" || name.starts_with("ft_") {
        Some("frankentorch")
    } else if name == "fsci" || name.starts_with("fsci_") {
        Some("frankenscipy")
    } else if name == "fp" || name.starts_with("fp_") {
        Some("frankenpandas")
    } else {
        None
    }
}

fn initial_usage() -> BTreeMap<&'static str, Usage> {
    POLICIES
        .iter()
        .map(|policy| (policy.lib, Usage::default()))
        .collect()
}

fn source_is_test(path: &str) -> bool {
    path.contains("/tests/")
        || path.contains("/benches/")
        || path.contains("/examples/")
        || path.ends_with("/build.rs")
}

fn scan_source_references(
    path: &str,
    source: &str,
    usage: &mut BTreeMap<&'static str, Usage>,
) -> Result<(), String> {
    let tokens = super::casual_rust_tokens(source)
        .map_err(|error| format!("cannot tokenize {path} for constellation usage: {error}"))?;
    for (index, token) in tokens.iter().enumerate() {
        let Some(lib) = repo_for_source_identifier(token.text) else {
            continue;
        };
        if tokens.get(index + 1).is_none_or(|token| token.text != ":")
            || tokens.get(index + 2).is_none_or(|token| token.text != ":")
        {
            continue;
        }
        let entry = usage
            .get_mut(lib)
            .expect("every source prefix has a sibling policy");
        if source_is_test(path) {
            entry.test_references = entry
                .test_references
                .checked_add(1)
                .ok_or_else(|| format!("test reference count overflow for {lib}"))?;
        } else {
            entry.production_references = entry
                .production_references
                .checked_add(1)
                .ok_or_else(|| format!("production reference count overflow for {lib}"))?;
        }
        if entry.samples.len() < MAX_SAMPLES_PER_SIBLING {
            let line = source[..token.start]
                .bytes()
                .filter(|byte| *byte == b'\n')
                .count()
                + 1;
            entry.samples.push(Sample {
                path: path.to_string(),
                line,
                symbol: token.text.to_string(),
            });
        }
    }
    Ok(())
}

fn derive_usage(root: &Path) -> Result<BTreeMap<&'static str, Usage>, String> {
    let mut usage = initial_usage();
    for manifest in super::load_workspace(root)? {
        for dependency in &manifest.runtime_deps {
            if let Some(lib) = repo_for_dependency(dependency) {
                usage
                    .get_mut(lib)
                    .expect("every dependency prefix has a sibling policy")
                    .runtime_consumers
                    .insert(manifest.name.clone());
            }
        }
        for dependency in &manifest.dev_deps {
            if let Some(lib) = repo_for_dependency(dependency) {
                usage
                    .get_mut(lib)
                    .expect("every dependency prefix has a sibling policy")
                    .dev_consumers
                    .insert(manifest.name.clone());
            }
        }
    }
    for (path, source) in super::workspace_rust_sources(root)? {
        scan_source_references(&path, &source, &mut usage)?;
    }
    Ok(usage)
}

fn build_assessment(root: &Path) -> Result<Assessment, String> {
    let lock_text = super::read_constellation_lock(&root.join("constellation.lock"))?;
    let (lock_hash, rows) = super::parse_lock_rows(&lock_text)?;
    let mut policies: BTreeMap<&str, Policy> = POLICIES
        .iter()
        .copied()
        .map(|policy| (policy.lib, policy))
        .collect();
    let mut usage = derive_usage(root)?;
    let mut siblings = Vec::with_capacity(rows.len());
    for row in rows {
        let policy = policies
            .remove(row.lib.as_str())
            .ok_or_else(|| format!("no trust policy for locked sibling {}", row.lib))?;
        let measured = usage
            .remove(policy.lib)
            .ok_or_else(|| format!("no usage bucket for locked sibling {}", row.lib))?;
        siblings.push(Sibling {
            policy,
            lock: row,
            usage: measured,
        });
    }
    if !policies.is_empty() || !usage.is_empty() {
        return Err(
            "trust policy or usage schema names a sibling absent from the lock".to_string(),
        );
    }
    let assessment = Assessment {
        lock_hash,
        siblings,
    };
    validate(&assessment)?;
    Ok(assessment)
}

fn usage_state(usage: &Usage) -> &'static str {
    if usage.production_references + usage.test_references > 0 {
        "active"
    } else if !usage.runtime_consumers.is_empty() || !usage.dev_consumers.is_empty() {
        "declared-only"
    } else {
        "pinned-unused"
    }
}

fn validate(assessment: &Assessment) -> Result<(), String> {
    if assessment.siblings.len() != POLICIES.len() {
        return Err(format!(
            "assessment has {} siblings; expected {}",
            assessment.siblings.len(),
            POLICIES.len()
        ));
    }
    let mut names = BTreeSet::new();
    for sibling in &assessment.siblings {
        let policy = sibling.policy;
        if !names.insert(policy.lib) {
            return Err(format!("duplicate trust assessment for {}", policy.lib));
        }
        if policy.role.is_empty()
            || policy.layers.is_empty()
            || policy.claims.is_empty()
            || policy.security_surface.is_empty()
            || policy.review_status.is_empty()
            || policy.verification.is_empty()
            || policy.gaps.is_empty()
            || policy.review_targets.is_empty()
            || !(1..=4).contains(&policy.priority)
        {
            return Err(format!(
                "incomplete trust assessment axes for {}",
                policy.lib
            ));
        }
        for risk in [policy.correctness_risk, policy.availability_risk] {
            if !matches!(risk, "critical" | "high" | "medium" | "low") {
                return Err(format!(
                    "unsupported risk class {risk:?} for {}",
                    policy.lib
                ));
            }
        }
        let references = sibling.usage.production_references + sibling.usage.test_references;
        if (references == 0) != sibling.usage.samples.is_empty() {
            return Err(format!(
                "{} must retain a sample iff measured source references exist",
                policy.lib
            ));
        }
        if sibling.usage.samples.len() > MAX_SAMPLES_PER_SIBLING {
            return Err(format!("{} exceeds the sample retention bound", policy.lib));
        }
    }
    Ok(())
}

fn json_string(output: &mut String, value: &str) {
    output.push('"');
    output.push_str(&super::json_escape(value));
    output.push('"');
}

fn json_array<'a>(output: &mut String, values: impl IntoIterator<Item = &'a str>) {
    output.push('[');
    for (index, value) in values.into_iter().enumerate() {
        if index > 0 {
            output.push_str(", ");
        }
        json_string(output, value);
    }
    output.push(']');
}

fn render_json(assessment: &Assessment) -> String {
    let mut output = String::new();
    let _ = writeln!(output, "{{");
    let _ = write!(output, "  \"schema\": ");
    json_string(&mut output, SCHEMA);
    let _ = write!(output, ",\n  \"identity_domain\": ");
    json_string(&mut output, IDENTITY_DOMAIN);
    let _ = writeln!(output, ",\n  \"identity_version\": 1,");
    let _ = write!(output, "  \"bead_id\": ");
    json_string(&mut output, BEAD_ID);
    let _ = write!(output, ",\n  \"constellation_lock_hash\": ");
    json_string(&mut output, &assessment.lock_hash);
    let _ = write!(output, ",\n  \"usage_authority\": ");
    json_string(
        &mut output,
        "xtask derives dependency consumers from Cargo manifests and API references from tokenized workspace Rust sources; governance fields are curated",
    );
    let _ = writeln!(output, ",\n  \"siblings\": [");
    for (sibling_index, sibling) in assessment.siblings.iter().enumerate() {
        let policy = sibling.policy;
        let usage = &sibling.usage;
        let _ = writeln!(output, "    {{");
        let _ = write!(output, "      \"lib\": ");
        json_string(&mut output, policy.lib);
        let _ = write!(output, ",\n      \"version\": ");
        json_string(&mut output, &sibling.lock.version);
        let _ = write!(output, ",\n      \"git_head\": ");
        json_string(&mut output, &sibling.lock.git_head);
        let _ = write!(output, ",\n      \"role\": ");
        json_string(&mut output, policy.role);
        let _ = write!(output, ",\n      \"layers\": ");
        json_array(&mut output, policy.layers.iter().copied());
        let _ = write!(output, ",\n      \"load_bearing_claims\": ");
        json_array(&mut output, policy.claims.iter().copied());
        let _ = write!(output, ",\n      \"risk\": {{\"correctness\": ");
        json_string(&mut output, policy.correctness_risk);
        let _ = write!(output, ", \"availability\": ");
        json_string(&mut output, policy.availability_risk);
        let _ = write!(output, ", \"security_surface\": ");
        json_string(&mut output, policy.security_surface);
        let _ = write!(output, "}},\n      \"review_status\": ");
        json_string(&mut output, policy.review_status);
        let _ = write!(output, ",\n      \"verification\": ");
        json_array(&mut output, policy.verification.iter().copied());
        let _ = write!(output, ",\n      \"gaps\": ");
        json_array(&mut output, policy.gaps.iter().copied());
        let _ = write!(
            output,
            ",\n      \"review_priority\": {},\n      \"review_targets\": ",
            policy.priority
        );
        json_array(&mut output, policy.review_targets.iter().copied());
        let _ = write!(output, ",\n      \"usage\": {{\"state\": ");
        json_string(&mut output, usage_state(usage));
        let _ = write!(output, ", \"runtime_consumers\": ");
        json_array(
            &mut output,
            usage.runtime_consumers.iter().map(String::as_str),
        );
        let _ = write!(output, ", \"dev_consumers\": ");
        json_array(&mut output, usage.dev_consumers.iter().map(String::as_str));
        let _ = write!(
            output,
            ", \"production_api_references\": {}, \"test_api_references\": {}, \"samples\": [",
            usage.production_references, usage.test_references
        );
        for (sample_index, sample) in usage.samples.iter().enumerate() {
            if sample_index > 0 {
                output.push_str(", ");
            }
            output.push_str("{\"path\": ");
            json_string(&mut output, &sample.path);
            let _ = write!(output, ", \"line\": {}, \"symbol\": ", sample.line);
            json_string(&mut output, &sample.symbol);
            output.push('}');
        }
        let _ = write!(output, "]}}\n    }}");
        if sibling_index + 1 != assessment.siblings.len() {
            output.push(',');
        }
        output.push('\n');
    }
    output.push_str("  ]\n}\n");
    output
}

fn markdown_list(values: &[&str]) -> String {
    values.join("; ")
}

fn consumer_list(values: &BTreeSet<String>) -> String {
    if values.is_empty() {
        "none".to_string()
    } else {
        values.iter().cloned().collect::<Vec<_>>().join(", ")
    }
}

fn render_markdown(assessment: &Assessment) -> String {
    let mut output = String::new();
    let _ = writeln!(output, "# FrankenSim Constellation Trust-Cone Assessment\n");
    let _ = writeln!(
        output,
        "> Generated by `cargo run -p xtask -- generate-constellation-assessment`; do not edit by hand. The pinned input is `constellation.lock` hash `{}`.\n",
        assessment.lock_hash
    );
    output.push_str("This document separates two authorities: manifest consumers and Rust API references are measured from the live workspace, while roles, risk classes, verification gaps, and review priorities are explicit governance judgments. A pin proves content identity and cleanliness; it does not prove sibling correctness.\n\n");
    output.push_str("## Measured usage summary\n\n");
    output.push_str("| Sibling | State | Runtime consumers | Dev consumers | Production refs | Test refs | Correctness | Availability | Review priority |\n");
    output.push_str("| --- | --- | ---: | ---: | ---: | ---: | --- | --- | ---: |\n");
    for sibling in &assessment.siblings {
        let usage = &sibling.usage;
        let _ = writeln!(
            output,
            "| `{}` | {} | {} | {} | {} | {} | {} | {} | P{} |",
            sibling.policy.lib,
            usage_state(usage),
            usage.runtime_consumers.len(),
            usage.dev_consumers.len(),
            usage.production_references,
            usage.test_references,
            sibling.policy.correctness_risk,
            sibling.policy.availability_risk,
            sibling.policy.priority
        );
    }
    output.push_str("\n`pinned-unused` means the sibling is governed by the lock but no manifest dependency or tokenized `crate_name::...` reference was measured. It is an explicit no-claim state, not evidence that the planned role exists.\n\n");
    output.push_str("## Prioritized findings\n\n");
    output.push_str("1. **Independent review starts with asupersync cancellation and FrankenSQLite durability (`f85xj.13.5`).** Both sit directly beneath cross-layer correctness claims and remain only partially exercised from FrankenSim.\n");
    output.push_str("2. **The compatibility suite follows measured surfaces (`f85xj.13.4`).** Active dependency and API-reference rows define the first release-train matrix; planned-only surfaces cannot silently enter it as implemented claims.\n");
    output.push_str("3. **The SBOM/source manifest binds all seven pins (`f85xj.13.2`).** It must retain pinned-but-unused siblings so absence of use is visible rather than omitted.\n");
    output.push_str("4. **Governance covers availability as well as hashes (`f85xj.13.6`).** Release cadence, incident response, archival, and maintainer continuity are unresolved for the critical siblings.\n\n");
    output.push_str("## Per-sibling assessment\n\n");
    for sibling in &assessment.siblings {
        let policy = sibling.policy;
        let usage = &sibling.usage;
        let _ = writeln!(output, "### `{}`\n", policy.lib);
        let _ = writeln!(
            output,
            "- **Pinned identity:** `{}` at `{}`",
            sibling.lock.version, sibling.lock.git_head
        );
        let _ = writeln!(output, "- **Role:** {}", policy.role);
        let _ = writeln!(output, "- **Layers:** {}", markdown_list(policy.layers));
        let _ = writeln!(
            output,
            "- **Load-bearing claims:** {}",
            markdown_list(policy.claims)
        );
        let _ = writeln!(
            output,
            "- **Measured usage:** `{}`; runtime consumers: {}; dev consumers: {}; production/test references: {}/{}",
            usage_state(usage),
            consumer_list(&usage.runtime_consumers),
            consumer_list(&usage.dev_consumers),
            usage.production_references,
            usage.test_references
        );
        if usage.samples.is_empty() {
            output.push_str("- **Sampled API references:** none measured\n");
        } else {
            output.push_str("- **Sampled API references:**\n");
            for sample in &usage.samples {
                let _ = writeln!(
                    output,
                    "  - `{}` at `{}`:{}",
                    sample.symbol, sample.path, sample.line
                );
            }
        }
        let _ = writeln!(
            output,
            "- **Risk:** correctness `{}`, availability `{}`; security surface: {}",
            policy.correctness_risk, policy.availability_risk, policy.security_surface
        );
        let _ = writeln!(output, "- **Review status:** {}", policy.review_status);
        let _ = writeln!(
            output,
            "- **Current verification:** {}",
            markdown_list(policy.verification)
        );
        let _ = writeln!(
            output,
            "- **Gaps/no-claim boundary:** {}",
            markdown_list(policy.gaps)
        );
        let _ = writeln!(
            output,
            "- **Review priority P{}:** {}\n",
            policy.priority,
            markdown_list(policy.review_targets)
        );
    }
    output.push_str("## Regeneration and checking\n\n");
    output.push_str("```bash\ncargo run -p xtask -- generate-constellation-assessment\ncargo run -p xtask -- check-constellation-assessment\n```\n\n");
    output.push_str("The check fails if the seven-axis schema is incomplete, the lock identity changes, measured consumers/references drift, or either retained artifact differs from deterministic rendering. It does not establish that a reference executes in a default feature set or that a sibling is correct.\n");
    output
}

fn expected_artifacts(root: &Path) -> Result<(String, String), String> {
    let assessment = build_assessment(root)?;
    Ok((render_json(&assessment), render_markdown(&assessment)))
}

pub(crate) fn generate(root: &Path) -> Result<(), String> {
    let (json, markdown) = expected_artifacts(root)?;
    std::fs::write(root.join(DATA_PATH), json)
        .map_err(|error| format!("cannot write {DATA_PATH}: {error}"))?;
    std::fs::write(root.join(DOC_PATH), markdown)
        .map_err(|error| format!("cannot write {DOC_PATH}: {error}"))?;
    Ok(())
}

pub(crate) fn check(root: &Path) -> Vec<Violation> {
    let (json, markdown) = match expected_artifacts(root) {
        Ok(artifacts) => artifacts,
        Err(detail) => {
            return vec![Violation {
                check: CHECK,
                crate_name: "<repo>".to_string(),
                detail,
            }];
        }
    };
    [(DATA_PATH, json), (DOC_PATH, markdown)]
        .into_iter()
        .filter_map(|(path, expected)| match std::fs::read_to_string(root.join(path)) {
            Ok(actual) if actual == expected => None,
            Ok(_) => Some(Violation {
                check: CHECK,
                crate_name: path.to_string(),
                detail: format!(
                    "tracked assessment is stale; run cargo run -p xtask -- generate-constellation-assessment"
                ),
            }),
            Err(error) => Some(Violation {
                check: CHECK,
                crate_name: path.to_string(),
                detail: format!("cannot read retained assessment artifact: {error}"),
            }),
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn usage_map_finds_sampled_calls_and_preserves_unused_state() {
        let mut usage = initial_usage();
        scan_source_references(
            "crates/fs-exec/src/lib.rs",
            "use asupersync::Cx; fn run() { let _ = asupersync::Cx::for_testing(); }",
            &mut usage,
        )
        .expect("production fixture scans");
        scan_source_references(
            "crates/fs-ledger/tests/e2e.rs",
            "let db = fsqlite::Connection::open(path);",
            &mut usage,
        )
        .expect("test fixture scans");
        assert_eq!(usage["asupersync"].production_references, 2);
        assert_eq!(usage["frankensqlite"].test_references, 1);
        assert_eq!(usage_state(&usage["frankenpandas"]), "pinned-unused");
        assert!(usage["frankenpandas"].samples.is_empty());
    }

    #[test]
    fn tracked_workspace_assessment_is_complete_and_renders_deterministically() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .expect("xtask has a workspace parent");
        let assessment = build_assessment(root).expect("live assessment builds");
        validate(&assessment).expect("all seven assessment axes are complete");
        assert_eq!(render_json(&assessment), render_json(&assessment));
        assert_eq!(render_markdown(&assessment), render_markdown(&assessment));
        let pandas = assessment
            .siblings
            .iter()
            .find(|sibling| sibling.policy.lib == "frankenpandas")
            .expect("FrankenPandas is retained even while unused");
        assert_eq!(usage_state(&pandas.usage), "pinned-unused");
    }
}
