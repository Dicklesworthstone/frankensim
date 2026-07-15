//! Proposed new-domain Cargo graph admission.
//!
//! The fixture is intentionally ahead of implementation: `present` edges must
//! already exist, while `proposed` edges reserve the only edges that may appear
//! as the expansion lands.  Cargo metadata supplies the observed graph; this
//! module refuses undeclared edges, wrong layers, same-layer cycles/order drift,
//! and duplicate ownership of the program's load-bearing public types.

use super::depgraph::{JsonParser, JsonValue};
use super::{Layer, PolicyNote, Violation};
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Component, Path, PathBuf};
use std::process::Command;

const CHECK: &str = "manifest-fixture";
const FIXTURE_PATH: &str = "proposed-manifest-fixture.json";
const FIXTURE_SCHEMA: &str = "frankensim-proposed-manifest-fixture-v1";
const AUDITED_SOURCE_COMMIT: &str = "2acc39f13b49724cf94cb6bb88afbb9c2ed3ad0b";
const MAX_FIXTURE_BYTES: u64 = 1 << 20;
const MAX_METADATA_BYTES: usize = 32 << 20;
const MAX_PACKAGES: usize = 256;
const MAX_DEPENDENCIES_PER_PACKAGE: usize = 128;
const MAX_TYPE_RULES: usize = 64;

const EXPECTED_PROPOSED_PACKAGES: &[&str] = &[
    "fs-acoustics",
    "fs-circuit",
    "fs-contact",
    "fs-control",
    "fs-em",
    "fs-gas",
    "fs-kinematics",
    "fs-machine",
    "fs-matdb",
    "fs-mbd",
    "fs-motion",
    "fs-power",
    "fs-thermal",
    "fs-thermochem",
    "fs-tribo",
];

const REQUIRED_POLICY_EDGES: &[(&str, &str)] = &[
    ("fs-couple", "fs-iface"),
    ("fs-contact", "fs-motion"),
    ("fs-contact", "fs-query"),
    ("fs-contact", "fs-ivl"),
    ("fs-contact", "fs-solver"),
    ("fs-contact", "fs-tribo"),
    ("fs-contact", "fs-matdb"),
    ("fs-contact", "fs-couple"),
    ("fs-contact", "fs-iface"),
    ("fs-mbd", "fs-kinematics"),
    ("fs-mbd", "fs-contact"),
];

const FORBIDDEN_POLICY_EDGES: &[(&str, &str)] = &[
    ("fs-couple", "fs-contact"),
    ("fs-couple", "fs-mbd"),
    ("fs-couple", "fs-machine"),
    ("fs-couple", "fs-em"),
    ("fs-couple", "fs-circuit"),
    ("fs-couple", "fs-power"),
    ("fs-couple", "fs-control"),
    ("fs-couple", "fs-thermal"),
    ("fs-couple", "fs-gas"),
    ("fs-couple", "fs-tribo"),
    ("fs-couple", "fs-acoustics"),
    ("fs-couple", "fs-solid"),
    ("fs-couple", "fs-flux"),
    ("fs-contact", "fs-solid"),
    ("fs-contact", "fs-mbd"),
    ("fs-contact", "fs-machine"),
    ("fs-mbd", "fs-solid"),
    ("fs-mbd", "fs-machine"),
    ("fs-circuit", "fs-em"),
    ("fs-gas", "fs-flux"),
    ("fs-machine", "fs-em"),
    ("fs-machine", "fs-circuit"),
    ("fs-machine", "fs-power"),
    ("fs-machine", "fs-control"),
    ("fs-machine", "fs-thermal"),
    ("fs-machine", "fs-gas"),
    ("fs-machine", "fs-acoustics"),
    ("fs-machine", "fs-solid"),
    ("fs-machine", "fs-flux"),
    ("fs-motion", "fs-scenario"),
];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Presence {
    Required,
    Proposed,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum EdgeState {
    Present,
    Proposed,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TypeState {
    Present,
    Proposed,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct EdgeKey {
    to: String,
    kind: String,
    target: Option<String>,
}

#[derive(Debug, Clone)]
struct EdgeSpec {
    key: EdgeKey,
    state: EdgeState,
}

#[derive(Debug, Clone)]
struct PackageSpec {
    name: String,
    manifest_path: String,
    layer: Layer,
    presence: Presence,
    dependencies: Vec<EdgeSpec>,
    forbidden_dependencies: Vec<String>,
}

#[derive(Debug, Clone)]
struct WorkspaceSpec {
    id: String,
    manifest_path: String,
    packages: Vec<PackageSpec>,
    same_layer_topological_order: Vec<String>,
    minimal_compile_roots: Vec<String>,
}

#[derive(Debug, Clone)]
struct UniqueTypeRule {
    name: String,
    owner: String,
    state: TypeState,
    scan_crates: Vec<String>,
}

#[derive(Debug, Clone)]
struct Fixture {
    schema: String,
    source_commit: String,
    source_plan: String,
    bead_id: String,
    workspaces: Vec<WorkspaceSpec>,
    unique_types: Vec<UniqueTypeRule>,
}

#[derive(Debug, Clone)]
struct MetadataPackage {
    name: String,
    manifest_path: String,
    layer: Layer,
    dependencies: BTreeSet<EdgeKey>,
    source_dir: PathBuf,
}

#[derive(Debug, Default)]
pub(super) struct FixtureReport {
    pub(super) violations: Vec<Violation>,
    pub(super) decisions: Vec<PolicyNote>,
}

impl FixtureReport {
    fn admit(&mut self, crate_name: impl Into<String>, detail: impl Into<String>) {
        self.decisions.push(PolicyNote {
            check: CHECK,
            crate_name: crate_name.into(),
            verdict: "admit",
            detail: detail.into(),
        });
    }

    fn reject(&mut self, crate_name: impl Into<String>, detail: impl Into<String>) {
        let crate_name = crate_name.into();
        let detail = detail.into();
        self.decisions.push(PolicyNote {
            check: CHECK,
            crate_name: crate_name.clone(),
            verdict: "reject",
            detail: detail.clone(),
        });
        self.violations.push(Violation {
            check: CHECK,
            crate_name,
            detail,
        });
    }

    fn merge(&mut self, mut other: FixtureReport) {
        self.violations.append(&mut other.violations);
        self.decisions.append(&mut other.decisions);
    }
}

fn json_object<'a>(
    value: &'a JsonValue,
    context: &str,
) -> Result<&'a BTreeMap<String, JsonValue>, String> {
    if let JsonValue::Object(object) = value {
        Ok(object)
    } else {
        Err(format!("{context} must be an object"))
    }
}

fn json_array<'a>(value: &'a JsonValue, context: &str) -> Result<&'a [JsonValue], String> {
    if let JsonValue::Array(values) = value {
        Ok(values)
    } else {
        Err(format!("{context} must be an array"))
    }
}

fn json_string<'a>(value: &'a JsonValue, context: &str) -> Result<&'a str, String> {
    if let JsonValue::String(value) = value {
        Ok(value)
    } else {
        Err(format!("{context} must be a string"))
    }
}

fn required_field<'a>(
    object: &'a BTreeMap<String, JsonValue>,
    key: &str,
    context: &str,
) -> Result<&'a JsonValue, String> {
    object
        .get(key)
        .ok_or_else(|| format!("{context} is missing field {key:?}"))
}

fn exact_fields(
    object: &BTreeMap<String, JsonValue>,
    expected: &[&str],
    context: &str,
) -> Result<(), String> {
    let expected: BTreeSet<&str> = expected.iter().copied().collect();
    for field in object.keys() {
        if !expected.contains(field.as_str()) {
            return Err(format!("{context} has unknown field {field:?}"));
        }
    }
    for field in expected {
        if !object.contains_key(field) {
            return Err(format!("{context} is missing field {field:?}"));
        }
    }
    Ok(())
}

fn parse_string_array(value: &JsonValue, context: &str, max: usize) -> Result<Vec<String>, String> {
    let values = json_array(value, context)?;
    if values.len() > max {
        return Err(format!("{context} exceeds the {max}-entry bound"));
    }
    values
        .iter()
        .enumerate()
        .map(|(index, value)| {
            json_string(value, &format!("{context}[{index}]")).map(str::to_string)
        })
        .collect()
}

fn valid_package_name(name: &str) -> bool {
    !name.is_empty()
        && !name.contains('/')
        && !name.contains('\\')
        && name
            .bytes()
            .all(|byte| byte.is_ascii_lowercase() || byte.is_ascii_digit() || byte == b'-')
}

fn valid_type_name(name: &str) -> bool {
    let mut bytes = name.bytes();
    bytes
        .next()
        .is_some_and(|byte| byte.is_ascii_alphabetic() || byte == b'_')
        && bytes.all(|byte| byte.is_ascii_alphanumeric() || byte == b'_')
}

fn validate_relative_path(path: &str, context: &str) -> Result<(), String> {
    if path.is_empty() || path.contains('\\') {
        return Err(format!(
            "{context} must be a non-empty slash-normalized relative path"
        ));
    }
    for component in Path::new(path).components() {
        if !matches!(component, Component::Normal(_)) {
            return Err(format!(
                "{context} is not a safe repository-relative path: {path:?}"
            ));
        }
    }
    Ok(())
}

fn parse_edge(value: &JsonValue, context: &str) -> Result<EdgeSpec, String> {
    let object = json_object(value, context)?;
    exact_fields(object, &["to", "kind", "target", "state"], context)?;
    let to = json_string(
        required_field(object, "to", context)?,
        &format!("{context}.to"),
    )?;
    if !valid_package_name(to) {
        return Err(format!(
            "{context}.to must be one explicit package name; slash/shared shorthand is forbidden"
        ));
    }
    let kind = json_string(
        required_field(object, "kind", context)?,
        &format!("{context}.kind"),
    )?;
    if !matches!(kind, "normal" | "build") {
        return Err(format!("{context}.kind must be normal or build"));
    }
    let target = match required_field(object, "target", context)? {
        JsonValue::Null => None,
        value => Some(json_string(value, &format!("{context}.target"))?.to_string()),
    };
    let state = match json_string(
        required_field(object, "state", context)?,
        &format!("{context}.state"),
    )? {
        "present" => EdgeState::Present,
        "proposed" => EdgeState::Proposed,
        other => return Err(format!("{context}.state has unknown value {other:?}")),
    };
    Ok(EdgeSpec {
        key: EdgeKey {
            to: to.to_string(),
            kind: kind.to_string(),
            target,
        },
        state,
    })
}

fn parse_package(value: &JsonValue, context: &str) -> Result<PackageSpec, String> {
    let object = json_object(value, context)?;
    exact_fields(
        object,
        &[
            "name",
            "manifest_path",
            "layer",
            "presence",
            "dependencies",
            "forbidden_dependencies",
        ],
        context,
    )?;
    let name = json_string(
        required_field(object, "name", context)?,
        &format!("{context}.name"),
    )?;
    if !valid_package_name(name) {
        return Err(format!(
            "{context}.name must be one explicit package name; slash/shared shorthand is forbidden"
        ));
    }
    let manifest_path = json_string(
        required_field(object, "manifest_path", context)?,
        &format!("{context}.manifest_path"),
    )?;
    validate_relative_path(manifest_path, &format!("{context}.manifest_path"))?;
    let layer_text = json_string(
        required_field(object, "layer", context)?,
        &format!("{context}.layer"),
    )?;
    let layer = Layer::parse(layer_text)
        .ok_or_else(|| format!("{context}.layer has unknown value {layer_text:?}"))?;
    let presence = match json_string(
        required_field(object, "presence", context)?,
        &format!("{context}.presence"),
    )? {
        "required" => Presence::Required,
        "proposed" => Presence::Proposed,
        other => return Err(format!("{context}.presence has unknown value {other:?}")),
    };
    let dependencies_json = json_array(
        required_field(object, "dependencies", context)?,
        &format!("{context}.dependencies"),
    )?;
    if dependencies_json.len() > MAX_DEPENDENCIES_PER_PACKAGE {
        return Err(format!(
            "{context}.dependencies exceeds the {MAX_DEPENDENCIES_PER_PACKAGE}-edge bound"
        ));
    }
    let dependencies = dependencies_json
        .iter()
        .enumerate()
        .map(|(index, edge)| parse_edge(edge, &format!("{context}.dependencies[{index}]")))
        .collect::<Result<Vec<_>, _>>()?;
    let forbidden_dependencies = parse_string_array(
        required_field(object, "forbidden_dependencies", context)?,
        &format!("{context}.forbidden_dependencies"),
        MAX_DEPENDENCIES_PER_PACKAGE,
    )?;
    for dependency in &forbidden_dependencies {
        if !valid_package_name(dependency) {
            return Err(format!(
                "{context}.forbidden_dependencies contains non-explicit package {dependency:?}"
            ));
        }
    }
    Ok(PackageSpec {
        name: name.to_string(),
        manifest_path: manifest_path.to_string(),
        layer,
        presence,
        dependencies,
        forbidden_dependencies,
    })
}

fn parse_workspace(value: &JsonValue, context: &str) -> Result<WorkspaceSpec, String> {
    let object = json_object(value, context)?;
    exact_fields(
        object,
        &[
            "id",
            "manifest_path",
            "packages",
            "same_layer_topological_order",
            "minimal_compile_roots",
        ],
        context,
    )?;
    let id = json_string(
        required_field(object, "id", context)?,
        &format!("{context}.id"),
    )?;
    let manifest_path = json_string(
        required_field(object, "manifest_path", context)?,
        &format!("{context}.manifest_path"),
    )?;
    validate_relative_path(manifest_path, &format!("{context}.manifest_path"))?;
    let packages_json = json_array(
        required_field(object, "packages", context)?,
        &format!("{context}.packages"),
    )?;
    if packages_json.len() > MAX_PACKAGES {
        return Err(format!(
            "{context}.packages exceeds the {MAX_PACKAGES}-package bound"
        ));
    }
    let packages = packages_json
        .iter()
        .enumerate()
        .map(|(index, package)| parse_package(package, &format!("{context}.packages[{index}]")))
        .collect::<Result<Vec<_>, _>>()?;
    let same_layer_topological_order = parse_string_array(
        required_field(object, "same_layer_topological_order", context)?,
        &format!("{context}.same_layer_topological_order"),
        MAX_PACKAGES,
    )?;
    let minimal_compile_roots = parse_string_array(
        required_field(object, "minimal_compile_roots", context)?,
        &format!("{context}.minimal_compile_roots"),
        MAX_PACKAGES,
    )?;
    Ok(WorkspaceSpec {
        id: id.to_string(),
        manifest_path: manifest_path.to_string(),
        packages,
        same_layer_topological_order,
        minimal_compile_roots,
    })
}

fn parse_type_rule(value: &JsonValue, context: &str) -> Result<UniqueTypeRule, String> {
    let object = json_object(value, context)?;
    exact_fields(object, &["name", "owner", "state", "scan_crates"], context)?;
    let name = json_string(
        required_field(object, "name", context)?,
        &format!("{context}.name"),
    )?;
    if !valid_type_name(name) {
        return Err(format!("{context}.name is not a Rust type identifier"));
    }
    let owner = json_string(
        required_field(object, "owner", context)?,
        &format!("{context}.owner"),
    )?;
    if !valid_package_name(owner) {
        return Err(format!("{context}.owner is not one explicit package name"));
    }
    let state = match json_string(
        required_field(object, "state", context)?,
        &format!("{context}.state"),
    )? {
        "present" => TypeState::Present,
        "proposed" => TypeState::Proposed,
        other => return Err(format!("{context}.state has unknown value {other:?}")),
    };
    let scan_crates = parse_string_array(
        required_field(object, "scan_crates", context)?,
        &format!("{context}.scan_crates"),
        MAX_PACKAGES,
    )?;
    if !scan_crates.iter().all(|name| valid_package_name(name)) {
        return Err(format!(
            "{context}.scan_crates contains shorthand or an invalid name"
        ));
    }
    Ok(UniqueTypeRule {
        name: name.to_string(),
        owner: owner.to_string(),
        state,
        scan_crates,
    })
}

fn parse_fixture(text: &str) -> Result<Fixture, String> {
    let value = JsonParser::new(text).finish()?;
    let object = json_object(&value, "fixture")?;
    exact_fields(
        object,
        &[
            "schema",
            "source_commit",
            "source_plan",
            "bead_id",
            "workspaces",
            "unique_types",
        ],
        "fixture",
    )?;
    let workspaces_json = json_array(
        required_field(object, "workspaces", "fixture")?,
        "fixture.workspaces",
    )?;
    if workspaces_json.len() > 8 {
        return Err("fixture.workspaces exceeds the 8-workspace bound".to_string());
    }
    let workspaces = workspaces_json
        .iter()
        .enumerate()
        .map(|(index, workspace)| {
            parse_workspace(workspace, &format!("fixture.workspaces[{index}]"))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let types_json = json_array(
        required_field(object, "unique_types", "fixture")?,
        "fixture.unique_types",
    )?;
    if types_json.len() > MAX_TYPE_RULES {
        return Err(format!(
            "fixture.unique_types exceeds the {MAX_TYPE_RULES}-rule bound"
        ));
    }
    let unique_types = types_json
        .iter()
        .enumerate()
        .map(|(index, rule)| parse_type_rule(rule, &format!("fixture.unique_types[{index}]")))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(Fixture {
        schema: json_string(
            required_field(object, "schema", "fixture")?,
            "fixture.schema",
        )?
        .to_string(),
        source_commit: json_string(
            required_field(object, "source_commit", "fixture")?,
            "fixture.source_commit",
        )?
        .to_string(),
        source_plan: json_string(
            required_field(object, "source_plan", "fixture")?,
            "fixture.source_plan",
        )?
        .to_string(),
        bead_id: json_string(
            required_field(object, "bead_id", "fixture")?,
            "fixture.bead_id",
        )?
        .to_string(),
        workspaces,
        unique_types,
    })
}

fn package_map(workspace: &WorkspaceSpec) -> BTreeMap<&str, &PackageSpec> {
    workspace
        .packages
        .iter()
        .map(|package| (package.name.as_str(), package))
        .collect()
}

fn declared_edges(package: &PackageSpec) -> BTreeSet<EdgeKey> {
    package
        .dependencies
        .iter()
        .map(|edge| edge.key.clone())
        .collect()
}

fn find_same_layer_cycle(workspace: &WorkspaceSpec) -> Option<Vec<String>> {
    fn visit(
        name: &str,
        packages: &BTreeMap<&str, &PackageSpec>,
        state: &mut BTreeMap<String, u8>,
        stack: &mut Vec<String>,
    ) -> Option<Vec<String>> {
        state.insert(name.to_string(), 1);
        stack.push(name.to_string());
        let package = packages[name];
        let mut dependencies: Vec<&str> = package
            .dependencies
            .iter()
            .filter_map(|edge| {
                packages
                    .get(edge.key.to.as_str())
                    .filter(|dependency| dependency.layer == package.layer)
                    .map(|_| edge.key.to.as_str())
            })
            .collect();
        dependencies.sort_unstable();
        dependencies.dedup();
        for dependency in dependencies {
            match state.get(dependency).copied().unwrap_or(0) {
                0 => {
                    if let Some(cycle) = visit(dependency, packages, state, stack) {
                        return Some(cycle);
                    }
                }
                1 => {
                    let start = stack.iter().position(|entry| entry == dependency)?;
                    let mut cycle = stack[start..].to_vec();
                    cycle.push(dependency.to_string());
                    return Some(cycle);
                }
                _ => {}
            }
        }
        stack.pop();
        state.insert(name.to_string(), 2);
        None
    }

    let packages = package_map(workspace);
    let mut state = BTreeMap::new();
    let mut stack = Vec::new();
    for name in packages.keys().copied() {
        if state.get(name).copied().unwrap_or(0) == 0
            && let Some(cycle) = visit(name, &packages, &mut state, &mut stack)
        {
            return Some(cycle);
        }
    }
    None
}

fn fixture_reachable<'a>(
    workspace: &'a WorkspaceSpec,
    roots: impl IntoIterator<Item = &'a str>,
) -> BTreeSet<&'a str> {
    let packages = package_map(workspace);
    let mut reachable = BTreeSet::new();
    let mut stack: Vec<&str> = roots.into_iter().collect();
    while let Some(name) = stack.pop() {
        if !reachable.insert(name) {
            continue;
        }
        if let Some(package) = packages.get(name) {
            for edge in &package.dependencies {
                if packages.contains_key(edge.key.to.as_str()) {
                    stack.push(edge.key.to.as_str());
                }
            }
        }
    }
    reachable
}

#[allow(clippy::too_many_lines)] // The edge, order, cycle, and compile-root decisions share one report.
fn validate_workspace_static(workspace: &WorkspaceSpec) -> FixtureReport {
    let mut report = FixtureReport::default();
    let mut names = BTreeSet::new();
    for package in &workspace.packages {
        if !names.insert(package.name.clone()) {
            report.reject(
                &package.name,
                format!("workspace {} declares package more than once", workspace.id),
            );
        }
        let expected_manifest = if workspace.id == "root" {
            format!("crates/{}/Cargo.toml", package.name)
        } else if workspace.id == "fs-wasm" && package.name == "fs-wasm" {
            "crates/fs-wasm/Cargo.toml".to_string()
        } else {
            package.manifest_path.clone()
        };
        if package.manifest_path != expected_manifest {
            report.reject(
                &package.name,
                format!(
                    "manifest path {} does not match explicit package path {expected_manifest}",
                    package.manifest_path
                ),
            );
        }
        let mut edges = BTreeSet::new();
        for edge in &package.dependencies {
            if !edges.insert(edge.key.clone()) {
                report.reject(
                    &package.name,
                    format!(
                        "duplicate declared edge {} --{}--> {} target={:?}",
                        package.name, edge.key.kind, edge.key.to, edge.key.target
                    ),
                );
            }
            if package.forbidden_dependencies.contains(&edge.key.to) {
                report.reject(
                    &package.name,
                    format!("declared edge to forbidden dependency {}", edge.key.to),
                );
            } else {
                report.admit(
                    &package.name,
                    format!(
                        "declared edge {} --{}--> {} target={:?} state={:?}",
                        package.name, edge.key.kind, edge.key.to, edge.key.target, edge.state
                    ),
                );
            }
        }
        let mut forbidden = BTreeSet::new();
        for dependency in &package.forbidden_dependencies {
            if !forbidden.insert(dependency) {
                report.reject(
                    &package.name,
                    format!("forbidden dependency {dependency} is listed more than once"),
                );
            }
        }
    }

    let order: BTreeMap<&str, usize> = workspace
        .same_layer_topological_order
        .iter()
        .enumerate()
        .map(|(index, name)| (name.as_str(), index))
        .collect();
    if order.len() != workspace.same_layer_topological_order.len() {
        report.reject(
            format!("<{}>", workspace.id),
            "same-layer topological order contains a duplicate package",
        );
    }
    let order_names: BTreeSet<&str> = order.keys().copied().collect();
    let package_names: BTreeSet<&str> = workspace
        .packages
        .iter()
        .map(|package| package.name.as_str())
        .collect();
    if order_names != package_names {
        report.reject(
            format!("<{}>", workspace.id),
            "same-layer topological order must list every fixture package exactly once",
        );
    }

    let packages = package_map(workspace);
    for package in &workspace.packages {
        for edge in &package.dependencies {
            let Some(dependency) = packages.get(edge.key.to.as_str()) else {
                continue;
            };
            if !package.layer.may_depend_on(dependency.layer) {
                report.reject(
                    &package.name,
                    format!(
                        "fixture layer violation: {} ({}) must not depend on {} ({})",
                        package.name,
                        package.layer.name(),
                        dependency.name,
                        dependency.layer.name()
                    ),
                );
            } else if package.layer == dependency.layer {
                let dependency_position = order.get(dependency.name.as_str()).copied();
                let package_position = order.get(package.name.as_str()).copied();
                if dependency_position >= package_position {
                    report.reject(
                        &package.name,
                        format!(
                            "same-layer order violation: dependency {} must precede {}",
                            dependency.name, package.name
                        ),
                    );
                } else {
                    report.admit(
                        &package.name,
                        format!(
                            "same-layer order admits {} before {} in {}",
                            dependency.name,
                            package.name,
                            package.layer.name()
                        ),
                    );
                }
            }
        }
    }
    if let Some(cycle) = find_same_layer_cycle(workspace) {
        report.reject(
            format!("<{}>", workspace.id),
            format!("same-layer dependency cycle: {}", cycle.join(" -> ")),
        );
    }

    let compile_roots: BTreeSet<&str> = workspace
        .minimal_compile_roots
        .iter()
        .map(String::as_str)
        .collect();
    if compile_roots.len() != workspace.minimal_compile_roots.len() {
        report.reject(
            format!("<{}>", workspace.id),
            "minimal compile roots contain a duplicate package",
        );
    }
    let proposed: BTreeSet<&str> = workspace
        .packages
        .iter()
        .filter(|package| package.presence == Presence::Proposed)
        .map(|package| package.name.as_str())
        .collect();
    if !compile_roots.is_subset(&proposed) {
        report.reject(
            format!("<{}>", workspace.id),
            "minimal compile roots may name only proposed packages",
        );
    } else {
        let reachable = fixture_reachable(workspace, compile_roots.iter().copied());
        if !proposed.is_subset(&reachable) {
            let missing = proposed
                .difference(&reachable)
                .copied()
                .collect::<Vec<_>>()
                .join(", ");
            report.reject(
                format!("<{}>", workspace.id),
                format!("minimal compile roots do not reach proposed packages: {missing}"),
            );
        } else {
            for package in &proposed {
                report.admit(*package, "minimal compile target reaches proposed package");
            }
        }
        for root in &compile_roots {
            let reduced = compile_roots
                .iter()
                .copied()
                .filter(|candidate| candidate != root);
            if proposed.is_subset(&fixture_reachable(workspace, reduced)) {
                report.reject(
                    *root,
                    "minimal compile root is redundant; a smaller root set covers the same proposed graph",
                );
            }
        }
    }
    report
}

#[allow(clippy::too_many_lines)] // The program-set and load-bearing ownership pins are one contract.
fn validate_fixture_static(fixture: &Fixture) -> FixtureReport {
    let mut report = FixtureReport::default();
    if fixture.schema != FIXTURE_SCHEMA {
        report.reject(
            "<fixture>",
            format!(
                "schema must be {FIXTURE_SCHEMA:?}, got {:?}",
                fixture.schema
            ),
        );
    }
    if fixture.source_commit != AUDITED_SOURCE_COMMIT {
        report.reject(
            "<fixture>",
            format!(
                "source_commit must retain audited snapshot {AUDITED_SOURCE_COMMIT}, got {}",
                fixture.source_commit
            ),
        );
    }
    if fixture.source_plan != "COMPREHENSIVE_PLAN_TO_EXTEND_FRANKENSIM_TO_NEW_DOMAINS.md" {
        report.reject(
            "<fixture>",
            "source_plan is not the ratified new-domain charter",
        );
    }
    if fixture.bead_id != "frankensim-ext-manifest-fixture-r56j" {
        report.reject("<fixture>", "bead_id does not bind the owning Bead");
    }
    let workspace_ids: BTreeSet<&str> = fixture
        .workspaces
        .iter()
        .map(|workspace| workspace.id.as_str())
        .collect();
    if workspace_ids != BTreeSet::from(["root", "fs-wasm"]) {
        report.reject(
            "<fixture>",
            "workspaces must explicitly cover exactly root and standalone fs-wasm",
        );
    }
    for workspace in &fixture.workspaces {
        let expected_manifest = match workspace.id.as_str() {
            "root" => "Cargo.toml",
            "fs-wasm" => "crates/fs-wasm/Cargo.toml",
            _ => "",
        };
        if workspace.manifest_path != expected_manifest {
            report.reject(
                format!("<{}>", workspace.id),
                format!(
                    "workspace manifest {} must be {expected_manifest}",
                    workspace.manifest_path
                ),
            );
        }
        report.merge(validate_workspace_static(workspace));
    }

    if let Some(root) = fixture
        .workspaces
        .iter()
        .find(|workspace| workspace.id == "root")
    {
        let proposed: BTreeSet<&str> = root
            .packages
            .iter()
            .filter(|package| package.presence == Presence::Proposed)
            .map(|package| package.name.as_str())
            .collect();
        let expected: BTreeSet<&str> = EXPECTED_PROPOSED_PACKAGES.iter().copied().collect();
        if proposed != expected {
            report.reject(
                "<root>",
                "root fixture must enumerate the exact 15-crate expansion without shorthand or omissions",
            );
        }
        let packages = package_map(root);
        for &(from, to) in REQUIRED_POLICY_EDGES {
            let declared = packages
                .get(from)
                .is_some_and(|package| package.dependencies.iter().any(|edge| edge.key.to == to));
            if !declared {
                report.reject(from, format!("load-bearing edge to {to} is missing"));
            }
        }
        for &(from, to) in FORBIDDEN_POLICY_EDGES {
            let Some(package) = packages.get(from) else {
                report.reject(from, "load-bearing package is missing from root fixture");
                continue;
            };
            if !package
                .forbidden_dependencies
                .iter()
                .any(|dependency| dependency == to)
            {
                report.reject(
                    from,
                    format!("load-bearing forbidden edge to {to} is not pinned"),
                );
            }
            if package.dependencies.iter().any(|edge| edge.key.to == to) {
                report.reject(
                    from,
                    format!("load-bearing forbidden edge to {to} is declared"),
                );
            }
        }
    }

    let mut type_names = BTreeSet::new();
    for rule in &fixture.unique_types {
        if !type_names.insert(rule.name.as_str()) {
            report.reject(
                &rule.owner,
                format!("duplicate unique-type rule for {}", rule.name),
            );
        }
        let scan: BTreeSet<&str> = rule.scan_crates.iter().map(String::as_str).collect();
        if scan.len() != rule.scan_crates.len() {
            report.reject(
                &rule.owner,
                format!("{} scan_crates contains a duplicate", rule.name),
            );
        }
        if !scan.contains(rule.owner.as_str()) {
            report.reject(
                &rule.owner,
                format!("{} ownership rule does not scan its owner", rule.name),
            );
        }
    }
    report
}

fn metadata_string<'a>(
    object: &'a BTreeMap<String, JsonValue>,
    key: &str,
    context: &str,
) -> Result<&'a str, String> {
    json_string(
        required_field(object, key, context)?,
        &format!("{context}.{key}"),
    )
}

fn metadata_target(value: &JsonValue, context: &str) -> Result<Option<String>, String> {
    match value {
        JsonValue::Null => Ok(None),
        value => Ok(Some(json_string(value, context)?.to_string())),
    }
}

fn parse_metadata(root: &Path, text: &str) -> Result<BTreeMap<String, MetadataPackage>, String> {
    let value = JsonParser::new(text).finish()?;
    let object = json_object(&value, "cargo metadata")?;
    let packages = json_array(
        required_field(object, "packages", "cargo metadata")?,
        "cargo metadata.packages",
    )?;
    if packages.len() > 2048 {
        return Err("Cargo metadata exceeds the 2048-package bound".to_string());
    }
    let mut result = BTreeMap::new();
    for (index, package_value) in packages.iter().enumerate() {
        let context = format!("cargo metadata.packages[{index}]");
        let package = json_object(package_value, &context)?;
        let name = metadata_string(package, "name", &context)?.to_string();
        let manifest_absolute = PathBuf::from(metadata_string(package, "manifest_path", &context)?);
        let manifest_relative = manifest_absolute
            .strip_prefix(root)
            .map_err(|_| format!("{context}.manifest_path escaped repository root"))?
            .to_string_lossy()
            .replace('\\', "/");
        let metadata = json_object(
            required_field(package, "metadata", &context)?,
            &format!("{context}.metadata"),
        )?;
        let frankensim = json_object(
            required_field(metadata, "frankensim", &format!("{context}.metadata"))?,
            &format!("{context}.metadata.frankensim"),
        )?;
        let layer_text = metadata_string(
            frankensim,
            "layer",
            &format!("{context}.metadata.frankensim"),
        )?;
        let layer = Layer::parse(layer_text)
            .ok_or_else(|| format!("{context} has unknown FrankenSim layer {layer_text:?}"))?;
        let dependencies = json_array(
            required_field(package, "dependencies", &context)?,
            &format!("{context}.dependencies"),
        )?;
        let mut edges = BTreeSet::new();
        for (dependency_index, dependency_value) in dependencies.iter().enumerate() {
            let dependency_context = format!("{context}.dependencies[{dependency_index}]");
            let dependency = json_object(dependency_value, &dependency_context)?;
            let kind = match required_field(dependency, "kind", &dependency_context)? {
                JsonValue::Null => "normal".to_string(),
                value => json_string(value, &format!("{dependency_context}.kind"))?.to_string(),
            };
            if kind == "dev" {
                continue;
            }
            if !matches!(kind.as_str(), "normal" | "build") {
                return Err(format!(
                    "{dependency_context}.kind has unsupported value {kind:?}"
                ));
            }
            let edge = EdgeKey {
                to: metadata_string(dependency, "name", &dependency_context)?.to_string(),
                kind,
                target: metadata_target(
                    required_field(dependency, "target", &dependency_context)?,
                    &format!("{dependency_context}.target"),
                )?,
            };
            if !edges.insert(edge.clone()) {
                return Err(format!(
                    "{dependency_context} duplicates edge to {} kind={} target={:?}",
                    edge.to, edge.kind, edge.target
                ));
            }
        }
        let source_dir = manifest_absolute
            .parent()
            .ok_or_else(|| format!("{context}.manifest_path has no parent"))?
            .to_path_buf();
        let metadata_package = MetadataPackage {
            name: name.clone(),
            manifest_path: manifest_relative,
            layer,
            dependencies: edges,
            source_dir,
        };
        if result.insert(name.clone(), metadata_package).is_some() {
            return Err(format!(
                "Cargo metadata contains duplicate package name {name:?}"
            ));
        }
    }
    Ok(result)
}

fn run_metadata(
    root: &Path,
    workspace: &WorkspaceSpec,
) -> Result<BTreeMap<String, MetadataPackage>, String> {
    let cargo = std::env::var_os("CARGO").unwrap_or_else(|| "cargo".into());
    let output = Command::new(cargo)
        .current_dir(root)
        .args([
            "metadata",
            "--format-version",
            "1",
            "--locked",
            "--no-deps",
            "--manifest-path",
        ])
        .arg(root.join(&workspace.manifest_path))
        .output()
        .map_err(|error| format!("cannot run Cargo metadata for {}: {error}", workspace.id))?;
    if output.stdout.len() > MAX_METADATA_BYTES || output.stderr.len() > MAX_METADATA_BYTES {
        return Err(format!(
            "Cargo metadata for {} exceeds the 32 MiB output bound",
            workspace.id
        ));
    }
    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        let detail: String = stderr.chars().take(4096).collect();
        return Err(format!(
            "cargo metadata --locked --no-deps failed for {}: {detail}",
            workspace.id
        ));
    }
    let text = String::from_utf8(output.stdout)
        .map_err(|error| format!("Cargo metadata for {} is not UTF-8: {error}", workspace.id))?;
    parse_metadata(root, &text)
}

#[allow(clippy::too_many_lines)] // Each observed/declared edge emits its own structured decision.
fn compare_workspace_metadata(
    root: &Path,
    workspace: &WorkspaceSpec,
    metadata: &BTreeMap<String, MetadataPackage>,
) -> FixtureReport {
    let mut report = FixtureReport::default();
    let specifications = package_map(workspace);
    let order: BTreeMap<&str, usize> = workspace
        .same_layer_topological_order
        .iter()
        .enumerate()
        .map(|(index, name)| (name.as_str(), index))
        .collect();
    for package in &workspace.packages {
        for edge in &package.dependencies {
            if let Some(dependency) = metadata.get(&edge.key.to)
                && !package.layer.may_depend_on(dependency.layer)
            {
                report.reject(
                    &package.name,
                    format!(
                        "declared layer violation: {} ({}) must not depend on {} ({})",
                        package.name,
                        package.layer.name(),
                        dependency.name,
                        dependency.layer.name()
                    ),
                );
            }
        }
    }
    for package in &workspace.packages {
        let Some(actual) = metadata.get(&package.name) else {
            if package.presence == Presence::Required {
                report.reject(
                    &package.name,
                    "required package is absent from Cargo metadata",
                );
            } else {
                let manifest = root.join(&package.manifest_path);
                if manifest.is_file() {
                    report.reject(
                        &package.name,
                        format!(
                            "proposed manifest {} exists but the package is absent from Cargo metadata/workspace membership",
                            package.manifest_path
                        ),
                    );
                } else {
                    report.admit(
                        &package.name,
                        "proposed package is reserved but not yet present",
                    );
                }
            }
            continue;
        };
        if actual.name != package.name {
            report.reject(&package.name, "Cargo metadata package identity mismatch");
        }
        if actual.manifest_path != package.manifest_path {
            report.reject(
                &package.name,
                format!(
                    "Cargo metadata manifest {} differs from declared {}",
                    actual.manifest_path, package.manifest_path
                ),
            );
        }
        if actual.layer != package.layer {
            report.reject(
                &package.name,
                format!(
                    "Cargo metadata layer {} differs from declared {}",
                    actual.layer.name(),
                    package.layer.name()
                ),
            );
        }
        let declared = declared_edges(package);
        for edge in &actual.dependencies {
            if !declared.contains(edge) {
                report.reject(
                    &package.name,
                    format!(
                        "undeclared observed edge {} --{}--> {} target={:?}",
                        package.name, edge.kind, edge.to, edge.target
                    ),
                );
            } else if package.forbidden_dependencies.contains(&edge.to) {
                report.reject(
                    &package.name,
                    format!("Cargo metadata contains forbidden edge to {}", edge.to),
                );
            } else {
                report.admit(
                    &package.name,
                    format!(
                        "observed edge matches declaration: {} --{}--> {} target={:?}",
                        package.name, edge.kind, edge.to, edge.target
                    ),
                );
            }
            if let Some(dependency) = metadata.get(&edge.to) {
                if !actual.layer.may_depend_on(dependency.layer) {
                    report.reject(
                        &package.name,
                        format!(
                            "observed layer violation: {} ({}) depends on {} ({})",
                            package.name,
                            actual.layer.name(),
                            dependency.name,
                            dependency.layer.name()
                        ),
                    );
                }
            }
            if let Some(dependency) = specifications.get(edge.to.as_str())
                && package.layer == dependency.layer
                && order.get(dependency.name.as_str()) >= order.get(package.name.as_str())
            {
                report.reject(
                    &package.name,
                    format!(
                        "observed same-layer edge violates declared order: {} must precede {}",
                        dependency.name, package.name
                    ),
                );
            }
        }
        for edge in &package.dependencies {
            if edge.state == EdgeState::Present && !actual.dependencies.contains(&edge.key) {
                report.reject(
                    &package.name,
                    format!(
                        "declared-present edge is absent: {} --{}--> {} target={:?}",
                        package.name, edge.key.kind, edge.key.to, edge.key.target
                    ),
                );
            } else if edge.state == EdgeState::Proposed && !actual.dependencies.contains(&edge.key)
            {
                report.admit(
                    &package.name,
                    format!(
                        "proposed edge remains reserved: {} --{}--> {} target={:?}",
                        package.name, edge.key.kind, edge.key.to, edge.key.target
                    ),
                );
            }
        }
    }
    report
}

#[allow(clippy::too_many_lines)] // One bounded lexer keeps comments and strings out of ownership checks.
fn rust_identifiers(source: &str) -> Result<Vec<(String, usize)>, String> {
    let bytes = source.as_bytes();
    let mut output = Vec::new();
    let mut cursor = 0usize;
    let mut line = 1usize;
    while cursor < bytes.len() {
        if bytes.get(cursor..cursor + 2) == Some(b"//") {
            while cursor < bytes.len() && bytes[cursor] != b'\n' {
                cursor += 1;
            }
            continue;
        }
        if bytes.get(cursor..cursor + 2) == Some(b"/*") {
            let mut depth = 1usize;
            cursor += 2;
            while cursor < bytes.len() && depth > 0 {
                if bytes.get(cursor..cursor + 2) == Some(b"/*") {
                    depth += 1;
                    cursor += 2;
                } else if bytes.get(cursor..cursor + 2) == Some(b"*/") {
                    depth -= 1;
                    cursor += 2;
                } else {
                    if bytes[cursor] == b'\n' {
                        line += 1;
                    }
                    cursor += 1;
                }
            }
            if depth != 0 {
                return Err("unterminated Rust block comment".to_string());
            }
            continue;
        }
        if bytes[cursor] == b'"' {
            cursor += 1;
            let mut escaped = false;
            while cursor < bytes.len() {
                let byte = bytes[cursor];
                cursor += 1;
                if byte == b'\n' {
                    line += 1;
                }
                if escaped {
                    escaped = false;
                } else if byte == b'\\' {
                    escaped = true;
                } else if byte == b'"' {
                    break;
                }
            }
            continue;
        }
        let raw_r = if bytes[cursor] == b'r' {
            Some(cursor)
        } else if matches!(bytes[cursor], b'b' | b'c') && bytes.get(cursor + 1) == Some(&b'r') {
            Some(cursor + 1)
        } else {
            None
        };
        if let Some(raw_r) = raw_r {
            let mut marker = raw_r + 1;
            while bytes.get(marker) == Some(&b'#') {
                marker += 1;
            }
            if bytes.get(marker) == Some(&b'"') {
                let hashes = marker - raw_r - 1;
                cursor = marker + 1;
                loop {
                    if cursor >= bytes.len() {
                        return Err("unterminated Rust raw string".to_string());
                    }
                    if bytes[cursor] == b'\n' {
                        line += 1;
                    }
                    if bytes[cursor] == b'"'
                        && bytes
                            .get(cursor + 1..cursor + 1 + hashes)
                            .is_some_and(|suffix| suffix.iter().all(|byte| *byte == b'#'))
                    {
                        cursor += hashes + 1;
                        break;
                    }
                    cursor += 1;
                }
                continue;
            }
        }
        let byte = bytes[cursor];
        if byte == b'\n' {
            line += 1;
            cursor += 1;
            continue;
        }
        if byte.is_ascii_alphabetic() || byte == b'_' {
            let start = cursor;
            cursor += 1;
            while bytes
                .get(cursor)
                .is_some_and(|next| next.is_ascii_alphanumeric() || *next == b'_')
            {
                cursor += 1;
            }
            output.push((source[start..cursor].to_string(), line));
        } else {
            cursor += 1;
        }
    }
    Ok(output)
}

fn declares_public_type(source: &str, type_name: &str) -> Result<Vec<usize>, String> {
    let tokens = rust_identifiers(source)?;
    let mut lines = Vec::new();
    for window in tokens.windows(3) {
        if window[0].0 == "pub"
            && matches!(window[1].0.as_str(), "struct" | "enum" | "trait" | "type")
            && window[2].0 == type_name
        {
            lines.push(window[0].1);
        }
    }
    Ok(lines)
}

fn rust_sources(directory: &Path) -> Result<Vec<PathBuf>, String> {
    let mut stack = vec![directory.to_path_buf()];
    let mut files = Vec::new();
    while let Some(path) = stack.pop() {
        let entries = std::fs::read_dir(&path).map_err(|error| {
            format!(
                "cannot read Rust source directory {}: {error}",
                path.display()
            )
        })?;
        for entry in entries {
            let entry =
                entry.map_err(|error| format!("cannot enumerate {}: {error}", path.display()))?;
            let entry_path = entry.path();
            let file_type = entry
                .file_type()
                .map_err(|error| format!("cannot inspect {}: {error}", entry_path.display()))?;
            if file_type.is_dir() {
                stack.push(entry_path);
            } else if file_type.is_file()
                && entry_path
                    .extension()
                    .is_some_and(|extension| extension == "rs")
            {
                files.push(entry_path);
            }
        }
    }
    files.sort();
    Ok(files)
}

fn check_unique_types(
    rules: &[UniqueTypeRule],
    metadata: &BTreeMap<String, MetadataPackage>,
) -> FixtureReport {
    let mut report = FixtureReport::default();
    for rule in rules {
        let mut occurrences = Vec::new();
        for crate_name in &rule.scan_crates {
            let Some(package) = metadata.get(crate_name) else {
                continue;
            };
            let source_root = package.source_dir.join("src");
            if !source_root.is_dir() {
                continue;
            }
            match rust_sources(&source_root) {
                Err(detail) => report.reject(crate_name, detail),
                Ok(files) => {
                    for file in files {
                        match std::fs::read_to_string(&file) {
                            Err(error) => report.reject(
                                crate_name,
                                format!("cannot read {}: {error}", file.display()),
                            ),
                            Ok(source) => match declares_public_type(&source, &rule.name) {
                                Err(detail) => report.reject(crate_name, detail),
                                Ok(lines) => {
                                    for line in lines {
                                        occurrences.push((crate_name.clone(), file.clone(), line));
                                    }
                                }
                            },
                        }
                    }
                }
            }
        }
        if occurrences.is_empty() {
            if rule.state == TypeState::Present {
                report.reject(
                    &rule.owner,
                    format!(
                        "required unique type {} is absent from its scan scope",
                        rule.name
                    ),
                );
            } else {
                report.admit(
                    &rule.owner,
                    format!("proposed unique type {} is not implemented yet", rule.name),
                );
            }
            continue;
        }
        let owners: BTreeSet<&str> = occurrences
            .iter()
            .map(|(owner, _, _)| owner.as_str())
            .collect();
        if owners != BTreeSet::from([rule.owner.as_str()]) || occurrences.len() != 1 {
            let locations = occurrences
                .iter()
                .map(|(owner, path, line)| format!("{owner}:{}:{line}", path.display()))
                .collect::<Vec<_>>()
                .join(", ");
            report.reject(
                &rule.owner,
                format!(
                    "unique type {} must have exactly one public declaration in owner {}; observed {locations}",
                    rule.name, rule.owner
                ),
            );
        } else {
            let (_, path, line) = &occurrences[0];
            report.admit(
                &rule.owner,
                format!(
                    "unique type {} is owned exactly once at {}:{line}",
                    rule.name,
                    path.display()
                ),
            );
        }
    }
    report
}

pub(super) fn check_manifest_fixture(root: &Path) -> FixtureReport {
    let path = root.join(FIXTURE_PATH);
    match std::fs::metadata(&path) {
        Ok(metadata) if metadata.is_file() && metadata.len() <= MAX_FIXTURE_BYTES => {}
        Ok(metadata) if metadata.len() > MAX_FIXTURE_BYTES => {
            let mut report = FixtureReport::default();
            report.reject(
                "<fixture>",
                format!("{FIXTURE_PATH} exceeds the 1 MiB processing bound"),
            );
            return report;
        }
        Ok(_) => {
            let mut report = FixtureReport::default();
            report.reject("<fixture>", format!("{FIXTURE_PATH} is not a regular file"));
            return report;
        }
        Err(error) => {
            let mut report = FixtureReport::default();
            report.reject(
                "<fixture>",
                format!("cannot inspect {FIXTURE_PATH}: {error}"),
            );
            return report;
        }
    }
    let text = match std::fs::read_to_string(&path) {
        Ok(text) => text,
        Err(error) => {
            let mut report = FixtureReport::default();
            report.reject("<fixture>", format!("cannot read {FIXTURE_PATH}: {error}"));
            return report;
        }
    };
    let fixture = match parse_fixture(&text) {
        Ok(fixture) => fixture,
        Err(detail) => {
            let mut report = FixtureReport::default();
            report.reject("<fixture>", detail);
            return report;
        }
    };
    let mut report = validate_fixture_static(&fixture);
    let mut root_metadata = None;
    for workspace in &fixture.workspaces {
        match run_metadata(root, workspace) {
            Ok(metadata) => {
                report.admit(
                    format!("<{}>", workspace.id),
                    format!("Cargo metadata resolved {} packages", metadata.len()),
                );
                report.merge(compare_workspace_metadata(root, workspace, &metadata));
                if workspace.id == "root" {
                    root_metadata = Some(metadata);
                }
            }
            Err(detail) => report.reject(format!("<{}>", workspace.id), detail),
        }
    }
    if let Some(metadata) = root_metadata {
        report.merge(check_unique_types(&fixture.unique_types, &metadata));
    }
    report
}

#[cfg(test)]
mod tests {
    use super::*;

    fn committed_fixture() -> Fixture {
        parse_fixture(include_str!("../../proposed-manifest-fixture.json"))
            .expect("committed fixture parses")
    }

    #[test]
    fn committed_fixture_is_a_complete_minimal_compile_target() {
        let report = validate_fixture_static(&committed_fixture());
        assert!(
            report.violations.is_empty(),
            "committed fixture static contract failed: {:?}",
            report.violations
        );
        let compile_decisions = report
            .decisions
            .iter()
            .filter(|decision| decision.detail == "minimal compile target reaches proposed package")
            .count();
        assert_eq!(compile_decisions, EXPECTED_PROPOSED_PACKAGES.len());
    }

    #[test]
    fn seeded_same_layer_cycle_is_rejected_with_structured_edge_context() {
        let mut fixture = committed_fixture();
        let root = fixture
            .workspaces
            .iter_mut()
            .find(|workspace| workspace.id == "root")
            .expect("root workspace");
        let contact = root
            .packages
            .iter_mut()
            .find(|package| package.name == "fs-contact")
            .expect("contact package");
        contact.dependencies.push(EdgeSpec {
            key: EdgeKey {
                to: "fs-machine".to_string(),
                kind: "normal".to_string(),
                target: None,
            },
            state: EdgeState::Proposed,
        });
        let report = validate_fixture_static(&fixture);
        assert!(
            report
                .violations
                .iter()
                .any(|violation| violation.detail.contains("same-layer dependency cycle")),
            "seeded cycle was not rejected: {:?}",
            report.violations
        );
        assert!(
            report
                .decisions
                .iter()
                .any(|decision| decision.verdict == "reject"
                    && decision.detail.contains("fs-contact")),
            "rejection must retain structured edge/package context"
        );
    }

    #[test]
    fn observed_undeclared_edge_is_rejected_and_logged() {
        let fixture = committed_fixture();
        let workspace = fixture
            .workspaces
            .iter()
            .find(|workspace| workspace.id == "root")
            .expect("root workspace");
        let package = workspace
            .packages
            .iter()
            .find(|package| package.name == "fs-couple")
            .expect("couple package");
        let mut metadata = BTreeMap::new();
        metadata.insert(
            package.name.clone(),
            MetadataPackage {
                name: package.name.clone(),
                manifest_path: package.manifest_path.clone(),
                layer: package.layer,
                dependencies: BTreeSet::from([EdgeKey {
                    to: "fs-solid".to_string(),
                    kind: "normal".to_string(),
                    target: None,
                }]),
                source_dir: PathBuf::new(),
            },
        );
        let report = compare_workspace_metadata(Path::new("/nonexistent"), workspace, &metadata);
        assert!(
            report
                .violations
                .iter()
                .any(|violation| violation.detail.contains("undeclared observed edge")),
            "undeclared edge was not rejected: {:?}",
            report.violations
        );
        assert!(report.decisions.iter().any(|decision| {
            decision.verdict == "reject"
                && decision.crate_name == "fs-couple"
                && decision.detail.contains("fs-solid")
        }));
    }

    #[test]
    fn duplicate_public_type_owners_are_detected_without_matching_comments_or_strings() {
        let source = r#"
            // pub struct PortSchema;
            const NOTE: &str = "pub struct PortSchema";
            pub struct PortSchema;
        "#;
        assert_eq!(declares_public_type(source, "PortSchema").unwrap(), vec![4]);
        let second = "pub enum PortSchema { V1 }";
        assert_eq!(declares_public_type(second, "PortSchema").unwrap(), vec![1]);
    }

    #[test]
    fn slash_notation_is_refused_in_executable_package_fields() {
        let text = include_str!("../../proposed-manifest-fixture.json").replacen(
            "\"name\": \"fs-motion\"",
            "\"name\": \"fs-motion/fs-query\"",
            1,
        );
        let error = parse_fixture(&text).expect_err("slash shorthand must fail");
        assert!(error.contains("slash/shared shorthand"), "{error}");
    }
}
