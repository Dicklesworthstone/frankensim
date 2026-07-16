//! FrankenScript Machine graph/behavior codec conformance (Gauntlet G0/G3).

use fs_blake3::identity::StrongIdentity;

use fs_ir::machine::codec::{
    MAX_MACHINE_BEHAVIOR_AST_NODES, MAX_MACHINE_GRAPH_AST_NODES, MachineBehaviorAstAdmissionError,
    MachineBehaviorCodecRule, MachineGraphAstAdmissionError, MachineGraphCodecRule,
    admit_machine_behavior_ast_v1, admit_machine_graph_ast_v1, parse_machine_behavior_program_v1,
    parse_machine_behavior_v1, parse_machine_graph_program_v1, parse_machine_graph_v1,
    write_machine_behavior_program_v1, write_machine_behavior_v1, write_machine_graph_program_v1,
    write_machine_graph_v1,
};
use fs_ir::machine::semantics::{
    MAX_MACHINE_BEHAVIOR_EVENTS, MAX_MACHINE_BEHAVIOR_NESTED_REFERENCES, MachineBehaviorRule,
};
use fs_ir::machine::{
    MAX_MACHINE_GRAPH_CLOCKS, MAX_MACHINE_GRAPH_OWNED_ELEMENTS, MachineGraphRule,
};
use fs_ir::{Node, NodeKind, VersionedProgram, json, sexpr};

const SOURCE_MODEL_DIGEST_BYTE: u8 = 0xab;

fn digest(byte: u8) -> String {
    format!("{byte:02x}").repeat(32)
}

fn valid_source(version: &str) -> String {
    let source_model = digest(SOURCE_MODEL_DIGEST_BYTE);
    let load_model = digest(2);
    let source_material = digest(3);
    let load_material = digest(4);
    let interface = digest(5);
    format!(
        r#"(machine-graph-v1
  (clocks
    (clock "clock/mechanical" (periodic "1000000" "0")))
  (subsystems
    (subsystem "subsystem/source"
      (ref "models/source" "{version}" "{source_model}")
      (bodies "body/source")
      (surface-patches)
      (contact-features)
      (state-slots))
    (subsystem "subsystem/load"
      (ref "models/load" "{version}" "{load_model}")
      (bodies "body/load")
      (surface-patches)
      (contact-features)
      (state-slots)))
  (terminals
    (terminal "terminal/source-effort" "subsystem/source"
      (semantic (pressure) static)
      (scalar) output "clock/mechanical" (frame "world/mechanical" preserving))
    (terminal "terminal/source-flow" "subsystem/source"
      (dims 3 0 -1 0 0 0)
      (scalar) output "clock/mechanical" (frame "world/mechanical" preserving))
    (terminal "terminal/load-effort" "subsystem/load"
      (semantic (pressure) static)
      (scalar) input "clock/mechanical" (frame "world/mechanical" preserving))
    (terminal "terminal/load-flow" "subsystem/load"
      (dims 3 0 -1 0 0 0)
      (scalar) input "clock/mechanical" (frame "world/mechanical" preserving)))
  (ports
    (port "port/source" "subsystem/source"
      "terminal/source-effort" "terminal/source-flow" out-of-subsystem)
    (port "port/load" "subsystem/load"
      "terminal/load-effort" "terminal/load-flow" into-subsystem))
  (relations
    (relation "relation/effort" "terminal/source-effort" "terminal/load-effort"
      (algebraic))
    (relation "relation/flow" "terminal/source-flow" "terminal/load-flow"
      (algebraic)))
  (materials
    (material (body "body/source")
      (ref "materials/source" "{version}" "{source_material}"))
    (material (body "body/load")
      (ref "materials/load" "{version}" "{load_material}")))
  (interfaces
    (interface "interface/source-load" "port/source" "port/load"
      (ref "interfaces/hydraulic" "{version}" "{interface}") aligned)))"#
    )
}

fn behavior_graph_source(model_byte: u8) -> String {
    let model = digest(model_byte);
    let material = digest(62);
    format!(
        r#"(machine-graph-v1
  (clocks
    (clock "clock/continuous" (continuous))
    (clock "clock/events" (event-driven)))
  (subsystems
    (subsystem "subsystem/plant"
      (ref "models/behavior-fixture" "1" "{model}")
      (bodies "body/plant")
      (surface-patches)
      (contact-features)
      (state-slots "state/position" "state/velocity")))
  (terminals
    (terminal "terminal/state-a-source" "subsystem/plant"
      (dims 0 0 0 0 0 0) (scalar) output "clock/continuous"
      (frame "world/mechanical" preserving))
    (terminal "terminal/state-a-sink" "subsystem/plant"
      (dims 0 0 0 0 0 0) (scalar) input "clock/continuous"
      (frame "world/mechanical" preserving))
    (terminal "terminal/state-b-source" "subsystem/plant"
      (dims 0 0 0 0 0 0) (scalar) output "clock/continuous"
      (frame "world/mechanical" preserving))
    (terminal "terminal/state-b-sink" "subsystem/plant"
      (dims 0 0 0 0 0 0) (scalar) input "clock/continuous"
      (frame "world/mechanical" preserving))
    (terminal "terminal/external-command" "subsystem/plant"
      (dims 0 0 0 0 0 0) (scalar) external-input "clock/continuous"
      (frame "world/mechanical" preserving))
    (terminal "terminal/guard-observation" "subsystem/plant"
      (dims 0 0 0 0 0 0) (scalar) output "clock/continuous"
      (frame "world/mechanical" preserving)))
  (ports)
  (relations
    (relation "relation/state-a" "terminal/state-a-source" "terminal/state-a-sink"
      (stateful "state/position"))
    (relation "relation/state-b" "terminal/state-b-source" "terminal/state-b-sink"
      (stateful "state/velocity")))
  (materials
    (material (body "body/plant")
      (ref "materials/behavior-fixture" "1" "{material}")))
  (interfaces))"#
    )
}

fn admitted_behavior_graph(model_byte: u8) -> fs_ir::machine::AdmittedMachineGraph {
    admit_machine_graph_ast_v1(
        &sexpr::parse(&behavior_graph_source(model_byte)).expect("behavior graph syntax parses"),
    )
    .expect("behavior graph admits")
}

fn identity_hex(bytes: [u8; 32]) -> String {
    bytes.iter().map(|byte| format!("{byte:02x}")).collect()
}

fn synthetic_symbol(value: &str) -> Node {
    Node::synthetic(NodeKind::Symbol(value.to_string()))
}

fn synthetic_string(value: &str) -> Node {
    Node::synthetic(NodeKind::Str(value.to_string()))
}

fn synthetic_form(head: &str, args: Vec<Node>) -> Node {
    let mut items = Vec::with_capacity(args.len() + 1);
    items.push(synthetic_symbol(head));
    items.extend(args);
    Node::synthetic(NodeKind::List(items))
}

fn synthetic_ref(namespace: &str, byte: u8) -> Node {
    synthetic_form(
        "ref",
        vec![
            synthetic_string(namespace),
            synthetic_string("1"),
            synthetic_string(&digest(byte)),
        ],
    )
}

fn synthetic_dims() -> Node {
    synthetic_form(
        "dims",
        (0..6).map(|_| Node::synthetic(NodeKind::Int(0))).collect(),
    )
}

fn generated_bounded_tolerance(index: usize) -> Node {
    synthetic_form(
        "tolerance",
        vec![
            synthetic_string(&format!("tolerance/generated-{index}")),
            synthetic_form(
                "element",
                vec![synthetic_form("body", vec![synthetic_string("body/plant")])],
            ),
            synthetic_ref(&format!("parameters/generated-{index}"), 71),
            synthetic_dims(),
            synthetic_form("scalar", Vec::new()),
            synthetic_form(
                "bounded",
                vec![
                    Node::synthetic(NodeKind::Float(0.01)),
                    Node::synthetic(NodeKind::Float(0.02)),
                    synthetic_ref("tolerances/generated", 72),
                ],
            ),
        ],
    )
}

fn generated_event(index: usize) -> Node {
    synthetic_form(
        "event",
        vec![
            synthetic_string(&format!("event/generated-{index}")),
            synthetic_string("clock/events"),
            synthetic_ref("guards/generated", 73),
            synthetic_symbol("negative-to-positive"),
            synthetic_form("unknown", vec![synthetic_ref("claims/generated", 74)]),
            synthetic_form(
                "dependencies",
                vec![
                    synthetic_form("state", vec![synthetic_string("state/position")]),
                    synthetic_form(
                        "terminal",
                        vec![synthetic_string("terminal/guard-observation")],
                    ),
                ],
            ),
            synthetic_form(
                "deterministic",
                vec![
                    synthetic_ref("resets/generated", 75),
                    synthetic_form(
                        "writes",
                        vec![
                            synthetic_string("state/position"),
                            synthetic_string("state/velocity"),
                        ],
                    ),
                ],
            ),
            synthetic_form("total-priority", vec![synthetic_string(&index.to_string())]),
        ],
    )
}

fn ast_node_count(node: &Node) -> usize {
    match &node.kind {
        NodeKind::List(items) => 1 + items.iter().map(ast_node_count).sum::<usize>(),
        NodeKind::Int(_)
        | NodeKind::Float(_)
        | NodeKind::Qty { .. }
        | NodeKind::Count { .. }
        | NodeKind::Seed(_)
        | NodeKind::Str(_)
        | NodeKind::Symbol(_)
        | NodeKind::Keyword(_) => 1,
    }
}

fn behavior_source(base_graph: &str) -> String {
    let value = digest(10);
    let distribution_condition = digest(11);
    let history = digest(12);
    let guard = digest(20);
    let no_claim = digest(21);
    let reset = digest(22);
    let parameter_position = digest(30);
    let tolerance_position = digest(31);
    let marginal_position = digest(32);
    let parameter_velocity = digest(33);
    let tolerance_velocity = digest(34);
    let marginal_velocity = digest(35);
    let parameter_clearance = digest(36);
    let tolerance_clearance = digest(37);
    let correlation = digest(40);
    format!(
        r#"(machine-behavior-v1
  (base-graph "{base_graph}")
  (state-contracts
    (state-contract "state/position" "subsystem/plant"
      (dims 0 0 0 0 0 0) (scalar) "clock/continuous"
      (frame "world/mechanical" preserving))
    (state-contract "state/velocity" "subsystem/plant"
      (dims 0 0 0 0 0 0) (scalar) "clock/continuous"
      (frame "world/mechanical" preserving)))
  (conditions
    (condition (initial "state/position")
      (dims 0 0 0 0 0 0) (scalar) "clock/continuous"
      (frame "world/mechanical" preserving)
      (fixed (ref "values/condition" "1" "{value}")))
    (condition (initial "state/velocity")
      (dims 0 0 0 0 0 0) (scalar) "clock/continuous"
      (frame "world/mechanical" preserving)
      (distribution (ref "distributions/marginal" "1" "{distribution_condition}")))
    (condition (boundary "terminal/external-command")
      (dims 0 0 0 0 0 0) (scalar) "clock/continuous"
      (frame "world/mechanical" preserving)
      (history (ref "signals/history" "1" "{history}")
        (reset-at-events "event/contact"))))
  (motions
    (motion "body/plant" "clock/continuous"
      (frame "world/mechanical" preserving) (static)))
  (events
    (event "event/contact" "clock/events"
      (ref "guards/state-threshold" "1" "{guard}")
      negative-to-positive
      (unknown (ref "claims/crossing-unknown" "1" "{no_claim}"))
      (dependencies
        (state "state/position")
        (terminal "terminal/guard-observation"))
      (deterministic (ref "resets/state-map" "1" "{reset}")
        (writes "state/position" "state/velocity"))
      (total-priority "0")))
  (tolerances
    (tolerance "tolerance/position"
      (element (state-slot "state/position"))
      (ref "parameters/position" "1" "{parameter_position}")
      (dims 0 0 0 0 0 0) (scalar)
      (random 0.1
        (ref "tolerances/additive" "1" "{tolerance_position}")
        (ref "distributions/marginal" "1" "{marginal_position}")))
    (tolerance "tolerance/velocity"
      (element (state-slot "state/velocity"))
      (ref "parameters/velocity" "1" "{parameter_velocity}")
      (dims 0 0 0 0 0 0) (scalar)
      (random 0.2
        (ref "tolerances/additive" "1" "{tolerance_velocity}")
        (ref "distributions/marginal" "1" "{marginal_velocity}")))
    (tolerance "tolerance/body-clearance"
      (element (body "body/plant"))
      (ref "parameters/body-clearance" "1" "{parameter_clearance}")
      (dims 0 0 0 0 0 0) (scalar)
      (bounded 0.01 0.02
        (ref "tolerances/additive" "1" "{tolerance_clearance}"))))
  (dependences
    (dependence
      (members
        (condition (initial "state/velocity"))
        (tolerance "tolerance/position")
        (tolerance "tolerance/velocity"))
      (correlated (ref "correlations/joint-law" "1" "{correlation}")))))"#
    )
}

fn parse_valid_behavior(graph: &fs_ir::machine::AdmittedMachineGraph) -> Node {
    let graph_id = identity_hex(*graph.identity().as_bytes());
    sexpr::parse(&behavior_source(&graph_id)).expect("valid Machine behavior literal")
}

fn parse_valid() -> Node {
    sexpr::parse(&valid_source("1")).expect("valid Machine graph literal")
}

fn root_items(node: &mut Node) -> &mut Vec<Node> {
    let NodeKind::List(items) = &mut node.kind else {
        panic!("fixture root is a list")
    };
    items
}

fn section_items(node: &mut Node) -> &mut Vec<Node> {
    let NodeKind::List(items) = &mut node.kind else {
        panic!("fixture section is a list")
    };
    items
}

#[test]
fn g0_literal_sexpr_and_json_publish_the_same_machine_identity() {
    let literal = parse_valid();
    let admitted = admit_machine_graph_ast_v1(&literal).expect("literal graph admits");

    let program = write_machine_graph_program_v1(&admitted).expect("admitted graph encodes");
    let sexpr_bytes = program
        .print_sexpr_checked()
        .expect("canonical s-expression");
    let json_bytes = program.print_json_checked().expect("canonical JSON");
    let from_sexpr = VersionedProgram::parse_sexpr(&sexpr_bytes).expect("s-expression reparses");
    let from_json = VersionedProgram::parse_json(&json_bytes).expect("JSON reparses");
    let admitted_sexpr = parse_machine_graph_program_v1(&from_sexpr)
        .expect("s-expression graph decodes")
        .admit()
        .expect("s-expression graph admits");
    let admitted_json = parse_machine_graph_program_v1(&from_json)
        .expect("JSON graph decodes")
        .admit()
        .expect("JSON graph admits");

    assert_eq!(admitted.identity(), admitted_sexpr.identity());
    assert_eq!(admitted.identity(), admitted_json.identity());
    let canonical_node = write_machine_graph_v1(&admitted).expect("graph writes");
    assert!(canonical_node.same_shape(from_sexpr.program()));
    assert!(canonical_node.same_shape(from_json.program()));
    assert_eq!(
        sexpr::print(&canonical_node).expect("canonical node prints"),
        sexpr::print(from_json.program()).expect("JSON-derived node prints")
    );
}

#[test]
fn g3_source_row_permutations_do_not_move_admitted_identity() {
    let baseline = admit_machine_graph_ast_v1(&parse_valid())
        .expect("baseline admits")
        .identity();
    let mut permuted = parse_valid();
    for section in &mut root_items(&mut permuted)[1..] {
        section_items(section)[1..].reverse();
    }
    let moved = admit_machine_graph_ast_v1(&permuted)
        .expect("permuted graph admits")
        .identity();
    assert_eq!(baseline, moved);
}

#[test]
fn g0_semantic_mutation_moves_identity_and_full_u64_versions_round_trip() {
    let baseline_source = valid_source("1");
    let baseline = admit_machine_graph_ast_v1(
        &sexpr::parse(&baseline_source).expect("baseline syntax parses"),
    )
    .expect("baseline admits");

    let changed_source = baseline_source.replacen(&digest(SOURCE_MODEL_DIGEST_BYTE), &digest(9), 1);
    let changed =
        admit_machine_graph_ast_v1(&sexpr::parse(&changed_source).expect("mutated syntax parses"))
            .expect("mutated graph admits");
    assert_ne!(baseline.identity(), changed.identity());

    let max_source = valid_source("18446744073709551615");
    let max = admit_machine_graph_ast_v1(&sexpr::parse(&max_source).expect("u64 syntax parses"))
        .expect("u64::MAX reference versions admit");
    let max_program = write_machine_graph_program_v1(&max).expect("max versions encode");
    assert!(max_program.print_sexpr().contains("18446744073709551615"));
    let reparsed =
        VersionedProgram::parse_json(&max_program.print_json()).expect("max-version JSON reparses");
    let readmitted = admit_machine_graph_ast_v1(reparsed.program()).expect("reparsed graph admits");
    assert_eq!(max.identity(), readmitted.identity());
}

#[test]
fn g0_codec_refusals_retain_rule_span_path_and_hint() {
    let source = valid_source("1");
    let uppercase = source.replacen(
        &digest(SOURCE_MODEL_DIGEST_BYTE),
        &digest(SOURCE_MODEL_DIGEST_BYTE).to_uppercase(),
        1,
    );
    let uppercase_node = sexpr::parse(&uppercase).expect("generic syntax remains valid");
    let error = parse_machine_graph_v1(&uppercase_node).expect_err("uppercase digest must refuse");
    assert_eq!(error.rule(), MachineGraphCodecRule::InvalidReference);
    assert_eq!(error.code(), "MachineGraphCodecInvalidReference");
    assert!(error.span().end > error.span().start);
    assert_eq!(error.path(), "$[2][1][2][3]");
    assert!(!error.detail().is_empty());
    assert!(!error.hint().is_empty());

    let leading_zero = source.replacen("\"1\"", "\"01\"", 1);
    let error =
        parse_machine_graph_v1(&sexpr::parse(&leading_zero).expect("generic syntax remains valid"))
            .expect_err("noncanonical u64 must refuse");
    assert_eq!(error.rule(), MachineGraphCodecRule::InvalidNumber);

    let zero_digest = source.replacen(&digest(SOURCE_MODEL_DIGEST_BYTE), &"0".repeat(64), 1);
    let error =
        parse_machine_graph_v1(&sexpr::parse(&zero_digest).expect("generic syntax remains valid"))
            .expect_err("zero semantic digest must refuse");
    assert_eq!(error.rule(), MachineGraphCodecRule::InvalidReference);
    assert_eq!(error.path(), "$[2][1][2][3]");

    let invalid_namespace = source.replacen("models/source", "Models Source", 1);
    let error = parse_machine_graph_v1(
        &sexpr::parse(&invalid_namespace).expect("generic syntax remains valid"),
    )
    .expect_err("invalid reference namespace must refuse");
    assert_eq!(error.rule(), MachineGraphCodecRule::InvalidReference);
    assert_eq!(error.path(), "$[2][1][2][1]");

    let wrong_order = source.replacen("(clocks", "(terminals", 1);
    let error =
        parse_machine_graph_v1(&sexpr::parse(&wrong_order).expect("generic syntax remains valid"))
            .expect_err("out-of-order section must refuse");
    assert_eq!(error.rule(), MachineGraphCodecRule::UnexpectedForm);
    assert_eq!(error.path(), "$[1]");

    let unknown_orientation = source.replacen("preserving", "sideways", 1);
    let error = parse_machine_graph_v1(
        &sexpr::parse(&unknown_orientation).expect("generic syntax remains valid"),
    )
    .expect_err("unknown nested orientation must refuse");
    assert_eq!(error.rule(), MachineGraphCodecRule::UnknownTag);
    assert_eq!(error.path(), "$[3][1][7][2]");
}

#[test]
fn g0_codec_resource_preflight_refuses_before_entry_decode() {
    assert!(MAX_MACHINE_GRAPH_AST_NODES > MAX_MACHINE_GRAPH_CLOCKS);
    let mut oversized = sexpr::parse(
        "(machine-graph-v1 (clocks) (subsystems) (terminals) (ports) (relations) (materials) (interfaces))",
    )
    .expect("empty graph syntax parses");
    let root = root_items(&mut oversized);
    let clocks = section_items(&mut root[1]);
    clocks.push(Node::synthetic(NodeKind::Float(f64::NAN)));
    clocks
        .extend((0..MAX_MACHINE_GRAPH_CLOCKS).map(|_| Node::synthetic(NodeKind::List(Vec::new()))));
    let error = parse_machine_graph_v1(&oversized).expect_err("clock cap must refuse");
    assert_eq!(error.rule(), MachineGraphCodecRule::ResourceLimit);
    assert_eq!(error.path(), "$[1]");

    let mut oversized_owned = parse_valid();
    let root = root_items(&mut oversized_owned);
    let subsystems = section_items(&mut root[2]);
    let subsystem = root_items(&mut subsystems[1]);
    let bodies = section_items(&mut subsystem[3]);
    bodies.push(Node::synthetic(NodeKind::Float(f64::NAN)));
    bodies.extend(
        (1..MAX_MACHINE_GRAPH_OWNED_ELEMENTS)
            .map(|_| Node::synthetic(NodeKind::Str("body/repeated".to_string()))),
    );
    let error = parse_machine_graph_v1(&oversized_owned)
        .expect_err("aggregate ownership cap must precede recursive AST validation");
    assert_eq!(error.rule(), MachineGraphCodecRule::ResourceLimit);
    assert_eq!(error.path(), "$[2][1][3]");
}

#[test]
fn g0_syntax_success_cannot_bypass_semantic_graph_refusal() {
    let mut unclosed = parse_valid();
    let root = root_items(&mut unclosed);
    let relations = section_items(&mut root[5]);
    relations.remove(1);

    let error = admit_machine_graph_ast_v1(&unclosed).expect_err("missing source must refuse");
    let MachineGraphAstAdmissionError::Graph(refusal) = error else {
        panic!("valid syntax must reach semantic graph refusal")
    };
    assert!(
        refusal
            .findings()
            .iter()
            .any(|finding| finding.rule() == MachineGraphRule::MissingSourceClosure)
    );
}

#[test]
fn g3_generic_json_shape_is_not_a_second_machine_grammar() {
    let node = parse_valid();
    let json_bytes = json::print(&node).expect("generic JSON prints");
    let from_json = json::parse(&json_bytes).expect("generic JSON reparses");
    let left = admit_machine_graph_ast_v1(&node).expect("s-expression AST admits");
    let right = admit_machine_graph_ast_v1(&from_json).expect("JSON AST admits");
    assert_eq!(left.identity(), right.identity());
}

#[test]
fn g0_behavior_literal_program_and_json_preserve_graph_and_identity() {
    let graph = admitted_behavior_graph(61);
    let literal = parse_valid_behavior(&graph);
    let decoded = parse_machine_behavior_v1(&literal).expect("behavior syntax decodes");
    assert_eq!(decoded.base_graph(), graph.identity());
    let admitted = decoded.admit_against(&graph).expect("behavior admits");

    let program = write_machine_behavior_program_v1(&admitted).expect("behavior encodes");
    let sexpr_bytes = program
        .print_sexpr_checked()
        .expect("canonical behavior s-expression");
    let json_bytes = program
        .print_json_checked()
        .expect("canonical behavior JSON");
    let from_sexpr = VersionedProgram::parse_sexpr(&sexpr_bytes).expect("s-expression reparses");
    let from_json = VersionedProgram::parse_json(&json_bytes).expect("JSON reparses");
    let admitted_sexpr = parse_machine_behavior_program_v1(&from_sexpr)
        .expect("s-expression behavior decodes")
        .admit_against(&graph)
        .expect("s-expression behavior admits");
    let admitted_json = parse_machine_behavior_program_v1(&from_json)
        .expect("JSON behavior decodes")
        .admit_against(&graph)
        .expect("JSON behavior admits");

    assert_eq!(admitted.identity(), admitted_sexpr.identity());
    assert_eq!(admitted.identity(), admitted_json.identity());
    let canonical = write_machine_behavior_v1(&admitted).expect("behavior writes");
    assert!(canonical.same_shape(from_sexpr.program()));
    assert!(canonical.same_shape(from_json.program()));

    let max_version = behavior_source(&identity_hex(*graph.identity().as_bytes())).replacen(
        "(ref \"values/condition\" \"1\"",
        "(ref \"values/condition\" \"18446744073709551615\"",
        1,
    );
    let max = admit_machine_behavior_ast_v1(
        &sexpr::parse(&max_version).expect("u64::MAX behavior syntax parses"),
        &graph,
    )
    .expect("u64::MAX behavior reference admits");
    assert!(
        write_machine_behavior_program_v1(&max)
            .expect("max-version behavior encodes")
            .print_sexpr()
            .contains("18446744073709551615")
    );

    let max_microstep = behavior_source(&identity_hex(*graph.identity().as_bytes())).replacen(
        "(total-priority \"0\")",
        "(total-priority \"4294967295\")",
        1,
    );
    let max_microstep = admit_machine_behavior_ast_v1(
        &sexpr::parse(&max_microstep).expect("u32::MAX microstep syntax parses"),
        &graph,
    )
    .expect("u32::MAX microstep admits");
    assert!(
        write_machine_behavior_program_v1(&max_microstep)
            .expect("max-microstep behavior encodes")
            .print_sexpr()
            .contains("4294967295")
    );

    let signed_zero = behavior_source(&identity_hex(*graph.identity().as_bytes())).replacen(
        "(bounded 0.01 0.02",
        "(bounded -0.0 0.02",
        1,
    );
    let signed_zero = admit_machine_behavior_ast_v1(
        &sexpr::parse(&signed_zero).expect("signed-zero float syntax parses"),
        &graph,
    )
    .expect("signed zero canonicalizes during semantic admission");
    assert!(
        write_machine_behavior_program_v1(&signed_zero)
            .expect("canonicalized zero behavior encodes")
            .print_sexpr()
            .contains("(bounded 0.0 0.02")
    );
}

#[test]
fn g0_high_cardinality_admitted_behavior_writer_remains_decodable() {
    let graph = admitted_behavior_graph(61);
    let mut high_cardinality = parse_valid_behavior(&graph);
    let root = root_items(&mut high_cardinality);
    {
        let events = section_items(&mut root[5]);
        events.extend((1..1_024).map(generated_event));
        assert_eq!(events.len() - 1, 1_024);
    }
    {
        let tolerances = section_items(&mut root[6]);
        tolerances.extend((3..8_192).map(generated_bounded_tolerance));
        assert_eq!(tolerances.len() - 1, 8_192);
    }
    assert!(ast_node_count(&high_cardinality) > 262_144);
    assert!(ast_node_count(&high_cardinality) < MAX_MACHINE_BEHAVIOR_AST_NODES);

    let admitted = admit_machine_behavior_ast_v1(&high_cardinality, &graph)
        .expect("high-cardinality behavior remains inside semantic bounds");
    let canonical = write_machine_behavior_v1(&admitted).expect("admitted behavior writes");
    assert!(ast_node_count(&canonical) <= MAX_MACHINE_BEHAVIOR_AST_NODES);
    let readmitted = parse_machine_behavior_v1(&canonical)
        .expect("canonical high-cardinality behavior decodes")
        .admit_against(&graph)
        .expect("canonical high-cardinality behavior readmits");
    assert_eq!(admitted.identity(), readmitted.identity());
}

#[test]
fn g3_behavior_source_and_nested_permutations_do_not_move_identity() {
    let graph = admitted_behavior_graph(61);
    let baseline = admit_machine_behavior_ast_v1(&parse_valid_behavior(&graph), &graph)
        .expect("baseline behavior admits")
        .identity();
    let mut permuted = parse_valid_behavior(&graph);
    let root = root_items(&mut permuted);
    for section_index in [2_usize, 3, 4, 5, 6, 7] {
        section_items(&mut root[section_index])[1..].reverse();
    }
    let events = section_items(&mut root[5]);
    let event = root_items(&mut events[1]);
    section_items(&mut event[6])[1..].reverse();
    let reset = root_items(&mut event[7]);
    section_items(&mut reset[2])[1..].reverse();
    let dependences = section_items(&mut root[7]);
    let dependence = root_items(&mut dependences[1]);
    section_items(&mut dependence[1])[1..].reverse();

    let readmitted = admit_machine_behavior_ast_v1(&permuted, &graph)
        .expect("permuted behavior admits")
        .identity();
    assert_eq!(baseline, readmitted);
}

#[test]
fn g0_behavior_codec_refusals_retain_exact_paths() {
    let graph = admitted_behavior_graph(61);
    let source = behavior_source(&identity_hex(*graph.identity().as_bytes()));

    let graph_id = identity_hex(*graph.identity().as_bytes());
    let uppercase_graph_id = format!("A{}", &graph_id[1..]);
    let uppercase_graph = source.replacen(&graph_id, &uppercase_graph_id, 1);
    let error = parse_machine_behavior_v1(
        &sexpr::parse(&uppercase_graph).expect("generic syntax remains valid"),
    )
    .expect_err("uppercase graph identity must refuse");
    assert_eq!(error.rule(), MachineBehaviorCodecRule::InvalidReference);
    assert_eq!(error.code(), "MachineBehaviorCodecInvalidReference");
    assert_eq!(error.path(), "$[1][1]");
    assert!(error.span().end > error.span().start);
    assert!(!error.detail().is_empty());
    assert!(!error.hint().is_empty());

    let guard_digest = digest(20);
    let uppercase_guard_digest = format!("A{}", &guard_digest[1..]);
    let uppercase_guard = source.replacen(&guard_digest, &uppercase_guard_digest, 1);
    let error = parse_machine_behavior_v1(
        &sexpr::parse(&uppercase_guard).expect("generic syntax remains valid"),
    )
    .expect_err("uppercase guard digest must refuse");
    assert_eq!(error.rule(), MachineBehaviorCodecRule::InvalidReference);
    assert_eq!(error.path(), "$[5][1][3][3]");

    let unknown_orientation = source.replacen("negative-to-positive", "sideways", 1);
    let error = parse_machine_behavior_v1(
        &sexpr::parse(&unknown_orientation).expect("generic syntax remains valid"),
    )
    .expect_err("unknown guard orientation must refuse");
    assert_eq!(error.rule(), MachineBehaviorCodecRule::UnknownTag);
    assert_eq!(error.path(), "$[5][1][4]");

    let integer_scale = source.replacen("(random 0.1", "(random 1", 1);
    let error = parse_machine_behavior_v1(
        &sexpr::parse(&integer_scale).expect("generic syntax remains valid"),
    )
    .expect_err("integer is not the canonical finite-float shape");
    assert_eq!(error.rule(), MachineBehaviorCodecRule::UnexpectedForm);
    assert_eq!(error.path(), "$[6][1][6][1]");

    let leading_zero = source.replacen("(total-priority \"0\")", "(total-priority \"00\")", 1);
    let error = parse_machine_behavior_v1(
        &sexpr::parse(&leading_zero).expect("generic syntax remains valid"),
    )
    .expect_err("noncanonical microstep must refuse");
    assert_eq!(error.rule(), MachineBehaviorCodecRule::InvalidNumber);
    assert_eq!(error.path(), "$[5][1][8][1]");

    let mut forged = parse_valid_behavior(&graph);
    let root = root_items(&mut forged);
    let events = section_items(&mut root[5]);
    let event = root_items(&mut events[1]);
    let guard = root_items(&mut event[3]);
    guard[3] = Node::synthetic(NodeKind::Float(f64::NAN));
    let error = parse_machine_behavior_v1(&forged).expect_err("forged nested AST must refuse");
    assert_eq!(error.rule(), MachineBehaviorCodecRule::InvalidAst);
    assert_eq!(error.path(), "$[5][1][3][3]");
}

#[test]
fn g0_behavior_resource_preflight_wins_before_recursive_validation() {
    assert!(MAX_MACHINE_BEHAVIOR_AST_NODES > MAX_MACHINE_BEHAVIOR_NESTED_REFERENCES);
    let graph = admitted_behavior_graph(61);
    let mut oversized_events = parse_valid_behavior(&graph);
    let root = root_items(&mut oversized_events);
    let events = section_items(&mut root[5]);
    events.push(Node::synthetic(NodeKind::Float(f64::NAN)));
    events.extend(
        (0..MAX_MACHINE_BEHAVIOR_EVENTS).map(|_| Node::synthetic(NodeKind::List(Vec::new()))),
    );
    let error = parse_machine_behavior_v1(&oversized_events)
        .expect_err("event section cap must precede recursive AST validation");
    assert_eq!(error.rule(), MachineBehaviorCodecRule::ResourceLimit);
    assert_eq!(error.path(), "$[5]");

    let mut aggregate_overflow = parse_valid_behavior(&graph);
    let root = root_items(&mut aggregate_overflow);
    {
        let conditions = section_items(&mut root[3]);
        let condition = root_items(&mut conditions[3]);
        let source = root_items(&mut condition[6]);
        let reset_events = root_items(&mut source[2]);
        let missing = 8_192 - (reset_events.len() - 1);
        reset_events.extend((0..missing).map(|_| synthetic_string("event/contact")));
    }
    {
        let events = section_items(&mut root[5]);
        let event = root_items(&mut events[1]);
        let dependencies = section_items(&mut event[6]);
        let missing = 8_192 - (dependencies.len() - 1);
        dependencies.extend((0..missing).map(|_| Node::synthetic(NodeKind::List(Vec::new()))));
        let reset = root_items(&mut event[7]);
        let writes = section_items(&mut reset[2]);
        let missing = 8_192 - (writes.len() - 1);
        writes.extend((0..missing).map(|_| synthetic_string("state/position")));
    }
    {
        let dependences = section_items(&mut root[7]);
        let dependence = root_items(&mut dependences[1]);
        let members = section_items(&mut dependence[1]);
        members.push(Node::synthetic(NodeKind::Float(f64::NAN)));
        let missing = 8_193 - (members.len() - 1);
        members.extend((0..missing).map(|_| Node::synthetic(NodeKind::List(Vec::new()))));
    }
    let error = parse_machine_behavior_v1(&aggregate_overflow)
        .expect_err("aggregate nested cap must precede recursive AST validation");
    assert_eq!(error.rule(), MachineBehaviorCodecRule::ResourceLimit);
    assert_eq!(error.path(), "$[7][1][1]");

    for (set_valued, expected_path) in [(false, "$[5][1][7][2]"), (true, "$[5][1][7][3]")] {
        let mut reset_overflow = parse_valid_behavior(&graph);
        let root = root_items(&mut reset_overflow);
        let events = section_items(&mut root[5]);
        let event = root_items(&mut events[1]);
        {
            let dependencies = section_items(&mut event[6]);
            let missing = 16_384 - (dependencies.len() - 1);
            dependencies.extend((0..missing).map(|_| Node::synthetic(NodeKind::List(Vec::new()))));
        }
        let reset = root_items(&mut event[7]);
        let writes_index = if set_valued {
            reset[0] = synthetic_symbol("set-valued");
            reset.insert(2, synthetic_ref("resets/generated-outcomes", 76));
            3
        } else {
            2
        };
        let writes = section_items(&mut reset[writes_index]);
        writes.push(Node::synthetic(NodeKind::Float(f64::NAN)));
        let missing = 16_385 - (writes.len() - 1);
        writes.extend((0..missing).map(|_| synthetic_string("state/position")));

        let error = parse_machine_behavior_v1(&reset_overflow)
            .expect_err("aggregate cap must report the exact reset-writes container");
        assert_eq!(error.rule(), MachineBehaviorCodecRule::ResourceLimit);
        assert_eq!(error.path(), expected_path);
    }
}

#[test]
fn g0_behavior_artifact_cannot_rebind_or_bypass_semantic_admission() {
    let graph_a = admitted_behavior_graph(61);
    let graph_b = admitted_behavior_graph(63);
    assert_ne!(graph_a.identity(), graph_b.identity());
    let mut incomplete = parse_valid_behavior(&graph_a);
    let root = root_items(&mut incomplete);
    let conditions = section_items(&mut root[3]);
    conditions.remove(1);

    let error = admit_machine_behavior_ast_v1(&incomplete, &graph_b)
        .expect_err("graph mismatch must precede inspection of invalid behavior semantics");
    assert_eq!(error.code(), "MachineBehaviorBaseGraphMismatch");
    let MachineBehaviorAstAdmissionError::BaseGraphMismatch { declared, provided } = error else {
        panic!("wrong graph must fail at the explicit binding boundary")
    };
    assert_eq!(declared, graph_a.identity());
    assert_eq!(provided, graph_b.identity());

    let error = admit_machine_behavior_ast_v1(&incomplete, &graph_a)
        .expect_err("missing initial condition must refuse semantically");
    let MachineBehaviorAstAdmissionError::Behavior(refusal) = error else {
        panic!("valid syntax must preserve behavior semantic refusal")
    };
    assert!(
        refusal
            .findings()
            .iter()
            .any(|finding| finding.rule() == MachineBehaviorRule::MissingInitialCondition)
    );
}

#[test]
fn g3_generic_json_shape_is_not_a_second_behavior_grammar() {
    let graph = admitted_behavior_graph(61);
    let node = parse_valid_behavior(&graph);
    let json_bytes = json::print(&node).expect("generic JSON prints");
    let from_json = json::parse(&json_bytes).expect("generic JSON reparses");
    let left = admit_machine_behavior_ast_v1(&node, &graph).expect("s-expression AST admits");
    let right = admit_machine_behavior_ast_v1(&from_json, &graph).expect("JSON AST admits");
    assert_eq!(left.identity(), right.identity());
}
