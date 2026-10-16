//! G0/G3 coverage for the bounded Part-21 syntax/graph kernel.

use fs_io::{
    IoError, StepEntity, StepLimits, StepProfileHint, StepValue, parse_step,
    parse_step_with_limits, write_step, write_step_with_limits,
};

fn prefix(schema: &str) -> String {
    format!(
        "ISO-10303-21;\n\
         HEADER;\n\
         FILE_DESCRIPTION(('bounded syntax fixture'),'2;1');\n\
         FILE_NAME('fixture.step','2026-07-16T00:00:00',('fs-io'),('FrankenSim'),\
         'fs-io','FrankenSim','');\n\
         FILE_SCHEMA(('{schema}'));\n\
         ENDSEC;\n\
         DATA;\n"
    )
}

fn document(schema: &str, data: &str) -> String {
    format!("{}{data}\nENDSEC;\nEND-ISO-10303-21;\n", prefix(schema))
}

#[test]
fn step_001_forward_references_complex_entities_and_strings_round_trip() {
    let source = document(
        "CONFIG_CONTROL_DESIGN",
        "#20=NEXT(#10,(1,.T.,LENGTH_MEASURE(2.5E+0)),'O''Brien');\n\
         #10=ORIGIN($,*);\n\
         #30=(ALPHA(#10)BETA((#20,#30)));",
    );
    let parsed = parse_step(source.as_bytes()).expect("bounded Part-21 parses");

    assert_eq!(parsed.receipt().profile_hint(), StepProfileHint::Ap203);
    parsed
        .document()
        .require_declared_schema("config_control_design")
        .expect("schema admission is exact apart from ASCII case");
    assert!(matches!(
        parsed
            .document()
            .require_declared_schema("AUTOMOTIVE_DESIGN"),
        Err(IoError::Unsupported { .. })
    ));
    assert_eq!(parsed.receipt().instance_count(), 3);
    assert_eq!(parsed.receipt().reference_count(), 4);
    assert_eq!(parsed.receipt().max_instance_id(), 30);
    assert_eq!(
        parsed.receipt().schema_identifiers(),
        &["CONFIG_CONTROL_DESIGN".to_string()]
    );

    let canonical = String::from_utf8(write_step(parsed.document()).expect("writes"))
        .expect("canonical output is ASCII");
    let at_10 = canonical.find("#10=ORIGIN").expect("#10 emitted");
    let at_20 = canonical.find("#20=NEXT").expect("#20 emitted");
    let at_30 = canonical.find("#30=(ALPHA").expect("#30 emitted");
    assert!(at_10 < at_20 && at_20 < at_30, "instances sorted by ID");
    assert!(canonical.contains("2.5E+0"), "numeric spelling preserved");
    assert!(canonical.contains("'O''Brien'"), "apostrophe doubled");

    let reparsed = parse_step(canonical.as_bytes()).expect("canonical bytes reparse");
    assert_eq!(
        write_step(reparsed.document()).expect("rewrites"),
        canonical.as_bytes()
    );
    let json = parsed.receipt().to_json();
    assert!(json.contains("\"authority\":\"syntax-only\""));
    assert!(json.contains("\"syntax_version\":\"part21-syntax-v1\""));
    assert!(json.contains("\"source_fingerprint_fnv1a64\""));
    assert!(json.contains("\"limits\""));
    assert!(json.contains("\"profile_hint\":\"ap203\""));
    assert!(json.contains("no EXPRESS-schema"));
}

#[test]
fn step_002_canonical_data_order_is_permutation_invariant() {
    let first = document("AUTOMOTIVE_DESIGN", "#2=SECOND(#1);\n#1=FIRST('same');");
    let second = document("AUTOMOTIVE_DESIGN", "#1=FIRST('same');\n#2=SECOND(#1);");
    let first = parse_step(first.as_bytes()).expect("first order parses");
    let second = parse_step(second.as_bytes()).expect("second order parses");

    assert_eq!(first.receipt().profile_hint(), StepProfileHint::Ap214);
    assert_eq!(write_step(first.document()), write_step(second.document()));
    assert_ne!(
        first.receipt().source_fingerprint(),
        second.receipt().source_fingerprint(),
        "fixture source fingerprints retain declaration order"
    );
    assert_eq!(
        first.receipt().canonical_layout_fingerprint(),
        second.receipt().canonical_layout_fingerprint(),
        "canonical-layout fingerprint ignores DATA declaration order"
    );
}

#[test]
fn step_003_duplicate_and_dangling_instance_graphs_refuse() {
    let duplicate = document("CONFIG_CONTROL_DESIGN", "#1=ONE();\n#1=AGAIN();");
    let error = parse_step(duplicate.as_bytes()).expect_err("duplicate must refuse");
    assert!(
        matches!(&error, IoError::Malformed { what, .. } if what.contains("duplicate")),
        "{error:?}"
    );

    let dangling = document("CONFIG_CONTROL_DESIGN", "#1=ONE(#99);");
    let error = parse_step(dangling.as_bytes()).expect_err("dangling must refuse");
    assert!(
        matches!(&error, IoError::Malformed { what, .. } if what.contains("dangling") && what.contains("#99")),
        "{error:?}"
    );

    let zero = document("CONFIG_CONTROL_DESIGN", "#0=ZERO();");
    assert!(matches!(
        parse_step(zero.as_bytes()),
        Err(IoError::Malformed { .. })
    ));
}

#[test]
fn step_004_envelope_and_lexical_truncations_fail_closed() {
    let valid = document("CONFIG_CONTROL_DESIGN", "#1=ONE('ok');");
    for malformed in [
        valid.replace("ISO-10303-21;", "ISO-10303-21"),
        valid.replace("HEADER;", "HEAD;"),
        valid.replace("FILE_NAME", "FILE_NAM"),
        valid.replace("ENDSEC;\nDATA;", "DATA;"),
        valid.replace("#1=ONE", "#1 ONE"),
        valid.replace("END-ISO-10303-21;", ""),
        format!("{}/* never closed", valid),
        document("CONFIG_CONTROL_DESIGN", "#1=ONE('never closes);"),
    ] {
        assert!(
            parse_step(malformed.as_bytes()).is_err(),
            "malformed input unexpectedly parsed: {malformed:?}"
        );
    }

    for cut in (0..valid.len()).step_by(11) {
        assert!(
            parse_step(&valid.as_bytes()[..cut]).is_err(),
            "truncated prefix {cut} unexpectedly parsed"
        );
    }
}

#[test]
fn step_005_comments_profile_hints_and_uppercase_grammar_are_deterministic() {
    let source = "/*lead*/ISO-10303-21;HEADER;\
        FILE_DESCRIPTION(/*a*/('x'),'2;1');\
        FILE_NAME('f','t',('a'),('o'),'p','s','');\
        FILE_SCHEMA(('AP203_DEMO','AP214_DEMO'));ENDSEC;DATA;\
        #2=THING(#1);#1=BASE();ENDSEC;END-ISO-10303-21;/*tail*/";
    let parsed = parse_step(source.as_bytes()).expect("commented uppercase syntax parses");
    assert_eq!(parsed.receipt().profile_hint(), StepProfileHint::Ambiguous);
    let canonical = String::from_utf8(write_step(parsed.document()).expect("writes"))
        .expect("canonical output is ASCII");
    assert!(canonical.starts_with("ISO-10303-21;\nHEADER;\n"));
    assert!(canonical.contains("#1=BASE();\n#2=THING(#1);"));
    assert!(!canonical.contains("/*"));

    let lowercase_keyword = source.replace("HEADER", "header");
    assert!(matches!(
        parse_step(lowercase_keyword.as_bytes()),
        Err(IoError::Malformed { .. })
    ));
}

#[test]
fn step_006_resource_caps_refuse_before_partial_admission() {
    let source = document(
        "CONFIG_CONTROL_DESIGN",
        "#1=ONE(((1)));\n#2=TWO('long string');",
    );

    let limits = StepLimits {
        max_bytes: source.len() - 1,
        ..StepLimits::default()
    };
    assert!(matches!(
        parse_step_with_limits(source.as_bytes(), limits),
        Err(IoError::ResourceBound { .. })
    ));

    let limits = StepLimits {
        max_instances: 1,
        ..StepLimits::default()
    };
    assert!(matches!(
        parse_step_with_limits(source.as_bytes(), limits),
        Err(IoError::ResourceBound { .. })
    ));

    let limits = StepLimits {
        max_nesting: 2,
        ..StepLimits::default()
    };
    assert!(matches!(
        parse_step_with_limits(source.as_bytes(), limits),
        Err(IoError::ResourceBound { .. })
    ));

    let limits = StepLimits {
        max_nesting: usize::MAX,
        ..StepLimits::default()
    };
    assert!(matches!(
        parse_step_with_limits(source.as_bytes(), limits),
        Err(IoError::ResourceBound { .. })
    ));

    let limits = StepLimits {
        max_string_bytes: 8,
        ..StepLimits::default()
    };
    assert!(matches!(
        parse_step_with_limits(source.as_bytes(), limits),
        Err(IoError::ResourceBound { .. })
    ));

    let number_source = document("CONFIG_CONTROL_DESIGN", "#1=ONE(1234);");
    let limits = StepLimits {
        max_number_bytes: 3,
        ..StepLimits::default()
    };
    assert!(matches!(
        parse_step_with_limits(number_source.as_bytes(), limits),
        Err(IoError::ResourceBound { .. })
    ));

    let component_source = document("CONFIG_CONTROL_DESIGN", "#1=(ONE()TWO());");
    let limits = StepLimits {
        max_components_per_instance: 1,
        ..StepLimits::default()
    };
    assert!(matches!(
        parse_step_with_limits(component_source.as_bytes(), limits),
        Err(IoError::ResourceBound { .. })
    ));

    let limits = StepLimits {
        max_tokens: 12,
        ..StepLimits::default()
    };
    assert!(matches!(
        parse_step_with_limits(source.as_bytes(), limits),
        Err(IoError::ResourceBound { .. })
    ));

    let limits = StepLimits {
        max_values: 4,
        ..StepLimits::default()
    };
    assert!(matches!(
        parse_step_with_limits(source.as_bytes(), limits),
        Err(IoError::ResourceBound { .. })
    ));
}

#[test]
fn step_007_staged_encodings_and_binary_literals_refuse_explicitly() {
    let encoded = document("CONFIG_CONTROL_DESIGN", r"#1=ONE('a\X2\0041\X0\');");
    assert!(matches!(
        parse_step(encoded.as_bytes()),
        Err(IoError::Unsupported { .. })
    ));

    let binary = document("CONFIG_CONTROL_DESIGN", "#1=ONE(\"0F\");");
    assert!(matches!(
        parse_step(binary.as_bytes()),
        Err(IoError::Unsupported { .. })
    ));

    let mut non_ascii = document("CONFIG_CONTROL_DESIGN", "#1=ONE('ok');").into_bytes();
    let position = non_ascii
        .iter()
        .position(|byte| *byte == b'o')
        .expect("fixture has o");
    non_ascii[position] = 0xff;
    assert!(matches!(
        parse_step(&non_ascii),
        Err(IoError::Unsupported { .. })
    ));
}

#[test]
fn step_008_writer_revalidates_caller_constructed_documents() {
    let source = document("CONFIG_CONTROL_DESIGN", "#1=ONE();\n#2=TWO(#1);");
    let parsed = parse_step(source.as_bytes()).expect("fixture parses");

    let mut duplicate = parsed.document().clone();
    duplicate.instances.push(duplicate.instances[0].clone());
    assert!(matches!(
        write_step(&duplicate),
        Err(IoError::Malformed { .. })
    ));

    let mut dangling = parsed.document().clone();
    dangling.instances[1].components[0]
        .parameters
        .push(StepValue::Reference(99));
    assert!(matches!(
        write_step(&dangling),
        Err(IoError::Malformed { .. })
    ));

    let mut bad_number = parsed.document().clone();
    bad_number.instances[0].components[0]
        .parameters
        .push(StepValue::Number("1e".to_string()));
    assert!(matches!(
        write_step(&bad_number),
        Err(IoError::Malformed { .. })
    ));

    let mut duplicate_component = parsed.document().clone();
    duplicate_component.instances[0].components = vec![
        StepEntity {
            name: "A".to_string(),
            parameters: vec![],
        },
        StepEntity {
            name: "a".to_string(),
            parameters: vec![],
        },
    ];
    assert!(matches!(
        write_step(&duplicate_component),
        Err(IoError::Malformed { .. })
    ));

    let limits = StepLimits {
        max_bytes: 32,
        ..StepLimits::default()
    };
    assert!(matches!(
        write_step_with_limits(parsed.document(), limits),
        Err(IoError::ResourceBound { .. })
    ));

    let limits = StepLimits {
        max_tokens: 8,
        ..StepLimits::default()
    };
    assert!(matches!(
        write_step_with_limits(parsed.document(), limits),
        Err(IoError::ResourceBound { .. })
    ));
}

#[test]
fn step_009_file_schema_shape_and_value_grammar_refuse() {
    let wrong_schema_shape = document("CONFIG_CONTROL_DESIGN", "#1=ONE();").replace(
        "FILE_SCHEMA(('CONFIG_CONTROL_DESIGN'))",
        "FILE_SCHEMA('CONFIG_CONTROL_DESIGN')",
    );
    assert!(matches!(
        parse_step(wrong_schema_shape.as_bytes()),
        Err(IoError::Malformed { .. })
    ));

    let bad_description = document("CONFIG_CONTROL_DESIGN", "#1=ONE();").replace(
        "FILE_DESCRIPTION(('bounded syntax fixture'),'2;1')",
        "FILE_DESCRIPTION() ",
    );
    assert!(matches!(
        parse_step(bad_description.as_bytes()),
        Err(IoError::Malformed { .. })
    ));

    let bad_implementation_level =
        document("CONFIG_CONTROL_DESIGN", "#1=ONE();").replace("'2;1'", "'edition-two'");
    assert!(matches!(
        parse_step(bad_implementation_level.as_bytes()),
        Err(IoError::Malformed { .. })
    ));

    let bad_name = document("CONFIG_CONTROL_DESIGN", "#1=ONE();").replace(
        "FILE_NAME('fixture.step','2026-07-16T00:00:00',('fs-io'),('FrankenSim'),'fs-io','FrankenSim','')",
        "FILE_NAME(1)",
    );
    assert!(matches!(
        parse_step(bad_name.as_bytes()),
        Err(IoError::Malformed { .. })
    ));

    for value in [
        "1E",
        "1E05",
        ".5",
        "+",
        ".BAD",
        "1..2",
        "(#1,)",
        "TYPE()",
        "TYPE(1,2)",
        "1.0e+2",
    ] {
        let source = document("CONFIG_CONTROL_DESIGN", &format!("#1=ONE({value});"));
        assert!(
            parse_step(source.as_bytes()).is_err(),
            "invalid value unexpectedly parsed: {value}"
        );
    }
}
