use super::*;

const STUDY: &str = r#"(fsim-uncertainty-study
    :version 1 :project "cooling-reference.fsim" :samples 8 :seed 29
    :wall-time 120s :method monte-carlo :correlation independent :qoi "temperature-max"
    :geometry ((mesh :role "enclosure" :path "plate.stl" :unit "m" :max-hole-edges 0))
    :materials ("aa6061.fsmcdpk") :interfaces ()
    :parameters (
        (uniform :name "power" :target power :entity "air" :low 4W :high 6W)
        (uniform :name "ambient" :target convection-temperature :entity "air" :low 294K :high 300K)))"#;

fn base() -> ProjectSpec {
    crate::parse_sexpr_migrating(include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
        "/../../data/reference-project/cooling-reference.fsim"))).unwrap().decoded.spec
}

#[test]
fn probability_study_roundtrips_through_existing_typed_ast() {
    let study = UncertaintyStudy::parse(STUDY).unwrap();
    assert_eq!(UncertaintyStudy::parse(study.canonical()).unwrap(), study);
    assert_eq!(study.samples(), 8);
    assert_eq!(study.seed(), 29);
    assert_eq!(study.parameters()[0].target.unit(), "W");
    assert_eq!(study.parameters()[0].low, 4.0);
    assert_eq!(study.geometry()[0].max_hole_edges, 0);
    assert_eq!(UncertaintyStudy::parse(&format!("; study comment\n{STUDY}")).unwrap(), study);
}

#[test]
fn unknown_repeated_implicit_and_out_of_range_declarations_refuse() {
    for (from, to) in [
        (":version 1", ":version 2"),
        (":samples 8", ":samples 1"),
        (":samples 8", ":samples 257"),
        (":seed 29", ":seed -1"),
        (":seed 29", ":seed 29 :seed 30"),
        (":wall-time 120s", ":wall-time 120"),
        (":correlation independent", ":correlation unknown"),
        (":method monte-carlo", ":method quasi-monte-carlo"),
        (":low 4W", ":low 4K"),
        (":low 4W", ":low -4W"),
        (":high 6W", ":high 3W"),
        ("uniform :name", "interval :name"),
        (":entity \"air\" :low 4W", ":entity \"air\" :surrogate true :low 4W"),
    ] {
        assert!(STUDY.contains(from));
        assert!(UncertaintyStudy::parse(&STUDY.replace(from, to)).is_err(), "admitted {to}");
    }
    assert!(UncertaintyStudy::parse(&" ".repeat(MAX_SOURCE_BYTES + 1)).is_err());
}

#[test]
fn native_samples_change_only_bound_inputs_and_do_not_accumulate() {
    let original = base();
    let before = crate::print_sexpr(&original);
    let bound = UncertaintyStudy::parse(STUDY).unwrap().bind(&original).unwrap();
    let first = bound.sample_project(&[4.5, 297.0]).unwrap();
    let second = bound.sample_project(&[5.5, 299.0]).unwrap();
    assert_eq!(first.power.as_ref().unwrap()[0].watts.value, 4.5);
    assert_eq!(second.power.as_ref().unwrap()[0].watts.value, 5.5);
    assert_eq!(first.power.as_ref().unwrap()[0].duty, original.power.as_ref().unwrap()[0].duty);
    let ThermalBoundaryCondition::Convection { reference_temperature, .. } =
        &first.cooling.as_ref().unwrap().conduction.as_ref().unwrap().boundaries[0].condition else { panic!() };
    assert_eq!(reference_temperature.value, 297.0);
    assert_eq!(first.geometry, original.geometry);
    assert_eq!(first.materials, original.materials);
    assert_eq!(first.requirements, original.requirements);
    assert_eq!(first.seeds, original.seeds);
    assert_eq!(first.budgets, original.budgets);
    assert_eq!(first.solver, original.solver);
    assert_eq!(first, bound.sample_project(&[4.5, 297.0]).unwrap());
    assert_eq!(crate::print_sexpr(&original), before);
    assert_eq!(bound.threshold_k(), 348.15);
    let decoded = crate::parse_sexpr(&crate::print_sexpr(&first)).unwrap();
    assert!(decoded.findings().is_empty());
    assert_eq!(decoded.spec, first);
}

#[test]
fn probability_support_cannot_expand_project_domains_or_select_missing_fields() {
    for changed in [
        STUDY.replace(":high 300K", ":high 400K"),
        STUDY.replace(":target power", ":target heat-flux"),
        STUDY.replace(":entity \"air\"", ":entity \"missing\""),
        STUDY.replace(":role \"enclosure\"", ":role \"unknown\""),
    ] {
        assert!(UncertaintyStudy::parse(&changed).and_then(|s| s.bind(&base())).is_err());
    }
    let bound = UncertaintyStudy::parse(STUDY).unwrap().bind(&base()).unwrap();
    for sample in [vec![], vec![4.0], vec![3.99, 295.0], vec![4.0, 301.0], vec![f64::NAN, 295.0]] {
        assert!(bound.sample_project(&sample).is_err());
    }
    let mut signoff = base();
    signoff.metadata.as_mut().unwrap().decision_gate = DecisionGate::ComplianceSignoff;
    assert!(UncertaintyStudy::parse(STUDY).unwrap().bind(&signoff).is_err());
}

#[test]
fn duplicate_physical_targets_refuse_even_under_distinct_statistical_names() {
    let second = "(uniform :name \"alias\" :target power :entity \"air\" :low 4W :high 6W)";
    let old = "(uniform :name \"ambient\" :target convection-temperature :entity \"air\" :low 294K :high 300K)";
    assert!(UncertaintyStudy::parse(&STUDY.replace(old, second)).is_err());
    assert!(UncertaintyStudy::parse(&STUDY.replace(":name \"ambient\"", ":name \"power\"")).is_err());
}

#[test]
fn air_inlet_binding_moves_the_whole_named_branch_but_not_other_branches() {
    let mut project = base();
    let setup = project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    let law = |branch: &str, order| ThermalBoundaryCondition::AirflowConvection {
        branch: branch.into(), order,
        inlet_temperature: QtyAny::new(295.0, crate::spec::dims::TEMPERATURE),
        hydraulic_diameter: QtyAny::new(0.01, crate::spec::dims::LENGTH),
        flow_area: QtyAny::new(0.001, crate::spec::dims::AREA),
        channel_length: QtyAny::new(0.1, crate::spec::dims::LENGTH),
        correlation: "convection.circular-duct-laminar-cwt".into(),
    };
    setup.boundaries = vec![
        crate::ThermalBoundary { target: "first".into(), condition: law("a", 0) },
        crate::ThermalBoundary { target: "second".into(), condition: law("a", 1) },
        crate::ThermalBoundary { target: "other".into(), condition: law("b", 0) },
    ];
    // Direct field-binding test; these isolated boundary declarations are not
    // presented as an admitted geometric model or an executed airflow solve.
    let parameter = UniformParameter { name: "inlet".into(), target: Target::AirInletTemperature,
        entity: "a".into(), low: 294.0, high: 300.0 };
    apply(&mut project, &parameter, 298.0).unwrap();
    let found = project.cooling.unwrap().conduction.unwrap().boundaries.iter().map(|row| {
        let ThermalBoundaryCondition::AirflowConvection { inlet_temperature, .. } = &row.condition else { panic!() };
        inlet_temperature.value
    }).collect::<Vec<_>>();
    assert_eq!(found, [298.0, 298.0, 295.0]);
}
