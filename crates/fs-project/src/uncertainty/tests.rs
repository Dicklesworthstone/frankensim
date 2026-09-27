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
    assert_eq!(study.compliance(), None);
    assert_eq!(study.seed(), 29);
    assert_eq!(study.parameters()[0].target.unit(), "W");
    assert_eq!(study.parameters()[0].low, 4.0);
    assert_eq!(study.geometry()[0].max_hole_edges, 0);
    assert_eq!(UncertaintyStudy::parse(&format!("; study comment\n{STUDY}")).unwrap(), study);
}

#[test]
fn unknown_repeated_implicit_and_out_of_range_declarations_refuse() {
    for (from, to) in [
        (":version 1", ":version 3"),
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
    let decoded = crate::parse_sexpr(&crate::print_sexpr(&first).unwrap()).unwrap();
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

fn fan_study() -> String {
    STUDY.replace(
        "(uniform :name \"power\" :target power :entity \"air\" :low 4W :high 6W)",
        "(uniform :name \"speed\" :target fan-speed-ratio :entity \"fixture-bank\" :low 0.6 :high 1.4)",
    )
}

#[test]
fn fan_speed_is_an_absolute_ratio_on_only_the_named_retained_bank() {
    let mut original = base();
    let system = original.cooling.as_mut().unwrap().fan_system.as_mut().unwrap();
    system.banks[0].speed_ratio = 0.8;
    let mut other = system.banks[0].clone();
    other.bank_id = "other-bank".into();
    other.speed_ratio = 1.2;
    system.banks.push(other);
    system.topology = crate::fansystem::FanSystemTopology::Series(
        vec!["fixture-bank".into(), "other-bank".into()]);
    let canonical = crate::print_sexpr(&original).unwrap();
    let study = UncertaintyStudy::parse(&fan_study()).unwrap();
    assert_eq!(study.parameters()[0].target, Target::FanSpeedRatio);
    assert_eq!(study.parameters()[0].target.unit(), "1");
    assert_eq!(UncertaintyStudy::parse(study.canonical()).unwrap(), study);
    let bound = study.bind(&original).unwrap();
    for ratio in [0.6, 1.0, 1.4, 0.6] {
        let sample = bound.sample_project(&[ratio, 297.0]).unwrap();
        let mut expected = original.clone();
        let cooling = expected.cooling.as_mut().unwrap();
        cooling.fan_system.as_mut().unwrap().banks[0].speed_ratio = ratio;
        let ThermalBoundaryCondition::Convection { reference_temperature, .. } =
            &mut cooling.conduction.as_mut().unwrap().boundaries[0].condition else { panic!() };
        reference_temperature.value = 297.0;
        assert_eq!(sample, expected, "only the bank speed and declared ambient may change");
        let wire = crate::print_sexpr(&sample).unwrap();
        assert_eq!(crate::parse_sexpr(&wire).unwrap().spec, sample);
        crate::fansystem::lower_fan_system(sample.cooling.as_ref().unwrap().fan_system.as_ref().unwrap())
            .expect("each sampled bank lowers through the ordinary production fan laws");
    }
    assert_eq!(crate::print_sexpr(&original).unwrap(), canonical);
    assert_eq!(bound.base(), &original);
}

#[test]
fn fan_speed_support_cannot_expand_domains_repair_bad_banks_or_infer_legacy_fans() {
    let source = fan_study();
    for (from, to) in [
        (":low 0.6", ":low 0.49"),
        (":high 1.4", ":high 2.01"),
        (":entity \"fixture-bank\"", ":entity \"unknown-bank\""),
    ] {
        assert!(UncertaintyStudy::parse(&source.replace(from, to)).unwrap().bind(&base()).is_err());
    }
    for defect in 0..4 {
        let mut original = base();
        let cooling = original.cooling.as_mut().unwrap();
        if defect == 0 {
            cooling.fan_system = None;
        } else {
            let system = cooling.fan_system.as_mut().unwrap();
            match defect {
                1 => system.banks.push(system.banks[0].clone()),
                2 => system.banks[0].speed_ratio_domain = (0.0, 2.0),
                _ => system.banks[0].speed_ratio = 3.0,
            }
        }
        assert!(UncertaintyStudy::parse(&source).unwrap().bind(&original).is_err(), "defect {defect}");
    }
    let bound = UncertaintyStudy::parse(&source).unwrap().bind(&base()).unwrap();
    for ratio in [0.59, 1.41, f64::NAN, f64::INFINITY] {
        assert!(bound.sample_project(&[ratio, 297.0]).is_err());
    }
}

#[test]
fn fan_speed_requires_positive_finite_dimensionless_support_and_unique_target() {
    let source = fan_study();
    for value in ["0", "-0.1", "0.6K", "0.6W", "\"0.6\"", "NaN", "1e999"] {
        assert!(UncertaintyStudy::parse(&source.replace(":low 0.6", &format!(":low {value}"))).is_err(),
            "admitted speed {value}");
    }
    let fixed = source.replace(":low 0.6 :high 1.4", ":low 1 :high 1");
    assert!(UncertaintyStudy::parse(&fixed).unwrap().bind(&base()).is_ok());
    let duplicate = source.replace(
        "(uniform :name \"ambient\" :target convection-temperature :entity \"air\" :low 294K :high 300K)",
        "(uniform :name \"second-speed\" :target fan-speed-ratio :entity \"fixture-bank\" :low 0.8 :high 1.2)",
    );
    assert!(UncertaintyStudy::parse(&duplicate).is_err());
}

fn qmc_study() -> String {
    STUDY.replace(":method monte-carlo",
        ":method quasi-monte-carlo :qmc (owen-scrambled-sobol :replicates 2 :samples-per-replicate 4)")
}

#[test]
fn randomized_qmc_layout_roundtrips_without_changing_fixed_count_or_native_binding() {
    let study = UncertaintyStudy::parse(&qmc_study()).unwrap();
    assert_eq!(study.qmc(), Some(QmcLayout { replicates: 2, samples_per_replicate: 4 }));
    assert_eq!(study.compliance(), None);
    assert_eq!(UncertaintyStudy::parse(study.canonical()).unwrap(), study);
    let mc = UncertaintyStudy::parse(STUDY).unwrap().bind(&base()).unwrap();
    let qmc = study.bind(&base()).unwrap();
    assert_eq!(mc.sample_project(&[5.0, 297.0]).unwrap(), qmc.sample_project(&[5.0, 297.0]).unwrap());
    let three = qmc_study().replace(":samples 8", ":samples 12").replace(":replicates 2", ":replicates 3");
    assert_eq!(UncertaintyStudy::parse(&three).unwrap().qmc().unwrap().replicates, 3,
        "only the point count, not the independent replicate count, must be a power of two");
}

#[test]
fn qmc_refuses_missing_extraneous_incompatible_or_unbalanced_layouts() {
    let source = qmc_study();
    for (from, to) in [
        (":method quasi-monte-carlo", ":method monte-carlo"),
        (":replicates 2", ":replicates 1"),
        (":replicates 2", ":replicates 3"),
        (":replicates 2", ":replicates 257"),
        (":samples-per-replicate 4", ":samples-per-replicate 3"),
        (":samples-per-replicate 4", ":samples-per-replicate 0"),
        (":samples-per-replicate 4", ":samples-per-replicate 4 :samples-per-replicate 4"),
        (":samples-per-replicate 4", ":samples-per-replicate 4 :skip-origin true"),
        ("owen-scrambled-sobol", "sobol"),
        (":qmc (owen-scrambled-sobol :replicates 2 :samples-per-replicate 4)", ""),
        (":version 1", ":version 2 :compliance (bernoulli-mixture :required-probability 0.9 :alpha 0.05 :min-samples 2)"),
    ] {
        assert!(source.contains(from));
        assert!(UncertaintyStudy::parse(&source.replace(from, to)).is_err(), "admitted {to}");
    }
    let eleven = (0..11).map(|index| format!(
        "(uniform :name \"p{index}\" :target power :entity \"r{index}\" :low 1W :high 2W)"
    )).collect::<Vec<_>>().join(" ");
    let start = source.find(":parameters (").unwrap() + ":parameters (".len();
    let excessive = format!("{}{eleven}))", &source[..start]);
    assert!(UncertaintyStudy::parse(&excessive).is_err(), "no undeclared Monte Carlo fallback after Sobol dimension ten");
    let fixed = source.replace(":low 4W :high 6W", ":low 5W :high 5W")
        .replace(":low 294K :high 300K", ":low 297K :high 297K");
    assert!(UncertaintyStudy::parse(&fixed).unwrap().bind(&base()).is_ok());
}
