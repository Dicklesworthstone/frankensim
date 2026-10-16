//! Independent ambient-fluid and radiative-reservoir inputs of the native model.
use super::*;
use crate::spec::dims;

const STUDY: &str = r#"(fsim-uncertainty-study
 :version 3 :project "cooling-radiation.fsim" :samples 4 :seed 29
 :wall-time 120s :method monte-carlo :correlation independent :qoi "temperature-max"
 :geometry ((mesh :role "enclosure" :path "plate.stl" :unit "m" :max-hole-edges 0))
 :materials ("aa6061.fsmcdpk" "gray-surface.fsmcdpk") :interfaces ()
 :mean-control (nominal-adjoint :max-solves 1)
 :parameters (
  (uniform :name "ambient" :target natural-convection-ambient :entity "air" :low 294K :high 300K)
  (uniform :name "reservoir" :target radiation-reservoir-temperature :entity "radiator" :low 280K :high 288K)))"#;

fn base() -> ProjectSpec {
    let mut project = crate::parse_sexpr_migrating(include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"), "/../../data/reference-project/cooling-radiation.fsim"
    ))).unwrap().decoded.spec;
    let envelope = project.envelope.as_mut().unwrap();
    envelope.ambient_lo.value = 290.0;
    envelope.ambient_hi.value = 310.0;
    let setup = project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.boundaries[0].condition = ThermalBoundaryCondition::NaturalConvection {
        characteristic_length: QtyAny::new(0.06, dims::LENGTH),
        ambient_temperature: QtyAny::new(297.0, dims::TEMPERATURE),
        correlation: "convection.churchill-chu-vertical-plate".into(),
    };
    let surface = &mut setup.radiation.as_mut().unwrap().surfaces[0];
    surface.name = "radiator".into();
    surface.reservoir_temperature.value = 285.0;
    assert!(project.validate().is_empty());
    project
}

#[test]
fn boundary_temperature_samples_change_only_the_two_explicit_physical_inputs() {
    let original = base();
    let source = crate::print_sexpr(&original).unwrap();
    let study = UncertaintyStudy::parse(STUDY).unwrap();
    assert_eq!(UncertaintyStudy::parse(study.canonical()).unwrap(), study);
    assert_eq!(study.parameters()[0].target.unit(), "K");
    assert_eq!(study.parameters()[1].target.unit(), "K");
    let policy = study.mean_control().unwrap();
    assert_eq!(policy.probe_count(study.parameters()), 1);
    assert_eq!(policy.probe(study.parameters(), 0).unwrap(), [297.0, 284.0]);
    let bound = study.bind(&original).unwrap();
    for values in [[295.0, 287.0], [299.0, 281.0], [295.0, 287.0]] {
        let mut expected = original.clone();
        let setup = expected.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
        let ThermalBoundaryCondition::NaturalConvection { ambient_temperature, .. }
            = &mut setup.boundaries[0].condition else { panic!("natural boundary"); };
        ambient_temperature.value = values[0];
        setup.radiation.as_mut().unwrap().surfaces[0].reservoir_temperature.value = values[1];
        let sample = bound.sample_project(&values).unwrap();
        assert_eq!(sample, expected, "no envelope/card/query/ownership mutation or accumulated samples");
        let roundtrip = crate::parse_sexpr(&crate::print_sexpr(&sample).unwrap()).unwrap();
        assert_eq!(roundtrip.spec, sample);
    }
    assert_eq!(crate::print_sexpr(bound.base()).unwrap(), source);
    assert_eq!(crate::print_sexpr(&original).unwrap(), source);
    // Reservoir support below the ambient-fluid envelope is intentional: they
    // are different physical temperatures, not aliases of one control.
    assert!(bound.sample_project(&[297.0, 280.0]).is_ok());
    assert!(bound.sample_project(&[297.0, 279.0]).is_err());
    assert!(bound.sample_project(&[f64::NAN, 285.0]).is_err());
}

#[test]
fn boundary_temperature_inputs_require_units_domains_and_exact_law_ownership() {
    let original = base();
    for (old, new) in [(":low 294K", ":low 294"), (":low 280K", ":low 280W"),
        (":low 280K", ":low 0K"), (":high 288K", ":high 279K")] {
        assert!(UncertaintyStudy::parse(&STUDY.replace(old, new)).is_err());
    }
    for source in [STUDY.replace(":low 294K", ":low 280K"),
        STUDY.replace("\"radiator\"", "\"air\""),
        STUDY.replace("natural-convection-ambient", "convection-temperature")] {
        assert!(UncertaintyStudy::parse(&source).unwrap().bind(&original).is_err());
    }
    let mut no_radiation = original.clone();
    no_radiation.cooling.as_mut().unwrap().conduction.as_mut().unwrap().radiation = None;
    assert!(UncertaintyStudy::parse(STUDY).unwrap().bind(&no_radiation).is_err());
    let duplicate = STUDY.replace(":name \"reservoir\" :target radiation-reservoir-temperature :entity \"radiator\"",
        ":name \"second-ambient\" :target natural-convection-ambient :entity \"air\"");
    assert!(UncertaintyStudy::parse(&duplicate).is_err());
    let fixed = STUDY.replace(":high 300K", ":high 294K").replace(":high 288K", ":high 280K")
        .replace(":max-solves 1", ":max-solves 0");
    let study = UncertaintyStudy::parse(&fixed).unwrap();
    assert_eq!(study.mean_control().unwrap().probe_count(study.parameters()), 0);
    assert!(study.bind(&original).is_ok());
}
