//! A prescribed solid temperature is not the fluid ambient or a material edit.
use super::*;

const STUDY: &str = r#"(fsim-uncertainty-study
 :version 3 :project "contact-pair.fsim" :samples 4 :seed 59
 :wall-time 120s :method monte-carlo :correlation independent :qoi "temperature-max"
 :geometry ((mesh :role "cold-body" :path "cold-body.stl" :unit "m" :max-hole-edges 0)
            (mesh :role "hot-body" :path "hot-body.stl" :unit "m" :max-hole-edges 0))
 :materials ("solid.fsmcdpk") :interfaces ("contact.fsintpk")
 :mean-control (nominal-adjoint :max-solves 1)
 :parameters (
  (uniform :name "wall" :target fixed-temperature :entity "cold" :low 280K :high 290K)
  (uniform :name "power" :target power :entity "hot" :low 4W :high 6W)))"#;

fn base() -> ProjectSpec {
    crate::parse_sexpr_migrating(include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
        "/../../examples/contact-pair/contact-pair.fsim"))).unwrap().decoded.spec
}

#[test]
fn fixed_temperature_sampling_preserves_the_base_and_every_unrelated_field() {
    let original = base();
    let study = UncertaintyStudy::parse(STUDY).unwrap();
    assert_eq!(UncertaintyStudy::parse(study.canonical()).unwrap(), study);
    assert_eq!(study.parameters()[0].target, Target::FixedTemperature);
    assert_eq!(study.parameters()[0].target.unit(), "K");
    let policy = study.mean_control().unwrap();
    assert_eq!(policy.probe_count(study.parameters()), 1);
    assert_eq!(policy.probe(study.parameters(), 0).unwrap(), [285.0, 5.0]);
    let bound = study.bind(&original).unwrap();
    for values in [[281.0,4.5], [289.0,5.5], [281.0,4.5]] {
        let mut expected = original.clone();
        let setup = expected.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
        let ThermalBoundaryCondition::FixedTemperature {temperature} = &mut setup.boundaries[0].condition
            else { panic!("contact-pair fixed boundary"); };
        temperature.value = values[0];
        expected.power.as_mut().unwrap()[0].watts.value = values[1];
        let actual = bound.sample_project(&values).unwrap();
        assert_eq!(actual, expected, "no ambient, material, contact, duty, output or geometry mutation");
        assert_eq!(crate::parse_sexpr(&crate::print_sexpr(&actual).unwrap()).unwrap().spec, actual);
    }
    assert_eq!(bound.base(), &original);
    // The wall support is below the ambient-fluid envelope, by declaration.
    assert!(bound.sample_project(&[280.0,5.0]).is_ok());
    assert!(bound.sample_project(&[279.0,5.0]).is_err());
    assert!(bound.sample_project(&[f64::NAN,5.0]).is_err());
    assert!(bound.sample_project(&[285.0]).is_err());
}

#[test]
fn fixed_temperature_requires_explicit_units_positive_support_and_exact_boundary_law() {
    for (old,new) in [(":low 280K",":low 280"), (":low 280K",":low 280W"),
        (":low 280K",":low 0K"), (":high 290K",":high 279K")] {
        assert!(UncertaintyStudy::parse(&STUDY.replace(old,new)).is_err());
    }
    for source in [STUDY.replace("fixed-temperature", "convection-temperature"),
        STUDY.replace(":entity \"cold\"", ":entity \"hot\""),
        STUDY.replace(":entity \"cold\"", ":entity \"missing\"")] {
        assert!(UncertaintyStudy::parse(&source).unwrap().bind(&base()).is_err());
    }
    let mut changed = base();
    changed.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries[0].condition =
        ThermalBoundaryCondition::Convection {
            coefficient: QtyAny::new(10.0,crate::spec::dims::HEAT_TRANSFER_COEFFICIENT),
            reference_temperature: QtyAny::new(293.15,crate::spec::dims::TEMPERATURE),
        };
    assert!(UncertaintyStudy::parse(STUDY).unwrap().bind(&changed).is_err());
    let deterministic = STUDY.replace(":high 290K",":high 280K")
        .replace(":high 6W",":high 4W").replace(":max-solves 1",":max-solves 0");
    let study = UncertaintyStudy::parse(&deterministic).unwrap();
    assert_eq!(study.mean_control().unwrap().probe_count(study.parameters()),0);
    assert!(study.bind(&base()).is_ok());
    let secant = STUDY.replace("nominal-adjoint :max-solves 1", "coordinate-secant :max-solves 4");
    assert!(UncertaintyStudy::parse(&secant).unwrap().bind(&base()).is_ok());
}
