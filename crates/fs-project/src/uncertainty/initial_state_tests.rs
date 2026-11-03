//! Startup uncertainty changes the original stored state, not a thermal bath.
use super::*;
use crate::{TransientInitialTemperature, spec::dims};

fn base(enthalpy: bool) -> ProjectSpec {
    let source = if enthalpy {
        include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../data/reference-project/cooling-enthalpy.fsim"))
    } else {
        include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../data/reference-project/cooling-transient.fsim"))
    };
    let mut project = crate::parse_sexpr_migrating(source).unwrap().decoded.spec;
    if enthalpy {
        let time = project.cooling.as_mut().unwrap().conduction.as_mut().unwrap().transient.as_mut().unwrap();
        time.initial_temperature = TransientInitialTemperature::FromRegionalEnthalpy;
        let material = &mut time.enthalpy.as_mut().unwrap().materials[0];
        // An arbitrary energy origin must not prohibit negative valid h.
        material.initial_specific_enthalpy.value -= 100.0;
        for knot in &mut material.knots { knot.specific_enthalpy.value -= 100.0; }
    }
    assert!(project.validate().is_empty());
    project
}

fn study(enthalpy: bool) -> String {
    let initial = if enthalpy {
        "initial-specific-enthalpy :entity \"solid\" :low -90J/kg :high -30J/kg"
    } else {
        "initial-temperature :entity \"initial-state\" :low 320K :high 340K"
    };
    format!(r#"(fsim-uncertainty-study
        :version 3 :mean-control (nominal-adjoint :max-solves 1)
        :project "cooling.fsim" :samples 4 :seed 29 :wall-time 120s
        :method monte-carlo :correlation independent :qoi "temperature-peak"
        :geometry ((mesh :role "enclosure" :path "plate.stl" :unit "m" :max-hole-edges 0))
        :materials ("aa6061.fsmcdpk") :interfaces ()
        :parameters (
            (uniform :name "power" :target power :entity "solid" :low 5W :high 5W)
            (uniform :name "startup" :target {initial})))"#)
}

#[test]
fn initial_state_samples_preserve_storage_sources_and_do_not_accumulate() {
    for enthalpy in [false, true] {
        let original = base(enthalpy);
        let parsed = UncertaintyStudy::parse(&study(enthalpy)).unwrap();
        assert_eq!(UncertaintyStudy::parse(parsed.canonical()).unwrap(), parsed);
        assert_eq!(parsed.parameters()[1].target.unit(), if enthalpy { "J/kg" } else { "K" });
        let bound = parsed.bind(&original).unwrap();
        let values = if enthalpy { [-85.0, -40.0, -85.0] } else { [325.0, 335.0, 325.0] };
        for value in values {
            let sample = bound.sample_project(&[5.0, value]).unwrap();
            let mut expected = original.clone();
            let time = expected.cooling.as_mut().unwrap().conduction.as_mut().unwrap().transient.as_mut().unwrap();
            if enthalpy {
                time.enthalpy.as_mut().unwrap().materials[0].initial_specific_enthalpy.value = value;
            } else {
                time.initial_temperature = TransientInitialTemperature::Uniform(QtyAny::new(value, dims::TEMPERATURE));
            }
            assert_eq!(sample, expected, "charts, density, source cards, boundaries and time grid stay fixed");
            assert_eq!(crate::parse_sexpr(&crate::print_sexpr(&sample).unwrap()).unwrap().spec, sample);
        }
        assert_eq!(bound.base(), &original);
        assert!(bound.sample_project(&[5.0, f64::NAN]).is_err());
    }
}

#[test]
fn initial_state_binding_refuses_wrong_owners_units_and_original_chart_domain() {
    let capacity = base(false);
    let enthalpy = base(true);
    for (source, project) in [
        (study(false).replace("initial-state", "solid"), &capacity),
        (study(true).replace("initial-specific-enthalpy :entity \"solid\"", "initial-specific-enthalpy :entity \"initial-state\""), &enthalpy),
        (study(false), &enthalpy), (study(true), &capacity),
        (study(true).replace("-90J/kg", "-110J/kg"), &enthalpy),
        (study(true).replace("-30J/kg", "90J/kg"), &enthalpy),
    ] {
        assert!(UncertaintyStudy::parse(&source).unwrap().bind(project).is_err());
    }
    for source in [study(false).replace("320K", "0K"),
        study(false).replace("320K", "320J/kg"), study(true).replace("-90J/kg", "290K")]
    {
        assert!(UncertaintyStudy::parse(&source).is_err());
    }
    let mut uniform_enthalpy = enthalpy.clone();
    uniform_enthalpy.cooling.as_mut().unwrap().conduction.as_mut().unwrap().transient.as_mut().unwrap()
        .initial_temperature = TransientInitialTemperature::Uniform(QtyAny::new(293.15, dims::TEMPERATURE));
    assert!(uniform_enthalpy.validate().is_empty());
    assert!(UncertaintyStudy::parse(&study(true)).unwrap().bind(&uniform_enthalpy).is_err(),
        "sampling must not silently replace an authored uniform-temperature assertion");
    let mut steady = capacity;
    steady.cooling.as_mut().unwrap().conduction.as_mut().unwrap().transient = None;
    assert!(UncertaintyStudy::parse(&study(false).replace("temperature-peak", "temperature-max"))
        .unwrap().bind(&steady).is_err());
}

#[test]
fn operating_cases_retain_signed_regional_initial_enthalpy() {
    let source = r#"(fsim-native-design-study
        :version 1 :project "cooling.fsim" :evaluations 32 :seed 29
        :wall-time 120s :qoi "temperature-max"
        :geometry ((mesh :role "enclosure" :path "plate.stl" :unit "m" :max-hole-edges 0))
        :materials ("aa6061.fsmcdpk") :interfaces ()
        :control (regional-power :entity "solid" :minimum 0W :maximum 20W :tolerance 0.01W)
        :temperature-tolerance 0.01K :assumption nondecreasing-temperature-with-power
        :operating-cases (
            (case :name "cold" :inputs ((initial-specific-enthalpy :entity "solid" :value -90J/kg)))
            (case :name "warm" :inputs ((initial-specific-enthalpy :entity "solid" :value -30J/kg)))))"#;
    let original = base(true);
    let bound = UncertaintyStudy::parse(source).unwrap().bind(&original).unwrap();
    let sample = bound.sample_project(&[10.0, -90.0]).unwrap();
    let material = &sample.cooling.as_ref().unwrap().conduction.as_ref().unwrap()
        .transient.as_ref().unwrap().enthalpy.as_ref().unwrap().materials[0];
    assert_eq!(material.initial_specific_enthalpy.value, -90.0);
    assert!(bound.sample_project(&[10.0, -60.0]).is_err(), "only the authored finite cases are admitted");
}
