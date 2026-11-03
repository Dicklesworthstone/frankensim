//! Authored storage coordinates retain their original engineering declarations.
use super::*;

fn base(enthalpy: bool) -> ProjectSpec {
    let source = if enthalpy {
        include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../data/reference-project/cooling-enthalpy.fsim"
        ))
    } else {
        include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../data/reference-project/cooling-transient.fsim"
        ))
    };
    crate::parse_sexpr_migrating(source).unwrap().decoded.spec
}

fn study(enthalpy: bool) -> String {
    let coordinate = if enthalpy {
        "enthalpy-reference-density :entity \"solid\" :low 5kg/m^3 :high 20kg/m^3"
    } else {
        "thermal-capacity-scale :entity \"thermal-storage\" :low 0.5 :high 2"
    };
    format!(
        r#"(fsim-uncertainty-study
        :version 3 :mean-control (nominal-adjoint :max-solves 1)
        :project "cooling.fsim" :samples 4 :seed 29 :wall-time 120s
        :method monte-carlo :correlation independent :qoi "temperature-max"
        :geometry ((mesh :role "enclosure" :path "plate.stl" :unit "m" :max-hole-edges 0))
        :materials ("aa6061.fsmcdpk") :interfaces ()
        :parameters (
            (uniform :name "power" :target power :entity "solid" :low 5W :high 5W)
            (uniform :name "storage" :target {coordinate})))"#
    )
}

#[test]
fn storage_samples_preserve_sources_charts_and_original_scale() {
    for enthalpy in [false, true] {
        let original = base(enthalpy);
        let parsed = UncertaintyStudy::parse(&study(enthalpy)).unwrap();
        assert_eq!(UncertaintyStudy::parse(parsed.canonical()).unwrap(), parsed);
        assert_eq!(
            parsed.parameters()[1].target.unit(),
            if enthalpy { "kg/m^3" } else { "1" }
        );
        let bound = parsed.bind(&original).unwrap();
        let values = if enthalpy {
            [8.0, 18.0, 8.0]
        } else {
            [0.75, 1.5, 0.75]
        };
        for value in values {
            let sample = bound.sample_project(&[5.0, value]).unwrap();
            let mut expected = original.clone();
            let time = expected
                .cooling
                .as_mut()
                .unwrap()
                .conduction
                .as_mut()
                .unwrap()
                .transient
                .as_mut()
                .unwrap();
            if enthalpy {
                time.enthalpy.as_mut().unwrap().materials[0]
                    .reference_density
                    .value = value;
            } else {
                for row in &mut time.capacities {
                    row.volumetric_heat_capacity.value *= value;
                }
            }
            assert_eq!(
                sample, expected,
                "the h/T chart, initial state, sources, cards, geometry and time grid stay fixed"
            );
            assert_eq!(
                crate::parse_sexpr(&crate::print_sexpr(&sample).unwrap())
                    .unwrap()
                    .spec,
                sample
            );
        }
        assert_eq!(bound.base(), &original);
        assert!(bound.sample_project(&[5.0, f64::NAN]).is_err());
        assert!(bound.sample_project(&[5.0, 0.0]).is_err());
    }
}

#[test]
fn storage_coordinates_refuse_wrong_owners_units_and_nonpositive_support() {
    let capacity = base(false);
    let enthalpy = base(true);
    for (source, project) in [
        (study(false).replace("thermal-storage", "solid"), &capacity),
        (
            study(true).replace("density :entity \"solid\"", "density :entity \"unknown\""),
            &enthalpy,
        ),
        (study(false), &enthalpy),
        (study(true), &capacity),
    ] {
        assert!(
            UncertaintyStudy::parse(&source)
                .unwrap()
                .bind(project)
                .is_err()
        );
    }
    for source in [
        study(false).replace(":low 0.5", ":low 0"),
        study(false).replace(":low 0.5", ":low 0.5K"),
        study(true).replace("5kg/m^3", "0kg/m^3"),
        study(true).replace("5kg/m^3", "-5kg/m^3"),
        study(true).replace("5kg/m^3", "5J/kg"),
    ] {
        assert!(UncertaintyStudy::parse(&source).is_err());
    }
    let mut steady = capacity.clone();
    steady
        .cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap()
        .transient = None;
    assert!(
        UncertaintyStudy::parse(&study(false))
            .unwrap()
            .bind(&steady)
            .is_err()
    );
    let parameter = UncertaintyStudy::parse(&study(false)).unwrap().parameters()[1].clone();
    let mut overflow = capacity.clone();
    assert!(apply(&mut overflow, &parameter, f64::MAX).is_err());
    assert_eq!(
        overflow, capacity,
        "an unrepresentable storage scale cannot partly mutate the model"
    );
}
