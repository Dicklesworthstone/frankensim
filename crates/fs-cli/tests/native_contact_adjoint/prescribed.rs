//! Reuse the actual native contact-pair import/card/solve fixture.
use super::*;
use fs_project::{ThermalBoundaryCondition as B, spec::dims};
use fs_qty::QtyAny;

const BOUNDARY_OUTPUT: &str = "temperature-max-boundary-adjoint";

fn fixed_value(f: &mut Fixture, target: &str, value: f64) {
    let setup = f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    let boundary = setup.boundaries.iter_mut().find(|b| b.target == target).unwrap();
    let B::FixedTemperature {temperature} = &mut boundary.condition else {panic!("fixed boundary");};
    temperature.value = value;
}

#[test]
fn native_prescribed_controls_include_contact_material_and_cooling_feedback() {
    for nonlinear in [false,true] { for cooling in 0..3 {
        let mut f = Fixture::new(nonlinear,cooling != 0);
        let (target,value) = if cooling == 0 {("cold",293.15)} else {
            let setup = f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
            setup.boundaries[1].condition = B::FixedTemperature {temperature:QtyAny::new(330.0,dims::TEMPERATURE)};
            if cooling == 2 {
                setup.boundaries[0].condition = B::NaturalConvection {
                    characteristic_length:QtyAny::new(0.06,dims::LENGTH),
                    ambient_temperature:QtyAny::new(297.0,dims::TEMPERATURE),
                    correlation:"convection.churchill-chu-vertical-plate".into(),
                };
            }
            f.project.requirements.as_mut().unwrap()[0].region = "cold".into();
            ("hot",330.0)
        };
        let (plain,original,_) = f.solve(0.13,None,0);
        let (extended,field,run) = f.solve(0.13,Some(BOUNDARY_OUTPUT),1);
        let (legacy,legacy_field,legacy_run) = f.solve(0.13,Some("temperature-max-adjoint"),2);
        assert_ne!(run,legacy_run,"opt-in output has its own native project/run identity");
        assert_eq!(original.get("temperature"),field.get("temperature"));
        assert_eq!(legacy_field.get("temperature"),field.get("temperature"));
        for key in ["energy","interfaces","conjugate","radiation","natural_convection"] {
            assert_eq!(plain.get(key),extended.get(key),"{key}");
        }
        let report = extended.get("nominal_adjoint").unwrap();
        let legacy = legacy.get("nominal_adjoint").unwrap();
        assert_eq!(report.str_field("output"),Some(BOUNDARY_OUTPUT));
        let parameters = report.get("parameters").unwrap().as_array().unwrap();
        let fixed:Vec<_> = parameters.iter().filter(|r| r.str_field("target")==Some("fixed-temperature")).collect();
        assert_eq!(fixed.len(),1);
        assert_eq!(fixed[0].str_field("entity"),Some(target));
        assert_eq!(fixed[0].str_field("parameter_unit"),Some("K"));
        let actual = fixed[0].f64_field("derivative").unwrap();
        let selected = report.f64_field("selected_vertex").unwrap() as usize;
        let unchanged:Vec<_> = parameters.iter().filter(|r| r.str_field("target")!=Some("fixed-temperature")).collect();
        assert_eq!(unchanged,legacy.get("parameters").unwrap().as_array().unwrap().iter().collect::<Vec<_>>());
        for key in ["mode","selected_vertex","value_k","dual_iterations","true_relative_residual"] {
            assert_eq!(report.get(key),legacy.get(key),"the original complete adjoint must be reused");
        }
        assert!(report.get("unsupported").unwrap().as_array().unwrap().iter()
            .all(|r| r.str_field("target")!=Some("fixed-temperature")));
        let step = 0.002;
        let mut values = Vec::new();
        for (side,sign) in [-1.0,1.0].into_iter().enumerate() {
            fixed_value(&mut f,target,value+sign*step);
            let (_,field,_) = f.solve(0.13,None,3+side);
            values.push(field.get("temperature").unwrap().as_array().unwrap()[selected].as_f64().unwrap());
        }
        let expected = (values[1]-values[0])/(2.0*step);
        assert!(actual>0.0 && expected>0.0);
        assert!((actual-expected).abs()<3e-4*expected.abs().max(0.01),
            "nonlinear={nonlinear} cooling={cooling}: {actual:e} != {expected:e}");
    } }
}

#[path = "uncertainty.rs"]
mod uncertainty;

#[test]
fn native_cooled_walls_keep_inward_heat_and_complete_boundary_sensitivities() {
    for nonlinear in [false, true] {
        let mut f = Fixture::new(nonlinear, false);
        f.project.power.as_mut().unwrap()[0].watts.value = 0.1;
        let setup = f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
        setup.boundaries[0].condition = B::NaturalConvection {
            characteristic_length: QtyAny::new(0.06, dims::LENGTH),
            ambient_temperature: QtyAny::new(310.0, dims::TEMPERATURE),
            correlation: "convection.churchill-chu-vertical-plate".into(),
        };
        setup.boundaries[1].condition = B::FixedTemperature {
            temperature: QtyAny::new(285.0, dims::TEMPERATURE),
        };
        f.project.requirements.as_mut().unwrap()[0].region = "cold".into();
        let (plain, original, _) = f.solve(0.13, None, 20);
        let (nominal, field, _) = f.solve(0.13, Some(BOUNDARY_OUTPUT), 21);
        assert_eq!(original.get("temperature"), field.get("temperature"));
        for key in ["energy", "interfaces", "natural_convection"] {
            assert_eq!(plain.get(key), nominal.get(key), "{key}");
        }
        let law = &nominal.get("natural_convection").unwrap()
            .get("laws").unwrap().as_array().unwrap()[0];
        assert!(law.f64_field("mean_wall_k").unwrap() < 310.0);
        assert!(law.f64_field("delta_t_k").unwrap() < 0.0);
        assert!(law.f64_field("heat_rate_w").unwrap() < 0.0,
            "a colder wall absorbs heat from ambient; do not take |flux|");
        assert!(law.f64_field("htc_w_m2_k").unwrap() > 0.0);
        assert_eq!(nominal.get("solver_control").unwrap().str_field("status"),
            Some("unsupported-model"));
        let report = nominal.get("nominal_adjoint").unwrap();
        let selected = report.f64_field("selected_vertex").unwrap() as usize;
        let parameters = report.get("parameters").unwrap().as_array().unwrap();
        let step = 0.002;
        for (case, (target, entity)) in [("fixed-temperature", "hot"),
            ("natural-convection-ambient", "cold")].into_iter().enumerate()
        {
            let actual = parameters.iter().find(|r| r.str_field("target") == Some(target)
                && r.str_field("entity") == Some(entity)).unwrap()
                .f64_field("derivative").unwrap();
            let mut values = Vec::new();
            for (side, sign) in [-1.0, 1.0].into_iter().enumerate() {
                fixed_value(&mut f, "hot", 285.0 + if case == 0 { sign * step } else { 0.0 });
                let setup = f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
                let B::NaturalConvection { ambient_temperature, .. } = &mut setup.boundaries[0].condition
                    else { panic!("cooled-wall law"); };
                ambient_temperature.value = 310.0 + if case == 1 { sign * step } else { 0.0 };
                let (_, field, _) = f.solve(0.13, None, 22 + case * 2 + side);
                values.push(field.get("temperature").unwrap().as_array().unwrap()[selected].as_f64().unwrap());
            }
            let expected = (values[1] - values[0]) / (2.0 * step);
            assert!(actual > 0.0 && expected > 0.0);
            assert!((actual - expected).abs() < 3e-4 * expected.abs().max(0.01),
                "nonlinear={nonlinear} {target}: {actual:e} != {expected:e}");
        }
    }
}

#[path = "combined.rs"]
mod combined;
