//! Both control families reuse the same physical state and complete dual.
use super::*;

const COMBINED: &str = "temperature-max-contact-boundary-adjoint";
const CONTACT: &str = "contact-resistance-multiplier";
const FIXED: &str = "fixed-temperature";

#[test]
fn combined_controls_equal_the_separate_reports_and_joint_physical_resolves() {
    for nonlinear in [false, true] { for cooling in 0..3 {
        let mut f = Fixture::new(nonlinear, cooling != 0);
        let (target, wall) = if cooling == 0 { ("cold", 293.15) } else {
            let setup = f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
            setup.boundaries[1].condition = B::FixedTemperature {
                temperature: QtyAny::new(330.0, dims::TEMPERATURE),
            };
            if cooling == 2 {
                setup.boundaries[0].condition = B::NaturalConvection {
                    characteristic_length: QtyAny::new(0.06, dims::LENGTH),
                    ambient_temperature: QtyAny::new(297.0, dims::TEMPERATURE),
                    correlation: "convection.churchill-chu-vertical-plate".into(),
                };
            }
            f.project.requirements.as_mut().unwrap()[0].region = "cold".into();
            ("hot", 330.0)
        };
        let (plain, original, _) = f.solve(0.13, None, 100);
        let (combined, field, combined_run) = f.solve(0.13, Some(COMBINED), 101);
        let report = combined.get("nominal_adjoint").unwrap();
        assert_eq!(report.str_field("output"), Some(COMBINED));
        assert_eq!(report.str_field("authority"), Some("Estimated"));
        assert_eq!(original.get("temperature"), field.get("temperature"));
        for key in ["energy", "interfaces", "conjugate", "radiation", "natural_convection"] {
            assert_eq!(plain.get(key), combined.get(key), "{key}");
        }
        let parameters = report.get("parameters").unwrap().as_array().unwrap();
        for (i, (output, omit)) in [(OUTPUT, FIXED), (BOUNDARY_OUTPUT, CONTACT)].into_iter().enumerate() {
            let (single, single_field, run) = f.solve(0.13, Some(output), 102 + i);
            assert_ne!(combined_run, run, "the new explicit output has a distinct cache identity");
            assert_eq!(single_field.get("temperature"), field.get("temperature"));
            let single = single.get("nominal_adjoint").unwrap();
            let shared: Vec<_> = parameters.iter().filter(|r| r.str_field("target") != Some(omit)).collect();
            assert_eq!(shared, single.get("parameters").unwrap().as_array().unwrap().iter().collect::<Vec<_>>());
            for key in ["selected_vertex", "value_k", "mode", "dual_iterations", "true_relative_residual",
                "stability_iterations", "response_iterations"] {
                assert_eq!(report.get(key), single.get(key), "no second dual or changed numerical work: {key}");
            }
        }
        assert!(report.get("unsupported").unwrap().as_array().unwrap().iter()
            .all(|r| !matches!(r.str_field("target"), Some(FIXED | CONTACT))));
        let derivative = |kind, entity| {
            let rows: Vec<_> = parameters.iter().filter(|r|
                r.str_field("target") == Some(kind) && r.str_field("entity") == Some(entity)).collect();
            assert_eq!(rows.len(), 1);
            rows[0].f64_field("derivative").unwrap()
        };
        let contact = derivative(CONTACT, "cold-hot-joint");
        let fixed = derivative(FIXED, target);
        assert!(contact.abs() > 1e-6 && fixed > 0.0, "both controls must act on the observed region");
        let selected = report.f64_field("selected_vertex").unwrap() as usize;
        // Independent and simultaneous perturbations: the mixed case detects
        // dropping either coordinate or replacing the complete feedback dual.
        for (case, (dr, dt)) in [(1.0, 0.0), (0.0, 1.0), (0.25, 0.5)].into_iter().enumerate() {
            let step = 0.002_f64;
            let mut values = Vec::new();
            for (side, sign) in [-1.0, 1.0].into_iter().enumerate() {
                fixed_value(&mut f, target, wall + sign * step * dt);
                let (_, perturbed, _) = f.solve(0.13 * (sign * step * dr).exp(), None, 104 + case * 2 + side);
                values.push(perturbed.get("temperature").unwrap().as_array().unwrap()[selected].as_f64().unwrap());
            }
            let expected = (values[1] - values[0]) / (2.0 * step);
            let actual = dr * contact + dt * fixed;
            assert!((actual - expected).abs() < 5e-4 * expected.abs().max(0.01),
                "nonlinear={nonlinear} cooling={cooling} direction=({dr},{dt}): {actual:e} != {expected:e}");
        }
    } }
}
