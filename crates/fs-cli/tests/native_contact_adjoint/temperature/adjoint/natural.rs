//! G3: source-contact and signed natural-convection derivatives in real projects.
use super::*;

fn fixture(cooled: bool) -> Fixture {
    let mut f = Fixture::new(true, false);
    let ambient = if cooled { 310.0 } else { 293.15 };
    let setup = f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.boundaries[0].condition = B::NaturalConvection {
        characteristic_length: QtyAny::new(1.0, dims::LENGTH),
        ambient_temperature: QtyAny::new(ambient, dims::TEMPERATURE),
        correlation: "convection.churchill-chu-vertical-plate".into(),
    };
    if cooled {
        setup.boundaries[1].condition = B::FixedTemperature {
            temperature: QtyAny::new(285.0, dims::TEMPERATURE),
        };
        f.project.power.as_mut().unwrap()[0].watts.value = 0.1;
        f.project.requirements.as_mut().unwrap()[0].region = "cold".into();
    }
    let envelope = f.project.envelope.as_mut().unwrap();
    envelope.ambient_lo.value = ambient;
    envelope.ambient_hi.value = ambient;
    f
}

#[test]
fn heated_and_cooled_walls_keep_contact_and_natural_feedback_in_one_adjoint() {
    // The heated project uses the Newton curve producer; the cooled one uses
    // the original source-surface Picard producer. Both retain nonlinear K(T).
    for (cooled, surface) in [(false, false), (true, true)] {
        let mut f = fixture(cooled);
        let power = f.project.power.as_ref().unwrap()[0].watts.value;
        let ambient = if cooled { 310.0 } else { 293.15 };
        let (nominal, field) = solve(&f, &[1.0], surface, Some(OUTPUT), 100);
        let report = nominal.get("nominal_adjoint").unwrap();
        assert_eq!(report.str_field("mode"), Some("nonlinear-contact-natural-convection-feedback"));
        assert_eq!(report.str_field("authority"), Some("Estimated"));
        assert!(report.f64_field("true_relative_residual").unwrap() <= 1e-12);
        let selected = report.f64_field("selected_vertex").unwrap() as usize;
        assert_eq!(report.f64_field("value_k"), Some(value(&field, selected)));
        let natural = &nominal.get("natural_convection").unwrap().get("laws").unwrap().as_array().unwrap()[0];
        assert_eq!(natural.f64_field("delta_t_k").unwrap() < 0.0, cooled);
        assert_eq!(natural.f64_field("heat_rate_w").unwrap() < 0.0, cooled);
        assert!(derivative(report, TARGET, "cold-hot-joint") > 1e-4);
        assert!(derivative(report, "natural-convection-ambient", "cold") > 0.0);
        let contact = nominal.get("temperature_contact").unwrap();
        assert!(contact.f64_field("residual_norm_w").unwrap() <= contact.f64_field("residual_tolerance_w").unwrap());
        assert!(contact.f64_field("energy_closure_w").unwrap().abs() <= contact.f64_field("energy_tolerance_w").unwrap());
        // Both report forms use exactly the same accepted field and single
        // complete transpose. Extra resistance rows add only contractions.
        let (legacy, legacy_field) = solve(&f, &[1.0], surface, Some("temperature-max-adjoint"), 101);
        assert_eq!(legacy_field.get("temperature"), field.get("temperature"));
        let legacy = legacy.get("nominal_adjoint").unwrap();
        let shared: Vec<_> = report.get("parameters").unwrap().as_array().unwrap().iter()
            .filter(|row| row.str_field("target") != Some(TARGET)).collect();
        assert_eq!(shared, legacy.get("parameters").unwrap().as_array().unwrap().iter().collect::<Vec<_>>());
        for key in ["mode", "selected_vertex", "value_k", "dual_iterations", "true_relative_residual"] {
            assert_eq!(report.get(key), legacy.get(key));
        }
        for (case, (target, entity, step)) in [(TARGET, "cold-hot-joint", 0.002),
            ("natural-convection-ambient", "cold", 0.002), ("power", "hot", 0.01)]
            .into_iter().enumerate()
        {
            // The cooled hot-body boundary fixes the entire contact trace,
            // so its interior source cannot influence the observed cold body.
            if cooled && target == "power" { continue; }
            let actual = derivative(report, target, entity);
            let mut values = Vec::new();
            for (side, sign) in [-1.0, 1.0].into_iter().enumerate() {
                f.project.power.as_mut().unwrap()[0].watts.value = power
                    + if target == "power" { sign * step } else { 0.0 };
                let B::NaturalConvection { ambient_temperature, .. } = &mut f.project.cooling.as_mut()
                    .unwrap().conduction.as_mut().unwrap().boundaries[0].condition else { panic!("natural law"); };
                let reservoir = ambient
                    + if target == "natural-convection-ambient" { sign * step } else { 0.0 };
                ambient_temperature.value = reservoir;
                let envelope = f.project.envelope.as_mut().unwrap();
                envelope.ambient_lo.value = reservoir;
                envelope.ambient_hi.value = reservoir;
                let scale = 1.0 + if target == TARGET { sign * step } else { 0.0 };
                let (_, changed) = solve(&f, &[scale], surface, None, 102 + case * 2 + side);
                values.push(value(&changed, selected));
            }
            agrees(actual, (values[1] - values[0]) / (2.0 * step));
        }
    }
}
