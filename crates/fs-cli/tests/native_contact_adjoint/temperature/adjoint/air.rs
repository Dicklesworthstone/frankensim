//! G3: source contact, stream-wise air and natural walls share one total dual.
use super::*;

fn fixture(mixed: bool) -> Fixture {
    let mut f = Fixture::new(true, true);
    let setup = f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.radiation = None;
    if mixed {
        setup.boundaries[1].condition = B::NaturalConvection {
            characteristic_length: QtyAny::new(1.0, dims::LENGTH),
            ambient_temperature: QtyAny::new(297.0, dims::TEMPERATURE),
            correlation: "convection.churchill-chu-vertical-plate".into(),
        };
    } else {
        let mut downstream = setup.boundaries[0].condition.clone();
        let B::AirflowConvection { order, .. } = &mut downstream else { panic!("air law"); };
        *order = 1;
        setup.boundaries[1].condition = downstream;
    }
    // Keep density and the hydraulic operating point fixed while perturbing
    // the independent inlet or natural reservoir. Every trial stays admitted.
    let envelope = f.project.envelope.as_mut().unwrap();
    envelope.ambient_lo.value = 296.0;
    envelope.ambient_hi.value = 298.0;
    f
}

#[test]
fn live_contact_air_and_mixed_controls_match_complete_physical_resolves() {
    for mixed in [false, true] {
        let mut f = fixture(mixed);
        let surface = mixed;
        let (nominal, field) = solve(&f, &[1.0], surface, Some(OUTPUT), 120);
        let report = nominal.get("nominal_adjoint").unwrap();
        assert_eq!(report.str_field("mode"), Some(if mixed {
            "nonlinear-contact-natural-full-air-feedback"
        } else { "nonlinear-contact-full-air-feedback" }));
        assert!(report.f64_field("true_relative_residual").unwrap() <= 1e-12);
        let selected = report.f64_field("selected_vertex").unwrap() as usize;
        assert_eq!(report.f64_field("value_k"), Some(value(&field, selected)));
        let inlet: Vec<_> = report.get("parameters").unwrap().as_array().unwrap().iter()
            .filter(|row| row.str_field("target") == Some("air-inlet-temperature")).collect();
        assert_eq!(inlet.len(), 1, "all ordered segments share one physical branch inlet");
        assert_eq!(inlet[0].str_field("entity"), Some("cold"));
        assert!(inlet[0].f64_field("derivative").unwrap() > 0.0);
        assert!(derivative(report, "fan-speed-ratio", "fixture-bank").abs() > 1e-6);
        assert!(derivative(report, TARGET, "cold-hot-joint").abs() > 1e-6);
        let contact = nominal.get("temperature_contact").unwrap();
        assert!(contact.f64_field("residual_norm_w").unwrap() <= contact.f64_field("residual_tolerance_w").unwrap());
        assert!(contact.f64_field("energy_closure_w").unwrap().abs() <= contact.f64_field("energy_tolerance_w").unwrap());
        let air = nominal.get("conjugate").unwrap();
        assert!((air.f64_field("solid_total_w").unwrap() - air.f64_field("air_total_w").unwrap()).abs()
            <= air.f64_field("balance_tolerance_w").unwrap());
        if mixed {
            let natural: f64 = nominal.get("natural_convection").unwrap().get("laws").unwrap()
                .as_array().unwrap().iter().map(|row| row.f64_field("heat_rate_w").unwrap()).sum();
            assert!(natural > 0.01);
            assert!((air.f64_field("solid_total_w").unwrap() + natural - 4.44).abs() < 1e-7);
        } else {
            assert_eq!(air.get("segments").unwrap().as_array().unwrap().len(), 2);
        }
        let (legacy, legacy_field) = solve(&f, &[1.0], surface, Some("temperature-max-adjoint"), 121);
        assert_eq!(legacy_field.get("temperature"), field.get("temperature"));
        let legacy = legacy.get("nominal_adjoint").unwrap();
        let shared: Vec<_> = report.get("parameters").unwrap().as_array().unwrap().iter()
            .filter(|row| row.str_field("target") != Some(TARGET)).collect();
        assert_eq!(shared, legacy.get("parameters").unwrap().as_array().unwrap().iter().collect::<Vec<_>>());
        for key in ["mode", "selected_vertex", "value_k", "dual_iterations", "true_relative_residual"] {
            assert_eq!(report.get(key), legacy.get(key));
        }
        for (case, (target, entity, step)) in [(TARGET, "cold-hot-joint", 0.002),
            ("air-inlet-temperature", "cold", 0.002), ("power", "hot", 0.01),
            ("fan-speed-ratio", "fixture-bank", 0.002), ("natural-convection-ambient", "hot", 0.002)]
            .into_iter().enumerate()
        {
            if !mixed && target == "natural-convection-ambient" { continue; }
            let actual = derivative(report, target, entity);
            let mut values = Vec::new();
            for (side, sign) in [-1.0, 1.0].into_iter().enumerate() {
                f.project.power.as_mut().unwrap()[0].watts.value = 12.0
                    + if target == "power" { sign * step } else { 0.0 };
                let cooling = f.project.cooling.as_mut().unwrap();
                cooling.fan_system.as_mut().unwrap().banks[0].speed_ratio = 1.0
                    + if target == "fan-speed-ratio" { sign * step } else { 0.0 };
                for boundary in &mut cooling.conduction.as_mut().unwrap().boundaries {
                    match &mut boundary.condition {
                        B::AirflowConvection { inlet_temperature, .. } => {
                            inlet_temperature.value = 297.0
                                + if target == "air-inlet-temperature" { sign * step } else { 0.0 };
                        }
                        B::NaturalConvection { ambient_temperature, .. } => {
                            ambient_temperature.value = 297.0
                                + if target == "natural-convection-ambient" { sign * step } else { 0.0 };
                        }
                        _ => {}
                    }
                }
                let scale = 1.0 + if target == TARGET { sign * step } else { 0.0 };
                let (_, changed) = solve(&f, &[scale], surface, None, 122 + 2 * case + side);
                values.push(value(&changed, selected));
            }
            agrees(actual, (values[1] - values[0]) / (2.0 * step));
        }
    }
}
