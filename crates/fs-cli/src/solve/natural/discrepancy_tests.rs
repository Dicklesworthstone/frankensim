//! Native natural-law admission and explicit discrepancy vertices.
use super::*;

fn project() -> fs_project::ProjectSpec {
    let mut spec = fs_project::parse_sexpr_migrating(include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"), "/../../data/reference-project/cooling-reference.fsim"
    ))).unwrap().decoded.spec;
    spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries[0].condition =
        ThermalBoundaryCondition::NaturalConvection {
            characteristic_length: fs_qty::QtyAny::new(0.06, fs_project::spec::dims::LENGTH),
            ambient_temperature: fs_qty::QtyAny::new(300.0, fs_project::spec::dims::TEMPERATURE),
            correlation: "convection.churchill-chu-vertical-plate".into(),
        };
    spec
}

#[test]
fn discrepancy_multiplier_travels_with_every_temperature_evaluation() {
    let spec = project();
    let setup = spec.cooling.as_ref().unwrap().conduction.as_ref().unwrap();
    let nominal = natural_laws(setup).unwrap().remove(0);
    for scale in [0.75, 1.0, 1.25] {
        let scaled = natural_laws_scaled(setup, scale).unwrap().remove(0);
        let initial = initial_coefficients(std::slice::from_ref(&scaled), 101_325.0).unwrap();
        assert_eq!(initial[&scaled.target].to_bits(),
            (scale * coefficient(&nominal, INITIAL_DELTA_T_K, 101_325.0).unwrap().htc_w_m2_k).to_bits());
        let mut previous = None;
        for delta in [-40.0, -10.0, -1.0, -0.2, 0.2, 1.0, 10.0, 40.0] {
            let raw = coefficient(&nominal, delta, 101_325.0).unwrap();
            let actual = coefficient(&scaled, delta, 101_325.0).unwrap();
            assert_eq!(actual.htc_w_m2_k.to_bits(), (scale * raw.htc_w_m2_k).to_bits());
            // Scaling h is not a change to the card's dimensionless inputs or
            // an instruction to freeze the original/initial wall temperature.
            assert_eq!(actual.rayleigh.to_bits(), raw.rayleigh.to_bits());
            assert_eq!(actual.nusselt.to_bits(), raw.nusselt.to_bits());
            if let Some(previous) = previous { assert_ne!(actual.htc_w_m2_k, previous); }
            previous = Some(actual.htc_w_m2_k);
        }
    }
}

#[test]
fn invalid_discrepancy_vertices_refuse_without_clipping_or_domain_extension() {
    let spec = project();
    let setup = spec.cooling.as_ref().unwrap().conduction.as_ref().unwrap();
    for scale in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(natural_laws_scaled(setup, scale).unwrap_err().code,
            "cli-solve-conduction-natural-scale");
    }
    let huge = natural_laws_scaled(setup, f64::MAX).unwrap().remove(0);
    assert!(coefficient(&huge, 40.0, 101_325.0).is_err(), "overflow is not a finite coefficient");
    let scaled = natural_laws_scaled(setup, 0.75).unwrap().remove(0);
    for delta in [0.0, -0.0, -300.0, f64::NAN] {
        assert!(coefficient(&scaled, delta, 101_325.0).is_err());
    }
    let mut outside = scaled;
    outside.length_m = 1.0e9;
    assert!(coefficient(&outside, 10.0, 101_325.0).is_err(), "a multiplier cannot extend a Rayleigh domain");
}

#[test]
fn scaled_receipt_distinguishes_applied_htc_from_the_unmodified_card() {
    let spec = project();
    let setup = spec.cooling.as_ref().unwrap().conduction.as_ref().unwrap();
    for scale in [1.0, 0.75] {
        let law = natural_laws_scaled(setup, scale).unwrap().remove(0);
        let row = Converged {
            coefficient: coefficient(&law, 10.0, 101_325.0).unwrap(),
            law, mean_wall_k: 310.0, heat_rate_w: 4.0,
        };
        let text = receipt_fragment(&[row], 3).unwrap();
        let value = crate::json_read::JsonValue::parse(&text).unwrap();
        let law = &value.get("laws").unwrap().as_array().unwrap()[0];
        if scale == 1.0 {
            assert!(law.get("htc_multiplier").is_none(), "preserve ordinary nominal receipt bytes");
        } else {
            assert_eq!(law.f64_field("htc_multiplier"), Some(scale));
        }
        let expected = scale * law.f64_field("nusselt").unwrap()
            * AIR_THERMAL_CONDUCTIVITY_W_M_K / law.f64_field("characteristic_length_m").unwrap();
        let actual = law.f64_field("htc_w_m2_k").unwrap();
        assert!((actual - expected).abs() <= 4.0 * f64::EPSILON * actual.abs());
    }
}

#[test]
fn adaptive_admission_cannot_freeze_the_law() {
    let mut spec = project();
    let setup = spec.cooling.as_ref().unwrap().conduction.as_ref().unwrap();
    assert_eq!(super::super::coolest_declared_temperature(setup), Some(300.0));
    assert!(admit_fidelity(&spec, setup).is_ok());
    spec.solver.as_mut().unwrap().fidelity = "ladder".into();
    assert!(admit_fidelity(&spec, spec.cooling.as_ref().unwrap().conduction.as_ref().unwrap()).is_ok());
    spec.solver.as_mut().unwrap().fidelity = "adaptive".into();
    let setup = spec.cooling.as_ref().unwrap().conduction.as_ref().unwrap();
    assert!(admit_fidelity(&spec, setup).is_ok(),
        "the admitted vertical-plate law now has a complete enriched-goal tangent");
    let setup = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    let mut cold = setup.boundaries[0].clone();
    cold.target = "cold-support".into();
    cold.condition = ThermalBoundaryCondition::FixedTemperature {
        temperature: fs_qty::QtyAny::new(280.0, fs_project::spec::dims::TEMPERATURE),
    };
    setup.boundaries.push(cold);
    assert_eq!(super::super::coolest_declared_temperature(setup), Some(280.0));
    // Unrelated affine adaptive projects keep their existing admission.
    setup.boundaries.remove(0);
    assert!(admit_fidelity(&spec, spec.cooling.as_ref().unwrap().conduction.as_ref().unwrap()).is_ok());
}
