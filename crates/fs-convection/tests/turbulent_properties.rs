//! G1/G3: original tube entrance, explicit wall properties and complete tangents.

use fs_convection::developing::{
    TurbulentCircularDuctSection, TurbulentTubePropertyCorrection as Correction,
    evaluate_turbulent_circular_duct_section,
};
use fs_convection::{CorrelationError, ThermalConductivity};
use fs_qty::{Length, Temperature};

fn section(start: f64, end: f64, pr: f64, correction: Correction) -> TurbulentCircularDuctSection {
    evaluate_turbulent_circular_duct_section(100_000.0, pr, start, end)
        .unwrap().with_property_correction(correction).unwrap()
}

fn air(bulk: f64, wall: f64) -> Correction {
    Correction::Air {
        bulk_temperature: Temperature::new(bulk),
        wall_temperature: Temperature::new(wall),
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "{actual:.17e} != {expected:.17e}");
}

#[test]
fn source_property_corrections_change_actual_mean_and_finite_section_heat() {
    // Independent Decimal(80) evaluation of Eqs. (3)/(4). At Re=1e5,
    // Pr=1 the corrected Darcy equation gives Nu_infinity=220 exactly;
    // L/D=1000 gives mean Nu=222.2 before the property correction.
    for (correction, expected) in [
        (Correction::Liquid { wall_prandtl: 0.5 }, 239.804_514_151_283_72),
        (Correction::Liquid { wall_prandtl: 2.0 }, 205.887_867_352_040_4),
        (air(300.0, 400.0), 195.218_792_930_244_72),
        (air(500.0, 400.0), 245.670_811_534_528_46),
    ] {
        let mean = section(0.0, 1000.0, 1.0, correction);
        close(mean.nusselt(), expected, 3e-13);
        close(mean.downstream().evidence().value, 222.2, 3e-13);
        assert_eq!(mean.property_correction(), Some(correction));
        let h = mean.heat_transfer_coefficient(
            ThermalConductivity::new(0.6), Length::new(0.02),
        ).unwrap();
        close(h.qoi, expected * 30.0, 3e-13);
        assert!(h.model.discrepancy_rel.is_infinite());
        assert!(h.model.cards.contains(&correction.model_card().name));
        assert!(h.model.validity.bounds().contains_key(match correction {
            Correction::Liquid { .. } => "Pr_over_Pr_wall",
            Correction::Air { .. } => "T_bulk_over_T_wall",
        }));
        assert!(h.sensitivity.d_qoi.is_empty());
        let robin = mean.robin_boundary(
            ThermalConductivity::new(0.6), Length::new(0.02), Temperature::new(300.0),
        ).unwrap();
        assert_eq!(robin.coefficient(), &h);
        assert_eq!(robin.boundary_condition(), &fs_conduction::ThermalBc::robin(h.qoi, 300.0).unwrap());

        // cbrt(27)=3, cbrt(125)=5: the section factor is 1+1/49.
        let finite = section(27.0, 125.0, 1.0, correction);
        close(finite.nusselt(), (220.0 * 50.0 / 49.0) * (expected / 222.2), 3e-13);
    }
}

#[test]
fn corrected_subdivision_conserves_heat_and_replacement_never_compounds_factors() {
    for correction in [Correction::Liquid { wall_prandtl: 3.5 }, air(300.0, 450.0)] {
        let whole = section(0.0, 1000.0, 7.0, correction);
        let integrated: f64 = [0.0, 10.0, 27.0, 64.0, 125.0, 500.0, 1000.0]
            .windows(2).map(|pair| {
                section(pair[0], pair[1], 7.0, correction).nusselt() * (pair[1] - pair[0])
            }).sum();
        close(integrated, 1000.0 * whole.nusselt(), 5e-13);
        let unchanged = whole.with_property_correction(correction).unwrap();
        assert_eq!(unchanged.nusselt(), section(0.0, 1000.0, 7.0, correction).nusselt());
        let changed = unchanged.with_property_correction(Correction::Liquid { wall_prandtl: 14.0 }).unwrap();
        assert_eq!(changed.nusselt(), section(0.0, 1000.0, 7.0,
            Correction::Liquid { wall_prandtl: 14.0 }).nusselt());
    }
    let baseline = evaluate_turbulent_circular_duct_section(100_000.0, 7.0, 10.0, 100.0).unwrap();
    let original = baseline.nusselt();
    assert_eq!(baseline.with_property_correction(Correction::Liquid { wall_prandtl: 7.0 })
        .unwrap().nusselt(), original);
}

#[test]
fn liquid_bulk_and_actual_wall_tangents_match_the_same_corrected_heat_law() {
    for start in [0.0, 27.0] {
        let value = |re, pr, wall_pr, k, d| {
            evaluate_turbulent_circular_duct_section(re, pr, start, 125.0).unwrap()
                .with_property_correction(Correction::Liquid { wall_prandtl: wall_pr }).unwrap()
                .heat_transfer_coefficient(ThermalConductivity::new(k), Length::new(d)).unwrap().qoi
        };
        let inputs = [100_000.0, 7.0, 3.5, 0.6, 0.02];
        let evaluation = section(start, 125.0, 7.0, Correction::Liquid { wall_prandtl: 3.5 })
            .heat_transfer_coefficient_with_sensitivities(
                ThermalConductivity::new(0.6), Length::new(0.02),
            ).unwrap();
        for (index, axis) in ["Re", "Pr", "Pr_wall", "k", "L"].into_iter().enumerate() {
            let mut plus = inputs;
            let mut minus = inputs;
            let delta = 1e-5 * inputs[index];
            plus[index] += delta;
            minus[index] -= delta;
            let run = |x: [f64; 5]| value(x[0], x[1], x[2], x[3], x[4]);
            close(evaluation.sensitivity.d_qoi[axis], (run(plus) - run(minus)) / (2.0 * delta), 3e-8);
        }
        // A viscosity-only or held-wall correction misses these coupled terms.
        close(evaluation.sensitivity.d_qoi["Pr_wall"], -0.11 * evaluation.qoi / 3.5, 1e-14);
        let baseline = evaluate_turbulent_circular_duct_section(100_000.0, 7.0, start, 125.0).unwrap()
            .heat_transfer_coefficient_with_sensitivities(
                ThermalConductivity::new(0.6), Length::new(0.02),
            ).unwrap();
        close(evaluation.sensitivity.d_qoi["Pr"] * 7.0 / evaluation.qoi,
            baseline.sensitivity.d_qoi["Pr"] * 7.0 / baseline.qoi + 0.11, 1e-13);
    }
}

#[test]
fn air_tangents_use_absolute_temperatures_and_retain_actual_coordinate_provenance() {
    let heat = |bulk, wall| section(27.0, 125.0, 0.71, air(bulk, wall))
        .heat_transfer_coefficient_with_sensitivities(
            ThermalConductivity::new(0.026), Length::new(0.02),
        ).unwrap();
    let nominal = heat(300.0, 400.0);
    close(nominal.sensitivity.d_qoi["T_bulk"], (heat(300.001, 400.0).qoi - heat(299.999, 400.0).qoi) / 0.002, 3e-8);
    close(nominal.sensitivity.d_qoi["T_wall"], (heat(300.0, 400.001).qoi - heat(300.0, 399.999).qoi) / 0.002, 3e-8);
    close(nominal.sensitivity.d_qoi["T_bulk"] * 300.0
        + nominal.sensitivity.d_qoi["T_wall"] * 400.0, 0.0, 2e-13);
    let scaled = heat(600.0, 800.0);
    assert_eq!(nominal.qoi, scaled.qoi);
    assert_ne!(nominal.provenance, scaled.provenance);
    assert_eq!(nominal.sensitivity.d_qoi["T_wall"], 2.0 * scaled.sensitivity.d_qoi["T_wall"]);
    let equal = heat(300.0, 300.0);
    assert!(equal.sensitivity.d_qoi["T_wall"] < 0.0);
    assert!(equal.sensitivity.d_qoi["T_bulk"] > 0.0);
}

#[test]
fn property_admissions_are_explicit_and_only_moving_endpoint_partials_are_omitted() {
    let base = || evaluate_turbulent_circular_duct_section(100_000.0, 1.0, 0.0, 100.0).unwrap();
    for correction in [Correction::Liquid { wall_prandtl: 20.0 }, Correction::Liquid { wall_prandtl: 0.05 }] {
        let h = base().with_property_correction(correction).unwrap()
            .heat_transfer_coefficient_with_sensitivities(ThermalConductivity::new(1.0), Length::new(1.0)).unwrap();
        assert!(!h.sensitivity.d_qoi.contains_key("Pr"));
        assert!(!h.sensitivity.d_qoi.contains_key("Pr_wall"));
        assert!(h.sensitivity.d_qoi.contains_key("Re"));
        assert!(h.sensitivity.d_qoi.contains_key("k"));
    }
    for correction in [air(200.0, 400.0), air(600.0, 400.0)] {
        let h = base().with_property_correction(correction).unwrap()
            .heat_transfer_coefficient_with_sensitivities(ThermalConductivity::new(1.0), Length::new(1.0)).unwrap();
        assert!(!h.sensitivity.d_qoi.contains_key("T_bulk"));
        assert!(!h.sensitivity.d_qoi.contains_key("T_wall"));
        assert!(h.sensitivity.d_qoi.contains_key("Pr"));
    }
    for correction in [Correction::Liquid { wall_prandtl: 20.01 }, Correction::Liquid { wall_prandtl: 0.049 },
        air(199.0, 400.0), air(601.0, 400.0)] {
        assert!(matches!(base().with_property_correction(correction), Err(CorrelationError::OutOfDomain { .. })));
    }
    for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(base().with_property_correction(Correction::Liquid { wall_prandtl: invalid }).is_err());
        assert!(base().with_property_correction(air(invalid, 300.0)).is_err());
        assert!(base().with_property_correction(air(300.0, invalid)).is_err());
    }
    assert!(evaluate_turbulent_circular_duct_section(9999.0, 1.0, 0.0, 100.0).is_err());
    assert!(evaluate_turbulent_circular_duct_section(100_000.0, 1.0, 1.0, 100.0).is_err());
}
