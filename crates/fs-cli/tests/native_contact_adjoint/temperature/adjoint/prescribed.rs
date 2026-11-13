//! Prescribed values and regional means reuse the complete live-contact dual.
use super::*;

const BOUNDARY: &str = "temperature-max-boundary-adjoint";
const COMBINED: &str = "temperature-max-contact-boundary-adjoint";
const MEAN: &str = "temperature-volume-mean-adjoint";

fn fixed_temperature(f: &mut Fixture, value: f64) {
    let B::FixedTemperature { temperature } = &mut f.project.cooling.as_mut().unwrap()
        .conduction.as_mut().unwrap().boundaries[0].condition else { panic!("original cold reservoir"); };
    temperature.value = value;
}

#[test]
fn live_contact_fixed_temperature_and_volume_mean_controls_match_physical_resolves() {
    let f = Fixture::new(true, false);
    let (plain, original) = solve(&f, &[1.0], false, None, 150);
    let (combined, field) = solve(&f, &[1.0], false, Some(COMBINED), 151);
    let (boundary, boundary_field) = solve(&f, &[1.0], false, Some(BOUNDARY), 152);
    assert_eq!(original.get("temperature"), field.get("temperature"));
    assert_eq!(boundary_field.get("temperature"), field.get("temperature"));
    for receipt in [&combined, &boundary] {
        for key in ["energy", "interfaces", "temperature_contact", "solver"] {
            assert_eq!(
                forward_section(&plain, key),
                forward_section(receipt, key),
                "adjoint changed primal {key}"
            );
        }
    }
    let report = combined.get("nominal_adjoint").unwrap();
    let boundary = boundary.get("nominal_adjoint").unwrap();
    assert_eq!(report.str_field("output"), Some(COMBINED));
    assert_eq!(boundary.str_field("output"), Some(BOUNDARY));
    let shared: Vec<_> = report.get("parameters").unwrap().as_array().unwrap().iter()
        .filter(|row| row.str_field("target") != Some(TARGET)).collect();
    assert_eq!(shared, boundary.get("parameters").unwrap().as_array().unwrap().iter().collect::<Vec<_>>());
    for key in ["selected_vertex", "value_k", "mode", "dual_iterations", "true_relative_residual"] {
        assert_eq!(report.get(key), boundary.get(key), "both reports reuse the same dual: {key}");
    }
    let selected = report.f64_field("selected_vertex").unwrap() as usize;
    let actual = derivative(report, "fixed-temperature", "cold");
    assert!(actual > 0.0 && derivative(report, TARGET, "cold-hot-joint") > 0.0);
    let step = 0.01;
    let mut values = Vec::new();
    for (side, sign) in [-1.0, 1.0].into_iter().enumerate() {
        let mut changed = Fixture { dir: f.dir.clone(), project: f.project.clone() };
        fixed_temperature(&mut changed, 293.15 + sign * step);
        let (_, field) = solve(&changed, &[1.0], false, None, 153 + side);
        values.push(value(&field, selected));
    }
    agrees(actual, (values[1] - values[0]) / (2.0 * step));

    // The same thermal/contact operator also accepts the existing P1 regional
    // volume load. Only the objective changes; no maximum node is substituted.
    let (mean, mean_field) = solve(&f, &[1.0], false, Some(MEAN), 155);
    assert_eq!(mean_field.get("temperature"), field.get("temperature"));
    let mean = mean.get("nominal_adjoint").unwrap();
    assert_eq!(mean.str_field("output"), Some(MEAN));
    assert_eq!(mean.str_field("functional"), Some("region-volume-mean-temperature"));
    assert_eq!(mean.get("selected_vertex"), Some(&J::Null));
    assert!(mean.f64_field("value_k").unwrap() > 293.15);
    assert!(mean.f64_field("value_k").unwrap() <= value(&field, selected));
    for (case, (target, entity, step)) in [("fixed-temperature", "cold", 0.01),
        (TARGET, "cold-hot-joint", 0.002)].into_iter().enumerate()
    {
        let actual = derivative(mean, target, entity);
        assert!(actual > 0.0);
        let mut values = Vec::new();
        for (side, sign) in [-1.0, 1.0].into_iter().enumerate() {
            let mut changed = Fixture { dir: f.dir.clone(), project: f.project.clone() };
            let scale = if target == TARGET { 1.0 + sign * step } else { 1.0 };
            if target == "fixed-temperature" { fixed_temperature(&mut changed, 293.15 + sign * step); }
            let (changed, _) = solve(&changed, &[scale], false, Some(MEAN), 156 + 2 * case + side);
            // Read the forward volume functional of each independently solved
            // project, never a derivative or a maximum of the perturbed field.
            values.push(changed.get("nominal_adjoint").unwrap().f64_field("value_k").unwrap());
        }
        agrees(actual, (values[1] - values[0]) / (2.0 * step));
    }
}

#[test]
fn a_selected_prescribed_temperature_keeps_its_direct_objective_term() {
    let mut f = Fixture::new(true, false);
    // The unpowered hot body loses heat to a colder reservoir, so the cold
    // body's hottest node is on its original prescribed-temperature boundary.
    f.project.power.as_mut().unwrap()[0].watts.value = 0.0;
    f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries[1].condition =
        B::Convection { coefficient: QtyAny::new(8.0, dims::HEAT_TRANSFER_COEFFICIENT),
            reference_temperature: QtyAny::new(280.0, dims::TEMPERATURE) };
    f.project.requirements.as_mut().unwrap()[0].region = "cold".into();
    let (receipt, _) = solve(&f, &[1.0], false, Some(COMBINED), 160);
    let report = receipt.get("nominal_adjoint").unwrap();
    assert_eq!(report.f64_field("value_k"), Some(293.15));
    assert_eq!(report.f64_field("dual_iterations"), Some(0.0));
    assert_eq!(derivative(report, "fixed-temperature", "cold"), 1.0);
    assert_eq!(derivative(report, TARGET, "cold-hot-joint"), 0.0);
    assert_eq!(derivative(report, "power", "hot"), 0.0);
}
