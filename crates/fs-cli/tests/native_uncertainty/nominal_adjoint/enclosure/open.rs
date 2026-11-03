//! Native open-enclosure controls traverse the same accepted solid/air state.
use super::*;
use fs_project::{ThermalBoundaryCondition as B, spec::dims};
use fs_qty::QtyAny;

fn open(project: &mut fs_project::ProjectSpec) {
    let declaration = project.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
        .enclosure_radiation.as_mut().unwrap();
    declaration.ambient_temperature = Some(QtyAny::new(280.0, dims::TEMPERATURE));
    declaration.view_factors = vec![vec![0.0, 0.25], vec![0.25, 0.0]];
}

fn move_bath(project: &mut fs_project::ProjectSpec, delta: f64) {
    project.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
        .enclosure_radiation.as_mut().unwrap().ambient_temperature.as_mut().unwrap().value += delta;
}

fn unchanged(f: &Fixture, baseline: &JsonValue, derivative: &JsonValue) -> usize {
    assert_eq!(field(f, baseline), field(f, derivative));
    for key in ["energy", "radiation", "conjugate", "natural_convection"] {
        assert_eq!(baseline.get(key), derivative.get(key), "{key}");
    }
    let report = derivative.get("nominal_adjoint").unwrap();
    assert_eq!(report.str_field("mode"), Some("open-enclosure-full-radiosity-feedback"));
    let parameters = report.get("parameters").unwrap().as_array().unwrap();
    assert_eq!(parameters.iter().filter(|row|
        row.str_field("target") == Some("radiation-emissivity")).count(), 2);
    let bath: Vec<_> = parameters.iter().filter(|row|
        row.str_field("target") == Some("radiation-reservoir-temperature")).collect();
    assert_eq!(bath.len(), 1, "one shared bath, independent of surface emissivity controls");
    assert_eq!(bath[0].str_field("entity"), Some("enclosure-environment"));
    assert!(bath[0].f64_field("derivative").unwrap() > 0.0);
    report.f64_field("selected_vertex").unwrap() as usize
}

#[test]
fn native_open_vacuum_enclosure_prescribed_lift_retains_black_bath_feedback() {
    let mut f = fixture();
    open(&mut f.project);
    // Keep the selected maximum on a free node. Otherwise the cold external
    // bath would make the prescribed wall the maximum and its bath partial zero.
    f.project.power.as_mut().unwrap()[0].watts.value = 100.0;
    f.project.solver.as_mut().unwrap().fidelity = "ladder".into();
    for (boundary, flux) in f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
        .boundaries.iter_mut().zip([-8.0, 2.0]) {
        boundary.condition = B::HeatFlux { outward_flux: QtyAny::new(flux, dims::HEAT_FLUX) };
    }
    extra_y_face(&mut f.project, "controlled-wall", B::FixedTemperature {
        temperature: QtyAny::new(290.0, dims::TEMPERATURE),
    });
    let baseline = native(&f, &f.project, 220);
    let mut requested = f.project.clone();
    requested.outputs.as_mut().unwrap().push(fs_project::spec::OutputRequest {
        name: "temperature-max-boundary-adjoint".into(), kind: "report".into(), region: None,
    });
    let receipt = native(&f, &requested, 221);
    let vertex = unchanged(&f, &baseline, &receipt);
    let parameters = receipt.get("nominal_adjoint").unwrap().get("parameters").unwrap()
        .as_array().unwrap();
    assert!(!parameters.iter().any(|row| matches!(row.str_field("target"),
        Some("convection-coefficient" | "convection-temperature"))));
    for (case, target) in [(0, "fixed-temperature"), (1, "radiation-reservoir-temperature")] {
        let mut values = Vec::new();
        for (side, sign) in [-1.0, 1.0].into_iter().enumerate() {
            let mut changed = f.project.clone();
            if case == 0 {
                let B::FixedTemperature { temperature } = &mut changed.cooling.as_mut().unwrap()
                    .conduction.as_mut().unwrap().boundaries.last_mut().unwrap().condition
                    else { panic!("prescribed bath"); };
                temperature.value += sign * 0.01;
            } else { move_bath(&mut changed, sign * 0.01); }
            values.push(at_vertex(&f, &native(&f, &changed, 222 + 2 * case + side), vertex));
        }
        close_control(coefficient(&receipt, target), (values[1] - values[0]) / 0.02, target);
    }
}

#[test]
fn native_open_enclosure_mixed_air_and_natural_feedback_uses_the_same_complete_dual() {
    let mut f = mixed_fixture();
    open(&mut f.project);
    let baseline = native(&f, &f.project, 240);
    let mut requested = f.project.clone();
    request(&mut requested);
    let receipt = native(&f, &requested, 241);
    let vertex = unchanged(&f, &baseline, &receipt);
    for (case, target) in [(0, "radiation-reservoir-temperature"),
        (1, "air-inlet-temperature"), (2, "natural-convection-ambient")]
    {
        let mut values = Vec::new();
        for (side, sign) in [-1.0, 1.0].into_iter().enumerate() {
            let mut changed = f.project.clone();
            if case == 0 { move_bath(&mut changed, sign * 0.01); }
            else {
                for boundary in &mut changed.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries {
                    match &mut boundary.condition {
                        B::AirflowConvection { inlet_temperature, .. } if case == 1 =>
                            inlet_temperature.value += sign * 0.01,
                        B::NaturalConvection { ambient_temperature, .. } if case == 2 =>
                            ambient_temperature.value += sign * 0.01,
                        _ => {},
                    }
                }
            }
            values.push(at_vertex(&f, &native(&f, &changed, 242 + 2 * case + side), vertex));
        }
        close_control(coefficient(&receipt, target), (values[1] - values[0]) / 0.02, target);
    }
}
