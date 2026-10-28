//! Native steady mixed-flux lowering and the empty-Robin enclosure dual.
use super::*;
use fs_project::{ThermalBoundaryCondition as B, spec::dims};
use fs_qty::QtyAny;

fn flux(value: f64) -> B {
    B::HeatFlux { outward_flux: QtyAny::new(value, dims::HEAT_FLUX) }
}

fn heat(receipt: &JsonValue, external_neumann_w: f64) {
    let energy = receipt.get("energy").unwrap();
    let value = |name| energy.f64_field(name).unwrap();
    assert!((value("neumann_out_w") - external_neumann_w).abs() < 1e-12,
        "the external flux must exclude reflected enclosure exchange");
    let balance = value("source_w") + value("dirichlet_in_w")
        - value("neumann_out_w") - value("robin_out_w");
    assert!((balance - value("closure_w")).abs() < 1e-8);
    assert!(balance.abs() < 1e-6, "physical steady heat balance {balance:e}");
}

#[test]
fn native_mixed_flux_enclosure_keeps_external_heat_and_complete_power_dual() {
    let mut f = fixture();
    f.project.power.as_mut().unwrap()[0].duty = 0.37;
    let setup = f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.boundaries[0].condition = flux(-8.0);
    let flux_target = setup.boundaries[0].target.clone();
    let baseline = native(&f, &f.project, 100);
    heat(&baseline, -4.0); // Explicit -8 W/m2 over a 0.5 m2 tetra face.
    let mut requested = f.project.clone();
    request(&mut requested);
    let receipt = native(&f, &requested, 101);
    assert_eq!(field(&f, &baseline), field(&f, &receipt));
    assert_eq!(baseline.get("energy"), receipt.get("energy"));
    heat(&receipt, -4.0);
    let report = receipt.get("nominal_adjoint").unwrap();
    assert_eq!(report.str_field("mode"), Some("enclosure-full-radiosity-feedback"));
    let parameters = report.get("parameters").unwrap().as_array().unwrap();
    assert!(!parameters.iter().any(|row|
        row.str_field("entity") == Some(flux_target.as_str())
            && matches!(row.str_field("target"),
                Some("convection-coefficient" | "convection-temperature"))),
        "a Neumann radiation patch cannot acquire an invented convection control");
    let radiation = receipt.get("radiation").unwrap();
    let patches = radiation.get("surfaces").unwrap().as_array().unwrap();
    assert!(patches[0].f64_field("nonlinear_heat_w").unwrap().abs() > 1e-6,
        "this mixed-boundary test must transfer internal heat through the flux patch");
    let selected = report.f64_field("selected_vertex").unwrap() as usize;
    let step = 0.02;
    let mut values = Vec::new();
    for (side, sign) in [-1.0, 1.0].into_iter().enumerate() {
        let mut changed = f.project.clone();
        changed.power.as_mut().unwrap()[0].watts.value += sign * step;
        values.push(at_vertex(&f, &native(&f, &changed, 102 + side), selected));
    }
    close_control(coefficient(&receipt, "power"), (values[1] - values[0]) / (2.0 * step),
        "mixed-flux power");
}

#[test]
fn native_all_flux_enclosure_uses_a_real_bath_and_prescribed_lift_with_no_robin_controls() {
    let mut f = fixture();
    // Refine the tetra so the radiating means are independent of the bath.
    f.project.solver.as_mut().unwrap().fidelity = "ladder".into();
    let setup = f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.boundaries[0].condition = flux(-8.0);
    setup.boundaries[1].condition = flux(2.0);
    // This separately declared y=0 bath shares vertices with both radiating
    // faces. It supplies the absolute steady anchor; the enclosure does not.
    extra_y_face(&mut f.project, "controlled-wall", B::FixedTemperature {
        temperature: QtyAny::new(280.0, dims::TEMPERATURE),
    });
    let baseline = native(&f, &f.project, 110);
    heat(&baseline, -3.0);
    assert!(baseline.get("mesh").unwrap().f64_field("elements").unwrap() > 1.0,
        "the bath-lift fixture requires a refined mesh");
    let patches = baseline.get("radiation").unwrap().get("surfaces")
        .unwrap().as_array().unwrap();
    assert!(patches.iter().any(|patch|
        patch.f64_field("nonlinear_heat_w").unwrap().abs() > 1e-6),
        "the refined vacuum fixture must exercise internal radiation");
    let mut requested = f.project.clone();
    requested.outputs.as_mut().unwrap().push(fs_project::spec::OutputRequest {
        name: "temperature-max-boundary-adjoint".into(), kind: "report".into(), region: None,
    });
    let receipt = native(&f, &requested, 111);
    assert_eq!(field(&f, &baseline), field(&f, &receipt));
    assert_eq!(baseline.get("energy"), receipt.get("energy"));
    heat(&receipt, -3.0);
    let report = receipt.get("nominal_adjoint").unwrap();
    let parameters = report.get("parameters").unwrap().as_array().unwrap();
    assert!(!parameters.iter().any(|row|
        matches!(row.str_field("target"),
            Some("convection-coefficient" | "convection-temperature"))));
    assert_eq!(parameters.iter().filter(|row|
        row.str_field("target") == Some("radiation-emissivity")).count(), 2);
    let selected = report.f64_field("selected_vertex").unwrap() as usize;
    let step = 0.01;
    let mut values = Vec::new();
    for (side, sign) in [-1.0, 1.0].into_iter().enumerate() {
        let mut changed = f.project.clone();
        let B::FixedTemperature { temperature } = &mut changed.cooling.as_mut().unwrap()
            .conduction.as_mut().unwrap().boundaries.last_mut().unwrap().condition
            else { panic!("explicit prescribed bath"); };
        temperature.value += sign * step;
        values.push(at_vertex(&f, &native(&f, &changed, 112 + side), selected));
    }
    close_control(coefficient(&receipt, "fixed-temperature"),
        (values[1] - values[0]) / (2.0 * step), "all-flux prescribed bath");
}

