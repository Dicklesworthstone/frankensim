//! G1 signed reservoir balances and G3 composition of native mixed convection.
//! These manufactured boundary laws exercise the real import/solve pipeline.

use super::*;

#[allow(dead_code)]
#[path = "../../src/json_read.rs"]
mod json;
use json::JsonValue as J;

#[path = "mixed_convection/adjoint.rs"]
mod adjoint;

fn n(row: &J, key: &str) -> f64 {
    row.f64_field(key).unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual}, expected={expected}, tolerance={tolerance}"
    );
}

fn project(bytes: &[u8], ambient_k: f64) -> ProjectSpec {
    let mut spec = conjugate_fixture_project(7, bytes, "convection.gnielinski");
    spec.solver.as_mut().unwrap().tolerance_rel = 1e-9;
    spec.assembly.as_mut().unwrap().push(EntityDecl::Surface {
        parent: "enclosure".to_string(),
        name: "natural-wall".to_string(),
        display: "Manufactured natural-convection boundary".to_string(),
        expect_id: None,
    });
    let assignments = spec.assignments.as_mut().unwrap();
    assignments[0].allow_overlap = true;
    assignments.push(GeometryAssignment {
        artifact: "enclosure".to_string(),
        target: "natural-wall".to_string(),
        length_unit: "m".to_string(),
        selector: MeshSelector::HalfSpace {
            normal: [0.0, 0.0, 1.0],
            offset: 0.0,
            side: HalfSpaceSide::AtMost,
            tolerance: 0.0,
        },
        allow_overlap: true,
    });
    spec.cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap()
        .boundaries
        .push(ThermalBoundary {
            target: "natural-wall".to_string(),
            condition: ThermalBoundaryCondition::NaturalConvection {
                characteristic_length: QtyAny::new(0.1, fs_project::spec::dims::LENGTH),
                ambient_temperature: QtyAny::new(
                    ambient_k,
                    fs_project::spec::dims::TEMPERATURE,
                ),
                correlation: "convection.churchill-chu-vertical-plate".to_string(),
            },
        });
    assert!(spec.validate().is_empty(), "{:?}", spec.validate());
    spec
}

fn run(spec: &ProjectSpec, bytes: &[u8], cards: &CardPackSet) -> (J, J) {
    let ledger = Ledger::open(":memory:").unwrap();
    import_fixture(&ledger, spec, bytes.to_vec());
    let (_, _, receipt, field) =
        run_conjugate_to_completion(&ledger, &decode(spec), cards);
    (J::parse(&receipt).unwrap(), J::parse(&field).unwrap())
}

// Independent direct evaluation of the declared manufactured Churchill-Chu
// law. This does not call fs-convection or its deterministic power routine.
fn physical_htc(wall_k: f64, ambient_k: f64) -> f64 {
    let film_k = (wall_k + ambient_k) * 0.5;
    let density = 101_325.0 / (287.05 * film_k);
    let nu = 1.846e-5 / density;
    let ra = 9.80665 * (wall_k - ambient_k).abs() * 0.1_f64.powi(3)
        / (film_k * nu * (nu / 0.707));
    let nusselt =
        (0.825 + 0.387 * ra.powf(1.0 / 6.0)
            / (1.0 + (0.492_f64 / 0.707).powf(9.0 / 16.0)).powf(8.0 / 27.0))
        .powi(2);
    nusselt * 0.0263 / 0.1
}

fn check_balance(receipt: &J, ambient_k: f64) -> f64 {
    let natural = receipt.get("natural_convection").unwrap();
    let law = &natural.get("laws").unwrap().as_array().unwrap()[0];
    let htc = physical_htc(n(law, "mean_wall_k"), ambient_k);
    close(n(law, "htc_w_m2_k"), htc, 1e-9);
    let natural_heat = n(law, "heat_rate_w");
    // The selected z=0 face is a right triangle with legs 1 m.
    close(
        natural_heat,
        0.5 * htc * (n(law, "mean_wall_k") - ambient_k),
        1e-7,
    );
    let air = receipt.get("conjugate").unwrap();
    let air_heat = n(air, "air_total_w");
    close(
        air_heat,
        n(air, "mass_flow_kg_s") * 1007.0 * (n(air, "outlet_k") - n(air, "inlet_k")),
        1e-7,
    );
    close(
        n(air, "solid_total_w"),
        air_heat,
        n(air, "balance_tolerance_w"),
    );
    let radiation_heat = receipt
        .get("radiation")
        .and_then(|radiation| radiation.f64_field("radiative_out_w"))
        .unwrap_or(0.0);
    close(air_heat + natural_heat + radiation_heat, 5.0, 1e-6);
    close(n(receipt.get("energy").unwrap(), "robin_out_w"), 5.0, 1e-6);
    let coupling = natural.get("coupling").unwrap();
    assert!(n(coupling, "max_natural_heat_gate_ratio") <= 1.0);
    assert!(n(coupling, "max_relative_coefficient_change") <= 1e-10);
    close(
        n(coupling, "total_solid_responses"),
        n(coupling, "total_air_iterations") + n(natural, "iterations"),
        0.0,
    );
    assert_eq!(
        receipt.get("solver_control").unwrap().str_field("status"),
        Some("unsupported-model")
    );
    natural_heat
}

#[test]
fn g1_mixed_steady_convection_preserves_each_reservoir_and_signed_heat() {
    let bytes = tetra_stl();
    for ambient_k in [293.15, 300.0] {
        let spec = project(&bytes, ambient_k);
        let (receipt, field) = run(&spec, &bytes, &fixture_cards());
        let natural_heat = check_balance(&receipt, ambient_k);
        if ambient_k == 293.15 {
            assert!(natural_heat > 0.0);
        } else {
            assert!(natural_heat < 0.0);
            assert!(n(receipt.get("conjugate").unwrap(), "air_total_w") > 5.0);
        }
        let (replayed, replayed_field) = run(&spec, &bytes, &fixture_cards());
        assert_eq!(field, replayed_field);
        assert_eq!(receipt.get("natural_convection"), replayed.get("natural_convection"));
        assert_eq!(receipt.get("conjugate"), replayed.get("conjugate"));
    }
}

#[test]
fn g3_mixed_convection_keeps_radiation_out_of_both_convective_heat_accounts() {
    let bytes = tetra_stl();
    let mut spec = project(&bytes, 300.0);
    let (_, baseline) = run(&spec, &bytes, &fixture_cards());
    let (cards, card) = radiation_product::emissivity_cards(&fixture_cards(), 0.85);
    radiation_product::declare(&mut spec, &card, "natural-wall", 293.15);
    let (receipt, field) = run(&spec, &bytes, &cards);
    let natural_heat = check_balance(&receipt, 300.0);
    assert_ne!(field.get("temperature"), baseline.get("temperature"));
    let radiation = receipt.get("radiation").unwrap();
    assert!(n(radiation, "radiative_out_w") > 0.0);
    close(
        n(radiation, "convective_out_w"),
        n(receipt.get("conjugate").unwrap(), "solid_total_w") + natural_heat,
        1e-6,
    );
}

#[test]
fn g0_mixed_adaptive_goal_refuses_before_solving() {
    let bytes = tetra_stl();
    for (mode, expected_code) in [
        ("adaptive", "cli-solve-conduction-mixed-adaptive"),
    ] {
        let mut spec = project(&bytes, 300.0);
        match mode {
            "adaptive" => spec.solver.as_mut().unwrap().fidelity = "adaptive".to_string(),
            _ => unreachable!(),
        }
        let ledger = Ledger::open(":memory:").unwrap();
        import_fixture(&ledger, &spec, bytes.clone());
        let refusal = run_solve(
            &ledger,
            &CancelGate::new_clock_free(),
            &mut benign_clock(),
            &decode(&spec),
            &fixture_cards(),
            &mut Vec::new(),
        )
        .expect_err("an unimplemented coupled adaptive goal must refuse");
        assert_eq!(refusal.code, expected_code, "{refusal:?}");
    }
}
