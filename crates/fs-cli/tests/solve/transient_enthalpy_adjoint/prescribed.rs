//! Prescribed baths retain both their direct objective term and their complete
//! free-node heat-transfer lift. A one-free-node tetrahedron supplies a scalar
//! backward-Euler oracle independent of the production derivative operators.

use super::*;

const DBATH: &str = "dtemperature_dprescribed_temperature_k_per_k";
const DK: &str = "dtemperature_dconductivity_multiplier_k";

fn baths(run: &ThermalRun) -> &[J] {
    run.adjoint().path(&["prescribed_temperature", "boundaries"])
        .unwrap().as_array().unwrap()
}

fn set_bath(spec: &mut ProjectSpec, target: &str, value: f64) {
    let setup = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.boundaries.iter_mut().find(|row| row.target == target).unwrap().condition =
        ThermalBoundaryCondition::FixedTemperature {
            temperature: QtyAny::new(value, dims::TEMPERATURE),
        };
}

fn partial_project(delta: f64) -> (ProjectSpec, Vec<u8>, CardPackSet) {
    let bytes = tetra_stl();
    let cards = transient_product::transient_conductivity_cards("constant-scalar");
    let mut spec = transient_product::driven_face_project(&bytes);
    install_storage(&mut spec, vec![material("air", 10.0 + 0.2 * delta, 300.0 + 2.0 * delta, false)]);
    rebind_to(&mut spec, &cards);
    time(&mut spec).enthalpy.as_mut().unwrap().materials[0].material_card_identity =
        fs_blake3::ContentHash::from_hex(&cards.materials()[0].card().to_hex()).unwrap();
    set_bath(&mut spec, "driven-face", 295.0 + 0.5 * delta);
    spec.power.as_mut().unwrap()[0].watts.value = 4.0 + 3.0 * delta;
    (spec, bytes, cards)
}

#[test]
fn g1_native_enthalpy_prescribed_adjoint_keeps_direct_fixed_objective_and_masks_loads() {
    let (mut spec, bytes) = uniform_project(10.0, 0.0, 40.0, false);
    set_bath(&mut spec, "air", 330.0);
    let solved = run(&spec, Some(&bytes), &fixture_cards(), Some(FINAL));
    check_final_adjoint(&solved);
    for &temperature in &solved.temperature { close(temperature, 330.0, 1e-10); }
    let rows = baths(&solved);
    assert_eq!(rows.len(), 1);
    assert_eq!(rows[0].str_field("target"), Some("air"));
    assert_eq!(rows[0].str_field("status"), Some("available"));
    close(n(&rows[0], "mapped_vertices"), 4.0, 0.0);
    close(n(&rows[0], DBATH), 1.0, 1e-10);
    let gradient = solved.adjoint();
    close(n(gradient, "dtemperature_duniform_initial_specific_enthalpy_k_kg_per_j"), 0.0, 1e-12);
    close(n(gradient, "dtemperature_dreference_density_multiplier_k"), 0.0, 1e-12);
    for value in values(gradient, "dtemperature_dinitial_specific_enthalpy_k_kg_per_j") {
        close(value, 0.0, 1e-12);
    }
    close(n(region(gradient, "regional_power", "air"),
        "dtemperature_duniform_added_watts_k_per_w"), 0.0, 1e-12);
    close(n(region(gradient, "regional_storage", "air"),
        "dtemperature_dreference_density_k_m3_per_kg"), 0.0, 1e-12);
    let conductivity = gradient.path(&["material_controls", "conductivity_regions"])
        .unwrap().as_array().unwrap();
    close(n(&conductivity[0], DK), 0.0, 1e-12);
    for interval in gradient.get("intervals").unwrap().as_array().unwrap() {
        close(n(interval, "dtemperature_dsource_multiplier_k"), 0.0, 1e-12);
    }
}

#[test]
fn g1_native_enthalpy_prescribed_adjoint_matches_free_apex_balance_and_fresh_resolves() {
    let (spec, bytes, cards) = partial_project(0.0);
    let solved = run(&spec, Some(&bytes), &cards, Some(FINAL));
    check_final_adjoint(&solved);
    let gradient = solved.adjoint();
    close(n(gradient, "active_ties"), 1.0, 0.0);
    close(n(&baths(&solved)[0], "mapped_vertices"), 3.0, 0.0);
    // Reference tetrahedron volume is 1/6. Its free origin has lumped
    // capacity rho*Cp/24, conductance k/2 to the opposite prescribed face,
    // and receives exactly one quarter of the total uniform source power.
    let capacity = 10.0 * CP / 24.0;
    let conductance = 2.0 / 2.0;
    let bath = 295.0;
    let mut temperature = 330.0;
    let (mut db, mut dh, mut drho, mut dp, mut dk) = (0.0, 1.0 / CP, 0.0, 0.0, 0.0);
    for interval in gradient.get("intervals").unwrap().as_array().unwrap() {
        let dt = n(interval, "end_s") - n(interval, "start_s");
        let denominator = capacity + dt * conductance;
        let a = capacity / denominator;
        let next = (capacity * temperature + dt * (conductance * bath + 4.0 / 4.0))
            / denominator;
        db = a * db + dt * conductance / denominator;
        dh *= a;
        drho = a * drho + (capacity / 10.0) * (temperature - next) / denominator;
        dp = a * dp + dt / (4.0 * denominator);
        dk = a * dk - dt * conductance * (next - bath) / denominator;
        temperature = next;
    }
    close(solved.final_temperature(), temperature, 2e-8);
    close(n(&baths(&solved)[0], DBATH), db, 2e-9);
    let storage = region(gradient, "regional_storage", "air");
    close(initial_gradient(storage), dh, 2e-9);
    close(density_gradient(storage), drho, 2e-9);
    let power = n(region(gradient, "regional_power", "air"),
        "dtemperature_duniform_added_watts_k_per_w");
    close(power, dp, 2e-9);
    let conductivity = gradient.path(&["material_controls", "conductivity_regions"])
        .unwrap().as_array().unwrap();
    close(n(&conductivity[0], DK), dk, 2e-8);
    let active = n(gradient, "active_vertex") as usize;
    for (vertex, value) in values(gradient,
        "dtemperature_dinitial_specific_enthalpy_k_kg_per_j").into_iter().enumerate()
    {
        close(value, if vertex == active { dh } else { 0.0 }, 2e-9);
    }
    let evaluate = |delta| {
        let (changed, bytes, cards) = partial_project(delta);
        run(&changed, Some(&bytes), &cards, None).final_temperature()
    };
    let finite_difference = (evaluate(0.01) - evaluate(-0.01)) / 0.02;
    let derivative = 0.5 * n(&baths(&solved)[0], DBATH)
        + 2.0 * initial_gradient(storage) + 0.2 * density_gradient(storage) + 3.0 * power;
    close(derivative, finite_difference, 2e-6);
}

#[test]
fn g1_native_enthalpy_prescribed_adjoint_initial_peak_has_zero_future_bath_effect() {
    let (spec, bytes, cards) = partial_project(0.0);
    let solved = run(&spec, Some(&bytes), &cards, Some(PEAK));
    let gradient = solved.adjoint();
    close(n(gradient, "value_k"), 330.0, 1e-12);
    close(n(gradient, "state_index"), 0.0, 0.0);
    close(n(gradient, "adjoint_krylov_iterations"), 0.0, 0.0);
    close(n(&baths(&solved)[0], DBATH), 0.0, 0.0);
    close(n(gradient.get("prescribed_temperature").unwrap(), "control_steps"), 0.0, 0.0);
    close(n(gradient, "dtemperature_duniform_initial_specific_enthalpy_k_kg_per_j"),
        1.0 / CP, 1e-12);
    assert!(solved.final_temperature() < 330.0);
}

#[test]
fn g1_native_enthalpy_prescribed_adjoint_shared_vertices_refuse_independent_bath_coordinates() {
    let (mut spec, bytes, cards) = partial_project(0.0);
    spec.assembly.as_mut().unwrap().push(EntityDecl::Surface {
        parent: "enclosure".into(), name: "second-bath".into(),
        display: "Second prescribed face".into(), expect_id: None,
    });
    spec.assignments.as_mut().unwrap().push(GeometryAssignment {
        artifact: "enclosure".into(), target: "second-bath".into(), length_unit: "m".into(),
        selector: MeshSelector::HalfSpace {
            normal: [1.0, 0.0, 0.0], offset: 0.0,
            side: HalfSpaceSide::AtMost, tolerance: 0.0,
        },
        allow_overlap: true,
    });
    spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries.push(ThermalBoundary {
        target: "second-bath".into(),
        condition: ThermalBoundaryCondition::FixedTemperature {
            temperature: QtyAny::new(295.0, dims::TEMPERATURE),
        },
    });
    let solved = run(&spec, Some(&bytes), &cards, Some(FINAL));
    assert_eq!(baths(&solved).len(), 2);
    for row in baths(&solved) {
        assert_eq!(row.str_field("status"), Some("unsupported"));
        assert!(row.str_field("reason").unwrap().contains("shares prescribed vertices"));
        assert!(row.get(DBATH).is_none(), "conflicting controls must not publish a numeric zero");
    }
}
