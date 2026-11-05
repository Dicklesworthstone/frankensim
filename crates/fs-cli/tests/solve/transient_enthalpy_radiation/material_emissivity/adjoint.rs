//! Physical source-ordinate twins rebuild the original immutable material card.
//! No perturbed run consumes the reported adjoint or a frozen emissivity proxy.

use super::*;

const ORDINATE_DERIVATIVE: &str = "dtemperature_demissivity_ordinate_k";
const BATH_DERIVATIVE: &str = "dtemperature_dreservoir_temperature_k_per_k";

fn request(spec: &mut ProjectSpec, goal: &str) {
    spec.outputs.as_mut().unwrap().push(OutputRequest {
        name: goal.into(),
        kind: "report".into(),
        region: None,
    });
}

fn derivative_close(actual: f64, expected: f64) {
    close(actual, expected, 3e-7 + 3e-5 * expected.abs());
}

fn sources<'a>(solved: &'a ThermalRun, material: &MaterialFixture, shared: usize) -> &'a [J] {
    let gradient = solved.transient().get("adjoint").unwrap();
    let radiation = gradient.get("radiation").unwrap();
    assert_eq!(
        radiation.str_field("method"),
        Some("complete-pointwise-p1-material-ambient-enthalpy-tangent")
    );
    assert_eq!(radiation.str_field("spatial_model"), Some("pointwise"));
    assert_eq!(
        radiation.str_field("emissivity_response"),
        Some("local-temperature")
    );
    close(
        n(radiation, "control_steps"),
        n(gradient, "state_index"),
        0.0,
    );
    let rows = radiation
        .get("emissivity_sources")
        .unwrap()
        .as_array()
        .unwrap();
    assert_eq!(
        rows.len(),
        1,
        "one independent control vector per original card and claim"
    );
    let source = &rows[0];
    assert_eq!(source.str_field("card"), Some(material.card.as_str()));
    assert_eq!(source.str_field("claim"), Some(material.claim.as_str()));
    assert_eq!(
        source.str_field("material_state"),
        Some(material.state.as_str())
    );
    assert_eq!(source.str_field("model"), Some("temperature-curve"));
    assert_eq!(source.str_field("axis"), Some("T"));
    close(n(source, "shared_surface_count"), shared as f64, 0.0);
    let knots = source.get("knots").unwrap().as_array().unwrap();
    assert_eq!(
        knots.len(),
        SOURCE_KNOTS.len(),
        "retain all original source ordinates"
    );
    for (index, (row, &(temperature, value))) in knots.iter().zip(&SOURCE_KNOTS).enumerate() {
        close(n(row, "index"), index as f64, 0.0);
        close(n(row, "temperature_k"), temperature, 0.0);
        close(n(row, "emissivity"), value, 0.0);
    }
    let surfaces = radiation.get("surfaces").unwrap().as_array().unwrap();
    assert_eq!(surfaces.len(), shared);
    for row in surfaces {
        assert_eq!(row.str_field("card"), Some(material.card.as_str()));
        assert_eq!(row.str_field("claim"), Some(material.claim.as_str()));
        close(n(row, "source_index"), 0.0, 0.0);
        assert!(
            row.get("dtemperature_demissivity_k").is_none(),
            "a frozen per-patch coefficient is not an original shared source control"
        );
    }
    knots
}

/// Differentiate the scalar backward-Euler balance, including epsilon'(T).
/// This recurrence uses no material query, FEM operator or production dual.
fn recurrence(solved: &ThermalRun, reservoir: f64) -> ([f64; 3], f64, f64, f64) {
    let area = area_m2();
    let mut old = INITIAL_K;
    let mut ordinates = [0.0; 3];
    let mut bath = 0.0;
    let mut initial_h = 1.0 / CP;
    let mut added_watts = 0.0;
    for row in solved.rows("fine") {
        let dt = n(row, "dt_s");
        let temperature = step(old, dt, reservoir, &SOURCE_KNOTS, true);
        // These independently checked histories stay inside the last original
        // segment. The first authored ordinate remains a genuine zero control.
        assert!(temperature > 300.0 && temperature < 450.0);
        let fraction = (temperature - 300.0) / 150.0;
        let value = epsilon(temperature, &SOURCE_KNOTS);
        let slope = (SOURCE_KNOTS[2].1 - SOURCE_KNOTS[1].1) / 150.0;
        let black = SIGMA * area * (temperature.powi(4) - reservoir.powi(4));
        let storage = mass_kg(DENSITY) * CP / dt;
        let denominator = storage
            + 10.0 * area
            + 4.0 * SIGMA * area * value * temperature.powi(3)
            + slope * black;
        for (derivative, weight) in ordinates.iter_mut().zip([0.0, 1.0 - fraction, fraction]) {
            *derivative = (storage * *derivative - black * weight) / denominator;
        }
        bath = (storage * bath + 4.0 * SIGMA * area * value * reservoir.powi(3)) / denominator;
        initial_h = storage * initial_h / denominator;
        added_watts = (storage * added_watts + 1.0) / denominator;
        old = temperature;
    }
    (ordinates, bath, initial_h, added_watts)
}

#[test]
fn g1_original_emissivity_ordinate_adjoint_matches_scalar_recurrence_and_full_card_resolves() {
    let original = material(&SOURCE_KNOTS, "original-specimen");
    for reservoir in [270.0, 400.0] {
        let spec = fixture(&original, reservoir);
        let baseline = run(&spec, &original.cards);
        let mut differentiated = spec.clone();
        request(&mut differentiated, "temperature-final-adjoint");
        let solved = run(&differentiated, &original.cards);
        assert_eq!(solved.temperature, baseline.temperature);
        assert_eq!(solved.receipt.get("energy"), baseline.receipt.get("energy"));
        assert_eq!(
            solved.receipt.get("radiation"),
            baseline.receipt.get("radiation")
        );
        let gradient = solved.transient().get("adjoint").unwrap();
        assert!(n(gradient, "maximum_dual_relative_residual") <= 1e-9);
        close(
            n(gradient, "state_index"),
            n(solved.transient(), "fine_steps"),
            0.0,
        );
        let vertex = n(gradient, "active_vertex") as usize;
        close(n(gradient, "value_k"), solved.temperature[vertex], 0.0);
        let knots = sources(&solved, &original, 1);
        let (expected, expected_bath, expected_initial, expected_power) =
            recurrence(&solved, reservoir);
        for (index, row) in knots.iter().enumerate() {
            let actual = n(row, ORDINATE_DERIVATIVE);
            derivative_close(actual, expected[index]);
            let perturb = |delta: f64| {
                let mut shifted = SOURCE_KNOTS;
                shifted[index].1 += delta;
                let changed = material(&shifted, "original-specimen");
                assert_ne!(changed.card, original.card);
                assert_ne!(changed.claim, original.claim);
                // fixture() rebinds bulk material, exact caloric-card identity
                // and selected radiation claim together. Caloric knots, rho,
                // initial h and source abscissae remain physically unchanged.
                run(&fixture(&changed, reservoir), &changed.cards).temperature[vertex]
            };
            let h = 0.0002;
            derivative_close(actual, (perturb(h) - perturb(-h)) / (2.0 * h));
            if index == 0 {
                close(actual, 0.0, 0.0);
            } else {
                assert!(
                    actual.abs() > 0.1,
                    "source ordinate must alter the actual solution"
                );
                assert_eq!(actual.is_sign_positive(), reservoir > INITIAL_K);
            }
        }
        let surface = &gradient
            .get("radiation")
            .unwrap()
            .get("surfaces")
            .unwrap()
            .as_array()
            .unwrap()[0];
        let actual_bath = n(surface, BATH_DERIVATIVE);
        derivative_close(actual_bath, expected_bath);
        let perturb = |delta: f64| {
            let mut changed = spec.clone();
            changed
                .cooling
                .as_mut()
                .unwrap()
                .conduction
                .as_mut()
                .unwrap()
                .radiation
                .as_mut()
                .unwrap()
                .surfaces[0]
                .reservoir_temperature
                .value += delta;
            run(&changed, &original.cards).temperature[vertex]
        };
        let h = 0.002;
        derivative_close(actual_bath, (perturb(h) - perturb(-h)) / (2.0 * h));
        let actual_initial = n(
            gradient,
            "dtemperature_duniform_initial_specific_enthalpy_k_kg_per_j",
        );
        let power = &gradient.get("regional_power").unwrap().as_array().unwrap()[0];
        let actual_power = n(power, "dtemperature_duniform_added_watts_k_per_w");
        for (coordinate, actual, expected, h) in [
            (0, actual_initial, expected_initial, 0.02),
            (1, actual_power, expected_power, 0.05),
        ] {
            derivative_close(actual, expected);
            assert!(actual > 0.0);
            let perturb = |delta: f64| {
                let mut changed = spec.clone();
                let policy = time(&mut changed);
                if coordinate == 0 {
                    policy.enthalpy.as_mut().unwrap().materials[0]
                        .initial_specific_enthalpy
                        .value += delta;
                } else {
                    policy.power_schedules[0].steps[0].watts.value += delta;
                }
                run(&changed, &original.cards).temperature[vertex]
            };
            derivative_close(actual, (perturb(h) - perturb(-h)) / (2.0 * h));
        }
    }
}

#[test]
fn g1_initial_peak_keeps_original_emissivity_ordinates_and_future_baths_zero() {
    let original = material(&SOURCE_KNOTS, "original-specimen");
    let mut spec = fixture(&original, 270.0);
    request(&mut spec, "temperature-sampled-peak-adjoint");
    let solved = run(&spec, &original.cards);
    let gradient = solved.transient().get("adjoint").unwrap();
    close(n(gradient, "state_index"), 0.0, 0.0);
    close(n(gradient, "adjoint_krylov_iterations"), 0.0, 0.0);
    close(n(gradient, "value_k"), INITIAL_K, 1e-12);
    close(
        n(
            gradient,
            "dtemperature_duniform_initial_specific_enthalpy_k_kg_per_j",
        ),
        1.0 / CP,
        1e-14,
    );
    for knot in sources(&solved, &original, 1) {
        close(n(knot, ORDINATE_DERIVATIVE), 0.0, 0.0);
    }
    for surface in gradient
        .get("radiation")
        .unwrap()
        .get("surfaces")
        .unwrap()
        .as_array()
        .unwrap()
    {
        close(n(surface, BATH_DERIVATIVE), 0.0, 0.0);
    }
}

fn four_faces(material: &MaterialFixture, reservoir: f64) -> ProjectSpec {
    let mut spec = fixture(material, reservoir);
    spec.assignments.as_mut().unwrap()[0].allow_overlap = true;
    let setup = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.boundaries.clear();
    let prototype = setup.radiation.as_mut().unwrap().surfaces.remove(0);
    // Each exact plane selects the face opposite one regular-tetra vertex.
    // All four traces share one card and source claim, while each bath remains
    // an independent physical control. Identical baseline baths preserve the
    // original uniform scalar mode and total exterior area.
    for (index, normal) in [
        [1.0, 1.0, 1.0],
        [1.0, -1.0, -1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, -1.0, 1.0],
    ]
    .into_iter()
    .enumerate()
    {
        let target = format!("emitting-face-{index}");
        spec.assembly.as_mut().unwrap().push(EntityDecl::Surface {
            parent: "enclosure".into(),
            name: target.clone(),
            display: target.clone(),
            expect_id: None,
        });
        spec.assignments.as_mut().unwrap().push(GeometryAssignment {
            artifact: "enclosure".into(),
            target: target.clone(),
            length_unit: "m".into(),
            selector: MeshSelector::HalfSpace {
                normal,
                offset: -1.0,
                side: HalfSpaceSide::AtMost,
                tolerance: 0.0,
            },
            allow_overlap: true,
        });
        let setup = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
        setup.boundaries.push(ThermalBoundary {
            target: target.clone(),
            condition: ThermalBoundaryCondition::Convection {
                coefficient: QtyAny::new(10.0, dims::HEAT_TRANSFER_COEFFICIENT),
                reference_temperature: QtyAny::new(AMBIENT_K, dims::TEMPERATURE),
            },
        });
        let mut surface = prototype.clone();
        surface.name = format!("original-gray-face-{index}");
        surface.target = target;
        setup.radiation.as_mut().unwrap().surfaces.push(surface);
    }
    spec
}

#[test]
fn g1_shared_original_emissivity_sums_physical_faces_but_keeps_bath_controls_independent() {
    const RESERVOIR: f64 = 400.0;
    let original = material(&SOURCE_KNOTS, "original-specimen");
    let mut spec = four_faces(&original, RESERVOIR);
    request(&mut spec, "temperature-final-adjoint");
    let solved = run(&spec, &original.cards);
    let gradient = solved.transient().get("adjoint").unwrap();
    let vertex = n(gradient, "active_vertex") as usize;
    let knots = sources(&solved, &original, 4);
    let (expected, expected_bath, _, _) = recurrence(&solved, RESERVOIR);
    for (row, expected) in knots.iter().zip(expected) {
        derivative_close(n(row, ORDINATE_DERIVATIVE), expected);
    }
    let mut oracle = INITIAL_K;
    for row in solved.rows("fine") {
        oracle = step(oracle, n(row, "dt_s"), RESERVOIR, &SOURCE_KNOTS, true);
    }
    for &temperature in &solved.temperature {
        close(temperature, oracle, 2e-8);
    }
    let perturb_source = |delta: f64| {
        let mut shifted = SOURCE_KNOTS;
        shifted[2].1 += delta;
        let changed = material(&shifted, "original-specimen");
        run(&four_faces(&changed, RESERVOIR), &changed.cards).temperature[vertex]
    };
    let h = 0.0002;
    derivative_close(
        n(&knots[2], ORDINATE_DERIVATIVE),
        (perturb_source(h) - perturb_source(-h)) / (2.0 * h),
    );
    let surfaces = gradient
        .get("radiation")
        .unwrap()
        .get("surfaces")
        .unwrap()
        .as_array()
        .unwrap();
    let mut shared_bath = 0.0;
    for (index, surface) in surfaces.iter().enumerate() {
        let actual = n(surface, BATH_DERIVATIVE);
        assert!(actual > 0.0);
        shared_bath += actual;
        assert_eq!(
            surface.str_field("target"),
            Some(format!("emitting-face-{index}").as_str())
        );
        let perturb = |delta: f64| {
            let mut changed = four_faces(&original, RESERVOIR);
            changed
                .cooling
                .as_mut()
                .unwrap()
                .conduction
                .as_mut()
                .unwrap()
                .radiation
                .as_mut()
                .unwrap()
                .surfaces[index]
                .reservoir_temperature
                .value += delta;
            // Independent bath probes break uniformity. Compare the retained
            // active node, never the perturbed run's re-selected maximum.
            run(&changed, &original.cards).temperature[vertex]
        };
        let h = 0.002;
        derivative_close(actual, (perturb(h) - perturb(-h)) / (2.0 * h));
    }
    derivative_close(shared_bath, expected_bath);
}
