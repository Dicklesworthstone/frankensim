//! Whole original-card perturbations include changing reflected irradiation.
//! Surface aliases share source controls; an open cavity has one black bath.

use super::super::adjoint::{derivative_close, gradient};
use super::*;

const ORDINATE: &str = "dtemperature_demissivity_ordinate_k";
const BATH: &str = "dtemperature_dambient_temperature_k_per_k";

fn controls<'a>(receipt: &'a J, material: &MaterialFixture, open: bool) -> &'a J {
    let g = gradient(receipt);
    let radiation = g.get("radiation").unwrap();
    assert_eq!(
        radiation.str_field("method"),
        Some(if open {
            "complete-pointwise-p1-material-open-enclosure-enthalpy-tangent"
        } else {
            "complete-pointwise-p1-material-closed-enclosure-enthalpy-tangent"
        })
    );
    assert_eq!(radiation.str_field("spatial_model"), Some("pointwise"));
    assert_eq!(
        radiation.str_field("emissivity_response"),
        Some("local-temperature")
    );
    let steps = n(g, "state_index");
    close(n(radiation, "control_steps"), steps, 0.0);
    close(n(radiation, "parameter_linear_solves"), steps, 0.0);
    close(
        n(radiation, "parameter_radiosity_transpose_backsolves"),
        steps,
        0.0,
    );
    close(n(radiation, "parameter_thermal_linear_solves"), 0.0, 0.0);
    close(n(radiation, "parameter_krylov_iterations"), 0.0, 0.0);
    assert!(
        fs_blake3::ContentHash::from_hex(radiation.str_field("view_factor_identity").unwrap())
            .is_some()
    );
    let sources = radiation
        .get("emissivity_sources")
        .unwrap()
        .as_array()
        .unwrap();
    assert_eq!(
        sources.len(),
        1,
        "the same source on two surfaces is one original control vector"
    );
    let source = &sources[0];
    assert_eq!(source.str_field("card"), Some(material.card.as_str()));
    assert_eq!(source.str_field("claim"), Some(material.claim.as_str()));
    assert_eq!(
        source.str_field("material_state"),
        Some(material.state.as_str())
    );
    assert_eq!(source.str_field("model"), Some("temperature-curve"));
    assert_eq!(source.str_field("axis"), Some("T"));
    close(n(source, "shared_surface_count"), 2.0, 0.0);
    let knots = source.get("knots").unwrap().as_array().unwrap();
    assert_eq!(knots.len(), SOURCE_KNOTS.len());
    for (index, (row, &(temperature, emissivity))) in knots.iter().zip(&SOURCE_KNOTS).enumerate() {
        close(n(row, "index"), index as f64, 0.0);
        close(n(row, "temperature_k"), temperature, 0.0);
        close(n(row, "emissivity"), emissivity, 0.0);
    }
    let surfaces = radiation.get("surfaces").unwrap().as_array().unwrap();
    assert_eq!(surfaces.len(), 2);
    for (surface, target) in surfaces.iter().zip(["wall-x", "wall-z"]) {
        assert_eq!(surface.str_field("target"), Some(target));
        assert_eq!(surface.str_field("card"), Some(material.card.as_str()));
        assert_eq!(surface.str_field("claim"), Some(material.claim.as_str()));
        close(n(surface, "source_index"), 0.0, 0.0);
        for key in [
            "dtemperature_demissivity_k",
            "reservoir_temperature_k",
            "dtemperature_dreservoir_temperature_k_per_k",
        ] {
            assert!(
                surface.get(key).is_none(),
                "enclosure surfaces must not invent {key}"
            );
        }
    }
    if open {
        let bath = radiation.get("black_environment").unwrap();
        close(
            n(bath, "temperature_k"),
            n(receipt.get("radiation").unwrap(), "ambient_temperature_k"),
            0.0,
        );
    } else {
        assert!(radiation.get("black_environment").is_none());
        assert_eq!(
            radiation.str_field("view_factor_identity"),
            receipt
                .get("radiation")
                .unwrap()
                .str_field("view_factor_identity")
        );
    }
    radiation
}

#[test]
fn g1_original_enclosure_ordinate_adjoint_matches_whole_card_changes_through_reflection() {
    let original = material(&SOURCE_KNOTS, "original-enclosure-specimen");
    for open in [false, true] {
        let spec = fixture(&original, open);
        let (forward, baseline) = evaluate(&spec, &original.cards, None);
        let (receipt, field) = evaluate(&spec, &original.cards, Some("temperature-final-adjoint"));
        assert_eq!(field, baseline);
        for key in ["energy", "radiation"] {
            assert_eq!(receipt.get(key), forward.get(key), "{key}");
        }
        let g = gradient(&receipt);
        assert!(n(g, "maximum_dual_relative_residual") < 1e-9);
        assert!(n(g, "adjoint_krylov_iterations") <= n(g, "max_adjoint_krylov_iterations"));
        let history = receipt.get("transient").unwrap();
        close(n(g, "state_index"), n(history, "fine_steps"), 0.0);
        let vertex = n(g, "active_vertex") as usize;
        close(n(g, "value_k"), field[vertex], 0.0);
        let radiation = controls(&receipt, &original, open);
        let source = &radiation
            .get("emissivity_sources")
            .unwrap()
            .as_array()
            .unwrap()[0];
        let knots = source.get("knots").unwrap().as_array().unwrap();
        for (index, knot) in knots.iter().enumerate() {
            let derivative = n(knot, ORDINATE);
            let probe = |delta: f64| {
                let mut shifted = SOURCE_KNOTS;
                shifted[index].1 += delta;
                let changed = material(&shifted, "original-enclosure-specimen");
                assert_ne!(changed.card, original.card);
                assert_ne!(changed.claim, original.claim);
                // Both physical patches, conductivity binding and original
                // caloric-card identity follow the rebuilt material source.
                // Geometry, view factors, source abscissae, caloric law and
                // initialization remain fixed. No perturbed solve requests a dual.
                evaluate(&fixture(&changed, open), &changed.cards, None).1[vertex]
            };
            let h = 0.0005;
            derivative_close(derivative, (probe(h) - probe(-h)) / (2.0 * h));
            if index == 0 {
                close(derivative, 0.0, 0.0);
            } else {
                assert!(
                    derivative.abs() > 1e-7,
                    "original source ordinate must affect reflected heat"
                );
            }
        }
        if open {
            let derivative = n(radiation.get("black_environment").unwrap(), BATH);
            assert!(derivative > 0.0);
            let probe = |delta: f64| {
                let mut changed = spec.clone();
                changed
                    .cooling
                    .as_mut()
                    .unwrap()
                    .conduction
                    .as_mut()
                    .unwrap()
                    .enclosure_radiation
                    .as_mut()
                    .unwrap()
                    .ambient_temperature
                    .as_mut()
                    .unwrap()
                    .value += delta;
                evaluate(&changed, &original.cards, None).1[vertex]
            };
            let h = 0.02;
            derivative_close(derivative, (probe(h) - probe(-h)) / (2.0 * h));
        }
        // These non-radiation perturbations force the thermal transpose to
        // include the changing source slope AND reflected matrix, beyond the
        // direct original-ordinate parameter contractions above.
        let initial = n(
            g,
            "dtemperature_duniform_initial_specific_enthalpy_k_kg_per_j",
        );
        let pulse: f64 = g
            .get("intervals")
            .unwrap()
            .as_array()
            .unwrap()
            .iter()
            .filter(|row| n(row, "end_s") <= 0.3)
            .map(|row| {
                n(
                    &row.get("regional_power").unwrap().as_array().unwrap()[0],
                    "dtemperature_ddelivered_watts_k_per_w",
                )
            })
            .sum();
        for (coordinate, derivative) in [(0, initial), (1, pulse)] {
            assert!(derivative > 0.0);
            let probe = |delta: f64| {
                let mut changed = spec.clone();
                let policy = time(&mut changed);
                if coordinate == 0 {
                    policy.enthalpy.as_mut().unwrap().materials[0]
                        .initial_specific_enthalpy
                        .value += delta;
                } else {
                    policy.power_schedules[0].steps[0].watts.value += delta;
                }
                evaluate(&changed, &original.cards, None).1[vertex]
            };
            let h = 0.02;
            derivative_close(derivative, (probe(h) - probe(-h)) / (2.0 * h));
        }
    }
}

#[test]
fn g1_open_material_enclosure_initial_peak_keeps_original_sources_and_shared_bath_zero() {
    let original = material(&SOURCE_KNOTS, "original-enclosure-specimen");
    let mut spec = fixture(&original, true);
    spec.cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap()
        .enclosure_radiation
        .as_mut()
        .unwrap()
        .ambient_temperature
        .as_mut()
        .unwrap()
        .value = 280.0;
    let policy = time(&mut spec);
    policy.enthalpy.as_mut().unwrap().materials[0]
        .initial_specific_enthalpy
        .value = 600.0;
    for step in &mut policy.power_schedules[0].steps {
        step.watts.value = 0.0;
    }
    let (receipt, field) = evaluate(
        &spec,
        &original.cards,
        Some("temperature-sampled-peak-adjoint"),
    );
    let g = gradient(&receipt);
    close(n(g, "state_index"), 0.0, 0.0);
    close(n(g, "value_k"), INITIAL_K + 600.0 / CP, 1e-12);
    assert!(
        field
            .iter()
            .all(|&temperature| temperature < n(g, "value_k"))
    );
    close(n(g, "adjoint_krylov_iterations"), 0.0, 0.0);
    close(
        n(
            g,
            "dtemperature_duniform_initial_specific_enthalpy_k_kg_per_j",
        ),
        1.0 / CP,
        1e-14,
    );
    let radiation = controls(&receipt, &original, true);
    for source in radiation
        .get("emissivity_sources")
        .unwrap()
        .as_array()
        .unwrap()
    {
        for knot in source.get("knots").unwrap().as_array().unwrap() {
            close(n(knot, ORDINATE), 0.0, 0.0);
        }
    }
    close(
        n(radiation.get("black_environment").unwrap(), BATH),
        0.0,
        0.0,
    );
    for row in g.get("regional_power").unwrap().as_array().unwrap() {
        close(
            n(row, "dtemperature_duniform_added_watts_k_per_w"),
            0.0,
            0.0,
        );
    }
}
