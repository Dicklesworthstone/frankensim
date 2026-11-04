//! Imported bodies compose live source contacts, reflected heat and air transport.
use super::*;
use fs_project::{ConductionEnclosure, EnclosureRadiatingSurface, EnclosureViewFactorEvidence};
use fs_project::spec::dims;
use fs_qty::QtyAny;

fn fixture() -> Fixture {
    let mut f = super::fixture(false);
    let setup = f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    let ambient = setup.radiation.take().unwrap();
    setup.enclosure_radiation = Some(ConductionEnclosure {
        spatial_model: None, ambient_temperature: None,
        surfaces: ["cold", "hot"].into_iter().map(|target| EnclosureRadiatingSurface {
            name: format!("{target}-radiator"), target: target.into(),
            card: ambient.surfaces[0].card.clone(), claim: ambient.surfaces[0].claim.clone(),
            query_temperature: QtyAny::new(300.0, dims::TEMPERATURE),
        }).collect(),
        view_factors: vec![vec![0.0, 1.0], vec![1.0, 0.0]],
        evidence: EnclosureViewFactorEvidence::Analytic {
            geometry: "synthetic equal-area exchange; not a visibility reconstruction".into(),
        },
        closure_tolerance: 1e-12, reciprocity_relative_tolerance: 1e-12,
        max_iterations: 256, temperature_tolerance: QtyAny::new(1e-10, dims::TEMPERATURE),
        heat_tolerance: QtyAny::new(1e-7, dims::POWER), relaxation: 0.05,
    });
    f
}

#[test]
fn source_curve_and_surface_contact_share_closed_enclosure_and_original_air_heat() {
    let f = fixture();
    let curve = receipt(&f, 190);
    let surface_pack = temperature_pack_from_value(PropertyValue::Surface2D {
        axes: ["T".into(), "normal_pressure".into()],
        axis_dims: [dims::TEMPERATURE, dims::PRESSURE],
        samples: vec![([200.0, 0.8e6], 0.18), ([450.0, 0.8e6], 0.43),
            ([200.0, 1.2e6], 0.22), ([450.0, 1.2e6], 0.47)],
        triangles: vec![[0, 1, 3], [0, 3, 2]],
        dims: fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
    });
    let selected = surface_pack.2.clone();
    let surface = receipt_with_pack(&f, 191, surface_pack);
    let field = |r: &J| artifact(&f.dir.join("temperature.db"),
        r.str_field("solution_artifact").unwrap()).get("temperature").unwrap().clone();
    for (a, b) in field(&curve).as_array().unwrap().iter().zip(field(&surface).as_array().unwrap()) {
        assert!((a.as_f64().unwrap() - b.as_f64().unwrap()).abs() < 2e-7);
    }
    for (result, picard) in [(&curve, false), (&surface, true)] {
        let contact = result.get("temperature_contact").unwrap();
        let composition = contact.get("enclosure_radiation").unwrap();
        let radiation = result.get("radiation").unwrap();
        let air = result.get("conjugate").unwrap();
        assert_eq!(radiation.str_field("model"), Some("closed-gray-diffuse-enclosure"));
        assert_eq!(contact.str_field("algorithm"), Some(if picard {
            "bounded-contact-picard" } else { "coupled-newton" }));
        assert_eq!(composition.get("physical_endpoint_reassembled"), Some(&J::Bool(true)));
        assert_eq!(number(composition, "thermal_responses"), number(radiation, "solid_solves"));
        assert_eq!(number(composition, "radiation_responses"), number(radiation, "surface_iterations") + 1.0);
        if picard {
            assert!(number(composition, "thermal_responses") > number(composition, "radiation_responses"));
        } else {
            assert_eq!(number(composition, "thermal_responses"), number(composition, "radiation_responses"));
        }
        assert!(number(contact, "total_thermal_responses") >= number(composition, "thermal_responses"));
        assert!(number(contact, "total_nonlinear_iterations") >= number(radiation, "solid_iterations"));
        assert!(number(contact, "total_krylov_iterations") >= number(radiation, "krylov_iterations"));
        assert!(number(contact, "residual_norm_w") <= number(contact, "residual_tolerance_w"));
        assert!(number(contact, "energy_closure_w").abs() <= number(contact, "energy_tolerance_w"));
        assert_eq!(number(result.get("energy").unwrap(), "closure_w"), number(contact, "energy_closure_w"));
        for state in contact.get("face_states").unwrap().as_array().unwrap() {
            let mean = 0.5 * number(state, "side_a_temperature_k") + 0.5 * number(state, "side_b_temperature_k");
            assert_eq!(number(state, "query_temperature_k"), mean);
            assert!((number(state, "resistance_m2_k_w") - 0.001 * mean).abs() < 1e-14);
            if picard { assert_eq!(state.str_field("selected_claim"), Some(selected.as_str())); }
        }
        let patches = radiation.get("surfaces").unwrap().as_array().unwrap();
        assert!(patches.iter().any(|p| number(p, "nonlinear_heat_w") > 0.1));
        assert!(patches.iter().any(|p| number(p, "nonlinear_heat_w") < -0.1));
        assert!(number(radiation, "nonlinear_radiative_out_w").abs() < 1e-7);
        let convection = number(radiation, "convective_out_w");
        assert!((convection - 4.44).abs() < 1e-7, "closed exchange remains internal");
        assert!((number(air, "solid_total_w") - convection).abs() < 1e-7);
        assert!((number(air, "air_total_w") - convection).abs() <= number(air, "balance_tolerance_w"));
        assert!(number(radiation, "decomposition_residual_w").abs() < 1e-7);
    }
}
