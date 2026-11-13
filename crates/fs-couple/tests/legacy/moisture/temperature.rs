//! Prescribed temperature histories retain the original source and inventories.
//! All coefficients are authored verification data, not measured wood evidence.
use super::*;
use fs_blake3::ContentHash;
use fs_couple::moisture::cartesian::{MoistureCartesian, MoistureCartesianGeometry};
use fs_material::state_point::integrate_orthotropic_hygroscopic_expansion;
use fs_qty::Temperature;

#[path = "temperature/radial.rs"]
mod radial;

const COLD: f64 = 293.15;
const HOT: f64 = 313.15;

fn thermal_card(degenerate: bool, hot_moisture_limit: f64) -> MaterialCard {
    let original = card([0.01; 2], 0.3);
    let mut claims = ClaimSet::new();
    for (_, original_claim) in original.claims().claims_ordered() {
        let mut claim = original_claim.clone();
        // Density is deliberately valid ONLY at the original cold datum.
        // Changing the temperature must never re-resolve that density.
        if claim.key.name() != "density" {
            claim.validity = claim.validity.with_quantity(
                "temperature", qty(QuantityKind::AbsoluteTemperature), COLD, HOT,
            );
        }
        if claim.key.name() == fs_material::state_point::HYGROSCOPIC_REFERENCE_MOISTURE_PROPERTY {
            claim.interpolation = InterpolationPolicy::ConstantWithinValidity;
        }
        if claim.key.name() == "moisture_diffusivity" {
            let value = |m: f64, hot: bool| {
                if degenerate { m * if hot { 0.3 } else { 0.1 } }
                else { 0.1 * m + if hot { 0.03 } else { 0.01 } }
            };
            claim.value = PropertyValue::Surface2D {
                axes: ["moisture_content".into(), "temperature".into()],
                axis_dims: [Dims::NONE, Temperature::DIMS],
                samples: vec![
                    ([0.0, COLD], value(0.0, false)),
                    ([0.3, COLD], value(0.3, false)),
                    ([0.0, HOT], value(0.0, true)),
                    ([hot_moisture_limit, HOT], value(hot_moisture_limit, true)),
                ],
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims: MOISTURE_DIFFUSIVITY_DIMS,
            };
            for (name, factor) in [("diffusion_x", 1.0), ("diffusion_y", 2.0), ("diffusion_z", 3.0)] {
                let mut directional = claim.clone();
                directional.key = PropertyKey::new(name, MOISTURE_DIFFUSIVITY_DIMS);
                if let PropertyValue::Surface2D { samples, .. } = &mut directional.value {
                    for (_, value) in samples { *value *= factor; }
                }
                claims.insert_claim(directional).unwrap();
            }
        }
        claims.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(original.id().clone(), claims, original.models().to_vec()).unwrap()
}

fn nonlinear() -> MoistureNonlinearConfig {
    MoistureNonlinearConfig {
        relative_tolerance: 1e-12,
        maximum_iterations: 128,
        relaxation: 1.0,
    }
}

fn grid(material: MaterialCard, initial: Vec<f64>, cx: &Cx<'_>) -> MoistureCartesian {
    MoistureCartesian::new(
        material, point(0.1), ["diffusion_x", "diffusion_y", "diffusion_z"],
        MaterialPropertySelection::SingleClaimOnly,
        MoistureCartesianGeometry {
            lengths_m: [1.0, 0.1, 0.01],
            cells: [initial.len(), 1, 1],
            frame: ContentHash([4; 32]),
            maximum_moisture: 0.3,
            maximum_cells: 256,
        },
        initial, cx,
    ).unwrap()
}

#[test]
fn temperature_event_preserves_history_and_requeries_the_original_diffusivity() {
    with_cx(|cx| {
        let start = slab(thermal_card(false, 0.3), 1.0, vec![0.1], cx);
        let (accepted, _) = start.advance(0.5, &faces(0.2), config(), cx).unwrap();
        let before = accepted.moisture()[0];
        let warm = accepted.at_prescribed_temperature(HOT, cx).unwrap();
        assert_eq!(warm.moisture(), accepted.moisture());
        assert_eq!(warm.time_s(), accepted.time_s());
        assert_eq!(warm.dry_mass_kg(), accepted.dry_mass_kg());
        assert_eq!(warm.water_mass_kg(), accepted.water_mass_kg());
        assert_eq!(warm.wet_mass_kg(), accepted.wet_mass_kg());
        assert_eq!(warm.definition_identity(), accepted.definition_identity());
        assert_eq!(warm.reference().identity(), accepted.reference().identity());
        assert_eq!(warm.temperature_kelvin(), HOT);
        assert_eq!(warm.reference().state_coordinate("temperature"), Some(COLD));
        let original_pin = accepted.reference().property("moisture_diffusivity")
            .unwrap().answer().receipt.selected;
        assert_eq!(warm.transport_context().property("moisture_diffusivity")
            .unwrap().answer().receipt.selected, original_pin);

        let (next, report) = warm.advance(1.0, &faces(0.2), config(), cx).unwrap();
        let diffusion = 0.03 + 0.1 * before;
        near(next.moisture()[0], (before + 4.0 * diffusion * 0.2) / (1.0 + 4.0 * diffusion), 2e-13);
        near(report.transport_start[0].property("moisture_diffusivity").unwrap().value_si(), diffusion, 2e-15);
        assert_eq!(report.transport_end[0].state_coordinate("temperature"), Some(HOT));
        let cold_step = accepted.advance(1.0, &faces(0.2), config(), cx).unwrap().0;
        assert!(next.moisture()[0] > cold_step.moisture()[0]);

        let returned = next.at_prescribed_temperature(COLD, cx).unwrap();
        assert_eq!(returned.moisture(), next.moisture());
        assert_eq!(returned.time_s(), 1.5);
        assert_eq!(returned.reference().identity(), start.reference().identity());
        let source = returned.transport_context().property("moisture_diffusivity").unwrap();
        assert_eq!(source.answer().receipt.selected, original_pin);
        near(source.value_si(), 0.01 + 0.1 * next.moisture()[0], 2e-15);
    });
}

#[test]
fn cartesian_temperature_event_updates_all_principal_sources_and_endpoint_transport() {
    with_cx(|cx| {
        let old = grid(thermal_card(false, 0.3), vec![0.1], cx);
        let warm = old.at_prescribed_temperature(HOT, cx).unwrap();
        assert_eq!(warm.moisture(), old.moisture());
        assert_eq!(warm.water_mass_kg(), old.water_mass_kg());
        assert_eq!(warm.dry_mass_kg(), old.dry_mass_kg());
        assert_eq!(warm.definition_identity(), old.definition_identity());
        assert_eq!(warm.reference().identity(), old.reference().identity());
        assert_eq!(warm.temperature_kelvin(), HOT);
        for (name, factor) in [("diffusion_x", 1.0), ("diffusion_y", 2.0), ("diffusion_z", 3.0)] {
            let response = warm.transport_context().property(name).unwrap();
            near(response.value_si(), 0.04 * factor, 2e-15);
            assert_eq!(response.answer().receipt.selected,
                old.reference().property(name).unwrap().answer().receipt.selected);
        }
        let mut boundaries = std::array::from_fn(|_| MoistureFace::Sealed);
        boundaries[0] = MoistureFace::PrescribedDryBasis(0.2);
        boundaries[1] = MoistureFace::PrescribedDryBasis(0.2);
        let (next, report) = warm.advance(1.0, &boundaries, config(), nonlinear(), cx).unwrap();
        // Backward Euler with D_x=.03+.1M: .4M²+1.04M-.124=0.
        let expected = 0.248 / (1.04 + 1.28_f64.sqrt());
        near(next.moisture()[0], expected, 3e-12);
        near(report.water_balance_residual_kg, 0.0, 2e-12);
        for state in &report.transport_end {
            assert_eq!(state.state_coordinate("temperature"), Some(HOT));
        }
        let returned = next.at_prescribed_temperature(COLD, cx).unwrap();
        assert_eq!(returned.moisture(), next.moisture());
        assert_eq!(returned.time_s(), next.time_s());
        near(returned.transport_context().property("diffusion_z").unwrap().value_si(),
            3.0 * (0.01 + 0.1 * expected), 2e-12);
    });
}

#[test]
fn temperature_events_rebuild_kirchhoff_paths_without_sealing_a_dry_front() {
    with_cx(|cx| {
        let old = MoistureSlab::new_kirchhoff(
            thermal_card(true, 0.3), point(0.1), "moisture_diffusivity",
            MaterialPropertySelection::SingleClaimOnly, geometry(1.0), vec![0.0],
            nonlinear(), cx,
        ).unwrap();
        let warm = old.at_prescribed_temperature(HOT, cx).unwrap();
        assert_eq!(warm.moisture(), &[0.0]);
        let (next, report) = warm.advance(10.0, &faces(0.2), config(), cx).unwrap();
        // Hot D=.3M, Phi=.15M²: 6M²+M-.24=0. Use the stable root.
        near(next.moisture()[0], 0.48 / (1.0 + 6.76_f64.sqrt()), 3e-12);
        near(report.water_balance_residual_kg, 0.0, 3e-12);
        assert!(report.kirchhoff_source.is_some());
        let cold = old.advance(10.0, &faces(0.2), config(), cx).unwrap();
        assert!(next.moisture()[0] > cold.0.moisture()[0]);
        assert_ne!(report.kirchhoff_source, cold.1.kirchhoff_source);
    });
}

#[test]
fn temperature_rebinding_checks_every_cell_and_preserves_the_accepted_state_on_refusal() {
    with_cx(|cx| {
        let old = slab(thermal_card(false, 0.15), 1.0, vec![0.1, 0.25], cx);
        let spatial = grid(thermal_card(false, 0.15), vec![0.1, 0.25], cx);
        // The first cell is supported at HOT; the second lies outside the actual
        // source triangulation even though the enclosing validity box admits it.
        assert!(old.at_prescribed_temperature(HOT, cx).is_err());
        assert!(spatial.at_prescribed_temperature(HOT, cx).is_err());
        for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY, HOT + 1.0] {
            assert!(old.at_prescribed_temperature(invalid, cx).is_err());
            assert!(spatial.at_prescribed_temperature(invalid, cx).is_err());
        }
        assert_eq!(old.moisture(), &[0.1, 0.25]);
        assert_eq!(spatial.moisture(), &[0.1, 0.25]);
        assert_eq!(old.temperature_kelvin(), COLD);
        assert_eq!(spatial.temperature_kelvin(), COLD);
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(spatial.time_s(), 0.0);
    });
}

#[test]
fn source_boundary_and_both_sdirk_paths_use_the_current_temperature() {
    with_cx(|cx| {
        let warm = slab(thermal_card(false, 0.3), 1.0, vec![0.1], cx)
            .at_prescribed_temperature(HOT, cx).unwrap();
        let spatial = grid(thermal_card(false, 0.3), vec![0.1], cx)
            .at_prescribed_temperature(HOT, cx).unwrap();
        let stale = MoistureFace::SourceSorption {
            law: atmosphere::sorption(0.25), relative_humidity: 0.8,
        };
        assert!(warm.advance(0.1, &[stale.clone(), stale.clone()], config(), cx).is_err());
        assert!(warm.advance_sdirk2(0.1, &[stale.clone(), stale.clone()], config(), nonlinear(), cx).is_err());
        let stale_spatial = std::array::from_fn(|_| stale.clone());
        assert!(spatial.advance(0.1, &stale_spatial, config(), nonlinear(), cx).is_err());
        assert!(spatial.advance_sdirk2(0.1, &stale_spatial, config(), nonlinear(), cx).is_err());
        let (_, report) = warm.advance_sdirk2(0.1, &faces(0.2), config(), nonlinear(), cx).unwrap();
        assert_eq!(report.endpoint.transport_end[0].state_coordinate("temperature"), Some(HOT));
        let sealed = std::array::from_fn(|_| MoistureFace::Sealed);
        let (next, _) = spatial.advance_sdirk2(0.1, &sealed, config(), nonlinear(), cx).unwrap();
        assert_eq!(next.moisture(), spatial.moisture());
        assert_eq!(next.temperature_kelvin(), HOT);
    });
}

#[test]
fn finite_atmosphere_matching_uses_current_temperature_not_the_density_datum() {
    with_cx(|cx| {
        let warm = slab(thermal_card(false, 0.3), 1.0, vec![0.1], cx)
            .at_prescribed_temperature(HOT, cx).unwrap();
        let cold_air = atmosphere::atmosphere(0.5, 1.0);
        assert!(cold_air.advance(&warm, 0.1, &[None, None], config(), atmosphere::controls(), cx).is_err());
        let (heated_gas, _) = cold_air.gas().at_temperature(HOT).unwrap();
        let warm_air = fs_couple::moisture::atmosphere::FiniteMoistureAtmosphere::new(heated_gas).unwrap();
        let next = warm_air.advance(&warm, 0.1, &[None, None], config(), atmosphere::controls(), cx).unwrap();
        assert_eq!(next.moisture().moisture(), warm.moisture());
        assert_eq!(next.moisture().temperature_kelvin(), HOT);
        near(next.moisture().water_mass_kg(), warm.water_mass_kg(), 0.0);
    });
}

#[test]
fn temperature_events_keep_the_original_unqualified_thermal_strain_refusal() {
    with_cx(|cx| {
        let old = slab(thermal_card(false, 0.3), 1.0, vec![0.2], cx);
        let warm = old.at_prescribed_temperature(HOT, cx).unwrap();
        // Both source endpoints exist. Their differing non-moisture coordinate
        // still forbids treating thermal deformation as a hygroscopic integral.
        let cold_path = integrate_orthotropic_hygroscopic_expansion(
            old.card(), &old.reference().query(), &point(0.2),
            MaterialPropertySelection::SingleClaimOnly,
        ).unwrap();
        assert!(cold_path.free_principal_strains()[0] > 0.0);
        let unqualified = integrate_orthotropic_hygroscopic_expansion(
            warm.card(), &warm.reference().query(), &warm.transport_context().query(),
            MaterialPropertySelection::SingleClaimOnly,
        );
        assert!(matches!(unqualified, Err(fs_material::state_point::MaterialStatePointError::InvalidDerived { .. })));
        assert_eq!(warm.reference().identity(), old.reference().identity());
        assert_eq!(warm.moisture(), old.moisture());
    });
}
