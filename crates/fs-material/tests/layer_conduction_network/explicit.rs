use super::*;
use fs_material::thermal::conduction::network::LayerConductionObserver;

fn local_layer<'a>(
    card: &'a MaterialCard,
    caloric: &'a EquilibriumEnthalpyPhaseCurve,
    context: &'a QueryPoint,
    range: [f64; 2],
) -> ResolvedConductionLayer<'a> {
    ResolvedConductionLayer::new(
        caloric,
        ThermalMaterialPoint {
            material: card,
            context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        &PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        SolidLiquidPhase::Solid,
        range,
        || Ok(()),
    )
    .unwrap()
}

#[test]
fn g1_disjoint_contact_sources_advance_explicitly_with_signed_heat_and_entropy() {
    let card = source();
    let caloric = curve(&card);
    let context = QueryPoint::new().with("T", 300.0).unwrap();
    let hot = local_layer(&card, &caloric, &context, [500.0, 800.0]);
    let cold = local_layer(&card, &caloric, &context, [250.0, 450.0]);
    // Independent k(T)=2+0.01T integrals: hot 600->520 K gives 608 W/m;
    // cold 300->400 K gives 550 W/m. j=30400 W/m² and a 120 K contact jump.
    let q = 30_400.0 * 0.004;
    let dt = 1.0 / q;
    for reversed in [false, true] {
        let mut old = initial(&caloric);
        if reversed {
            let patch = old.contact_patch(7).unwrap();
            let mut replacement = patch;
            std::mem::swap(&mut replacement.link.a, &mut replacement.link.b);
            old = old
                .replace_contact_patches(
                    &[ContactPatchUpdate {
                        id: 7,
                        expected: Some(patch),
                        replacement: Some(replacement),
                    }],
                    || Ok(()),
                )
                .unwrap();
        }
        let mut layers = [&hot, &cold];
        let mut distances_m = [0.02, 550.0 / 30_400.0];
        if reversed {
            layers.reverse();
            distances_m.reverse();
        }
        let prepared = PreparedLayerConduction::new(
            &old,
            vec![LayerLinkBinding {
                link: 7,
                layers,
                distances_m,
                contact_area_resistance_m2_k_w: 120.0 / 30_400.0,
            }],
            1,
            || Ok(()),
        )
        .unwrap();
        let observer = LayerConductionObserver::new(&prepared, &old, config(), || Ok(())).unwrap();
        let (next, used, energy, entropy) = old
            .advance_observed_with_entropy(
                &observer,
                dt,
                StepConfig {
                    max_dt_s: dt,
                    energy_tolerance_j: 1e-8,
                },
                || Ok(()),
            )
            .unwrap();
        assert_eq!(energy.substeps, 1);
        assert!((next.temperature_k(CellId(1)).unwrap() - 599.0).abs() < 1e-8);
        assert!((next.temperature_k(CellId(2)).unwrap() - 301.0).abs() < 1e-8);
        assert_eq!(energy.external_heat_j, 0.0);
        assert!(energy.residual_j.abs() < 1e-8);
        // The diagnostic remains the actual starting-state flux, not the
        // endpoint revalidation used to admit the accepted state.
        let expected = if reversed { -q } else { q };
        assert!((used.last_samples()[0].applied_power_a_to_b_w - expected).abs() < 1e-8);
        assert!((entropy.conduction_generation_j_k - (1.0 / 300.0 - 1.0 / 600.0)).abs() < 1e-10);
        assert!(entropy.thermal_storage_change_j_k > 0.0);
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(next.links(), old.links());
    }
}

#[test]
fn g0_supported_bulk_states_do_not_hide_a_final_contact_face_source_exit() {
    let card = source();
    let caloric = curve(&card);
    let context = QueryPoint::new().with("T", 300.0).unwrap();
    let hot = local_layer(&card, &caloric, &context, [500.0, 800.0]);
    let cold = local_layer(&card, &caloric, &context, [250.0, 450.0]);
    let original = initial(&caloric);
    let q = 30_400.0 * 0.004;
    let dt = 1.0 / q;
    // One explicit interval would reach [599,440] K, both admitted BULK
    // temperatures. The cold contact face would exceed its 450 K support.
    let mut events: Vec<_> = original
        .cells()
        .values()
        .cloned()
        .map(|mut cell| {
            if cell.id == CellId(2) {
                cell.boundary.source_w = 139.0 * q;
            }
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(cell),
            }
        })
        .collect();
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(original.links()[&7]),
    });
    let old = DepositionSimulation::new(
        events,
        &[(MaterialId(1), &caloric)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 10,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: dt,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let prepared = PreparedLayerConduction::new(
        &old,
        vec![LayerLinkBinding {
            link: 7,
            layers: [&hot, &cold],
            distances_m: [0.02, 550.0 / 30_400.0],
            contact_area_resistance_m2_k_w: 120.0 / 30_400.0,
        }],
        1,
        || Ok(()),
    )
    .unwrap();
    let observer = LayerConductionObserver::new(&prepared, &old, config(), || Ok(())).unwrap();
    let samples = observer.last_samples().to_vec();
    assert!(hot.temperature_at_enthalpy(39_900.0).is_ok());
    assert!(cold.temperature_at_enthalpy(24_000.0).is_ok());
    assert!(matches!(
        old.advance_observed(
            &observer,
            dt,
            StepConfig {
                max_dt_s: dt,
                energy_tolerance_j: 1e-8
            },
            || Ok(()),
        ),
        Err(LayerSourceError::Thermal(DepositionError::MaterialDomain))
    ));
    assert_eq!(old.time_s(), 0.0);
    assert_eq!(old.temperature_k(CellId(2)).unwrap(), 300.0);
    assert_eq!(observer.last_samples(), samples);
}
