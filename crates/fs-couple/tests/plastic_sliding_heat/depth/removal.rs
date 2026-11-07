//! Native worn inventory, carried enthalpy and subsequent depth-hardness heat.

use super::*;
use fs_thermal::enthalpy::removal::MaterialRemovalTolerances;

fn export_tolerances() -> MaterialRemovalTolerances {
    MaterialRemovalTolerances {
        mass_kg: 1e-20,
        volume_m3: 1e-23,
    }
}

// Declared uniform prismatic receiver approximation for the one loaded cell.
// Original spatial recession and permanent indentation stay in native contact.
pub(super) fn export_geometry(
    pending_m3: [f64; 2],
    thermal: &DepositionSimulation<'_>,
) -> ([CellGeometry; 2], [ContactPatchUpdate; 1]) {
    let geometry = std::array::from_fn(|side| {
        let mut g = thermal.cells()[&CELLS[side]].geometry;
        if pending_m3[side] > 0.0 {
            g.volume_m3 -= pending_m3[side];
            let length = g.volume_m3 / AREA;
            if side == 0 {
                g.start_m[0] = -length;
            } else {
                g.end_m[0] = length;
            }
            g.surface_area_m2 = 2.0 * AREA + 4.0 * SIDE * length;
        }
        g
    });
    let old = thermal.contact_patch(LINK).unwrap();
    let distances = geometry.map(|g| 0.5 * (g.end_m[0] - g.start_m[0]));
    let replacement = ContactPatch {
        link: ThermalLink {
            conductance_w_k: series_conductance_w_k(AREA, distances, [CONDUCTIVITY; 2]).unwrap(),
            ..old.link
        },
        ..old
    };
    (
        geometry,
        [ContactPatchUpdate {
            id: LINK,
            expected: Some(old),
            replacement: Some(replacement),
        }],
    )
}

fn contact<'a>(
    data: &'a DepthFixture,
    thermal: &DepositionSimulation<'_>,
    physical_limits: PlasticSlidingHeatLimits,
    cx: &Cx<'_>,
) -> DepthSlidingHeatContact<'a> {
    let route = FrictionHeatRoute::from_link(thermal, LINK, [CONDUCTIVITY; 2]).unwrap();
    bind_depth(data, thermal, &route, depth_loaded(data, cx), physical_limits, cx)
}

#[test]
fn g1_depth_export_retains_history_and_carries_signed_enthalpy_into_next_heat_balance() {
    with_cx(|cx, _| {
        for datum in [0.0, -1e6] {
            let mut data = depth_fixture(320.0, cx);
            data.caloric = data.cards.each_ref().map(|card| {
                ResolvedHeatCapacityCurve::try_new_exact(
                    card,
                    HeatCapacitySourceQuery {
                        heat_capacity: &PropertyKey::new(
                            "specific_heat_capacity",
                            SPECIFIC_HEAT_CAPACITY_DIMS,
                        ),
                        density: &PropertyKey::new("density", Density::DIMS),
                        temperature_axis: "T",
                        lower: &point(280.0),
                        upper: &point(320.0),
                        selection: SelectionPolicy::SingleClaimOnly,
                    },
                    SolidLiquidPhase::Solid,
                    datum,
                )
                .unwrap()
            });
            let thermal = network(&data.caloric, T0);
            let bound = contact(&data, &thermal, limits(), cx);
            let step = bound
                .step_candidate(
                    &thermal, 0.1, &motion(0.01), PlasticSlidingHeatBoundary::HeldApproach,
                    controls(), &[0.0; 2], T0, transport, cx,
                )
                .unwrap();
            let before = step.thermal();
            let pending = step.contact().unexported_wear_volume_m3();
            assert!(pending.into_iter().all(|volume| volume > 0.0));
            assert_ne!(pending[0], pending[1]);
            let (geometry, patches) = export_geometry(pending, before);
            let exported = step.contact()
                .export_wear_candidate(
                    before, geometry, &patches, [CONDUCTIVITY; 2], export_tolerances(), cx,
                )
                .unwrap();
            let after = exported.thermal();
            assert_eq!(exported.contact().source(), step.contact().source());
            assert_eq!(exported.contact().unexported_wear_volume_m3(), [0.0; 2]);
            assert_eq!(after.time_s(), before.time_s());
            assert_eq!(after.pending_events(), before.pending_events());
            for side in 0..2 {
                let old = &before.cells()[&CELLS[side]];
                let new = &after.cells()[&CELLS[side]];
                let receipt = &exported.report().cells[side];
                assert_eq!(receipt.requested_volume_m3, pending[side]);
                close(old.mass_kg - new.mass_kg, RHO * pending[side], 1e-10, 1e-21);
                close(
                    receipt.outgoing_enthalpy_j,
                    (old.mass_kg - new.mass_kg) * old.specific_h_j_kg,
                    0.0, 1e-22,
                );
                assert_eq!(new.specific_h_j_kg, old.specific_h_j_kg);
                assert_eq!(after.temperature_k(CELLS[side]), before.temperature_k(CELLS[side]));
                assert!(new.geometry.volume_m3 < old.geometry.volume_m3);
                if datum < 0.0 {
                    assert!(receipt.outgoing_enthalpy_j < 0.0);
                    assert!(receipt.retained_energy_change_j > 0.0);
                }
            }
            close(after.energy_j() + exported.report().outgoing_enthalpy_j, before.energy_j(), 1e-14, 1e-18);
            close(after.mass_kg() + exported.report().removed_mass_kg, before.mass_kg(), 1e-14, 1e-22);
            assert_eq!(exported.report().energy_residual_j, 0.0);

            // A second export has no material to remove, including at negative h.
            let (geometry, patches) = export_geometry([0.0; 2], after);
            let repeated = exported.contact()
                .export_wear_candidate(
                    after, geometry, &patches, [CONDUCTIVITY; 2], export_tolerances(), cx,
                )
                .unwrap();
            assert_eq!(repeated.report().removed_mass_kg, 0.0);
            assert_eq!(repeated.report().outgoing_enthalpy_j, 0.0);
            assert_eq!(repeated.thermal().cells(), after.cells());
            assert_eq!(repeated.contact().source(), exported.contact().source());

            // Independent two-capacity backward-Euler reference uses REMAINING
            // mass and the shortened original conduction distance.
            let power_w = 1e-4;
            let dt = 0.1;
            let capacity = CELLS.map(|id| after.cells()[&id].mass_kg * CP);
            let initial = CELLS.map(|id| after.temperature_k(id).unwrap());
            let gdt = after.links()[&LINK].conductance_w_k * dt;
            let determinant = capacity[0] * capacity[1] + gdt * (capacity[0] + capacity[1]);
            let difference = initial[0] - initial[1];
            let expected = [
                initial[0] + (power_w * dt * (capacity[1] + gdt) - gdt * capacity[1] * difference) / determinant,
                initial[1] + (power_w * dt * gdt + gdt * capacity[0] * difference) / determinant,
            ];
            let next = exported.contact()
                .step_candidate(
                    after, 0.2, &motion(0.0), PlasticSlidingHeatBoundary::HeldApproach,
                    controls(), &[0.0; 2], step.report().heat.interface_temperature_k,
                    |state, temperature| {
                        let mut source = transport(state, temperature)?;
                        source.transport.sources_w.insert(CELLS[0], power_w);
                        Ok::<_, DepositionError>(source)
                    }, cx,
                )
                .unwrap();
            for side in 0..2 {
                close(next.thermal().temperature_k(CELLS[side]).unwrap(), expected[side], 0.0, 2e-6);
                close(
                    next.report().removed_volume_fraction[side],
                    next.report().removed_volume_m3[side] / (AREA * LENGTH),
                    0.0, 1e-15,
                );
            }
            assert_eq!(next.report().frictional_work_j, 0.0);
            assert_eq!(next.contact().unexported_wear_volume_m3(), [0.0; 2]);
        }
    });
}

#[test]
fn g3_depth_export_preserves_prebinding_inventory_and_original_cumulative_limit() {
    with_cx(|cx, _| {
        let data = depth_fixture(320.0, cx);
        let thermal = network(&data.caloric, T0);
        let preworn = depth_loaded(&data, cx)
            .advance(&motion(0.001), 0.1, controls().friction, HardnessSlidingControl::Approach(5e-5), cx)
            .unwrap().candidate;
        let baseline = preworn.wear_state().volumes_m3();
        let route = FrictionHeatRoute::from_link(&thermal, LINK, [CONDUCTIVITY; 2]).unwrap();
        let bounded = PlasticSlidingHeatLimits {
            maximum_removed_volume_fraction: [5e-6; 2],
            ..limits()
        };
        let bound = bind_depth(&data, &thermal, &route, preworn, bounded, cx);
        assert_eq!(bound.unexported_wear_volume_m3(), [0.0; 2]);
        let step = bound
            .step_candidate(
                &thermal, 0.1, &motion(0.001), PlasticSlidingHeatBoundary::HeldApproach,
                controls(), &[0.0; 2], T0, transport, cx,
            )
            .unwrap();
        let cumulative = step.contact().source().wear_state().volumes_m3();
        let pending = step.contact().unexported_wear_volume_m3();
        for side in 0..2 {
            assert_eq!(pending[side], cumulative[side] - baseline[side]);
            assert!(pending[side] > 0.0 && pending[side] < cumulative[side]);
        }
        let (geometry, patches) = export_geometry(pending, step.thermal());
        let exported = step.contact()
            .export_wear_candidate(
                step.thermal(), geometry, &patches, [CONDUCTIVITY; 2], export_tolerances(), cx,
            )
            .unwrap();
        assert_eq!(exported.contact().source().wear_state().volumes_m3(), cumulative);
        for side in 0..2 {
            assert_eq!(exported.report().cells[side].requested_volume_m3, pending[side]);
        }
        let saved_source = exported.contact().source().clone();
        let saved_cells = exported.thermal().cells().clone();
        // Export does not renew the original cumulative 5e-6 wear envelope.
        assert!(matches!(
            exported.contact().step_candidate(
                exported.thermal(), 0.2, &motion(0.005), PlasticSlidingHeatBoundary::HeldApproach,
                controls(), &[0.0; 2], step.report().heat.interface_temperature_k, transport, cx,
            ),
            Err(PlasticSlidingHeatError::Invalid(
                "cumulative removal exceeds original thermal wear envelope"
            ))
        ));
        assert_eq!(exported.contact().source(), &saved_source);
        assert_eq!(exported.thermal().cells(), &saved_cells);
        assert_eq!(exported.contact().unexported_wear_volume_m3(), [0.0; 2]);
    });
}

#[test]
fn g0_depth_export_refuses_invalid_geometry_or_cancellation_without_losing_pending_wear() {
    with_cx(|cx, gate| {
        let data = depth_fixture(320.0, cx);
        let thermal = network(&data.caloric, T0);
        let bound = contact(&data, &thermal, limits(), cx);
        let step = bound
            .step_candidate(
                &thermal, 0.1, &motion(0.002), PlasticSlidingHeatBoundary::HeldApproach,
                controls(), &[0.0; 2], T0, transport, cx,
            )
            .unwrap();
        let before = step.thermal();
        let owner = step.contact();
        let pending = owner.unexported_wear_volume_m3();
        let saved_source = owner.source().clone();
        let saved_cells = before.cells().clone();
        let (geometry, patches) = export_geometry(pending, before);
        assert!(owner.export_wear_candidate(
            before, geometry, &[], [CONDUCTIVITY; 2], export_tolerances(), cx,
        ).is_err());
        let mut shifted = geometry;
        for g in &mut shifted {
            g.start_m[1] += SIDE;
            g.end_m[1] += SIDE;
        }
        assert!(matches!(owner.export_wear_candidate(
            before, shifted, &patches, [CONDUCTIVITY; 2], export_tolerances(), cx,
        ), Err(PlasticSlidingHeatError::Invalid(_))));
        let mut wrong_volume = geometry;
        wrong_volume[1].volume_m3 = before.cells()[&CELLS[1]].geometry.volume_m3;
        assert!(owner.export_wear_candidate(
            before, wrong_volume, &patches, [CONDUCTIVITY; 2], export_tolerances(), cx,
        ).is_err());
        gate.request();
        assert!(owner.export_wear_candidate(
            before, geometry, &patches, [CONDUCTIVITY; 2], export_tolerances(), cx,
        ).is_err());
        assert_eq!(owner.source(), &saved_source);
        assert_eq!(before.cells(), &saved_cells);
        assert_eq!(owner.unexported_wear_volume_m3(), pending);
    });
}
