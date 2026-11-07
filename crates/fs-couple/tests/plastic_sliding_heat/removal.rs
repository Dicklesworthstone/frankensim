use super::*;
use fs_thermal::enthalpy::removal::MaterialRemovalTolerances;

fn export_tolerances() -> MaterialRemovalTolerances {
    MaterialRemovalTolerances {
        mass_kg: 1e-20,
        volume_m3: 1e-23,
    }
}

// Authored uniform-cell surrogate for this synthetic one-loaded-cell fixture.
// The actual normal recession remains in the native, unchanged contact history.
pub(super) fn export_geometry(
    contact: &PlasticSlidingHeatContact<'_>,
    thermal: &DepositionSimulation<'_>,
) -> ([CellGeometry; 2], [ContactPatchUpdate; 1]) {
    let volumes = contact.unexported_wear_volume_m3();
    let geometry = std::array::from_fn(|side| {
        let mut g = thermal.cells()[&CELLS[side]].geometry;
        if volumes[side] > 0.0 {
            g.volume_m3 -= volumes[side];
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

fn worn<'a>(
    data: &'a Fixture,
    thermal: &DepositionSimulation<'a>,
    cx: &Cx<'_>,
) -> PlasticSlidingHeatCandidate<'a, 'a, [f64; 3]> {
    let route = FrictionHeatRoute::from_link(thermal, LINK, [CONDUCTIVITY; 2]).unwrap();
    PlasticSlidingHeatContact::new(
        thermal,
        data.source.clone(),
        loaded(data, cx),
        [&data.cards[0], &data.cards[1]],
        [&data.caloric[0], &data.caloric[1]],
        &route,
        limits(),
        cx,
    )
    .unwrap()
    .step_candidate(
        thermal,
        0.1,
        &motion(0.01),
        PlasticSlidingHeatBoundary::HeldApproach,
        controls(),
        &[0.0; 2],
        T0,
        transport,
        cx,
    )
    .unwrap()
}

#[test]
fn g1_worn_mass_and_carried_enthalpy_leave_before_the_next_finite_heat_solve() {
    with_cx(|cx, _| {
        let data = fixture(cx);
        let thermal = network(&data.caloric, T0);
        let step = worn(&data, &thermal, cx);
        let before = step.thermal();
        let pending = step.contact().unexported_wear_volume_m3();
        assert!(pending[0] > 0.0 && pending[1] > 0.0);
        assert_ne!(pending[0], pending[1]);
        let (geometry, patches) = export_geometry(step.contact(), before);
        let exported = step
            .contact()
            .export_wear_candidate(
                before,
                geometry,
                &patches,
                [CONDUCTIVITY; 2],
                export_tolerances(),
                cx,
            )
            .unwrap();
        let after = exported.thermal();
        assert_eq!(exported.contact().state(), step.contact().state());
        assert_eq!(exported.contact().source(), step.contact().source());
        assert_eq!(exported.contact().limits(), step.contact().limits());
        assert_eq!(exported.contact().time_s(), before.time_s());
        assert_eq!(after.time_s(), before.time_s());
        assert_eq!(after.pending_events(), before.pending_events());
        assert_eq!(exported.contact().unexported_wear_volume_m3(), [0.0; 2]);
        for side in 0..2 {
            let old = &before.cells()[&CELLS[side]];
            let new = &after.cells()[&CELLS[side]];
            let receipt = &exported.report().cells[side];
            assert_eq!(receipt.cell_id, CELLS[side]);
            assert_eq!(receipt.requested_volume_m3, pending[side]);
            close(old.mass_kg - new.mass_kg, RHO * pending[side], 1e-10, 1e-21);
            close(
                receipt.outgoing_enthalpy_j,
                (old.mass_kg - new.mass_kg) * old.specific_h_j_kg,
                0.0,
                1e-22,
            );
            assert_eq!(new.specific_h_j_kg, old.specific_h_j_kg);
            assert_eq!(after.temperature_k(CELLS[side]), before.temperature_k(CELLS[side]));
            assert!(
                exported.contact().route().distances_m()[side]
                    < step.contact().route().distances_m()[side]
            );
        }
        close(
            after.energy_j() + exported.report().outgoing_enthalpy_j,
            before.energy_j(),
            1e-14,
            1e-18,
        );
        close(
            after.mass_kg() + exported.report().removed_mass_kg,
            before.mass_kg(),
            1e-14,
            1e-22,
        );
        assert_eq!(exported.report().energy_residual_j, 0.0);

        // Independent two-capacity backward-Euler oracle after material export.
        // With zero sliding and default zero plastic heat, only this authored
        // source and the shortened native series conductance heat the bodies.
        let power_w = 1e-4;
        let dt = 0.1;
        let capacity = CELLS.map(|id| after.cells()[&id].mass_kg * CP);
        let initial = CELLS.map(|id| after.temperature_k(id).unwrap());
        let gdt = after.links()[&LINK].conductance_w_k * dt;
        let determinant = capacity[0] * capacity[1] + gdt * (capacity[0] + capacity[1]);
        let difference = initial[0] - initial[1];
        let expected = [
            initial[0]
                + (power_w * dt * (capacity[1] + gdt) - gdt * capacity[1] * difference)
                    / determinant,
            initial[1]
                + (power_w * dt * gdt + gdt * capacity[0] * difference) / determinant,
        ];
        let next = exported
            .contact()
            .step_candidate(
                after,
                0.2,
                &motion(0.0),
                PlasticSlidingHeatBoundary::HeldApproach,
                controls(),
                &[0.0; 2],
                step.report().heat.interface_temperature_k,
                |state, temperature| {
                    let mut source = transport(state, temperature)?;
                    source.transport.sources_w.insert(CELLS[0], power_w);
                    Ok::<_, DepositionError>(source)
                },
                cx,
            )
            .unwrap();
        for side in 0..2 {
            close(next.thermal().temperature_k(CELLS[side]).unwrap(), expected[side], 0.0, 2e-6);
            close(
                next.report().removed_volume_fraction[side],
                next.report().removed_volume_m3[side] / (AREA * LENGTH),
                0.0,
                1e-15,
            );
        }
        assert_eq!(next.report().frictional_work_j, 0.0);
        assert_eq!(next.contact().unexported_wear_volume_m3(), [0.0; 2]);
    });
}

#[test]
fn g3_wear_export_is_not_repeated_and_does_not_remove_prebinding_history() {
    with_cx(|cx, _| {
        let data = fixture(cx);
        let thermal = network(&data.caloric, T0);
        let step = worn(&data, &thermal, cx);
        let (geometry, patches) = export_geometry(step.contact(), step.thermal());
        let exported = step
            .contact()
            .export_wear_candidate(
                step.thermal(), geometry, &patches, [CONDUCTIVITY; 2], export_tolerances(), cx,
            )
            .unwrap();
        let (geometry, patches) = export_geometry(exported.contact(), exported.thermal());
        let repeated = exported
            .contact()
            .export_wear_candidate(
                exported.thermal(), geometry, &patches, [CONDUCTIVITY; 2], export_tolerances(), cx,
            )
            .unwrap();
        assert_eq!(repeated.report().removed_mass_kg, 0.0);
        assert_eq!(repeated.report().outgoing_enthalpy_j, 0.0);
        assert_eq!(repeated.thermal().cells(), exported.thermal().cells());
        assert_eq!(repeated.contact().state(), exported.contact().state());

        let later = repeated
            .contact()
            .step_candidate(
                repeated.thermal(),
                0.2,
                &motion(0.002),
                PlasticSlidingHeatBoundary::HeldApproach,
                controls(),
                &[0.0; 2],
                step.report().heat.interface_temperature_k,
                transport,
                cx,
            )
            .unwrap();
        let old_total = exported.contact().state().wear_state().volumes_m3();
        let total = later.contact().state().wear_state().volumes_m3();
        let pending = later.contact().unexported_wear_volume_m3();
        for side in 0..2 {
            assert_eq!(pending[side], total[side] - old_total[side]);
            assert!(pending[side] > 0.0 && pending[side] < total[side]);
        }
        let (geometry, patches) = export_geometry(later.contact(), later.thermal());
        let later_export = later
            .contact()
            .export_wear_candidate(
                later.thermal(), geometry, &patches, [CONDUCTIVITY; 2], export_tolerances(), cx,
            )
            .unwrap();
        for side in 0..2 {
            assert_eq!(later_export.report().cells[side].requested_volume_m3, pending[side]);
        }

        // A new binding takes the supplied physical inventory as its baseline.
        // It must not remove historical wear that preceded that binding.
        let rebound = PlasticSlidingHeatContact::new(
            later.thermal(),
            later.contact().source().clone(),
            later.contact().state().clone(),
            [&data.cards[0], &data.cards[1]],
            [&data.caloric[0], &data.caloric[1]],
            later.contact().route(),
            limits(),
            cx,
        )
        .unwrap();
        assert_eq!(rebound.unexported_wear_volume_m3(), [0.0; 2]);
    });
}

#[test]
fn g0_wear_export_refuses_stale_geometry_and_preserves_both_histories() {
    with_cx(|cx, gate| {
        let data = fixture(cx);
        let thermal = network(&data.caloric, T0);
        let step = worn(&data, &thermal, cx);
        let contact = step.contact();
        let before = step.thermal();
        let saved_cells = before.cells().clone();
        let saved_state = contact.state().clone();
        let saved_pending = contact.unexported_wear_volume_m3();
        let (geometry, patches) = export_geometry(contact, before);
        assert!(
            contact.export_wear_candidate(
                before, geometry, &[], [CONDUCTIVITY; 2], export_tolerances(), cx,
            ).is_err()
        );
        let mut stale = patches;
        stale[0].expected = None;
        assert!(
            contact.export_wear_candidate(
                before, geometry, &stale, [CONDUCTIVITY; 2], export_tolerances(), cx,
            ).is_err()
        );
        let mut shifted = geometry;
        for g in &mut shifted {
            g.start_m[1] += SIDE;
            g.end_m[1] += SIDE;
        }
        assert!(matches!(
            contact.export_wear_candidate(
                before, shifted, &patches, [CONDUCTIVITY; 2], export_tolerances(), cx,
            ),
            Err(PlasticSlidingHeatError::Invalid(_))
        ));
        let mut rotated = geometry;
        for g in &mut rotated {
            g.start_m.swap(0, 1);
            g.end_m.swap(0, 1);
        }
        assert!(matches!(
            contact.export_wear_candidate(
                before, rotated, &patches, [CONDUCTIVITY; 2], export_tolerances(), cx,
            ),
            Err(PlasticSlidingHeatError::Invalid(_))
        ));
        let mut widened = geometry;
        for g in &mut widened {
            g.start_m[0] *= 0.5;
            g.end_m[0] *= 0.5;
        }
        let mut wider_patches = patches;
        wider_patches[0].replacement.as_mut().unwrap().link.conductance_w_k *= 2.0;
        assert!(matches!(
            contact.export_wear_candidate(
                before, widened, &wider_patches, [CONDUCTIVITY; 2], export_tolerances(), cx,
            ),
            Err(PlasticSlidingHeatError::Invalid(_))
        ));
        gate.request();
        assert!(
            contact.export_wear_candidate(
                before, geometry, &patches, [CONDUCTIVITY; 2], export_tolerances(), cx,
            ).is_err()
        );
        assert_eq!(before.cells(), &saved_cells);
        assert_eq!(contact.state(), &saved_state);
        assert_eq!(contact.unexported_wear_volume_m3(), saved_pending);
    });
}
