use super::*;
use crate::resistive_film::FilmChannelReport;
use crate::resistive_film::pressure::PressureDrivenFilm;

fn supply(pressure_pa: f64, liquid_fraction: f64) -> Option<CavitationReservoir> {
    Some(CavitationReservoir {
        pressure_pa,
        liquid_fraction,
    })
}

fn supplied_film(alpha: f64) -> CavitatingFilm {
    let geometry = ResistiveFilm::new(
        vec![FilmCell {
            area_m2: 0.001,
            gap: gap(0.001, 1.0),
        }],
        vec![
            FilmChannel {
                from: 0,
                to: None,
                width_m: 0.0012,
                length_m: 0.1,
                gap: gap(0.001, 0.0),
            };
            2
        ],
        1,
        0.01,
        FilmLimits {
            minimum_cell_gap_m: 1e-8,
            maximum_gap_m: 0.1,
            maximum_pressure_pa: 1e9,
        },
    )
    .unwrap()
    .with_barus_viscosity(alpha)
    .unwrap();
    CavitatingFilm::new(geometry, liquid(), config()).unwrap()
}

#[test]
fn supplied_full_film_matches_native_pressure_flow_and_independent_barus_integral() {
    for alpha in [0.0, 1e-6] {
        for (left, right) in [(300_000.0, 100_000.0), (100_000.0, 300_000.0)] {
            let base = supplied_film(alpha);
            let film = base
                .with_reservoirs(&[supply(left, 1.0), supply(right, 1.0)])
                .unwrap();
            // Rebinding accepts the SAME retained inventory, including a state
            // created before the supplies were attached.
            let old = base.initial_state(&[0.0]).unwrap();
            let before = old.clone();
            let step = film
                .advance(&old, &[0.0], &[0.0; 2], 0.01, || false)
                .unwrap();
            assert_eq!(old, before);
            assert_eq!(step.state.liquid_fraction(), &[1.0]);
            let expected_pressure = if alpha == 0.0 {
                0.5 * (left + right)
            } else {
                let reference = liquid().ambient_pressure_pa;
                reference
                    - (0.5
                        * ((-alpha * (left - reference)).exp()
                            + (-alpha * (right - reference)).exp()))
                    .ln()
                        / alpha
            };
            close(step.state.pressure_pa()[0], expected_pressure);
            close(
                step.generalized_force[0],
                0.001 * (expected_pressure - 100_000.0),
            );
            close(
                step.channel_volume_flow_m3_s[0],
                -step.channel_volume_flow_m3_s[1],
            );
            assert!(step.report.mass_residual_kg.abs() < 1e-16);
            assert!(step.report.maximum_cell_mass_residual_kg < 1e-16);

            let native = PressureDrivenFilm::new(
                base.geometry().clone(),
                &[Some(left - 100_000.0), Some(right - 100_000.0)],
            )
            .unwrap();
            let mut force = [0.0];
            let mut pressure = [0.0];
            let mut channels = [FilmChannelReport::default(); 2];
            native
                .evaluate_into(&[0.0], &[0.0], &mut force, &mut pressure, &mut channels)
                .unwrap();
            close(step.generalized_force[0], force[0]);
            close(step.state.pressure_pa()[0] - 100_000.0, pressure[0]);
            for (actual, expected) in step.channel_volume_flow_m3_s.iter().zip(channels) {
                close(*actual, expected.volume_flow_m3_s);
            }
        }
    }
}

#[test]
fn independent_supply_starvation_refills_then_reforms_without_resetting_mass() {
    let film = supplied_film(0.0)
        .with_reservoirs(&[supply(200_000.0, 1.0), supply(2000.0, 0.25)])
        .unwrap();
    let old = film.initialize(&[0.0], &[0.2]).unwrap();
    let before = old.clone();
    // One pressure-fed inlet and one partially filled Couette inlet. Both
    // contributions are independently known while the cell remains at p_c.
    let first = film
        .advance(&old, &[0.0], &[0.0, -2e-6], 0.01, || false)
        .unwrap();
    close(first.channel_volume_flow_m3_s[0], -1e-10 * 198_000.0);
    close(first.channel_volume_flow_m3_s[1], -2e-6 * 0.25);
    close(first.state.liquid_fraction()[0], 0.403);
    assert_eq!(first.state.pressure_pa(), &[2000.0]);
    close(first.report.liquid_mass_before_kg, 0.0002);
    close(first.report.liquid_mass_after_kg, 0.000403);
    close(first.report.boundary_mass_out_kg, -0.000203);
    assert_eq!(old, before);

    let second = film
        .advance(&first.state, &[0.0], &[0.0, -2e-6], 0.05, || false)
        .unwrap();
    assert_eq!(second.state.liquid_fraction(), &[1.0]);
    close(second.state.pressure_pa()[0], 43_800.0);
    close(
        second.report.liquid_mass_before_kg,
        first.report.liquid_mass_after_kg,
    );
    close(second.report.liquid_mass_after_kg, 0.001);
    close(second.report.boundary_mass_out_kg, -0.000597);

    // Change which reservoir is pressurized without changing the retained
    // state identity. The next full-film through-flow reverses immediately.
    let reversed = film
        .with_reservoirs(&[supply(100_000.0, 1.0), supply(200_000.0, 1.0)])
        .unwrap();
    let last = reversed
        .advance(&second.state, &[0.0], &[0.0; 2], 0.01, || false)
        .unwrap();
    close(last.state.pressure_pa()[0], 150_000.0);
    close(last.channel_volume_flow_m3_s[0], 5e-6);
    close(last.channel_volume_flow_m3_s[1], -5e-6);
    close(
        last.report.liquid_mass_after_kg,
        second.report.liquid_mass_after_kg,
    );
    assert!(last.report.mass_residual_kg.abs() < 1e-16);
    assert!(
        reversed
            .advance(&last.state, &[0.0], &[0.0; 2], 0.01, || true)
            .is_err()
    );
}

#[test]
fn uniform_explicit_supplies_preserve_the_original_trajectory_exactly() {
    for alpha in [0.0, 1e-6] {
        let base = supplied_film(alpha);
        let conditions: Vec<_> = base.reservoir_conditions().collect();
        assert_eq!(conditions, vec![supply(100_000.0, 1.0); 2]);
        let bound = base.with_reservoirs(&conditions).unwrap();
        let mut state = base.initial_state(&[0.0]).unwrap();
        for q in [-0.001, -0.0004, 0.0] {
            let implicit = base
                .advance(&state, &[q], &[0.0; 2], 0.01, || false)
                .unwrap();
            let explicit = bound
                .advance(&state, &[q], &[0.0; 2], 0.01, || false)
                .unwrap();
            assert_eq!(implicit, explicit);
            state = implicit.state;
        }
        let changed = bound.with_viscosity(0.02).unwrap();
        assert_eq!(
            changed.reservoir_conditions().collect::<Vec<_>>(),
            conditions
        );
        assert!(
            changed
                .advance(&state, &[0.0], &[0.0; 2], 0.01, || false)
                .is_ok()
        );
        let changed = bound.with_channel_viscosities(&[0.02, 0.04]).unwrap();
        assert_eq!(
            changed.reservoir_conditions().collect::<Vec<_>>(),
            conditions
        );
        assert!(
            changed
                .advance(&state, &[0.0], &[0.0; 2], 0.01, || false)
                .is_ok()
        );
    }
}

#[test]
fn material_refresh_retains_independent_supplies_and_the_original_inventory() {
    let conditions = [supply(200_000.0, 1.0), supply(100_000.0, 1.0)];
    let film = supplied_film(0.0).with_reservoirs(&conditions).unwrap();
    let old = film.initialize(&[0.0], &[0.2]).unwrap();
    let before = old.clone();
    let mut next_liquid = liquid();
    next_liquid.density_kg_m3 = 800.0;
    next_liquid.cavitation_pressure_pa = 3000.0;
    next_liquid.ambient_pressure_pa = 110_000.0;
    let continuation = film
        .advance_with_liquid(
            &old,
            &[0.0],
            &[0.0; 2],
            0.001,
            next_liquid,
            &[0.01; 2],
            || false,
        )
        .unwrap();
    assert_eq!(
        continuation.film.reservoir_conditions().collect::<Vec<_>>(),
        conditions
    );
    let step = &continuation.step;
    close(step.channel_volume_flow_m3_s[0], -1e-10 * 197_000.0);
    close(step.channel_volume_flow_m3_s[1], -1e-10 * 97_000.0);
    close(step.report.liquid_mass_before_kg, 0.0002);
    close(
        step.report.boundary_mass_out_kg,
        -800.0 * 0.001 * 1e-10 * 294_000.0,
    );
    close(step.report.liquid_mass_after_kg, 0.00022352);
    close(step.state.liquid_fraction()[0], 0.2794);
    assert_eq!(step.state.pressure_pa(), &[3000.0]);
    assert!(step.report.mass_residual_kg.abs() < 1e-16);
    assert_eq!(old, before);

    // Raising the rupture pressure above an authored supply must refuse,
    // rather than replacing that supply with the new ambient reference.
    next_liquid.cavitation_pressure_pa = 150_000.0;
    next_liquid.ambient_pressure_pa = 200_000.0;
    assert!(
        film.advance_with_liquid(
            &old,
            &[0.0],
            &[0.0; 2],
            0.001,
            next_liquid,
            &[0.01; 2],
            || false,
        )
        .is_err()
    );
    assert_eq!(old, before);
}

#[test]
fn reservoir_topology_material_domain_and_barus_barrier_refuse_transactionally() {
    let film = supplied_film(1e-5);
    let old = film.initialize(&[0.0], &[0.5]).unwrap();
    let before = old.clone();
    assert!(film.with_reservoirs(&[]).is_err());
    assert!(
        film.with_reservoirs(&[None, supply(100_000.0, 1.0)])
            .is_err()
    );
    for invalid in [
        supply(f64::NAN, 1.0),
        supply(f64::INFINITY, 1.0),
        supply(1000.0, 1.0),
        supply(100_000.0, 0.9),
        supply(2000.0, -0.1),
        supply(2000.0, 1.1),
        supply(2000.0, f64::NAN),
        supply(2e9, 1.0),
        supply(4_002_000.0, 1.0),
    ] {
        assert!(
            film.with_reservoirs(&[invalid, supply(100_000.0, 1.0)])
                .is_err()
        );
    }
    let chain = two_cells();
    assert!(
        chain
            .with_reservoirs(&[
                supply(2000.0, 1.0),
                supply(2000.0, 1.0),
                supply(2000.0, 1.0),
            ])
            .is_err()
    );
    let bound = chain
        .with_reservoirs(&[supply(2000.0, 0.2), None, supply(2000.0, 0.8)])
        .unwrap();
    assert_eq!(
        bound.reservoir_conditions().collect::<Vec<_>>(),
        vec![supply(2000.0, 0.2), None, supply(2000.0, 0.8)]
    );
    assert_eq!(old, before);
    assert!(
        film.advance(&old, &[0.0], &[0.0; 2], 0.001, || false)
            .is_ok()
    );
}
