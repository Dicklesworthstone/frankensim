use super::*;
use crate::resistive_film::{FilmCell, FilmChannel, FilmLimits, GapPort};

fn gap(reference: f64, closure: f64) -> GapPort {
    GapPort {
        reference_m: reference,
        closure: vec![closure],
    }
}
fn liquid() -> CavitationLiquid {
    CavitationLiquid {
        density_kg_m3: 1000.0,
        cavitation_pressure_pa: 2000.0,
        ambient_pressure_pa: 100_000.0,
        ambient_liquid_fraction: 1.0,
    }
}
fn config() -> CavitationSolveConfig {
    CavitationSolveConfig {
        maximum_iterations: 128,
        absolute_mass_tolerance_kg: 1e-14,
        relative_mass_tolerance: 1e-11,
    }
}
fn geometry(closed: bool) -> ResistiveFilm {
    ResistiveFilm::new(
        vec![FilmCell {
            area_m2: 0.001,
            gap: gap(0.001, 1.0),
        }],
        vec![FilmChannel {
            from: 0,
            to: None,
            width_m: 0.0012,
            length_m: 0.1,
            gap: gap(if closed { 0.0 } else { 0.001 }, 0.0),
        }],
        1,
        0.01,
        FilmLimits {
            minimum_cell_gap_m: 1e-8,
            maximum_gap_m: 0.1,
            maximum_pressure_pa: 1e9,
        },
    )
    .unwrap()
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-11 * expected.abs().max(1e-12),
        "{actual} != {expected}"
    );
}

#[test]
fn g1_full_squeeze_is_the_existing_reynolds_pressure_and_force() {
    let original = geometry(false);
    let film = CavitatingFilm::new(original.clone(), liquid(), config()).unwrap();
    let old = film.initial_state(&[0.0]).unwrap();
    let step = film
        .advance(&old, &[0.0001], &[0.0], 0.01, || false)
        .unwrap();
    let mut force = [0.0];
    let mut pressure = [0.0];
    original
        .evaluate_into(&[0.0001], &[0.01], &mut force, &mut pressure)
        .unwrap();
    close(step.state.pressure_pa()[0], 200_000.0);
    close(
        step.state.pressure_pa()[0] - liquid().ambient_pressure_pa,
        pressure[0],
    );
    close(step.generalized_force[0], force[0]);
    assert_eq!(step.state.liquid_fraction(), &[1.0]);
    close(step.report.boundary_mass_out_kg, 0.0001);
    assert!(step.report.mass_residual_kg.abs() < 1e-16);
    assert_eq!(old.configuration(), &[0.0]);
}

#[test]
fn g1_separation_creates_void_without_manufacturing_liquid_then_refills() {
    let film = CavitatingFilm::new(geometry(false), liquid(), config()).unwrap();
    let old = film.initial_state(&[0.0]).unwrap();
    let step = film
        .advance(&old, &[-0.001], &[0.0], 0.01, || false)
        .unwrap();
    // G=1e-10 m^3/(Pa s). Cavitated inflow is G*(100000-2000).
    let inflow = 1e-10 * 98_000.0;
    close(step.channel_volume_flow_m3_s[0], -inflow);
    close(
        step.state.liquid_fraction()[0],
        (1e-6 + 0.01 * inflow) / 2e-6,
    );
    assert_eq!(step.state.pressure_pa(), &[2000.0]);
    close(step.generalized_force[0], -98.0);
    close(
        step.report.liquid_mass_after_kg,
        step.report.liquid_mass_before_kg - step.report.boundary_mass_out_kg,
    );
    let refill = film
        .advance(&step.state, &[0.0], &[0.0], 0.01, || false)
        .unwrap();
    assert_eq!(refill.state.liquid_fraction(), &[1.0]);
    close(refill.state.pressure_pa()[0], 198_000.0);
    assert!(refill.report.mass_residual_kg.abs() < 1e-16);
}

#[test]
fn g1_sealed_expansion_and_reformation_conserve_the_original_liquid() {
    let film = CavitatingFilm::new(geometry(true), liquid(), config()).unwrap();
    let old = film.initial_state(&[0.0]).unwrap();
    let first = film
        .advance(&old, &[-0.001], &[0.0], 0.1, || false)
        .unwrap();
    close(first.state.liquid_fraction()[0], 0.5);
    assert_eq!(first.state.pressure_pa(), &[2000.0]);
    assert_eq!(first.report.boundary_mass_out_kg, 0.0);
    close(
        first.report.liquid_mass_after_kg,
        first.report.liquid_mass_before_kg,
    );
    let second = film
        .advance(&first.state, &[-0.0005], &[0.0], 0.1, || false)
        .unwrap();
    close(second.state.liquid_fraction()[0], 2.0 / 3.0);
    let full = film
        .advance(&second.state, &[0.0], &[0.0], 0.1, || false)
        .unwrap();
    close(full.state.liquid_fraction()[0], 1.0);
    assert!(
        film.advance(&full.state, &[0.0001], &[0.0], 0.1, || false)
            .is_err()
    );
    assert_eq!(full.state.configuration(), &[0.0]);
}

fn two_cells() -> CavitatingFilm {
    let geometry = ResistiveFilm::new(
        vec![
            FilmCell {
                area_m2: 0.001,
                gap: gap(0.001, 0.0)
            };
            2
        ],
        vec![
            FilmChannel {
                from: 0,
                to: None,
                width_m: 0.0012,
                length_m: 0.1,
                gap: gap(0.001, 0.0),
            },
            FilmChannel {
                from: 0,
                to: Some(1),
                width_m: 0.0012,
                length_m: 0.1,
                gap: gap(0.001, 0.0),
            },
            FilmChannel {
                from: 1,
                to: None,
                width_m: 0.0012,
                length_m: 0.1,
                gap: gap(0.001, 0.0),
            },
        ],
        1,
        0.01,
        FilmLimits {
            minimum_cell_gap_m: 1e-8,
            maximum_gap_m: 0.1,
            maximum_pressure_pa: 1e9,
        },
    )
    .unwrap();
    let mut medium = liquid();
    medium.ambient_pressure_pa = medium.cavitation_pressure_pa;
    medium.ambient_liquid_fraction = 0.4;
    CavitatingFilm::new(geometry, medium, config()).unwrap()
}

#[test]
fn g1_upwind_couette_uses_inlet_filling_and_closes_every_face() {
    let film = two_cells();
    let old = film.initialize(&[0.0], &[0.2, 0.1]).unwrap();
    // dt*Q/V = 1. The implicit upwind recurrence is
    // theta0=(0.2+0.4)/2=0.3, theta1=(0.1+theta0)/2=0.2.
    let step = film
        .advance(&old, &[0.0], &[-1e-5, 1e-5, 1e-5], 0.1, || false)
        .unwrap();
    close(step.state.liquid_fraction()[0], 0.3);
    close(step.state.liquid_fraction()[1], 0.2);
    close(step.channel_volume_flow_m3_s[0], -4e-6);
    close(step.channel_volume_flow_m3_s[1], 3e-6);
    close(step.channel_volume_flow_m3_s[2], 2e-6);
    close(step.report.boundary_mass_out_kg, -0.0002);
    assert!(step.report.maximum_cell_mass_residual_kg < 1e-16);
}

#[test]
fn g3_reversing_entrainment_mirrors_the_same_transport_solution() {
    let film = two_cells();
    let forward = film.initialize(&[0.0], &[0.2, 0.1]).unwrap();
    let reverse = film.initialize(&[0.0], &[0.1, 0.2]).unwrap();
    let a = film
        .advance(&forward, &[0.0], &[-1e-5, 1e-5, 1e-5], 0.1, || false)
        .unwrap();
    let b = film
        .advance(&reverse, &[0.0], &[1e-5, -1e-5, -1e-5], 0.1, || false)
        .unwrap();
    close(a.state.liquid_fraction()[0], b.state.liquid_fraction()[1]);
    close(a.state.liquid_fraction()[1], b.state.liquid_fraction()[0]);
    close(a.report.liquid_mass_after_kg, b.report.liquid_mass_after_kg);
    close(a.report.boundary_mass_out_kg, b.report.boundary_mass_out_kg);
}

#[test]
fn g3_starved_content_is_nonnegative_even_at_large_transport_steps() {
    let film = two_cells();
    let old = film.initialize(&[0.0], &[0.0, 0.0]).unwrap();
    let step = film
        .advance(&old, &[0.0], &[-1e-5, 1e-5, 1e-5], 100.0, || false)
        .unwrap();
    for theta in step.state.liquid_fraction() {
        assert!((0.0..=1.0).contains(theta));
    }
    assert!(step.state.liquid_fraction()[0] > 0.399);
    assert!(step.report.mass_residual_kg.abs() < 1e-14);
}

#[test]
fn g0_absolute_pressure_budget_singular_and_cancelled_steps_are_transactional() {
    let film = CavitatingFilm::new(geometry(false), liquid(), config()).unwrap();
    let old = film.initial_state(&[0.0]).unwrap();
    let saved = old.clone();
    let mut polls = 0;
    let result = film.advance(&old, &[-0.001], &[0.0], 0.01, || {
        polls += 1;
        polls >= 3
    });
    assert_eq!(
        result.unwrap_err(),
        FilmError("cavitating-film solve cancelled")
    );
    assert_eq!(old, saved);
    assert!(
        film.advance(&old, &[0.001], &[0.0], 0.01, || false)
            .is_err()
    );
    assert!(
        film.advance(&old, &[0.0], &[f64::NAN], 0.01, || false)
            .is_err()
    );
    let mut invalid = liquid();
    invalid.ambient_pressure_pa = 0.0;
    assert!(CavitatingFilm::new(geometry(false), invalid, config()).is_err());
    invalid = liquid();
    invalid.ambient_liquid_fraction = 0.5;
    assert!(CavitatingFilm::new(geometry(false), invalid, config()).is_err());
    let sealed = CavitatingFilm::new(geometry(true), liquid(), config()).unwrap();
    let sealed_old = sealed.initial_state(&[0.0]).unwrap();
    assert!(
        sealed
            .advance(&sealed_old, &[0.0001], &[0.0], 0.1, || false)
            .is_err()
    );
    assert!(
        sealed
            .advance(&sealed_old, &[-0.001], &[1.0], 0.1, || false)
            .is_err()
    );
    assert_eq!(old, saved);
}

#[test]
fn g0_another_liquid_or_geometry_cannot_reinterpret_retained_mass() {
    let film = CavitatingFilm::new(geometry(false), liquid(), config()).unwrap();
    let old = film.initial_state(&[0.0]).unwrap();
    let mut other_liquid = liquid();
    other_liquid.density_kg_m3 = 500.0;
    let other = CavitatingFilm::new(geometry(false), other_liquid, config()).unwrap();
    assert!(other.advance(&old, &[0.0], &[0.0], 0.1, || false).is_err());
    let mut other_geometry = geometry(false);
    std::sync::Arc::make_mut(&mut other_geometry.cells)[0].area_m2 *= 2.0;
    let other = CavitatingFilm::new(other_geometry, liquid(), config()).unwrap();
    assert!(other.advance(&old, &[0.0], &[0.0], 0.1, || false).is_err());
    // A clone is the SAME immutable physical owner and can continue its state.
    assert!(
        film.clone()
            .advance(&old, &[0.0], &[0.0], 0.1, || false)
            .is_ok()
    );
}

#[test]
fn g5_identical_inputs_replay_the_same_candidate_and_receipts() {
    let film = two_cells();
    let old = film.initialize(&[0.0], &[0.2, 0.1]).unwrap();
    let a = film
        .advance(&old, &[0.0], &[-1e-5, 1e-5, 1e-5], 0.1, || false)
        .unwrap();
    let b = film
        .advance(&old, &[0.0], &[-1e-5, 1e-5, 1e-5], 0.1, || false)
        .unwrap();
    assert_eq!(a, b);
}

#[test]
fn g1_squeezed_full_cell_feeds_adjacent_cavity_with_one_conserved_flux() {
    let geometry = ResistiveFilm::new(
        vec![
            FilmCell {
                area_m2: 0.001,
                gap: gap(0.001, 1.0),
            },
            FilmCell {
                area_m2: 0.001,
                gap: gap(0.001, 0.0),
            },
        ],
        vec![FilmChannel {
            from: 0,
            to: Some(1),
            width_m: 0.0012,
            length_m: 0.1,
            gap: gap(0.001, 0.0),
        }],
        1,
        0.01,
        FilmLimits {
            minimum_cell_gap_m: 1e-8,
            maximum_gap_m: 0.1,
            maximum_pressure_pa: 1e9,
        },
    )
    .unwrap();
    let mut medium = liquid();
    medium.ambient_pressure_pa = medium.cavitation_pressure_pa;
    let film = CavitatingFilm::new(geometry, medium, config()).unwrap();
    let old = film.initialize(&[0.0], &[1.0, 0.2]).unwrap();
    let step = film
        .advance(&old, &[0.0005], &[0.0], 0.01, || false)
        .unwrap();
    // The squeezed cell loses 0.5e-6 m^3; the closed second cell gains exactly
    // that liquid. Its remaining void fixes p1=p_c, hence p0-p1=Q/G.
    close(step.state.liquid_fraction()[0], 1.0);
    close(step.state.liquid_fraction()[1], 0.7);
    close(step.channel_volume_flow_m3_s[0], 5e-5);
    close(
        step.state.pressure_pa()[0] - step.state.pressure_pa()[1],
        500_000.0,
    );
    close(step.generalized_force[0], 500.0);
    close(
        step.report.liquid_mass_after_kg,
        step.report.liquid_mass_before_kg,
    );
    assert_eq!(step.report.boundary_mass_out_kg, 0.0);
    // Shrinking total capacity below the retained liquid requires compression.
    assert!(
        film.advance(&step.state, &[0.0009], &[0.0], 0.01, || false)
            .is_err()
    );
}

#[test]
fn g1_g3_viscosity_continuation_preserves_cavity_and_changes_refill_flux() {
    let original = CavitatingFilm::new(geometry(false), liquid(), config()).unwrap();
    let full = original.initial_state(&[0.0]).unwrap();
    let cavity = original
        .advance(&full, &[-0.001], &[0.0], 0.01, || false)
        .unwrap();
    let saved = cavity.state.clone();
    assert!(saved.liquid_fraction()[0] < 1.0);
    let continued = original.with_viscosity(0.02).unwrap();
    assert_eq!(continued.liquid(), original.liquid());
    assert_eq!(original.geometry().viscosity_pa_s(), Some(0.01));
    let baseline = original
        .advance(&saved, &[-0.001], &[0.0], 0.01, || false)
        .unwrap();
    let slower = continued
        .advance(&saved, &[-0.001], &[0.0], 0.01, || false)
        .unwrap();
    assert_eq!(
        slower.state.pressure_pa(),
        &[liquid().cavitation_pressure_pa]
    );
    close(
        slower.channel_volume_flow_m3_s[0],
        baseline.channel_volume_flow_m3_s[0] / 2.0,
    );
    close(
        slower.report.liquid_mass_before_kg,
        cavity.report.liquid_mass_after_kg,
    );
    close(
        slower.report.liquid_mass_after_kg,
        cavity.report.liquid_mass_after_kg - slower.report.boundary_mass_out_kg,
    );
    close(
        slower.state.liquid_fraction()[0],
        saved.liquid_fraction()[0] + 0.01 * 0.5e-10 * 98_000.0 / 2e-6,
    );
    assert!(slower.report.mass_residual_kg.abs() < 1e-16);
    assert_eq!(cavity.state, saved);
    let reformed = continued
        .advance(&slower.state, &[0.0], &[0.0], 0.01, || false)
        .unwrap();
    assert_eq!(reformed.state.liquid_fraction(), &[1.0]);
    close(
        reformed.report.liquid_mass_before_kg,
        slower.report.liquid_mass_after_kg,
    );
    close(
        reformed.report.liquid_mass_after_kg,
        reformed.report.liquid_mass_before_kg - reformed.report.boundary_mass_out_kg,
    );
}

#[test]
fn g1_g3_viscosity_continuation_scales_squeeze_pressure_without_changing_mass() {
    let original = CavitatingFilm::new(geometry(false), liquid(), config()).unwrap();
    let old = original.initial_state(&[0.0]).unwrap();
    let continued = original.with_viscosity(0.03).unwrap();
    let a = original
        .advance(&old, &[0.0001], &[0.0], 0.01, || false)
        .unwrap();
    let b = continued
        .advance(&old, &[0.0001], &[0.0], 0.01, || false)
        .unwrap();
    let ambient = liquid().ambient_pressure_pa;
    close(
        b.state.pressure_pa()[0] - ambient,
        3.0 * (a.state.pressure_pa()[0] - ambient),
    );
    close(b.generalized_force[0], 3.0 * a.generalized_force[0]);
    close(b.channel_volume_flow_m3_s[0], a.channel_volume_flow_m3_s[0]);
    close(
        b.report.liquid_mass_before_kg,
        a.report.liquid_mass_before_kg,
    );
    close(b.report.liquid_mass_after_kg, a.report.liquid_mass_after_kg);
    assert_eq!(old.configuration(), &[0.0]);
    let sealed = CavitatingFilm::new(geometry(true), liquid(), config()).unwrap();
    let full = sealed.initial_state(&[0.0]).unwrap();
    let cavity = sealed
        .advance(&full, &[-0.001], &[0.0], 0.1, || false)
        .unwrap();
    let next = sealed
        .with_viscosity(0.03)
        .unwrap()
        .advance(&cavity.state, &[-0.0005], &[0.0], 0.1, || false)
        .unwrap();
    close(next.state.liquid_fraction()[0], 2.0 / 3.0);
    close(
        next.report.liquid_mass_after_kg,
        cavity.report.liquid_mass_after_kg,
    );
    assert_eq!(next.report.boundary_mass_out_kg, 0.0);
}

#[test]
fn g0_viscosity_continuation_keeps_state_ownership_and_refuses_invalid_values() {
    let original = CavitatingFilm::new(geometry(false), liquid(), config()).unwrap();
    let state = original.initialize(&[0.0], &[0.4]).unwrap();
    let saved = state.clone();
    for viscosity in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(original.with_viscosity(viscosity).is_err());
    }
    let continued = original.with_viscosity(0.02).unwrap();
    assert!(
        continued
            .advance(&state, &[0.0], &[0.0], 0.001, || false)
            .is_ok()
    );
    let foreign = CavitatingFilm::new(geometry(false), liquid(), config())
        .unwrap()
        .with_viscosity(0.02)
        .unwrap();
    assert!(
        foreign
            .advance(&state, &[0.0], &[0.0], 0.001, || false)
            .is_err()
    );
    assert_eq!(state, saved);
    let first = original
        .advance(&state, &[0.0], &[0.0], 0.001, || false)
        .unwrap();
    let retry = original
        .advance(&state, &[0.0], &[0.0], 0.001, || false)
        .unwrap();
    assert_eq!(first, retry);
}

#[test]
fn g1_g3_channel_viscosities_keep_retained_cavities_and_conserve_separate_refill_fluxes() {
    let graph = geometry(false);
    let mut channels = graph.channels.to_vec();
    channels.push(channels[0].clone());
    // Construct the complete original graph before any state acquires ownership.
    let graph = ResistiveFilm::new(graph.cells.to_vec(), channels, 1, 0.01, graph.limits).unwrap();
    let original = CavitatingFilm::new(graph, liquid(), config()).unwrap();
    let old = original.initialize(&[-0.001], &[0.2]).unwrap();
    let saved = old.clone();
    let current = original.with_channel_viscosities(&[0.005, 0.02]).unwrap();
    let step = current.advance(&old, &[-0.001], &[0.0; 2], 0.01, || false).unwrap();
    let rates = [-2e-10 * 98_000.0, -0.5e-10 * 98_000.0];
    for (actual, expected) in step.channel_volume_flow_m3_s.iter().zip(rates) {
        close(*actual, expected);
    }
    assert_eq!(step.state.pressure_pa(), &[2000.0]);
    close(step.state.liquid_fraction()[0], 0.2 - 0.01 * (rates[0] + rates[1]) / 2e-6);
    close(step.report.liquid_mass_before_kg, 0.0004);
    close(step.report.boundary_mass_out_kg, 1000.0 * 0.01 * (rates[0] + rates[1]));
    close(step.report.liquid_mass_after_kg,
        step.report.liquid_mass_before_kg - step.report.boundary_mass_out_kg);
    assert!(step.report.mass_residual_kg.abs() < 1e-16);
    assert_eq!(old, saved);

    // Once reformed, both different resistances share ONE pressure and their
    // actual fluxes still sum to the volume expelled by the moving wall.
    let full = original.initial_state(&[0.0]).unwrap();
    let squeeze = current.advance(&full, &[0.0001], &[0.0; 2], 0.01, || false).unwrap();
    let volume_rate = 0.001 * 0.0001 / 0.01;
    let pressure = volume_rate / 2.5e-10;
    close(squeeze.state.pressure_pa()[0] - liquid().ambient_pressure_pa, pressure);
    close(squeeze.channel_volume_flow_m3_s[0], 2e-10 * pressure);
    close(squeeze.channel_volume_flow_m3_s[1], 0.5e-10 * pressure);
    close(squeeze.generalized_force[0], 0.001 * pressure);
    assert!(squeeze.report.mass_residual_kg.abs() < 1e-16);

    let uniform = original.with_channel_viscosities(&[0.01; 2]).unwrap();
    let baseline = original.advance(&old, &[-0.001], &[0.0; 2], 0.01, || false).unwrap();
    let repeated = uniform.advance(&old, &[-0.001], &[0.0; 2], 0.01, || false).unwrap();
    assert_eq!(baseline, repeated);
    assert!(original.with_channel_viscosities(&[0.01]).is_err());
    assert!(original.with_channel_viscosities(&[0.01, -0.01]).is_err());
    assert_eq!(old, saved);
}
