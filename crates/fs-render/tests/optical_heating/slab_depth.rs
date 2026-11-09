//! G1/G3/G4: depth-localized glass heat on the original conservative clock.

use super::*;
use fs_blake3::ContentHash;
use fs_render::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};
use fs_render::thermal_surface::heating::slab_depth::{
    SlabDepthIlluminationSource, SlabDepthPowerResponse,
};
use fs_render::thermal_surface::heating::{
    SurfaceHeatingBudget, SurfaceHeatingError, SurfaceHeatingReceiver, SurfaceHeatingTolerance,
    SurfaceOpticalHeating,
};
use fs_render::volumes::heating::illumination::schedule::{
    SpectralPowerSchedule, SpectralPowerSegment,
};

fn glass(alpha: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(1.5).unwrap(),
        BeerLambertAbsorption::try_constant(alpha).unwrap(),
        GlassProvenance::Custom,
    )
}
fn response(alpha: f64, edges: &[f64], reverse: bool) -> SlabDepthPowerResponse<1> {
    SlabDepthPowerResponse::try_new(
        &glass(alpha),
        [550.0],
        edges,
        [1.0],
        [1.0],
        1.0,
        reverse,
        [0.5],
        ContentHash([244; 32]),
        8,
        || Ok(()),
    )
    .unwrap()
}
fn source(response: SlabDepthPowerResponse<1>, ids: &[u64]) -> SlabDepthIlluminationSource<1> {
    SlabDepthIlluminationSource {
        source_id: 1,
        depth_receivers: ids
            .iter()
            .map(|&id| SurfaceHeatingReceiver {
                thermal_cell: CellId(id),
                material_id: MaterialId(7),
                temperature_bounds_k: (300.0, 600.0),
            })
            .collect(),
        response,
        power: SpectralPowerSchedule::try_new(
            &[SpectralPowerSegment {
                start_s: 0.0,
                end_s: 1.0,
                start_power_w: [20.0],
                end_power_w: [20.0],
            }],
            1,
        )
        .unwrap(),
    }
}
fn budget() -> SurfaceHeatingBudget {
    SurfaceHeatingBudget {
        max_sources: 4,
        max_total_power_segments: 4,
        max_source_evaluations: 1000,
    }
}
fn tolerance() -> SurfaceHeatingTolerance {
    SurfaceHeatingTolerance {
        power_roundoff_w: 1e-11,
    }
}
fn heating(
    native: &DepositionSimulation<'_>,
    sources: &[SlabDepthIlluminationSource<1>],
) -> SurfaceOpticalHeating<1> {
    SurfaceOpticalHeating::try_new_slab_depths(
        sources,
        native,
        budget(),
        16,
        tolerance(),
        || Ok(()),
    )
    .unwrap()
}

#[test]
fn g1_skin_and_interior_heat_use_exact_bounce_integrals_and_one_optical_budget() {
    let alpha: f64 = 60.0;
    let d: f64 = 0.02;
    let edges = [0.0, 0.002, 0.01, d];
    let native = thermal(
        vec![
            cell(11, 1.0, 1.0, 0.0),
            cell(12, 2.0, 1.0, 0.0),
            cell(13, 3.0, 1.0, 0.0),
        ],
        &SENSIBLE,
        vec![],
    );
    let compiled = response(alpha, &edges, false);
    let input = heating(&native, &[source(compiled.clone(), &[11, 12, 13])]);
    let (heated, accepted, audit) = native
        .advance_observed(&input, 0.5, config(0.11), || Ok(()))
        .unwrap();
    // Independent incoherent ray sum in a plate with n=1.5 (R=0.04).
    let r: f64 = 0.04;
    let a = (-alpha * d).exp();
    let entry = (1.0 - r) / (1.0 - r * r * a * a);
    let mut total = 0.0;
    for (i, interval) in edges.windows(2).enumerate() {
        let [z0, z1] = [interval[0], interval[1]];
        let fraction = entry
            * ((-alpha * z0).exp() - (-alpha * z1).exp()
                + r * ((-alpha * (2.0 * d - z1)).exp() - (-alpha * (2.0 * d - z0)).exp()));
        let expected = 10.0 * fraction;
        let id = CellId(11 + i as u64);
        close(
            heated.cells()[&id].mass_kg * heated.cells()[&id].specific_h_j_kg,
            expected,
        );
        total += expected;
    }
    close(accepted.energy().incident_j[0], 10.0);
    close(accepted.energy().absorbed_j[0], total);
    close(accepted.energy().thermal_input_j, total);
    close(audit.external_heat_j, total);
    close(heated.energy_j() - native.energy_j(), total);
    close(
        accepted.energy().reflected_j[0],
        10.0 * compiled.total().reflected_fraction()[0],
    );
    close(
        accepted.energy().transmitted_j[0],
        10.0 * compiled.total().transmitted_fraction()[0],
    );
    assert!(accepted.energy().residual_j[0].abs() < 1e-13);
    assert_eq!(
        accepted.response_evaluations(),
        0,
        "fixed depth optics must not be rebuilt in time loop"
    );
}

#[test]
fn g3_reverse_illumination_moves_heat_to_the_authored_rear_cells() {
    let edges = [0.0, 0.005, 0.01, 0.015, 0.02];
    let forward = response(200.0, &edges, false);
    let backward = response(200.0, &edges, true);
    assert_eq!(
        forward.total().absorbed_fraction(),
        backward.total().absorbed_fraction()
    );
    assert_ne!(forward.identity(), backward.identity());
    for i in 0..4 {
        close(
            forward.depth_absorbed_fractions()[i][0],
            backward.depth_absorbed_fractions()[3 - i][0],
        );
    }
    assert!(
        forward.depth_absorbed_fractions()[0][0] > 10.0 * forward.depth_absorbed_fractions()[3][0]
    );
    assert_eq!(forward.depth_edges_m(), &edges);
}

#[test]
fn g3_reverse_illumination_keeps_heat_in_an_ultrathin_authored_cell() {
    let alpha: f64 = 1e-20;
    let thickness: f64 = 1e10;
    let width: f64 = 1e-10;
    assert_eq!(thickness - width, thickness);
    let compiled = response(alpha, &[0.0, width, thickness], true);
    let r: f64 = 0.04;
    let a = (-alpha * thickness).exp();
    let entry = (1.0 - r) / (1.0 - r * r * a * a);
    let expected_fraction = entry * a * (1.0 + r) * -(-alpha * width).exp_m1();
    let fraction = compiled.depth_absorbed_fractions()[0][0];
    assert!(fraction > 0.0);
    assert!((fraction / expected_fraction - 1.0).abs() < 2e-14);

    let native = thermal(
        vec![cell(11, 1.0, 1.0, 0.0), cell(12, 1.0, 1.0, 0.0)],
        &SENSIBLE,
        vec![],
    );
    let input = heating(&native, &[source(compiled, &[11, 12])]);
    let (heated, _, _) = native
        .advance_observed(&input, 1.0, config(0.25), || Ok(()))
        .unwrap();
    let heat = heated.cells()[&CellId(11)].specific_h_j_kg;
    assert!(heat > 0.0);
    assert!((heat / (20.0 * expected_fraction) - 1.0).abs() < 2e-14);
}

#[test]
fn g3_shared_receiver_and_scheduled_pulse_conserve_the_same_total_heat() {
    let native = thermal(vec![cell(11, 2.0, 1.0, 0.0)], &SENSIBLE, vec![]);
    let compiled = response(50.0, &[0.0, 0.002, 0.015, 0.02], false);
    let mut illumination = source(compiled.clone(), &[11, 11, 11]);
    illumination.power = SpectralPowerSchedule::try_new(
        &[SpectralPowerSegment {
            start_s: 0.2,
            end_s: 0.7,
            start_power_w: [0.0],
            end_power_w: [20.0],
        }],
        1,
    )
    .unwrap();
    let input = heating(&native, &[illumination]);
    let (heated, accepted, audit) = native
        .advance_observed(&input, 1.5, config(0.6), || Ok(()))
        .unwrap();
    let expected = 5.0 * compiled.total().absorbed_fraction()[0];
    close(accepted.energy().incident_j[0], 5.0);
    close(accepted.energy().thermal_input_j, expected);
    close(heated.energy_j(), expected);
    close(audit.external_heat_j, expected);
    assert_eq!(heated.cells().len(), 1);
}

#[test]
fn g1_weak_absorption_reaches_each_original_enthalpy_without_a_residual_floor() {
    let alpha = 1e-24;
    let native = thermal(
        vec![cell(11, 1.0, 1.0, 0.0), cell(12, 2.0, 1.0, 0.0)],
        &SENSIBLE,
        vec![],
    );
    let input = heating(
        &native,
        &[source(
            response(alpha, &[0.0, 0.005, 0.02], false),
            &[11, 12],
        )],
    );
    let (heated, accepted, _) = native
        .advance_observed(&input, 1.0, config(0.25), || Ok(()))
        .unwrap();
    // With equal exterior media at normal incidence, the zero-loss limiting
    // internal bounce intensity is exactly one; local loss is alpha*dz.
    for (id, mass, length) in [(11, 1.0, 0.005), (12, 2.0, 0.015)] {
        let expected = 20.0 * alpha * length;
        let actual = mass * heated.cells()[&CellId(id)].specific_h_j_kg;
        assert!(actual > 0.0);
        assert!((actual / expected - 1.0).abs() < 2e-14);
    }
    assert!(accepted.energy().absorbed_j[0] > 0.0);
    assert_eq!((-alpha * 0.02_f64).exp(), 1.0);
}

#[test]
fn g0_g4_partition_source_binding_budget_and_cancellation_refuse_atomically() {
    for edges in [
        &[][..],
        &[0.0][..],
        &[0.001, 0.02][..],
        &[0.0, 0.01, 0.01][..],
        &[0.0, f64::INFINITY][..],
    ] {
        assert!(
            SlabDepthPowerResponse::try_new(
                &glass(50.0),
                [550.0],
                edges,
                [1.0],
                [1.0],
                1.0,
                false,
                [0.5],
                ContentHash([244; 32]),
                8,
                || Ok(()),
            )
            .is_err()
        );
    }
    assert!(matches!(
        SlabDepthPowerResponse::try_new(
            &glass(50.0),
            [550.0],
            &[0.0, 0.01, 0.02],
            [1.0],
            [1.0],
            1.0,
            false,
            [0.5],
            ContentHash([244; 32]),
            1,
            || Ok(()),
        ),
        Err(SurfaceHeatingError::Budget)
    ));
    let native = thermal(
        vec![cell(11, 1.0, 1.0, 0.0), cell(12, 1.0, 1.0, 0.0)],
        &SENSIBLE,
        vec![],
    );
    let valid = source(response(50.0, &[0.0, 0.01, 0.02], false), &[11, 12]);
    let mut missing = valid.clone();
    missing.depth_receivers.pop();
    assert!(
        SurfaceOpticalHeating::try_new_slab_depths(
            &[missing],
            &native,
            budget(),
            8,
            tolerance(),
            || Ok(())
        )
        .is_err()
    );
    assert!(matches!(
        SurfaceOpticalHeating::try_new_slab_depths(
            core::slice::from_ref(&valid),
            &native,
            budget(),
            1,
            tolerance(),
            || Ok(()),
        ),
        Err(SurfaceHeatingError::Budget)
    ));
    let mut wrong = valid.clone();
    wrong.depth_receivers[0].material_id = MaterialId(99);
    assert!(
        SurfaceOpticalHeating::try_new_slab_depths(
            &[wrong],
            &native,
            budget(),
            8,
            tolerance(),
            || Ok(())
        )
        .is_err()
    );
    let input = heating(&native, &[valid]);
    let mut calls = 0;
    assert!(matches!(
        native.advance_observed(&input, 1.0, config(0.25), || {
            calls += 1;
            if calls > 12 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(SurfaceHeatingError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(native.time_s(), 0.0);
    assert_eq!(native.energy_j(), 0.0);
    assert_eq!(input.time_s(), 0.0);
    assert_eq!(input.energy().incident_j, [0.0]);
}
