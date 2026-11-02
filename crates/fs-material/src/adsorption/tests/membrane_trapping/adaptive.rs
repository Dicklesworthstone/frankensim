//! Original-source checks participate in the actual adaptive kinetic trajectory.

use super::*;
use crate::thermochemical_gas::energy::parcel::transfer::species::network::transport::membrane::transient::slab::{
    MembraneDiffusionSolve,
    adaptive::{MembraneDiffusionAdvanceControls, MembraneDiffusionAdvanceError},
    trapping::adaptive::TrappingMembraneAdvanceControls,
};

fn controls() -> TrappingMembraneAdvanceControls {
    TrappingMembraneAdvanceControls {
        temporal: MembraneDiffusionAdvanceControls {
            step: MembraneDiffusionSolve {
                absolute_amount_tolerance_mol: 1e-16,
                relative_amount_tolerance: 1e-11,
            },
            initial_step_s: 2.0,
            minimum_step_s: 1e-8,
            maximum_step_s: 5.0,
            maximum_evaluations: 20_000,
            absolute_pressure_tolerance_pa: 20.0,
            absolute_amount_tolerance_mol: 1e-10,
            relative_tolerance: 0.0,
            amount_balance_tolerance_mol: 1e-13,
        },
        maximum_iterations: 64,
        maximum_backtracks: 32,
    }
}

#[test]
fn original_sources_follow_every_accepted_adaptive_step_without_changing_native_physics() {
    let catalog = pure_catalog();
    let trap = law(&catalog, vacuum_source(PRESSURE));
    let card = membrane_card();
    let context = point(T0, PRESSURE);
    let state = ResolvedTrappingMembraneDiffusion::new(
        mobile(&card, &context),
        &[Some(LangmuirMembraneTrapSource {
            adsorbent: trap.clone(),
            dry_skeleton_density_kg_m3: 4.0,
        })],
        &[0.0; 8],
        [PRESSURE, 0.0],
        || false,
    )
    .unwrap();
    let mut source_work = state.adaptive_workspace().unwrap();
    let candidate = state
        .advance_candidate([PRESSURE, 0.0], 20.0, controls(), &mut source_work, || {
            false
        })
        .unwrap();
    let mut native_work = state.state().adaptive_workspace().unwrap();
    let native = state
        .state()
        .advance_candidate([PRESSURE, 0.0], 20.0, controls(), &mut native_work, || {
            false
        })
        .unwrap();
    assert_eq!(candidate.state().state(), native.state());
    assert_eq!(candidate.report(), native.report());
    assert!(candidate.report().transport.accepted_intervals > 1);
    assert!(candidate.report().trapped_change_mol > 0.0);
    assert_eq!(candidate.state().state().mobile().elapsed_s(), 20.0);
    for support in candidate.state().trap_source_support()[0].as_ref().unwrap() {
        assert_eq!(
            support.source.card_identity(),
            trap.source.material.content_hash()
        );
        assert_eq!(support.affinity_pa_inverse, 2e-5);
        assert_eq!(support.desorption_rate_s_inverse, 0.2);
    }
}

#[test]
fn adaptive_source_refusal_and_late_cancellation_preserve_loading_time_and_receipts() {
    let catalog = pure_catalog();
    let trap = law(&catalog, vacuum_source(PRESSURE));
    let card = membrane_card();
    let context = point(T0, PRESSURE);
    let state = ResolvedTrappingMembraneDiffusion::new(
        mobile(&card, &context),
        &[Some(LangmuirMembraneTrapSource {
            adsorbent: trap,
            dry_skeleton_density_kg_m3: 4.0,
        })],
        &[0.0; 8],
        [PRESSURE, 0.0],
        || false,
    )
    .unwrap();
    let saved = state.state().clone();
    let receipts = state.trap_source_support().to_vec();
    let mut work = state.adaptive_workspace().unwrap();
    assert!(matches!(
        state.advance_candidate([2.0 * PRESSURE, 0.0], 5.0, controls(), &mut work, || false,),
        Err(MembraneDiffusionAdvanceError::Native(
            ResolvedTrappingMembraneError::Adsorption(AdsorptionError::Material(_))
        ))
    ));
    let mut polls = 0;
    let accepted = state
        .advance_candidate([PRESSURE, 0.0], 5.0, controls(), &mut work, || {
            polls += 1;
            false
        })
        .unwrap();
    let mut calls = 0;
    assert!(matches!(
        state.advance_candidate([PRESSURE, 0.0], 5.0, controls(), &mut work, || {
            calls += 1;
            calls == polls
        },),
        Err(MembraneDiffusionAdvanceError::Cancelled)
    ));
    assert_eq!(calls, polls);
    assert_eq!(state.state(), &saved);
    assert_eq!(state.trap_source_support(), receipts);
    let repeated = state
        .advance_candidate([PRESSURE, 0.0], 5.0, controls(), &mut work, || false)
        .unwrap();
    assert_eq!(accepted.state().state(), repeated.state().state());
    assert_eq!(accepted.report(), repeated.report());
}

#[test]
fn supported_coarse_endpoint_cannot_hide_an_unsupported_fine_release_probe() {
    let catalog = pure_catalog();
    let mut source = vacuum_source(3_500.0);
    source.point = point(T0, 0.0);
    let trap = law(&catalog, source);
    let card = membrane_card();
    let context = point(T0, 0.0);
    let layer = MembraneDiffusionLayerSource {
        thickness_m: 1e-4,
        cells: 1,
        material: MembranePermeabilityMaterial {
            species: SpeciesId::new("A").unwrap(),
            material: &card,
            context: &context,
            selection: MaterialPropertySelection::SingleClaimOnly,
        },
    };
    let mobile =
        ResolvedMembraneDiffusion::new(T0, 0.01, &[layer], &[], &[0.0], [0.0; 2], || false)
            .unwrap();
    let state = ResolvedTrappingMembraneDiffusion::new(
        mobile,
        &[Some(LangmuirMembraneTrapSource {
            adsorbent: trap,
            dry_skeleton_density_kg_m3: 4.0,
        })],
        &[(0.01 * 1e-4) * 2.0],
        [0.0; 2],
        || false,
    )
    .unwrap();
    // A long coarse endpoint has already released much of the initial load.
    // Its first half-step still has a higher internal mobile pressure.
    let mut native_work = state.workspace().unwrap();
    let coarse = state
        .step_candidate(
            [0.0; 2],
            1_000.0,
            TrappingMembraneSolve::default(),
            &mut native_work,
            || false,
        )
        .unwrap();
    let fine = state
        .state()
        .step_candidate(
            [0.0; 2],
            500.0,
            TrappingMembraneSolve::default(),
            &mut native_work,
            || false,
        )
        .unwrap();
    assert!(coarse.state().state().mobile().partial_pressures_pa()[0] < 3_500.0);
    assert!(fine.state().mobile().partial_pressures_pa()[0] > 3_500.0);
    let mut c = controls();
    c.temporal.initial_step_s = 1_000.0;
    c.temporal.maximum_step_s = 1_000.0;
    c.temporal.absolute_pressure_tolerance_pa = 1e9;
    c.temporal.absolute_amount_tolerance_mol = 1.0;
    let saved = state.state().clone();
    let receipts = state.trap_source_support().to_vec();
    let mut work = state.adaptive_workspace().unwrap();
    assert!(matches!(
        state.advance_candidate([0.0; 2], 1_000.0, c, &mut work, || false,),
        Err(MembraneDiffusionAdvanceError::Native(
            ResolvedTrappingMembraneError::Adsorption(AdsorptionError::Material(_))
        ))
    ));
    assert_eq!(state.state(), &saved);
    assert_eq!(state.trap_source_support(), receipts);
}

#[test]
fn finite_source_adaptation_preserves_reservoir_mobile_and_trapped_inventories() {
    use crate::thermochemical_gas::energy::parcel::transfer::species::network::transport::membrane::transient::{
        slab::finite::MolarPressureCapacity,
        trapping::ResolvedFiniteTrappingMembraneDiffusion,
    };
    let catalog = pure_catalog();
    let trap = law(&catalog, vacuum_source(PRESSURE));
    let card = membrane_card();
    let context = point(T0, PRESSURE);
    let membrane = ResolvedTrappingMembraneDiffusion::new(
        mobile(&card, &context),
        &[Some(LangmuirMembraneTrapSource {
            adsorbent: trap,
            dry_skeleton_density_kg_m3: 4.0,
        })],
        &[0.0; 8],
        [PRESSURE, 0.0],
        || false,
    )
    .unwrap();
    let capacity = MolarPressureCapacity::try_new(2e-11).unwrap();
    let state = ResolvedFiniteTrappingMembraneDiffusion::new(
        membrane,
        [capacity; 2],
        [PRESSURE, 0.0],
        || false,
    )
    .unwrap();
    let total = state.state().total_amount_mol().unwrap();
    let mut work = state.adaptive_workspace().unwrap();
    let accepted = state
        .advance_candidate(20.0, controls(), &mut work, || false)
        .unwrap();
    let mut native_work = state.state().adaptive_workspace().unwrap();
    let native = state
        .state()
        .advance_candidate(20.0, controls(), &mut native_work, || false)
        .unwrap();
    assert_eq!(accepted.state().state(), native.state());
    assert_eq!(accepted.report(), native.report());
    assert!(accepted.state().reservoir_pressures_pa()[0] < PRESSURE);
    assert!(accepted.state().reservoir_pressures_pa()[1] > 0.0);
    assert!(accepted.report().trapped_after_mol > 0.0);
    assert_eq!(accepted.state().membrane().mobile().elapsed_s(), 20.0);
    assert!((accepted.state().state().total_amount_mol().unwrap() - total).abs() < 1e-15);
    let report = accepted.report();
    let reservoirs = report.transport.finite_reservoirs.as_ref().unwrap();
    let closed = reservoirs.amount_changes_mol[0]
        + report.mobile_change_mol
        + report.trapped_change_mol
        + reservoirs.amount_changes_mol[1];
    assert!(closed.abs() < 1e-15);
    assert!(
        reservoirs.closed_amount_residual_mol.abs()
            <= controls().temporal.amount_balance_tolerance_mol
    );
}
