//! Evolving source kinetics share the original conservative transport operators.

use super::*;
use fs_flux::miscibility::{
    BinaryMiscibilityError as Error, adaptive::BinaryMiscibilityTimeControl,
};
use fs_material::liquid::solution::diffusion::{
    BINARY_LIQUID_DIFFUSIVITY_DIMS as DIMS, MaterialBinaryDiffusivity,
};
use std::sync::Arc;

#[path = "diffusivity/nonisothermal.rs"]
mod nonisothermal;

fn kinetics(axis: &str, knots: Vec<(f64, f64)>, maximum_x: f64) -> MaterialBinaryDiffusivity {
    let property = PropertyKey::new("binary_ms_diffusivity", DIMS);
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: property.clone(),
            value: PropertyValue::Curve {
                abscissa: axis.into(),
                abscissa_dims: if axis == "T" {
                    Temperature::DIMS
                } else {
                    Dimensionless::DIMS
                },
                knots,
                dims: DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "T",
                    quantity(QuantityKind::AbsoluteTemperature),
                    300.0,
                    600.0,
                )
                .with_quantity("x_a", composition(), 0.0, maximum_x)
                .with_quantity(
                    "pressure",
                    quantity(QuantityKind::Pressure),
                    100_000.0,
                    100_000.0,
                ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured binary molecular diffusivity".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic A+B".into(),
            phase: "liquid".into(),
            process: "kinetic fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let reference = QueryPoint::new()
        .with_quantity("T", quantity(QuantityKind::AbsoluteTemperature), 300.0)
        .unwrap()
        .with_quantity("x_a", composition(), 0.5)
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), 100_000.0)
        .unwrap();
    MaterialBinaryDiffusivity::new(
        pair(0.0).species().clone(),
        Arc::new(card),
        property,
        MaterialPropertySelection::SingleClaimOnly,
        reference,
        "T",
        "x_a",
    )
    .unwrap()
}

fn unequal(diffusivity: f64) -> BinaryMiscibilityNetwork {
    BinaryMiscibilityNetwork::new(
        gradient(300.0, 10_000.0, 1e-6, 0.0),
        vec![
            BinaryMiscibilityCell { volume_m3: 1e-4 },
            BinaryMiscibilityCell { volume_m3: 2e-4 },
        ],
        vec![BinaryMiscibilityFace {
            cells: [0, 1],
            area_m2: 1e-2,
            distance_m: 1e-2,
            diffusivity_m2_s: diffusivity,
        }],
        1 << 20,
        || false,
    )
    .unwrap()
}

#[test]
fn constant_source_reproduces_authored_diffusion_and_complete_material_transport() {
    let bare = network(model(), 1e-5);
    let initial = bare
        .state_from_mole_fractions(vec![[0.4, 0.6], [0.6, 0.4]], || false)
        .unwrap();
    let source = kinetics("T", vec![(300.0, 1e-5), (600.0, 1e-5)], 1.0);
    let (net, bound) = bare
        .bind_diffusivity_source(&initial, source.clone(), || false)
        .unwrap()
        .into_parts();
    assert_eq!(net.diffusivity_source(), Some(&source));
    assert_eq!(bound.component_moles(), initial.component_moles());
    let expected = bare
        .step_candidate(&initial, 0.1, BinaryMiscibilitySolve::default(), || false)
        .unwrap();
    let actual = net
        .step_candidate(&bound, 0.1, BinaryMiscibilitySolve::default(), || false)
        .unwrap();
    assert_eq!(
        actual.state().component_moles(),
        expected.state().component_moles()
    );
    assert_eq!(
        actual.state().mole_fractions(),
        expected.state().mole_fractions()
    );
    assert_eq!(actual.state().elapsed_s(), expected.state().elapsed_s());
    assert_eq!(actual.step(), expected.step());

    let flow = bare
        .prepare_advection(
            vec![BinaryMiscibilityAdvectionFace {
                forward_m3_s: 1e-5,
                reverse_m3_s: 1e-5,
            }],
            0.0,
            || false,
        )
        .unwrap();
    let rebound = net.rebind_advection(&flow, || false).unwrap();
    let expected = bare
        .transport_candidate(
            &initial,
            &flow,
            0.1,
            BinaryMiscibilityAdvectionControls::default(),
            BinaryMiscibilitySolve::default(),
            || false,
        )
        .unwrap();
    let card = card(1.0, 590.0);
    let liquid = liquid(&card);
    let liquid_source = BinaryLiquidSquareGradientLaw::try_new(&liquid, model(), 1e-12).unwrap();
    let field = net
        .resolve_liquid_field(&bound, &liquid_source, 2, || false)
        .unwrap();
    let actual = field
        .transport_candidate(
            &rebound,
            0.1,
            BinaryMiscibilityAdvectionControls::default(),
            BinaryMiscibilitySolve::default(),
            || false,
        )
        .unwrap();
    assert_eq!(
        actual.field().state().component_moles(),
        expected.state().component_moles()
    );
    assert_eq!(
        actual.field().state().mole_fractions(),
        expected.state().mole_fractions()
    );
    assert_eq!(actual.step(), expected.step());
    assert_eq!(field.state(), &bound);
}

#[test]
fn temperature_cycles_replay_original_diffusivity_and_rebind_changed_flow_snapshots() {
    let bare = unequal(0.0);
    let initial = bare
        .state_from_mole_fractions(vec![[0.2, 0.8], [0.8, 0.2]], || false)
        .unwrap();
    let source = kinetics("T", vec![(300.0, 1e-5), (600.0, 4e-5)], 1.0);
    let (net, initial) = bare
        .bind_diffusivity_source(&initial, source.clone(), || false)
        .unwrap()
        .into_parts();
    let (hot, hot_state, _) = net
        .temperature_candidate(&initial, Temperature::new(600.0), || false)
        .unwrap()
        .into_parts();
    near(hot.faces().next().unwrap().diffusivity_m2_s, 4e-5);
    assert_eq!(hot.diffusivity_source(), Some(&source));
    assert_eq!(hot_state.component_moles(), initial.component_moles());
    let cold_step = net
        .step_candidate(&initial, 0.5, BinaryMiscibilitySolve::default(), || false)
        .unwrap();
    let hot_step = hot
        .step_candidate(&hot_state, 0.5, BinaryMiscibilitySolve::default(), || false)
        .unwrap();
    assert!(
        hot_step.step().transferred_a_mol[0].abs() > cold_step.step().transferred_a_mol[0].abs()
    );
    let cycled = hot
        .temperature_candidate(&hot_state, Temperature::new(300.0), || false)
        .unwrap();
    assert_eq!(cycled.network(), &net);
    assert_eq!(cycled.state(), &initial);

    let flow = net
        .prepare_advection(
            vec![BinaryMiscibilityAdvectionFace {
                forward_m3_s: 1e-5,
                reverse_m3_s: 1e-5,
            }],
            0.0,
            || false,
        )
        .unwrap();
    let rebound = hot.rebind_advection(&flow, || false).unwrap();
    let transported = hot
        .transport_candidate(
            &hot_state,
            &rebound,
            0.1,
            BinaryMiscibilityAdvectionControls::default(),
            BinaryMiscibilitySolve::default(),
            || false,
        )
        .unwrap();
    assert_eq!(transported.state().elapsed_s(), 0.1);
    assert!(
        transported
            .step()
            .component_residual_mol
            .iter()
            .all(|r| r.abs() < 1e-12)
    );
    assert_eq!(hot.diffusivity_source(), Some(&source));
}

#[test]
fn unequal_cell_composition_refresh_matches_independent_stage_kinetics_and_adaptive_fine_steps() {
    let bare = unequal(0.0);
    let initial = bare
        .state_from_mole_fractions(vec![[0.2, 0.8], [0.8, 0.2]], || false)
        .unwrap();
    // D(x)=1e-5*(1+4*x). Unequal capacities ensure the arithmetic face x
    // changes under conservative diffusion, exposing a frozen-source defect.
    let source = kinetics("x_a", vec![(0.0, 1e-5), (1.0, 5e-5)], 1.0);
    let (net, initial) = bare
        .bind_diffusivity_source(&initial, source, || false)
        .unwrap()
        .into_parts();
    let totals = [0, 1].map(|s| initial.component_moles().iter().map(|n| n[s]).sum::<f64>());
    let mut state = initial.clone();
    for _ in 0..4 {
        let mean = 0.5 * (state.mole_fractions()[0][0] + state.mole_fractions()[1][0]);
        let oracle = unequal(1e-5 * (1.0 + 4.0 * mean));
        let seed = oracle
            .state_from_component_moles(state.component_moles().to_vec(), || false)
            .unwrap();
        let expected = oracle
            .step_candidate(&seed, 0.5, BinaryMiscibilitySolve::default(), || false)
            .unwrap();
        let actual = net
            .step_candidate(&state, 0.5, BinaryMiscibilitySolve::default(), || false)
            .unwrap();
        for (&a, &b) in actual
            .state()
            .component_moles()
            .iter()
            .flatten()
            .zip(expected.state().component_moles().iter().flatten())
        {
            near(a, b);
        }
        assert!(actual.step().free_energy_change_j <= 0.0);
        assert!(actual.step().convex_split_dissipation_j >= 0.0);
        for s in 0..2 {
            near(
                actual.state().component_moles().iter().map(|n| n[s]).sum(),
                totals[s],
            );
        }
        state = actual.into_parts().0;
    }
    assert!(0.5 * (state.mole_fractions()[0][0] + state.mole_fractions()[1][0]) > 0.501);
    near(net.faces().next().unwrap().diffusivity_m2_s, 3e-5); // immutable binding snapshot
    let half = net
        .step_candidate(&initial, 0.25, BinaryMiscibilitySolve::default(), || false)
        .unwrap();
    let fine = net
        .step_candidate(
            half.state(),
            0.25,
            BinaryMiscibilitySolve::default(),
            || false,
        )
        .unwrap();
    let adaptive = net
        .advance_candidate(
            &initial,
            0.5,
            BinaryMiscibilityTimeControl {
                initial_step_s: 0.5,
                minimum_step_s: 0.5,
                maximum_step_s: 0.5,
                absolute_mole_fraction_tolerance: 1.0,
                relative_mole_fraction_tolerance: 0.0,
                maximum_attempts: 1,
            },
            BinaryMiscibilitySolve::default(),
            || false,
        )
        .unwrap();
    assert_eq!(adaptive.advance().accepted_intervals, 1);
    assert_eq!(adaptive.state(), fine.state());
}

#[test]
fn final_source_domain_negative_endpoint_and_late_cancellation_leave_original_state_intact() {
    let bare = network(model(), 1e-5);
    let initial = bare
        .state_from_mole_fractions(vec![[0.4, 0.6], [0.6, 0.4]], || false)
        .unwrap();
    let numerical = bare
        .step_candidate(&initial, 0.1, BinaryMiscibilitySolve::default(), || false)
        .unwrap();
    assert!(numerical.state().mole_fractions()[1][0] > 0.600001);
    for source in [
        kinetics("x_a", vec![(0.0, 1e-5), (1.0, 1e-5)], 0.6),
        kinetics(
            "x_a",
            vec![(0.0, 1e-5), (0.6, 1e-5), (0.600001, -1e-5), (1.0, -1e-5)],
            1.0,
        ),
    ] {
        let (net, state) = bare
            .bind_diffusivity_source(&initial, source, || false)
            .unwrap()
            .into_parts();
        let before = state.clone();
        assert!(matches!(
            net.step_candidate(&state, 0.1, BinaryMiscibilitySolve::default(), || false),
            Err(Error::DiffusivitySource(_))
        ));
        assert_eq!(state, before);
        assert_eq!(state.component_moles(), initial.component_moles());
    }
    let (net, state) = bare
        .bind_diffusivity_source(
            &initial,
            kinetics("T", vec![(300.0, 1e-5), (600.0, 1e-5)], 1.0),
            || false,
        )
        .unwrap()
        .into_parts();
    let before = state.clone();
    let mut polls = 0;
    net.step_candidate(&state, 0.1, BinaryMiscibilitySolve::default(), || {
        polls += 1;
        false
    })
    .unwrap();
    let mut count = 0;
    assert!(matches!(
        net.step_candidate(&state, 0.1, BinaryMiscibilitySolve::default(), || {
            count += 1;
            count == polls
        }),
        Err(Error::Cancelled)
    ));
    assert_eq!(state, before);
}
