//! G1/G3: original temperature-dependent hardness changes an existing dent.

use super::*;
use fs_contact::normal_patch::plastic_source::temperature::SourcePlasticTemperatureError;
use fs_contact::normal_patch::{
    FiniteGapPlasticContinuationKind, FiniteGapPlasticDrive, FiniteGapPlasticLoadControls,
    FiniteGapPlasticLoadError,
};
use fs_matdb::ClaimSet;
use fs_qty::Temperature;

const COLD: f64 = 293.15;
const HOT: f64 = 313.15;

fn thermal_card(
    name: &str,
    cold_hardness: f64,
    hot_hardness: f64,
    hot_young: f64,
) -> (MaterialCard, PropertyKey) {
    let (base, key) = material(name, cold_hardness, "temperature-continuation");
    let mut claims = ClaimSet::new();
    for id in base.claims().observation_ids() {
        claims
            .register_observation(base.claims().observation(id).unwrap().clone())
            .unwrap();
    }
    for (_, original) in base.claims().claims_ordered() {
        let mut claim = original.clone();
        let values = match claim.key.name() {
            "indentation-hardness" => Some((cold_hardness, hot_hardness)),
            "young_modulus" => Some((2.0e7, hot_young)),
            _ => None,
        };
        if let Some((cold, hot)) = values {
            claim.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(COLD, cold), (HOT, hot)],
                dims: Pressure::DIMS,
            };
            claim.validity = ValidityDomain::unconstrained().with("T", COLD, HOT);
            claim.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(claim).unwrap();
    }
    (
        MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap(),
        key,
    )
}

fn thermal_source(
    cold_hardness: [f64; 2],
    hot_hardness: [f64; 2],
    hot_young: [f64; 2],
    envelope: PlasticContactEnvelope,
    cx: &Cx<'_>,
) -> ([MaterialCard; 2], SourcePlasticContact) {
    let a = thermal_card("A", cold_hardness[0], hot_hardness[0], hot_young[0]);
    let b = thermal_card("B", cold_hardness[1], hot_hardness[1], hot_young[1]);
    let resolve = |(card, key): &(MaterialCard, PropertyKey)| {
        let state = resolve_material_state_point(
            card,
            &common::query(),
            &[
                ScalarPropertyRequirement::try_with_key(key, ScalarAdmissibility::StrictlyPositive)
                    .unwrap(),
            ],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        IndentationHardnessStatePoint::try_from_resolved(state, key.clone()).unwrap()
    };
    let hardness = [resolve(&a), resolve(&b)];
    let normal = common::bind(&a.0, &b.0, None, 0.0).normal_state().clone();
    let law = build(normal, hardness, envelope, cx).unwrap();
    ([a.0, b.0], law)
}

fn self_compliance(law: &SourcePlasticContact) -> f64 {
    // Independent exact integral over one square constant-pressure cell.
    let side = law.model().config().grid.cell_width_m;
    4.0 * side * (1.0 + 2.0_f64.sqrt()).ln()
        / (core::f64::consts::PI * law.normal_source().elastic().reduced_modulus_pa())
}

fn thermal_load_controls() -> FiniteGapPlasticLoadControls {
    FiniteGapPlasticLoadControls {
        minimum_approach_m: -1.0e-4,
        maximum_approach_m: 1.0e-4,
        absolute_force_tolerance_n: 1.0e-14,
        relative_force_tolerance: 1.0e-12,
        approach_tolerance_m: 1.0e-18,
        maximum_equilibrium_evaluations: 128,
    }
}

#[test]
fn g1_maintained_load_thermal_modulus_change_retains_dent_and_separates_work() {
    with_cx(|cx, _| {
        let (cards, law) = thermal_source(
            [1.0e6, 2.0e6],
            [0.75e6, 2.0e6],
            [1.0e7, 2.0e7],
            envelope(),
            cx,
        );
        let dented = law.advance(&law.initial_state(), 5.0e-5, cx).unwrap();
        let force = 0.005;
        let loaded = law
            .advance_at_load(&dented.candidate, force, thermal_load_controls(), cx)
            .unwrap();
        let before = loaded.candidate;
        let saved = before.clone();
        let old = before.contact_state();
        let old_energy = loaded.receipt.normal.reversible_energy_j;
        assert!(old.plastic_depth_m()[12] > 0.0);

        let unchanged = law
            .refresh_at_temperatures_at_load(
                &before,
                cards.each_ref(),
                [COLD; 2],
                COLD,
                force,
                thermal_load_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(unchanged.source(), &law);
        assert_eq!(unchanged.candidate(), &before);
        assert_eq!(unchanged.load().equilibrium_evaluations, 0);
        assert_eq!(unchanged.receipt().normal_work_j, 0.0);
        assert_eq!(unchanged.receipt().parameter_work_j, 0.0);

        let hot = law
            .refresh_at_temperatures_at_load(
                &before,
                cards.each_ref(),
                [HOT, COLD],
                COLD,
                force,
                thermal_load_controls(),
                cx,
            )
            .unwrap();
        let next = hot.candidate().contact_state();
        let account = hot.receipt();
        assert_eq!(
            account.kind,
            FiniteGapPlasticContinuationKind::MaterialRefreshAtLoad
        );
        assert_eq!(hot.original_source_identity(), law.identity());
        assert_eq!(account.predecessor_state_identity, old.identity());
        assert_eq!(hot.load().predecessor_state_identity, old.identity());
        assert_eq!(hot.load().response_identity, account.response.identity);
        common::close(hot.load().achieved_force_n, force);
        common::close(next.pressure_pa()[12], force / 1.0e-8);
        common::close(
            next.approach_m() - old.approach_m(),
            0.5 * self_compliance(&law) * force / 1.0e-8,
        );
        common::close(account.elastic_energy_j, 1.5 * old_energy);
        common::close(account.elastic_energy_change_j, 0.5 * old_energy);
        common::close(account.normal_work_j, old_energy);
        common::close(account.parameter_work_j, -0.5 * old_energy);
        assert_eq!(account.response.plastic_work_j, 0.0);
        assert_eq!(next.plastic_depth_m(), old.plastic_depth_m());
        assert_eq!(next.recession_m(), old.recession_m());
        assert_eq!(
            next.cumulative_plastic_work_j(),
            old.cumulative_plastic_work_j()
        );
        common::close(hot.source().hardness_sources()[0].hardness_pa(), 0.75e6);
        assert_eq!(
            hot.source().hardness_sources()[1],
            law.hardness_sources()[1]
        );
        for side in 0..2 {
            assert_eq!(
                hot.source().hardness_sources()[side].key(),
                law.hardness_sources()[side].key()
            );
            assert_eq!(
                hot.source().hardness_sources()[side]
                    .property()
                    .answer()
                    .receipt
                    .selected,
                law.hardness_sources()[side]
                    .property()
                    .answer()
                    .receipt
                    .selected,
            );
        }
        assert!(law.advance(hot.candidate(), next.approach_m(), cx).is_err());
        let cooled = hot
            .source()
            .refresh_at_temperatures_at_load(
                hot.candidate(),
                cards.each_ref(),
                [COLD; 2],
                COLD,
                force,
                thermal_load_controls(),
                cx,
            )
            .unwrap();
        common::close(
            cooled.candidate().contact_state().approach_m(),
            old.approach_m(),
        );
        common::close(cooled.receipt().elastic_energy_j, old_energy);
        assert_eq!(
            cooled.candidate().contact_state().plastic_depth_m(),
            old.plastic_depth_m()
        );
        let (source, state, _, load) = cooled.into_parts();
        common::close(load.achieved_force_n, force);
        assert!(
            source
                .advance_at_load(&state, force, thermal_load_controls(), cx)
                .is_ok()
        );
        assert_eq!(before, saved);
    });
}

#[test]
fn g3_direct_thermal_load_solve_avoids_a_fictitious_held_approach_dent() {
    with_cx(|cx, _| {
        // This synthetic modulus/hardness path deliberately makes a held-
        // approach refresh yield although the maintained-force endpoint is elastic.
        let (cards, law) =
            thermal_source([1.0e6, 2.0e6], [0.75e6, 1.5e6], [2.0e8; 2], envelope(), cx);
        let force = 0.005;
        let before = law
            .advance_at_load(&law.initial_state(), force, thermal_load_controls(), cx)
            .unwrap()
            .candidate;
        let held = law
            .refresh_at_temperatures(&before, cards.each_ref(), [HOT; 2], COLD, cx)
            .unwrap();
        assert!(held.receipt().response.plastic_work_j > 0.0);
        assert!(held.candidate().contact_state().plastic_depth_m()[12] > 0.0);

        let direct = law
            .refresh_at_temperatures_at_load(
                &before,
                cards.each_ref(),
                [HOT; 2],
                COLD,
                force,
                thermal_load_controls(),
                cx,
            )
            .unwrap();
        common::close(direct.load().achieved_force_n, force);
        common::close(
            direct.candidate().contact_state().approach_m(),
            0.1 * before.contact_state().approach_m(),
        );
        assert_eq!(
            direct.candidate().contact_state().plastic_depth_m(),
            &[0.0; 25]
        );
        assert_eq!(direct.receipt().response.plastic_work_j, 0.0);
        assert!(direct.receipt().normal_work_j < 0.0);
        assert!(direct.receipt().parameter_work_j > 0.0);
    });
}

#[test]
fn g0_thermal_load_refusals_preserve_original_sources_and_history() {
    with_cx(|cx, gate| {
        let (cards, law) =
            thermal_source([1.0e6, 2.0e6], [0.25e6, 2.0e6], [2.0e7; 2], envelope(), cx);
        let force = 0.005;
        let before = law
            .advance_at_load(&law.initial_state(), force, thermal_load_controls(), cx)
            .unwrap()
            .candidate;
        let saved = before.clone();
        assert!(matches!(
            law.refresh_at_temperatures_at_load(
                &before,
                cards.each_ref(),
                [HOT, COLD],
                COLD,
                force,
                thermal_load_controls(),
                cx,
            ),
            Err(SourcePlasticTemperatureError::Load(
                FiniteGapPlasticLoadError::TargetOutsideApproachBounds { .. }
            ))
        ));
        assert!(matches!(
            law.refresh_at_temperatures_at_load(
                &before,
                cards.each_ref(),
                [COLD; 2],
                COLD,
                0.004,
                thermal_load_controls(),
                cx,
            ),
            Err(SourcePlasticTemperatureError::Load(
                FiniteGapPlasticLoadError::PredecessorLoadMismatch { .. }
            ))
        ));
        let outside = FiniteGapPlasticLoadControls {
            maximum_approach_m: 1.01e-4,
            ..thermal_load_controls()
        };
        assert!(matches!(
            law.refresh_at_temperatures_at_load(
                &before,
                cards.each_ref(),
                [COLD; 2],
                COLD,
                force,
                outside,
                cx,
            ),
            Err(SourcePlasticTemperatureError::Contact(
                SourcePlasticContactError::OutsideApplicability { .. }
            ))
        ));
        let revised =
            MaterialCard::supersede(&cards[0], cards[0].claims().clone(), vec![]).unwrap();
        assert!(
            law.refresh_at_temperatures_at_load(
                &before,
                [&revised, &cards[1]],
                [COLD; 2],
                COLD,
                force,
                thermal_load_controls(),
                cx,
            )
            .is_err()
        );
        assert!(
            law.refresh_at_temperatures_at_load(
                &before,
                [&cards[1], &cards[0]],
                [COLD; 2],
                COLD,
                force,
                thermal_load_controls(),
                cx,
            )
            .is_err()
        );
        assert!(
            law.refresh_at_temperatures_at_load(
                &before,
                cards.each_ref(),
                [330.0, COLD],
                COLD,
                force,
                thermal_load_controls(),
                cx,
            )
            .is_err()
        );
        gate.request();
        assert!(
            law.refresh_at_temperatures_at_load(
                &before,
                cards.each_ref(),
                [COLD; 2],
                COLD,
                force,
                thermal_load_controls(),
                cx,
            )
            .is_err()
        );
        assert_eq!(before, saved);
    });
}

#[test]
fn g1_thermal_softening_adds_indentation_and_cooling_cannot_heal_it() {
    with_cx(|cx, _| {
        let (cards, law) =
            thermal_source([1.0e6, 2.0e6], [0.5e6, 2.0e6], [2.0e7; 2], envelope(), cx);
        let loaded = law.advance(&law.initial_state(), 5.0e-5, cx).unwrap();
        let before = loaded.candidate.clone();
        let c = self_compliance(&law);
        let area = 1.0e-8;
        let hot = law
            .refresh_at_temperatures(&before, [&cards[0], &cards[1]], [HOT, COLD], COLD, cx)
            .unwrap();
        let account = hot.receipt();
        common::close(account.response.normal.normal_force_n, 0.5e6 * area);
        common::close(account.response.plastic_increment_m[12], c * 0.5e6);
        common::close(account.response.plastic_work_j, area * c * 0.5e6 * 0.5e6);
        common::close(
            hot.candidate().contact_state().plastic_depth_m()[12],
            5.0e-5 - c * 0.5e6,
        );
        common::close(
            account.parameter_work_j,
            account.response.normal.reversible_energy_j - loaded.receipt.normal.reversible_energy_j
                + account.response.plastic_work_j,
        );
        assert!(account.parameter_work_j < 0.0);
        assert_eq!(before, loaded.candidate);
        assert_eq!(hot.original_source_identity(), law.identity());
        assert_eq!(
            hot.source().hardness_sources()[0].key(),
            law.hardness_sources()[0].key()
        );
        assert_eq!(
            hot.source().hardness_sources()[0]
                .property()
                .answer()
                .receipt
                .selected,
            law.hardness_sources()[0]
                .property()
                .answer()
                .receipt
                .selected
        );
        assert_eq!(
            hot.source().hardness_sources()[1],
            law.hardness_sources()[1]
        );

        let cooled = hot
            .source()
            .refresh_at_temperatures(hot.candidate(), [&cards[0], &cards[1]], [COLD; 2], COLD, cx)
            .unwrap();
        assert_eq!(
            cooled.candidate().contact_state().plastic_depth_m(),
            hot.candidate().contact_state().plastic_depth_m()
        );
        assert_eq!(cooled.receipt().response.plastic_work_j, 0.0);
        common::close(cooled.receipt().response.normal.normal_force_n, 0.005);
        assert_eq!(
            cooled
                .candidate()
                .contact_state()
                .cumulative_plastic_work_j(),
            hot.candidate().contact_state().cumulative_plastic_work_j()
        );
    });
}

#[test]
fn g3_separate_modulus_temperatures_change_storage_without_resetting_history() {
    with_cx(|cx, _| {
        let (cards, law) = thermal_source([1.0e8; 2], [1.0e8; 2], [1.0e7, 2.0e7], envelope(), cx);
        let loaded = law.advance(&law.initial_state(), 1.0e-6, cx).unwrap();
        let unchanged = law
            .refresh_at_temperatures(
                &loaded.candidate,
                [&cards[0], &cards[1]],
                [COLD; 2],
                COLD,
                cx,
            )
            .unwrap();
        assert_eq!(unchanged.source(), &law);
        assert_eq!(unchanged.candidate(), &loaded.candidate);
        assert_eq!(unchanged.receipt().parameter_work_j, 0.0);

        let hot = law
            .refresh_at_temperatures(
                &loaded.candidate,
                [&cards[0], &cards[1]],
                [HOT, COLD],
                COLD,
                cx,
            )
            .unwrap();
        let ratio = 2.0 / 3.0;
        common::close(
            hot.receipt().response.normal.normal_force_n,
            ratio * loaded.receipt.normal.normal_force_n,
        );
        common::close(
            hot.receipt().response.normal.reversible_energy_j,
            ratio * loaded.receipt.normal.reversible_energy_j,
        );
        assert_eq!(hot.receipt().response.plastic_work_j, 0.0);
        assert_eq!(
            hot.candidate().contact_state().plastic_depth_m(),
            loaded.candidate.contact_state().plastic_depth_m()
        );
        common::close(
            hot.receipt().parameter_work_j,
            (ratio - 1.0) * loaded.receipt.normal.reversible_energy_j,
        );
        assert!(law.advance(hot.candidate(), 1.0e-6, cx).is_err());
        let (source, state, _) = hot.into_parts();
        assert!(source.advance(&state, 1.1e-6, cx).is_ok());
    });
}

#[test]
fn g0_thermal_source_coverage_and_late_depth_refusals_leave_originals_unchanged() {
    with_cx(|cx, gate| {
        let mut bounds = envelope();
        bounds.maximum_plastic_depth_m = 4.2e-5;
        let (cards, law) = thermal_source([1.0e6, 2.0e6], [0.5e6, 2.0e6], [2.0e7; 2], bounds, cx);
        let loaded = law.advance(&law.initial_state(), 5.0e-5, cx).unwrap();
        let saved = loaded.candidate.clone();
        assert!(matches!(
            law.refresh_at_temperatures(&saved, [&cards[0], &cards[1]], [HOT, COLD], COLD, cx,),
            Err(SourcePlasticTemperatureError::Contact(
                SourcePlasticContactError::OutsideApplicability {
                    quantity: "permanent relative indentation [m]",
                    ..
                }
            ))
        ));
        for temperature in [0.0, f64::NAN, 330.0] {
            assert!(
                law.refresh_at_temperatures(
                    &saved,
                    [&cards[0], &cards[1]],
                    [temperature, COLD],
                    COLD,
                    cx,
                )
                .is_err()
            );
        }
        let changed =
            MaterialCard::supersede(&cards[0], cards[0].claims().clone(), vec![]).unwrap();
        assert!(
            law.refresh_at_temperatures(&saved, [&changed, &cards[1]], [COLD; 2], COLD, cx,)
                .is_err()
        );
        assert!(
            law.refresh_at_temperatures(&saved, [&cards[1], &cards[0]], [COLD; 2], COLD, cx,)
                .is_err()
        );
        gate.request();
        assert!(
            law.refresh_at_temperatures(&saved, [&cards[0], &cards[1]], [COLD; 2], COLD, cx,)
                .is_err()
        );
        assert_eq!(saved, loaded.candidate);
    });
}

#[test]
fn g1_softening_at_a_held_device_releases_spring_energy_and_cooling_retains_the_dent() {
    with_cx(|cx, _| {
        let (cards, law) =
            thermal_source([1.0e6, 2.0e6], [0.5e6, 2.0e6], [2.0e7; 2], envelope(), cx);
        let c = self_compliance(&law);
        let h = 1.0e6;
        let area = 1.0e-8;
        let drive = FiniteGapPlasticDrive {
            stiffness_n_per_m: area / c,
            position_m: 3.0 * c * h,
        };
        let loaded = law
            .advance_with_drive(&law.initial_state(), drive, thermal_load_controls(), cx)
            .unwrap();
        let before = loaded.candidate;
        let saved = before.clone();
        common::close(before.contact_state().plastic_depth_m()[12], c * h);
        let unchanged = law
            .refresh_at_temperatures_with_drive(
                &before,
                cards.each_ref(),
                [COLD; 2],
                COLD,
                drive,
                thermal_load_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(unchanged.source(), &law);
        assert_eq!(unchanged.candidate(), &before);
        assert_eq!(unchanged.drive().equilibrium_evaluations, 0);
        assert_eq!(unchanged.receipt().normal_work_j, 0.0);
        assert_eq!(unchanged.receipt().parameter_work_j, 0.0);

        let hot = law
            .refresh_at_temperatures_with_drive(
                &before,
                cards.each_ref(),
                [HOT, COLD],
                COLD,
                drive,
                thermal_load_controls(),
                cx,
            )
            .unwrap();
        let next = hot.candidate().contact_state();
        let account = hot.receipt();
        let scale = area * c * h * h;
        assert_eq!(
            account.kind,
            FiniteGapPlasticContinuationKind::MaterialRefreshWithDrive
        );
        assert_eq!(hot.original_source_identity(), law.identity());
        assert_eq!(hot.drive().boundary, drive);
        assert_eq!(
            hot.drive().predecessor_state_identity,
            before.contact_state().identity()
        );
        assert_eq!(hot.drive().response_identity, account.response.identity);
        common::close(hot.drive().achieved_force_n, 0.5 * area * h);
        common::close(next.approach_m(), 2.5 * c * h);
        common::close(next.plastic_depth_m()[12], 2.0 * c * h);
        common::close(account.response.plastic_work_j, 0.5 * scale);
        common::close(account.elastic_energy_change_j, -0.375 * scale);
        common::close(account.normal_work_j, 0.375 * scale);
        common::close(account.parameter_work_j, -0.25 * scale);
        common::close(hot.drive().spring_energy_j, 0.125 * scale);
        common::close(
            account.normal_work_j,
            loaded.drive.spring_energy_j - hot.drive().spring_energy_j,
        );
        assert_eq!(
            hot.source().hardness_sources()[1],
            law.hardness_sources()[1]
        );
        assert_eq!(
            hot.source().hardness_sources()[0].key(),
            law.hardness_sources()[0].key()
        );
        assert!(
            law.advance_with_drive(hot.candidate(), drive, thermal_load_controls(), cx)
                .is_err()
        );

        let cooled = hot
            .source()
            .refresh_at_temperatures_with_drive(
                hot.candidate(),
                cards.each_ref(),
                [COLD; 2],
                COLD,
                drive,
                thermal_load_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(
            cooled.candidate().contact_state().plastic_depth_m(),
            next.plastic_depth_m()
        );
        assert_eq!(
            cooled
                .candidate()
                .contact_state()
                .cumulative_plastic_work_j(),
            next.cumulative_plastic_work_j()
        );
        assert_eq!(cooled.receipt().response.plastic_work_j, 0.0);
        common::close(cooled.drive().achieved_force_n, 0.5 * area * h);
        common::close(
            cooled.candidate().contact_state().approach_m(),
            next.approach_m(),
        );
        let (source, state, _, receipt) = cooled.into_parts();
        assert_eq!(receipt.boundary, drive);
        assert!(
            source
                .advance_with_drive(&state, drive, thermal_load_controls(), cx)
                .is_ok()
        );
        assert_eq!(before, saved);
    });
}

#[test]
fn g1_independent_body_temperature_follows_the_elastic_series_spring_reference() {
    with_cx(|cx, _| {
        let (cards, law) = thermal_source([1.0e8; 2], [1.0e8; 2], [1.0e7, 2.0e7], envelope(), cx);
        let c = self_compliance(&law);
        let area = 1.0e-8;
        let drive = FiniteGapPlasticDrive {
            stiffness_n_per_m: area / c,
            position_m: 1.0e-6,
        };
        let loaded = law
            .advance_with_drive(&law.initial_state(), drive, thermal_load_controls(), cx)
            .unwrap();
        let before = loaded.candidate;
        let old_pressure = drive.position_m / (2.0 * c);
        common::close(before.contact_state().pressure_pa()[12], old_pressure);
        common::close(before.contact_state().approach_m(), 0.5 * drive.position_m);

        // Halving only E_A increases the two-body compliance by 3/2. The
        // series-spring equation p=(u-z-r)/(C+A/K) gives p_new=4*p_old/5.
        let hot = law
            .refresh_at_temperatures_with_drive(
                &before,
                cards.each_ref(),
                [HOT, COLD],
                COLD,
                drive,
                thermal_load_controls(),
                cx,
            )
            .unwrap();
        let next = hot.candidate().contact_state();
        let scale = drive.stiffness_n_per_m * drive.position_m * drive.position_m;
        common::close(self_compliance(hot.source()), 1.5 * c);
        common::close(next.pressure_pa()[12], 0.8 * old_pressure);
        common::close(next.approach_m(), 0.6 * drive.position_m);
        common::close(hot.drive().achieved_force_n, 0.8 * area * old_pressure);
        common::close(hot.drive().spring_energy_j, 0.08 * scale);
        common::close(hot.receipt().elastic_energy_j, 0.12 * scale);
        common::close(hot.receipt().elastic_energy_change_j, -0.005 * scale);
        common::close(hot.receipt().normal_work_j, 0.045 * scale);
        common::close(hot.receipt().parameter_work_j, -0.05 * scale);
        assert_eq!(hot.receipt().response.plastic_work_j, 0.0);
        assert_eq!(
            next.plastic_depth_m(),
            before.contact_state().plastic_depth_m()
        );
        assert_eq!(
            hot.source().hardness_sources()[1],
            law.hardness_sources()[1]
        );
    });
}
