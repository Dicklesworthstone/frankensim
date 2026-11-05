//! G1/G3: original temperature-dependent hardness changes an existing dent.

use super::*;
use fs_contact::normal_patch::plastic_source::temperature::SourcePlasticTemperatureError;
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
