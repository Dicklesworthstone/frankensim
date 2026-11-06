//! Original source-bound local fretting, Archard limits and joint acceptance.
use super::*;
use fs_contact::interface_binding::mindlin::fretting::{
    SourceFrettingError, SourceMindlinSurfaceWearAdvance, SourceMindlinWear,
};
use fs_matdb::{
    HardnessLoadStep, HardnessTestContext, MaterialCard, ObservationDataset, ObservationId,
    QueryPoint,
};
use fs_material::hardness::INDENTATION_HARDNESS_QUANTITY;
use fs_material::state_point::{ResolvedMaterialStatePoint, resolve_material_state_point};
use fs_qty::QtyAny;
use fs_tribo::WearState;

#[path = "fretting/loading.rs"]
mod loading;

const COEFFICIENTS: [f64; 2] = [2e-6, 7e-6];
const HARDNESS: [f64; 2] = [1e9, 4e8];
fn coefficient_keys() -> [PropertyKey; 2] {
    ["wear-a", "wear-b"].map(|name| PropertyKey::new(name, Dims::NONE))
}
fn hardness_key() -> PropertyKey {
    PropertyKey::with_quantity("indentation-hardness", INDENTATION_HARDNESS_QUANTITY)
        .with_hardness_test(
            HardnessTestContext::new(
                "synthetic spherical diamond indenter, radius 1 mm",
                vec![
                    HardnessLoadStep::new(
                        QtyAny::new(10.0, Dims([1, 1, -2, 0, 0, 0])),
                        QtyAny::new(5.0, Dims([0, 0, 1, 0, 0, 0])),
                    )
                    .unwrap(),
                ],
                "synthetic projected-contact-area protocol revision 1",
                ObservationId(hardness_observation().content_hash()),
            )
            .unwrap(),
        )
        .unwrap()
}
fn hardness_observation() -> ObservationDataset {
    ObservationDataset {
        specimen: "synthetic original fretting surfaces".into(),
        method: "synthetic indentation pressure".into(),
        artifact: fs_blake3::hash_domain("fixture", b"synthetic fretting hardness"),
        caveats: "numerical fixture only; no measured calibration".into(),
        provenance: Provenance {
            source: "synthetic source fretting fixture; not measured data".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}
fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 3e-10 * a.abs().max(b.abs()).max(1e-300),
        "{a:e} != {b:e}"
    );
}
fn claim(key: PropertyKey, value: f64) -> PropertyClaim {
    let observations = key
        .hardness_test()
        .map(HardnessTestContext::observation)
        .into_iter()
        .collect();
    PropertyClaim {
        value: PropertyValue::Scalar {
            value,
            dims: key.quantity().dims(),
        },
        key,
        validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations,
        provenance: Provenance {
            source: "synthetic source fretting fixture; not measured data".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}
struct Fixture {
    cards: [MaterialCard; 2],
    wear: SourceMindlinWear,
}
fn fixture(mu: f64, coefficients: [f64; 2], hardness: [f64; 2]) -> Fixture {
    // Hardness is part of each ORIGINAL full card before any elastic binding.
    let cards = [0, 1].map(|i| {
        let original = common::solid(["A", "B"][i], [1000.0, 2000.0][i], [1e7, 2e7][i]);
        let mut claims = original.claims().clone();
        claims.register_observation(hardness_observation()).unwrap();
        claims
            .insert_claim(claim(hardness_key(), hardness[i]))
            .unwrap();
        MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap()
    });
    let base = common::bind(&cards[0], &cards[1], None, 0.0);
    let original = base.normal_state().elastic().interface().card();
    let mut claims = original.claims().clone();
    claims
        .insert_claim(claim(PropertyKey::new(MU, Dims::NONE), mu))
        .unwrap();
    for (key, value) in coefficient_keys().into_iter().zip(coefficients) {
        claims.insert_claim(claim(key, value)).unwrap();
    }
    let card = InterfaceSystemCard::assemble(
        original.surface_a().clone(),
        original.surface_b().clone(),
        original.context().clone(),
        claims,
        original.models().to_vec(),
    )
    .unwrap();
    let point = common::query();
    let selection = MaterialPropertySelection::SingleClaimOnly;
    let solid = cards
        .each_ref()
        .map(|card| resolve_isotropic_solid_state_point(card, &point, selection.clone()).unwrap());
    let dry = bind_dry_interface_system_card(&card, InputAuthority::SyntheticFixture).unwrap();
    let elastic = bind_isotropic_elastic_interface(&dry, &solid[0], &solid[1]).unwrap();
    let adhesion = resolve_interface_state_point(
        &card,
        &point,
        &[ScalarPropertyRequirement::try_new(
            ADHESION_ENERGY_PROPERTY,
            ADHESION_ENERGY_DIMS,
            ScalarAdmissibility::NonNegative,
        )
        .unwrap()],
        selection.clone(),
    )
    .unwrap();
    let state = bind_normal_interface_state(&elastic, &adhesion).unwrap();
    let model =
        bind_normal_contact_model(&state, NormalContactModelSelection::SingleSupported).unwrap();
    let contact = resolve(
        SourcePointContact::new(
            model,
            PointContactGeometry::SpherePlane {
                effective_radius_m: 0.1,
            },
            common::envelope(),
        )
        .unwrap(),
    )
    .unwrap();
    let requirements = coefficient_keys().each_ref().map(|key| {
        ScalarPropertyRequirement::try_with_key(key, ScalarAdmissibility::NonNegative).unwrap()
    });
    let coefficients =
        resolve_interface_state_point(&card, &point, &requirements, selection).unwrap();
    let hardness = cards.each_ref().map(|card| hardness_at(card, &point));
    let wear = SourceMindlinWear::new(
        contact,
        coefficients,
        coefficient_keys(),
        hardness,
        [hardness_key(), hardness_key()],
    )
    .unwrap();
    Fixture { cards, wear }
}
fn hardness_at(card: &MaterialCard, point: &QueryPoint) -> ResolvedMaterialStatePoint {
    resolve_material_state_point(
        card,
        point,
        &[ScalarPropertyRequirement::try_with_key(
            &hardness_key(),
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

#[test]
fn small_reversal_fretting_uses_local_dissipative_sliding_and_preserves_sticking_centre() {
    let f = fixture(0.6, COEFFICIENTS, HARDNESS);
    let law = f.wear.contact();
    let limit = law.law().gross_slip_displacement_m();
    let mut state = f.wear.initial_state();
    let (mut measure, mut macro_travel) = (0.0, 0.0);
    for fraction in [0.04, -0.08, 0.08, -0.08, 0.08] {
        let input = step(law, fraction * limit);
        let native = law
            .advance_with_surface_sliding(state.contact_state(), &input)
            .unwrap();
        assert_eq!(native.surface_sliding(0.0).unwrap().sliding_distance_m, 0.0);
        let expected = native.tangential_step().dissipated_j / 0.6;
        let candidate = f
            .wear
            .advance_candidate(&state, &input, &mut || false)
            .unwrap();
        assert!(candidate.pressure_weighted_sliding_j > 0.0);
        near(candidate.pressure_weighted_sliding_j, expected);
        let contact = native.into_tangential_step();
        assert_eq!(candidate.tangential.state.contact_state(), &contact.state);
        assert_eq!(candidate.tangential.dissipated_j, contact.dissipated_j);
        assert_eq!(
            candidate.tangential.stored_change_j,
            contact.stored_change_j
        );
        for i in 0..2 {
            near(
                candidate.wear[i].increment_m3,
                COEFFICIENTS[i] * expected / HARDNESS[i],
            );
        }
        measure += expected;
        macro_travel += candidate.tangential.slip_increment_m[0].abs();
        state = candidate.tangential.state;
    }
    assert!(measure < 0.01 * law.normal_response().normal_force_n * macro_travel);
    for i in 0..2 {
        near(
            state.volumes_m3()[i],
            COEFFICIENTS[i] * measure / HARDNESS[i],
        );
    }
    assert!(state.contact_state().history_intervals() > 0);
}

#[test]
fn established_gross_sliding_recovers_archard_and_independent_ordered_surface_rates() {
    let f = fixture(0.6, COEFFICIENTS, HARDNESS);
    let law = f.wear.contact();
    let limit = law.law().gross_slip_displacement_m();
    let loaded = f
        .wear
        .advance_candidate(
            &f.wear.initial_state(),
            &step(law, 2.0 * limit),
            &mut || false,
        )
        .unwrap()
        .tangential
        .state;
    let distance = 0.75 * limit;
    let candidate = f
        .wear
        .advance_candidate(&loaded, &step(law, distance), &mut || false)
        .unwrap();
    let expected = law.normal_response().normal_force_n * distance;
    near(candidate.pressure_weighted_sliding_j, expected);
    near(candidate.tangential.dissipated_j, 0.6 * expected);
    assert!(candidate.tangential.stored_change_j.abs() < 1e-13 * expected);
    for i in 0..2 {
        near(
            candidate.wear[i].increment_m3,
            COEFFICIENTS[i] * expected / HARDNESS[i],
        );
        near(
            candidate.tangential.state.volumes_m3()[i],
            loaded.volumes_m3()[i] + candidate.wear[i].increment_m3,
        );
    }
    near(
        candidate.wear[1].increment_m3 / candidate.wear[0].increment_m3,
        (COEFFICIENTS[1] / HARDNESS[1]) / (COEFFICIENTS[0] / HARDNESS[0]),
    );
    // Zero local traction does not manufacture elasticity or erase declared k.
    let free = fixture(0.0, COEFFICIENTS, HARDNESS);
    let input = step(free.wear.contact(), -1e-6);
    let answer = free
        .wear
        .advance_candidate(&free.wear.initial_state(), &input, &mut || false)
        .unwrap();
    assert_eq!(answer.tangential.dissipated_j, 0.0);
    near(answer.pressure_weighted_sliding_j, input.normal_n * 1e-6);
    for i in 0..2 {
        near(
            answer.wear[i].increment_m3,
            COEFFICIENTS[i] * input.normal_n * 1e-6 / HARDNESS[i],
        );
    }
}

#[test]
fn zero_released_drive_uses_actual_compliance_motion_and_retains_independent_inventories() {
    let f = fixture(0.6, COEFFICIENTS, HARDNESS);
    let law = f.wear.contact();
    let loaded = law
        .advance(
            &law.initial_state(),
            &step(law, 0.4 * law.law().gross_slip_displacement_m()),
        )
        .unwrap();
    let old_volumes = [1e-15, 3e-15];
    let state = f
        .wear
        .state_from_contact(
            loaded.state,
            old_volumes.map(|v| WearState::new(v).unwrap()),
        )
        .unwrap();
    let mut input = step(law, 0.0);
    input.compliance_m_per_n = 0.5 / law.law().initial_stiffness_n_per_m();
    let observer = law
        .advance_with_surface_sliding(state.contact_state(), &input)
        .unwrap();
    let expected = observer.pressure_weighted_sliding_j();
    let result = f
        .wear
        .advance_candidate(&state, &input, &mut || false)
        .unwrap();
    assert!(result.tangential.slip_increment_m[0] < 0.0);
    assert!(expected > 0.0);
    assert!(result.tangential.stored_change_j < 0.0);
    for i in 0..2 {
        near(
            result.wear[i].increment_m3,
            COEFFICIENTS[i] * expected / HARDNESS[i],
        );
        near(
            result.tangential.state.volumes_m3()[i],
            old_volumes[i] + result.wear[i].increment_m3,
        );
    }
    assert_eq!(state.volumes_m3(), old_volumes);
}

#[test]
fn original_sources_fixed_temperature_and_late_wear_refusals_preserve_exact_retry() {
    let f = fixture(0.6, COEFFICIENTS, HARDNESS);
    let h = f.wear.hardness_sources();
    let bind = |coefficients, hardness| {
        SourceMindlinWear::new(
            f.wear.contact().clone(),
            coefficients,
            coefficient_keys(),
            hardness,
            [hardness_key(), hardness_key()],
        )
    };
    assert!(
        bind(
            f.wear.coefficient_source().clone(),
            [h[1].clone(), h[0].clone()]
        )
        .is_err()
    );
    let warmer = hardness_at(&f.cards[0], &common::query().with("T", 300.0).unwrap());
    assert!(bind(f.wear.coefficient_source().clone(), [warmer, h[1].clone()]).is_err());
    let mut changed = f.cards[0].claims().clone();
    changed
        .insert_claim(claim(
            PropertyKey::new("other-original-property", Dims::NONE),
            1.0,
        ))
        .unwrap();
    let replacement = MaterialCard::assemble(f.cards[0].id().clone(), changed, vec![]).unwrap();
    assert!(
        bind(
            f.wear.coefficient_source().clone(),
            [hardness_at(&replacement, &common::query()), h[1].clone()]
        )
        .is_err()
    );
    let other = fixture(0.6, [3e-6, COEFFICIENTS[1]], HARDNESS);
    assert!(bind(other.wear.coefficient_source().clone(), h.clone()).is_err());
    let law = f.wear.contact();
    let loaded = f
        .wear
        .advance_candidate(
            &f.wear.initial_state(),
            &step(law, 0.3 * law.law().gross_slip_displacement_m()),
            &mut || false,
        )
        .unwrap()
        .tangential
        .state;
    let saved = loaded.clone();
    let input = step(law, -0.1 * law.law().gross_slip_displacement_m());
    let expected = f
        .wear
        .advance_candidate(&saved, &input, &mut || false)
        .unwrap();
    let mut invalid = input;
    invalid.contact_temperature_k = Some(law.temperature_kelvin() + 1.0);
    assert!(
        f.wear
            .advance_candidate(&saved, &invalid, &mut || false)
            .is_err()
    );
    assert!(
        other
            .wear
            .advance_candidate(&saved, &input, &mut || false)
            .is_err()
    );
    let mut polls = 0;
    let cancelled = f.wear.advance_candidate(&saved, &input, &mut || {
        polls += 1;
        polls == 3
    });
    assert_eq!(cancelled.unwrap_err(), SourceFrettingError::Cancelled);
    assert_eq!(polls, 3);
    assert_eq!(loaded, saved);
    assert_eq!(
        f.wear
            .advance_candidate(&saved, &input, &mut || false)
            .unwrap(),
        expected
    );
    // Surface A's finite candidate is computed before B overflows; neither may publish.
    let huge = fixture(0.6, [COEFFICIENTS[0], 1e308], [HARDNESS[0], 1e-300]);
    let state = huge.wear.initial_state();
    let before = state.clone();
    let input = step(
        huge.wear.contact(),
        0.2 * huge.wear.contact().law().gross_slip_displacement_m(),
    );
    assert!(
        huge.wear
            .advance_candidate(&state, &input, &mut || false)
            .is_err()
    );
    assert_eq!(state, before);
    let rest = step(huge.wear.contact(), 0.0);
    let zero = huge
        .wear
        .advance_candidate(&state, &rest, &mut || false)
        .unwrap();
    assert_eq!(zero.tangential.state, before);
    assert_eq!(zero.wear.map(|w| w.increment_m3), [0.0; 2]);
}

#[test]
fn g1_g3_local_wear_retains_sticking_core_and_each_original_surface_hardness() {
    let f = fixture(0.6, COEFFICIENTS, HARDNESS);
    let harder = fixture(0.6, COEFFICIENTS, [2.0 * HARDNESS[0], 4.0 * HARDNESS[1]]);
    let initial = f.wear.initial_state();
    let harder_initial = harder.wear.initial_state();
    let law = f.wear.contact();
    let radius = law.law().parameters().contact_radius_m;
    let input = step(law, 0.4 * law.law().gross_slip_displacement_m());
    let observed = f
        .wear
        .advance_with_surface_wear(&initial, &input, &mut || false)
        .unwrap();
    let hard = harder
        .wear
        .advance_with_surface_wear(&harder_initial, &input, &mut || false)
        .unwrap();
    // Virgin Cattaneo--Mindlin loading has c/a = sqrt(1 - du/d_gross).
    // Local removal must preserve that sticking disk, not smear total volume
    // over the nominal patch or charge recoverable displacement as sliding.
    for normalized in [0.0, 0.3, 0.7] {
        let local = observed.surface_wear(normalized * radius).unwrap();
        assert!(local.sliding.normal_pressure_pa > 0.0);
        assert_eq!(local.depth_increment_m, [0.0; 2]);
    }
    for normalized in [0.8, 0.9, 0.99] {
        let local = observed.surface_wear(normalized * radius).unwrap();
        let resistant = hard.surface_wear(normalized * radius).unwrap();
        assert!(local.depth_increment_m.into_iter().all(|depth| depth > 0.0));
        assert_eq!(local.sliding, resistant.sliding);
        near(
            resistant.depth_increment_m[0],
            0.5 * local.depth_increment_m[0],
        );
        near(
            resistant.depth_increment_m[1],
            0.25 * local.depth_increment_m[1],
        );
        near(
            local.depth_increment_m[1] / local.depth_increment_m[0],
            (COEFFICIENTS[1] / HARDNESS[1]) / (COEFFICIENTS[0] / HARDNESS[0]),
        );
    }
    let edge = observed.surface_wear(radius).unwrap();
    assert!(edge.sliding.sliding_distance_m > 0.0);
    assert_eq!(edge.depth_increment_m, [0.0; 2]);
    assert_eq!(
        observed.tangential_step().dissipated_j,
        hard.tangential_step().dissipated_j
    );
    assert_eq!(
        observed.into_candidate(),
        f.wear
            .advance_candidate(&initial, &input, &mut || false)
            .unwrap()
    );
}

fn integrated_local_depth(
    observed: &SourceMindlinSurfaceWearAdvance<'_>,
    radius_m: f64,
    cells: usize,
) -> [f64; 2] {
    // Independent equal-area annuli: integrate actual returned depth, not the
    // native pressure/sliding integral used by the total-volume calculation.
    let mut volume = [0.0; 2];
    for i in 0..cells {
        let r = radius_m * ((i as f64 + 0.5) / cells as f64).sqrt();
        let local = observed.surface_wear(r).unwrap();
        for (total, depth) in volume.iter_mut().zip(local.depth_increment_m) {
            *total += depth;
        }
    }
    volume.map(|sum| core::f64::consts::PI * radius_m * radius_m * sum / cells as f64)
}

#[test]
fn g1_local_depth_integrals_recover_analytic_volumes_through_reversal_and_gross_slip() {
    let f = fixture(0.6, COEFFICIENTS, HARDNESS);
    let law = f.wear.contact();
    let radius = law.law().parameters().contact_radius_m;
    let limit = law.law().gross_slip_displacement_m();
    let mut state = f.wear.initial_state();
    for fraction in [0.7, -0.5, 1.6, 0.3, -3.0] {
        let input = step(law, fraction * limit);
        let native = law.advance(state.contact_state(), &input).unwrap();
        let observed = f
            .wear
            .advance_with_surface_wear(&state, &input, &mut || false)
            .unwrap();
        let exact = observed.wear().each_ref().map(|wear| wear.increment_m3);
        let volumes =
            [512, 2048, 8192].map(|cells| integrated_local_depth(&observed, radius, cells));
        for side in 0..2 {
            let errors = volumes.map(|volume| (volume[side] - exact[side]).abs());
            assert!(errors[2] < 0.2 * errors[0], "{fraction}: {errors:?}");
            assert!(errors[2] < 3e-5 * exact[side], "{fraction}: {errors:?}");
        }
        near(
            observed.pressure_weighted_sliding_j(),
            native.dissipated_j / 0.6,
        );
        let candidate = observed.into_candidate();
        assert_eq!(candidate.tangential.state.contact_state(), &native.state);
        assert_eq!(candidate.tangential.dissipated_j, native.dissipated_j);
        assert_eq!(candidate.tangential.stored_change_j, native.stored_change_j);
        state = candidate.tangential.state;
    }
}

#[test]
fn g0_g1_local_wear_has_frictionless_hertz_limit_and_refuses_invalid_samples() {
    let f = fixture(0.0, [0.0, COEFFICIENTS[1]], HARDNESS);
    let law = f.wear.contact();
    let radius = law.law().parameters().contact_radius_m;
    let normal_n = law.normal_response().normal_force_n;
    let initial = f.wear.initial_state();
    let distance = 1e-6;
    let observed = f
        .wear
        .advance_with_surface_wear(&initial, &step(law, -distance), &mut || false)
        .unwrap();
    assert_eq!(observed.tangential_step().dissipated_j, 0.0);
    for normalized in [0.0_f64, 0.5, 0.9, 1.0] {
        let local = observed.surface_wear(normalized * radius).unwrap();
        let pressure = 1.5 * normal_n / (core::f64::consts::PI * radius * radius)
            * (1.0 - normalized * normalized).sqrt();
        assert_eq!(local.depth_increment_m[0], 0.0);
        near(local.sliding.signed_sliding_increment_m, -distance);
        near(
            local.depth_increment_m[1],
            COEFFICIENTS[1] * pressure * distance / HARDNESS[1],
        );
    }
    let centre = observed.surface_wear(0.0).unwrap();
    for invalid_radius in [-radius, 1.1 * radius, f64::NAN, f64::INFINITY] {
        assert!(observed.surface_wear(invalid_radius).is_err());
    }
    assert_eq!(observed.surface_wear(0.0).unwrap(), centre);
    let volume = observed.into_candidate();
    near(
        volume.wear[1].increment_m3,
        COEFFICIENTS[1] * normal_n * distance / HARDNESS[1],
    );
    // A finite area-integrated volume does not guarantee a representable local
    // depth. Refusing one sample leaves the already validated candidate intact.
    let extreme = fixture(0.0, [0.0, 1e10], [HARDNESS[0], 1e-300]);
    let initial = extreme.wear.initial_state();
    let observed = extreme
        .wear
        .advance_with_surface_wear(
            &initial,
            &step(extreme.wear.contact(), distance),
            &mut || false,
        )
        .unwrap();
    assert!(observed.wear()[1].increment_m3.is_finite());
    assert!(observed.surface_wear(0.0).is_err());
    assert!(observed.into_candidate().wear[1].increment_m3.is_finite());
}
