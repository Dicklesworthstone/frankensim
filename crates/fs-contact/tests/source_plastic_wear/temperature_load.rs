//! G1/G3: original thermal k/H changes future wear without resetting its history.

use super::*;
use fs_contact::normal_patch::plastic_source::temperature::SourcePlasticTemperatureError;
use fs_contact::normal_patch::{
    FiniteGapPlasticContinuationKind, FiniteGapPlasticDrive, FiniteGapPlasticDriveError,
    FiniteGapPlasticLoadControls,
};
use fs_matdb::ClaimSet;
use fs_qty::{Pressure, Temperature};

const COLD: f64 = 293.15;
const HOT: f64 = 313.15;

fn load_controls() -> FiniteGapPlasticLoadControls {
    FiniteGapPlasticLoadControls {
        minimum_approach_m: -1.0e-4,
        maximum_approach_m: 1.0e-4,
        absolute_force_tolerance_n: 1.0e-14,
        relative_force_tolerance: 1.0e-12,
        approach_tolerance_m: 1.0e-18,
        maximum_equilibrium_evaluations: 128,
    }
}

fn thermal_material(
    name: &str,
    cold_h: f64,
    hot_h: f64,
    hot_e: f64,
) -> (MaterialCard, PropertyKey) {
    let (base, key) = material(name, cold_h);
    let mut claims = ClaimSet::new();
    for id in base.claims().observation_ids() {
        claims
            .register_observation(base.claims().observation(id).unwrap().clone())
            .unwrap();
    }
    for (_, original) in base.claims().claims_ordered() {
        let mut claim = original.clone();
        let values = match claim.key.name() {
            "indentation-hardness" => Some((cold_h, hot_h)),
            "young_modulus" => Some((2.0e7, hot_e)),
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

fn thermal_fixture(cx: &Cx<'_>) -> ([MaterialCard; 2], SourcePlasticArchardWear) {
    let cards = [
        thermal_material("A", 1.0e6, 0.75e6, 1.0e7),
        thermal_material("B", 2.0e6, 1.5e6, 2.0e7),
    ];
    let keys = ["k_a", "k_b"].map(|name| PropertyKey::new(name, Dims::NONE));
    let base = common::bind(&cards[0].0, &cards[1].0, None, 0.0);
    let original = base.normal_state().elastic().interface().card();
    let mut claims = original.claims().clone();
    for (key, cold) in keys.clone().into_iter().zip([0.01, 0.03]) {
        let mut coefficient = claim(key, cold);
        coefficient.value = PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(COLD, cold), (HOT, 2.0 * cold)],
            dims: Dims::NONE,
        };
        coefficient.validity = ValidityDomain::unconstrained().with("T", COLD, HOT);
        coefficient.interpolation = InterpolationPolicy::LinearInside;
        claims.insert_claim(coefficient).unwrap();
    }
    let interface = InterfaceSystemCard::assemble(
        original.surface_a().clone(),
        original.surface_b().clone(),
        original.context().clone(),
        claims,
        original.models().to_vec(),
    )
    .unwrap();
    let point = common::query();
    let select = MaterialPropertySelection::SingleClaimOnly;
    let dry = bind_dry_interface_system_card(&interface, InputAuthority::SyntheticFixture).unwrap();
    let elastic = bind_isotropic_elastic_interface(
        &dry,
        base.normal_state().elastic().surface_a_state(),
        base.normal_state().elastic().surface_b_state(),
    )
    .unwrap();
    let adhesion = resolve_interface_state_point(
        &interface,
        &point,
        &[ScalarPropertyRequirement::try_new(
            ADHESION_ENERGY_PROPERTY,
            ADHESION_ENERGY_DIMS,
            ScalarAdmissibility::NonNegative,
        )
        .unwrap()],
        select.clone(),
    )
    .unwrap();
    let normal = bind_normal_interface_state(&elastic, &adhesion).unwrap();
    let hardness = cards.each_ref().map(|(card, key)| {
        let source = resolve_material_state_point(
            card,
            &point,
            &[
                ScalarPropertyRequirement::try_with_key(key, ScalarAdmissibility::StrictlyPositive)
                    .unwrap(),
            ],
            select.clone(),
        )
        .unwrap();
        IndentationHardnessStatePoint::try_from_resolved(source, key.clone()).unwrap()
    });
    let coefficients = resolve_interface_state_point(
        &interface,
        &point,
        &keys.each_ref().map(|key| {
            ScalarPropertyRequirement::try_with_key(key, ScalarAdmissibility::NonNegative).unwrap()
        }),
        select,
    )
    .unwrap();
    let mut gap = vec![1.0e-3; 25];
    gap[12] = 0.0;
    let contact = SourcePlasticContact::new(
        normal,
        hardness.clone(),
        FiniteGapGrid {
            cells_x: 5,
            cells_y: 5,
            cell_width_m: 1.0e-4,
            cell_depth_m: 1.0e-4,
        },
        gap,
        PlasticContactControls {
            maximum_active_set_iterations: 100,
            complementarity_tolerance_m: 1.0e-12,
            boundary_clearance_cells: 1,
        },
        PlasticContactEnvelope {
            maximum_approach_m: 1.0e-4,
            maximum_plastic_depth_m: 1.0e-4,
            half_space_depth_m: [1.0e-2; 2],
            maximum_patch_to_depth: [0.1; 2],
        },
        cx,
    )
    .unwrap();
    let wear =
        SourcePlasticArchardWear::new(contact, coefficients, keys, hardness, [1.0e-4; 2]).unwrap();
    (cards.map(|(card, _)| card), wear)
}

#[test]
fn g1_hot_body_and_interface_change_future_wear_at_the_maintained_force() {
    with_cx(|cx, _| {
        let (cards, law) = thermal_fixture(cx);
        let dented = law
            .advance_normal(&law.initial_state(), 5.0e-5, cx)
            .unwrap()
            .candidate;
        let force = 0.005;
        let loaded = law
            .advance_normal_at_load(&dented, force, load_controls(), cx)
            .unwrap()
            .candidate;
        let first = law
            .advance_wear_at_load(&loaded, &[1.0e-4; 25], force, load_controls(), cx)
            .unwrap();
        let before = first.wear.candidate;
        let saved = before.clone();
        let old = before.contact_state().contact_state();
        let old_energy = first.wear.receipt.elastic_energy_j;
        let unchanged = law
            .refresh_at_temperatures_at_load(
                &before,
                cards.each_ref(),
                [COLD; 2],
                COLD,
                force,
                load_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(unchanged.source(), &law);
        assert_eq!(unchanged.candidate(), &before);
        assert_eq!(unchanged.load().equilibrium_evaluations, 0);
        assert_eq!(unchanged.receipt().parameter_work_j, 0.0);

        let hot = law
            .refresh_at_temperatures_at_load(
                &before,
                cards.each_ref(),
                [HOT, COLD],
                0.5 * (HOT + COLD),
                force,
                load_controls(),
                cx,
            )
            .unwrap();
        let next = hot.candidate().contact_state().contact_state();
        assert_eq!(
            hot.receipt().kind,
            FiniteGapPlasticContinuationKind::MaterialRefreshAtLoad
        );
        close(hot.load().achieved_force_n, force);
        close(next.pressure_pa()[12], 0.5e6);
        close(hot.receipt().elastic_energy_j, 1.5 * old_energy);
        close(hot.receipt().normal_work_j, old_energy);
        close(hot.receipt().parameter_work_j, -0.5 * old_energy);
        assert_eq!(hot.candidate().volumes_m3(), before.volumes_m3());
        assert_eq!(hot.candidate().recession_m(), before.recession_m());
        assert_eq!(next.plastic_depth_m(), old.plastic_depth_m());
        assert_eq!(next.recession_m(), old.recession_m());
        assert_eq!(
            next.cumulative_plastic_work_j(),
            old.cumulative_plastic_work_j()
        );
        assert_eq!(hot.receipt().response.plastic_work_j, 0.0);
        for side in 0..2 {
            let new_h = &hot.source().hardness_sources()[side];
            let old_h = &law.hardness_sources()[side];
            close(new_h.hardness_pa(), [0.75e6, 2.0e6][side]);
            assert_eq!(new_h.key(), old_h.key());
            assert_eq!(
                new_h.resolved().card_identity(),
                old_h.resolved().card_identity()
            );
            assert_eq!(
                new_h.property().answer().receipt.selected,
                old_h.property().answer().receipt.selected
            );
            let key = &law.coefficient_keys()[side];
            let new_k = hot
                .source()
                .coefficient_source()
                .property_by_key(key)
                .unwrap();
            close(new_k.value_si(), [0.015, 0.045][side]);
            assert_eq!(
                new_k.answer().receipt.selected,
                law.coefficient_source()
                    .property_by_key(key)
                    .unwrap()
                    .answer()
                    .receipt
                    .selected
            );
        }
        assert_eq!(
            hot.source().hardness_sources()[1],
            law.hardness_sources()[1]
        );
        assert!(
            law.advance_wear_at_load(hot.candidate(), &[1.0e-4; 25], force, load_controls(), cx)
                .is_err()
        );

        // Independently evaluated original k(T_interface)/H(T_body) for each body.
        let continued = hot
            .source()
            .advance_wear_at_load(hot.candidate(), &[1.0e-4; 25], force, load_controls(), cx)
            .unwrap();
        for side in 0..2 {
            let expected = [1.0e-14, 1.125e-14][side];
            close(continued.wear.wear[side].increment_m3, expected);
            close(
                continued.wear.candidate.volumes_m3()[side],
                before.volumes_m3()[side] + expected,
            );
        }
        let cooled = hot
            .source()
            .refresh_at_temperatures_at_load(
                &continued.wear.candidate,
                cards.each_ref(),
                [COLD; 2],
                COLD,
                force,
                load_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(
            cooled.candidate().volumes_m3(),
            continued.wear.candidate.volumes_m3()
        );
        assert_eq!(
            cooled.candidate().recession_m(),
            continued.wear.candidate.recession_m()
        );
        assert_eq!(
            cooled
                .candidate()
                .contact_state()
                .contact_state()
                .plastic_depth_m(),
            old.plastic_depth_m()
        );
        let (source, candidate, _, load) = cooled.into_parts();
        close(load.achieved_force_n, force);
        assert!(
            source
                .advance_wear_at_load(&candidate, &[0.0; 25], force, load_controls(), cx)
                .is_ok()
        );
        assert_eq!(before, saved);
    });
}

#[test]
fn g0_temperature_load_and_late_coefficient_refusals_leave_wear_unchanged() {
    with_cx(|cx, gate| {
        let (cards, law) = thermal_fixture(cx);
        let force = 0.005;
        let loaded = law
            .advance_normal_at_load(&law.initial_state(), force, load_controls(), cx)
            .unwrap()
            .candidate;
        let before = law
            .advance_wear_at_load(&loaded, &[1.0e-4; 25], force, load_controls(), cx)
            .unwrap()
            .wear
            .candidate;
        let saved = before.clone();
        // 318 K is supported by the original adhesion/normal source but lies
        // outside both original coefficient curves; refuse after contact refresh.
        assert!(matches!(
            law.refresh_at_temperatures_at_load(
                &before,
                cards.each_ref(),
                [COLD; 2],
                318.0,
                force,
                load_controls(),
                cx
            ),
            Err(SourcePlasticWearError::Material(_))
        ));
        assert!(
            law.refresh_at_temperatures_at_load(
                &before,
                cards.each_ref(),
                [HOT, COLD],
                COLD,
                0.004,
                load_controls(),
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
                load_controls(),
                cx,
            )
            .is_err()
        );
        let revised =
            MaterialCard::supersede(&cards[0], cards[0].claims().clone(), vec![]).unwrap();
        assert!(
            law.refresh_at_temperatures_at_load(
                &before,
                [&revised, &cards[1]],
                [COLD; 2],
                COLD,
                force,
                load_controls(),
                cx,
            )
            .is_err()
        );
        let outside = FiniteGapPlasticLoadControls {
            maximum_approach_m: 1.01e-4,
            ..load_controls()
        };
        assert!(
            law.refresh_at_temperatures_at_load(
                &before,
                cards.each_ref(),
                [COLD; 2],
                COLD,
                force,
                outside,
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
                load_controls(),
                cx,
            )
            .is_err()
        );
        assert_eq!(before, saved);
    });
}

#[test]
fn g1_worn_contact_heats_and_continues_original_archard_wear_with_the_device_held() {
    with_cx(|cx, _| {
        let (cards, law) = thermal_fixture(cx);
        let area = 1.0e-8;
        let width = law.contact().model().config().grid.cell_width_m;
        let modulus = law.contact().normal_source().elastic().reduced_modulus_pa();
        let c = 4.0 * width * (1.0 + 2.0_f64.sqrt()).ln() / (core::f64::consts::PI * modulus);
        let h = 1.0e6;
        let drive = FiniteGapPlasticDrive {
            stiffness_n_per_m: area / c,
            position_m: 3.0 * c * h,
        };
        let loaded = law
            .advance_normal_with_drive(&law.initial_state(), drive, load_controls(), cx)
            .unwrap()
            .candidate;
        let first = law
            .advance_wear_with_drive(&loaded, &[1.0e-5; 25], drive, load_controls(), cx)
            .unwrap();
        let before = first.wear.candidate;
        let saved = before.clone();
        let old = before.contact_state().contact_state();
        let recession = old.recession_m()[12];
        close(recession, 2.5e-7);
        close(old.plastic_depth_m()[12], c * h);
        let unchanged = law
            .refresh_at_temperatures_with_drive(
                &before,
                cards.each_ref(),
                [COLD; 2],
                COLD,
                drive,
                load_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(unchanged.source(), &law);
        assert_eq!(unchanged.candidate(), &before);
        assert_eq!(unchanged.drive().equilibrium_evaluations, 0);

        let hot = law
            .refresh_at_temperatures_with_drive(
                &before,
                cards.each_ref(),
                [HOT, COLD],
                0.5 * (HOT + COLD),
                drive,
                load_controls(),
                cx,
            )
            .unwrap();
        let next = hot.candidate().contact_state().contact_state();
        assert_eq!(
            hot.receipt().kind,
            FiniteGapPlasticContinuationKind::MaterialRefreshWithDrive
        );
        assert_eq!(hot.drive().boundary, drive);
        close(next.pressure_pa()[12], 0.75 * h);
        close(next.approach_m(), 2.25 * c * h);
        close(next.plastic_depth_m()[12], 1.125 * c * h - recession);
        close(
            hot.receipt().response.plastic_work_j,
            area * 0.75 * h * (0.125 * c * h - recession),
        );
        close(
            hot.receipt().normal_work_j,
            first.drive.spring_energy_j - hot.drive().spring_energy_j,
        );
        assert_eq!(hot.candidate().volumes_m3(), before.volumes_m3());
        assert_eq!(hot.candidate().recession_m(), before.recession_m());
        assert_eq!(next.recession_m(), old.recession_m());
        assert!(next.plastic_depth_m()[12] > old.plastic_depth_m()[12]);
        for side in 0..2 {
            close(
                hot.source().hardness_sources()[side].hardness_pa(),
                [0.75e6, 2.0e6][side],
            );
            assert_eq!(
                hot.source().hardness_sources()[side].key(),
                law.hardness_sources()[side].key()
            );
            assert_eq!(
                hot.source().coefficient_keys()[side],
                law.coefficient_keys()[side]
            );
        }
        assert!(
            law.advance_wear_with_drive(hot.candidate(), &[1.0e-5; 25], drive, load_controls(), cx)
                .is_err()
        );

        // At the new 0.75 MPa pressure, the independently prescribed interface
        // k values 0.015/0.045 and body H values 0.75/2 MPa set different removals.
        let continued = hot
            .source()
            .advance_wear_with_drive(hot.candidate(), &[1.0e-5; 25], drive, load_controls(), cx)
            .unwrap();
        for side in 0..2 {
            let volume = [1.5e-15, 1.6875e-15][side];
            close(continued.wear.wear[side].increment_m3, volume);
            close(
                continued.wear.candidate.volumes_m3()[side],
                before.volumes_m3()[side] + volume,
            );
        }
        close(
            continued
                .wear
                .candidate
                .contact_state()
                .contact_state()
                .pressure_pa()[12],
            0.75 * h - 3.1875e-7 / (2.5 * c),
        );
        assert_eq!(continued.wear.receipt.response.plastic_work_j, 0.0);
        assert_eq!(
            continued
                .wear
                .candidate
                .contact_state()
                .contact_state()
                .plastic_depth_m(),
            next.plastic_depth_m()
        );
        let (source, state, _, receipt) = hot.into_parts();
        assert_eq!(receipt.boundary, drive);
        assert!(
            source
                .advance_wear_with_drive(&state, &[0.0; 25], drive, load_controls(), cx)
                .is_ok()
        );
        assert_eq!(before, saved);
    });
}

#[test]
fn g0_thermal_device_and_late_wear_source_refusals_are_atomic() {
    with_cx(|cx, gate| {
        let (cards, law) = thermal_fixture(cx);
        let drive = FiniteGapPlasticDrive {
            stiffness_n_per_m: 1000.0,
            position_m: 1.0e-6,
        };
        let loaded = law
            .advance_normal_with_drive(&law.initial_state(), drive, load_controls(), cx)
            .unwrap()
            .candidate;
        let before = law
            .advance_wear_with_drive(&loaded, &[1.0e-5; 25], drive, load_controls(), cx)
            .unwrap()
            .wear
            .candidate;
        let saved = before.clone();
        assert!(matches!(
            law.refresh_at_temperatures_with_drive(
                &before,
                cards.each_ref(),
                [COLD; 2],
                318.0,
                drive,
                load_controls(),
                cx,
            ),
            Err(SourcePlasticWearError::Material(_))
        ));
        let moved = FiniteGapPlasticDrive {
            position_m: 2.0e-6,
            ..drive
        };
        assert!(matches!(
            law.refresh_at_temperatures_with_drive(
                &before,
                cards.each_ref(),
                [HOT, COLD],
                COLD,
                moved,
                load_controls(),
                cx,
            ),
            Err(SourcePlasticWearError::Temperature(
                SourcePlasticTemperatureError::Drive(
                    FiniteGapPlasticDriveError::PredecessorDriveMismatch { .. }
                )
            ))
        ));
        let revised =
            MaterialCard::supersede(&cards[0], cards[0].claims().clone(), vec![]).unwrap();
        assert!(
            law.refresh_at_temperatures_with_drive(
                &before,
                [&revised, &cards[1]],
                [COLD; 2],
                COLD,
                drive,
                load_controls(),
                cx,
            )
            .is_err()
        );
        assert!(
            law.refresh_at_temperatures_with_drive(
                &before,
                cards.each_ref(),
                [330.0, COLD],
                COLD,
                drive,
                load_controls(),
                cx,
            )
            .is_err()
        );
        let outside = FiniteGapPlasticLoadControls {
            maximum_approach_m: 1.01e-4,
            ..load_controls()
        };
        assert!(
            law.refresh_at_temperatures_with_drive(
                &before,
                cards.each_ref(),
                [COLD; 2],
                COLD,
                drive,
                outside,
                cx,
            )
            .is_err()
        );
        gate.request();
        assert!(
            law.refresh_at_temperatures_with_drive(
                &before,
                cards.each_ref(),
                [COLD; 2],
                COLD,
                drive,
                load_controls(),
                cx,
            )
            .is_err()
        );
        assert_eq!(before, saved);
    });
}
