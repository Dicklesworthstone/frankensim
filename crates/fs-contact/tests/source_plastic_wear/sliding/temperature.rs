//! Independent temperatures change real H/k/mu without erasing accepted shear.

use super::*;
use fs_contact::normal_patch::FiniteGapPlasticLinearSolver;
use fs_matdb::ClaimSet;
use fs_qty::Temperature;

const COLD: f64 = 293.15;
const HOT: f64 = 313.15;

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

fn thermal_source(
    pressure_domain_pa: [f64; 2],
    cx: &Cx<'_>,
) -> ([MaterialCard; 2], SourcePlasticSlidingWear) {
    let cards = [
        thermal_material("A", 1.0e6, 0.5e6, 1.0e7),
        thermal_material("B", 2.0e6, 1.5e6, 2.0e7),
    ];
    let keys = ["k_a", "k_b"].map(|name| PropertyKey::new(name, Dims::NONE));
    let base = common::bind(&cards[0].0, &cards[1].0, None, 0.0);
    let original = base.normal_state().elastic().interface().card();
    let mut claims = original.claims().clone();
    let pressure_quantity =
        QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static));
    for (name, cold, hot) in [
        ("k_a", 0.01, 0.02),
        ("k_b", 0.03, 0.06),
        ("mu_static", 0.6, 0.3),
        ("mu_kinetic", 0.4, 0.2),
    ] {
        let mut coefficient = claim(PropertyKey::new(name, Dims::NONE), cold);
        coefficient.value = PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(COLD, cold), (HOT, hot)],
            dims: Dims::NONE,
        };
        coefficient.validity = ValidityDomain::unconstrained().with("T", COLD, HOT);
        if name.starts_with("mu_") {
            coefficient.validity = coefficient.validity.with_quantity(
                "p",
                pressure_quantity,
                pressure_domain_pa[0],
                pressure_domain_pa[1],
            );
        }
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
        select.clone(),
    )
    .unwrap();
    let mut gap = vec![1.0e-3; 25];
    gap[12] = 0.0;
    let contact = SourcePlasticContact::new_with_linear_solver(
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
        FiniteGapPlasticLinearSolver::Fft {
            maximum_iterations: 64,
            maximum_corrections: 1,
        },
        cx,
    )
    .unwrap();
    let wear =
        SourcePlasticArchardWear::new(contact, coefficients, keys, hardness, [1.0e-4; 2]).unwrap();
    let requirement = |name| {
        ScalarPropertyRequirement::try_new(name, Dims::NONE, ScalarAdmissibility::NonNegative)
            .unwrap()
    };
    let friction = SourceThermalCoulomb::resolve(
        dry,
        point
            .with_quantity(
                "p",
                pressure_quantity,
                0.5 * pressure_domain_pa[0] + 0.5 * pressure_domain_pa[1],
            )
            .unwrap(),
        requirement("mu_static"),
        requirement("mu_kinetic"),
        select,
        ThermalFrictionAxes {
            temperature: "T".into(),
            nominal_pressure: Some("p".into()),
        },
    )
    .unwrap();
    (
        cards.map(|(card, _)| card),
        SourcePlasticSlidingWear::new(wear, friction).unwrap(),
    )
}

#[test]
fn independent_body_and_interface_temperatures_continue_original_shear_and_future_wear() {
    with_cx(|cx, _| {
        let (cards, law) = thermal_source([1.0e5, 2.0e6], cx);
        let dented = law
            .advance_normal(&law.initial_state(), 5.0e-5, cx)
            .unwrap();
        let force = 0.003;
        let loaded = law
            .advance_normal_at_load(&dented.candidate, force, normal_controls(), cx)
            .unwrap();
        let mut input = motion();
        input[12].released_increment_m = [1.0e-4, 0.0];
        let before = law
            .advance_at_load(
                &loaded.candidate,
                &input,
                0.1,
                friction_controls(),
                force,
                normal_controls(),
                cx,
            )
            .unwrap()
            .candidate;
        let saved = before.clone();
        let unchanged = law
            .refresh_at_temperatures_at_load(
                &before,
                cards.each_ref(),
                [COLD; 2],
                COLD,
                force,
                normal_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(unchanged.source(), &law);
        assert_eq!(unchanged.candidate(), &before);
        assert_eq!(unchanged.receipt().boundary.equilibrium_evaluations, 0);
        let hot = law
            .refresh_at_temperatures_at_load(
                &before,
                cards.each_ref(),
                [HOT, COLD],
                HOT,
                force,
                normal_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(before, saved);
        assert_eq!(
            hot.candidate().tangential_states(),
            before.tangential_states()
        );
        assert_eq!(
            hot.candidate().wear_state().volumes_m3(),
            before.wear_state().volumes_m3()
        );
        assert_eq!(
            hot.candidate().wear_state().recession_m(),
            before.wear_state().recession_m()
        );
        assert_eq!(hot.source().interface_temperature_k(), HOT);
        assert_eq!(
            hot.candidate().tangential_states()[12].contact_temperature_k(),
            Some(COLD)
        );
        assert_eq!(hot.source().friction_source(), law.friction_source());
        assert_eq!(hot.receipt().original_source_identity, law.identity());
        assert_ne!(hot.receipt().source_identity, law.identity());
        assert_eq!(hot.receipt().source_identity, hot.source().identity());
        assert_eq!(
            hot.receipt().normal.kind,
            FiniteGapPlasticContinuationKind::MaterialRefreshAtLoad
        );
        close(hot.receipt().boundary.achieved_force_n, force);
        assert_eq!(hot.receipt().normal.response.plastic_work_j, 0.0);
        assert_eq!(
            hot.source().wear_source().contact().model().linear_solver(),
            law.wear_source().contact().model().linear_solver()
        );
        for side in 0..2 {
            let refreshed = &hot.source().wear_source().hardness_sources()[side];
            let original = &law.wear_source().hardness_sources()[side];
            close(refreshed.hardness_pa(), [0.5e6, 2.0e6][side]);
            assert_eq!(
                refreshed.property().answer().receipt.selected,
                original.property().answer().receipt.selected
            );
            assert_eq!(refreshed.key(), original.key());
        }
        assert_eq!(hot.receipt().friction_sources.len(), 1);
        assert_eq!(hot.receipt().cell_source_indices[0], None);
        assert_eq!(hot.receipt().cell_source_indices[12], Some(0));
        let selected = &hot.receipt().friction_sources[0];
        assert_eq!(selected.state_coordinate("T"), Some(HOT));
        close(selected.state_coordinate("p").unwrap(), force / 1.0e-8);
        close(selected.property("mu_kinetic").unwrap().value_si(), 0.2);
        assert_eq!(
            selected
                .property("mu_kinetic")
                .unwrap()
                .answer()
                .receipt
                .selected,
            law.friction_source()
                .initial_source()
                .property("mu_kinetic")
                .unwrap()
                .answer()
                .receipt
                .selected
        );
        assert!(law.advance_normal(hot.candidate(), 5.0e-5, cx).is_err());

        input[12].released_increment_m = [0.0, 1.0e-5];
        let (continued, candidate, _) = hot.into_parts();
        let next = continued
            .advance_at_load(
                &candidate,
                &input,
                0.1,
                friction_controls(),
                force,
                normal_controls(),
                cx,
            )
            .unwrap();
        close(next.cells[12].response.traction_n[1], 0.2 * force);
        close(next.frictional_work_j, 0.2 * force * 1.0e-5);
        assert_eq!(
            next.candidate.tangential_states()[12].slip_m(),
            [1.0e-4, 1.0e-5]
        );
        assert_eq!(
            next.candidate.tangential_states()[12].contact_temperature_k(),
            Some(HOT)
        );
        for side in 0..2 {
            let expected = [0.02, 0.06][side] * force * 1.0e-5 / [0.5e6, 2.0e6][side];
            close(next.wear[side].increment_m3, expected);
            close(
                next.candidate.wear_state().volumes_m3()[side],
                before.wear_state().volumes_m3()[side] + expected,
            );
        }
    });
}

#[test]
fn held_approach_and_device_softening_add_normal_indentation_without_fictitious_sliding() {
    with_cx(|cx, _| {
        let (cards, law) = thermal_source([1.0e5, 2.0e6], cx);
        let loaded = law
            .advance_normal(&law.initial_state(), 5.0e-5, cx)
            .unwrap();
        let mut input = motion();
        input[12].released_increment_m = [1.0e-5, -2.0e-5];
        let first = law
            .advance(&loaded.candidate, &input, 0.1, friction_controls(), cx)
            .unwrap();
        let before = first.candidate;
        let old = before.wear_state().contact_state().contact_state();
        let approach = law
            .refresh_at_temperatures(&before, cards.each_ref(), [HOT, COLD], HOT, cx)
            .unwrap();
        let drive = FiniteGapPlasticDrive {
            stiffness_n_per_m: 1000.0,
            position_m: old.approach_m() + first.normal.response.normal.normal_force_n / 1000.0,
        };
        let driven = law
            .refresh_at_temperatures_with_drive(
                &before,
                cards.each_ref(),
                [HOT, COLD],
                HOT,
                drive,
                normal_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(
            approach.receipt().normal.kind,
            FiniteGapPlasticContinuationKind::MaterialRefresh
        );
        assert_eq!(
            driven.receipt().normal.kind,
            FiniteGapPlasticContinuationKind::MaterialRefreshWithDrive
        );
        assert_eq!(driven.receipt().boundary.boundary, drive);
        close(approach.receipt().normal.approach_m, old.approach_m());
        assert!(driven.receipt().normal.approach_m > old.approach_m());
        for (source, next, receipt) in [
            (
                approach.source(),
                approach.candidate(),
                &approach.receipt().normal,
            ),
            (
                driven.source(),
                driven.candidate(),
                &driven.receipt().normal,
            ),
        ] {
            assert_eq!(next.tangential_states(), before.tangential_states());
            assert_eq!(
                next.wear_state().volumes_m3(),
                before.wear_state().volumes_m3()
            );
            assert_eq!(
                next.wear_state().recession_m(),
                before.wear_state().recession_m()
            );
            let native = next.wear_state().contact_state().contact_state();
            assert!(native.plastic_depth_m()[12] > old.plastic_depth_m()[12]);
            assert!(receipt.response.plastic_work_j > 0.0);
            close(native.pressure_pa()[12], 0.5e6);
            assert_eq!(source.interface_temperature_k(), HOT);
        }
        close(driven.receipt().boundary.drive_force_n, 0.5e6 * 1.0e-8);
    });
}

#[test]
fn detached_temperature_change_binds_the_owner_without_inventing_a_pressure_sample() {
    with_cx(|cx, _| {
        let (cards, law) = thermal_source([1.0e5, 2.0e6], cx);
        let before = law.initial_state();
        let next = law
            .refresh_at_temperatures(&before, cards.each_ref(), [HOT, COLD], HOT, cx)
            .unwrap();
        assert_eq!(next.source().interface_temperature_k(), HOT);
        assert_ne!(next.source().identity(), law.identity());
        assert!(next.receipt().friction_sources.is_empty());
        assert!(
            next.receipt()
                .cell_source_indices
                .iter()
                .all(Option::is_none)
        );
        assert_eq!(
            next.candidate().tangential_states(),
            before.tangential_states()
        );
        assert_eq!(next.candidate().wear_state().volumes_m3(), [0.0; 2]);
        assert!(matches!(
            law.refresh_at_temperatures(&before, cards.each_ref(), [COLD; 2], f64::NAN, cx),
            Err(SourcePlasticSlidingError::Source(
                SourceThermalCoulombError::Invalid(_)
            ))
        ));
    });
}

#[test]
fn late_successor_pressure_domain_refusal_and_cancellation_leave_all_histories_intact() {
    with_cx(|cx, gate| {
        let (cards, law) = thermal_source([0.75e6, 2.0e6], cx);
        let loaded = law
            .advance_normal(&law.initial_state(), 5.0e-5, cx)
            .unwrap();
        let mut input = motion();
        input[12].released_increment_m = [1.0e-5, 0.0];
        let before = law
            .advance(&loaded.candidate, &input, 0.1, friction_controls(), cx)
            .unwrap()
            .candidate;
        let saved = before.clone();
        assert!(matches!(
            law.refresh_at_temperatures(&before, cards.each_ref(), [HOT, COLD], HOT, cx),
            Err(SourcePlasticSlidingError::Cell {
                index: 12,
                error: SourceThermalCoulombError::Material(_)
            })
        ));
        assert_eq!(before, saved);
        let retry = law
            .refresh_at_temperatures(&before, cards.each_ref(), [COLD; 2], HOT, cx)
            .unwrap();
        close(
            retry.receipt().friction_sources[0]
                .state_coordinate("p")
                .unwrap(),
            before
                .wear_state()
                .contact_state()
                .contact_state()
                .pressure_pa()[12],
        );
        assert_eq!(
            retry.candidate().tangential_states(),
            before.tangential_states()
        );
        gate.request();
        assert!(
            law.refresh_at_temperatures(&before, cards.each_ref(), [COLD; 2], HOT, cx)
                .is_err()
        );
        assert_eq!(before, saved);
    });
}
