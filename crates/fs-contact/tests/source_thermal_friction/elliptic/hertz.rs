//! G1/G3/G4: true source-owned pressure axes and directional friction together.
use super::*;
use fs_contact::interface_binding::NormalContactModelSelection;
use fs_contact::interface_binding::thermal_friction::elliptic::{
    SourceFrictionConditions, SourceThermalEllipticState,
    hertz::{SourceEllipticHertzError, SourceEllipticHertzInput},
};
use fs_contact::interface_binding::thermal_friction::{FrictionQueryAxis, FrictionQueryAxisKind};
use fs_contact::normal_patch::source::state::PointContactMaterialCards;
use fs_contact::normal_patch::source::{PointContactGeometry, SourcePointContact};
use fs_tribo::sliding_patch::elliptic::{
    EllipticHertzPatchParameters, EllipticHertzQuadrature, EllipticHertzSlidingPatch,
};
use fs_tribo::tangential::elliptic::EllipticCoulomb;
use fs_tribo::tangential::transport::{ObjectiveTangentialLaw, TangentialHistoryRotation};

#[path = "../../source_point/common.rs"]
mod normal_fixture;

fn fixture(
    bulk_temperatures: [f64; 2],
    texture: [f64; 2],
    speed_dependent: bool,
    rank_one: bool,
) -> (SourceThermalEllipticCoulomb, SourcePointContact) {
    let a = normal_fixture::solid("A", 1000.0, 1e7);
    let b = normal_fixture::solid("B", 2000.0, 2e7);
    let base = normal_fixture::bind(&a, &b, Some(0.4), 0.0);
    let original = base.normal_state().elastic().interface().card();
    let mut claims = original.claims().clone();
    for (name, axis, values) in [
        (STATIC, "Ta", [0.8, 1.0]),
        (STATIC_Y, "Tb", [0.9, 1.1]),
        (
            KINETIC,
            if speed_dependent { "speed" } else { "Ta" },
            [0.2, 0.6],
        ),
        (
            KINETIC_Y,
            "Tb",
            if rank_one { [0.0; 2] } else { [0.1, 0.5] },
        ),
    ] {
        let mut property = curve(name, values);
        let (dims, abscissas) = if axis == "speed" {
            (Dims([1, 0, -1, 0, 0, 0]), [0.0, 2.0])
        } else {
            (Temperature::DIMS, [280.0, 320.0])
        };
        property.value = PropertyValue::Curve {
            abscissa: axis.into(),
            abscissa_dims: dims,
            knots: vec![(abscissas[0], values[0]), (abscissas[1], values[1])],
            dims: Dims::NONE,
        };
        property.validity = property.validity.with_quantity(
            axis,
            QuantitySpec::dimensional(dims),
            abscissas[0],
            abscissas[1],
        );
        claims.insert_claim(property).unwrap();
    }
    let card = InterfaceSystemCard::assemble(
        original.surface_a().clone(),
        original.surface_b().clone(),
        original.context().clone(),
        claims,
        original.models().to_vec(),
    )
    .unwrap();
    let p = source_point(speed_dependent);
    let principal = [[STATIC, KINETIC], [STATIC_Y, KINETIC_Y]].map(|names| {
        let mut owned = vec![
            FrictionQueryAxis {
                name: "Ta".into(),
                kind: FrictionQueryAxisKind::Temperature,
            },
            FrictionQueryAxis {
                name: "Tb".into(),
                kind: FrictionQueryAxisKind::Temperature,
            },
            FrictionQueryAxis {
                name: "path".into(),
                kind: FrictionQueryAxisKind::AccumulatedDistance,
            },
        ];
        if speed_dependent {
            owned.push(FrictionQueryAxis {
                name: "speed".into(),
                kind: FrictionQueryAxisKind::SlidingSpeed,
            });
        }
        SourceThermalCoulomb::resolve_with_query_axes(
            bind_dry_interface_system_card(&card, InputAuthority::SyntheticFixture).unwrap(),
            p.clone(),
            requirement(names[0]),
            requirement(names[1]),
            MaterialPropertySelection::SingleClaimOnly,
            axes(),
            owned,
        )
        .unwrap()
    });
    let friction =
        SourceThermalEllipticCoulomb::from_principal_sources(principal, texture).unwrap();
    let cards = PointContactMaterialCards {
        surface_a: &a,
        surface_b: &b,
        interface: &card,
        bulk_selection: [
            MaterialPropertySelection::SingleClaimOnly,
            MaterialPropertySelection::SingleClaimOnly,
        ],
        interface_selection: MaterialPropertySelection::SingleClaimOnly,
        model_selection: NormalContactModelSelection::SingleSupported,
        authority: InputAuthority::SyntheticFixture,
    };
    let normal = cards
        .resolve_at_temperatures(
            [&query(bulk_temperatures[0]), &query(bulk_temperatures[1])],
            &query(300.0),
            PointContactGeometry::EllipticParaboloid {
                maximum_curvature_m_inverse: 20.0,
                minimum_curvature_m_inverse: 5.0,
            },
            normal_fixture::envelope(),
        )
        .unwrap();
    (friction, normal)
}

fn source_point(speed_dependent: bool) -> QueryPoint {
    let p = query(300.0)
        .with_quantity("Ta", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
        .unwrap()
        .with_quantity("Tb", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
        .unwrap()
        .with_quantity(
            "path",
            QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0])),
            0.25,
        )
        .unwrap();
    if speed_dependent {
        p.with_quantity(
            "speed",
            QuantitySpec::dimensional(Dims([1, 0, -1, 0, 0, 0])),
            0.5,
        )
        .unwrap()
    } else {
        p
    }
}

fn conditions(speed_dependent: bool) -> SourceFrictionConditions {
    let mut p = source_point(speed_dependent)
        .with("Ta", 999.0)
        .unwrap()
        .with("Tb", 999.0)
        .unwrap()
        .with("path", 99.0)
        .unwrap();
    if speed_dependent {
        p = p.with("speed", 99.0).unwrap();
    }
    SourceFrictionConditions {
        query: p,
        contact_temperature_k: Some(300.0),
        nominal_pressure_pa: None,
    }
}

fn motion() -> SourceEllipticHertzInput<'static> {
    SourceEllipticHertzInput {
        normal_state_id: "elliptic-friction-sample",
        indentation_m: 2e-5,
        indentation_rate_m_s: 0.05,
        sample_duration_s: 1e-4,
        footprint_axis: [0.8, -0.6],
        slip_velocity_m_s: [0.03, -0.02],
        spin_rad_s: 12.0,
        surface_temperature_axes: [Some("Ta"), Some("Tb")],
    }
}

fn quadrature() -> EllipticHertzQuadrature {
    EllipticHertzQuadrature { angular_panels: 32 }
}

fn verify_sources(
    law: &SourceThermalEllipticCoulomb,
    state: &SourceThermalEllipticState,
    temperatures: [f64; 2],
) {
    for (index, source) in state.friction_sources().iter().enumerate() {
        assert_eq!(source.state_coordinate("Ta"), Some(temperatures[0]));
        assert_eq!(source.state_coordinate("Tb"), Some(temperatures[1]));
        assert_eq!(source.state_coordinate("T"), Some(300.0));
        assert_eq!(source.state_coordinate("path"), Some(0.25));
        for property in source.properties() {
            let original = law.principal_sources()[index]
                .initial_source()
                .property(property.requirement().name())
                .unwrap();
            assert_eq!(
                property.answer().receipt.selected,
                original.answer().receipt.selected
            );
            law.principal_sources()[index]
                .interface()
                .card()
                .claims()
                .verify_receipt(&property.answer().receipt)
                .unwrap();
        }
    }
}

#[test]
fn original_normal_axes_and_ordered_temperatures_drive_the_actual_wrench() {
    let texture = [0.6, 0.8];
    for (temperatures, mu) in [([290.0, 310.0], [0.3, 0.4]), ([310.0, 290.0], [0.5, 0.2])] {
        let (law, normal) = fixture(temperatures, texture, false, false);
        let state = law.initial_state();
        let input = motion();
        let report = law
            .evaluate_hertz_patch(
                &state,
                &normal,
                input,
                &conditions(false),
                quadrature(),
                &mut || false,
            )
            .unwrap();
        let expected_normal = normal
            .evaluate(
                input.normal_state_id,
                input.indentation_m,
                input.indentation_rate_m_s,
                input.sample_duration_s,
            )
            .unwrap();
        assert_eq!(report.normal, expected_normal);
        assert_eq!(report.normal_model_identity, normal.model().identity());
        let axes = expected_normal.elliptic_patch_axes.unwrap();
        assert!(axes.semi_major_axis_m > axes.semi_minor_axis_m);
        let native = EllipticHertzSlidingPatch::new(
            EllipticHertzPatchParameters {
                normal_force_n: expected_normal.normal_force_n,
                semi_axes_m: [axes.semi_major_axis_m, axes.semi_minor_axis_m],
                footprint_axis: input.footprint_axis,
            },
            EllipticCoulomb::new([1.1; 2], mu, texture).unwrap(),
        )
        .unwrap()
        .evaluate(input.slip_velocity_m_s, input.spin_rad_s, quadrature())
        .unwrap();
        for i in 0..2 {
            close(report.wrench.force_n[i], native.force_n[i]);
        }
        close(report.wrench.spin_moment_nm, native.spin_moment_nm);
        close(report.wrench.dissipated_power_w, native.dissipated_power_w);
        close(report.wrench.power_residual_w, 0.0);
        verify_sources(&law, &report.state, temperatures);
        assert_eq!(
            report.state.constitutive_state().slip_m(),
            state.constitutive_state().slip_m()
        );
        assert_eq!(
            report.state.constitutive_state().principal_axis(),
            state.constitutive_state().principal_axis()
        );
        assert_eq!(
            report.state.constitutive_state().contact_temperature_k(),
            state.constitutive_state().contact_temperature_k()
        );
        assert_eq!(state, law.initial_state());
        // Normal hysteretic loss and tangential heat have separate owners.
        assert!(report.normal.dissipated_power_w > 0.0);
        close(
            report.normal.irreversible_work_j,
            report.normal.dissipated_power_w * input.sample_duration_s,
        );
    }
}

#[test]
fn source_rank_one_spin_uses_true_minor_axis_and_hertz_load_scaling() {
    let (law, normal) = fixture([290.0, 310.0], [1.0, 0.0], false, true);
    let state = law.initial_state();
    let mut input = motion();
    input.footprint_axis = [1.0, 0.0];
    input.slip_velocity_m_s = [0.0; 2];
    input.indentation_rate_m_s = 0.0;
    let mut previous = None;
    for indentation in [1e-5, 4e-5] {
        input.indentation_m = indentation;
        let r = law
            .evaluate_hertz_patch(
                &state,
                &normal,
                input,
                &conditions(false),
                quadrature(),
                &mut || false,
            )
            .unwrap();
        let b = r.normal.elliptic_patch_axes.unwrap().semi_minor_axis_m;
        let expected = -3.0 * r.normal.normal_force_n * 0.3 * b / 8.0;
        normal_fixture::close(r.wrench.spin_moment_nm, expected);
        close(r.wrench.force_n[0], 0.0);
        close(r.wrench.force_n[1], 0.0);
        normal_fixture::close(r.wrench.dissipated_power_w, -expected * input.spin_rad_s);
        assert!(
            (expected / (-3.0 * r.normal.normal_force_n * 0.3 * r.normal.patch_radius_m / 8.0)
                - 1.0)
                .abs()
                > 0.1
        );
        if let Some((force, moment)) = previous {
            normal_fixture::close(r.normal.normal_force_n, 8.0 * force);
            normal_fixture::close(r.wrench.spin_moment_nm, 16.0 * moment);
        }
        previous = Some((r.normal.normal_force_n, r.wrench.spin_moment_nm));
    }
}

#[test]
fn pure_translation_queries_actual_speed_and_retains_original_path_without_fake_spin_sampling() {
    let (law, normal) = fixture([290.0, 310.0], [1.0, 0.0], true, false);
    let state = law.initial_state();
    let saved = state.clone();
    let mut input = motion();
    input.spin_rad_s = 0.0;
    input.slip_velocity_m_s = [0.6, 0.8];
    let r = law
        .evaluate_hertz_patch(
            &state,
            &normal,
            input,
            &conditions(true),
            quadrature(),
            &mut || false,
        )
        .unwrap();
    // mu_x=0.2+0.2*|u|=0.4 and mu_y(Tb=310)=0.4, independently.
    close(r.wrench.force_n[0], -0.4 * r.normal.normal_force_n * 0.6);
    close(r.wrench.force_n[1], -0.4 * r.normal.normal_force_n * 0.8);
    assert_eq!(r.wrench.spin_moment_nm, 0.0);
    close(r.wrench.dissipated_power_w, 0.4 * r.normal.normal_force_n);
    for source in r.state.friction_sources() {
        assert_eq!(source.state_coordinate("speed"), Some(1.0));
    }
    verify_sources(&law, &r.state, [290.0, 310.0]);
    assert_eq!(
        r.state.constitutive_state().slip_m(),
        state.constitutive_state().slip_m()
    );
    assert_eq!(
        r.state.constitutive_state().principal_axis(),
        state.constitutive_state().principal_axis()
    );
    assert_eq!(
        r.state.constitutive_state().contact_temperature_k(),
        state.constitutive_state().contact_temperature_k()
    );
    input.spin_rad_s = 1e-200;
    assert!(matches!(
        law.evaluate_hertz_patch(
            &state,
            &normal,
            input,
            &conditions(true),
            quadrature(),
            &mut || false
        ),
        Err(SourceEllipticHertzError::Invalid(
            "spinning speed-dependent friction requires resolved local speed queries"
        ))
    ));
    assert_eq!(state, saved);
}

#[test]
fn transported_texture_and_independent_footprint_rotate_the_whole_native_wrench() {
    let (law, normal) = fixture([290.0, 310.0], [1.0, 0.0], false, false);
    let initial = law.initial_state();
    let input = motion();
    let reference = law
        .evaluate_hertz_patch(
            &initial,
            &normal,
            input,
            &conditions(false),
            quadrature(),
            &mut || false,
        )
        .unwrap();
    let old = TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap();
    let new = old.rotated(0.7).unwrap();
    let rotation = TangentialHistoryRotation::between(old, new, 0.0).unwrap();
    let transported = law.transport_state(&initial, rotation).unwrap();
    let rotate = |v: [f64; 2]| {
        let c = 0.7_f64.cos();
        let s = 0.7_f64.sin();
        [c * v[0] + s * v[1], -s * v[0] + c * v[1]]
    };
    let mut moved = input;
    moved.footprint_axis = rotate(input.footprint_axis);
    moved.slip_velocity_m_s = rotate(input.slip_velocity_m_s);
    let report = law
        .evaluate_hertz_patch(
            &transported,
            &normal,
            moved,
            &conditions(false),
            quadrature(),
            &mut || false,
        )
        .unwrap();
    let expected = rotate(reference.wrench.force_n);
    for i in 0..2 {
        close(report.wrench.force_n[i], expected[i]);
    }
    close(
        report.wrench.spin_moment_nm,
        reference.wrench.spin_moment_nm,
    );
    close(
        report.wrench.dissipated_power_w,
        reference.wrench.dissipated_power_w,
    );
    assert_eq!(
        report.state.constitutive_state().principal_axis(),
        transported.constitutive_state().principal_axis()
    );
    assert_ne!(
        report.state.constitutive_state().principal_axis(),
        initial.constitutive_state().principal_axis()
    );
}

#[test]
fn inactive_contact_keeps_zero_footprint_and_validates_sources_before_returning() {
    let (law, normal) = fixture([290.0, 310.0], [1.0, 0.0], false, false);
    let state = law.initial_state();
    let mut input = motion();
    input.indentation_m = 0.0;
    let report = law
        .evaluate_hertz_patch(
            &state,
            &normal,
            input,
            &conditions(false),
            quadrature(),
            &mut || false,
        )
        .unwrap();
    assert!(report.patch.is_none());
    assert_eq!(report.normal.normal_force_n, 0.0);
    assert_eq!(report.wrench.force_n, [0.0; 2]);
    assert_eq!(report.wrench.spin_moment_nm, 0.0);
    assert_eq!(report.wrench.dissipated_power_w, 0.0);
    assert!(!report.wrench.at_rest);
    verify_sources(&law, &report.state, [290.0, 310.0]);
    assert!(
        law.evaluate_hertz_patch(
            &state,
            &normal,
            input,
            &conditions(false),
            EllipticHertzQuadrature { angular_panels: 0 },
            &mut || false
        )
        .is_err()
    );
}

#[test]
fn source_geometry_and_late_cancellation_refusals_leave_the_original_checkpoint_retryable() {
    let (law, normal) = fixture([290.0, 310.0], [1.0, 0.0], false, false);
    let state = law.initial_state();
    let saved = state.clone();
    let baseline = law
        .evaluate_hertz_patch(
            &state,
            &normal,
            motion(),
            &conditions(false),
            quadrature(),
            &mut || false,
        )
        .unwrap();
    let (foreign_law, foreign_normal) = fixture([290.0, 310.0], [0.0, 1.0], false, true);
    assert!(
        law.evaluate_hertz_patch(
            &state,
            &foreign_normal,
            motion(),
            &conditions(false),
            quadrature(),
            &mut || false,
        )
        .is_err()
    );
    assert!(
        law.evaluate_hertz_patch(
            &foreign_law.initial_state(),
            &normal,
            motion(),
            &conditions(false),
            quadrature(),
            &mut || false,
        )
        .is_err()
    );
    for (indentation, names) in [
        (1.01e-4, [Some("Ta"), Some("Tb")]),
        (2e-5, [Some("Ta"), Some("Ta")]),
        (2e-5, [Some("Ta"), None]),
    ] {
        let mut input = motion();
        input.indentation_m = indentation;
        input.surface_temperature_axes = names;
        assert!(
            law.evaluate_hertz_patch(
                &state,
                &normal,
                input,
                &conditions(false),
                quadrature(),
                &mut || false
            )
            .is_err()
        );
    }
    let mut wrong = conditions(false);
    wrong.contact_temperature_k = Some(301.0);
    assert!(
        law.evaluate_hertz_patch(&state, &normal, motion(), &wrong, quadrature(), &mut || {
            false
        })
        .is_err()
    );
    // Missing Ta cannot be concealed by filling its normal-source temperature.
    wrong = conditions(false);
    wrong.query = query(300.0)
        .with_quantity("Tb", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
        .unwrap()
        .with_quantity(
            "path",
            QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0])),
            0.25,
        )
        .unwrap();
    assert!(
        law.evaluate_hertz_patch(&state, &normal, motion(), &wrong, quadrature(), &mut || {
            false
        })
        .is_err()
    );
    for limit in [1, 4, 12] {
        let mut polls = 0;
        let result = law.evaluate_hertz_patch(
            &state,
            &normal,
            motion(),
            &conditions(false),
            quadrature(),
            &mut || {
                polls += 1;
                polls >= limit
            },
        );
        assert!(matches!(result, Err(SourceEllipticHertzError::Cancelled)));
    }
    assert_eq!(state, saved);
    assert_eq!(
        law.evaluate_hertz_patch(
            &state,
            &normal,
            motion(),
            &conditions(false),
            quadrature(),
            &mut || false
        )
        .unwrap(),
        baseline
    );
}
