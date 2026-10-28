//! Synthetic source/VOF integration, not contact-line calibration observations.

use super::*;
use fs_lbm::freesurface::{
    ADVANCING_CONTACT_ANGLE_PROPERTY, ContactLineRegime2, DYNAMIC_WETTING_SPEED_AXIS,
    DYNAMIC_WETTING_SPEED_DIMS, RECEDING_CONTACT_ANGLE_PROPERTY,
};
use fs_lbm::material::free_surface::wetting::{
    FreeSurfaceWettingControls, MaterialFreeSurfaceWetting3,
};

fn controls() -> FreeSurfaceWettingControls {
    FreeSurfaceWettingControls {
        pinning_speed_m_s: 0.0,
        maximum_wall_cells: 256,
    }
}

fn dynamic_card(
    bulk: &ResolvedLiquidSurfaceState,
    receding: [f64; 2],
    advancing: [f64; 2],
    maximum_speed: f64,
    include_fallback: bool,
) -> InterfaceSystemCard {
    let base = interface_card(bulk, None, Dims::NONE, false, None);
    let mut claims = ClaimSet::new();
    for (name, angles) in [
        (RECEDING_CONTACT_ANGLE_PROPERTY, receding),
        (ADVANCING_CONTACT_ANGLE_PROPERTY, advancing),
    ] {
        let mut source = claim(name, angles[0], Dims::NONE);
        source.value = PropertyValue::Curve {
            abscissa: DYNAMIC_WETTING_SPEED_AXIS.into(),
            abscissa_dims: DYNAMIC_WETTING_SPEED_DIMS,
            knots: vec![(0.0, angles[0]), (maximum_speed, angles[1])],
            dims: Dims::NONE,
        };
        source.interpolation = InterpolationPolicy::LinearInside;
        source.validity = source.validity.with_quantity(
            DYNAMIC_WETTING_SPEED_AXIS,
            QuantitySpec::dimensional(DYNAMIC_WETTING_SPEED_DIMS),
            0.0,
            maximum_speed,
        );
        claims.insert_claim(source).unwrap();
        if include_fallback {
            let mut fallback = claim(name, angles[0], Dims::NONE);
            fallback.validity = fallback.validity.with_quantity(
                DYNAMIC_WETTING_SPEED_AXIS,
                QuantitySpec::dimensional(DYNAMIC_WETTING_SPEED_DIMS),
                maximum_speed * 2.0,
                100.0,
            );
            claims.insert_claim(fallback).unwrap();
        }
    }
    InterfaceSystemCard::assemble(
        base.surface_a().clone(),
        base.surface_b().clone(),
        base.context().clone(),
        claims,
        vec![],
    )
    .unwrap()
}

fn dynamic_flow<'s, 'm>(
    bulk: &'s ResolvedLiquidSurfaceState,
    card: &'m InterfaceSystemCard,
) -> MaterialFreeSurfaceWetting3<'s, 'm> {
    plan(bulk)
        .free_surface3(domain(), drop_shape)
        .unwrap()
        .with_dynamic_planar_wetting(
            PlanarWall3::ZMin,
            card,
            MaterialPropertySelection::SingleClaimOnly,
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap()
}

#[test]
fn g1_collapsed_hysteresis_reuses_actual_equilibrium_capillary_motion() {
    let bulk = liquid(0, 300.0, Some(0.05));
    let angle = PI / 3.0;
    let card = dynamic_card(&bulk, [angle; 2], [angle; 2], 100.0, false);
    let source = wetting(&bulk, angle);
    let mut dynamic = dynamic_flow(&bulk, &card);
    let mut equilibrium = plan(&bulk)
        .with_planar_wetting(PlanarWall3::ZMin, &source)
        .unwrap()
        .free_surface3(domain(), drop_shape)
        .unwrap();
    let initial_mass = dynamic.flow().ledger_mass_kg().unwrap();
    let gate = CancelGate::new_clock_free();
    assert!(
        dynamic
            .samples()
            .values()
            .all(|s| s.regime == ContactLineRegime2::Pinned)
    );
    for _ in 0..2 {
        dynamic.step_with_gate(&gate).unwrap();
        equilibrium.step_with_gate(&gate).unwrap();
        assert_eq!(
            dynamic.flow().lattice().grid().state_bits(),
            equilibrium.lattice().grid().state_bits()
        );
        assert_eq!(
            dynamic.flow().ledger_mass_kg().unwrap().to_bits(),
            equilibrium.ledger_mass_kg().unwrap().to_bits()
        );
    }
    assert!((dynamic.flow().ledger_mass_kg().unwrap() - initial_mass).abs() < initial_mass * 1e-10);
    assert!(
        dynamic
            .samples()
            .values()
            .any(|s| s.signed_speed_m_s != 0.0)
    );
}

#[test]
fn g3_speed_dependent_interface_data_changes_native_populations_with_one_step_lag() {
    let bulk = liquid(0, 300.0, Some(0.05));
    let start = PI / 3.0;
    let constant = dynamic_card(&bulk, [start; 2], [start; 2], 100.0, false);
    let dependent = dynamic_card(
        &bulk,
        [start, PI / 6.0],
        [start, 5.0 * PI / 6.0],
        100.0,
        false,
    );
    let mut a = dynamic_flow(&bulk, &constant);
    let mut b = dynamic_flow(&bulk, &dependent);
    let gate = CancelGate::new_clock_free();
    a.step_with_gate(&gate).unwrap();
    b.step_with_gate(&gate).unwrap();
    // Same initial boundary and evolved speed; the rate-law responses differ.
    assert_eq!(
        a.flow().lattice().grid().state_bits(),
        b.flow().lattice().grid().state_bits()
    );
    assert!(
        b.samples()
            .values()
            .any(|s| s.applied_angle_radians != start)
    );
    for sample in b.samples().values() {
        assert_eq!(sample.material.card_identity(), dependent.content_hash());
        for property in sample.material.properties() {
            assert_eq!(
                Some(property.answer().receipt.selected),
                b.source()
                    .property(property.requirement().name())
                    .map(|p| p.answer().receipt.selected)
            );
        }
    }
    a.step_with_gate(&gate).unwrap();
    b.step_with_gate(&gate).unwrap();
    assert_ne!(
        a.flow().lattice().grid().state_bits(),
        b.flow().lattice().grid().state_bits()
    );
}

#[test]
fn g3_omitted_capillarity_keeps_wetting_data_out_of_native_momentum() {
    let bulk = liquid(0, 300.0, None);
    let a_card = dynamic_card(&bulk, [PI / 3.0; 2], [PI / 3.0; 2], 100.0, false);
    let b_card = dynamic_card(
        &bulk,
        [2.0 * PI / 3.0; 2],
        [2.0 * PI / 3.0; 2],
        100.0,
        false,
    );
    let mut a = dynamic_flow(&bulk, &a_card);
    let mut b = dynamic_flow(&bulk, &b_card);
    let gate = CancelGate::new_clock_free();
    a.step_with_gate(&gate).unwrap();
    b.step_with_gate(&gate).unwrap();
    assert_eq!(
        a.flow().lattice().grid().state_bits(),
        b.flow().lattice().grid().state_bits()
    );
}

#[test]
fn g3_typed_held_pressure_changes_sourced_angles_and_actual_capillary_motion() {
    let low = liquid_at_pressure(0, 300.0, Some(0.05), 95_000.0);
    let high = liquid_at_pressure(0, 300.0, Some(0.05), 105_000.0);
    assert_eq!(
        low.material().card_identity(),
        high.material().card_identity()
    );
    let base = interface_card(&low, None, Dims::NONE, false, None);
    let mut claims = ClaimSet::new();
    for name in [
        RECEDING_CONTACT_ANGLE_PROPERTY,
        ADVANCING_CONTACT_ANGLE_PROPERTY,
    ] {
        let mut source = claim(name, PI / 3.0, Dims::NONE);
        source.value = PropertyValue::Curve {
            abscissa: "pressure".into(),
            abscissa_dims: fs_qty::Pressure::DIMS,
            knots: vec![(90_000.0, PI / 3.0), (110_000.0, 2.0 * PI / 3.0)],
            dims: Dims::NONE,
        };
        source.interpolation = InterpolationPolicy::LinearInside;
        source.validity = source.validity.with_quantity(
            DYNAMIC_WETTING_SPEED_AXIS,
            QuantitySpec::dimensional(DYNAMIC_WETTING_SPEED_DIMS),
            0.0,
            100.0,
        );
        claims.insert_claim(source).unwrap();
    }
    let card = InterfaceSystemCard::assemble(
        base.surface_a().clone(),
        base.surface_b().clone(),
        base.context().clone(),
        claims,
        vec![],
    )
    .unwrap();
    let mut a = dynamic_flow(&low, &card);
    let mut b = dynamic_flow(&high, &card);
    for (flow, pressure, expected) in [
        (&a, 95_000.0, 5.0 * PI / 12.0),
        (&b, 105_000.0, 7.0 * PI / 12.0),
    ] {
        for sample in flow.samples().values() {
            assert!((sample.applied_angle_radians - expected).abs() < 2e-15);
            assert_eq!(sample.material.state_coordinate("pressure"), Some(pressure));
            assert_eq!(
                sample.material.axis_quantities().get("pressure"),
                Some(&QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::Pressure,
                    ValueForm::Static
                )))
            );
            assert_eq!(sample.material.card_identity(), card.content_hash());
        }
    }
    let gate = CancelGate::new_clock_free();
    a.step_with_gate(&gate).unwrap();
    b.step_with_gate(&gate).unwrap();
    assert_ne!(
        a.flow().lattice().grid().state_bits(),
        b.flow().lattice().grid().state_bits()
    );
    for (flow, pressure) in [(&a, 95_000.0), (&b, 105_000.0)] {
        assert!(
            flow.samples()
                .values()
                .all(|sample| sample.material.state_coordinate("pressure") == Some(pressure))
        );
    }
}

#[test]
fn g4_source_exit_and_cancellation_preserve_populations_and_contact_history() {
    let bulk = liquid(0, 300.0, Some(0.05));
    let angle = PI / 3.0;
    // A different claim DOES cover moving speeds; original selections must
    // still refuse after the actual native motion leaves their narrow support.
    let card = dynamic_card(&bulk, [angle; 2], [angle; 2], 1e-12, true);
    let mut flow = dynamic_flow(&bulk, &card);
    let bits = flow.flow().lattice().grid().state_bits();
    let mass = flow.flow().ledger_mass_kg().unwrap().to_bits();
    let conversions = flow.flow().lattice().conversions();
    let samples = flow.samples().clone();
    assert!(flow.step_with_gate(&CancelGate::new_clock_free()).is_err());
    assert_eq!(flow.flow().lattice().grid().state_bits(), bits);
    assert_eq!(flow.flow().ledger_mass_kg().unwrap().to_bits(), mass);
    assert_eq!(flow.flow().lattice().conversions(), conversions);
    assert_eq!(flow.samples(), &samples);

    let supported = dynamic_card(&bulk, [angle; 2], [angle; 2], 100.0, false);
    let mut retry = dynamic_flow(&bulk, &supported);
    let mut reference = dynamic_flow(&bulk, &supported);
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(retry.step_with_gate(&cancelled).is_err());
    assert_eq!(retry.samples(), reference.samples());
    let gate = CancelGate::new_clock_free();
    retry.step_with_gate(&gate).unwrap();
    reference.step_with_gate(&gate).unwrap();
    assert_eq!(
        retry.flow().lattice().grid().state_bits(),
        reference.flow().lattice().grid().state_bits()
    );
    assert_eq!(retry.samples(), reference.samples());
}

fn cox_controls(
    maximum_capillary_number: f64,
) -> fs_lbm::material::free_surface::wetting::CoxVoinovControls {
    fs_lbm::material::free_surface::wetting::CoxVoinovControls {
        microscopic_cutoff_m: 1e-9,
        observation_length_m: 0.001,
        maximum_capillary_number,
        minimum_angle_radians: 0.01,
        maximum_angle_radians: 3.0,
    }
}

fn static_card(bulk: &ResolvedLiquidSurfaceState, angle: f64) -> InterfaceSystemCard {
    let base = interface_card(bulk, None, Dims::NONE, false, None);
    let mut claims = ClaimSet::new();
    for name in [
        RECEDING_CONTACT_ANGLE_PROPERTY,
        ADVANCING_CONTACT_ANGLE_PROPERTY,
    ] {
        let mut endpoint = claim(name, angle, Dims::NONE);
        // Measurement support is stationary ONLY. Predicted dynamic angles
        // must not relabel these source receipts as measurements at speed U.
        endpoint.validity = endpoint.validity.with_quantity(
            DYNAMIC_WETTING_SPEED_AXIS,
            QuantitySpec::dimensional(DYNAMIC_WETTING_SPEED_DIMS),
            0.0,
            0.0,
        );
        claims.insert_claim(endpoint).unwrap();
    }
    InterfaceSystemCard::assemble(
        base.surface_a().clone(),
        base.surface_b().clone(),
        base.context().clone(),
        claims,
        vec![],
    )
    .unwrap()
}

fn slow_static_flow<'s, 'm>(
    bulk: &'s ResolvedLiquidSurfaceState,
    card: &'m InterfaceSystemCard,
) -> MaterialFreeSurfaceWetting3<'s, 'm> {
    plan_free_surface_from_material(bulk, 0.001, 0.00002, 0.1)
        .unwrap()
        .free_surface3(domain(), drop_shape)
        .unwrap()
        .with_dynamic_planar_wetting(
            PlanarWall3::ZMin,
            card,
            MaterialPropertySelection::SingleClaimOnly,
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap()
}

#[test]
fn g1_g3_cox_wetting_drives_native_flow_from_static_source_endpoints() {
    use fs_lbm::material::free_surface::wetting::CoxVoinovLaw;
    let bulk = liquid(0, 300.0, Some(0.05));
    let angle = 1.2;
    let static_data = static_card(&bulk, angle);
    let measured_rate_independent = dynamic_card(&bulk, [angle; 2], [angle; 2], 100.0, false);
    let gate = CancelGate::new_clock_free();
    let mut cox = slow_static_flow(&bulk, &static_data)
        .with_cox_voinov(
            CoxVoinovLaw::try_new(cox_controls(0.01), &gate).unwrap(),
            &gate,
        )
        .unwrap();
    let mut unchanged_angle = slow_static_flow(&bulk, &measured_rate_independent);
    let mass = cox.flow().ledger_mass_kg().unwrap();
    cox.step_with_gate(&gate).unwrap();
    unchanged_angle.step_with_gate(&gate).unwrap();
    assert_eq!(
        cox.flow().lattice().grid().state_bits(),
        unchanged_angle.flow().lattice().grid().state_bits()
    );
    assert!(
        cox.samples()
            .values()
            .any(|sample| sample.signed_speed_m_s != 0.0)
    );
    for sample in cox.samples().values() {
        assert_eq!(
            sample.material.state_coordinate(DYNAMIC_WETTING_SPEED_AXIS),
            Some(0.0)
        );
        let prediction = sample.hydrodynamic.unwrap();
        assert_eq!(
            prediction.dynamic_viscosity_pa_s,
            bulk.dynamic_viscosity_pa_s()
        );
        assert_eq!(prediction.surface_tension_n_m, 0.05);
        if sample.regime != ContactLineRegime2::Pinned {
            let expected = cox
                .cox_voinov_law()
                .unwrap()
                .apparent_angle_radians(angle, prediction.signed_capillary_number)
                .unwrap();
            assert_eq!(sample.applied_angle_radians, expected);
        }
    }
    cox.step_with_gate(&gate).unwrap();
    unchanged_angle.step_with_gate(&gate).unwrap();
    assert_ne!(
        cox.flow().lattice().grid().state_bits(),
        unchanged_angle.flow().lattice().grid().state_bits()
    );
    assert!((cox.flow().ledger_mass_kg().unwrap() - mass).abs() < mass * 1e-10);
}

#[test]
fn g0_g4_cox_refuses_double_rate_counting_and_preserves_state_on_model_exit() {
    use fs_lbm::material::free_surface::wetting::CoxVoinovLaw;
    let bulk = liquid(0, 300.0, Some(0.05));
    let gate = CancelGate::new_clock_free();
    let dynamic = dynamic_card(&bulk, [1.2, 1.0], [1.2, 1.4], 100.0, false);
    assert!(
        slow_static_flow(&bulk, &dynamic)
            .with_cox_voinov(
                CoxVoinovLaw::try_new(cox_controls(0.01), &gate).unwrap(),
                &gate
            )
            .is_err()
    );
    let static_data = static_card(&bulk, 1.2);
    let mut flow = slow_static_flow(&bulk, &static_data)
        .with_cox_voinov(
            CoxVoinovLaw::try_new(cox_controls(1e-18), &gate).unwrap(),
            &gate,
        )
        .unwrap();
    let bits = flow.flow().lattice().grid().state_bits();
    let mass = flow.flow().ledger_mass_kg().unwrap();
    let samples = flow.samples().clone();
    let conversions = flow.flow().lattice().conversions();
    assert!(flow.step_with_gate(&gate).is_err());
    assert_eq!(flow.flow().lattice().grid().state_bits(), bits);
    assert_eq!(flow.flow().ledger_mass_kg().unwrap(), mass);
    assert_eq!(flow.samples(), &samples);
    assert_eq!(flow.flow().lattice().conversions(), conversions);

    let neglected = liquid(0, 300.0, None);
    let no_sigma = static_card(&neglected, 1.2);
    assert!(
        slow_static_flow(&neglected, &no_sigma)
            .with_cox_voinov(
                CoxVoinovLaw::try_new(cox_controls(0.01), &gate).unwrap(),
                &gate
            )
            .is_err()
    );
}
