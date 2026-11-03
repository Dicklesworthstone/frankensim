//! Native dense wetting with the same source law and SI conventions as 3-D.
use super::*;
use fs_lbm::freesurface::{ContactModel, FreeSurface, PlanarWall2};
use fs_lbm::material::free_surface::wetting::CoxVoinovLaw;
use fs_lbm::material::free_surface::{MaterialFreeSurfaceDomain2, MaterialFreeSurfaceWetting2};

fn domain2() -> MaterialFreeSurfaceDomain2 {
    MaterialFreeSurfaceDomain2 {
        cells: [20; 2],
        gravity_m_s2: [0.0; 2],
        contact: ContactModel::Neutral,
        max_cells: 400,
    }
}

fn shape(wall: PlanarWall2, x: usize, y: usize) -> bool {
    let (normal, tangent) = match wall {
        PlanarWall2::XMin => (x, y),
        PlanarWall2::XMax => (19 - x, y),
        PlanarWall2::YMin => (y, x),
        PlanarWall2::YMax => (19 - y, x),
    };
    (normal * normal) as f64 + (tangent as f64 - 10.0).powi(2) < 36.0
}

fn flow2<'s, 'm>(
    bulk: &'s ResolvedLiquidSurfaceState,
    card: &'m InterfaceSystemCard,
    wall: PlanarWall2,
    dt: f64,
) -> MaterialFreeSurfaceWetting2<'s, 'm> {
    plan_free_surface_from_material(bulk, 0.001, dt, 0.1)
        .unwrap()
        .free_surface2(domain2(), |x, y| shape(wall, x, y))
        .unwrap()
        .with_dynamic_planar_wetting(
            wall,
            card,
            MaterialPropertySelection::SingleClaimOnly,
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap()
}

fn same_state(a: &FreeSurface, b: &FreeSurface) {
    assert_eq!(a.grid.f, b.grid.f);
    assert_eq!(a.grid.flags, b.grid.flags);
    assert_eq!(a.mass, b.mass);
    assert_eq!(a.carry.to_bits(), b.carry.to_bits());
    assert_eq!(a.conversions.to_fluid, b.conversions.to_fluid);
    assert_eq!(a.conversions.to_gas, b.conversions.to_gas);
    assert_eq!(
        a.conversions.gas_to_interface,
        b.conversions.gas_to_interface
    );
    assert_eq!(
        a.conversions.fluid_to_interface,
        b.conversions.fluid_to_interface
    );
}

#[test]
fn dense_collapsed_bracket_replays_equilibrium_at_every_wall() {
    let bulk = liquid(0, 300.0, Some(0.05));
    let gate = CancelGate::new_clock_free();
    for (wall, wall3) in [
        (PlanarWall2::XMin, PlanarWall3::XMin),
        (PlanarWall2::XMax, PlanarWall3::XMax),
        (PlanarWall2::YMin, PlanarWall3::YMin),
        (PlanarWall2::YMax, PlanarWall3::YMax),
    ] {
        for angle in [PI / 3.0, FRAC_PI_2, 2.0 * PI / 3.0] {
            let card = dynamic_card(&bulk, [angle; 2], [angle; 2], 100.0, false);
            let boundary = wetting(&bulk, angle);
            let mut dynamic = flow2(&bulk, &card, wall, 0.0002);
            let mut equilibrium = plan(&bulk)
                .with_planar_wetting(wall3, &boundary)
                .unwrap()
                .free_surface2(domain2(), |x, y| shape(wall, x, y))
                .unwrap();
            let mass = dynamic.flow().ledger_mass_kg_per_m().unwrap();
            assert_eq!(dynamic.samples().len(), 18);
            for _ in 0..3 {
                dynamic.step_with_gate(&gate).unwrap();
                equilibrium.step();
                same_state(dynamic.flow().lattice(), equilibrium.lattice());
            }
            assert!((dynamic.flow().ledger_mass_kg_per_m().unwrap() - mass).abs() < mass * 1e-10);
        }
    }
}

#[test]
fn dense_source_speed_changes_real_capillary_momentum_with_one_step_lag() {
    let gate = CancelGate::new_clock_free();
    for sigma in [Some(0.05), None] {
        let bulk = liquid(0, 300.0, sigma);
        let angle = PI / 3.0;
        let constant = dynamic_card(&bulk, [angle; 2], [angle; 2], 100.0, false);
        let varying = dynamic_card(
            &bulk,
            [angle, PI / 6.0],
            [angle, 5.0 * PI / 6.0],
            100.0,
            false,
        );
        let mut a = flow2(&bulk, &constant, PlanarWall2::YMin, 0.0002);
        let mut b = flow2(&bulk, &varying, PlanarWall2::YMin, 0.0002);
        a.step_with_gate(&gate).unwrap();
        b.step_with_gate(&gate).unwrap();
        same_state(a.flow().lattice(), b.flow().lattice());
        if sigma.is_some() {
            assert!(
                b.samples()
                    .values()
                    .any(|s| s.applied_angle_radians != angle)
            );
        }
        for sample in b.samples().values() {
            assert_eq!(sample.material.card_identity(), varying.content_hash());
            assert_eq!(
                sample.material.state_coordinate("pressure"),
                Some(100_000.0)
            );
            assert_eq!(
                sample.material.state_coordinate(DYNAMIC_WETTING_SPEED_AXIS),
                Some(sample.signed_speed_m_s.abs())
            );
            for property in sample.material.properties() {
                assert_eq!(
                    property.answer().receipt.selected,
                    b.source()
                        .property(property.requirement().name())
                        .unwrap()
                        .answer()
                        .receipt
                        .selected
                );
            }
        }
        a.step_with_gate(&gate).unwrap();
        b.step_with_gate(&gate).unwrap();
        if sigma.is_some() {
            assert_ne!(a.flow().lattice().grid.f, b.flow().lattice().grid.f);
        } else {
            same_state(a.flow().lattice(), b.flow().lattice());
        }
    }
}

#[test]
fn dense_cox_prediction_uses_static_receipts_and_changes_the_native_operator() {
    let bulk = liquid(0, 300.0, Some(0.05));
    let gate = CancelGate::new_clock_free();
    let angle = 1.2;
    let source = static_card(&bulk, angle);
    let constant = dynamic_card(&bulk, [angle; 2], [angle; 2], 100.0, false);
    let mut cox = flow2(&bulk, &source, PlanarWall2::YMin, 0.00002)
        .with_cox_voinov(
            CoxVoinovLaw::try_new_with_viscosity_ratio(cox_controls(0.01), 0.1, &gate).unwrap(),
            &gate,
        )
        .unwrap();
    let mut a = flow2(&bulk, &constant, PlanarWall2::YMin, 0.00002);
    cox.step_with_gate(&gate).unwrap();
    a.step_with_gate(&gate).unwrap();
    same_state(cox.flow().lattice(), a.flow().lattice());
    assert!(cox.samples().values().any(|s| s.signed_speed_m_s != 0.0));
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
        assert_eq!(prediction.outer_to_liquid_viscosity_ratio, 0.1);
        if sample.regime != ContactLineRegime2::Pinned {
            assert_eq!(
                sample.applied_angle_radians,
                cox.cox_voinov_law()
                    .unwrap()
                    .apparent_angle_radians(angle, prediction.signed_capillary_number)
                    .unwrap()
            );
        }
    }
    cox.step_with_gate(&gate).unwrap();
    a.step_with_gate(&gate).unwrap();
    assert_ne!(cox.flow().lattice().grid.f, a.flow().lattice().grid.f);
    assert!(
        flow2(&bulk, &constant, PlanarWall2::YMin, 0.00002)
            .with_cox_voinov(
                CoxVoinovLaw::try_new(cox_controls(0.01), &gate).unwrap(),
                &gate
            )
            .is_err()
    );
}

#[test]
fn dense_source_exit_cox_limit_and_cancel_roll_back_complete_flow_and_history() {
    let bulk = liquid(0, 300.0, Some(0.05));
    let gate = CancelGate::new_clock_free();
    let angle = PI / 3.0;
    // A different source supports subsequent speed. The selected claim must
    // still fail, rather than silently adopting that unrelated measurement.
    let narrow = dynamic_card(&bulk, [angle; 2], [angle; 2], 1e-12, true);
    let mut flow = flow2(&bulk, &narrow, PlanarWall2::YMin, 0.0002);
    let before = flow.flow().lattice().clone();
    let samples = flow.samples().clone();
    assert!(flow.step_with_gate(&gate).is_err());
    same_state(flow.flow().lattice(), &before);
    assert_eq!(flow.samples(), &samples);

    let source = static_card(&bulk, angle);
    let mut limited = flow2(&bulk, &source, PlanarWall2::YMin, 0.00002)
        .with_cox_voinov(
            CoxVoinovLaw::try_new(cox_controls(1e-18), &gate).unwrap(),
            &gate,
        )
        .unwrap();
    let before = limited.flow().lattice().clone();
    let samples = limited.samples().clone();
    assert!(limited.step_with_gate(&gate).is_err());
    same_state(limited.flow().lattice(), &before);
    assert_eq!(limited.samples(), &samples);

    let supported = dynamic_card(&bulk, [angle; 2], [angle; 2], 100.0, false);
    let mut retry = flow2(&bulk, &supported, PlanarWall2::YMin, 0.0002);
    let mut reference = flow2(&bulk, &supported, PlanarWall2::YMin, 0.0002);
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(retry.step_with_gate(&cancelled).is_err());
    same_state(retry.flow().lattice(), reference.flow().lattice());
    assert_eq!(retry.samples(), reference.samples());
    retry.step_with_gate(&gate).unwrap();
    reference.step_with_gate(&gate).unwrap();
    same_state(retry.flow().lattice(), reference.flow().lattice());
    assert_eq!(retry.samples(), reference.samples());
}

#[test]
fn g1_dense_material_momentum_step_replays_native_dynamics_in_physical_units() {
    let bulk = liquid(0, 300.0, None);
    let gate = CancelGate::new_clock_free();
    for (dx, dt) in [(0.001, 0.0002), (0.002, 0.0005)] {
        let mut domain = domain2();
        domain.gravity_m_s2 = [0.5, -1.25];
        let mut physical = plan_free_surface_from_material(&bulk, dx, dt, 0.1)
            .unwrap()
            .free_surface2(domain, |x, y| shape(PlanarWall2::YMin, x, y))
            .unwrap();
        let mut native = physical.lattice().clone();
        let mass_unit = bulk.density_kg_m3() * dx * dx;
        let velocity_unit = dx / dt;
        for _ in 0..3 {
            let expected = native.step_with_momentum(&gate).unwrap();
            let actual = physical.step_with_momentum(&gate).unwrap();
            assert_eq!(actual.native, expected);
            same_state(physical.lattice(), &native);
            assert_eq!(
                actual.redistributed_mass_kg_per_m,
                expected.redistributed_mass * mass_unit
            );
            assert_eq!(
                actual.absolute_mass_defect_kg_per_m,
                expected.absolute_mass_defect * mass_unit
            );
            for d in 0..2 {
                assert_eq!(
                    actual.absolute_momentum_defect_kg_per_s[d],
                    (expected.absolute_momentum_defect[d] * mass_unit) * velocity_unit
                );
            }
            assert_eq!(
                actual.kinetic_energy_change_j_per_m,
                ((expected.kinetic_energy_change * mass_unit) * velocity_unit) * velocity_unit
            );
            assert_eq!(physical.carry_mass_kg_per_m().unwrap(), 0.0);
        }
        let before = physical.lattice().clone();
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(matches!(
            physical.step_with_momentum(&cancelled),
            Err(
                fs_lbm::material::free_surface::MaterialFreeSurfaceMomentumError2::Native(
                    fs_lbm::freesurface::FreeSurfaceMomentumError2::Cancelled
                )
            )
        ));
        same_state(physical.lattice(), &before);
    }
}

#[test]
fn g4_dense_momentum_wetting_keeps_source_refusal_and_retry_atomic() {
    use fs_lbm::material::free_surface::MaterialFreeSurfaceMomentumError2;
    let bulk = liquid(0, 300.0, Some(0.05));
    let gate = CancelGate::new_clock_free();
    let angle = PI / 3.0;
    let narrow = dynamic_card(&bulk, [angle; 2], [angle; 2], 1e-12, true);
    let mut refused = flow2(&bulk, &narrow, PlanarWall2::YMin, 0.0002);
    let before = refused.flow().lattice().clone();
    let samples = refused.samples().clone();
    // The native proposal is valid; the selected material's speed domain is
    // what refuses. No alternate source or partly advanced liquid may escape.
    let mut native = before.clone();
    native.step_with_momentum(&gate).unwrap();
    assert!(matches!(
        refused.step_with_momentum(&gate),
        Err(MaterialFreeSurfaceMomentumError2::Wetting(_))
    ));
    same_state(refused.flow().lattice(), &before);
    assert_eq!(refused.samples(), &samples);
    assert!(refused.step_with_gate(&gate).is_err());
    same_state(refused.flow().lattice(), &before);
    assert_eq!(refused.samples(), &samples);

    let supported = dynamic_card(&bulk, [angle; 2], [angle; 2], 100.0, false);
    for wall in [
        PlanarWall2::XMin,
        PlanarWall2::XMax,
        PlanarWall2::YMin,
        PlanarWall2::YMax,
    ] {
        let mut retry = flow2(&bulk, &supported, wall, 0.0002);
        let mut reference = flow2(&bulk, &supported, wall, 0.0002);
        let mut ordinary = flow2(&bulk, &supported, wall, 0.0002);
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(matches!(
            retry.step_with_momentum(&cancelled),
            Err(MaterialFreeSurfaceMomentumError2::Native(
                fs_lbm::freesurface::FreeSurfaceMomentumError2::Cancelled
            ))
        ));
        same_state(retry.flow().lattice(), reference.flow().lattice());
        assert_eq!(retry.samples(), reference.samples());
        assert!(matches!(
            ordinary.step_with_gate(&cancelled),
            Err(
                fs_lbm::material::free_surface::wetting::FreeSurfaceWettingError::Flow(
                    fs_lbm::material::MaterialLbmError::FreeSurface(
                        fs_lbm::d3q19::freesurface3::FreeSurfaceError3::Cancelled
                    )
                )
            )
        ));
        same_state(ordinary.flow().lattice(), reference.flow().lattice());
        let a = retry.step_with_momentum(&gate).unwrap();
        let b = reference.step_with_momentum(&gate).unwrap();
        ordinary.step_with_gate(&gate).unwrap();
        assert_eq!(a, b);
        same_state(retry.flow().lattice(), reference.flow().lattice());
        assert_eq!(retry.samples(), reference.samples());
        same_state(ordinary.flow().lattice(), reference.flow().lattice());
        assert_eq!(ordinary.samples(), reference.samples());
    }
}

#[test]
fn g3_dense_default_step_carries_momentum_through_live_conversions() {
    let bulk = liquid(0, 300.0, None);
    let gate = CancelGate::new_clock_free();
    let make = || {
        let mut domain = domain2();
        domain.gravity_m_s2 = [0.0, -9.81];
        plan_free_surface_from_material(&bulk, 0.001, 0.0004, 0.2)
            .unwrap()
            .free_surface2(domain, |x, y| {
                (x as f64 - 10.0).powi(2) + (y as f64 - 12.0).powi(2) < 16.0
            })
            .unwrap()
    };
    let mut ordinary = make();
    let mut explicit = make();
    let mut routes = 0;
    for _ in 0..41 {
        routes += explicit.step_with_momentum(&gate).unwrap().native.routes;
        ordinary.step();
        same_state(ordinary.lattice(), explicit.lattice());
    }
    assert!(
        routes > 0,
        "fixture must exercise actual conversion momentum"
    );
    let momentum = ordinary.resolved_linear_momentum_kg_per_s().unwrap();
    let mut sampled = [0.0; 2];
    for y in 0..20 {
        for x in 0..20 {
            if let Some(cell) = ordinary.sample_at_si(x, y).unwrap() {
                let mass_per_depth = cell.density_kg_m3 * cell.fill_fraction * 0.001_f64.powi(2);
                for d in 0..2 {
                    sampled[d] += mass_per_depth * cell.velocity_m_s[d];
                }
            }
        }
    }
    assert!(momentum[1] < 0.0);
    for d in 0..2 {
        assert!((momentum[d] - sampled[d]).abs() < 2e-14);
    }
}
