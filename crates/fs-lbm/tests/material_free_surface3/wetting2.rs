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
