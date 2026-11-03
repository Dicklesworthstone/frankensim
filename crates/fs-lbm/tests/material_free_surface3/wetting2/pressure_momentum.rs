//! G1 source-unit/force-frame coupling and G4 flow/contact-history publication.
use super::*;
use fs_lbm::freesurface::{ConversionRedistribution2, FreeSurfacePressureError2};
use fs_lbm::material::free_surface::MaterialFreeSurfaceMomentumError2;

#[test]
fn g1_angular_pressure_conversion_uses_original_material_units_and_guo_velocity() {
    let bulk = liquid(0, 300.0, None);
    let gate = CancelGate::new_clock_free();
    for (dx, dt) in [(0.001, 0.0002), (0.002, 0.0005)] {
        let domain = MaterialFreeSurfaceDomain2 {
            gravity_m_s2: [0.5, -1.25],
            ..domain2()
        };
        let make = |reference| {
            let mut flow = plan_free_surface_from_material(&bulk, dx, dt, 0.1)
                .unwrap()
                .with_pressure_reference(reference)
                .unwrap()
                .free_surface2(domain, |x, y| shape(PlanarWall2::YMin, x, y))
                .unwrap();
            flow.set_conversion_redistribution(ConversionRedistribution2::ConservativeLocalAngular)
                .unwrap();
            flow
        };
        let mut physical = make(100_000.0);
        let mut native = physical.lattice().clone();
        let speed_unit = dx / dt;
        let pressure_unit = (bulk.density_kg_m3() * speed_unit) * speed_unit;
        let mass_unit = (bulk.density_kg_m3() * dx) * dx;
        let expected = native
            .step_with_gas_link_pressure_and_momentum(&gate, 0.05, |link| {
                Ok((0.125 * link.incoming_direction as f64) / pressure_unit)
            })
            .unwrap();
        let visits = std::cell::Cell::new(0);
        let actual = physical
            .step_with_gas_link_pressure_and_momentum(&gate, 0.05, |link| {
                visits.set(visits.get() + 1);
                assert_eq!(
                    link.liquid_center_m,
                    link.link.liquid_cell.map(|c| (c as f64 + 0.5) * dx)
                );
                assert_eq!(
                    link.gas_center_m,
                    link.link.gas_cell.map(|c| (c as f64 + 0.5) * dx)
                );
                for d in 0..2 {
                    assert_eq!(
                        link.midpoint_m[d],
                        f64::midpoint(link.liquid_center_m[d], link.gas_center_m[d])
                    );
                }
                100_000.0 + 0.125 * link.link.incoming_direction as f64
            })
            .unwrap();
        assert!(visits.get() > 0);
        assert_eq!(actual.native, expected);
        assert!(actual.native.absolute_angular_momentum_defect.is_some());
        assert_eq!(
            actual.absolute_angular_momentum_defect_kg_m_per_s,
            expected
                .absolute_angular_momentum_defect
                .map(|defect| ((defect * mass_unit) * speed_unit) * dx)
        );
        same_state(physical.lattice(), &native);
        assert_eq!(
            actual.redistributed_mass_kg_per_m,
            expected.redistributed_mass * mass_unit
        );
        assert_eq!(
            actual.kinetic_energy_change_j_per_m,
            ((expected.kinetic_energy_change * mass_unit) * speed_unit) * speed_unit
        );
        for d in 0..2 {
            assert_eq!(
                actual.absolute_momentum_defect_kg_per_s[d],
                (expected.absolute_momentum_defect[d] * mass_unit) * speed_unit
            );
        }
        let mut observed = 0;
        for y in 0..20 {
            for x in 0..20 {
                let Some(sample) = physical.sample_at_si(x, y).unwrap() else {
                    continue;
                };
                observed += 1;
                let i = native.grid.idx(x, y);
                let rho = native.grid.f[i].iter().sum::<f64>();
                for d in 0..2 {
                    let raw = native.grid.f[i]
                        .iter()
                        .zip(fs_lbm::E)
                        .map(|(f, e)| f * f64::from([e.0, e.1][d]))
                        .sum::<f64>();
                    let force = native.grid.g[d].mul_add(rho, native.grid.fext[i][d]);
                    let expected_velocity = speed_unit * ((raw + 0.5 * force) / rho);
                    assert!((sample.velocity_m_s[d] - expected_velocity).abs() < 2e-15);
                }
            }
        }
        assert!(observed > 0);

        let mut cell = make(100_000.0);
        let mut shifted_link = make(100_500.0);
        let a = cell
            .step_with_gas_pressure_and_momentum(&gate, 0.05, |_| 100_000.125)
            .unwrap();
        let b = shifted_link
            .step_with_gas_link_pressure_and_momentum(&gate, 0.05, |_| 100_500.125)
            .unwrap();
        assert_eq!(a, b);
        same_state(cell.lattice(), shifted_link.lattice());
    }
    let mut unbound = plan(&bulk)
        .free_surface2(domain2(), |x, y| shape(PlanarWall2::YMin, x, y))
        .unwrap();
    unbound
        .set_conversion_redistribution(ConversionRedistribution2::ConservativeLocalAngular)
        .unwrap();
    let before = unbound.lattice().clone();
    assert!(matches!(
        unbound.step_with_gas_pressure_and_momentum(&gate, 0.05, |_| 100_000.0),
        Err(MaterialFreeSurfaceMomentumError2::Units(_))
    ));
    same_state(unbound.lattice(), &before);
}

fn loaded_flow<'s, 'm>(
    bulk: &'s ResolvedLiquidSurfaceState,
    card: &'m InterfaceSystemCard,
) -> MaterialFreeSurfaceWetting2<'s, 'm> {
    let mut flow = plan_free_surface_from_material(bulk, 0.001, 0.0002, 0.1)
        .unwrap()
        .with_pressure_reference(100_000.0)
        .unwrap()
        .free_surface2(domain2(), |x, y| shape(PlanarWall2::YMin, x, y))
        .unwrap();
    flow.set_conversion_redistribution(ConversionRedistribution2::ConservativeLocalAngular)
        .unwrap();
    flow.with_dynamic_planar_wetting(
        PlanarWall2::YMin,
        card,
        MaterialPropertySelection::SingleClaimOnly,
        controls(),
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}

#[test]
fn g4_angular_pressure_conversion_and_sourced_wetting_commit_or_refuse_together() {
    let bulk = liquid(0, 300.0, Some(0.05));
    let gate = CancelGate::new_clock_free();
    let angle = PI / 3.0;
    let narrow = dynamic_card(&bulk, [angle; 2], [angle; 2], 1e-12, true);
    let mut refused = loaded_flow(&bulk, &narrow);
    assert_eq!(
        refused.flow().conversion_redistribution(),
        ConversionRedistribution2::ConservativeLocalAngular
    );
    let before = refused.flow().lattice().clone();
    let samples = refused.samples().clone();
    let mut native = before.clone();
    native
        .step_with_gas_pressure_and_momentum(&gate, 0.05, |_| Ok(0.125 / 25_000.0))
        .unwrap();
    assert!(matches!(
        refused.step_with_gas_pressure_and_momentum(&gate, 0.05, |_| 100_000.125),
        Err(MaterialFreeSurfaceMomentumError2::Wetting(_))
    ));
    same_state(refused.flow().lattice(), &before);
    assert_eq!(refused.samples(), &samples);

    let supported = dynamic_card(&bulk, [angle; 2], [angle; 2], 100.0, false);
    let mut retry = loaded_flow(&bulk, &supported);
    let mut reference = loaded_flow(&bulk, &supported);
    let before = retry.flow().lattice().clone();
    let samples = retry.samples().clone();
    let visits = std::cell::Cell::new(0);
    assert!(matches!(
        retry.step_with_gas_pressure_and_momentum(&gate, 0.05, |_| {
            visits.set(visits.get() + 1);
            if visits.get() == 2 {
                f64::NAN
            } else {
                100_000.125
            }
        }),
        Err(MaterialFreeSurfaceMomentumError2::Pressure(
            FreeSurfacePressureError2::Invalid { .. }
        ))
    ));
    assert_eq!(visits.get(), 2);
    same_state(retry.flow().lattice(), &before);
    assert_eq!(retry.samples(), &samples);
    let cancelled = CancelGate::new_clock_free();
    assert!(matches!(
        retry.step_with_gas_link_pressure_and_momentum(&cancelled, 0.05, |_| {
            cancelled.request();
            100_000.125
        }),
        Err(MaterialFreeSurfaceMomentumError2::Pressure(
            FreeSurfacePressureError2::Cancelled
        ))
    ));
    same_state(retry.flow().lattice(), &before);
    assert_eq!(retry.samples(), &samples);
    let a = retry
        .step_with_gas_pressure_and_momentum(&gate, 0.05, |_| 100_000.125)
        .unwrap();
    let b = reference
        .step_with_gas_link_pressure_and_momentum(&gate, 0.05, |_| 100_000.125)
        .unwrap();
    assert_eq!(a, b);
    assert!(a.native.absolute_angular_momentum_defect.is_some());
    same_state(retry.flow().lattice(), reference.flow().lattice());
    assert_eq!(retry.samples(), reference.samples());
    assert!(retry.samples().values().any(|s| s.signed_speed_m_s != 0.0));
}
