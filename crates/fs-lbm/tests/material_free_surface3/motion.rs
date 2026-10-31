//! Sourced units and the real moving-wall VOF trajectory; synthetic liquid data.

use super::*;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

fn sourced_viscosity(mu: f64) -> ResolvedLiquidSurfaceState {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(claim(DENSITY_PROPERTY, 1000.0, Density::DIMS))
        .unwrap();
    claims
        .insert_claim(claim(
            DYNAMIC_VISCOSITY_PROPERTY,
            mu,
            DYNAMIC_VISCOSITY_DIMS,
        ))
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic moving-wall liquid".into(),
            phase: "liquid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_liquid_surface_state(
        &card,
        &point(300.0),
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

#[test]
fn g1_wall_motion_si_matches_native_impulse_work_mass_and_explicit_pressure_datum() {
    let source = liquid(0, 300.0, None);
    let setup = plan(&source).with_pressure_reference(100_000.0).unwrap();
    let mut fluid = setup.free_surface3(domain(), |_, _, z| z < 6).unwrap();
    let mut native = fluid.lattice().clone();
    let unit = setup.velocity_si(1.0).unwrap();
    let velocity = [0.03, 0.0, 0.0];
    let raw = native
        .step_with_planar_wall_motion(
            PlanarWall3::ZMin,
            velocity.map(|v| v / unit),
            1e-12,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    let report = fluid
        .step_with_planar_wall_motion_si(
            PlanarWall3::ZMin,
            velocity,
            1e-12,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(
        fluid.lattice().grid().state_bits(),
        native.grid().state_bits()
    );
    assert_eq!(fluid.scaling().material(), &source);
    let plane_cells = 256.0;
    let cell_mass = 1000.0 * setup.cell_size_m().powi(3);
    let expected_impulse = cell_mass * plane_cells * velocity[0] / 3.0;
    close(report.fluid_impulse_kg_m_s[0], expected_impulse, 2e-18);
    close(report.fluid_work_j, velocity[0] * expected_impulse, 2e-19);
    assert_eq!(report.wall_work_j, -report.fluid_work_j);
    for axis in 0..3 {
        assert_eq!(
            report.wall_impulse_kg_m_s[axis],
            -report.fluid_impulse_kg_m_s[axis]
        );
        close(
            report.fluid_impulse_kg_m_s[axis],
            raw.fluid_impulse[axis] * setup.mass_unit_kg() * unit,
            2e-17,
        );
    }
    close(report.pressure_perturbation_wall_force_n[2], 0.0, 2e-12);
    let area = plane_cells * setup.cell_size_m().powi(2);
    close(report.wetted_area_m2, area, 1e-18);
    close(
        report.absolute_wall_force_n.unwrap()[2],
        -100_000.0 * area,
        2e-12,
    );
    close(report.mass_before_kg, report.mass_after_kg, 2e-16);

    let mut shifted = plan(&source)
        .with_pressure_reference(101_000.0)
        .unwrap()
        .free_surface3(domain(), |_, _, z| z < 6)
        .unwrap();
    let other = shifted
        .step_with_planar_wall_motion_si(
            PlanarWall3::ZMin,
            velocity,
            1e-12,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(
        shifted.lattice().grid().state_bits(),
        fluid.lattice().grid().state_bits()
    );
    assert_eq!(other.fluid_work_j, report.fluid_work_j);
    close(
        other.absolute_wall_force_n.unwrap()[2] - report.absolute_wall_force_n.unwrap()[2],
        -1000.0 * area,
        2e-12,
    );
}

#[test]
fn g3_original_viscosity_changes_wall_driven_liquid_and_reaction_over_native_steps() {
    let sources = [sourced_viscosity(0.1), sourced_viscosity(0.2)];
    let mut flows: Vec<_> = sources
        .iter()
        .map(|source| {
            plan(source)
                .free_surface3(domain(), |_, _, z| z < 6)
                .unwrap()
        })
        .collect();
    let initial = flows[0].ledger_mass_kg().unwrap();
    let mut reactions = [0.0; 2];
    for _ in 0..4 {
        for (i, flow) in flows.iter_mut().enumerate() {
            let report = flow
                .step_with_planar_wall_motion_si(
                    PlanarWall3::ZMin,
                    [0.03, 0.0, 0.0],
                    1e-12,
                    &CancelGate::new_clock_free(),
                )
                .unwrap();
            assert!(report.absolute_wall_force_n.is_none());
            reactions[i] += report.wall_impulse_kg_m_s[0];
            close(report.mass_after_kg, initial, 1e-14);
        }
    }
    assert_ne!(
        flows[0].lattice().grid().state_bits(),
        flows[1].lattice().grid().state_bits()
    );
    assert!((reactions[0] - reactions[1]).abs() > 1e-10);
    for (flow, source) in flows.iter().zip(&sources) {
        assert_eq!(flow.scaling().material(), source);
        assert!(flow.sample_at_si(8, 8, 0).unwrap().unwrap().velocity_m_s[0] > 0.0);
    }
}

#[test]
fn g4_rebind_and_late_scenario_speed_refusals_preserve_accepted_material_flow() {
    let source = liquid(0, 300.0, None);
    let mut fluid = plan(&source)
        .free_surface3(domain(), |_, _, z| z < 6)
        .unwrap();
    let before = fluid.lattice().grid().state_bits();
    let mass = fluid.ledger_mass_kg().unwrap();
    assert!(
        fluid
            .step_with_planar_wall_motion_si(
                PlanarWall3::ZMin,
                [0.11, 0.0, 0.0],
                1e-12,
                &CancelGate::new_clock_free()
            )
            .is_err()
    );
    let gate = CancelGate::new_clock_free();
    gate.request();
    assert_eq!(
        fluid.step_with_planar_wall_motion_si(PlanarWall3::ZMin, [0.03, 0.0, 0.0], 1e-12, &gate),
        Err(MaterialLbmError::FreeSurface(FreeSurfaceError3::Cancelled))
    );
    assert_eq!(fluid.lattice().grid().state_bits(), before);
    assert_eq!(fluid.ledger_mass_kg().unwrap(), mass);

    // Initial Guo speed is 0.01 m/s; the first actual body-force step gives
    // interior speed about 0.03 m/s. Native Mach remains admissible, but the
    // original 0.02 m/s material scenario bound must reject AFTER the solve.
    let setup = plan_free_surface_from_material(&source, 1e-3, 2e-4, 0.02).unwrap();
    let mut fast = setup
        .free_surface3(
            MaterialFreeSurfaceDomain3 {
                gravity_m_s2: [0.0, 0.0, -100.0],
                ..domain()
            },
            |_, _, z| z < 6,
        )
        .unwrap();
    let before = fast.lattice().grid().state_bits();
    let mass = fast.ledger_mass_kg().unwrap();
    let counters = fast.lattice().conversions();
    assert!(
        matches!(fast.step_with_planar_wall_motion_si(PlanarWall3::ZMin,
        [0.005,0.0,0.0],1e-12,&CancelGate::new_clock_free()),
        Err(MaterialLbmError::Invalid(reason)) if reason.contains("liquid exceeds"))
    );
    assert_eq!(fast.lattice().grid().state_bits(), before);
    assert_eq!(fast.ledger_mass_kg().unwrap(), mass);
    assert_eq!(fast.lattice().conversions(), counters);
}
