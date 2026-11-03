//! Actual P1 liquid strain enters native stress without extra impulse or heat.

use super::*;

#[test]
fn g1_p1_translation_rotation_and_shear_retain_their_physical_strain() {
    let material = liquid(800.0);
    let (points, tets, spin) = cube([0.0; 3]);
    let translation = vec![[0.01, -0.005, 0.002]; points.len()];
    let shear: Vec<_> = points
        .iter()
        .map(|p| {
            let r = sub(*p, [0.004; 3]);
            [
                0.01 + r[0] + 2.0 * r[1] - r[2],
                -0.005 - 0.5 * r[1] + 1.5 * r[2],
                0.002 + 0.5 * r[0] - 0.5 * r[1] - 0.5 * r[2],
            ]
        })
        .collect();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        for gravity in [[0.0; 3], [0.3, -0.5, 0.7]] {
            let mut receiving = domain();
            receiving.gravity_m_s2 = gravity;
            for (velocity, expected_rate) in
                [(&translation, 0.0), (&spin, 0.0), (&shear, 8.25_f64.sqrt())]
            {
                let scaling =
                    plan_free_surface_from_material(&material, 0.001, 0.0005, 0.1).unwrap();
                let (mut flow, report) = scaling
                    .free_surface3_from_tetrahedra(
                        receiving,
                        [0.0; 3],
                        TetrahedralLiquidState {
                            positions_m: &points,
                            tetrahedra: &tets,
                            velocities_m_s: velocity,
                        },
                        limits(),
                        cx,
                    )
                    .unwrap();
                let mut checked = 0;
                for (location, mass) in flow.lattice.mass_snapshot() {
                    if mass <= 0.0 {
                        continue;
                    }
                    let rate = flow
                        .lattice
                        .transport_shear_rate(location, scaling.flow.tau)
                        .unwrap()
                        / scaling.timestep_s();
                    close(rate, expected_rate, 5e-10);
                    checked += 1;
                }
                assert_eq!(checked, report.wet_cells);
                close(report.target.mass_kg, report.source.mass_kg, 1e-17);
                for axis in 0..3 {
                    close(
                        report.target.momentum_kg_m_s[axis],
                        report.source.momentum_kg_m_s[axis],
                        1e-16,
                    );
                    close(
                        report.target.angular_momentum_kg_m2_s[axis],
                        report.source.angular_momentum_kg_m2_s[axis],
                        1e-18,
                    );
                }
                let mass = flow.ledger_mass_kg().unwrap();
                flow.step().unwrap();
                close(flow.ledger_mass_kg().unwrap(), mass, 1e-14);
            }
        }
    });
}

#[test]
fn g4_unrepresentable_viscous_ingress_returns_no_flow_and_retains_source_fields() {
    let material = liquid(800.0);
    let (points, tets, _) = cube([0.0; 3]);
    let velocity: Vec<_> = points
        .iter()
        .map(|p| [0.12 * (p[1] - 0.004), 0.0, 0.0])
        .collect();
    let original = (points.clone(), tets.clone(), velocity.clone());
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        // The original low-Mach equilibrium initializer can represent this
        // velocity. Its physical viscosity/strain at the deliberately coarse
        // timestep cannot be represented by positive D3Q19 populations.
        let scaling = plan_free_surface_from_material(&material, 0.001, 0.1, 0.0003).unwrap();
        let mut cells = Vec::new();
        for z in 2..6 {
            for y in 2..6 {
                for x in 2..6 {
                    cells.push(MaterialFreeSurfaceInitialCell3 {
                        coordinates: [x, y, z],
                        fill_fraction: 1.0,
                        velocity_m_s: [0.12 * ((y as f64 + 0.5) * 0.001 - 0.004), 0.0, 0.0],
                    });
                }
            }
        }
        assert!(
            scaling
                .free_surface3_from_cells(domain(), &cells, cx)
                .is_ok()
        );
        let result = scaling.free_surface3_from_tetrahedra(
            domain(),
            [0.0; 3],
            TetrahedralLiquidState {
                positions_m: &points,
                tetrahedra: &tets,
                velocities_m_s: &velocity,
            },
            limits(),
            cx,
        );
        assert!(matches!(result, Err(TetLiquidTransferError::Flow(_))));
        let error = result.err().unwrap().to_string();
        assert!(
            error.contains("initial viscous stress needs positive"),
            "{error}"
        );
        assert_eq!((points.clone(), tets.clone(), velocity.clone()), original);
    });
}
