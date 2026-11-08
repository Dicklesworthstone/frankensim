//! Physical planar torque uses original SI units per unit depth.

use super::*;
use fs_lbm::material::free_surface::{
    MaterialFreeSurfaceWrenchImpulse2, MaterialFreeSurfaceWrenchImpulseControls2,
};

const INDICES: [usize; 4] = [7, 8, 25, 26];

fn cross(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}

fn controls() -> MaterialFreeSurfaceWrenchImpulseControls2 {
    MaterialFreeSurfaceWrenchImpulseControls2 {
        impulse: controls2(),
        maximum_angular_momentum_defect_kg_m_s: 1e-15,
    }
}

fn scaled_flow(source: &ResolvedLiquidSurfaceState, dx: f64, dt: f64) -> MaterialFreeSurface2<'_> {
    plan_free_surface_from_material(source, dx, dt, 0.1)
        .unwrap()
        .free_surface2(
            MaterialFreeSurfaceDomain2 {
                cells: [6; 2],
                gravity_m_s2: [0.1, -0.2],
                contact: ContactModel::Neutral,
                max_cells: 36,
            },
            |x, _| x < 3,
        )
        .unwrap()
}

#[test]
fn planar_si_wrench_uses_partial_inertia_and_distinct_length_time_scales() {
    let source = liquid();
    let gate = CancelGate::new_clock_free();
    let coordinates = INDICES.map(|i| [i % 6, i / 6]);
    for (dx, dt) in [(0.001, 0.001), (0.001, 0.0005), (0.0005, 0.00025)] {
        let mut flow = scaled_flow(&source, dx, dt);
        let centers = coordinates.map(|p| p.map(|v| (v as f64 + 0.5) * dx));
        let before = coordinates.map(|[x, y]| flow.sample_at_si(x, y).unwrap().unwrap());
        assert_eq!(before.map(|s| s.fill_fraction), [1.0, 0.5, 1.0, 0.5]);
        let masses = before.map(|s| s.density_kg_m3 * s.fill_fraction * dx * dx);
        let total = masses.into_iter().sum::<f64>();
        let center: [f64; 2] = core::array::from_fn(|a| {
            (0..4).map(|i| masses[i] * centers[i][a]).sum::<f64>() / total
        });
        let jumps: [[f64; 2]; 4] = core::array::from_fn(|i| {
            [
                0.02 - 3.0 * (centers[i][1] - center[1]),
                -0.01 + 3.0 * (centers[i][0] - center[0]),
            ]
        });
        let mut wrench = MaterialFreeSurfaceWrenchImpulse2 {
            origin_m: [0.00025, -0.002],
            linear_impulse_kg_s: [0.0; 2],
            angular_impulse_kg_m_s: 0.0,
        };
        for i in 0..4 {
            let impulse = jumps[i].map(|u| u * masses[i]);
            let arm = core::array::from_fn(|a| centers[i][a] - wrench.origin_m[a]);
            wrench.angular_impulse_kg_m_s += cross(arm, impulse);
            for a in 0..2 {
                wrench.linear_impulse_kg_s[a] += impulse[a];
            }
        }
        let old_populations = flow.lattice().grid.f.clone();
        let old_mass = flow.lattice().mass.clone();
        let old_flags = flow.lattice().grid.flags.clone();
        let old_force = (flow.lattice().grid.g, flow.lattice().grid.fext.clone());
        let ledger = flow.ledger_mass_kg_per_m().unwrap().to_bits();
        let report = flow
            .apply_wrench_impulse(&INDICES, wrench, controls(), &gate)
            .unwrap();
        let after = coordinates.map(|[x, y]| flow.sample_at_si(x, y).unwrap().unwrap());
        let mut observed = [0.0; 2];
        let mut angular = 0.0;
        let mut kinetic = 0.0;
        for i in 0..4 {
            assert_eq!(before[i].cell, after[i].cell);
            assert_eq!(
                before[i].fill_fraction.to_bits(),
                after[i].fill_fraction.to_bits()
            );
            assert_eq!(
                before[i].density_kg_m3.to_bits(),
                after[i].density_kg_m3.to_bits()
            );
            let impulse: [f64; 2] = core::array::from_fn(|a| {
                masses[i] * (after[i].velocity_m_s[a] - before[i].velocity_m_s[a])
            });
            for a in 0..2 {
                near(
                    after[i].velocity_m_s[a] - before[i].velocity_m_s[a],
                    jumps[i][a],
                    2e-14,
                );
                observed[a] += impulse[a];
                kinetic += 0.5
                    * masses[i]
                    * (after[i].velocity_m_s[a].powi(2) - before[i].velocity_m_s[a].powi(2));
            }
            angular += cross(
                core::array::from_fn(|a| centers[i][a] - wrench.origin_m[a]),
                impulse,
            );
        }
        near(report.stencil_mass_kg_m, total, 2e-17);
        for a in 0..2 {
            near(report.stencil_center_of_mass_m[a], center[a], 2e-17);
            near(observed[a], wrench.linear_impulse_kg_s[a], 5e-17);
            near(report.impulse.momentum_change_kg_s[a], observed[a], 5e-17);
        }
        near(angular, wrench.angular_impulse_kg_m_s, 5e-20);
        near(report.angular_momentum_change_kg_m_s, angular, 5e-20);
        near(report.impulse.work_j_m, kinetic, 2e-18);
        near(report.impulse.kinetic_energy_change_j_m, kinetic, 2e-18);
        assert_eq!(report.impulse.mass_change_kg_m, 0.0);
        assert_eq!(report.impulse.affected_cells, 4);
        for (i, original) in old_populations.iter().enumerate() {
            if !INDICES.contains(&i) {
                assert_eq!(
                    original.map(f64::to_bits),
                    flow.lattice().grid.f[i].map(f64::to_bits)
                );
            }
        }
        assert_eq!(flow.lattice().mass, old_mass);
        assert_eq!(flow.lattice().grid.flags, old_flags);
        assert_eq!(flow.lattice().grid.g, old_force.0);
        assert_eq!(flow.lattice().grid.fext, old_force.1);
        assert_eq!(flow.ledger_mass_kg_per_m().unwrap().to_bits(), ledger);
        assert!(std::ptr::eq(flow.scaling().material(), &source));
    }
}

#[test]
fn planar_si_origin_covariance_and_reverse_wrench_preserve_signed_work() {
    let source = liquid();
    let gate = CancelGate::new_clock_free();
    let mut flow = flow2(&source, 0.0005);
    let mut shifted = flow2(&source, 0.0005);
    let original = flow.lattice().grid.f.clone();
    let wrench = MaterialFreeSurfaceWrenchImpulse2 {
        origin_m: [0.00025, -0.002],
        linear_impulse_kg_s: [0.00003, -0.000015],
        angular_impulse_kg_m_s: 0.00000004,
    };
    let shifted_origin = [0.010123456789, -0.020987654321];
    let shifted_wrench = MaterialFreeSurfaceWrenchImpulse2 {
        origin_m: shifted_origin,
        angular_impulse_kg_m_s: wrench.angular_impulse_kg_m_s
            + cross(
                core::array::from_fn(|a| wrench.origin_m[a] - shifted_origin[a]),
                wrench.linear_impulse_kg_s,
            ),
        ..wrench
    };
    let report = flow
        .apply_wrench_impulse(&INDICES, wrench, controls(), &gate)
        .unwrap();
    let shifted_report = shifted
        .apply_wrench_impulse(&INDICES, shifted_wrench, controls(), &gate)
        .unwrap();
    near(
        shifted_report.angular_momentum_change_kg_m_s,
        shifted_wrench.angular_impulse_kg_m_s,
        5e-20,
    );
    for i in 0..36 {
        for q in 0..fs_lbm::Q {
            near(
                flow.lattice().grid.f[i][q],
                shifted.lattice().grid.f[i][q],
                2e-16,
            );
        }
    }
    let reverse = MaterialFreeSurfaceWrenchImpulse2 {
        linear_impulse_kg_s: wrench.linear_impulse_kg_s.map(|value| -value),
        angular_impulse_kg_m_s: -wrench.angular_impulse_kg_m_s,
        ..wrench
    };
    let braking = flow
        .apply_wrench_impulse(&INDICES, reverse, controls(), &gate)
        .unwrap();
    assert!(braking.impulse.work_j_m < 0.0);
    near(
        report.impulse.work_j_m + braking.impulse.work_j_m,
        0.0,
        2e-18,
    );
    near(
        braking.impulse.work_j_m,
        braking.impulse.kinetic_energy_change_j_m,
        2e-18,
    );
    for (i, before) in original.iter().enumerate() {
        for q in 0..fs_lbm::Q {
            near(flow.lattice().grid.f[i][q], before[q], 3e-16);
        }
    }
}

#[test]
fn planar_si_native_and_physical_refusals_leave_the_original_owner_unchanged() {
    let source = liquid();
    let gate = CancelGate::new_clock_free();
    let mut flow = flow2(&source, 0.001);
    let before = flow.lattice().grid.f.clone();
    let mass = flow.lattice().mass.clone();
    let flags = flow.lattice().grid.flags.clone();
    let ledger = flow.ledger_mass_kg_per_m().unwrap().to_bits();
    let zero = MaterialFreeSurfaceWrenchImpulse2 {
        origin_m: [0.0; 2],
        linear_impulse_kg_s: [0.0; 2],
        angular_impulse_kg_m_s: 0.0,
    };
    let mut wide = controls();
    wide.impulse.maximum_speed_m_s = 0.11;
    assert!(matches!(
        flow.apply_wrench_impulse(&INDICES, zero, wide, &gate),
        Err(MaterialFreeSurfaceImpulseError::Units(_))
    ));
    let mut invalid_budget = controls();
    invalid_budget.maximum_angular_momentum_defect_kg_m_s = -1.0;
    assert!(matches!(
        flow.apply_wrench_impulse(&INDICES, zero, invalid_budget, &gate),
        Err(MaterialFreeSurfaceImpulseError::Units(_))
    ));
    for bad in [
        MaterialFreeSurfaceWrenchImpulse2 {
            origin_m: [f64::MAX, 0.0],
            ..zero
        },
        MaterialFreeSurfaceWrenchImpulse2 {
            angular_impulse_kg_m_s: f64::NAN,
            ..zero
        },
    ] {
        assert!(matches!(
            flow.apply_wrench_impulse(&INDICES, bad, controls(), &gate),
            Err(MaterialFreeSurfaceImpulseError::Units(_))
        ));
    }
    assert!(matches!(
        flow.apply_wrench_impulse(&[7, 35], zero, controls(), &gate),
        Err(MaterialFreeSurfaceImpulseError::Native(_))
    ));
    assert!(matches!(
        flow.apply_wrench_impulse(
            &[7],
            MaterialFreeSurfaceWrenchImpulse2 {
                angular_impulse_kg_m_s: 1e-8,
                ..zero
            },
            controls(),
            &gate,
        ),
        Err(MaterialFreeSurfaceImpulseError::Native(_))
    ));
    let identity = flow
        .apply_wrench_impulse(&INDICES, zero, controls(), &gate)
        .unwrap();
    assert_eq!(identity.impulse.affected_cells, 0);
    assert_eq!(identity.impulse.work_j_m, 0.0);
    gate.request();
    assert_eq!(
        flow.apply_wrench_impulse(&INDICES, zero, controls(), &gate),
        Err(MaterialFreeSurfaceImpulseError::Native(
            ImpulseError::Cancelled
        ))
    );
    assert_eq!(flow.lattice().grid.f, before);
    assert_eq!(flow.lattice().mass, mass);
    assert_eq!(flow.lattice().grid.flags, flags);
    assert_eq!(flow.ledger_mass_kg_per_m().unwrap().to_bits(), ledger);
}
