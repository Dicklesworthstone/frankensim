//! G1/G3: total physical wrench on the material owner's actual partial mass.

use super::*;
use fs_lbm::material::free_surface::{
    MaterialFreeSurfaceWrenchImpulse3, MaterialFreeSurfaceWrenchImpulseControls3,
};

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[test]
fn total_si_wrench_uses_actual_partial_liquid_inertia_at_two_timesteps() {
    let source = liquid();
    let indices = [0, 1, 12, 49];
    let positions = indices.map(|i| [i % 4, i / 4 % 4, i / 16]);
    let centers = positions.map(|p| p.map(|v| (v as f64 + 0.5) * 0.001));
    for dt in [0.001, 0.0005] {
        let mut flow = flow3(&source, dt);
        let before = positions.map(|[x, y, z]| flow.sample_at_si(x, y, z).unwrap().unwrap());
        let masses = before.map(|s| s.density_kg_m3 * s.fill_fraction * 1e-9);
        let total = masses.into_iter().sum::<f64>();
        let center: [f64; 3] = core::array::from_fn(|a| {
            (0..4).map(|i| masses[i] * centers[i][a]).sum::<f64>() / total
        });
        let jumps: [[f64; 3]; 4] = core::array::from_fn(|i| {
            let arm = core::array::from_fn(|a| centers[i][a] - center[a]);
            let rotation = cross([2.0, -1.0, 3.0], arm);
            core::array::from_fn(|a| [0.02, -0.01, 0.015][a] + rotation[a])
        });
        let mut wrench = MaterialFreeSurfaceWrenchImpulse3 {
            origin_m: [0.00025, -0.002, 0.001],
            linear_impulse_kg_m_s: [0.0; 3],
            angular_impulse_kg_m2_s: [0.0; 3],
        };
        for i in 0..4 {
            let impulse = jumps[i].map(|u| u * masses[i]);
            let arm = core::array::from_fn(|a| centers[i][a] - wrench.origin_m[a]);
            let moment = cross(arm, impulse);
            for a in 0..3 {
                wrench.linear_impulse_kg_m_s[a] += impulse[a];
                wrench.angular_impulse_kg_m2_s[a] += moment[a];
            }
        }
        let mass = flow.ledger_mass_kg().unwrap().to_bits();
        let report = flow
            .apply_wrench_impulse(
                &indices,
                wrench,
                MaterialFreeSurfaceWrenchImpulseControls3 {
                    impulse: controls3(),
                    maximum_angular_momentum_defect_kg_m2_s: 1e-19,
                },
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        let after = positions.map(|[x, y, z]| flow.sample_at_si(x, y, z).unwrap().unwrap());
        let mut kinetic = 0.0;
        for i in 0..4 {
            assert_eq!(
                before[i].fill_fraction.to_bits(),
                after[i].fill_fraction.to_bits()
            );
            assert_eq!(
                before[i].density_kg_m3.to_bits(),
                after[i].density_kg_m3.to_bits()
            );
            for a in 0..3 {
                near(
                    after[i].velocity_m_s[a] - before[i].velocity_m_s[a],
                    jumps[i][a],
                    2e-14,
                );
                kinetic += 0.5
                    * masses[i]
                    * (after[i].velocity_m_s[a].powi(2) - before[i].velocity_m_s[a].powi(2));
            }
        }
        near(report.stencil_mass_kg, total, 1e-20);
        near(report.impulse.work_j, kinetic, 1e-21);
        for a in 0..3 {
            near(report.stencil_center_of_mass_m[a], center[a], 1e-17);
            near(
                report.impulse.momentum_change_kg_m_s[a],
                wrench.linear_impulse_kg_m_s[a],
                1e-20,
            );
            near(
                report.angular_momentum_change_kg_m2_s[a],
                wrench.angular_impulse_kg_m2_s[a],
                1e-22,
            );
        }
        assert_eq!(flow.ledger_mass_kg().unwrap().to_bits(), mass);
        assert_eq!(flow.lattice().grid().steps(), 0);
        assert!(std::ptr::eq(flow.scaling().material(), &source));
    }
}
