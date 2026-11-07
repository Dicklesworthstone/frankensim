//! G1/G3: resolved held-bond heat below the spacing of stored-energy totals.
//! References integrate authored Maxwell spectra analytically; only the public
//! cohesive API supplies actual states, tractions and energy increments.

use fs_material::cohesive::relaxation::{IntactViscoelasticBond, ViscoelasticBondLimits};
use fs_material::visco::GeneralizedMaxwell;

const THICKNESS_M: f64 = 0.001;
const PRELOAD_S: f64 = 1e-4;
const GAP_M: [f64; 3] = [2e-5, -1e-5, 0.5e-5];
const NORMAL: [(f64, f64); 2] = [(2e6, 0.1), (3e5, 0.8)];
const SHEAR: [(f64, f64); 1] = [(8e5, 0.25)];

fn spectrum(axis: usize) -> (f64, &'static [(f64, f64)]) {
    if axis == 0 {
        (1e6, &NORMAL)
    } else {
        (4e5, &SHEAR)
    }
}

fn law() -> IntactViscoelasticBond {
    IntactViscoelasticBond::new(
        GeneralizedMaxwell::new(1e6, NORMAL.to_vec()).unwrap(),
        GeneralizedMaxwell::new(4e5, SHEAR.to_vec()).unwrap(),
        THICKNESS_M,
        ViscoelasticBondLimits {
            maximum_jump_m: 0.0001,
            maximum_rate_m_s: 1.0,
            maximum_elapsed_s: 10.0,
        },
    )
    .unwrap()
}

// Solve dq/dt + q/tau = strain_rate for a stress-free linear preload.
// This reference uses platform expm1, not the production PronyRamp owner.
fn preload_branch_strain(strain: f64, tau_s: f64) -> f64 {
    strain * (tau_s / PRELOAD_S) * -(-PRELOAD_S / tau_s).exp_m1()
}

#[test]
fn weak_held_bond_relaxation_matches_analytic_spring_energy_release() {
    let law = law();
    let loaded = law
        .bonded_point()
        .trial(&law.prepare(PRELOAD_S).unwrap(), GAP_M, 1e-10)
        .unwrap()
        .point;
    for (axis, jump) in GAP_M.into_iter().enumerate() {
        let (equilibrium, terms) = spectrum(axis);
        let strain = jump / THICKNESS_M;
        let expected = equilibrium * strain
            + terms
                .iter()
                .map(|&(modulus, tau)| modulus * preload_branch_strain(strain, tau))
                .sum::<f64>();
        assert!((loaded.traction_pa()[axis] - expected).abs() <= 1e-13 * expected.abs());
    }

    for rate in [1e-17, 1e-14, 1e-10, 1e-6, 0.1, 1.0, 10.0] {
        let duration = rate * NORMAL[0].1;
        // Held strain leaves the equilibrium spring unchanged. Each Maxwell
        // spring releases h*E*q0^2*(1-exp(-2*dt/tau))/2 as physical heat.
        let expected_heat = THICKNESS_M
            * GAP_M
                .into_iter()
                .enumerate()
                .map(|(axis, jump)| {
                    spectrum(axis)
                        .1
                        .iter()
                        .map(|&(modulus, tau)| {
                            let q0 = preload_branch_strain(jump / THICKNESS_M, tau);
                            0.5 * modulus * q0 * q0 * -(-2.0 * duration / tau).exp_m1()
                        })
                        .sum::<f64>()
                })
                .sum::<f64>();
        let tolerance = 1e-13 * expected_heat;
        let step = loaded
            .trial(&law.prepare(duration).unwrap(), loaded.jump_m(), tolerance)
            .expect("resolved physical heat must survive rounded endpoint storage");
        assert_eq!(step.work_j_m2.to_bits(), 0);
        assert!(step.stored_change_j_m2 < 0.0 && step.viscous_dissipated_j_m2 > 0.0);
        assert!((step.stored_change_j_m2 + expected_heat).abs() <= tolerance);
        assert!((step.viscous_dissipated_j_m2 - expected_heat).abs() <= tolerance);
        assert!(step.defect_j_m2.abs() <= tolerance);
        if rate == 1e-17 {
            assert_eq!(step.traction_pa, loaded.traction_pa());
            assert_eq!(step.point.stored_energy_j_m2(), loaded.stored_energy_j_m2());
        }
    }
}

#[test]
fn unstressed_bond_retains_exact_zero_work_storage_and_heat() {
    let law = law();
    let mut point = law.bonded_point();
    for duration in [1e-18, 1e-4, 0.1, 1.0] {
        let step = point
            .trial(&law.prepare(duration).unwrap(), [0.0; 3], 1e-30)
            .unwrap();
        assert_eq!(step.traction_pa, [0.0; 3]);
        assert_eq!(step.mean_traction_pa, [0.0; 3]);
        for value in [
            step.work_j_m2,
            step.stored_change_j_m2,
            step.viscous_dissipated_j_m2,
            step.defect_j_m2,
            step.point.stored_energy_j_m2(),
        ] {
            assert_eq!(value.to_bits(), 0);
        }
        point = step.point;
    }
}
