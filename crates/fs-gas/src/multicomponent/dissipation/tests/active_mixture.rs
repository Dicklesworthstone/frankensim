use super::*;
use crate::multicomponent::{NormalPrimitive, SpatialPrimitive};
use fs_thermochem::MassSpecificEnergyQuantityV1;

fn mixture(mass_factors: &[f64]) -> MixtureGas {
    let names = ["A", "B", "C", "D"];
    MixtureGas::new(
        bank(
            mass_factors
                .iter()
                .enumerate()
                .map(|(k, mass)| source(names[k], mass * R, 3.5, 0.0, 0.0))
                .collect(),
        ),
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
        || false,
    )
    .unwrap()
}

fn flow(mass_factors: &[f64]) -> MixtureEuler1d {
    field(mixture(mass_factors), 8, Boundaries::Periodic, |x| {
        let a = 0.5 + 0.1 * (2.0 * core::f64::consts::PI * x).sin();
        let mut densities = vec![0.0; mass_factors.len()];
        densities[0] = a;
        densities[1] = 2.0 * (1.0 - a);
        (densities, 0.0, 2.0)
    })
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * expected.abs().max(1e-15),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn inactive_masses_do_not_exhaust_engineering_or_conjugate_step_budgets() {
    for conjugate in [false, true] {
        let mut reference: Option<(Vec<f64>, f64)> = None;
        for masses in [[1.0, 2.0, 1.0, 1.0], [1.0, 2.0, 1048576.0, 1.0 / 1048576.0]] {
            let mut flow = flow(&masses);
            let initial = flow.states.clone();
            let mut op = operator(&flow, 0.0, 0.0, 0.0, 0.03);
            let controls = DissipationControl {
                maximum_steps: 1,
                ..Default::default()
            };
            // The active A/B ratio is two in both catalogs. The unused C/D
            // masses used to inflate the time-step estimate by about 2^39.
            let report = if conjugate {
                let (report, entropy) = op
                    .advance_entropy_consistent_for(
                        &mut flow,
                        Time::new(0.01),
                        controls,
                        MolecularEntropyControl::default(),
                        || false,
                    )
                    .unwrap();
                assert!(entropy.entropy_change_j_k > 0.0);
                report
            } else {
                op.advance_for(&mut flow, Time::new(0.01), controls, || false)
                    .unwrap()
            };
            assert_eq!(report.accepted_steps, 1);
            assert_eq!(report.rejected_steps, 0);
            assert_ne!(flow.states, initial);
            for i in 0..flow.cell_count() {
                assert_eq!(&flow.cell(i).unwrap().partial_densities[2..], &[0.0, 0.0]);
            }
            if let Some((states, scale)) = &reference {
                for (&actual, &expected) in flow.states.iter().zip(states) {
                    near(actual, expected);
                }
                near(report.maximum_parabolic_number, *scale);
            } else {
                reference = Some((flow.states.clone(), report.maximum_parabolic_number));
            }
        }
    }
}

#[test]
fn inactive_spectators_preserve_planar_and_spatial_fluxes_and_coupled_scales() {
    for thermal in [false, true] {
        let mut reference: Option<([f64; 5], f64)> = None;
        for masses in [
            &[1.0, 2.0][..],
            &[1.0, 2.0, 1.0, 1.0],
            &[1.0, 2.0, 1048576.0, 1.0 / 1048576.0],
        ] {
            let flow = flow(masses);
            let gas = flow.gas();
            let n = gas.species_count();
            let molecular = law(&flow, 0.0, 0.0, 0.01, 0.03);
            let authored = [0.5, -0.25, 7.0, -9.0];
            let heats = gas
                .source()
                .models()
                .iter()
                .enumerate()
                .map(|(k, model)| {
                    (
                        model.species().clone(),
                        MassSpecificEnergyQuantityV1::new(if thermal { authored[k] } else { 0.0 }),
                    )
                })
                .collect();
            let transport = SoretDufourTransport::new(
                molecular.clone(),
                "manufactured-active-mixture",
                1,
                molecular.domain(),
                heats,
            )
            .unwrap();
            let mut planar =
                PlanarTransport::new(gas.clone(), transport.clone(), || false).unwrap();
            let mut spatial = SpatialTransport::new(gas.clone(), transport, || false).unwrap();
            let left = [0.8, 0.4, 0.0, 0.0];
            let right = [0.2, 1.6, 0.0, 0.0];
            let mut flux = vec![0.0; n + 3];
            let report = planar
                .face_from_primitives(
                    NormalPrimitive {
                        partial_densities: &left[..n],
                        normal_velocity: 0.0,
                        transverse_velocity: 0.0,
                        pressure: 1.5,
                    },
                    NormalPrimitive {
                        partial_densities: &right[..n],
                        normal_velocity: 0.0,
                        transverse_velocity: 0.0,
                        pressure: 2.5,
                    },
                    0.5,
                    [0.0; 2],
                    MaxwellStefanSolveV1::default(),
                    &mut flux,
                    || false,
                )
                .unwrap();
            let mut flux3 = vec![0.0; n + 4];
            let report3 = spatial
                .face_from_primitives(
                    SpatialPrimitive {
                        partial_densities: &left[..n],
                        normal_velocity: 0.0,
                        transverse_velocity: [0.0; 2],
                        pressure: 1.5,
                    },
                    SpatialPrimitive {
                        partial_densities: &right[..n],
                        normal_velocity: 0.0,
                        transverse_velocity: [0.0; 2],
                        pressure: 2.5,
                    },
                    0.5,
                    [[0.0; 3]; 2],
                    MaxwellStefanSolveV1::default(),
                    &mut flux3,
                    || false,
                )
                .unwrap();
            assert_ne!(flux[0], 0.0);
            assert!(flux[2..n].iter().all(|&j| j == 0.0));
            assert_eq!(flux3[..n + 2], flux[..n + 2]);
            assert_eq!(flux3[n + 2], 0.0);
            assert_eq!(flux3[n + 3], flux[n + 2]);
            assert_eq!(report3.diffusivity_m2_s, report.diffusivity_m2_s);
            let active_flux = [flux[0], flux[1], flux[n], flux[n + 1], flux[n + 2]];
            if let Some((expected_flux, scale)) = &reference {
                for (actual, &expected) in active_flux.into_iter().zip(expected_flux) {
                    near(actual, expected);
                }
                near(report.diffusivity_m2_s, *scale);
            } else {
                reference = Some((active_flux, report.diffusivity_m2_s));
            }
        }
    }
}

#[test]
fn one_sided_and_trace_species_restore_their_full_molecular_mass_scale() {
    let flow = flow(&[1.0, 2.0, 64.0, 1048576.0]);
    let mut planar =
        PlanarTransport::new(flow.gas().clone(), law(&flow, 0.0, 0.0, 0.0, 0.03), || {
            false
        })
        .unwrap();
    for heavy_fraction in [0.25, 1e-200] {
        let a = [0.8, 0.4, 0.0, 0.0];
        let b = [0.5 - heavy_fraction, 1.0, 64.0 * heavy_fraction, 0.0];
        for [left, right] in [[a, b], [b, a]] {
            let mut flux = [0.0; 7];
            let report = planar
                .face_from_primitives(
                    NormalPrimitive {
                        partial_densities: &left,
                        normal_velocity: 0.0,
                        transverse_velocity: 0.0,
                        pressure: 2.0,
                    },
                    NormalPrimitive {
                        partial_densities: &right,
                        normal_velocity: 0.0,
                        transverse_velocity: 0.0,
                        pressure: 2.0,
                    },
                    0.5,
                    [0.0; 2],
                    MaxwellStefanSolveV1::default(),
                    &mut flux,
                    || false,
                )
                .unwrap();
            // All species have Cp/R=3.5, so gamma=1.4. C is present on one
            // trace and must retain M_C/M_A=64 even at x_C=1e-200. Only D
            // is irrelevant to the bound; no trace floor may erase C's flux.
            near(report.diffusivity_m2_s, 0.03 * 64.0 * 1.4);
            assert_ne!(flux[2], 0.0);
            assert_eq!(flux[3], 0.0);
            assert_eq!(flux[2].is_sign_positive(), right[2] < left[2]);
        }
    }
}
