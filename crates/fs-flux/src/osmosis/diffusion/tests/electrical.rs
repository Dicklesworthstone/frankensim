//! Independent salt transport and complete charged-network evolution.
use super::*;
use fs_material::osmosis::diffusion::electrical::ZeroCurrentDiffusionReport;
use fs_thermochem::{UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R, stoichiometry::FARADAY_C_PER_MOL};

fn ionic() -> ZeroCurrentDiffusion {
    ZeroCurrentDiffusion::new(vec![1, -1], 1e-15, 1e-16).unwrap()
}
fn salt() -> Arc<VirialSolutionModel> {
    source(1.0, vec![0.0; 4], vec![200.0; 2])
}
fn sample(
    model: &Arc<VirialSolutionModel>,
    left: &[f64],
    right: &[f64],
) -> ([f64; 2], ZeroCurrentDiffusionReport) {
    let law = law(model);
    let mut workspace = law.workspace().unwrap();
    let mut flux = [f64::NAN; 2];
    let report = law
        .flux_zero_current_into(
            cell(model, left).solution(),
            cell(model, right).solution(),
            0.001,
            physical().constitutive,
            &ionic(),
            &mut workspace,
            &mut flux,
            || false,
        )
        .unwrap();
    (flux, report)
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}

#[test]
fn binary_salt_matches_ambipolar_flux_and_floating_diffusion_potential() {
    for reference in [1e-100, 1.0, 1e100] {
        let model = source(reference, vec![0.0; 4], vec![200.0; 2]);
        for (a, b) in [(100.0_f64, 20.0_f64), (20.0, 100.0), (1e-10, 1e-11)] {
            let (flux, report) = sample(&model, &[a, a], &[b, b]);
            // Independent zero-current elimination for a 1:1 salt: D_salt is
            // the harmonic mean, and finite stationary solvent contributes 1/x0.
            let mean = (a - b) / (a / b).ln();
            let expected = 1.5e-9 * (a - b) * (1.0 + 2.0 * VBAR * mean) / 0.001;
            let potential = 0.5 * (b / a).ln(); // (D_minus-D_plus)/(D_minus+D_plus).
            for n in flux {
                close(n, expected, 2e-11 * expected.abs());
            }
            close(
                report.dimensionless_potential_right_minus_left.unwrap(),
                potential,
                3e-11,
            );
            close(
                report.potential_right_minus_left_v.unwrap(),
                potential * R * T / FARADAY_C_PER_MOL,
                1e-12,
            );
            assert!(report.charge_flux_residual_mol_m2_s.abs() <= 1e-15);
            assert_eq!(report.friction_evaluations, 3);
            close(
                report.molecular.chemical_entropy_production_w_m3_k.unwrap(),
                -2.0 * R * expected * (b / a).ln() / 0.001,
                1e-9,
            );
            assert!(report.molecular.friction_entropy_production_w_m3_k.unwrap() >= 0.0);
        }
    }
}

#[test]
fn nonideal_divalent_salt_and_neutral_cosolute_match_independent_two_by_two_equations() {
    let interactions = vec![
        0.001, 0.0005, 0.008, 0.0005, 0.001, 0.003, 0.008, 0.003, 0.1,
    ];
    let model = source(1.0, interactions.clone(), vec![200.0; 3]);
    let (left, right) = ([20.0_f64, 40.0, 100.0], [21.0_f64, 42.0, 1.0]);
    let binary = [
        0.0, 1e-9, 3e-9, 2e-9, 1e-9, 0.0, 5e-10, 8e-10, 3e-9, 5e-10, 0.0, 1.2e-9, 2e-9, 8e-10,
        1.2e-9, 0.0,
    ];
    let material = SolutionMaxwellStefan::new(
        model.clone(),
        (SpeciesId::new("water").unwrap(), MolarMass::new(0.018)),
        (0..3)
            .map(|i| {
                (
                    SpeciesId::new(format!("solute-{i}")).unwrap(),
                    MolarMass::new(0.04 + 0.02 * i as f64),
                )
            })
            .collect(),
        binary.to_vec(),
    )
    .unwrap();
    let policy = ZeroCurrentDiffusion::new(vec![2, -1, 0], 1e-15, 1e-16).unwrap();
    let mut flux = [0.0; 3];
    material
        .flux_zero_current_into(
            cell(&model, &left).solution(),
            cell(&model, &right).solution(),
            0.001,
            physical().constitutive,
            &policy,
            &mut material.workspace().unwrap(),
            &mut flux,
            || false,
        )
        .unwrap();
    // Eliminate both charge and electric force analytically: N_anion=2*N_cation.
    // Solve only salt extent and neutral-solute flux in the independent basis.
    let mean = std::array::from_fn::<_, 3, _>(|i| (right[i] - left[i]) / (right[i] / left[i]).ln());
    let total = 1.0 / VBAR + mean.iter().sum::<f64>();
    let x = [
        1.0 / VBAR / total,
        mean[0] / total,
        mean[1] / total,
        mean[2] / total,
    ];
    let force = std::array::from_fn::<_, 3, _>(|i| {
        ((right[i] / left[i]).ln()
            + (0..3)
                .map(|j| interactions[3 * i + j] * (right[j] - left[j]))
                .sum::<f64>())
            / 0.001
    });
    let mut friction = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            friction[i][j] = if i == j {
                (0..4)
                    .filter(|&k| k != i + 1)
                    .map(|k| x[k] / binary[4 * (i + 1) + k])
                    .sum::<f64>()
                    / (total * x[i + 1])
            } else {
                -1.0 / (total * binary[4 * (i + 1) + j + 1])
            };
        }
    }
    let a = friction[0][0] + 2.0 * friction[0][1] + 2.0 * friction[1][0] + 4.0 * friction[1][1];
    let b = friction[0][2] + 2.0 * friction[1][2];
    let c = friction[2][0] + 2.0 * friction[2][1];
    let d = friction[2][2];
    let salt_force = -(force[0] + 2.0 * force[1]);
    let q = (salt_force * d + b * force[2]) / (a * d - b * c);
    let neutral = (-a * force[2] - c * salt_force) / (a * d - b * c);
    close(flux[0], q, 2e-14);
    close(flux[1], 2.0 * q, 2e-14);
    close(flux[2], neutral, 2e-14);
    assert!(
        q > 0.0,
        "salt moves uphill in concentration under the original nonideal chemical force"
    );
}

#[test]
fn independent_salts_use_a_neutral_intermediate_face_and_conjugate_current_constraint() {
    let model = source(1.0, vec![0.0; 9], vec![200.0; 3]);
    let (left, right) = ([100.0_f64, 20.0, 120.0], [20.0_f64, 80.0, 100.0]);
    let binary = [
        0.0, 1e-9, 2e-9, 3e-9, 1e-9, 0.0, 4e-10, 7e-10, 2e-9, 4e-10, 0.0, 9e-10, 3e-9, 7e-10,
        9e-10, 0.0,
    ];
    let material = SolutionMaxwellStefan::new(
        model.clone(),
        (SpeciesId::new("water").unwrap(), MolarMass::new(0.018)),
        (0..3)
            .map(|i| {
                (
                    SpeciesId::new(format!("solute-{i}")).unwrap(),
                    MolarMass::new(0.04 + 0.02 * i as f64),
                )
            })
            .collect(),
        binary.to_vec(),
    )
    .unwrap();
    let policy = ZeroCurrentDiffusion::new(vec![1, 1, -1], 1e-15, 1e-16).unwrap();
    let mut flux = [0.0; 3];
    let report = material
        .flux_zero_current_into(
            cell(&model, &left).solution(),
            cell(&model, &right).solution(),
            0.001,
            physical().constitutive,
            &policy,
            &mut material.workspace().unwrap(),
            &mut flux,
            || false,
        )
        .unwrap();
    // Independent Simpson integration of the neutral straight path weighted
    // by 1/total ionic concentration; no production mean helper is used.
    let mut integral = [0.0; 4];
    for k in 0..=1000 {
        let t = k as f64 / 1000.0;
        let state = std::array::from_fn::<_, 3, _>(|i| left[i] + t * (right[i] - left[i]));
        let weight = if k == 0 || k == 1000 {
            1.0
        } else if k % 2 == 0 {
            2.0
        } else {
            4.0
        };
        let total = state.iter().sum::<f64>();
        integral[3] += weight / total;
        for i in 0..3 {
            integral[i] += weight * state[i] / total;
        }
    }
    let face = std::array::from_fn::<_, 3, _>(|i| integral[i] / integral[3]);
    close(face[0] + face[1] - face[2], 0.0, 1e-12);
    let wrong =
        std::array::from_fn::<_, 3, _>(|i| (right[i] - left[i]) / (right[i] / left[i]).ln());
    assert!((wrong[0] + wrong[1] - wrong[2]).abs() > 10.0);
    let total = 1.0 / VBAR + face.iter().sum::<f64>();
    let x = [
        1.0 / VBAR / total,
        face[0] / total,
        face[1] / total,
        face[2] / total,
    ];
    let force = std::array::from_fn::<_, 3, _>(|i| (right[i] / left[i]).ln() / 0.001);
    let mut f = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            f[i][j] = if i == j {
                (0..4)
                    .filter(|&k| k != i + 1)
                    .map(|k| x[k] / binary[4 * (i + 1) + k])
                    .sum::<f64>()
                    / (total * x[i + 1])
            } else {
                -1.0 / (total * binary[4 * (i + 1) + j + 1])
            };
        }
    }
    // NaCl and KCl component fluxes; N_Cl=N_Na+N_K removes current exactly.
    let a = f[0][0] + f[0][2] + f[2][0] + f[2][2];
    let b = f[0][1] + f[0][2] + f[2][1] + f[2][2];
    let c = f[1][0] + f[1][2] + f[2][0] + f[2][2];
    let d = f[1][1] + f[1][2] + f[2][1] + f[2][2];
    let r0 = -(force[0] + force[2]);
    let r1 = -(force[1] + force[2]);
    let na = (r0 * d - b * r1) / (a * d - b * c);
    let k = (a * r1 - c * r0) / (a * d - b * c);
    for (actual, expected) in flux.into_iter().zip([na, k, na + k]) {
        close(actual, expected, 2e-14);
    }
    let potential = -(force[0] + f[0][0] * na + f[0][1] * k + f[0][2] * (na + k)) * 0.001;
    close(
        report.dimensionless_potential_right_minus_left.unwrap(),
        potential,
        1e-10,
    );
    assert!(report.molecular.chemical_entropy_production_w_m3_k.unwrap() > 0.0);
}

#[test]
fn ionic_face_refuses_singular_empty_fronts_and_preserves_outputs_on_late_cancellation() {
    let model = salt();
    let material = law(&model);
    let policy = ionic();
    let mut workspace = material.workspace().unwrap();
    let (left, right) = (cell(&model, &[100.0, 100.0]), cell(&model, &[20.0, 20.0]));
    let mut output = [123.0, 456.0];
    let mut polls = 0;
    material
        .flux_zero_current_into(
            left.solution(),
            right.solution(),
            0.001,
            physical().constitutive,
            &policy,
            &mut workspace,
            &mut output,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    output = [123.0, 456.0];
    let mut calls = 0;
    assert_eq!(
        material
            .flux_zero_current_into(
                left.solution(),
                right.solution(),
                0.001,
                physical().constitutive,
                &policy,
                &mut workspace,
                &mut output,
                || {
                    calls += 1;
                    calls == polls
                }
            )
            .unwrap_err(),
        SolutionDiffusionError::Cancelled
    );
    assert_eq!(output, [123.0, 456.0]);
    assert_eq!(
        material
            .flux_zero_current_into(
                left.solution(),
                cell(&model, &[0.0, 0.0]).solution(),
                0.001,
                physical().constitutive,
                &policy,
                &mut workspace,
                &mut output,
                || false
            )
            .unwrap_err(),
        SolutionDiffusionError::SingularDiffusionPotential
    );
    assert_eq!(output, [123.0, 456.0]);
    let empty = cell(&model, &[0.0, 0.0]);
    let report = material
        .flux_zero_current_into(
            empty.solution(),
            empty.solution(),
            0.001,
            physical().constitutive,
            &policy,
            &mut workspace,
            &mut output,
            || false,
        )
        .unwrap();
    assert_eq!(output, [0.0, 0.0]);
    assert_eq!(report.potential_right_minus_left_v, None);
    assert!(matches!(
        policy.admit_inventories(&[1.0, 0.9]),
        Err(SolutionDiffusionError::InventoryCharge { .. })
    ));
}

fn pair() -> OsmoticNetwork {
    network(&salt(), &[&[100.0, 100.0], &[20.0, 20.0]], &[(0, 1, 0.01)])
        .with_zero_current_diffusion(ionic(), 1 << 22)
        .unwrap()
}

#[test]
fn charged_network_refines_to_independent_salt_ode_and_rolls_back_whole_windows() {
    let rate = |left: f64| {
        let right = 120.0 - left;
        let mean = (left - right) / (left / right).ln();
        -0.01 / VOLUME * 1.5e-9 * (left - right) * (1.0 + 2.0 * VBAR * mean) / 0.001
    };
    let mut reference = 100.0;
    let h = 50.0 / 10000.0;
    for _ in 0..10000 {
        let a = rate(reference);
        let b = rate(reference + 0.5 * h * a);
        let c = rate(reference + 0.5 * h * b);
        let d = rate(reference + h * c);
        reference += h * (a + 2.0 * b + 2.0 * c + d) / 6.0;
    }
    let mut errors = Vec::new();
    for tolerance in [1e-3, 6.25e-5] {
        let mut system = pair();
        let report = system
            .advance(50.0, numerical(tolerance), || false)
            .unwrap();
        assert!(report.diffusion_free_energy_change_j.unwrap() < 0.0);
        for compartment in system.compartments() {
            ionic()
                .admit_inventories(compartment.solution().solute_moles())
                .unwrap();
        }
        for i in 0..2 {
            assert!(report.solute_balance_residual_mol[i].abs() < 1e-16);
        }
        errors.push((concentration(&system, 0, 0) - reference).abs());
    }
    assert!(errors[1] < 0.4 * errors[0], "{errors:?}");
    let mut success = pair();
    let mut polls = 0;
    success
        .advance(1.0, numerical(1e-3), || {
            polls += 1;
            false
        })
        .unwrap();
    let mut failed = pair();
    let saved = bits(&failed);
    let mut calls = 0;
    assert!(matches!(
        failed.advance(1.0, numerical(1e-3), || {
            calls += 1;
            calls == polls
        }),
        Err(OsmoticNetworkError::Cancelled)
    ));
    assert_eq!(bits(&failed), saved);
    failed.advance(1.0, numerical(1e-3), || false).unwrap();
    assert_eq!(bits(&failed), bits(&success));
    let mut exhausted = pair();
    let saved = bits(&exhausted);
    assert!(matches!(
        exhausted.advance(
            50.0,
            OsmoticStepControls {
                maximum_steps: 1,
                ..numerical(1e-5)
            },
            || false
        ),
        Err(OsmoticNetworkError::WorkBudget)
    ));
    assert_eq!(bits(&exhausted), saved);
}

#[test]
fn ionic_selection_rejects_cycles_non_neutrality_and_independent_implicit_coordinates() {
    let model = salt();
    let cells = [&[100.0, 100.0][..], &[60.0, 60.0][..], &[20.0, 20.0][..]];
    assert!(
        network(&model, &cells, &[(0, 1, 0.01), (1, 2, 0.01)])
            .with_zero_current_diffusion(ionic(), 1 << 22)
            .is_ok()
    );
    assert!(
        network(&model, &cells, &[(0, 1, 0.01), (1, 2, 0.01), (0, 2, 0.01)])
            .with_zero_current_diffusion(ionic(), 1 << 22)
            .is_err()
    );
    assert!(
        network(&model, &cells, &[(0, 1, 0.01), (1, 2, 0.01), (0, 2, 0.0)])
            .with_zero_current_diffusion(ionic(), 1 << 22)
            .is_ok()
    );
    assert!(
        network(&model, &[&[100.0, 99.0], &[20.0, 20.0]], &[(0, 1, 0.01)])
            .with_zero_current_diffusion(ionic(), 1 << 22)
            .is_err()
    );
    assert!(
        pair()
            .with_implicit_multicomponent_diffusion(
                ImplicitMulticomponentDiffusionControls::default(),
                1 << 22
            )
            .is_err()
    );
    let implicit = network(&model, &cells, &[(0, 1, 0.01)])
        .with_implicit_multicomponent_diffusion(
            ImplicitMulticomponentDiffusionControls::default(),
            1 << 22,
        )
        .unwrap();
    assert!(
        implicit
            .with_zero_current_diffusion(ionic(), 1 << 22)
            .is_err()
    );
}

#[test]
fn retained_water_and_ionic_diffusion_share_one_neutral_finite_network() {
    use crate::osmosis::{OsmoticMembrane, coupled::CoupledSolutionControls};
    use fs_material::osmosis::{
        HydraulicPermeance, MembraneTransport, ReflectionCoefficient, SolutePermeance,
    };
    let model = salt();
    let membrane = OsmoticMembrane::new(
        0,
        1,
        0.01,
        MembraneTransport::new(
            HydraulicPermeance::new(5e-12).unwrap(),
            vec![
                (
                    ReflectionCoefficient::new(1.0).unwrap(),
                    SolutePermeance::new(0.0).unwrap()
                );
                2
            ],
        )
        .unwrap(),
    )
    .unwrap();
    let mut system = OsmoticNetwork::with_coupled_solution_transport(
        vec!["solute-0".into(), "solute-1".into()],
        vec![cell(&model, &[100.0, 100.0]), cell(&model, &[20.0, 20.0])],
        vec![membrane],
        vec![SolutionDiffusionLink::new(0, 1, 0.01, 0.001, law(&model)).unwrap()],
        CoupledSolutionControls {
            diffusion: physical(),
            water_balance_tolerance_m3: 1e-18,
        },
        1 << 22,
    )
    .unwrap()
    .with_zero_current_diffusion(ionic(), 1 << 22)
    .unwrap();
    let report = system.advance(5.0, numerical(1e-4), || false).unwrap();
    assert!(system.compartments()[0].solution().solvent_volume_m3() > VOLUME);
    assert!(system.compartments()[0].solution().solute_moles()[0] < 100.0 * VOLUME);
    assert!(
        report
            .coupled_potential_change
            .unwrap()
            .total_potential_change_j
            < 0.0
    );
    assert!(report.water_balance_residual_m3.abs() < 1e-18);
    for cell in system.compartments() {
        ionic()
            .admit_inventories(cell.solution().solute_moles())
            .unwrap();
    }
}
