//! Original bulk-law oracles and finite conservative osmotic consumers.
use super::*;
use crate::osmosis::{
    OsmoticCompartment, OsmoticMembrane, OsmoticNetwork, OsmoticStepControls,
    donnan::{DonnanOsmoticPair, DonnanStepControls},
    transient::{TransientFilm, TransientOsmoticPair, TransientStepControls},
};
use fs_material::{
    hyper::{
        Hyperelastic, HyperelasticModel, OgdenTerm, Volumetric,
        discrete_gradient::bulk::DiscreteBulkSettings, reference::ReferenceDeformation,
    },
    osmosis::{
        HydraulicPermeance, IdealDiluteSolution, MembraneTransport, ReflectionCoefficient,
        SolutePermeance, donnan::IdealDonnan,
    },
};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

const T: f64 = 298.15;
const RT: f64 = R * T;

fn settings() -> DiscreteBulkSettings {
    DiscreteBulkSettings {
        minimum_volume_ratio: 0.2,
        maximum_volume_ratio: 2.6,
        maximum_iterations: 96,
        absolute_pressure_tolerance_pa: 1e-10,
        relative_pressure_tolerance: 1e-13,
        absolute_work_tolerance_pa: 1e-10,
        relative_work_tolerance: 1e-13,
    }
}
fn material(volume: Volumetric) -> Hyperelastic {
    Hyperelastic::new(
        HyperelasticModel::Ogden {
            terms: vec![OgdenTerm {
                mu_pa: 3.0,
                alpha: 2.0,
            }],
            volumetric: volume,
        },
        4.0,
    )
    .unwrap()
}

// Independent original U, U' and U'' use platform ln, not material queries.
fn oracle(volume: Volumetric, j: f64) -> (f64, f64, f64) {
    let k = volume.kappa();
    match volume {
        Volumetric::Quadratic { .. } => (0.5 * k * (j - 1.0).powi(2), k * (j - 1.0), k),
        Volumetric::LogSquared { .. } => (
            0.5 * k * j.ln().powi(2),
            k * j.ln() / j,
            k * (1.0 - j.ln()) / (j * j),
        ),
        Volumetric::SimoTaylor { .. } => (
            0.25 * k * (j * j - 1.0 - 2.0 * j.ln()),
            0.5 * k * (j - 1.0 / j),
            0.5 * k * (1.0 + 1.0 / (j * j)),
        ),
    }
}
fn near(actual: f64, expected: f64, allowance: f64) {
    assert!(
        (actual - expected).abs() <= allowance,
        "actual={actual:e}, expected={expected:e}, allowance={allowance:e}"
    );
}
fn root(mut lo: f64, mut hi: f64, f: impl Fn(f64) -> f64) -> f64 {
    assert!(f(lo) < 0.0 && f(hi) > 0.0);
    for _ in 0..80 {
        let mid = lo + 0.5 * (hi - lo);
        if f(mid) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    lo + 0.5 * (hi - lo)
}
fn response(volume: Volumetric, limits: DiscreteBulkSettings) -> OsmoticPressureVolume {
    let point = DiscreteBulkPoint::try_new(&material(volume), 1.0, limits).unwrap();
    OsmoticPressureVolume::from_bulk_material(1.0, 0.0, 0.0, point).unwrap()
}
fn cell(amount: f64, pressure: OsmoticPressureVolume) -> OsmoticCompartment {
    OsmoticCompartment::new(
        IdealDiluteSolution::new(T, 18e-6, 1.0, vec![amount]).unwrap(),
        pressure,
        0.05,
    )
    .unwrap()
}
fn membrane() -> MembraneTransport {
    MembraneTransport::new(
        HydraulicPermeance::new(1.0 / RT).unwrap(),
        vec![(
            ReflectionCoefficient::new(1.0).unwrap(),
            SolutePermeance::new(0.0).unwrap(),
        )],
    )
    .unwrap()
}
fn charged_pair(pressure: OsmoticPressureVolume) -> DonnanOsmoticPair {
    DonnanOsmoticPair::new(
        [cell(1.0, pressure), cell(0.0, pressure)],
        IdealDonnan::new(vec![1]).unwrap(),
        [-1.0, 0.0],
        [0.1, 0.0],
        1.0,
        HydraulicPermeance::new(1.0 / RT).unwrap(),
        1 << 20,
    )
    .unwrap()
}

#[test]
fn all_original_bulk_laws_retain_pressure_tangent_and_exact_energy() {
    for volume in [
        Volumetric::Quadratic { kappa: 7e5 },
        Volumetric::LogSquared { kappa: 7e5 },
        Volumetric::SimoTaylor { kappa: 7e5 },
    ] {
        let point = DiscreteBulkPoint::try_new(&material(volume), 1.3, settings()).unwrap();
        let initial = 1.5e-6;
        let v0 = (initial + 2.5e-6) / 1.3;
        let pressure =
            OsmoticPressureVolume::from_bulk_material(initial, 2.5e-6, 1e5, point).unwrap();
        let shifted =
            OsmoticPressureVolume::from_bulk_material(initial, 2.5e-6, 1e24, point).unwrap();
        for v in [0.9e-6, initial, 2.1e-6] {
            let (_, p, tangent) = oracle(volume, 1.3 + (v - initial) / v0);
            near(pressure.pressure_pa(v).unwrap(), 1e5 + p, 2e-9);
            near(pressure.tangent_pa_per_m3(v).unwrap(), tangent / v0, 0.002);
        }
        let (a, b) = (0.9e-6, 2.1e-6);
        let expected = v0
            * (oracle(volume, 1.3 + (b - initial) / v0).0
                - oracle(volume, 1.3 + (a - initial) / v0).0);
        near(
            pressure.reference_cancelled_potential_change_j(a, b, 1e5).unwrap(),
            expected,
            2e-14,
        );
        near(
            shifted.reference_cancelled_potential_change_j(a, b, 1e24).unwrap(),
            expected,
            2e-14,
        );
        near(
            pressure.reference_cancelled_potential_change_j(b, a, 1e5).unwrap(),
            -expected,
            2e-14,
        );
        near(
            shifted.pressure_difference_pa(b, shifted, a).unwrap(),
            oracle(volume, 1.3 + (b - initial) / v0).1
                - oracle(volume, 1.3 + (a - initial) / v0).1,
            2e-9,
        );
    }
}

#[test]
fn dry_volume_original_prestrain_and_held_natural_reference_keep_their_chain_rule() {
    let volume = Volumetric::LogSquared { kappa: 8e5 };
    let reference =
        ReferenceDeformation::try_new([2.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]).unwrap();
    let point = DiscreteBulkPoint::try_new_with_reference(
        &material(volume),
        &reference,
        2.6,
        settings(),
    )
    .unwrap();
    let pressure = OsmoticPressureVolume::from_bulk_material(1e-6, 3e-6, 17.0, point).unwrap();
    let v0 = 4e-6 / 2.6;
    let (before, after) = (1e-6, 2e-6);
    let je_before = 1.3;
    let je_after = (2.6 + (after - before) / v0) / 2.0;
    assert!(pressure.pressure_pa(before).unwrap() > 1e5);
    near(
        pressure.pressure_pa(after).unwrap(),
        17.0 + oracle(volume, je_after).1,
        2e-9,
    );
    near(
        pressure.tangent_pa_per_m3(after).unwrap(),
        oracle(volume, je_after).2 / (2.0 * v0),
        0.002,
    );
    near(
        pressure.reference_cancelled_potential_change_j(before, after, 0.0).unwrap(),
        17.0 * (after - before)
            + v0 * 2.0 * (oracle(volume, je_after).0 - oracle(volume, je_before).0),
        2e-14,
    );
}

#[test]
fn donnan_nonlinear_swelling_solves_the_original_bulk_force_and_exact_work() {
    let volume = Volumetric::SimoTaylor { kappa: RT };
    let pressure = response(volume, settings());
    let mut pair = charged_pair(pressure);
    let dt = 0.5;
    let expected = root(-0.75, 0.0, |x| {
        x - dt / RT
            * (oracle(volume, 1.0 - x).1 - oracle(volume, 1.0 + x).1
                - RT * 1.1 / (1.0 - x))
    });
    let report = pair.step(dt, DonnanStepControls::default(), || false).unwrap();
    near(report.solvent_transfer_m3, expected, 2e-11);
    assert!(expected < -0.1);
    assert_eq!(pair.compartments()[0].solution().solute_moles(), &[1.0]);
    assert_eq!(pair.compartments()[1].solution().solute_moles(), &[0.0]);
    let vl = pair.compartments()[0].solution().solvent_volume_m3();
    let vr = pair.compartments()[1].solution().solvent_volume_m3();
    let exact_work = -oracle(volume, vl).0 - oracle(volume, vr).0;
    near(report.pressure_work_into_solution_j, exact_work, 2e-10);
    let trapezoid = 0.5 * (oracle(volume, vl).1 - oracle(volume, vr).1) * expected;
    assert!(
        (trapezoid - exact_work).abs() > 0.01,
        "nonlinear work must distinguish the original primitive from a trapezoid"
    );
    near(vl + vr, 2.0, 1e-14);
    assert!(RT * 1.1 * (-vl.ln()) - exact_work < 0.0);
}

#[test]
fn retaining_network_reaches_the_material_equilibrium_instead_of_a_tangent_spring() {
    let volume = Volumetric::SimoTaylor { kappa: RT };
    let pressure = response(volume, settings());
    let mut network = OsmoticNetwork::new(
        vec!["neutral".into()],
        vec![cell(1.0, pressure), cell(2.0, pressure)],
        vec![OsmoticMembrane::new(0, 1, 1.0, membrane()).unwrap()],
        1 << 20,
    )
    .unwrap();
    let balance = |x: f64, law| {
        -(oracle(law, 1.0 - x).1 - oracle(law, 1.0 + x).1)
            + RT * (1.0 / (1.0 - x) - 2.0 / (1.0 + x))
    };
    let expected = root(0.0, 0.6, |x| balance(x, volume));
    let linearized = root(0.0, 0.6, |x| balance(x, Volumetric::Quadratic { kappa: RT }));
    assert!((expected - linearized).abs() > 1e-4);
    let report = network
        .advance(
            6.0,
            OsmoticStepControls {
                maximum_step_s: 0.5,
                relative_tolerance: 1e-5,
                absolute_volume_tolerance_m3: 1e-12,
                absolute_amount_tolerance_mol: 1e-12,
                maximum_steps: 20_000,
                maximum_rejections: 128,
                ..Default::default()
            },
            || false,
        )
        .unwrap();
    near(network.compartments()[0].solution().solvent_volume_m3(), 1.0 - expected, 3e-5);
    near(report.water_balance_residual_m3, 0.0, 1e-12);
    assert_eq!(report.solute_balance_residual_mol[0], 0.0);
}

#[test]
fn finite_films_report_original_nonlinear_mechanical_storage_change() {
    let volume = Volumetric::SimoTaylor { kappa: RT };
    let pressure = response(volume, settings());
    let mut pair = TransientOsmoticPair::new(
        [cell(1.0, pressure), cell(2.0, pressure)],
        1.0,
        membrane(),
        Some(TransientFilm::new(0.05, 3, vec![0.2], vec![1.0; 3]).unwrap()),
        Some(TransientFilm::new(0.05, 3, vec![0.2], vec![2.0; 3]).unwrap()),
        1 << 20,
    )
    .unwrap();
    let before_water = pair.total_solvent_volume_m3();
    let report = pair.step(0.5, TransientStepControls::default(), || false).unwrap();
    let vl = pair.bulks()[0].solution().solvent_volume_m3();
    let vr = pair.bulks()[1].solution().solvent_volume_m3();
    let exact_work = -oracle(volume, vl).0 - oracle(volume, vr).0;
    assert!(report.solvent_transfer_m3 > 0.05);
    near(report.pressure_work_into_solution_j, exact_work, 2e-10);
    let trapezoid =
        0.5 * (oracle(volume, vl).1 - oracle(volume, vr).1) * report.solvent_transfer_m3;
    assert!((trapezoid - exact_work).abs() > 1e-3);
    near(pair.total_solvent_volume_m3(), before_water, 1e-14);
    near(report.solute_balance_residual_mol[0], 0.0, 1e-10);
}

#[test]
fn original_material_domain_refusal_preserves_donnan_inventories_and_clock() {
    let volume = Volumetric::LogSquared { kappa: RT };
    let pressure = response(
        volume,
        DiscreteBulkSettings {
            minimum_volume_ratio: 0.98,
            maximum_volume_ratio: 1.02,
            ..settings()
        },
    );
    assert!(matches!(pressure.pressure_pa(0.9), Err(OsmoticNetworkError::Mechanical(_))));
    assert!(matches!(pressure.tangent_pa_per_m3(1.1), Err(OsmoticNetworkError::Mechanical(_))));
    assert!(matches!(
        pressure.reference_cancelled_potential_change_j(1.0, 1.1, 0.0),
        Err(OsmoticNetworkError::Mechanical(_))
    ));
    let mut pair = charged_pair(pressure);
    let snapshot = |pair: &DonnanOsmoticPair| {
        let mut values = vec![pair.time_s().to_bits()];
        for c in pair.compartments() {
            values.push(c.solution().solvent_volume_m3().to_bits());
            values.extend(c.solution().solute_moles().iter().map(|x| x.to_bits()));
        }
        values.extend(pair.equilibrium().osmotic_pressure_pa.map(f64::to_bits));
        values.extend(pair.equilibrium().charge_residual_mol.map(f64::to_bits));
        values
    };
    let before = snapshot(&pair);
    assert!(pair.step(10.0, DonnanStepControls::default(), || false).is_err());
    assert_eq!(snapshot(&pair), before);
    let retry = pair.step(0.001, DonnanStepControls::default(), || false).unwrap();
    assert!(retry.solvent_transfer_m3 < 0.0);
    assert_eq!(pair.time_s(), 0.001);
}

#[test]
fn nonideal_total_potential_retains_original_bulk_storage_and_pressure_gauge() {
    use crate::osmosis::{
        coupled::CoupledSolutionControls, diffusion::SolutionDiffusionControls,
    };
    use fs_material::osmosis::{
        AqueousSolution, diffusion::SolutionDiffusionSolve, nonideal::VirialSolutionModel,
    };
    use fs_thermochem::transport::maxwell_stefan::MaxwellStefanMolarSolveV1;
    use std::sync::Arc;

    let volume = Volumetric::SimoTaylor { kappa: RT };
    let interaction = 0.2;
    let model = Arc::new(
        VirialSolutionModel::new(T, 18e-6, 1.0, vec![10.0], vec![interaction]).unwrap(),
    );
    let point = DiscreteBulkPoint::try_new(&material(volume), 1.0, settings()).unwrap();
    let mut reference = None;
    for gauge in [0.0, 1e24] {
        let pressure = OsmoticPressureVolume::from_bulk_material(1.0, 0.0, gauge, point).unwrap();
        let cells = [1.0, 2.0]
            .into_iter()
            .map(|amount| {
                OsmoticCompartment::new(
                    AqueousSolution::nonideal(model.clone(), 1.0, vec![amount]).unwrap(),
                    pressure,
                    0.05,
                )
                .unwrap()
            })
            .collect();
        let mut network = OsmoticNetwork::with_coupled_solution_transport(
            vec!["neutral".into()],
            cells,
            vec![OsmoticMembrane::new(0, 1, 1.0, membrane()).unwrap()],
            vec![],
            CoupledSolutionControls {
                diffusion: SolutionDiffusionControls {
                    constitutive: SolutionDiffusionSolve {
                        friction: MaxwellStefanMolarSolveV1::default(),
                        entropy_tolerance_w_m3_k: 1e-10,
                    },
                    maximum_species_fraction: 0.2,
                    amount_balance_tolerance_mol: 1e-12,
                    free_energy_tolerance_j: 1e-9,
                },
                water_balance_tolerance_m3: 1e-12,
            },
            1 << 20,
        )
        .unwrap();
        let report = network
            .advance(
                0.5,
                OsmoticStepControls {
                    maximum_step_s: 0.25,
                    relative_tolerance: 1e-5,
                    absolute_volume_tolerance_m3: 1e-12,
                    absolute_amount_tolerance_mol: 1e-12,
                    maximum_steps: 20_000,
                    maximum_rejections: 128,
                    ..Default::default()
                },
                || false,
            )
            .unwrap();
        let vl = network.compartments()[0].solution().solvent_volume_m3();
        let vr = network.compartments()[1].solution().solvent_volume_m3();
        let mechanical = oracle(volume, vl).0 + oracle(volume, vr).0;
        // Original F/RT = sum n(ln(n/V)-1) + B sum n^2/(2V).
        // Conserved n makes the reference terms cancel exactly.
        let mixing = RT
            * (-vl.ln() - 2.0 * vr.ln()
                + 0.5 * interaction * ((1.0 / vl - 1.0) + 4.0 * (1.0 / vr - 1.0)));
        let change = report.coupled_potential_change.unwrap();
        assert!(vl < 0.9 && mechanical > 1.0 && mixing < -mechanical);
        near(change.mechanical_potential_change_j, mechanical, 2e-10);
        near(change.mixing_free_energy_change_j, mixing, 2e-10);
        near(change.total_potential_change_j, mechanical + mixing, 3e-10);
        near(report.water_balance_residual_m3, 0.0, 1e-12);
        assert_eq!(report.solute_balance_residual_mol[0], 0.0);
        let bits = [
            vl.to_bits(),
            vr.to_bits(),
            change.mechanical_potential_change_j.to_bits(),
            change.mixing_free_energy_change_j.to_bits(),
            change.total_potential_change_j.to_bits(),
        ];
        if let Some(expected) = reference {
            assert_eq!(bits, expected);
        } else {
            reference = Some(bits);
        }
    }
}
