//! Actual crystal, solvent and molecular stocks share one native coupled clock.

use super::*;
use crate::osmosis::diffusion::ImplicitBinaryDiffusionControls;
use fs_material::osmosis::solubility::kinetics::KineticSolubilitySolution;
use fs_qty::{Amount, AmountConcentration, Temperature, Volume, VolumetricFlowRate};
use fs_thermochem::solution::solubility::{
    IdealSolubilityV1,
    kinetics::{PrecipitationPolicyV1, SolubilityKineticStateV1, SolubilityKineticsV1},
};

const SATURATION: f64 = 100.0;
const REFERENCE_SOLID: f64 = 200.0 * VL;
const CONDUCTANCE: f64 = 0.1 * VL;

fn ideal_model(reference: f64) -> Arc<VirialSolutionModel> {
    Arc::new(VirialSolutionModel::new(T, VBAR, reference, vec![500.0], vec![0.0]).unwrap())
}

fn phase_cell(
    concentration: f64,
    solid: f64,
    pressure: f64,
    particles: bool,
) -> OsmoticCompartment {
    let law = IdealSolubilityV1::new(
        id("water"),
        id("A"),
        Temperature::new(T),
        AmountConcentration::new(SATURATION),
        0.0,
        [Temperature::new(280.0), Temperature::new(340.0)],
    )
    .unwrap();
    let model = SolubilityKineticsV1::new(
        law,
        VolumetricFlowRate::new(CONDUCTANCE),
        PrecipitationPolicyV1::ExistingSolidOnly,
    )
    .unwrap();
    let model = if particles {
        model
            .with_particle_surface(Amount::new(REFERENCE_SOLID))
            .unwrap()
    } else {
        model
    };
    let phases = KineticSolubilitySolution::new(
        Temperature::new(T),
        VBAR,
        Volume::new(VL),
        vec![(
            model,
            SolubilityKineticStateV1::from_phases(
                Amount::new(concentration * VL),
                Amount::new(solid),
            )
            .unwrap(),
        )],
    )
    .unwrap();
    OsmoticCompartment::with_kinetic_solubility(
        phases,
        OsmoticPressureVolume::new(VL, pressure, K).unwrap(),
        0.1 * VL,
    )
    .unwrap()
}

fn prepared(reference: f64, particles: bool) -> OsmoticNetwork {
    let model = ideal_model(reference);
    OsmoticNetwork::with_kinetic_coupled_solution_transport(
        vec!["A".into()],
        vec![
            phase_cell(80.0, REFERENCE_SOLID, 1e5, particles),
            cell(&model, VR, 20.0, 1e5, K),
        ],
        vec![membrane(0, 1)],
        vec![SolutionDiffusionLink::new(0, 1, AREA, LENGTH, constant(&model)).unwrap()],
        model,
        physical(),
        1 << 22,
    )
    .unwrap()
}

fn complete_bits(system: &OsmoticNetwork) -> Vec<u64> {
    let mut result = bits(system);
    for cell in system.compartments() {
        if let Some(phases) = cell.kinetic_solubility_snapshot() {
            result.push(phases.solution().solvent_volume_m3().to_bits());
            for phase in phases.last_steps() {
                result.push(phase.state().solid().value().to_bits());
                result.push(phase.state().dissolved().value().to_bits());
                result.push(phase.state().total_solute().value().to_bits());
                result.push(phase.transferred_to_solution().value().to_bits());
            }
        }
    }
    result
}

fn phase(system: &OsmoticNetwork) -> KineticSolubilitySolution {
    system.compartments()[0]
        .kinetic_solubility_snapshot()
        .unwrap()
}

// Independent continuous equations for [left V, left dissolved n, right n,
// left solid n]. The native binary MS law in stationary solvent has
// J = D*(c_left-c_right)*(1+c_log*Vbar)/L. Water is a separate retaining path.
fn oracle(duration: f64, particles: bool) -> [f64; 4] {
    let rhs = |x: [f64; 4]| {
        let [vl, nl, nr, solid] = x;
        let vr = VL + VR - vl;
        let (cl, cr) = (nl / vl, nr / vr);
        let mean = if cl == cr {
            cl
        } else {
            (cl - cr) / (cl / cr).ln()
        };
        let molecular = AREA * 1e-9 * (cl - cr) * (1.0 + mean * VBAR) / LENGTH;
        let area = if particles {
            (solid / REFERENCE_SOLID).powf(2.0 / 3.0)
        } else {
            1.0
        };
        let dissolution = CONDUCTANCE * area * (SATURATION - cl);
        let pressure = K * (vl - VL) - K * (vr - VR);
        let water = -AREA * LP * (pressure - R * T * (cl - cr));
        [water, dissolution - molecular, molecular, -dissolution]
    };
    let add =
        |x: [f64; 4], rate: [f64; 4], scale: f64| core::array::from_fn(|i| x[i] + scale * rate[i]);
    let mut state = [VL, 80.0 * VL, 20.0 * VR, REFERENCE_SOLID];
    let h = duration / 16384.0;
    for _ in 0..16384 {
        let a = rhs(state);
        let b = rhs(add(state, a, h * 0.5));
        let c = rhs(add(state, b, h * 0.5));
        let d = rhs(add(state, c, h));
        for i in 0..4 {
            state[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    state
}

fn independent_chemical(volume: f64, dissolved: f64) -> f64 {
    if dissolved == 0.0 {
        0.0
    } else {
        R * T * dissolved * ((dissolved / (volume * SATURATION)).ln() - 1.0)
    }
}

#[test]
fn water_molecular_and_particle_exchange_refine_together_to_independent_ode() {
    for particles in [false, true] {
        let exact = oracle(5.0, particles);
        let mut errors = Vec::new();
        for tolerance in [1e-3, 6.25e-5] {
            let mut system = prepared(1.0, particles);
            let report = system.advance(5.0, temporal(tolerance), || false).unwrap();
            let phases = phase(&system);
            let left = system.compartments()[0].solution();
            let right = system.compartments()[1].solution();
            let actual = [
                left.solvent_volume_m3(),
                left.solute_moles()[0],
                right.solute_moles()[0],
                phases.last_steps()[0].state().solid().value(),
            ];
            errors.push(
                actual
                    .iter()
                    .zip(exact)
                    .zip([VL, 100.0 * VL, 100.0 * VR, REFERENCE_SOLID])
                    .map(|((&a, b), scale)| (a - b).abs() / scale)
                    .fold(0.0, f64::max),
            );
            assert!(actual[0] > VL && actual[2] > 20.0 * VR && actual[3] < REFERENCE_SOLID);
            assert_eq!(phases.solution().solvent_volume_m3(), actual[0]);
            assert_eq!(phases.solution().solute_moles(), left.solute_moles());
            assert!(report.phase_evaluations >= 6 * report.accepted_steps);
            assert!(
                report.water_balance_residual_m3.abs() <= physical().water_balance_tolerance_m3
            );
            assert!(
                report.solute_balance_residual_mol[0].abs()
                    <= physical().diffusion.amount_balance_tolerance_mol
            );
            near(
                actual[1] + actual[2] + actual[3],
                80.0 * VL + 20.0 * VR + REFERENCE_SOLID,
                1e-16,
            );
            let chemical = independent_chemical(actual[0], actual[1])
                + independent_chemical(right.solvent_volume_m3(), actual[2])
                - independent_chemical(VL, 80.0 * VL)
                - independent_chemical(VR, 20.0 * VR);
            let mechanical =
                0.5 * K * ((actual[0] - VL).powi(2) + (right.solvent_volume_m3() - VR).powi(2));
            let account = report.coupled_potential_change.unwrap();
            near(account.mixing_free_energy_change_j, chemical, 2e-12);
            near(account.mechanical_potential_change_j, mechanical, 2e-12);
            assert!(account.total_potential_change_j < 0.0);
            assert!(report.diffusion_free_energy_change_j.is_none());
            if particles {
                assert_eq!(
                    phases.models()[0].particle_reference_solid(),
                    Some(Amount::new(REFERENCE_SOLID))
                );
            }
        }
        assert!(
            errors[0] < 0.02 && errors[1] < 0.65 * errors[0],
            "{errors:?}"
        );
    }
}

#[test]
fn pressure_driven_water_removal_grows_crystals_and_retains_load_work() {
    let model = ideal_model(1.0);
    let mut system = OsmoticNetwork::with_kinetic_coupled_solution_transport(
        vec!["A".into()],
        vec![
            phase_cell(SATURATION, REFERENCE_SOLID, 1e6, true),
            cell(&model, VR, SATURATION, 1e5, K),
        ],
        vec![membrane(0, 1)],
        Vec::new(),
        model,
        physical(),
        1 << 22,
    )
    .unwrap();
    let report = system.advance(1.0, temporal(1e-5), || false).unwrap();
    let phases = phase(&system);
    let left = system.compartments()[0].solution();
    let right = system.compartments()[1].solution();
    assert!(left.solvent_volume_m3() < VL && right.solvent_volume_m3() > VR);
    assert!(phases.last_steps()[0].state().solid().value() > REFERENCE_SOLID);
    near(
        phases.last_steps()[0].state().total_solute().value(),
        SATURATION * VL + REFERENCE_SOLID,
        1e-17,
    );
    assert_eq!(right.solute_moles()[0], SATURATION * VR);
    assert_eq!(
        phases.solution().solvent_volume_m3(),
        left.solvent_volume_m3()
    );
    let chemical = independent_chemical(left.solvent_volume_m3(), left.solute_moles()[0])
        + independent_chemical(right.solvent_volume_m3(), right.solute_moles()[0])
        - independent_chemical(VL, SATURATION * VL)
        - independent_chemical(VR, SATURATION * VR);
    let mechanical = (1e6 - 1e5) * (left.solvent_volume_m3() - VL)
        + 0.5
            * K
            * ((left.solvent_volume_m3() - VL).powi(2) + (right.solvent_volume_m3() - VR).powi(2));
    let account = report.coupled_potential_change.unwrap();
    near(account.mixing_free_energy_change_j, chemical, 2e-12);
    near(account.mechanical_potential_change_j, mechanical, 2e-12);
    assert!(account.mixing_free_energy_change_j > 0.0);
    assert!(account.mechanical_potential_change_j < 0.0 && account.total_potential_change_j < 0.0);
}

#[test]
fn reference_changes_and_existing_implicit_stages_keep_original_phase_stocks() {
    let mut ordinary = prepared(1.0, true);
    let mut shifted = prepared(1e200, true);
    let a = ordinary.advance(2.0, temporal(1e-4), || false).unwrap();
    let b = shifted.advance(2.0, temporal(1e-4), || false).unwrap();
    assert_eq!(complete_bits(&ordinary), complete_bits(&shifted));
    assert_eq!(a.coupled_potential_change, b.coupled_potential_change);
    let mut implicit = prepared(1.0, true)
        .with_implicit_binary_diffusion(ImplicitBinaryDiffusionControls::default(), 1 << 22)
        .unwrap();
    let report = implicit.advance(2.0, temporal(1e-5), || false).unwrap();
    let exact = oracle(2.0, true);
    let phases = phase(&implicit);
    near(
        implicit.compartments()[0].solution().solvent_volume_m3(),
        exact[0],
        3e-10,
    );
    near(
        implicit.compartments()[0].solution().solute_moles()[0],
        exact[1],
        3e-8,
    );
    near(
        phases.last_steps()[0].state().solid().value(),
        exact[3],
        3e-8,
    );
    assert!(
        report
            .coupled_potential_change
            .unwrap()
            .total_potential_change_j
            < 0.0
    );
    assert!(
        report.solute_balance_residual_mol[0].abs()
            <= physical().diffusion.amount_balance_tolerance_mol
    );
}

#[test]
fn late_cancellation_and_molecular_work_refusal_restore_water_and_both_phases() {
    let mut system = prepared(1.0, true)
        .with_implicit_binary_diffusion(
            ImplicitBinaryDiffusionControls {
                maximum_evaluations: 1,
                ..ImplicitBinaryDiffusionControls::default()
            },
            1 << 22,
        )
        .unwrap();
    let before = complete_bits(&system);
    assert!(matches!(
        system.advance(1.0, temporal(1e-4), || false),
        Err(OsmoticNetworkError::WorkBudget)
    ));
    assert_eq!(complete_bits(&system), before);
    system = system
        .with_implicit_binary_diffusion(ImplicitBinaryDiffusionControls::default(), 1 << 22)
        .unwrap();
    let mut polls = 0;
    assert!(matches!(
        system.advance(2.0, temporal(1e-5), || {
            polls += 1;
            polls == 80
        }),
        Err(OsmoticNetworkError::Cancelled)
    ));
    assert_eq!(complete_bits(&system), before);
    let mut fresh = prepared(1.0, true)
        .with_implicit_binary_diffusion(ImplicitBinaryDiffusionControls::default(), 1 << 22)
        .unwrap();
    system.advance(1.0, temporal(1e-5), || false).unwrap();
    fresh.advance(1.0, temporal(1e-5), || false).unwrap();
    assert_eq!(complete_bits(&system), complete_bits(&fresh));
}

#[test]
fn unsupported_activity_projection_and_leaky_water_paths_refuse() {
    let nonideal = model(1.0);
    assert!(
        OsmoticNetwork::with_kinetic_coupled_solution_transport(
            vec!["A".into()],
            vec![phase_cell(80.0, REFERENCE_SOLID, 1e5, true)],
            Vec::new(),
            Vec::new(),
            nonideal,
            physical(),
            1 << 22,
        )
        .is_err()
    );
    let model = ideal_model(1.0);
    let leaky = OsmoticMembrane::new(
        0,
        1,
        AREA,
        MembraneTransport::new(
            HydraulicPermeance::new(LP).unwrap(),
            vec![(
                ReflectionCoefficient::new(0.5).unwrap(),
                SolutePermeance::new(1e-9).unwrap(),
            )],
        )
        .unwrap(),
    )
    .unwrap();
    assert!(
        OsmoticNetwork::with_kinetic_coupled_solution_transport(
            vec!["A".into()],
            vec![
                phase_cell(80.0, REFERENCE_SOLID, 1e5, true),
                cell(&model, VR, 20.0, 1e5, K)
            ],
            vec![leaky],
            Vec::new(),
            model,
            physical(),
            1 << 22,
        )
        .is_err()
    );
}

#[test]
fn explicitly_bound_nonideal_crystals_share_water_forces_and_complete_potential() {
    let model = Arc::new(
        VirialSolutionModel::new(
            T,
            VBAR,
            1.0,
            vec![300.0; 2],
            vec![0.002, 0.003, 0.003, 0.001],
        )
        .unwrap(),
    );
    let initial = [[80.0, 40.0], [20.0, 10.0]];
    let components = ["A", "B"]
        .into_iter()
        .enumerate()
        .map(|(i, name)| {
            let law = IdealSolubilityV1::new(
                id("water"),
                id(name),
                Temperature::new(T),
                AmountConcentration::new([100.0, 20.0][i]),
                0.0,
                [Temperature::new(280.0), Temperature::new(340.0)],
            )
            .unwrap();
            let stock = Amount::new([200.0, 60.0][i] * VL);
            let kinetics = SolubilityKineticsV1::new(
                law,
                VolumetricFlowRate::new([0.1, 0.05][i] * VL),
                PrecipitationPolicyV1::ExistingSolidOnly,
            )
            .unwrap()
            .with_particle_surface(stock)
            .unwrap();
            (
                kinetics,
                SolubilityKineticStateV1::from_phases(Amount::new(initial[0][i] * VL), stock)
                    .unwrap(),
            )
        })
        .collect();
    let original =
        KineticSolubilitySolution::new(Temperature::new(T), VBAR, Volume::new(VL), components)
            .unwrap()
            .with_nonideal_thermodynamics(model.clone())
            .unwrap();
    let left = OsmoticCompartment::with_kinetic_solubility(
        original,
        OsmoticPressureVolume::new(VL, 1e5, K).unwrap(),
        0.1 * VL,
    )
    .unwrap();
    let right = OsmoticCompartment::new(
        AqueousSolution::nonideal(model.clone(), VR, initial[1].map(|c| c * VR).to_vec()).unwrap(),
        OsmoticPressureVolume::new(VR, 1e5, K).unwrap(),
        0.1 * VR,
    )
    .unwrap();
    let water = OsmoticMembrane::new(
        0,
        1,
        AREA,
        MembraneTransport::new(
            HydraulicPermeance::new(LP).unwrap(),
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
    let molecular = SolutionMaxwellStefan::new(
        model.clone(),
        (id("water"), MolarMass::new(0.018)),
        vec![
            (id("A"), MolarMass::new(0.04)),
            (id("B"), MolarMass::new(0.06)),
        ],
        vec![0.0, 1e-9, 2e-9, 1e-9, 0.0, 3e-10, 2e-9, 3e-10, 0.0],
    )
    .unwrap();
    let mut system = OsmoticNetwork::with_kinetic_coupled_solution_transport(
        vec!["A".into(), "B".into()],
        vec![left, right],
        vec![water],
        vec![SolutionDiffusionLink::new(0, 1, AREA, LENGTH, molecular).unwrap()],
        model.clone(),
        physical(),
        1 << 24,
    )
    .unwrap();
    let report = system.advance(1.0, temporal(1e-4), || false).unwrap();
    let retained = phase(&system);
    assert_eq!(retained.solution().nonideal_model(), Some(model.as_ref()));
    assert!(retained.last_steps()[0].state().solid().value() < 200.0 * VL);
    assert!(retained.last_steps()[1].state().solid().value() > 60.0 * VL);
    let left = system.compartments()[0].solution();
    let right = system.compartments()[1].solution();
    assert!(left.solvent_volume_m3() > VL && right.solvent_volume_m3() < VR);
    assert_eq!(
        retained.solution().solvent_volume_m3(),
        left.solvent_volume_m3()
    );
    assert!(report.water_balance_residual_m3.abs() <= physical().water_balance_tolerance_m3);
    for (i, total) in [320.0 * VL, 120.0 * VL].into_iter().enumerate() {
        let actual = retained.last_steps()[i].state().solid().value()
            + left.solute_moles()[i]
            + right.solute_moles()[i];
        near(actual, total, 1e-16);
        assert!(
            report.solute_balance_residual_mol[i].abs()
                <= physical().diffusion.amount_balance_tolerance_mol
        );
    }
    let energy = |volume: f64, n: &[f64]| {
        let c = [n[0] / volume, n[1] / volume];
        let ideal = c
            .iter()
            .zip([100.0, 20.0])
            .zip([0.002, 0.001])
            .map(|((&c, saturation), diagonal)| {
                c * ((c / saturation).ln() - 1.0 - diagonal * saturation)
            })
            .sum::<f64>();
        R * T
            * volume
            * (ideal
                + 0.5 * (0.002 * c[0] * c[0] + 2.0 * 0.003 * c[0] * c[1] + 0.001 * c[1] * c[1]))
    };
    let chemical = energy(left.solvent_volume_m3(), left.solute_moles())
        + energy(right.solvent_volume_m3(), right.solute_moles())
        - energy(VL, &initial[0].map(|c| c * VL))
        - energy(VR, &initial[1].map(|c| c * VR));
    let mechanical = 0.5
        * K
        * ((left.solvent_volume_m3() - VL).powi(2) + (right.solvent_volume_m3() - VR).powi(2));
    let account = report.coupled_potential_change.unwrap();
    near(account.mixing_free_energy_change_j, chemical, 3e-12);
    near(account.mechanical_potential_change_j, mechanical, 3e-12);
    near(
        account.total_potential_change_j,
        chemical + mechanical,
        3e-12,
    );
    assert!(account.total_potential_change_j < 0.0);
}
