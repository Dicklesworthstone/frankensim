//! Independent coupled limits, actual phase conservation, and atomic retry.

use super::*;
use fs_material::osmosis::{HydraulicPermeance, ReflectionCoefficient, SolutePermeance};
use fs_qty::{Amount, AmountConcentration, Temperature, Volume, VolumetricFlowRate};
use fs_thermochem::{
    SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
    solution::solubility::{
        IdealSolubilityV1,
        kinetics::{PrecipitationPolicyV1, SolubilityKineticStateV1, SolubilityKineticsV1},
    },
};

const T: f64 = 298.15;

fn law(saturation: f64) -> IdealSolubilityV1 {
    IdealSolubilityV1::new(
        SpeciesId::new("water").unwrap(),
        SpeciesId::new("solute").unwrap(),
        Temperature::new(T),
        AmountConcentration::new(saturation),
        0.0,
        [Temperature::new(280.0), Temperature::new(350.0)],
    )
    .unwrap()
}

fn kinetic(
    dissolved: f64,
    solid: f64,
    conductance: f64,
    policy: PrecipitationPolicyV1,
) -> OsmoticCompartment {
    let phases = KineticSolubilitySolution::new(
        Temperature::new(T),
        18e-6,
        Volume::new(1.0),
        vec![(
            SolubilityKineticsV1::new(law(1.0), VolumetricFlowRate::new(conductance), policy)
                .unwrap(),
            SolubilityKineticStateV1::from_phases(Amount::new(dissolved), Amount::new(solid))
                .unwrap(),
        )],
    )
    .unwrap();
    OsmoticCompartment::with_kinetic_solubility(
        phases,
        OsmoticPressureVolume::new(1.0, 0.0, 0.0).unwrap(),
        1e-6,
    )
    .unwrap()
}

fn aqueous(dissolved: f64) -> OsmoticCompartment {
    OsmoticCompartment::new(
        IdealDiluteSolution::new(T, 18e-6, 1.0, vec![dissolved]).unwrap(),
        OsmoticPressureVolume::new(1.0, 0.0, 0.0).unwrap(),
        1e-6,
    )
    .unwrap()
}

fn membrane(hydraulic: f64, reflection: f64, permeance: f64) -> OsmoticMembrane {
    OsmoticMembrane::new(
        0,
        1,
        1.0,
        MembraneTransport::new(
            HydraulicPermeance::new(hydraulic).unwrap(),
            vec![(
                ReflectionCoefficient::new(reflection).unwrap(),
                SolutePermeance::new(permeance).unwrap(),
            )],
        )
        .unwrap(),
    )
    .unwrap()
}

fn network(cells: Vec<OsmoticCompartment>, membranes: Vec<OsmoticMembrane>) -> OsmoticNetwork {
    OsmoticNetwork::new(vec!["solute".into()], cells, membranes, 1 << 20).unwrap()
}

fn breakthrough(conductance: f64) -> OsmoticNetwork {
    network(
        vec![
            kinetic(
                0.0,
                10.0,
                conductance,
                PrecipitationPolicyV1::ExistingSolidOnly,
            ),
            aqueous(0.0),
        ],
        vec![membrane(0.0, 0.0, 0.4)],
    )
}

fn controls(step: f64) -> OsmoticStepControls {
    OsmoticStepControls {
        maximum_step_s: step,
        relative_tolerance: 1e-6,
        absolute_volume_tolerance_m3: 1e-11,
        absolute_amount_tolerance_mol: 1e-11,
        maximum_steps: 50_000,
        maximum_rejections: 128,
        ..OsmoticStepControls::default()
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.16e}, expected={expected:.16e}, tolerance={tolerance:.3e}"
    );
}

fn phase(cell: &OsmoticCompartment) -> SolubilityKineticStateV1 {
    *cell.kinetic_solubility.as_ref().unwrap().last_steps()[0].state()
}

// Independent continuous coupled ODE, valid while donor crystals remain:
// L'=g(1-L)-p(L-R), R'=p(L-R), L(0)=R(0)=0, V_L=V_R=c_sat=1.
fn analytic_breakthrough(time: f64) -> [f64; 2] {
    let g: f64 = 0.1;
    let p: f64 = 0.4;
    let gap = (g * g + 4.0 * p * p).sqrt();
    let slow = (-(g + 2.0 * p) + gap) / 2.0;
    let fast = (-(g + 2.0 * p) - gap) / 2.0;
    let right = 1.0 + (fast * (slow * time).exp() - slow * (fast * time).exp()) / gap;
    let right_rate = slow * fast * ((slow * time).exp() - (fast * time).exp()) / gap;
    [right + right_rate / p, right]
}

#[test]
fn dissolution_limited_receiver_breakthrough_matches_the_coupled_analytic_solution() {
    let mut system = breakthrough(0.1);
    let report = system.advance(4.0, controls(0.5), || false).unwrap();
    let expected = analytic_breakthrough(4.0);
    for (cell, amount) in system.cells.iter().zip(expected) {
        near(cell.solution().solute_moles()[0], amount, 5e-5);
    }
    let donor = phase(&system.cells[0]);
    assert!(donor.solid().value() > 9.0);
    near(
        donor.solid().value() + donor.dissolved().value() + expected[1],
        10.0,
        5e-5,
    );
    // An instantly saturated donor would deliver 1-exp(-p*t). This much
    // slower breakthrough must not be reconstructed from that target.
    assert!(expected[1] < 0.4 * (1.0 - (-0.4_f64 * 4.0).exp()));
    near(report.solute_balance_residual_mol[0], 0.0, 1e-10);
    assert_eq!(report.water_balance_residual_m3, 0.0);
    assert!(report.phase_evaluations >= 6 * report.accepted_steps);
    near(
        system.cells[0].solution().osmotic_pressure_pa().unwrap(),
        R * T * donor.dissolved().value(),
        1e-10,
    );

    let mut blocked = breakthrough(0.0);
    blocked.advance(4.0, controls(0.5), || false).unwrap();
    assert_eq!(blocked.cells[1].solution().solute_moles()[0], 0.0);
    assert_eq!(phase(&blocked.cells[0]).solid().value(), 10.0);
}

#[test]
fn isolated_native_phase_evolution_preserves_declared_precipitation_policy() {
    for policy in [
        PrecipitationPolicyV1::ExistingSolidOnly,
        PrecipitationPolicyV1::SeededSurface,
    ] {
        let mut system = network(vec![kinetic(2.0, 0.0, 0.5, policy)], vec![]);
        let report = system.advance(3.0, controls(3.0), || false).unwrap();
        let state = phase(&system.cells[0]);
        let expected = if policy == PrecipitationPolicyV1::SeededSurface {
            1.0 + (-1.5_f64).exp()
        } else {
            2.0
        };
        near(state.dissolved().value(), expected, 2e-14);
        near(state.solid().value(), 2.0 - expected, 2e-14);
        assert_eq!(state.total_solute().value(), 2.0);
        assert_eq!(report.accepted_steps, 1);
        assert_eq!(report.phase_evaluations, 6);
    }
}

#[test]
fn water_removal_and_dilution_evolve_actual_crystals_with_finite_lag() {
    let mut donor = kinetic(1.0, 0.1, 0.2, PrecipitationPolicyV1::ExistingSolidOnly);
    donor.pressure = OsmoticPressureVolume::new(1.0, 10_000.0, 0.0).unwrap();
    let mut system = network(vec![donor, aqueous(0.0)], vec![membrane(1e-5, 1.0, 0.0)]);
    let report = system.advance(2.0, controls(0.2), || false).unwrap();
    let concentrated = phase(&system.cells[0]);
    let volume = system.cells[0].solution().solvent_volume_m3();
    assert!(volume < 0.9);
    assert!(concentrated.solid().value() > 0.1);
    assert!(concentrated.dissolved().value() > volume); // still supersaturated
    assert_eq!(system.cells[1].solution().solute_moles()[0], 0.0);
    near(concentrated.total_solute().value(), 1.1, 1e-12);
    near(report.solute_balance_residual_mol[0], 0.0, 2e-12);
    near(report.water_balance_residual_m3, 0.0, 2e-12);

    // Prescribed mechanical load change retains the complete kinetic state.
    let mut cells = system.compartments().to_vec();
    for (index, cell) in cells.iter_mut().enumerate() {
        cell.pressure = OsmoticPressureVolume::new(
            cell.solution().solvent_volume_m3(),
            if index == 1 { 20_000.0 } else { 0.0 },
            0.0,
        )
        .unwrap();
    }
    let mut diluted = network(cells, vec![membrane(1e-5, 1.0, 0.0)]);
    let report = diluted.advance(2.0, controls(0.2), || false).unwrap();
    let final_state = phase(&diluted.cells[0]);
    assert!(diluted.cells[0].solution().solvent_volume_m3() > 1.1);
    assert!(final_state.solid().value() < concentrated.solid().value());
    assert!(final_state.solid().value() > 0.0); // dissolution has not caught up
    assert!(final_state.dissolved().value() < 1.1);
    near(final_state.total_solute().value(), 1.1, 2e-12);
    near(totals(&diluted.cells, 1).0, 2.0, 3e-12);
    near(report.solute_balance_residual_mol[0], 0.0, 3e-12);
}

#[test]
fn coupled_phase_transport_converges_under_native_step_refinement() {
    let expected = analytic_breakthrough(4.0)[1];
    let errors: Vec<_> = [1.0, 0.5, 0.25]
        .into_iter()
        .map(|step| {
            let mut system = breakthrough(0.1);
            let mut options = controls(step);
            // Deliberately loose local tolerance makes the caller's maximum
            // step control this independent temporal refinement experiment.
            options.relative_tolerance = 0.2;
            options.absolute_amount_tolerance_mol = 0.01;
            let report = system.advance(4.0, options, || false).unwrap();
            assert_eq!(report.rejected_steps, 0);
            (system.cells[1].solution().solute_moles()[0] - expected).abs()
        })
        .collect();
    assert!(errors[1] < 0.6 * errors[0], "errors={errors:?}");
    assert!(errors[2] < 0.6 * errors[1], "errors={errors:?}");
}

fn particle_compartment(dissolved: f64, solid: f64, conductance: f64) -> OsmoticCompartment {
    let model = SolubilityKineticsV1::new(
        law(1.0),
        VolumetricFlowRate::new(conductance),
        PrecipitationPolicyV1::ExistingSolidOnly,
    )
    .unwrap()
    .with_particle_surface(Amount::new(solid))
    .unwrap();
    let material = KineticSolubilitySolution::new(
        Temperature::new(T),
        18e-6,
        Volume::new(1.0),
        vec![(
            model,
            SolubilityKineticStateV1::from_phases(Amount::new(dissolved), Amount::new(solid))
                .unwrap(),
        )],
    )
    .unwrap();
    OsmoticCompartment::with_kinetic_solubility(
        material,
        OsmoticPressureVolume::new(1.0, 0.0, 0.0).unwrap(),
        1e-6,
    )
    .unwrap()
}

fn particle_breakthrough() -> OsmoticNetwork {
    network(
        vec![
            particle_compartment(0.0, 1.0, 0.4),
            particle_compartment(0.05, 0.25, 0.0),
        ],
        vec![membrane(0.0, 0.0, 0.4)],
    )
}

// Independent continuous three-stock ODE for the specific particle case above.
// The receiver's quarter mole of inert crystals never enters either flux.
// Classical RK4 is used only here; it calls no production kinetic or transport
// kernel and resolves the reference well below the native refinement errors.
fn particle_breakthrough_oracle(duration: f64, steps: usize) -> [f64; 3] {
    let rate = |[left, right, solid]: [f64; 3]| {
        let dissolution = 0.4 * solid.powf(2.0 / 3.0) * (1.0 - left);
        let transport = 0.4 * (left - right);
        [dissolution - transport, transport, -dissolution]
    };
    let h = duration / steps as f64;
    let mut state = [0.0, 0.05, 1.0];
    for _ in 0..steps {
        let k1 = rate(state);
        let k2 = rate(core::array::from_fn(|i| state[i] + 0.5 * h * k1[i]));
        let k3 = rate(core::array::from_fn(|i| state[i] + 0.5 * h * k2[i]));
        let k4 = rate(core::array::from_fn(|i| state[i] + h * k3[i]));
        state = core::array::from_fn(|i| {
            state[i] + h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0
        });
    }
    state
}

#[test]
fn native_membrane_coupling_evolves_particle_area_and_transports_only_dissolved_stock() {
    let expected = particle_breakthrough_oracle(3.0, 4000);
    let coarse_reference = particle_breakthrough_oracle(3.0, 2000);
    for (fine, coarse) in expected.into_iter().zip(coarse_reference) {
        near(fine, coarse, 1e-11);
    }
    let mut system = particle_breakthrough();
    let report = system.advance(3.0, controls(0.5), || false).unwrap();
    let donor = phase(&system.cells[0]);
    let receiver = phase(&system.cells[1]);
    for (actual, expected) in [
        donor.dissolved().value(),
        receiver.dissolved().value(),
        donor.solid().value(),
    ]
    .into_iter()
    .zip(expected)
    {
        near(actual, expected, 1e-4);
    }
    assert_eq!(receiver.solid().value(), 0.25);
    let original = &system.cells[0]
        .kinetic_solubility
        .as_ref()
        .unwrap()
        .models()[0];
    assert_eq!(original.particle_reference_solid(), Some(Amount::new(1.0)));
    assert!(original.conductance_at(donor.solid()).unwrap().value() < 0.4);
    near(
        donor.dissolved().value()
            + donor.solid().value()
            + receiver.dissolved().value()
            + receiver.solid().value(),
        1.3,
        2e-11,
    );
    near(report.solute_balance_residual_mol[0], 0.0, 2e-11);
    assert_eq!(report.water_balance_residual_m3, 0.0);
    near(
        system.cells[0].solution().osmotic_pressure_pa().unwrap(),
        R * T * donor.dissolved().value(),
        1e-10,
    );
    assert!(report.phase_evaluations >= 12 * report.accepted_steps);
}

#[test]
fn native_particle_transport_refines_toward_the_independent_coupled_ode() {
    let expected = particle_breakthrough_oracle(3.0, 4000)[1];
    let errors: Vec<_> = [0.5, 0.25, 0.125]
        .into_iter()
        .map(|step| {
            let mut system = particle_breakthrough();
            let mut options = controls(step);
            options.relative_tolerance = 0.2;
            options.absolute_amount_tolerance_mol = 0.01;
            let report = system.advance(3.0, options, || false).unwrap();
            assert_eq!(report.rejected_steps, 0);
            assert_eq!(phase(&system.cells[1]).solid().value(), 0.25);
            near(report.solute_balance_residual_mol[0], 0.0, 2e-12);
            (system.cells[1].solution().solute_moles()[0] - expected).abs()
        })
        .collect();
    assert!(errors[1] < 0.65 * errors[0], "errors={errors:?}");
    assert!(errors[2] < 0.65 * errors[1], "errors={errors:?}");
}

#[test]
fn late_cancellation_and_work_limit_preserve_every_phase_and_allow_exact_retry() {
    let original = breakthrough(0.1);
    let checkpoint = original.cells[0].kinetic_solubility_snapshot().unwrap();
    let mut expected = original.clone();
    let mut calls = 0;
    let accepted = expected
        .advance(1.0, controls(0.1), || {
            calls += 1;
            false
        })
        .unwrap();
    let mut cancelled = original.clone();
    let mut observed = 0;
    assert!(matches!(
        cancelled.advance(1.0, controls(0.1), || {
            observed += 1;
            observed == calls
        }),
        Err(OsmoticNetworkError::Cancelled)
    ));
    assert_eq!(cancelled.time_s(), 0.0);
    assert_eq!(
        cancelled.cells[0].kinetic_solubility_snapshot(),
        Some(checkpoint.clone())
    );
    assert_eq!(cancelled.cells[1].solution(), original.cells[1].solution());

    let mut limited = original.clone();
    let mut limits = controls(0.01);
    limits.maximum_steps = 1;
    assert!(matches!(
        limited.advance(1.0, limits, || false),
        Err(OsmoticNetworkError::WorkBudget)
    ));
    assert_eq!(limited.time_s(), 0.0);
    assert_eq!(
        limited.cells[0].kinetic_solubility_snapshot(),
        Some(checkpoint.clone())
    );
    for mut retry in [cancelled, limited] {
        let report = retry.advance(1.0, controls(0.1), || false).unwrap();
        assert_eq!(report.phase_evaluations, accepted.phase_evaluations);
        assert_eq!(retry.time_s(), expected.time_s());
        for (actual, expected) in retry.cells.iter().zip(&expected.cells) {
            assert_eq!(actual.solution(), expected.solution());
            assert_eq!(
                actual.kinetic_solubility_snapshot(),
                expected.kinetic_solubility_snapshot()
            );
        }
    }
    assert_eq!(checkpoint.solution().solute_moles()[0], 0.0);
    assert_eq!(checkpoint.last_steps()[0].state().solid().value(), 10.0);
}

#[test]
fn transport_and_error_control_retain_solid_tails_independently() {
    let mut system = network(
        vec![
            kinetic(1.0, 1e-20, 0.0, PrecipitationPolicyV1::ExistingSolidOnly),
            aqueous(0.0),
        ],
        vec![membrane(0.0, 0.0, 0.4)],
    );
    system.advance(1.0, controls(0.1), || false).unwrap();
    assert_eq!(phase(&system.cells[0]).solid().value(), 1e-20);
    assert!(system.cells[1].solution().solute_moles()[0] > 0.0);

    // These crystal errors are tiny relative to dissolved/total inventory,
    // but exceed their own requested absolute/relative amount allowance.
    let old = kinetic(1e6, 1e-4, 0.0, PrecipitationPolicyV1::ExistingSolidOnly);
    let coarse = kinetic(1e6, 2e-4, 0.0, PrecipitationPolicyV1::ExistingSolidOnly);
    let fine = old.clone();
    assert!(scaled_error(&[old], &[coarse], &[fine], controls(1.0)) > 1.0);
}

#[test]
fn equilibrium_donor_and_kinetic_receiver_share_the_existing_membrane_solve() {
    let donor = OsmoticCompartment::with_solubility(
        FiniteSolubilitySolution::new(
            Temperature::new(T),
            18e-6,
            Volume::new(1.0),
            vec![(law(1.0), Amount::new(2.0))],
        )
        .unwrap(),
        OsmoticPressureVolume::new(1.0, 0.0, 0.0).unwrap(),
        1e-6,
    )
    .unwrap();
    let receiver = kinetic(0.2, 0.5, 0.0, PrecipitationPolicyV1::ExistingSolidOnly);
    let mut system = network(vec![donor, receiver], vec![membrane(0.0, 0.0, 0.2)]);
    let report = system.advance(1.0, controls(0.1), || false).unwrap();
    // The saturated equilibrium donor supplies R'=p*(1-R); the receiver's
    // half mole of inert crystals must never enter that transport coordinate.
    let expected = 1.0 - 0.8 * (-0.2_f64).exp();
    near(system.cells[1].solution().solute_moles()[0], expected, 5e-5);
    assert_eq!(phase(&system.cells[1]).solid().value(), 0.5);
    near(system.cells[0].solution().solute_moles()[0], 1.0, 2e-14);
    near(totals(&system.cells, 1).1[0], 2.7, 1e-12);
    near(report.solute_balance_residual_mol[0], 0.0, 1e-12);
}

#[test]
fn unsupported_phase_pathways_and_wrong_species_are_refused_at_admission() {
    use fs_material::osmosis::{diffusion::SolutionDiffusionSolve, donnan::IdealDonnan};
    use fs_thermochem::transport::maxwell_stefan::MaxwellStefanMolarSolveV1;

    let cells = [
        kinetic(0.5, 0.5, 0.1, PrecipitationPolicyV1::ExistingSolidOnly),
        aqueous(0.0),
    ];
    let material = membrane(1e-6, 1.0, 0.0).membrane_law().clone();
    let polarized = OsmoticMembrane::new_polarized(
        0,
        1,
        1.0,
        PolarizedMembraneTransport::new(material.clone(), vec![None], vec![None]).unwrap(),
    )
    .unwrap();
    assert!(
        OsmoticNetwork::new(
            vec!["solute".into()],
            cells.to_vec(),
            vec![polarized],
            1 << 20,
        )
        .is_err()
    );
    assert!(
        transient::TransientOsmoticPair::new(cells.clone(), 1.0, material, None, None, 1 << 20,)
            .is_err()
    );
    assert!(
        donnan::DonnanOsmoticPair::new(
            cells.clone(),
            IdealDonnan::new(vec![1]).unwrap(),
            [0.0; 2],
            [0.0; 2],
            1.0,
            HydraulicPermeance::new(1e-6).unwrap(),
            1 << 20,
        )
        .is_err()
    );

    let physical = diffusion::SolutionDiffusionControls {
        constitutive: SolutionDiffusionSolve {
            friction: MaxwellStefanMolarSolveV1::default(),
            entropy_tolerance_w_m3_k: 1e-9,
        },
        maximum_species_fraction: 0.2,
        amount_balance_tolerance_mol: 1e-12,
        free_energy_tolerance_j: 1e-12,
    };
    for result in [
        OsmoticNetwork::with_solution_diffusion(
            vec!["solute".into()],
            cells.to_vec(),
            vec![],
            physical,
            1 << 20,
        ),
        OsmoticNetwork::with_coupled_solution_transport(
            vec!["solute".into()],
            cells.to_vec(),
            vec![],
            vec![],
            coupled::CoupledSolutionControls {
                diffusion: physical,
                water_balance_tolerance_m3: 1e-12,
            },
            1 << 20,
        ),
    ] {
        assert!(matches!(
            result,
            Err(OsmoticNetworkError::Invalid(
                "solution diffusion does not evolve solid phase inventories or their free energy"
            ))
        ));
    }
    assert!(
        OsmoticNetwork::new(
            vec!["different-solute".into()],
            cells.to_vec(),
            vec![],
            1 << 20,
        )
        .is_err()
    );
}
