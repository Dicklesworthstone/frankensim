//! Native open Euler steps against prescribed source-compatible exterior gas.
//! G0/G1: independent polynomial-caloric SSP and entropy boundary contractions.
//! G4: unsupported scopes and cancellation publish neither state nor clock.

use super::*;
use crate::multicomponent::tests::{gas, source};
use fs_matdb::LawParameter;
use fs_thermochem::{
    Nasa9RegionV1, Nasa9StandardStateModelV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
    required_nasa9_parameter_dims,
};

// [Cp/R constant, Cp/R linear coefficient, M/R], with source pressure p0=1.
const ROWS: [[f64; 3]; 2] = [[3.5, 0.2, 1.0], [4.0, 0.4, 2.0]];
const T0: f64 = 0.001;
const AREA: f64 = 0.7;

fn controls() -> StepControl {
    StepControl {
        reconstruction: Reconstruction::FirstOrder,
        ..StepControl::default()
    }
}

fn model(k: usize, formation: f64, entropy_constant: f64) -> Nasa9StandardStateModelV1 {
    let [a, b, mass] = ROWS[k];
    let original = source(if k == 0 { "A" } else { "B" }, mass * R, a, b, formation);
    let mut card = original.regions()[0].card().clone();
    let (name, dims) = required_nasa9_parameter_dims().into_iter().nth(8).unwrap();
    card.parameters.insert(name.to_owned(), LawParameter { value: entropy_constant, dims });
    Nasa9StandardStateModelV1::new(
        original.species().clone(),
        original.molar_mass(),
        original.convention().clone(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    ).unwrap()
}

fn field(
    gas: MixtureGas,
    boundaries: Boundaries,
    initial: impl Fn(usize) -> (Vec<f64>, f64, f64),
) -> MixtureEuler1d {
    let states = (0..4).map(|i| {
        let (rho, velocity, temperature) = initial(i);
        let concentration: f64 = rho.iter().zip(gas.source().models())
            .map(|(density, model)| density / model.molar_mass().value()).sum();
        gas.conserved(&rho, velocity, R * concentration * temperature).unwrap()
    }).collect();
    MixtureEuler1d::new(gas, Length::new(0.25), AREA, states, boundaries, || false).unwrap()
}

fn exterior(gas: &MixtureGas, rho: &[f64], velocity: f64, temperature: f64) -> Conserved {
    let concentration: f64 = rho.iter().zip(gas.source().models())
        .map(|(density, model)| density / model.molar_mass().value()).sum();
    gas.conserved(rho, velocity, R * concentration * temperature).unwrap()
}

fn fixture(formation: f64, constants: [f64; 2]) -> (MixtureEuler1d, [Conserved; 2]) {
    let gas = gas(vec![model(0, formation, constants[0]), model(1, -2.0 * formation, constants[1])]);
    let exteriors = [
        exterior(&gas, &[1.1, 0.15], 0.35, 3.0),
        exterior(&gas, &[0.25, 0.75], -0.2, 1.2),
    ];
    let flow = field(
        gas, Boundaries::Ends { left: Boundary::Transmissive, right: Boundary::Transmissive },
        |i| {
            let x = (i as f64 + 0.5) / 4.0;
            (vec![0.65 + 0.1 * x, 0.35 - 0.1 * x], 0.15 - 0.2 * x, 1.8 + 0.4 * x)
        },
    );
    (flow, exteriors)
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 3e-11 * expected.abs().max(1.0),
        "actual={actual:e}, expected={expected:e}",
    );
}

// Independently invert the fixture's quadratic caloric energy, then integrate
// its ideal partial-pressure entropy. No production entropy reducer is used.
fn source_entropy(flow: &MixtureEuler1d, constants: [f64; 2]) -> f64 {
    (0..flow.cell_count()).map(|i| {
        let q = flow.cell(i).unwrap();
        let rho: f64 = q.partial_densities.iter().sum();
        let a: f64 = ROWS.iter().enumerate().map(|(k, row)| {
            q.partial_densities[k] * (row[0] - 1.0) / row[2]
        }).sum();
        let b: f64 = ROWS.iter().enumerate().map(|(k, row)| {
            0.5 * q.partial_densities[k] * row[1] / row[2]
        }).sum();
        let target = q.energy - 0.5 * q.momentum * q.momentum / rho + a * T0 + b * T0 * T0;
        let t = 2.0 * target / (a + (a * a + 4.0 * b * target).sqrt());
        flow.mesh().volumes()[i] * ROWS.iter().enumerate().map(|(k, row)| {
            let [a, b, mass] = *row;
            let n = q.partial_densities[k] / mass;
            if n == 0.0 { 0.0 } else {
                n * ((a - 1.0) * t.ln() + b * t + constants[k] - n.ln())
            }
        }).sum::<f64>()
    }).sum()
}

// For eta=-S, v_i=e_i,ref/T + R/M_i - s_i - u²/(2T),
// v_m=u/T, v_E=-1/T, psi=u*sum_i rho_i R/M_i.
// Actual inward extensive transfers carry their sign already; only psi needs
// the explicit end orientation. This detects omitted kinetic, momentum or
// entropy-potential terms, including when the exterior has nonzero velocity.
#[allow(clippy::too_many_arguments)]
fn boundary_input(
    rho: [f64; 2],
    velocity: f64,
    temperature: f64,
    constants: [f64; 2],
    transfer: &Totals,
    sign: f64,
    dt: f64,
) -> f64 {
    let mut input = transfer.energy / temperature - velocity * transfer.momentum / temperature;
    let mut concentration = 0.0;
    for (k, [a, b, mass]) in ROWS.into_iter().enumerate() {
        let internal = ((a - 1.0) * (temperature - T0)
            + 0.5 * b * (temperature * temperature - T0 * T0)) / mass;
        let entropy = ((a - 1.0) * temperature.ln() + b * temperature
            + constants[k] - (rho[k] / mass).ln()) / mass;
        let effort = internal / temperature + 1.0 / mass - entropy
            - 0.5 * velocity * velocity / temperature;
        input -= effort * transfer.species_masses[k];
        concentration += rho[k] / mass;
    }
    input + sign * AREA * dt * velocity * concentration
}

#[test]
fn g0_uniform_moving_exterior_carries_physical_entropy_through_both_ends() {
    let gas = gas(vec![model(0, 0.0, 0.0), model(1, 0.0, 0.0)]);
    let outside = exterior(&gas, &[0.8, 0.2], 0.7, 2.0);
    let mut flow = field(
        gas, Boundaries::Ends { left: Boundary::Transmissive, right: Boundary::Transmissive },
        |_| (vec![0.8, 0.2], 0.7, 2.0),
    );
    let original = flow.states.clone();
    let original_entropy = source_entropy(&flow, [0.0; 2]);
    let (step, entropy) = flow.step_entropy_consistent_with_exterior_states(
        Time::new(0.01), controls(), MixtureEntropyControl::default(),
        [Some(outside.as_cell()), Some(outside.as_cell())], || false,
    ).unwrap();
    assert_eq!(step.elapsed.value(), 0.01);
    assert_eq!(step.halvings, 0);
    assert_eq!(flow.states, original);
    // Total volume is AREA * 1; an unchanged uniform flow transports S*u*A*dt.
    let expected = original_entropy * 0.7 * step.elapsed.value();
    assert!(expected > 0.0);
    near(entropy.boundary_entropy_input_j_k[0], expected);
    near(entropy.boundary_entropy_input_j_k[1], -expected);
    near(entropy.boundary_entropy_input_j_k[0], boundary_input(
        [0.8, 0.2], 0.7, 2.0, [0.0; 2], &step.left_boundary_change, 1.0, 0.01,
    ));
    near(entropy.source_entropy_change_j_k, 0.0);
    near(entropy.entropy_generation_j_k, 0.0);
    assert!(step.left_boundary_change.mass() > 0.0);
    assert!(step.right_boundary_change.mass() < 0.0);
}

#[test]
fn g1_variable_caloric_open_ssp_matches_independent_states_transfers_and_source_entropy() {
    // Independent polynomial EOS, logarithmic means and explicit SSP-RK2
    // evaluate the four-cell fixture at dt=.002. These values are not produced
    // by calling either production entropy face or the production time step.
    let expected = [
        [0.6677052758517895, 0.33582830214278103, 0.13355357261421108, 4.414355023262033],
        [0.6877880487350695, 0.3126405437930241, 0.07426977324983408, 4.649147852934464],
        [0.7127666162018477, 0.28762422207225347, 0.02415253069217061, 4.9565495790627425],
        [0.7346565232058466, 0.26616271805824265, -0.020644268983220596, 5.253713458215744],
    ];
    let expected_transfers = [
        [0.0009827316368473777, -0.00021996034525055547, 0.003806029950416452, 0.012100129587002606],
        [-0.0004723504378005628, 0.0006147229068532836, -0.0018229986251422877, -0.0026822747288805663],
    ];
    let mut reference = None;
    for (formation, constants) in [(0.0, [0.0; 2]), (1e18, [0.0; 2]), (1e18, [0.8, -0.4])] {
        let (mut flow, outside) = fixture(formation, constants);
        let before = source_entropy(&flow, constants);
        let (step, entropy) = flow.step_entropy_consistent_with_exterior_states(
            Time::new(0.002), controls(), MixtureEntropyControl::default(),
            [Some(outside[0].as_cell()), Some(outside[1].as_cell())], || false,
        ).unwrap();
        assert_eq!(step.elapsed.value(), 0.002);
        assert_eq!(step.halvings, 0);
        assert_eq!(step.hll_fallback_faces, 0);
        for (i, values) in expected.iter().enumerate() {
            let q = flow.cell(i).unwrap();
            for (actual, expected) in [
                q.partial_densities[0], q.partial_densities[1], q.momentum, q.energy,
            ].into_iter().zip(values) { near(actual, *expected); }
        }
        for (transfer, values) in [step.left_boundary_change, step.right_boundary_change]
            .into_iter().zip(expected_transfers) {
            for (actual, expected) in [
                transfer.species_masses[0], transfer.species_masses[1], transfer.momentum, transfer.energy,
            ].into_iter().zip(values) { near(actual, expected); }
        }
        let measured = source_entropy(&flow, constants) - before;
        let input = [
            boundary_input([1.1, 0.15], 0.35, 3.0, constants, &step.left_boundary_change, 1.0, 0.002),
            boundary_input([0.25, 0.75], -0.2, 1.2, constants, &step.right_boundary_change, -1.0, 0.002),
        ];
        for end in 0..2 { near(entropy.boundary_entropy_input_j_k[end], input[end]); }
        near(entropy.source_entropy_change_j_k, measured);
        near(entropy.source_entropy_generation_j_k, measured - input.iter().sum::<f64>());
        near(entropy.entropy_generation_j_k, entropy.source_entropy_generation_j_k
            - entropy.reference_entropy_residual_j_k);
        near(entropy.entropy_generation_j_k, 0.003542209728953078);
        near(entropy.minimum_step_entropy_generation_j_k, entropy.entropy_generation_j_k);
        assert!(entropy.maximum_relative_species_residual < 1e-12);
        assert!(entropy.entropy_generation_j_k > 0.0);
        for k in 0..2 {
            near(step.after.species_masses[k] - step.before.species_masses[k],
                step.boundary_change.species_masses[k]);
        }
        near(step.conservation_residual.energy, 0.0);
        near(step.conservation_residual.momentum, 0.0);
        if let Some((states, generation)) = &reference {
            assert_eq!(&flow.states, states);
            near(entropy.entropy_generation_j_k, *generation);
        } else {
            reference = Some((flow.states.clone(), entropy.entropy_generation_j_k));
            near(measured, 0.004589286221221567);
            near(input[0], 0.0033686468512359284);
            near(input[1], -0.0023215703589674398);
        }
    }
}

#[test]
fn g1_one_open_end_and_reflecting_wall_preserve_common_absent_species() {
    let gas = gas(vec![
        model(0, 0.0, 0.0), model(1, 0.0, 0.0), source("C", R, 5.0, 0.0, 0.0),
    ]);
    let outside = exterior(&gas, &[1.0, 0.2, 0.0], 0.15, 2.5);
    let mut flow = field(
        gas, Boundaries::Ends { left: Boundary::Transmissive, right: Boundary::Reflecting },
        |_| (vec![0.8, 0.2, 0.0], 0.0, 2.0),
    );
    let before = source_entropy(&flow, [0.0; 2]);
    let (step, entropy) = flow.step_entropy_consistent_with_exterior_states(
        Time::new(0.002), controls(), MixtureEntropyControl::default(),
        [Some(outside.as_cell()), None], || false,
    ).unwrap();
    assert_eq!(entropy.boundary_entropy_input_j_k[1], 0.0);
    assert_eq!(step.right_boundary_change.mass(), 0.0);
    assert_eq!(step.right_boundary_change.energy, 0.0);
    assert_ne!(step.right_boundary_change.momentum, 0.0);
    assert_eq!(step.boundary_change.species_masses[2], 0.0);
    for i in 0..flow.cell_count() { assert_eq!(flow.cell(i).unwrap().partial_densities[2], 0.0); }
    let input = boundary_input(
        [1.0, 0.2], 0.15, 2.5, [0.0; 2], &step.left_boundary_change, 1.0, step.elapsed.value(),
    );
    near(entropy.boundary_entropy_input_j_k[0], input);
    near(entropy.source_entropy_generation_j_k, source_entropy(&flow, [0.0; 2]) - before - input);
    assert!(entropy.entropy_generation_j_k >= -1e-10);
}

#[test]
fn g0_absent_exteriors_preserve_native_closed_step_bit_for_bit() {
    let (mut original, _) = fixture(0.0, [0.0; 2]);
    let (mut through_open_api, _) = fixture(0.0, [0.0; 2]);
    original.set_boundaries(Boundaries::Periodic, || false).unwrap();
    through_open_api.set_boundaries(Boundaries::Periodic, || false).unwrap();
    let (expected, closed) = original.step_entropy_consistent(
        Time::new(0.002), controls(), MixtureEntropyControl::default(), || false,
    ).unwrap();
    let (actual, open) = through_open_api.step_entropy_consistent_with_exterior_states(
        Time::new(0.002), controls(), MixtureEntropyControl::default(), [None; 2], || false,
    ).unwrap();
    assert_eq!(actual, expected);
    assert_eq!(through_open_api.states, original.states);
    assert_eq!(through_open_api.time(), original.time());
    assert_eq!(open.boundary_entropy_input_j_k, [0.0; 2]);
    assert_eq!(open.entropy_generation_j_k, closed.entropy_change_j_k);
    assert_eq!(open.source_entropy_change_j_k, closed.source_entropy_change_j_k);
    assert_eq!(open.minimum_step_entropy_generation_j_k, closed.minimum_step_entropy_change_j_k);
    assert_eq!(open.maximum_relative_species_residual, closed.maximum_relative_species_residual);
}

#[test]
fn g4_missing_exterior_species_support_dimensions_and_unsupported_scope_refuse_atomically() {
    let (mut flow, outside) = fixture(0.0, [0.0; 2]);
    let initial = flow.states.clone();
    let wrong_axis = Conserved { partial_densities: vec![1.0], momentum: 0.0, energy: 2.0 };
    let one_sided_zero = exterior(flow.gas(), &[1.0, 0.0], 0.0, 2.0);
    // A Cell carries the caller's canonical source axis, not a second source
    // owner. Dimension mismatch is detectable; equal-size foreign axes are not.
    for ends in [
        [None, Some(outside[1].as_cell())],
        [Some(wrong_axis.as_cell()), Some(outside[1].as_cell())],
        [Some(one_sided_zero.as_cell()), Some(outside[1].as_cell())],
    ] {
        assert!(flow.step_entropy_consistent_with_exterior_states(
            Time::new(0.002), controls(), MixtureEntropyControl::default(), ends, || false,
        ).is_err());
        assert_eq!(flow.states, initial);
        assert_eq!(flow.time().value(), 0.0);
    }
    let mut muscl = controls();
    muscl.reconstruction = Reconstruction::Muscl;
    assert!(flow.step_entropy_consistent_with_exterior_states(
        Time::new(0.002), muscl, MixtureEntropyControl::default(),
        [Some(outside[0].as_cell()), Some(outside[1].as_cell())], || false,
    ).is_err());
    flow.set_boundaries(Boundaries::Periodic, || false).unwrap();
    assert!(flow.step_entropy_consistent_with_exterior_states(
        Time::new(0.002), controls(), MixtureEntropyControl::default(),
        [Some(outside[0].as_cell()), Some(outside[1].as_cell())], || false,
    ).is_err());
    assert_eq!(flow.states, initial);
    assert_eq!(flow.time().value(), 0.0);
}

#[test]
fn g4_open_step_late_cancellation_rolls_back_and_the_owner_remains_reusable() {
    let (mut reference, outside) = fixture(0.0, [0.0; 2]);
    let mut polls = 0;
    let (_, expected) = reference.step_entropy_consistent_with_exterior_states(
        Time::new(0.002), controls(), MixtureEntropyControl::default(),
        [Some(outside[0].as_cell()), Some(outside[1].as_cell())],
        || { polls += 1; false },
    ).unwrap();
    let (mut flow, outside) = fixture(0.0, [0.0; 2]);
    let initial = flow.states.clone();
    let mut visits = 0;
    let error = flow.step_entropy_consistent_with_exterior_states(
        Time::new(0.002), controls(), MixtureEntropyControl::default(),
        [Some(outside[0].as_cell()), Some(outside[1].as_cell())],
        || { visits += 1; visits >= polls },
    ).unwrap_err();
    assert_eq!(error, Error::Cancelled);
    assert_eq!(flow.states, initial);
    assert_eq!(flow.time().value(), 0.0);
    let (_, actual) = flow.step_entropy_consistent_with_exterior_states(
        Time::new(0.002), controls(), MixtureEntropyControl::default(),
        [Some(outside[0].as_cell()), Some(outside[1].as_cell())], || false,
    ).unwrap();
    assert_eq!(flow.time(), reference.time());
    for (a, b) in flow.states.iter().zip(&reference.states) { near(*a, *b); }
    near(actual.entropy_generation_j_k, expected.entropy_generation_j_k);
    for end in 0..2 { near(actual.boundary_entropy_input_j_k[end], expected.boundary_entropy_input_j_k[end]); }
}
