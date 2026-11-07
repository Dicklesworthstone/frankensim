//! Physical-clock kinetics on the original Whitney A–phi and electrode owners.
//! Scalar return maps and exact six-tetrahedron reductions are independent of
//! the native vector inverse, curl-curl tangent and Newton/FGMRES kernel.

use super::super::preparation;
use super::{
    SpatialHystereticInductionError, SpatialMagneticError, VectorHysteresisError, assert_mode,
    common as data, grid, near,
};
use data::{CHI, MU, PIN, conductor, field_controls, initial, resistivity, stored};
use fs_em::spatial::vector_potential::induction::hysteretic::terminals::SpatialHystereticTerminalInduction;

#[allow(dead_code)]
#[path = "../terminals/hysteretic_common.rs"]
mod electrode_reference;

fn kinetic_play(old_z: f64, h: f64, eta: f64, dt: f64) -> f64 {
    let stiffness = CHI.recip() + eta / dt;
    let drive = h + eta / dt * old_z;
    old_z
        .clamp((drive - PIN) / stiffness, (drive + PIN) / stiffness)
        .clamp(-1.0, 1.0)
}

fn endpoint(old_b: f64, old_z: f64, sigma: f64, eta: f64, dt: f64) -> (f64, f64, f64) {
    let h = data::bracket(-29.0, 29.0, |h| {
        h + 3.0 * sigma / (160.0 * dt) * (MU * h + kinetic_play(old_z, h, eta, dt) - old_b)
    });
    let z = kinetic_play(old_z, h, eta, dt);
    (h, z, MU * h + z)
}

#[test]
fn timed_whitney_history_closes_distinct_dry_viscous_ohmic_and_numerical_transfers() {
    let eta = 0.6;
    let card = data::card_with_viscosity(500.0, eta);
    let (positions, tets) = grid(1);
    let problem = data::problem_with_viscosity(&card, positions, tets);
    let mut state = initial(&problem, Some(conductor(&card, 300.0)));
    let sigma = resistivity(300.0).recip();
    let mass = 3.0 * sigma / 40.0;
    let (mut old_h, mut old_z, mut old_b) = (0.0, 0.0, 0.0);
    let mut total_viscous = 0.0;
    let history = [
        (0.01, Some(6.0)),
        (0.1, Some(6.0)),
        (1.0, Some(12.0)),
        (0.01, Some(0.0)),
        (1.0, Some(0.0)),
        (1.0, Some(-12.0)),
        (0.1, None),
        (1.0, None),
        (10.0, None),
    ];
    for (index, (dt, field)) in history.into_iter().enumerate() {
        let (h, z, b, current) = if let Some(h) = field {
            let z = kinetic_play(old_z, h, eta, dt);
            let b = MU * h + z;
            (h, z, b, 4.0 * h + mass * (b - old_b) / dt)
        } else {
            let (h, z, b) = endpoint(old_b, old_z, sigma, eta, dt);
            (h, z, b, 0.0)
        };
        let before = (
            state.potential_wb().to_vec(),
            state.cell_states().to_vec(),
            state.time_s(),
        );
        let mut polls = 0;
        let candidate = state
            .step_candidate(&[[current; 3]; 6], dt, field_controls(), &mut || {
                polls += 1;
                false
            })
            .unwrap();
        assert_mode(candidate.state(), h, z, b);
        let step = candidate.step();
        let joule = mass * (0.5 * (b - old_b)).powi(2) / dt;
        let pinning = PIN * (z - old_z).abs();
        let viscous = eta * (z - old_z).powi(2) / dt;
        let change = stored(h, z) - stored(old_h, old_z);
        let numerical = h * (b - old_b) - change - pinning - viscous;
        near(step.joule_heat_j, joule, 3e-8);
        near(step.pinning_heat_j, pinning, 3e-8);
        near(step.viscous_heat_j, viscous, 3e-8);
        near(step.hysteresis_heat_j, pinning + viscous, 3e-8);
        near(step.heat_j, joule + pinning + viscous, 3e-8);
        near(step.stored_change_j, change, 3e-8);
        near(step.numerical_dissipation_j, numerical, 3e-8);
        near(step.impressed_work_j, 0.25 * current * (b - old_b), 3e-8);
        near(step.energy_residual_j, 0.0, 3e-8);
        assert!(step.maximum_charge_residual_a <= 2e-8);
        assert!(step.maximum_ampere_residual_a <= 2e-8);
        for cell in &step.magnetic_cells {
            near(cell.pinning_heat_j, pinning / 6.0, 3e-8);
            near(cell.viscous_heat_j, viscous / 6.0, 3e-8);
            near(cell.irreversible_heat_j, (pinning + viscous) / 6.0, 3e-8);
            near(cell.energy_change_j, change / 6.0, 3e-8);
            near(cell.numerical_dissipation_j, numerical / 6.0, 3e-8);
        }
        for cell in &step.cells {
            near(cell.heat_j, joule / 6.0, 3e-8);
        }
        if index == 1 {
            // A hold at the same H continues switching; accepted memory is
            // physical kinetics rather than a rate-independent return point.
            assert!(z > 5.0 * old_z);
            let mut seen = 0;
            assert!(matches!(
                state.step_candidate(&[[current; 3]; 6], dt, field_controls(), &mut || {
                    seen += 1;
                    seen >= polls - 1
                },),
                Err(
                    SpatialHystereticInductionError::Spatial(SpatialMagneticError::Cancelled)
                        | SpatialHystereticInductionError::Material(
                            VectorHysteresisError::Cancelled
                        )
                )
            ));
        }
        assert_eq!(state.potential_wb(), before.0);
        assert_eq!(state.cell_states(), before.1);
        assert_eq!(state.time_s(), before.2);
        total_viscous += viscous;
        state = candidate.into_parts().0;
        (old_h, old_z, old_b) = (h, z, b);
    }
    assert!(total_viscous > 0.0 && old_z < -0.39);
}

#[test]
fn zero_sourced_viscosity_preserves_the_3d_native_solution_bitwise() {
    let card = data::card_with_viscosity(500.0, 0.0);
    let (positions, tets) = grid(1);
    let baseline = data::problem(&card, positions.clone(), tets.clone());
    let zero = data::problem_with_viscosity(&card, positions, tets);
    let mut a = initial(&baseline, Some(conductor(&card, 300.0)));
    let mut b = initial(&zero, Some(conductor(&card, 300.0)));
    for (dt, current) in [(0.1, 50.0), (0.2, -25.0), (0.4, 0.0), (0.01, 40.0)] {
        let ca = a
            .step_candidate(&[[current; 3]; 6], dt, field_controls(), &mut || false)
            .unwrap();
        let cb = b
            .step_candidate(&[[current; 3]; 6], dt, field_controls(), &mut || false)
            .unwrap();
        assert_eq!(ca.state().potential_wb(), cb.state().potential_wb());
        for (a, b) in ca
            .state()
            .cell_states()
            .iter()
            .zip(cb.state().cell_states())
        {
            assert_eq!(a.field_a_m(), b.field_a_m());
            assert_eq!(a.flux_density_t(), b.flux_density_t());
            assert_eq!(a.polarizations_t(), b.polarizations_t());
        }
        let a_step = ca.step();
        let b_step = cb.step();
        for (a, b) in [
            (a_step.stored_change_j, b_step.stored_change_j),
            (a_step.pinning_heat_j, b_step.pinning_heat_j),
            (a_step.hysteresis_heat_j, b_step.hysteresis_heat_j),
            (a_step.joule_heat_j, b_step.joule_heat_j),
            (a_step.heat_j, b_step.heat_j),
            (a_step.impressed_work_j, b_step.impressed_work_j),
            (
                a_step.numerical_dissipation_j,
                b_step.numerical_dissipation_j,
            ),
        ] {
            assert_eq!(a.to_bits(), b.to_bits());
        }
        assert_eq!(b_step.viscous_heat_j, 0.0);
        a = ca.into_parts().0;
        b = cb.into_parts().0;
    }
}

#[test]
fn physical_voltage_and_current_electrodes_share_the_same_intrinsic_switching_clock() {
    use electrode_reference::{Endpoint, assert_memory, electrodes, joule};

    let eta = 0.6;
    let card = data::card_with_viscosity(500.0, eta);
    let (positions, tets) = grid(1);
    let problem = data::problem_with_viscosity(&card, positions, tets);
    let original = initial(&problem, Some(conductor(&card, 300.0)));
    let mut state = SpatialHystereticTerminalInduction::new(
        original,
        &electrodes(&problem),
        preparation(),
        &mut || false,
    )
    .unwrap();
    let sigma = resistivity(300.0).recip();
    let mut old = Endpoint::default();
    for (dt, h) in [
        (0.01, 6.0),
        (0.1, 6.0),
        (1.0, 12.0),
        (0.1, 0.0),
        (0.5, -8.0),
    ] {
        let z = kinetic_play(old.z, h, eta, dt);
        let b = MU * h + z;
        // Original physical-electrode identities, with the independent kinetic
        // polarization substituted for the static play response.
        let rate = (b - old.b) / (2.0 * dt);
        let voltage = 12.0 * h / sigma + 6.0 * rate / 5.0;
        let incoming = sigma * (voltage - rate / 6.0);
        let expected = Endpoint {
            h,
            z,
            b,
            voltage,
            incoming,
        };
        let voltage_trial = state
            .step_voltage_candidate(
                &[voltage, 0.0],
                &[[0.0; 3]; 6],
                dt,
                field_controls(),
                &mut || false,
            )
            .unwrap();
        let current_owner = state
            .clone()
            .into_current_driven(preparation(), &mut || false)
            .unwrap();
        let current_trial = current_owner
            .step_current_candidate(
                &[incoming, -incoming],
                &[[0.0; 3]; 6],
                dt,
                field_controls(),
                &mut || false,
            )
            .unwrap();
        let pinning = PIN * (z - old.z).abs();
        let viscous = eta * (z - old.z).powi(2) / dt;
        let ohmic = joule(old, expected, sigma, dt);
        let change = stored(h, z) - stored(old.h, old.z);
        let numerical = h * (b - old.b) - change - pinning - viscous;
        assert_memory(voltage_trial.state().cell_states(), expected);
        assert_memory(current_trial.state().cell_states(), expected);
        for step in [voltage_trial.step(), current_trial.step()] {
            near(step.pinning_heat_j, pinning, 3e-8);
            near(step.viscous_heat_j, viscous, 3e-8);
            near(step.joule_heat_j, ohmic, 4e-7);
            near(step.heat_j, pinning + viscous + ohmic, 4e-7);
            near(step.numerical_dissipation_j, numerical, 3e-8);
            near(step.terminal_work_j, dt * voltage * incoming, 4e-7);
            near(step.impressed_work_j, 0.0, 1e-14);
            near(step.energy_residual_j, 0.0, 3e-8);
            near(step.terminals[0].incoming_current_a, incoming, 4e-7);
        }
        state = current_trial
            .into_parts()
            .0
            .into_voltage_driven(preparation(), &mut || false)
            .unwrap();
        old = expected;
    }
}
