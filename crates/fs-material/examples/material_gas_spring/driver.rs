//! Two-way material gas / existing fs-mbd rigid-body integration.
//! Run with an explicit synthetic mixture file; see README.md.
mod source;
mod heat;

use heat::HeatDrive;

use fs_material::thermochemical_gas::energy::GasEnergySolve;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::ThermochemicalGasModel;
use fs_mbd::{Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, Vec3, Wrench};
use std::error::Error;

type Result<T> = std::result::Result<T, Box<dyn Error>>;
const AREA: f64 = 0.01; // m^2, flat piston normal to world +x
const VOLUME: f64 = 0.001; // m^3 at x=0
const EXTERIOR: f64 = 100000.0; // Pa absolute, held by a reservoir
const PISTON_MASS: f64 = 0.5; // kg
const FORCE_TOLERANCE: f64 = 1e-8; // N
const SOURCE_TOLERANCE: f64 = 1e-10; // J, separate from first-law accounting
const STEP_ENERGY_TOLERANCE: f64 = 1e-8; // absolute J, not relative to formation energy

#[derive(Clone, Debug, PartialEq)]
struct State { body: RigidBodyState, gas: ClosedGasParcel }
#[derive(Clone, Copy, Debug, Default)]
struct HeatStep {
    force_residual_n: f64,
    supplied_heat_j: f64,
    gas_work_j: f64,
    energy_residual_j: f64,
    entropy_change_j_k: f64,
    source_residual_j: f64,
}
fn properties() -> Result<MassProperties> {
    // Inertia is inactive: this centered one-axis load excites no rotation.
    Ok(MassProperties::new(PISTON_MASS, Vec3::ZERO, Vec3::new(0.001, 0.001, 0.001))?)
}
fn energy() -> GasEnergySolve {
    GasEnergySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1e-7, maximum_iterations: 128 }
}
fn initialize(model: ThermochemicalGasModel, velocity: f64) -> Result<State> {
    Ok(State {
        body: RigidBodyState::new(Pose::identity(), Vec3::new(PISTON_MASS*velocity, 0.0, 0.0), Vec3::ZERO)?,
        gas: ClosedGasParcel::new(model, source::point(400.0, EXTERIOR)?, VOLUME)?,
    })
}

/// Original no-source entry point, now using the same energy-conservative step.
fn advance(old: &State, dt: f64) -> Result<(State, f64)> {
    let (state, report) = advance_heated(old, dt, 0.0, || false)?;
    Ok((state, report.force_residual_n))
}

/// Solve body motion and the original gas law's heated midpoint TOGETHER.
/// Heat is a whole-step energy. The material owner computes pressure work and
/// the exact DISCRETE tangent; this driver never substitutes an isentropic
/// derivative, endpoint temperature average or Q/(m*cp) temperature increment.
fn advance_heated(
    old: &State, dt: f64, heat_j: f64, cancelled: impl FnMut() -> bool,
) -> Result<(State, HeatStep)> {
    advance_with_drive(old, dt, HeatDrive::PrescribedEnergy(heat_j), cancelled)
}

/// Close material-dependent heat exchange and mechanical force at one midpoint.
/// Pressure and its heat/volume derivatives stay with the material owner; only
/// the source boundary condition is new. Heat fixed-point iteration is bounded,
/// not an unconditional convergence claim for arbitrary thermal stiffness.
fn advance_with_drive(
    old: &State, dt: f64, drive: HeatDrive, mut cancelled: impl FnMut() -> bool,
) -> Result<(State, HeatStep)> {
    if !dt.is_finite() || dt <= 0.0 {
        return Err("positive finite dt required".into());
    }
    drive.validate()?;
    // A reservoir starts from an adiabatic trial, not an explicit heat impulse
    // that could cross the material domain before solving the coupled response.
    let mut heat_j = match drive {
        HeatDrive::ConductiveReservoir { .. } => 0.0,
        _ => drive.energy_j(dt, old.gas.state())?,
    };
    if cancelled() { return Err("cancelled before gas/body trial".into()); }
    let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
    let mass = properties()?;
    let gradient = [AREA];
    let mut force = old.gas.mechanical_load(old.gas.volume_m3(), &gradient, EXTERIOR)?
        .generalized_forces()?[0];
    for _ in 0..32 {
        if cancelled() { return Err("cancelled before gas/body iteration".into()); }
        let trial = integrator.step(old.body, mass,
            Wrench { force_world: Vec3::new(force, 0.0, 0.0), torque_body: Vec3::ZERO }, dt)?;
        let next_volume = VOLUME+AREA*trial.state_after.pose().position_world().x;
        let thermal = old.gas.midpoint_heat_volume(next_volume, heat_j, energy(), &mut cancelled)?;
        let residual = force-AREA*(thermal.midpoint().state().pressure-EXTERIOR);
        let target_heat = drive.energy_j(dt, thermal.midpoint())?;
        let source_residual = heat_j-target_heat;
        if residual.abs() <= FORCE_TOLERANCE && source_residual.abs() <= SOURCE_TOLERANCE {
            // Pair actual accepted mechanical work with the same pm*delta_V
            // used by gas storage. Do not hide a lagged-force defect in a large
            // baseline internal energy or publish a body without its gas state.
            let mechanical_change = trial.diagnostics_after.mechanical_energy-trial.diagnostics_before.mechanical_energy;
            let exterior_work = EXTERIOR*(next_volume-old.gas.volume_m3());
            let balance = (thermal.step().stored_energy_change_j+mechanical_change)+exterior_work-heat_j;
            if !balance.is_finite() || balance.abs() > STEP_ENERGY_TOLERANCE {
                return Err("coupled gas/body step failed its absolute joule balance".into());
            }
            let report = HeatStep { force_residual_n: residual, supplied_heat_j: heat_j,
                gas_work_j: thermal.step().work_by_gas_j, energy_residual_j: balance,
                entropy_change_j_k: thermal.step().entropy_change_j_k,
                source_residual_j: source_residual };
            if cancelled() { return Err("cancelled before gas/body publication".into()); }
            return Ok((State { body: trial.state_after, gas: thermal.into_parts().0 }, report));
        }
        // dV1/dF=AREA*dt^2/(2*m) for the ORIGINAL fs-mbd constant-force step.
        // The material derivative already includes the midpoint factor AND heat.
        let jacobian = 1.0-AREA*AREA*thermal.step().pressure_final_volume_derivative_pa_m3
            *dt*dt/(2.0*PISTON_MASS);
        if !jacobian.is_finite() || jacobian <= 0.0 || !residual.is_finite() {
            return Err("nonfinite or inadmissible coupled Newton tangent".into());
        }
        // Include pressure's response to the proposed source change; using the
        // static isentropic tangent here would solve a different coupled law.
        let heat_change = target_heat-heat_j;
        let correction = (residual-AREA*thermal.step().pressure_heat_derivative_pa_j*heat_change)/jacobian;
        if !source_residual.is_finite() || !correction.is_finite() {
            return Err("nonfinite coupled heat/force correction".into());
        }
        force -= correction;
        heat_j = target_heat;
    }
    Err("gas/rigid-body heat-source coupling exhausted its 32-iteration budget".into())
}
fn total_energy(state: &State) -> Result<f64> {
    let mechanical = RigidBodyIntegrator::new(Gravity::ZERO).diagnostics(state.body, properties()?)?.mechanical_energy;
    Ok(mechanical + state.gas.mass_kg()*state.gas.state().specific_internal_energy_j_kg()
        + EXTERIOR*state.gas.volume_m3())
}
fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.is_empty() || args.len() > 7 {
        return Err("usage: material_gas_spring MIXTURE.tsv [steps=500] [dt_s=0.0002] [heat_w=0 | wall temperature_K area_m2 distance_m]".into());
    }
    let steps: usize = args.get(1).map_or(Ok(500), |s| s.parse())?;
    let dt: f64 = args.get(2).map_or(Ok(0.0002), |s| s.parse())?;
    let drive = HeatDrive::parse(args.get(3..).unwrap_or(&[]))?;
    if !(1..=10000).contains(&steps) || !dt.is_finite() || dt <= 0.0 {
        return Err("steps must be 1..=10000 and dt must be finite and positive".into());
    }
    let text = std::fs::read_to_string(&args[0])?;
    let mut state = initialize(source::load(&text, &args[0])?, 0.2)?;
    // Admit power-to-energy conversion before printing a trajectory. Reservoir
    // conductivity is still re-evaluated at every actual coupled midpoint.
    drive.energy_j(dt, state.gas.state())?;
    let initial_energy = total_energy(&state)?;
    let mut supplied_heat = 0.0;
    let heated_output = args.get(3).is_some();
    let reservoir_output = matches!(drive, HeatDrive::ConductiveReservoir { .. });
    print!("time_s,x_m,v_m_s,temperature_k,pressure_pa,sound_speed_m_s,impedance_pa_s_m,viscosity_pa_s,conductivity_w_m_k,stiffness_n_m,total_energy_drift_j,force_residual_n");
    if heated_output { print!(",supplied_heat_j,step_gas_work_j,step_energy_residual_j,step_entropy_change_j_k"); }
    if reservoir_output { print!(",step_source_residual_j"); }
    println!();
    for i in 0..=steps {
        let report = if i == 0 { HeatStep::default() } else if heated_output {
            let (next, report) = advance_with_drive(&state, dt, drive, || false)?; state = next; report
        } else {
            let (next, residual) = advance(&state, dt)?; state = next;
            HeatStep { force_residual_n: residual, ..HeatStep::default() }
        };
        supplied_heat += report.supplied_heat_j;
        let medium = state.gas.state().state();
        let stiffness = -state.gas.mechanical_load(state.gas.volume_m3(), &[AREA], EXTERIOR)?
            .tangent_action(&[1.0], &[0.0])?[0];
        print!("{:.9},{:.12},{:.12},{:.9},{:.9},{:.9},{:.9},{:.12},{:.12},{:.9},{:.12},{:.12}",
            i as f64*dt, state.body.pose().position_world().x,
            state.body.center_of_mass_velocity_world(properties()?)?.x,
            medium.temperature, medium.pressure, medium.sound_speed, medium.characteristic_impedance,
            medium.dynamic_viscosity, medium.thermal_conductivity, stiffness,
            total_energy(&state)?-initial_energy-supplied_heat, report.force_residual_n);
        if heated_output { print!(",{:.12},{:.12},{:.12},{:.12}", supplied_heat,
            report.gas_work_j, report.energy_residual_j, report.entropy_change_j_k); }
        if reservoir_output { print!(",{:.12}", report.source_residual_j); }
        println!();
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    const A: &str = include_str!("mixture-a.tsv");
    const B: &str = include_str!("mixture-b.tsv");
    fn run(text: &str, n: usize) -> State {
        let mut state = initialize(source::load(text, "synthetic test data").unwrap(), 0.2).unwrap();
        let initial = state.clone();
        for _ in 0..n {
            let (next, residual) = advance(&state, 0.1/n as f64).unwrap();
            assert!(residual.abs() <= FORCE_TOLERANCE); state = next;
        }
        assert_eq!(state.gas.mass_kg(), initial.gas.mass_kg());
        assert_eq!(state.body.pose().position_world().y, 0.0);
        assert_eq!(state.body.angular_momentum_body(), Vec3::ZERO);
        assert!((total_energy(&state).unwrap()-total_energy(&initial).unwrap()).abs() < 1e-5);
        state
    }
    #[test]
    fn coupled_motion_converges_under_timestep_refinement() {
        let a = run(A, 100); let b = run(A, 200); let c = run(A, 400);
        let error = |x: &State, y: &State| (x.body.pose().position_world().x-y.body.pose().position_world().x).abs();
        assert!(error(&a, &b) > 3.0*error(&b, &c));
    }
    #[test]
    fn changing_only_mixture_data_changes_motion_and_sound() {
        let a = run(A, 200); let b = run(B, 200);
        assert!((a.body.pose().position_world().x-b.body.pose().position_world().x).abs() > 1e-4);
        assert!((a.gas.state().state().sound_speed-b.gas.state().state().sound_speed).abs() > 100.0);
    }
    #[test]
    fn failed_trial_does_not_modify_either_accepted_owner() {
        let state = initialize(source::load(A, "synthetic test").unwrap(), 0.2).unwrap();
        let before = state.clone();
        for dt in [0.0, -1.0, f64::NAN, f64::INFINITY, 1000.0] {
            assert!(advance(&state, dt).is_err()); assert_eq!(state, before);
        }
    }
    #[test]
    fn source_shutoff_preserves_heat_motion_and_material_history() {
        let mut state = initialize(source::load(A, "synthetic heat pulse").unwrap(), 0.2).unwrap();
        let initial = state.clone();
        let mut supplied = 0.0;
        for i in 0..600 {
            let heat = if i < 300 { 0.02 } else { 0.0 };
            let (next, report) = advance_heated(&state, 0.0001, heat, || false).unwrap();
            assert!(report.energy_residual_j.abs() <= STEP_ENERGY_TOLERANCE);
            state = next; supplied += report.supplied_heat_j;
            assert!((total_energy(&state).unwrap()-total_energy(&initial).unwrap()-supplied).abs() < 1e-7);
        }
        assert!((supplied-6.0).abs() < 1e-10);
        assert_eq!(state.gas.mass_kg(), initial.gas.mass_kg());
        assert!((state.gas.state().state().temperature-400.0).abs() > 0.1);
        assert!((state.body.pose().position_world().x).abs() > 1e-5);
    }
    fn heated_run(text: &str, count: usize) -> State {
        let mut state = initialize(source::load(text, "synthetic heated convergence").unwrap(), 0.2).unwrap();
        let dt = 0.06/count as f64;
        for _ in 0..count { state = advance_heated(&state, dt, 100.0*dt, || false).unwrap().0; }
        state
    }
    #[test]
    fn heated_trajectory_refines_and_composition_alone_changes_the_solution() {
        let a = heated_run(A, 64); let b = heated_run(A, 128); let c = heated_run(A, 256);
        let x = |s: &State| s.body.pose().position_world().x;
        assert!((x(&a)-x(&b)).abs() > 3.5*(x(&b)-x(&c)).abs());
        let other = heated_run(B, 256);
        assert!((x(&other)-x(&c)).abs() > 1e-4);
        assert!((other.gas.state().state().sound_speed-c.gas.state().state().sound_speed).abs() > 100.0);
    }
    #[test]
    fn zero_source_parity_and_late_cancellation_never_publish_one_owner_alone() {
        let old = initialize(source::load(A, "synthetic atomicity").unwrap(), 0.2).unwrap();
        assert_eq!(advance(&old, 0.0002).unwrap().0,
            advance_heated(&old, 0.0002, 0.0, || false).unwrap().0);
        let mut total = 0;
        advance_heated(&old, 0.0002, 0.02, || { total += 1; false }).unwrap();
        let saved = old.clone();
        let mut calls = 0;
        assert!(advance_heated(&old, 0.0002, 0.02, || { calls += 1; calls == total }).is_err());
        assert_eq!(old, saved);
        assert!(advance_heated(&old, 0.0002, 1e12, || false).is_err());
        assert_eq!(old, saved);
    }

    fn reservoir_run(text: &str, temperature_k: f64) -> (State, f64) {
        let mut state = initialize(source::load(text, "synthetic thermal reservoir").unwrap(), 0.0).unwrap();
        let initial = state.clone();
        let drive = HeatDrive::ConductiveReservoir { temperature_k, area_m2: 0.02, distance_m: 0.001 };
        let mut supplied = 0.0;
        for _ in 0..200 {
            let (next, report) = advance_with_drive(&state, 0.0005, drive, || false).unwrap();
            assert!(report.source_residual_j.abs() <= SOURCE_TOLERANCE);
            assert!(report.energy_residual_j.abs() <= STEP_ENERGY_TOLERANCE);
            supplied += report.supplied_heat_j;
            state = next;
            assert!((total_energy(&state).unwrap()-total_energy(&initial).unwrap()-supplied).abs() < 1e-7);
            assert_eq!(state.gas.mass_kg(), initial.gas.mass_kg());
        }
        (state, supplied)
    }
    #[test]
    fn conductivity_data_drives_heat_motion_and_sound_in_the_same_run() {
        // Only the two source conductivity columns change; same cp, viscosity,
        // species fractions, initial state, source temperature and geometry.
        let changed = A.replace("0.025", "0.050").replace("0.150", "0.300");
        let (original, original_heat) = reservoir_run(A, 500.0);
        let (doubled, doubled_heat) = reservoir_run(&changed, 500.0);
        assert!(original_heat > 0.0 && doubled_heat > original_heat);
        assert!(doubled.gas.state().state().temperature > original.gas.state().state().temperature);
        assert!(doubled.body.pose().position_world().x > original.body.pose().position_world().x);
        assert!(doubled.gas.state().state().sound_speed > original.gas.state().state().sound_speed);
        let (cooled, heat_out) = reservoir_run(A, 300.0);
        assert!(heat_out < 0.0);
        assert!(cooled.gas.state().state().temperature < 400.0);
        assert!(cooled.body.pose().position_world().x < 0.0);
    }
    #[test]
    fn isolated_source_parity_and_failed_reservoir_trials_preserve_both_owners() {
        let old = initialize(source::load(A, "synthetic thermal isolation").unwrap(), 0.2).unwrap();
        let isolated = HeatDrive::ConductiveReservoir { temperature_k: 500.0, area_m2: 0.0, distance_m: 0.001 };
        let (next, report) = advance_with_drive(&old, 0.0002, isolated, || false).unwrap();
        assert_eq!(next, advance(&old, 0.0002).unwrap().0);
        assert_eq!(report.supplied_heat_j, 0.0);
        let active = HeatDrive::ConductiveReservoir { temperature_k: 500.0, area_m2: 0.02, distance_m: 0.001 };
        let mut polls = 0;
        advance_with_drive(&old, 0.0002, active, || { polls += 1; false }).unwrap();
        let saved = old.clone();
        let mut actual = 0;
        assert!(advance_with_drive(&old, 0.0002, active, || { actual += 1; actual == polls }).is_err());
        assert_eq!(old, saved);
        let invalid = HeatDrive::ConductiveReservoir { temperature_k: 500.0, area_m2: 0.02, distance_m: 0.0 };
        assert!(advance_with_drive(&old, 0.0002, invalid, || false).is_err());
        assert_eq!(old, saved);
    }
}
