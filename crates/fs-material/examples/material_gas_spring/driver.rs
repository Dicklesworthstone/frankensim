//! Two-way material gas / existing fs-mbd rigid-body integration.
//! Run with an explicit synthetic mixture file; see README.md.
mod source;

use fs_material::thermochemical_gas::energy::entropy::GasEntropySolve;
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

#[derive(Clone, Debug, PartialEq)]
struct State { body: RigidBodyState, gas: ClosedGasParcel }
fn properties() -> Result<MassProperties> {
    // Inertia is inactive: this centered one-axis load excites no rotation.
    Ok(MassProperties::new(PISTON_MASS, Vec3::ZERO, Vec3::new(0.001, 0.001, 0.001))?)
}
fn entropy() -> GasEntropySolve {
    GasEntropySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 1800.0,
        entropy_tolerance_j_kg_k: 1e-9, maximum_iterations: 128 }
}
fn initialize(model: ThermochemicalGasModel, velocity: f64) -> Result<State> {
    Ok(State {
        body: RigidBodyState::new(Pose::identity(), Vec3::new(PISTON_MASS*velocity, 0.0, 0.0), Vec3::ZERO)?,
        gas: ClosedGasParcel::new(model, source::point(400.0, EXTERIOR)?, VOLUME)?,
    })
}

/// Solve the actual midpoint coupling, rather than lagging the gas pressure.
/// Trial rigid-body states come from fs-mbd; trial gas states come from the
/// source-backed parcel. Neither accepted state is changed on ANY failure.
/// The only example-specific geometry is V(x)=VOLUME+AREA*x.
fn advance(old: &State, dt: f64) -> Result<(State, f64)> {
    if !dt.is_finite() || dt <= 0.0 { return Err("positive finite dt required".into()); }
    let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
    let mass = properties()?;
    let gradient = [AREA];
    let mut force = old.gas.mechanical_load(old.gas.volume_m3(), &gradient, EXTERIOR)?
        .generalized_forces()?[0];
    for _ in 0..32 {
        let trial = integrator.step(old.body, mass,
            Wrench { force_world: Vec3::new(force, 0.0, 0.0), torque_body: Vec3::ZERO }, dt)?;
        let next_volume = VOLUME+AREA*trial.state_after.pose().position_world().x;
        let middle_volume = 0.5*(old.gas.volume_m3()+next_volume);
        let mut middle = old.gas.clone();
        middle.set_volume_isentropic(middle_volume, entropy(), || false)?;
        let load = middle.mechanical_load(middle_volume, &gradient, EXTERIOR)?;
        let residual = force-load.generalized_forces()?[0];
        if residual.abs() <= FORCE_TOLERANCE {
            let mut gas = old.gas.clone();
            gas.set_volume_isentropic(next_volume, entropy(), || false)?;
            return Ok((State { body: trial.state_after, gas }, residual));
        }
        // dx_mid/dF = dt^2/(4*m) for fs-mbd's centered constant-force update.
        let force_slope = load.tangent_action(&[1.0], &[0.0])?[0];
        let jacobian = 1.0-force_slope*dt*dt/(4.0*PISTON_MASS);
        if !jacobian.is_finite() || jacobian <= 0.0 || !residual.is_finite() {
            return Err("nonfinite or inadmissible coupled Newton tangent".into());
        }
        force -= residual/jacobian;
    }
    Err("gas/rigid-body midpoint coupling exhausted its 32-iteration budget".into())
}
fn total_energy(state: &State) -> Result<f64> {
    let mechanical = RigidBodyIntegrator::new(Gravity::ZERO).diagnostics(state.body, properties()?)?.mechanical_energy;
    Ok(mechanical + state.gas.mass_kg()*state.gas.state().specific_internal_energy_j_kg()
        + EXTERIOR*state.gas.volume_m3())
}
fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.is_empty() || args.len() > 3 {
        return Err("usage: material_gas_spring MIXTURE.tsv [steps=500] [dt_s=0.0002]".into());
    }
    let steps: usize = args.get(1).map_or(Ok(500), |s| s.parse())?;
    let dt: f64 = args.get(2).map_or(Ok(0.0002), |s| s.parse())?;
    if !(1..=10000).contains(&steps) || !dt.is_finite() || dt <= 0.0 {
        return Err("steps must be 1..=10000 and dt must be finite and positive".into());
    }
    let text = std::fs::read_to_string(&args[0])?;
    let mut state = initialize(source::load(&text, &args[0])?, 0.2)?;
    let initial_energy = total_energy(&state)?;
    println!("time_s,x_m,v_m_s,temperature_k,pressure_pa,sound_speed_m_s,impedance_pa_s_m,viscosity_pa_s,conductivity_w_m_k,stiffness_n_m,total_energy_drift_j,force_residual_n");
    for i in 0..=steps {
        let residual = if i == 0 { 0.0 } else {
            let (next, residual) = advance(&state, dt)?; state = next; residual
        };
        let medium = state.gas.state().state();
        let stiffness = -state.gas.mechanical_load(state.gas.volume_m3(), &[AREA], EXTERIOR)?
            .tangent_action(&[1.0], &[0.0])?[0];
        println!("{:.9},{:.12},{:.12},{:.9},{:.9},{:.9},{:.9},{:.12},{:.12},{:.9},{:.12},{:.12}",
            i as f64*dt, state.body.pose().position_world().x,
            state.body.center_of_mass_velocity_world(properties()?)?.x,
            medium.temperature, medium.pressure, medium.sound_speed, medium.characteristic_impedance,
            medium.dynamic_viscosity, medium.thermal_conductivity, stiffness,
            total_energy(&state)?-initial_energy, residual);
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
}
