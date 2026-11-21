//! A closed-loop device made entirely of reusable Components.
//!
//! Run `cargo run -p fs-system --example feedback_control` for proportional
//! feedback, or append `-- --sampled-pi` for a PI controller sampled every
//! 50 ms while the physical system steps at most 10 ms. The reference ramps
//! from +2 m/s to -1 m/s between 1 s and 1.1 s, so the actuator drives, brakes
//! and recovers work, then reverses the body. References, measurements,
//! controller memory and actuator work all belong to replayable Components.
//! CSV columns report accepted states, the held reference and the energy audit.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::PortKind;
use fs_system::SystemBuilder;
use fs_system::elements::{
    AffineSignal, ControlledSource, ControllerState, Mass, Motion, SampledController, SignalClock,
    SignalSource, SourceAccount, Waveform,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let sampled = match std::env::args().nth(1).as_deref() {
        None => false,
        Some("--sampled-pi") => true,
        Some(_) => return Err("usage: feedback_control [--sampled-pi]".into()),
    };
    let mut builder = SystemBuilder::new();
    let reference = builder.add_with_tangent(
        "reference",
        SignalSource::new(vec![Waveform::table(vec![
            (0.0, 2.0),
            (1.0, 2.0),
            (1.1, -1.0),
        ])?])?,
        SignalClock::default(),
    )?;
    let body = builder.add_with_tangent(
        "body",
        Mass::new(2.0)?.publishing(),
        Motion {
            position_m: 0.0,
            velocity_m_s: 0.0,
        },
    )?;
    let controller = if sampled {
        // z_next = z + Ts * (reference - velocity)
        // force = 4 * z + 8 * (reference - velocity), held until the next sample.
        let law = SampledController::new(
            0.05,
            1,
            2,
            1,
            vec![1.0],
            vec![0.05, -0.05],
            vec![4.0],
            vec![8.0, -8.0],
        )?;
        let initial = ControllerState::new(&law, vec![0.0])?;
        builder.add_with_tangent("controller", law, initial)?
    } else {
        // force = 8 * (reference - velocity), with gains in N s/m.
        builder.add_with_tangent(
            "controller",
            AffineSignal::new(2, 1, vec![8.0, -8.0], vec![0.0])?,
            (),
        )?
    };
    let actuator = builder.add_with_tangent(
        "actuator",
        ControlledSource::effort(PortKind::MechanicalForceVelocity),
        SourceAccount::default(),
    )?;
    let power = builder.port(actuator, "p")?;
    let reference_output = builder.port(reference, "output")?;
    builder.bond(builder.port(body, "p")?, power)?;
    builder.signal(reference_output, builder.port(controller, "input0")?)?;
    builder.signal(
        builder.port(body, "velocity")?,
        builder.port(controller, "input1")?,
    )?;
    builder.signal(
        builder.port(controller, "output")?,
        builder.port(actuator, "command")?,
    )?;
    let mut system = builder.seal()?;
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| -> Result<(), Box<dyn std::error::Error>> {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        println!("time_s,reference_at_step_start_m_s,velocity_m_s,position_m,force_n,source_energy_in_j,energy_residual_j");
        while system.time_s() < 2.0 {
            system.step(0.01_f64.min(2.0 - system.time_s()), &cx)?;
            let state = system.state::<Motion>(body).ok_or("body state missing")?;
            let source = system
                .state::<SourceAccount>(actuator)
                .ok_or("source state missing")?;
            println!(
                "{:.8},{:.8},{:.8},{:.8},{:.8},{:.12e},{:.12e}",
                system.time_s(),
                system.last_output(reference_output)[0],
                state.velocity_m_s,
                state.position_m,
                system.last_output(power)[0],
                source.energy_in_j,
                system.ledger().energy_residual_j(),
            );
        }
        Ok(())
    })
}
