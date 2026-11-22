//! Two RC charging branches on one thermal body, assembled from Components.
//!
//! Run `cargo run -p fs-system --example electrical_switching`. The switch
//! interrupts one branch between 0.20 and 0.35 s. The other branch's supply
//! reverses at 0.70 s; its diode prevents discharge into the reversed supply.
//! At 1.00 s an increased supply resumes charging. The two resistors heat a
//! shared thermal body. All states and work come from the composition kernel;
//! this example only describes, steps and reads the assembly.
//!
//! Both switching laws are ideal. They model neither semiconductor losses nor
//! impulses, and the authored events resolve the discontinuous commands.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::PortKind;
use fs_ports::component::LoadLine;
use fs_system::elements::{
    Admittance, ElementError, Force, HeatCapacity, IdealDiode, IdealSwitch, SignalClock,
    SignalSource, SourceAccount, Storage, Waveform,
};
use fs_system::{NodeId, PortRef, SolveConfig, SystemBuilder};

const RESISTANCE_OHM: f64 = 10.0;
const CAPACITANCE_F: f64 = 0.01;

// Return the source, capacitor, ideal device, and resistor's heat outlet.
fn branch<C>(
    builder: &mut SystemBuilder,
    name: &str,
    device: C,
    voltage: Waveform,
) -> Result<(NodeId, NodeId, NodeId, PortRef), Box<dyn std::error::Error>>
where
    C: LoadLine<State = (), Error = ElementError> + 'static,
{
    let electrical = PortKind::ElectricalVoltageCurrent;
    let source = builder.add_with_tangent(
        &format!("{name}_supply"),
        Force::in_domain(electrical, voltage),
        SourceAccount::default(),
    )?;
    let resistor = builder.add_with_tangent(
        &format!("{name}_resistor"),
        Admittance::in_domain(electrical, 1.0 / RESISTANCE_OHM)?,
        (),
    )?;
    let capacitor = builder.add_with_tangent(
        &format!("{name}_capacitor"),
        Storage::in_domain(electrical, CAPACITANCE_F)?,
        0.0,
    )?;
    let device = builder.add_on_line(name, device, ())?;
    builder.common_flow(
        builder.port(resistor, "a")?,
        &[
            builder.port(source, "p")?,
            builder.port(capacitor, "p")?,
            builder.port(device, "p")?,
        ],
    )?;
    builder.boundary(
        builder.port(resistor, "b")?,
        Box::new(|_, _, out| out[0] = 0.0),
    )?;
    Ok((source, capacitor, device, builder.port(resistor, "heat")?))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut builder = SystemBuilder::new();
    // Every series member receives positive charging current. The source's
    // voltage is negative in that orientation: it supplies, rather than
    // absorbs, the electrical work.
    let (diode_source, diode_capacitor, diode, diode_heat) = branch(
        &mut builder,
        "diode",
        IdealDiode::new(),
        Waveform::held_table(vec![(0.0, -12.0), (0.70, 6.0), (1.00, -18.0)])?,
    )?;
    let (switch_source, switch_capacitor, switch, switch_heat) = branch(
        &mut builder,
        "switch",
        IdealSwitch::new(),
        Waveform::constant(-12.0)?,
    )?;
    let thermal = builder.add_with_tangent("thermal_body", HeatCapacity::new(10.0)?, 300.0)?;
    builder.common_effort(builder.port(thermal, "heat")?, &[diode_heat, switch_heat])?;
    let command = builder.add_with_tangent(
        "switch_command",
        SignalSource::new(vec![Waveform::held_table(vec![
            (0.0, 1.0),
            (0.20, 0.0),
            (0.35, 1.0),
        ])?])?,
        SignalClock::default(),
    )?;
    let diode_port = builder.port(diode, "p")?;
    let switch_port = builder.port(switch, "p")?;
    builder.signal(
        builder.port(command, "output")?,
        builder.port(switch, "closed")?,
    )?;
    let mut system = builder.seal()?;
    // Temperatures inside a loop need a valid initial interface guess,
    // separately from the thermal component's accepted storage state.
    for heat_port in [diode_heat, switch_heat] {
        assert!(system.guess(heat_port, &[300.0]));
    }
    system.set_config(SolveConfig {
        absolute_tolerance: 1.0e-12,
        ..SolveConfig::default()
    });
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| -> Result<(), Box<dyn std::error::Error>> {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey { seed: 0, kernel_id: 1, tile: 0, iteration: 0 },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        println!("time_s,diode_capacitor_voltage_v,switch_capacitor_voltage_v,diode_current_a,switch_current_a,body_temperature_k,source_work_j,energy_residual_j");
        while system.time_s() < 1.4 {
            system.step(0.005_f64.min(1.4 - system.time_s()), &cx)?;
            let diode_charge = *system.state::<f64>(diode_capacitor).ok_or("missing diode capacitor")?;
            let switch_charge = *system.state::<f64>(switch_capacitor).ok_or("missing switch capacitor")?;
            let supplied_j = -system.state::<SourceAccount>(diode_source).ok_or("missing diode supply")?.energy_in_j
                -system.state::<SourceAccount>(switch_source).ok_or("missing switch supply")?.energy_in_j;
            println!(
                "{:.8},{:.8},{:.8},{:.8},{:.8},{:.12e},{:.12e},{:.12e}",
                system.time_s(),
                diode_charge / CAPACITANCE_F,
                switch_charge / CAPACITANCE_F,
                system.last_input(diode_port)[0],
                system.last_input(switch_port)[0],
                system.state::<f64>(thermal).ok_or("missing thermal body")?,
                supplied_j,
                system.ledger().energy_residual_j(),
            );
        }
        Ok(())
    })
}
