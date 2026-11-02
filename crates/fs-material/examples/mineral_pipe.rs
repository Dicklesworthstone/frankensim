//! Signed slurry pipe flow with explicitly authored Herschel–Bulkley parameters.
//! Run `cargo run -p fs-material --example mineral_pipe -- --help` for usage.

use fs_matdb::Provenance;
use fs_material::slurry::rheology::{
    HerschelBulkleyCard, HerschelBulkleyParameters, LaminarCircularPipe, SlurryConditions,
    SlurryDomain, SlurryTestContext,
};
use fs_qty::{Dimensionless, Frequency, Length, Stress, Temperature, VolumetricFlowRate};

const HELP: &str = "Usage: mineral_pipe [volume_flow_m3_s] [solids_mass_fraction]

Defaults: volume flow = 0.005 m³/s; dry-solids / total wet mass = 0.6.
Negative volume flow reverses the pressure-drop direction. Zero flow selects
zero pressure from the static yield interval. Solids mass fraction must be 0.4–0.7.

Authored example: inner radius 0.1 m, length 100 m, temperature 298.15 K,
yield stress 25 Pa, excess stress 8 Pa at 10/s, flow index 0.6.
The declared source domain is 0–500/s, 293–303 K and 0.4–0.7 solids mass fraction.
The caller separately declares a laminar wall-rate ceiling of 200/s.

These parameters are arithmetic scenario inputs, not measurements of mineral ore.
The calculation assumes homogeneous, fully developed, no-slip laminar flow.
It does not predict settling or transition to turbulence, or include elevation,
fittings, entrance losses or pump efficiency. Unsupported requests are refused.";

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.iter().any(|arg| arg == "--help" || arg == "-h") {
        println!("{HELP}");
        return Ok(());
    }
    if args.len() > 2 {
        return Err("expected at most two positional arguments; use --help".into());
    }
    let requested_flow: f64 = args.first().map_or(Ok(0.005), |value| value.parse())?;
    let solids_mass_fraction: f64 = args.get(1).map_or(Ok(0.6), |value| value.parse())?;
    let conditions = SlurryConditions {
        temperature: Temperature::new(298.15),
        solids_mass_fraction: Dimensionless::new(solids_mass_fraction),
    };
    let card = HerschelBulkleyCard::new(
        HerschelBulkleyParameters {
            yield_stress: Stress::new(25.0),
            reference_excess_stress: Stress::new(8.0),
            reference_shear_rate: Frequency::new(10.0),
            flow_index: Dimensionless::new(0.6),
        },
        SlurryDomain::new(
            [Frequency::new(0.0), Frequency::new(500.0)],
            [Temperature::new(293.0), Temperature::new(303.0)],
            [Dimensionless::new(0.4), Dimensionless::new(0.7)],
        )?,
        SlurryTestContext {
            specimen: "authored mineral-pipe scenario; no measured ore specimen".into(),
            method: "authored constant-parameter fit for an arithmetic example".into(),
            provenance: Provenance {
                source: "authored example parameters; law: NIST SP 946 (2001), section 6; \
                         pipe integral: Ashena et al. (2023), doi:10.1007/s13202-022-01566-2"
                    .into(),
                license: "project-authored example parameters".into(),
                artifact: None,
            },
        },
    )?;
    let pipe =
        LaminarCircularPipe::new(Length::new(0.1), Length::new(100.0), Frequency::new(200.0))?;
    let solution = pipe.pressure_for_flow(
        &card,
        conditions,
        VolumetricFlowRate::new(requested_flow),
        VolumetricFlowRate::new(1e-12),
        96,
    )?;
    let flow = solution.flow;

    println!("Authored mineral-pipe scenario; no measured ore or plant calibration.");
    println!("Homogeneous, fully developed, no-slip laminar flow; radius 0.1 m, length 100 m.");
    println!(
        "Temperature:                 {:.2} K",
        conditions.temperature.value()
    );
    println!("Dry-solids mass fraction:    {solids_mass_fraction:.6}");
    println!("Requested signed flow:      {requested_flow:.9} m³/s");
    println!(
        "Achieved signed flow:       {:.9} m³/s",
        flow.volumetric_flow.value()
    );
    println!(
        "Frictional pressure drop:   {:.6} Pa",
        flow.pressure_drop.value()
    );
    println!(
        "Signed wall shear stress:   {:.6} Pa",
        flow.wall_stress.value()
    );
    println!(
        "Signed wall shear rate:     {:.6} 1/s",
        flow.wall_shear_rate.value()
    );
    println!(
        "Unsheared core radius:      {:.9} m",
        flow.plug_radius.value()
    );
    println!(
        "Hydraulic dissipated power: {:.6} W",
        flow.dissipated_power.value()
    );
    println!(
        "Achieved flow residual:     {:.6e} m³/s",
        solution.flow_residual.value()
    );
    println!("Inverse iterations:         {}", solution.iterations);
    println!("Laminar wall-rate ceiling is caller-declared: 200/s.");
    println!("No settling/turbulence prediction; elevation, fittings and pump losses excluded.");
    Ok(())
}
