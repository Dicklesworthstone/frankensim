//! Fountain equilibrium and early superfluid flow toward a warmer helium-II channel.
//! Run: cargo run -p fs-duct --features superfluid --example superfluid_superleak
//! Both channels use the same sourced 1.8 K reference, including thermal expansion.

use fs_duct::superfluid::connection::{ConnectionDrive, ConnectionKind, TwoFluidConnection};
use fs_duct::superfluid::{LinearizationLimits, TwoFluidChannel};
use fs_material::superfluid::helium::{HeliumSvpReference, SOURCE_DOI};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    const CELLS: usize = 24;
    const AREA_M2: f64 = 1e-6;
    let reference = HeliumSvpReference::at_1_8_k()?;
    let material = *reference.state();
    let delta_t = 1e-6;
    // d(mu) = dp/rho - s*dT: zero superfluid drive requires this pressure rise.
    let fountain_pressure = material.density_kg_m3() * material.entropy_j_kg_k() * delta_t;
    let limits = LinearizationLimits {
        max_pressure_perturbation_pa: 1.0,
        max_temperature_perturbation_k: 1e-4,
        max_component_speed_m_s: 0.01,
    };
    println!(
        "# helium reference DOI {SOURCE_DOI}; T0=1.8 K; saturation pressure={:.9e} Pa",
        reference.saturation_pressure_pa()
    );
    println!(
        "# delta_T={delta_t:.9e} K; predicted fountain pressure rho*s*delta_T={fountain_pressure:.9e} Pa"
    );
    println!("# two 0.05 m channels, area={AREA_M2:.9e} m2, {CELLS} cells each; closed outer ends");
    println!("# fixed reference; limits |p'|<=1 Pa, |T'|<=1e-4 K, component speed<=0.01 m/s");
    println!(
        "# ideal junction normal clamp; no inferred porous resistance, critical velocity, or background-temperature evolution"
    );
    println!(
        "# positive superfluid current is cold-to-hot; undamped early transient, not a relaxation-to-equilibrium claim"
    );
    println!(
        "case,time_s,delta_pressure_pa,delta_temperature_k,delta_mu_j_kg,superfluid_current_after_m3_s,normal_velocity_mid_m_s,energy_defect_j,mass_defect_kg,entropy_defect_j_k"
    );
    for (case, pressure_factor) in [("fountain_equilibrium", 1.0), ("thermal_step", 0.0)] {
        let make_channel = || TwoFluidChannel::try_new(material, 0.05, AREA_M2, CELLS, limits);
        let mut cold = make_channel()?;
        let mut hot = make_channel()?;
        hot.initialize(
            &[pressure_factor * fountain_pressure; CELLS],
            &[delta_t; CELLS],
            &[0.0; CELLS + 1],
            &[0.0; CELLS + 1],
            || false,
        )?;
        let mut connection = TwoFluidConnection::try_new(&cold, &hot, ConnectionKind::Superleak)?;
        let dt = 0.4 * cold.maximum_step_s().min(hot.maximum_step_s());
        // Stop before outer-end reflections return to the junction.
        for step in 1..=32 {
            let report =
                connection.step(&mut cold, &mut hot, dt, ConnectionDrive::default(), || {
                    false
                })?;
            if step == 1 || step % 4 == 0 {
                let left = cold.cell(CELLS - 1).ok_or("missing cold endpoint")?;
                let right = hot.cell(0).ok_or("missing hot endpoint")?;
                let dp = right.pressure_perturbation_pa - left.pressure_perturbation_pa;
                let dtheta = right.temperature_perturbation_k - left.temperature_perturbation_k;
                let dmu = dp / material.density_kg_m3() - material.entropy_j_kg_k() * dtheta;
                // With vn=0, the instantaneous superfluid current is -A*(vn-vs).
                let superfluid_current = -connection.flow_m3_s()[1];
                let normal = cold
                    .face(CELLS)
                    .ok_or("missing junction face")?
                    .normal_velocity_m_s;
                println!(
                    "{case},{:.9e},{dp:.9e},{dtheta:.9e},{dmu:.9e},{superfluid_current:.9e},{normal:.9e},{:.9e},{:.9e},{:.9e}",
                    connection.time_s(),
                    report.combined_energy_defect_j,
                    report.mass_balance_defect_kg,
                    report.entropy_balance_defect_j_k,
                );
            }
        }
    }
    Ok(())
}
