//! cargo run -p fs-material --example nonideal_residence
//! Authored passive-constituent pulse, shutdown and flush; no plant calibration.
use fs_matdb::Provenance;
use fs_material::recovery::residence::nonideal::{
    NonidealTankCascade, ResidenceFlows, TankGeometry, ZoneInventory,
};
use fs_qty::{Density, Dimensionless, Mass, Time, Volume, VolumetricFlowRate};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let geometry = [(12.0, 3.0), (10.0, 2.0), (8.0, 2.0)]
        .into_iter()
        .map(|(a, s)| TankGeometry::new(Volume::new(a), Volume::new(s)))
        .collect::<Result<Vec<_>, _>>()?;
    let mut tanks = NonidealTankCascade::new(
        "authored passive dissolved constituent",
        geometry,
        vec![ZoneInventory { active: Mass::new(0.0), stagnant: Mass::new(0.0) }; 3],
        Provenance {
            source: "authored pulse/shutdown illustration; model families: Fogler ch18 and DOI 10.1016/j.apradiso.2011.09.003".into(),
            license: "original example; no measured plant parameter data".into(),
            artifact: None,
        },
    )?;
    let mut reports = Vec::new();
    for (name, seconds, flow, concentration) in [
        ("baseline", 60.0, 0.5, 0.2),
        ("feed pulse", 30.0, 0.5, 1.0),
        ("flow stop, mixing continues", 20.0, 0.0, 0.0),
        ("restart with lower feed", 60.0, 0.5, 0.1),
        ("clean-water flush", 120.0, 0.5, 0.0),
    ] {
        let flows = ResidenceFlows::new(
            VolumetricFlowRate::new(flow),
            Dimensionless::new(0.1),
            [0.04, 0.02, 0.03].map(VolumetricFlowRate::new).to_vec(),
        )?;
        reports.push((
            name,
            tanks.advance(Time::new(seconds), &flows, Density::new(concentration))?,
        ));
    }
    println!("Authored passive-transport example; geometry and flows are not fitted plant data.");
    println!(
        "Three unequal tanks, 30 m3 active + 7 m3 stagnant; 10% of total inlet flow bypasses the train."
    );
    println!("interval | inlet kg | outlet kg | active kg | stagnant kg | residual kg");
    let (mut total_in, mut total_out) = (0.0, 0.0);
    for (name, r) in &reports {
        total_in += r.inlet_mass.value();
        total_out += r.outlet_mass.value();
        println!(
            "{name} | {:.9} | {:.9} | {:.9} | {:.9} | {:.3e}",
            r.inlet_mass.value(),
            r.outlet_mass.value(),
            r.final_inventory.active.value(),
            r.final_inventory.stagnant.value(),
            r.conservation_residual.value()
        );
    }
    let final_inventory = tanks.total_inventory()?;
    println!(
        "total inlet={total_in:.9} kg; total outlet={total_out:.9} kg; retained={:.9} kg",
        final_inventory.active.value() + final_inventory.stagnant.value()
    );
    Ok(())
}
