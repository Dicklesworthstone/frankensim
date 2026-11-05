//! cargo run -p fs-material --example drop_weight_reference -- --material PACK.fsmatpk
//! Five published JKDWT fits compare t10 at 1 kWh/metric tonne and the energy
//! for 20 percent t10. These are explicit impact-model calculations, not mill
//! power, throughput or measured recoveries.
use fs_matdb::{NormalizedPack, QueryPoint, SelectionPolicy};
use fs_material::comminution::{
    SpecificGrindingEnergy,
    drop_weight::{ImpactEnergyDomain, resolve_jk_drop_weight},
};
use fs_qty::{Dimensionless, Dims, QuantitySpec};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() != 2 || args[0] != "--material" {
        return Err("usage: drop_weight_reference --material PACK.fsmatpk".into());
    }
    let bytes = std::fs::read(&args[1])?;
    let pack = NormalizedPack::from_bytes(&bytes)?;
    if pack.pack_id() != "bueno-jkdwt-2021" || pack.to_bytes() != bytes {
        return Err("expected canonical bueno-jkdwt-2021 pack".into());
    }
    let domain = ImpactEnergyDomain::new([
        SpecificGrindingEnergy::new(360.0),
        SpecificGrindingEnergy::new(9000.0),
    ])?;
    println!(
        "source_pack={} source_hash={}",
        pack.pack_id(),
        pack.content_hash()
    );
    println!("scenario: t10 at 1 kWh/metric tonne; inverse energy for t10=20 percent");
    println!("declared model interval: 0.1–2.5 kWh/metric tonne; no plant calibration");
    for ore in 1..=5 {
        let query = [
            ("source-bueno-ore", f64::from(ore)),
            ("source-test-temperature-known", 0.0),
            ("source-test-pressure-known", 0.0),
            ("operating-ore-qualified", 0.0),
            ("historical-testwork-reference", 1.0),
        ]
        .into_iter()
        .try_fold(QueryPoint::new(), |point, (axis, value)| {
            point.with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
        })?;
        let card = resolve_jk_drop_weight(&pack, &query, SelectionPolicy::SingleClaimOnly, domain)?;
        let passing = card
            .response(SpecificGrindingEnergy::new(3600.0))?
            .passing_one_tenth();
        let inverse = card.energy_for_t10(Dimensionless::new(0.2))?;
        println!(
            "specimen={:?} A={:.1}% b={:.2}/(kWh/metric-tonne) calculated_t10={:.6}% energy_for_20_percent_t10={:.6} kWh/metric-tonne",
            card.calibration_source().observation().specimen,
            card.asymptotic_t10().value() * 100.0,
            card.energy_sensitivity().value() * 3600.0,
            passing.value() * 100.0,
            inverse.value() / 3600.0
        );
    }
    Ok(())
}
