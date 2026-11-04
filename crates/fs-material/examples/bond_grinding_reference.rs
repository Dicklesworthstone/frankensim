//! cargo run -p fs-material --example bond_grinding_reference -- --material PACK.fsmatpk
//!
//! Five historical specimen-specific Bond tests drive the existing grinding law.
//! The example explicitly evaluates only each reported test-size pair. The
//! 1000 kW net-pinion power and 90% utilization are an authored scenario, not
//! measured plant inputs. Energy and dry throughput are empirical estimates.

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy};
use fs_material::comminution::{
    BondSizeDomain, WorkIndexMassBasis, source::resolve_bond_ball_mill,
};
use fs_qty::{Dimensionless, Dims, Length, Power, QuantitySpec};
use std::error::Error;

const PACK: &str = "joanna-hosco-sgs-bond-2011";
const USAGE: &str = "usage: bond_grinding_reference --material PACK.fsmatpk";

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.len() != 2 || args[0] != "--material" {
        return Err(USAGE.into());
    }
    let bytes = std::fs::read(&args[1])?;
    let pack = NormalizedPack::from_bytes(&bytes)?;
    if pack.pack_id() != PACK || pack.to_bytes() != bytes {
        return Err(format!("expected canonical material pack {PACK}").into());
    }
    println!(
        "source_pack={} source_hash={}",
        pack.pack_id(),
        pack.content_hash()
    );
    println!("authored scenario: net pinion power=1000 kW; utilization=0.9");
    println!("model size domain: each reported test point only; no plant calibration");
    for row in 1..=5 {
        let query = [
            ("source-joanna-composite", f64::from(row)),
            ("source-test-temperature-known", 0.0),
            ("source-test-pressure-known", 0.0),
            ("source-preparation-detail-known", 0.0),
            ("operating-ore-qualified", 0.0),
            ("historical-testwork-reference", 1.0),
        ]
        .into_iter()
        .try_fold(QueryPoint::new(), |point, (axis, value)| {
            point.with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
        })?;
        let size = |key| -> Result<Length, Box<dyn Error>> {
            let answer = pack.claims().query_typed(
                &PropertyKey::new(key, Length::DIMS),
                &query,
                SelectionPolicy::SingleClaimOnly,
            )?;
            Ok(Length::new(answer.evidence.value.value))
        };
        let feed = size("bond-test-feed-f80")?;
        let product = size("bond-test-product-p80")?;
        let card = resolve_bond_ball_mill(
            &pack,
            &query,
            WorkIndexMassBasis::MetricTonne,
            SelectionPolicy::SingleClaimOnly,
            BondSizeDomain::new([feed; 2], [product; 2])?,
        )?;
        let energy = card.specific_energy(feed, product)?;
        let flow = card.dry_throughput(
            feed,
            product,
            Power::new(1_000_000.0),
            Dimensionless::new(0.9),
        )?;
        println!(
            "specimen={:?} Wi={:.1} kWh/metric-tonne F80={:.0} um P80={:.0} um screen={:.0} um estimated_energy={:.6} kWh/metric-tonne estimated_dry_flow={:.6} tonne/h",
            card.test_context().specimen,
            card.work_index().source_value(),
            feed.value() * 1e6,
            product.value() * 1e6,
            card.test_context().closing_screen.value() * 1e6,
            energy.value() / 3600.0,
            flow.value() * 3.6,
        );
        let source = card.calibration_source().ok_or("missing retained source")?;
        for answer in source.answers() {
            println!(
                "  source_property={} selected_claim={} observation_backed={}",
                answer.receipt.property,
                answer.receipt.selected.0,
                answer.receipt.observation_backed,
            );
        }
    }
    Ok(())
}
