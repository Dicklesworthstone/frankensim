//! cargo run -p fs-material --example austin_grinding_reference -- --material PACK.fsmatpk
//! Uses original experimental Table 2 fits, not fitted plant performance.
use fs_matdb::{NormalizedPack, QueryPoint, SelectionPolicy};
use fs_material::comminution::{austin::resolve_austin_grinding, kinetics::UniformizationBudget};
use fs_qty::{Dimensionless, Frequency, Length, Mass, QuantitySpec, Time};

fn conditions(filling: f64) -> Result<QueryPoint, fs_matdb::MatDbError> {
    [
        ("ball-charge-volume-fraction", Dimensionless::DIMS, filling),
        ("powder-charge-volume-fraction", Dimensionless::DIMS, 0.04),
        ("mill-inside-diameter", Length::DIMS, 0.204),
        ("mill-inside-length", Length::DIMS, 0.166),
        ("grinding-ball-diameter", Length::DIMS, 0.0254),
        ("rotation-frequency", Frequency::DIMS, 1.1),
        ("dry-charge-mass", Mass::DIMS, 0.345),
        ("dry-grinding", Dimensionless::DIMS, 1.0),
        ("source-test-temperature-known", Dimensionless::DIMS, 0.0),
        ("source-test-pressure-known", Dimensionless::DIMS, 0.0),
        ("fitted-reference-model", Dimensionless::DIMS, 1.0),
        ("operating-ore-qualified", Dimensionless::DIMS, 0.0),
    ]
    .into_iter()
    .try_fold(QueryPoint::new(), |point, (axis, dims, value)| {
        point.with_quantity(axis, QuantitySpec::dimensional(dims), value)
    })
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() != 2 || args[0] != "--material" {
        return Err("usage: austin_grinding_reference --material PACK.fsmatpk".into());
    }
    let bytes = std::fs::read(&args[1])?;
    let pack = NormalizedPack::from_bytes(&bytes)?;
    if pack.pack_id() != "quartz-assiros-petrakis-2017-austin" || pack.to_bytes() != bytes {
        return Err("expected canonical quartz-assiros-petrakis-2017-austin pack".into());
    }
    println!(
        "source_pack={} source_hash={}",
        pack.pack_id(),
        pack.content_hash()
    );
    println!(
        "calculated PSD from a normalized initial 0.850 mm source class; dry laboratory reference"
    );
    println!(
        "Table 2 fits retain the source-reported PSD mismatch; fines sink has no resolved sub-grid distribution"
    );
    let budget = UniformizationBudget {
        maximum_l1_tail: Dimensionless::new(1e-13),
        max_poisson_terms: 256,
        max_coefficient_products: 1_000_000,
    };
    for filling in [0.2, 0.1, 0.067] {
        let query = conditions(filling)?;
        let model = resolve_austin_grinding(
            &pack,
            &query,
            Time::new(240.0),
            SelectionPolicy::SingleClaimOnly,
        )?;
        println!("\nball_filling_percent={:.1}", filling * 100.0);
        for seconds in [30.0, 60.0, 120.0, 240.0] {
            let pass = model.prepare(Time::new(seconds), &query, budget, || false)?;
            println!("  contact_time_s={seconds:.0}");
            let mut total = 0.0;
            for class in model.daughters().classes() {
                let fraction = pass.transition().fraction(&class.id, "size-000")?.value();
                total += fraction;
                println!(
                    "    class={} source_size_um={:.0} calculated_mass_percent={:.9}",
                    class.id,
                    class.diameter.value() * 1e6,
                    fraction * 100.0
                );
            }
            if (total - 1.0).abs() > 1e-12 {
                return Err("calculated distribution fails mass closure".into());
            }
            println!(
                "    mass_residual={:.3e} nominal_L1_tail_bound={:.3e}",
                total - 1.0,
                2.0 * pass.poisson_tail_bound.value()
            );
        }
    }
    Ok(())
}
