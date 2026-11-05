//! cargo run -p fs-material --example austin_grinding_reference -- --material PACK.fsmatpk
//! Uses original experimental Table 2 fits, not fitted plant performance.
use fs_matdb::{NormalizedPack, QueryPoint, SelectionPolicy};
use fs_material::comminution::{
    austin::resolve_austin_grinding,
    kinetics::{
        UniformizationBudget,
        target::{GrindTarget, GrindTargetSettings},
    },
};
use fs_material::separation::{ParticleClass, ParticleStream};
use fs_qty::{Dimensionless, Frequency, Length, Mass, MassFlowRate, Power, QuantitySpec, Time};
use std::collections::BTreeMap;

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
    if !matches!(args.len(), 2 | 4)
        || args[0] != "--material"
        || (args.len() == 4 && args[2] != "--net-power-w")
    {
        return Err("usage: austin_grinding_reference --material PACK.fsmatpk [--net-power-w DECLARED_WATTS]".into());
    }
    let net_power = if args.len() == 4 {
        Some(Power::new(args[3].parse()?))
    } else {
        None
    };
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
        let feed = ParticleStream::new(
            model
                .daughters()
                .classes()
                .iter()
                .map(|class| {
                    ParticleClass::new(
                        class.id.clone(),
                        class.diameter,
                        MassFlowRate::new(if class.id == "size-000" { 1.0 } else { 0.0 }),
                        BTreeMap::new(),
                    )
                })
                .collect::<Result<Vec<_>, _>>()?,
            MassFlowRate::new(0.0),
        )?;
        let target = model.time_to_target(
            &feed,
            GrindTarget {
                coarsest_class: "size-001".into(),
                minimum_fraction: Dimensionless::new(0.5),
            },
            &query,
            GrindTargetSettings {
                maximum_contact_time: Time::new(240.0),
                time_tolerance: Time::new(1e-6),
                fraction_roundoff_allowance: Dimensionless::new(1e-12),
                maximum_forward_evaluations: 64,
                forward_budget: budget,
            },
            || false,
        )?;
        let bracket = target
            .bracket
            .ok_or("expected a nonzero grind-target bracket")?;
        println!(
            "  target=50_percent_in_size-001_and_finer_source_bins included={:?}",
            target.included_classes
        );
        println!(
            "    contact_time_bracket_s=[{:.9},{:.9}] achieved_fraction={:.12} nominal_fraction_allowance={:.3e} forward_calls={}",
            bracket.below.contact_time.value(),
            bracket.meeting.contact_time.value(),
            target.point.mass_fraction.value(),
            target.point.nominal_fraction_error.value(),
            target.evaluations,
        );
        if let Some(power) = net_power {
            let energy = target.energy_at_constant_net_power(power, Mass::new(0.345))?;
            println!(
                "    declared_external_net_power_W={} source_dry_batch_kg={} energy_bracket_J=[{:.6},{:.6}] meeting_specific_energy_J_per_kg={:.6}",
                energy.net_power.value(),
                energy.dry_batch_mass.value(),
                energy.below_energy.value(),
                energy.meeting_energy.value(),
                energy.specific_energy.value(),
            );
        }
    }
    Ok(())
}
