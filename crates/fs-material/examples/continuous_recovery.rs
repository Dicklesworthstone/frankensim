//! Separately sourced continuous treatment laws, with actual hydraulic flow.
//!
//! `cargo run -p fs-material --example continuous_recovery -- [ACTUAL_FLOW_M3_S]`
//!
//! Default flow, fixed active volume, dry throughput and grade are authored
//! example inputs. A flow argument can come from `SlurryStream::state()`'s
//! `total_volume_flow`; the example does not infer or query density itself.

use std::{collections::BTreeMap, error::Error};

use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance, QueryPoint,
};
use fs_material::recovery::{
    MAXIMUM_RATE_PARAMETER, RATE_PARAMETER, ULTIMATE_RECOVERY_PARAMETER,
    continuous::{
        CONTINUOUS_SINGLE_RATE_LAW, CONTINUOUS_UNIFORM_RATE_LAW, ContinuousRecoveryCard,
        MEAN_RESIDENCE_TIME_AXIS, TANK_COUNT_AXIS,
    },
};
use fs_qty::{Dims, Frequency, MassFlowRate, QuantitySpec, Time, Volume, VolumetricFlowRate};

const HOUR: f64 = 3600.0;
const NOMINAL_FLOW: f64 = 0.19535657;
const USAGE: &str = "usage: continuous_recovery [ACTUAL_FLOW_M3_S, default 0.19535657]";

fn synthetic_card(uniform: bool) -> ConstitutiveModelCard {
    let (law, parameter, rate) = if uniform {
        (
            CONTINUOUS_UNIFORM_RATE_LAW,
            MAXIMUM_RATE_PARAMETER,
            2.0 / (6.0 * HOUR),
        )
    } else {
        (
            CONTINUOUS_SINGLE_RATE_LAW,
            RATE_PARAMETER,
            1.0 / (6.0 * HOUR),
        )
    };
    let declaration = format!(
        "Synthetic stationary local first-order law; no measured calibration. \
         law={law}; version=1; Rinf=0.9 fraction; {parameter}={rate:.17e} s^-1; \
         equal ideal full tanks; constant carrier flow and local rates; no hidden state. \
         Hydraulic mean residence domain=[43200,172800] s; tank count=[1,16]; pH=[9,11]. \
         This is an independently authored continuous law, not an extension of a finite batch test."
    );
    let artifact = hash_domain(
        "continuous-recovery-authored-local-law-v1",
        declaration.as_bytes(),
    );
    ConstitutiveModelCard {
        law: LawId(law.into()),
        law_version: 1,
        parameters: BTreeMap::from([
            (
                ULTIMATE_RECOVERY_PARAMETER.into(),
                LawParameter {
                    value: 0.9,
                    dims: Dims::NONE,
                },
            ),
            (
                parameter.into(),
                LawParameter {
                    value: rate,
                    dims: Frequency::DIMS,
                },
            ),
        ]),
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                MEAN_RESIDENCE_TIME_AXIS,
                QuantitySpec::dimensional(Time::DIMS),
                12.0 * HOUR,
                48.0 * HOUR,
            )
            .with_quantity(
                TANK_COUNT_AXIS,
                QuantitySpec::dimensional(Dims::NONE),
                1.0,
                16.0,
            )
            .with_quantity("pH", QuantitySpec::dimensional(Dims::NONE), 9.0, 11.0),
        sources: vec![artifact],
        provenance: Provenance {
            source: declaration,
            license: "project-owned synthetic example".into(),
            artifact: Some(artifact),
        },
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() == 1 && args[0] == "--help" {
        println!(
            "{USAGE}\nActive volume is fixed at nominal flow times 24 h. Source-declared mean residence must remain within 12..48 h; stopped or out-of-domain flow refuses."
        );
        return Ok(());
    }
    if args.len() > 1 {
        return Err(USAGE.into());
    }
    let flow = args
        .first()
        .map(|value| value.parse::<f64>())
        .transpose()?
        .unwrap_or(NOMINAL_FLOW);
    let actual_flow = VolumetricFlowRate::new(flow);
    let active_volume = Volume::new(NOMINAL_FLOW * 24.0 * HOUR);
    // Contained metal flow: 250 kg dry solids/s * 3 g/t * 1e-6 kg/kg per g/t.
    // This scalar balance does not predict a bulk concentrate mass or grade.
    let feed = MassFlowRate::new(250.0 * 3.0e-6);
    let conditions =
        QueryPoint::new().with_quantity("pH", QuantitySpec::dimensional(Dims::NONE), 10.0)?;
    let models = [
        (
            "single",
            ContinuousRecoveryCard::from_card(synthetic_card(false))?,
        ),
        (
            "uniform",
            ContinuousRecoveryCard::from_card(synthetic_card(true))?,
        ),
    ];
    // Admit every requested result before printing a partial output table.
    let mut results = Vec::new();
    for (label, model) in &models {
        for tanks in [1, 8] {
            results.push((
                *label,
                model.evaluate_train(tanks, active_volume, actual_flow, feed, &conditions)?,
            ));
        }
    }
    println!("# Authored continuous treatment illustration; no measured ore/plant calibration.");
    println!(
        "# Fixed active volume={} m3; actual flow={} m3/s; feed=250 kg dry solids/s at authored 3 g/t contained metal.",
        active_volume.value(),
        flow
    );
    println!(
        "# Both models have Rinf=0.9 and mean local rate 1/(6 h): single k=1/(6 h), uniform k in [0,2/(6 h)]."
    );
    println!(
        "# These cards declare stationary local kinetics, with mean-residence domain 12..48 h, 1..16 tanks and pH 9..11. They do not extend finite elapsed-time test evidence."
    );
    for (label, model) in &models {
        println!(
            "# Original {label} source card hash: {}",
            model.source().content_hash()
        );
    }
    println!(
        "model,tanks,actual_flow_m3_s,mean_residence_h,feed_metal_kg_h,recovery_pct,recovered_metal_kg_h,tail_metal_kg_h,conservation_residual_kg_s"
    );
    for (label, result) in results {
        println!(
            "{label},{},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.6e}",
            result.tank_count,
            result.actual_train_flow.value(),
            result.mean_residence_time.value() / HOUR,
            result.feed_metal_flow.value() * HOUR,
            result.recovery.value() * 100.0,
            result.recovered_metal_flow.value() * HOUR,
            result.tail_metal_flow.value() * HOUR,
            result.conservation_residual.value()
        );
    }
    Ok(())
}
