//! Continuous retained-carbon metal accounts under feed changes and pump stops.
//!
//! Run: `cargo run -p fs-material --example mineral_adsorption_flow`
//! All capacities, affinities, local rates and conditions below are explicitly
//! authored numerical examples. They are not measured gold/carbon data or a
//! calibrated CIP/CIL forecast. Carbon remains retained in one fixed-volume
//! ideal liquid tank; only dissolved analytical gold enters and leaves.

use std::error::Error;

use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance, QueryPoint,
};
use fs_material::aqueous_adsorption::{
    AQUEOUS_LANGMUIR_LAW, AdsorbentLoading, AqueousLangmuir, CAPACITY_PARAMETER,
    CONCENTRATION_AXIS, HALF_SATURATION_PARAMETER, TEMPERATURE_AXIS,
    flow::{
        AdsorptionCstr, AdsorptionFlowBudget, LOADING_AXIS, LOCAL_LDF_LAW, LocalLangmuirTransfer,
        RATE_PARAMETER, RetainedSorbent, VOLUME_AXIS,
    },
};
use fs_qty::{
    Amount, AmountConcentration, Dimensionless, Dims, Frequency, Mass, QuantitySpec, Temperature,
    Time, Volume, VolumetricFlowRate,
    chemistry::ElementId,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

fn source(
    law: &str,
    parameters: &[(&str, f64, Dims)],
    validity: ValidityDomain,
) -> ConstitutiveModelCard {
    let artifact = hash_domain(
        "authored-retained-carbon-flow-example-v1",
        b"synthetic original source declarations; no experimental fit or plant forecast",
    );
    ConstitutiveModelCard {
        law: LawId(law.into()), law_version: 1, state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        parameters: parameters.iter().map(|&(name, value, dims)|
            (name.into(), LawParameter { value, dims })).collect(),
        validity, sources: vec![artifact], provenance: Provenance {
            source: "Authored retained-carbon numerical illustration; not measured mineral-processing data".into(),
            license: "project example".into(), artifact: Some(artifact),
        },
    }
}

fn population(
    id: &str,
    mass: f64,
    loaded: f64,
    half: f64,
    rate: f64,
) -> Result<RetainedSorbent, Box<dyn Error>> {
    let domain = ValidityDomain::unconstrained()
        .with_quantity(
            CONCENTRATION_AXIS,
            QuantitySpec::dimensional(AmountConcentration::DIMS),
            0.0,
            1.0,
        )
        .with_quantity(
            TEMPERATURE_AXIS,
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            298.15,
            298.15,
        )
        .with_quantity("pH", QuantitySpec::dimensional(Dims::NONE), 11.0, 11.0);
    let equilibrium = AqueousLangmuir::from_card(
        ElementId::new("Au")?,
        source(
            AQUEOUS_LANGMUIR_LAW,
            &[
                (CAPACITY_PARAMETER, 0.02, AdsorbentLoading::DIMS),
                (HALF_SATURATION_PARAMETER, half, AmountConcentration::DIMS),
            ],
            domain.clone(),
        ),
    )?;
    let local = LocalLangmuirTransfer::from_card(
        equilibrium,
        source(
            LOCAL_LDF_LAW,
            &[(RATE_PARAMETER, rate, Frequency::DIMS)],
            domain
                .with_quantity(
                    LOADING_AXIS,
                    QuantitySpec::dimensional(AdsorbentLoading::DIMS),
                    0.0,
                    0.02,
                )
                .with_quantity(
                    VOLUME_AXIS,
                    QuantitySpec::dimensional(Volume::DIMS),
                    2.0,
                    2.0,
                ),
        ),
    )?;
    Ok(RetainedSorbent::new(
        id,
        local,
        Mass::new(mass),
        Amount::new(loaded),
        QueryPoint::new().with_quantity("pH", QuantitySpec::dimensional(Dims::NONE), 11.0)?,
    )?)
}

fn main() -> Result<(), Box<dyn Error>> {
    let mut tank = AdsorptionCstr::new(
        vec![
            population("fast-carbon", 2.0, 0.0, 0.0005, 0.02)?,
            population("slow-preloaded-carbon", 4.0, 0.04, 0.002, 0.005)?,
        ],
        Amount::new(0.004),
        Volume::new(2.0),
        Temperature::new(298.15),
    )?;
    let budget = AdsorptionFlowBudget {
        absolute_amount_error: Amount::new(1e-9),
        relative_amount_error: Dimensionless::new(1e-5),
        maximum_step: Time::new(10.0),
        maximum_work: 5_000_000,
        ..AdsorptionFlowBudget::default()
    };
    println!("Synthetic retained-carbon contactor: V=2 m3, T=298.15 K, declared fixed pH=11.");
    println!(
        "Capacities, affinities and rates are authored illustrations, not calibrated recovery data."
    );
    println!(
        "The sorbent stays in the tank. Solution outlet is exported dissolved gold, not loaded-carbon product."
    );
    for (label, flow, concentration, seconds) in [
        ("pregnant solution feed", 0.01, 0.002, 120.0),
        ("richer solution feed", 0.01, 0.01, 120.0),
        (
            "pump stop; retained phases continue exchanging metal",
            0.0,
            0.0,
            120.0,
        ),
        ("restart with metal-free wash", 0.02, 0.0, 300.0),
    ] {
        let report = tank.advance(
            VolumetricFlowRate::new(flow),
            AmountConcentration::new(concentration),
            Time::new(seconds),
            budget,
            || false,
        )?;
        println!("\n{label}: Q={flow} m3/s, inlet={concentration} mol/m3, duration={seconds} s");
        println!(
            "  interval incoming / solution outlet: {:.9} / {:.9} mol",
            report.inlet_amount.value(),
            report.outlet_amount.value()
        );
        println!(
            "  retained dissolved: {:.9} mol; concentration: {:.9} mol/m3",
            report.dissolved.value(),
            report.concentration.value()
        );
        for (state, (_, transfer)) in report.populations.iter().zip(&report.transfers) {
            println!(
                "  {}: {:.9} mol loaded, {:.9} mol vacant, interval transfer {:+.9} mol",
                state.id,
                state.adsorbed.value(),
                state.vacant_capacity.value(),
                transfer.value()
            );
        }
        println!(
            "  cumulative incoming / solution outlet: {:.9} / {:.9} mol",
            report.cumulative_inlet.value(),
            report.cumulative_outlet.value()
        );
        println!(
            "  interval / history metal-balance residual: {:.3e} / {:.3e} mol",
            report.interval_balance_residual.value(),
            report.cumulative_balance_residual.value()
        );
        println!(
            "  accepted/rejected steps: {}/{}; maximum local estimate: {:.3e} mol",
            report.accepted_steps,
            report.rejected_steps,
            report.maximum_local_error.value()
        );
    }
    println!(
        "\nFirst-order implicit integration; error estimates are not certified global bounds."
    );
    println!(
        "This tracks analytical Au only. No ion/charge removal, leach reaction, carbon motion or pH solve is inferred."
    );
    Ok(())
}
