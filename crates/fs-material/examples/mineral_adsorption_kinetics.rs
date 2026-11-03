//! Run: cargo run -p fs-material --example mineral_adsorption_kinetics
//!
//! Authored analytical-Au exchange between solution, activated carbon, and
//! natural carbon. Capacities, affinities, rates and operating boxes below are
//! synthetic: no measured gold/carbon data or calibrated CIP/CIL prediction.
//! The element identity does not determine its dissolved complex or a mechanism.

use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance, QueryPoint,
};
use fs_material::aqueous_adsorption::{
    AQUEOUS_LANGMUIR_LAW, AdsorbentLoading, AqueousLangmuir, CAPACITY_PARAMETER,
    CONCENTRATION_AXIS, HALF_SATURATION_PARAMETER, TEMPERATURE_AXIS,
    kinetics::{
        AdsorptionBatch, AdsorptionNumerics, ELAPSED_TIME_AXIS, LDF_LAW, LOADING_AXIS,
        LangmuirKinetics, RATE_PARAMETER, SorbentPopulation, VOLUME_AXIS,
    },
};
use fs_qty::{
    Amount, AmountConcentration, Dims, Frequency, Mass, QuantitySpec, Temperature, Time, Volume,
    chemistry::ElementId,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::{collections::BTreeMap, error::Error};

type Result<T> = std::result::Result<T, Box<dyn Error>>;

fn sorbent(
    id: &str,
    mass: f64,
    preloaded: f64,
    capacity: f64,
    half: f64,
    rate: f64,
) -> Result<SorbentPopulation> {
    let declaration = format!(
        "SYNTHETIC {id}: qmax={capacity} mol/kg, c_half={half} mol/m3, \
        independent LDF k={rate} per second; fixed pH=11 and T=298.15 K. No ore calibration."
    );
    let artifact = hash_domain(
        "mineral-adsorption-kinetics-example-v1",
        declaration.as_bytes(),
    );
    let provenance = Provenance {
        source: declaration.clone(),
        license: "project-authored example".into(),
        artifact: Some(artifact),
    };
    let domain = ValidityDomain::unconstrained()
        .with_quantity(
            CONCENTRATION_AXIS,
            QuantitySpec::dimensional(AmountConcentration::DIMS),
            0.0,
            0.02,
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
        ConstitutiveModelCard {
            law: LawId(AQUEOUS_LANGMUIR_LAW.into()),
            law_version: 1,
            state_schema_version: 0,
            initial_state: InitialStatePolicy::ZeroInternalState,
            parameters: BTreeMap::from([
                (
                    CAPACITY_PARAMETER.into(),
                    LawParameter {
                        value: capacity,
                        dims: AdsorbentLoading::DIMS,
                    },
                ),
                (
                    HALF_SATURATION_PARAMETER.into(),
                    LawParameter {
                        value: half,
                        dims: AmountConcentration::DIMS,
                    },
                ),
            ]),
            validity: domain.clone(),
            sources: vec![artifact],
            provenance: provenance.clone(),
        },
    )?;
    let kinetics = LangmuirKinetics::from_card(
        equilibrium,
        ConstitutiveModelCard {
            law: LawId(LDF_LAW.into()),
            law_version: 1,
            state_schema_version: 0,
            initial_state: InitialStatePolicy::ZeroInternalState,
            parameters: BTreeMap::from([(
                RATE_PARAMETER.into(),
                LawParameter {
                    value: rate,
                    dims: Frequency::DIMS,
                },
            )]),
            validity: domain
                .with_quantity(
                    ELAPSED_TIME_AXIS,
                    QuantitySpec::dimensional(Time::DIMS),
                    0.0,
                    1000.0,
                )
                .with_quantity(
                    LOADING_AXIS,
                    QuantitySpec::dimensional(AdsorbentLoading::DIMS),
                    0.0,
                    capacity,
                )
                .with_quantity(
                    VOLUME_AXIS,
                    QuantitySpec::dimensional(Volume::DIMS),
                    1.0,
                    1.0,
                ),
            sources: vec![artifact],
            provenance,
        },
    )?;
    println!("{declaration} source_hash={artifact}");
    Ok(SorbentPopulation::new(
        id,
        kinetics,
        Mass::new(mass),
        Amount::new(preloaded),
        QueryPoint::new().with_quantity("pH", QuantitySpec::dimensional(Dims::NONE), 11.0)?,
    )?)
}

fn run(label: &str, dissolved: f64, preloaded_natural: f64) -> Result<()> {
    println!("\nscenario={label}; solution=1 m3; total analytical Au=0.010 mol");
    let activated = sorbent("activated-carbon", 1.0, 0.0, 0.025, 0.002, 0.03)?;
    let natural = sorbent("natural-carbon", 2.0, preloaded_natural, 0.005, 0.001, 0.12)?;
    let mut bath = AdsorptionBatch::new(
        vec![natural, activated],
        Amount::new(dissolved),
        Volume::new(1.0),
        Temperature::new(298.15),
    )?;
    for seconds in [5.0, 20.0, 75.0] {
        let result = bath.advance(
            Time::new(seconds),
            AdsorptionNumerics {
                absolute_amount_error: Amount::new(1e-10),
                relative_amount_error: 1e-6,
                maximum_step: Time::new(5.0),
                max_work: 2_000_000,
                ..AdsorptionNumerics::default()
            },
            || false,
        )?;
        println!(
            "elapsed_s={:.1} previous_s={:.1} dissolved_mol={:.10} concentration_mol_m3={:.10}",
            result.elapsed.value(),
            result.previous_elapsed.value(),
            result.dissolved.value(),
            result.concentration.value()
        );
        for (state, (_, transfer)) in result.populations.iter().zip(&result.transfers) {
            println!(
                "  sorbent={} adsorbed_mol={:.10} vacant_mol={:.10} loading_mol_kg={:.10} interval_transfer_mol={:+.10} cumulative_transfer_mol={:+.10}",
                state.id,
                state.adsorbed.value(),
                state.vacant_capacity.value(),
                state.loading.value(),
                transfer.value(),
                state.cumulative_transfer.value()
            );
        }
        println!(
            "metal_balance_residual_mol={:.3e} max_local_L1_estimate_mol={:.3e} accumulated_local_estimates_mol={:.3e} accepted={} rejected={} work={}",
            result.conservation_residual.value(),
            result.maximum_local_error.value(),
            result.accumulated_error_estimate.value(),
            result.accepted_steps,
            result.rejected_steps,
            result.work
        );
        if result.conservation_residual.value().abs() > 1e-12 * bath.total_metal().value() {
            return Err("finite analytical-metal balance failed".into());
        }
    }
    Ok(())
}

fn main() -> Result<()> {
    println!("Authored reversible finite-bath LDF kinetics using original Langmuir source owners.");
    println!(
        "Every sorbent sees the same changing solution; negative signed transfer is desorption."
    );
    run("initially-unloaded competing carbons", 0.010, 0.0)?;
    run("preloaded natural carbon can release metal", 0.002, 0.008)?;
    println!(
        "\nClocks and site stocks persist across calls. Local step estimates are not certified global error bounds. No aqueous speciation, leaching, heat, flow or CIP/CIL calibration is implied."
    );
    Ok(())
}
