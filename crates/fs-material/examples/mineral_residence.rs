//! Authored recovery curves and forward mineral hold-up; no measured calibration.
//!
//! Run: `cargo run -p fs-material --example mineral_residence`
//!
//! The hydraulic example transports feed metal through eight ideal full tanks.
//! A separate, explicitly fixed 88% extraction is applied at the outlet. The
//! time-dependent recovery cards are evaluated separately, not coupled to the
//! hydraulic model or used to fit this example. A final, separately declared
//! transient example retains undissolved and dissolved inventories together;
//! its dissolved effluent is not final refining recovery.

use std::{collections::BTreeMap, error::Error};

use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance, QueryPoint,
};
use fs_material::recovery::{
    ELAPSED_TIME_AXIS, MAXIMUM_RATE_PARAMETER, RATE_PARAMETER, RecoveryKinetics, SINGLE_RATE_LAW,
    ULTIMATE_RECOVERY_PARAMETER, UNIFORM_RATE_LAW,
    reactive::{
        ConstituentInventory, LocalTransferCard, ReactiveTankCascade, TRANSFER_RATE_PARAMETER,
        TRANSIENT_TRANSFER_LAW,
    },
    residence::StirredTankCascade,
};
use fs_qty::{Density, Dims, Frequency, Mass, QuantitySpec, Time, Volume, VolumetricFlowRate};

const HOUR: f64 = 3600.0;
const DRY_SOLIDS_PER_VOLUME: f64 = 1000.0; // kg dry solids per m³ of slurry.
const FIXED_OUTLET_EXTRACTION: f64 = 0.88;
const SCENARIO: &str = "Authored GB41-shaped illustration, not observed plant data. \
Eight equal ideal full tanks; active volume=86400 m3; initial metal grade=3 g/dry metric tonne; \
dry solids per slurry volume=1000 kg/m3; actual train flow=1 m3/s; \
four 12-hour feed grades=9,3,3,3 g/dry metric tonne; then 6 hours stopped and \
12 hours at half flow with grade 3; separate downstream metal extraction=0.88. \
No segregation, reaction, slurry volume change or in-tank extraction.";

fn concentration(grade_g_per_tonne: f64) -> Density {
    // 1 g/t = 1e-6 kg metal/kg dry solids, not 1e-3.
    Density::new(grade_g_per_tonne * 1e-6 * DRY_SOLIDS_PER_VOLUME)
}

fn recovery_card(uniform: bool) -> ConstitutiveModelCard {
    let (law, parameter) = if uniform {
        (UNIFORM_RATE_LAW, MAXIMUM_RATE_PARAMETER)
    } else {
        (SINGLE_RATE_LAW, RATE_PARAMETER)
    };
    let rate = 1.0 / (6.0 * HOUR);
    let source = format!(
        "Synthetic authored curve; no ore test or plant fit. law={law}; version=1; \
         Rinf=0.9 fraction; {parameter}={rate:.17e} s^-1; no internal state; \
         elapsed-time validity=[0,172800] seconds; pH validity=[9,11]."
    );
    let artifact = hash_domain("mineral-residence-authored-recovery-v1", source.as_bytes());
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
                ELAPSED_TIME_AXIS,
                QuantitySpec::dimensional(Time::DIMS),
                0.0,
                48.0 * HOUR,
            )
            .with_quantity("pH", QuantitySpec::dimensional(Dims::NONE), 9.0, 11.0),
        sources: vec![artifact],
        provenance: Provenance {
            source,
            license: "project-owned synthetic example".into(),
            artifact: Some(artifact),
        },
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args_os().skip(1).collect();
    if args.len() == 1 && args[0] == "--help" {
        println!(
            "usage: mineral_residence\nPrint an authored forward hold-up CSV and independent source-declared recovery curves."
        );
        return Ok(());
    }
    if !args.is_empty() {
        return Err("usage: mineral_residence [--help]".into());
    }

    let artifact = hash_domain(
        "mineral-residence-authored-scenario-v1",
        SCENARIO.as_bytes(),
    );
    let mut tanks = StirredTankCascade::uniform(
        8,
        Volume::new(86_400.0),
        concentration(3.0),
        Provenance {
            source: SCENARIO.into(),
            license: "project-owned synthetic example".into(),
            artifact: Some(artifact),
        },
    )?;
    println!("# {SCENARIO}");
    println!("# Synthetic scenario source hash: {artifact}");
    println!(
        "# Baseline mean residence=24 h. Grade is g metal/dry metric tonne; inventories and residuals are kg metal."
    );
    println!(
        "# Extraction is applied downstream of hydraulic transport. Inventory and residual concern metal before that separate extraction."
    );
    println!("# A blank mean/recovery denotes no material leaving during a stopped interval.");
    println!(
        "interval,end_h,actual_flow_m3_s,inlet_grade_g_t,transported_mean_g_t,tail_grade_g_t,naive_paired_recovery_pct,inventory_change_kg,hydraulic_residual_kg"
    );
    let mut end_hours = 0.0;
    for (name, hours, flow, grade) in [
        ("pulse", 12.0, 1.0, 9.0),
        ("after_pulse_1", 12.0, 1.0, 3.0),
        ("after_pulse_2", 12.0, 1.0, 3.0),
        ("after_pulse_3", 12.0, 1.0, 3.0),
        ("stopped", 6.0, 0.0, 3.0),
        ("half_flow", 12.0, 0.5, 3.0),
    ] {
        let interval = tanks.advance(
            Time::new(hours * HOUR),
            VolumetricFlowRate::new(flow),
            concentration(grade),
        )?;
        end_hours += hours;
        let (transported, tail, naive) = match interval.outlet_mean {
            Some(mean) => {
                let transported_grade = mean.value() / DRY_SOLIDS_PER_VOLUME * 1e6;
                let tail_grade = transported_grade * (1.0 - FIXED_OUTLET_EXTRACTION);
                (
                    format!("{transported_grade:.12}"),
                    format!("{tail_grade:.12}"),
                    format!("{:.9}", 100.0 * (1.0 - tail_grade / grade)),
                )
            }
            None => (String::new(), String::new(), String::new()),
        };
        println!(
            "{name},{end_hours:.1},{flow:.3},{grade:.3},{transported},{tail},{naive},{:.12},{:.6e}",
            interval.final_inventory.value() - interval.initial_inventory.value(),
            interval.conservation_residual.value()
        );
    }

    let single = RecoveryKinetics::from_card(recovery_card(false))?;
    let uniform = RecoveryKinetics::from_card(recovery_card(true))?;
    let point =
        QueryPoint::new().with_quantity("pH", QuantitySpec::dimensional(Dims::NONE), 10.0)?;
    println!(
        "# Independent authored contact-time curves at pH 10; no automatic coupling to the hydraulic example or its fixed 88% extraction."
    );
    println!(
        "# Rinf=0.9; single k=1/(6 h); uniform rates span 0..1/(6 h). Both cards support 0..48 h and pH 9..11."
    );
    println!(
        "# Original complete single card hash: {}",
        single.source().content_hash()
    );
    println!(
        "# Original complete uniform card hash: {}",
        uniform.source().content_hash()
    );
    println!("contact_time_h,single_recovery_fraction,uniform_recovery_fraction");
    for hours in [0.0, 3.0, 6.0, 12.0, 24.0, 48.0] {
        let time = Time::new(hours * HOUR);
        println!(
            "{hours:.1},{:.12},{:.12}",
            single.evaluate(time, &point)?.value(),
            uniform.evaluate(time, &point)?.value()
        );
    }
    reactive_demo()
}

fn reactive_demo() -> Result<(), Box<dyn Error>> {
    let source = "Authored transient local dissolution of one conserved metal constituent: k=1/(6h), pH9..11; eight equal ideal full tanks, 86400m3 total; initial undissolved grade3g/t, zero dissolved metal; neither phase slips relative to carrier; no reagent limit or plant calibration.";
    let artifact = hash_domain("mineral-residence-authored-reactive-v1", source.as_bytes());
    let law = LocalTransferCard::from_card(ConstitutiveModelCard {
        law: LawId(TRANSIENT_TRANSFER_LAW.into()),
        law_version: 1,
        parameters: BTreeMap::from([(
            TRANSFER_RATE_PARAMETER.into(),
            LawParameter {
                value: 1.0 / (6.0 * HOUR),
                dims: Frequency::DIMS,
            },
        )]),
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with_quantity(
            "pH",
            QuantitySpec::dimensional(Dims::NONE),
            9.0,
            11.0,
        ),
        sources: vec![artifact],
        provenance: Provenance {
            source: source.into(),
            license: "project-owned synthetic example".into(),
            artifact: Some(artifact),
        },
    })?;
    let mut model = ReactiveTankCascade::new(
        "authored metal constituent",
        Volume::new(86_400.0),
        vec![law; 8],
        vec![
            ConstituentInventory {
                undissolved: Mass::new(concentration(3.0).value() * 86_400.0 / 8.0),
                dissolved: Mass::new(0.0),
            };
            8
        ],
    )?;
    let conditions =
        vec![
            QueryPoint::new().with_quantity("pH", QuantitySpec::dimensional(Dims::NONE), 10.0)?;
            8
        ];
    println!("# Separate coupled transient dissolution: {source}");
    println!("# Original local-law source hash: {artifact}");
    println!(
        "# Conversion transfers metal internally; both phases remain until advected. Dissolved outlet is potential downstream feed, not final recovered metal."
    );
    println!(
        "reactive_interval,flow_m3_s,undissolved_inventory_kg,dissolved_inventory_kg,converted_inside_kg,undissolved_outlet_kg,dissolved_outlet_kg,undissolved_residual_kg,dissolved_residual_kg"
    );
    for (name, hours, flow, grade) in [
        ("pulse", 12.0, 1.0, 9.0),
        ("after_pulse", 12.0, 1.0, 3.0),
        ("stopped", 6.0, 0.0, 3.0),
        ("restart_half_flow", 12.0, 0.5, 3.0),
    ] {
        let step = model.advance(
            Time::new(hours * HOUR),
            VolumetricFlowRate::new(flow),
            concentration(grade),
            Density::new(0.0),
            &conditions,
        )?;
        println!(
            "{name},{flow:.3},{:.9},{:.9},{:.9},{:.9},{:.9},{:.6e},{:.6e}",
            step.final_inventory.undissolved.value(),
            step.final_inventory.dissolved.value(),
            step.converted.value(),
            step.outlet.undissolved.value(),
            step.outlet.dissolved.value(),
            step.undissolved_residual.value(),
            step.dissolved_residual.value()
        );
    }
    Ok(())
}
