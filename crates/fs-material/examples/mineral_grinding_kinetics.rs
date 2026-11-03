//! Synthetic selection/breakage kinetics through the original particle stream.
//! Run: cargo run -p fs-material --example mineral_grinding_kinetics
//!
//! Every class receives the same declared contact time. This is a constant-rate
//! population-balance illustration, not a stirred-mill residence distribution,
//! measured cassiterite test, fitted plant model or mineral-liberation forecast.

use std::{collections::BTreeMap, error::Error};

use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance, QueryPoint,
};
use fs_material::{
    comminution::{
        breakage::{
            BreakageColumn, BreakageTestContext, BreakageTransfer, SizeClass,
            SizeRedistributionCard,
        },
        kinetics::{
            CONTACT_TIME_AXIS, RATE_PREFIX, SELECTION_BREAKAGE_LAW, SelectionBreakageModel,
            UniformizationBudget,
        },
    },
    separation::{ParticleClass, ParticleStream},
};
use fs_qty::{
    Dimensionless, Frequency, Length, MassFlowRate, QuantitySpec, Time, chemistry::ElementId,
};

fn main() -> Result<(), Box<dyn Error>> {
    let description = "SYNTHETIC example: 300/150/75 um classes; coarse breaks to medium \
        at 0.02/s; medium to fine at 0.01/s; fine is inert. These rates and \
        daughters are authored, not measured ore properties. Dry flows are \
        20/10/5 kg/s, Sn flows 0.4/0.1/0.025 kg/s, and water is 7 kg/s. \
        Constant kinetics are declared at 2 revolutions/s for 0..1000 s.";
    let artifact = hash_domain(
        "synthetic-mineral-grinding-kinetics-v1",
        description.as_bytes(),
    );
    let provenance = Provenance {
        source: description.into(),
        license: "project-authored example".into(),
        artifact: Some(artifact),
    };
    let quantity = QuantitySpec::dimensional;
    let operating_domain = ValidityDomain::unconstrained().with_quantity(
        "rotation-frequency",
        quantity(Frequency::DIMS),
        2.0,
        2.0,
    );
    let daughters = SizeRedistributionCard::new(
        [("coarse", 300e-6), ("medium", 150e-6), ("fine", 75e-6)]
            .into_iter()
            .map(|(id, size)| SizeClass {
                id: id.into(),
                diameter: Length::new(size),
            })
            .collect(),
        [("coarse", "medium"), ("medium", "fine"), ("fine", "fine")]
            .into_iter()
            .map(|(input, output)| BreakageColumn {
                input: input.into(),
                transfers: vec![BreakageTransfer {
                    output: output.into(),
                    mass_fraction: Dimensionless::new(1.0),
                }],
            })
            .collect(),
        BreakageTestContext {
            specimen: "authored tin-bearing three-class stream".into(),
            method: "synthetic congruent daughter chain".into(),
            provenance: provenance.clone(),
        },
        operating_domain.clone(),
    )?;
    let model = SelectionBreakageModel::from_card(
        daughters,
        ConstitutiveModelCard {
            law: LawId(SELECTION_BREAKAGE_LAW.into()),
            law_version: 1,
            state_schema_version: 0,
            initial_state: InitialStatePolicy::ZeroInternalState,
            parameters: [("coarse", 0.02), ("medium", 0.01), ("fine", 0.0)]
                .into_iter()
                .map(|(id, value)| {
                    (
                        format!("{RATE_PREFIX}{id}"),
                        LawParameter {
                            value,
                            dims: Frequency::DIMS,
                        },
                    )
                })
                .collect(),
            validity: operating_domain.with_quantity(
                CONTACT_TIME_AXIS,
                quantity(Time::DIMS),
                0.0,
                1000.0,
            ),
            sources: vec![artifact],
            provenance,
        },
    )?;
    let tin = ElementId::new("Sn")?;
    let inlet = ParticleStream::new(
        [
            ("coarse", 300e-6, 20.0, 0.4),
            ("medium", 150e-6, 10.0, 0.1),
            ("fine", 75e-6, 5.0, 0.025),
        ]
        .into_iter()
        .map(|(id, size, dry, sn)| {
            ParticleClass::new(
                id,
                Length::new(size),
                MassFlowRate::new(dry),
                BTreeMap::from([(tin.clone(), MassFlowRate::new(sn))]),
            )
        })
        .collect::<Result<Vec<_>, _>>()?,
        MassFlowRate::new(7.0),
    )?;
    let conditions =
        QueryPoint::new().with_quantity("rotation-frequency", quantity(Frequency::DIMS), 2.0)?;
    let budget = UniformizationBudget {
        maximum_l1_tail: Dimensionless::new(1e-12),
        max_poisson_terms: 256,
        max_coefficient_products: 1_000_000,
    };
    println!("{description}");
    for seconds in [30.0, 120.0] {
        let pass = model.prepare(Time::new(seconds), &conditions, budget, || false)?;
        let outlet = pass.transform(&inlet)?;
        println!("\ncontact_time_s={seconds:.1}");
        for class in outlet.classes() {
            println!(
                "class={} diameter_um={:.1} dry_kg_s={:.12} Sn_kg_s={:.12}",
                class.id(),
                class.diameter().value() * 1e6,
                class.dry_solids().value(),
                class.elements()[&tin].value()
            );
        }
        for (name, before, after) in [
            (
                "dry",
                inlet.dry_solids()?.value(),
                outlet.dry_solids()?.value(),
            ),
            (
                "Sn",
                inlet.element_flow(&tin)?.value(),
                outlet.element_flow(&tin)?.value(),
            ),
            ("water", inlet.water().value(), outlet.water().value()),
        ] {
            println!(
                "balance={name} inlet_kg_s={before:.12} outlet_kg_s={after:.12} residual_kg_s={:.3e}",
                after - before
            );
            if (after - before).abs() > 1e-12 * before.max(1.0) {
                return Err(format!(
                    "{name} conservation exceeds the example's numerical allowance"
                )
                .into());
            }
        }
        println!(
            "capped_events={} nominal_poisson_tail_bound={:.3e} nominal_L1_tail_bound={:.3e} coefficient_products={}",
            pass.capped_events,
            pass.poisson_tail_bound.value(),
            2.0 * pass.poisson_tail_bound.value(),
            pass.coefficient_products
        );
    }
    println!(
        "\nTail bounds exclude floating-point rounding and source uncertainty; no calibration is claimed."
    );
    Ok(())
}
