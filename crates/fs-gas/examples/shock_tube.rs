//! Run `cargo run -p fs-gas --features experimental-euler --example shock_tube -- 240`.
//! Emits CSV for a one-metre nitrogen shock tube, using the original NASA-9 owner.
//! This demonstrates an Estimated frozen ideal-gas simulation, not experimental validation.

use fs_evidence::ValidityDomain;
use fs_gas::dynamics::{Boundaries, Boundary, CaloricDomain, Euler1d, FrozenGas, StepControl};
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{Length, MolarMass, Pressure, Temperature, Time};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{
    Composition, CompositionBasis, ElementalReferenceIdV1, Nasa9RegionV1,
    Nasa9StandardStateModelV1, ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1,
    StandardStatePhaseV1, required_nasa9_parameter_dims,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cells = std::env::args()
        .nth(1)
        .map(|s| s.parse::<usize>())
        .transpose()?
        .unwrap_or(240);
    if !(2..=fs_gas::dynamics::MAX_EULER_CELLS).contains(&cells) {
        return Err("invalid cell count".into());
    }
    let gas = nitrogen()?;
    let pressure = gas.gas_constant() * 600.0;
    let duration = 0.2 / pressure.sqrt();
    let states = (0..cells)
        .map(|i| {
            if i < cells / 2 {
                gas.conserved(1.0, 0.0, pressure)
            } else {
                gas.conserved(0.125, 0.0, 0.1 * pressure)
            }
        })
        .collect::<Result<Vec<_>, _>>()?;
    let mut flow = Euler1d::new(
        gas,
        Length::new(1.0 / cells as f64),
        0.01,
        states,
        Boundaries::Ends {
            left: Boundary::Transmissive,
            right: Boundary::Transmissive,
        },
        || false,
    )?;
    let initial = flow.totals();
    let mut boundary_energy = 0.0;
    let mut boundary_mass = 0.0;
    let mut steps = 0;
    while flow.time().value() < duration {
        if steps >= 100_000 {
            return Err("step budget exhausted".into());
        }
        let step = flow.step(
            Time::new(duration - flow.time().value()),
            StepControl::default(),
            || false,
        )?;
        boundary_energy += step.boundary_change.energy;
        boundary_mass += step.boundary_change.mass;
        steps += 1;
    }
    let final_state = flow.totals();
    eprintln!(
        "N2 frozen NASA9 shock tube: cells={cells}, steps={steps}, time_s={:.12e}, mass_residual_kg={:.6e}, referenced_energy_residual_j={:.6e}",
        flow.time().value(),
        final_state.mass - initial.mass - boundary_mass,
        final_state.energy - initial.energy - boundary_energy
    );
    println!("x_m,density_kg_m3,velocity_m_s,pressure_pa,temperature_k");
    for (i, &q) in flow.states().iter().enumerate() {
        let p = flow.gas().primitive(q, || false)?;
        println!(
            "{:.9},{:.12e},{:.12e},{:.12e},{:.12e}",
            (i as f64 + 0.5) / cells as f64,
            p.density,
            p.velocity,
            p.pressure,
            p.temperature
        );
    }
    Ok(())
}

fn nitrogen() -> Result<FrozenGas, Box<dyn std::error::Error>> {
    // Same 200–1000 K nitrogen coefficients retained in fs-thermochem's source tests.
    let a = [
        2.210371497e4,
        -3.818461820e2,
        6.082738360,
        -8.530914410e-3,
        1.384646189e-5,
        -9.625793620e-9,
        2.519705809e-12,
        7.108460860e2,
        -1.076003744e1,
    ];
    let mut parameters = std::collections::BTreeMap::new();
    for ((name, dims), value) in required_nasa9_parameter_dims().into_iter().take(9).zip(a) {
        parameters.insert(name.to_owned(), LawParameter { value, dims });
    }
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: 100_000.0,
            dims: Pressure::DIMS,
        },
    );
    let card = ConstitutiveModelCard {
        law: LawId("nasa9-standard-state".into()),
        law_version: 1,
        parameters,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 200.0, 1000.0),
        sources: Vec::new(),
        provenance: Provenance {
            source: "NASA/TP-2002-211556 N2 example association".into(),
            license: "US government technical report".into(),
            artifact: None,
        },
    };
    let convention = StandardStateConventionV1::new(
        StandardStatePhaseV1::Gas,
        ReferenceEquationOfStateV1::IdealGas,
        Pressure::new(100_000.0),
        ElementalReferenceIdV1::new("nasa-glenn-2002-reference-elements")?,
    )?;
    let source = Nasa9StandardStateModelV1::new(
        SpeciesId::new("N2")?,
        MolarMass::new(0.0280134),
        convention,
        vec![Nasa9RegionV1::from_card(card)?],
    )?;
    let mixture = FrozenIdealGasMixtureModelV1::new(
        vec![source],
        Composition::new(CompositionBasis::MoleFraction, vec![1.0])?,
    )?;
    Ok(FrozenGas::new(
        mixture,
        CaloricDomain::new(Temperature::new(200.0), Temperature::new(1000.0)),
        || false,
    )?)
}
