//! Run `cargo run -p fs-duct --features nonlinear-euler --example nonlinear_horn -- 160`.
//! CSV snapshots show a finite-amplitude thermal pressure pulse splitting and
//! reflecting in a rigid cylinder/cone horn. The original NASA-9 nitrogen law
//! supplies changing heat capacity and sound speed on its 200–1000 K domain.
//! This inviscid, quasi-1D engineering example is not experimental validation.

use fs_duct::nonlinear::{
    Boundaries, Boundary, CaloricDomain, FrozenGas, InitialGasState, NonlinearDuct, StepControl,
};
use fs_duct::{Duct, Segment};
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_math::det;
use fs_qty::{MolarMass, Pressure, Temperature, Time};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{
    Composition, CompositionBasis, ElementalReferenceIdV1, Nasa9RegionV1,
    Nasa9StandardStateModelV1, ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1,
    StandardStatePhaseV1, required_nasa9_parameter_dims,
};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

fn main() -> Result<()> {
    let cells = std::env::args()
        .nth(1)
        .map(|s| s.parse::<usize>())
        .transpose()?
        .unwrap_or(160);
    if !(32..=4096).contains(&cells) {
        return Err("cell count must be in 32..=4096".into());
    }
    let gas = nitrogen()?;
    let reference_density = 101_325.0 / (gas.gas_constant() * 300.0);
    let geometry = Duct {
        segments: vec![
            Segment::Cylinder {
                radius: 0.01,
                length: 0.35,
            },
            Segment::Cone {
                inlet_radius: 0.01,
                outlet_radius: 0.025,
                length: 0.65,
            },
        ],
    };
    let mut horn = NonlinearDuct::new(
        &geometry,
        gas,
        &[cells / 3, cells - cells / 3],
        Boundaries::Ends {
            left: Boundary::Reflecting,
            right: Boundary::Reflecting,
        },
        |x| {
            Ok(InitialGasState {
                density_kg_m3: reference_density,
                velocity_m_s: 0.0,
                pressure_pa: 101_325.0 + 30_000.0 * det::exp(-0.5 * ((x - 0.55) / 0.035).powi(2)),
            })
        },
    )?;
    let initial = horn.solver().totals();
    let (mut boundary_mass, mut boundary_energy, mut steps) = (0.0, 0.0, 0usize);
    println!("time_s,x_m,area_m2,pressure_pa,temperature_k,velocity_m_s,volume_flow_m3_s");
    for snapshot_time in [0.0, 0.0008, 0.0024] {
        while horn.time().value() < snapshot_time {
            if steps >= 100_000 {
                return Err("step budget exhausted".into());
            }
            let step = horn.step(
                Time::new(snapshot_time - horn.time().value()),
                StepControl::default(),
                || false,
            )?;
            boundary_mass += step.boundary_change.mass;
            boundary_energy += step.boundary_change.energy;
            steps += 1;
        }
        for i in 0..horn.solver().mesh().cell_count() {
            let p = horn.sample(i)?;
            println!(
                "{:.9},{:.9},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}",
                horn.time().value(),
                p.axial_position_m,
                p.area_m2,
                p.pressure_pa,
                p.temperature_k,
                p.velocity_m_s,
                p.volume_velocity_m3_s
            );
        }
    }
    let final_state = horn.solver().totals();
    let mass_residual = final_state.mass - initial.mass - boundary_mass;
    let energy_residual = final_state.energy - initial.energy - boundary_energy;
    eprintln!(
        "N2 NASA9 rigid horn: cells={cells}, steps={steps}, time_s={:.9}, mass_residual_kg={mass_residual:.6e}, referenced_energy_residual_j={energy_residual:.6e}, physical_energy_residual_j={:.6e}",
        horn.time().value(),
        energy_residual + horn.solver().gas().energy_reference() * mass_residual
    );
    Ok(())
}

fn nitrogen() -> Result<FrozenGas> {
    // Same original 200–1000 K N2 coefficients/source association as fs-gas's
    // shock_tube example; no fitted sound speed or constant-gamma replacement.
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
