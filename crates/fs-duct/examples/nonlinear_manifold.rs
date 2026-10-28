//! Run `cargo run -p fs-duct --features nonlinear-euler --example nonlinear_manifold -- 64`.
//!
//! A finite-amplitude nitrogen pulse travels down a feeder duct, charges a rigid
//! 30 mL plenum, and enters two daughter ducts with different original geometry.
//! Exterior walls reflect the waves back through the common finite-volume owner.
//! The original NASA-9 nitrogen source supplies changing Cp, temperature and sound
//! speed; one coupled SSP RK2 step transfers each port's mass and energy once.
//!
//! This is an inviscid quasi-1D engineering example with an acoustically compact,
//! perfectly mixed plenum. It does not resolve a three-dimensional junction or
//! establish experimental shock accuracy. The declared source range is 200–1000 K.

use fs_duct::nonlinear::{
    Boundaries, Boundary, CaloricDomain, FrozenGas, InitialGasState, NonlinearDuct, StepControl,
};
use fs_duct::{Duct, Segment};
use fs_evidence::ValidityDomain;
use fs_gas::dynamics::network::{BranchEnd, GasNetwork1d, Plenum, PlenumConnection};
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

const REFERENCE_PRESSURE_PA: f64 = 101_325.0;
const REFERENCE_TEMPERATURE_K: f64 = 300.0;
const PLENUM_VOLUME_M3: f64 = 30.0e-6;
const SNAPSHOT_TIMES_S: [f64; 5] = [0.0, 0.0006, 0.0012, 0.0024, 0.0040];

fn main() -> Result<()> {
    let cells = std::env::args()
        .nth(1)
        .map(|s| s.parse::<usize>())
        .transpose()?
        .unwrap_or(64);
    if !(16..=2048).contains(&cells) {
        return Err("cell count per branch must be in 16..=2048".into());
    }
    let gas = nitrogen()?;
    let density = REFERENCE_PRESSURE_PA / (gas.gas_constant() * REFERENCE_TEMPERATURE_K);
    let geometry = [
        Duct {
            segments: vec![Segment::Cylinder {
                radius: 0.010,
                length: 0.6,
            }],
        },
        Duct {
            segments: vec![Segment::Cylinder {
                radius: 0.008,
                length: 0.5,
            }],
        },
        Duct {
            segments: vec![Segment::Cone {
                inlet_radius: 0.008,
                outlet_radius: 0.012,
                length: 0.7,
            }],
        },
    ];
    let mut branches = Vec::with_capacity(geometry.len());
    let mut coordinates = Vec::with_capacity(geometry.len());
    for (branch, duct) in geometry.iter().enumerate() {
        let image = NonlinearDuct::new(
            duct,
            gas.clone(),
            &[cells],
            Boundaries::Ends {
                left: Boundary::Reflecting,
                right: Boundary::Reflecting,
            },
            |x| {
                let overpressure = if branch == 0 {
                    30_000.0 * det::exp(-0.5 * ((x - 0.20) / 0.035).powi(2))
                } else {
                    0.0
                };
                Ok(InitialGasState {
                    density_kg_m3: density,
                    velocity_m_s: 0.0,
                    pressure_pa: REFERENCE_PRESSURE_PA + overpressure,
                })
            },
        )?;
        // Preserve the original geometric midpoint areas before transferring
        // ownership. Arithmetic interpolation of face areas is wrong for cones.
        coordinates.push(
            image
                .samples()?
                .into_iter()
                .map(|sample| (sample.axial_position_m, sample.area_m2))
                .collect::<Vec<_>>(),
        );
        branches.push(image.into_solver());
    }
    let plenum = Plenum::new(
        PLENUM_VOLUME_M3,
        gas.conserved(density, 0.0, REFERENCE_PRESSURE_PA)?,
    )?;
    let mut network = GasNetwork1d::new(
        branches,
        vec![plenum],
        vec![
            PlenumConnection {
                branch: 0,
                end: BranchEnd::Right,
                plenum: 0,
            },
            PlenumConnection {
                branch: 1,
                end: BranchEnd::Left,
                plenum: 0,
            },
            PlenumConnection {
                branch: 2,
                end: BranchEnd::Left,
                plenum: 0,
            },
        ],
        || false,
    )?;
    let initial = network.totals();
    let (mut boundary_mass, mut boundary_energy, mut steps) = (0.0, 0.0, 0usize);
    let mut maximum_cfl = 0.0_f64;
    let mut hll_fallback_faces = 0usize;
    let mut rejected_steps = 0usize;
    println!(
        "time_s,component,id,cell,x_m,area_m2,volume_m3,density_kg_m3,pressure_pa,temperature_k,velocity_m_s,volume_flow_m3_s"
    );
    for snapshot_time in SNAPSHOT_TIMES_S {
        while network.time().value() < snapshot_time {
            if steps >= 250_000 {
                return Err("coupled step budget exhausted".into());
            }
            let step = network.step(
                Time::new(snapshot_time - network.time().value()),
                StepControl::default(),
                || false,
            )?;
            boundary_mass += step.boundary_change.mass;
            boundary_energy += step.boundary_change.energy;
            maximum_cfl = maximum_cfl.max(step.maximum_cfl);
            hll_fallback_faces += step.hll_fallback_faces;
            rejected_steps += step.halvings;
            steps += 1;
        }
        for (branch_index, branch) in network.branches().iter().enumerate() {
            for (cell, &state) in branch.states().iter().enumerate() {
                let primitive = network.gas().primitive(state, || false)?;
                let (position, area) = coordinates[branch_index][cell];
                println!(
                    "{:.9},branch,{branch_index},{cell},{position:.9},{area:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}",
                    network.time().value(),
                    branch.mesh().volumes()[cell],
                    primitive.density,
                    primitive.pressure,
                    primitive.temperature,
                    primitive.velocity,
                    area * primitive.velocity,
                );
            }
        }
        for (plenum_index, plenum) in network.plenums().iter().enumerate() {
            let primitive = network.gas().primitive(plenum.state(), || false)?;
            // A lumped plenum has volume and a thermodynamic state, but no
            // meaningful axial coordinate, cross section or resolved velocity.
            println!(
                "{:.9},plenum,{plenum_index},,,,{:.12e},{:.12e},{:.12e},{:.12e},,",
                network.time().value(),
                plenum.volume_m3(),
                primitive.density,
                primitive.pressure,
                primitive.temperature,
            );
        }
    }
    let final_state = network.totals();
    let mass_residual = final_state.mass - initial.mass - boundary_mass;
    let energy_residual = final_state.energy - initial.energy - boundary_energy;
    let physical_energy_residual =
        energy_residual + network.gas().energy_reference() * mass_residual;
    eprintln!(
        "N2 NASA9 rigid manifold: branches=3, cells_per_branch={cells}, plenum_volume_m3={PLENUM_VOLUME_M3:.6e}, steps={steps}, time_s={:.9}, maximum_cfl={maximum_cfl:.6e}, hll_fallback_faces={hll_fallback_faces}, halvings={rejected_steps}, mass_residual_kg={mass_residual:.6e}, referenced_energy_residual_j={energy_residual:.6e}, physical_energy_residual_j={physical_energy_residual:.6e}",
        network.time().value(),
    );
    Ok(())
}

fn nitrogen() -> Result<FrozenGas> {
    // Same original 200–1000 K N2 coefficients/source association as the existing
    // nonlinear_horn and fs-gas shock_tube examples, preserving the source law.
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
