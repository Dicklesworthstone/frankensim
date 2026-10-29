//! Run `cargo run -p fs-gas --features experimental-euler --example mixture_manifold -- 64`.
//!
//! Pressurized nitrogen expands into an air-filled compact plenum and two daughter
//! passages. Compression waves, reflected waves and the composition front evolve
//! together; the closed system retains each species and physical total energy.
//! The 200--1000 K N2/O2 NASA-9 coefficients are retained from Cantera v3.2.0
//! data/airNASA9.yaml, blob 5975d1303819332c1bd77a783d9873b66f65c952:
//! https://github.com/Cantera/cantera/blob/v3.2.0/data/airNASA9.yaml
//!
//! This inviscid quasi-1D example uses a rigid perfectly mixed 30 mL plenum.
//! The air surrogate contains only N2 and O2; it adds no reaction or diffusion,
//! empirical junction loss, three-dimensional mixing or experimental validation.

use fs_evidence::ValidityDomain;
use fs_gas::multicomponent::network::{BranchEnd, MixtureNetwork1d, Plenum, PlenumConnection};
use fs_gas::multicomponent::{
    Boundaries, Boundary, CaloricDomain, Cell, Conserved, Mesh1d, MixtureEuler1d, MixtureGas,
    StepControl,
};
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{MolarMass, Pressure, Temperature, Time};
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::{
    ChargeVector, ElementId, ElementalMatrix, ElementalReferenceIdV1, Nasa9RegionV1,
    Nasa9StandardStateModelV1, ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1,
    StandardStatePhaseV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R, required_nasa9_parameter_dims,
};
use std::io::{BufWriter, Write};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
const AMBIENT_PRESSURE: f64 = 101_325.0;
const PLENUM_VOLUME: f64 = 30.0e-6;
const SNAPSHOTS: [f64; 5] = [0.0, 0.0005, 0.001, 0.002, 0.004];

fn main() -> Result<()> {
    let cells = std::env::args()
        .nth(1)
        .map(|value| value.parse::<usize>())
        .transpose()?
        .unwrap_or(64);
    if !(16..=2048).contains(&cells) {
        return Err("cell count per branch must be in 16..=2048".into());
    }
    let gas = nitrogen_oxygen()?;
    let air = state(&gas, [0.79, 0.21], 300.0, AMBIENT_PRESSURE)?;
    let nitrogen = state(&gas, [1.0, 0.0], 330.0, 180_000.0)?;
    let profiles = [
        (0.6, 0.010, 0.010),
        (0.5, 0.008, 0.008),
        (0.7, 0.008, 0.012),
    ];
    let mut branches = Vec::with_capacity(profiles.len());
    for (i, &(length, inlet_radius, outlet_radius)) in profiles.iter().enumerate() {
        let mesh = conical_mesh(cells, length, inlet_radius, outlet_radius)?;
        branches.push(MixtureEuler1d::new_on_mesh(
            gas.clone(),
            mesh,
            vec![
                if i == 0 {
                    nitrogen.clone()
                } else {
                    air.clone()
                };
                cells
            ],
            Boundaries::Ends {
                left: Boundary::Reflecting,
                right: Boundary::Reflecting,
            },
            || false,
        )?);
    }
    let mut network = MixtureNetwork1d::new(
        branches,
        vec![Plenum::new(PLENUM_VOLUME, air)?],
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
    let mut external_species = [0.0; 2];
    let mut external_energy = 0.0;
    let mut steps = 0_usize;
    let mut maximum_cfl = 0.0_f64;
    let mut fallbacks = 0_usize;
    let mut halvings = 0_usize;
    let mut output = BufWriter::new(std::io::stdout().lock());
    writeln!(
        output,
        "time_s,owner,index,x_m,pressure_pa,temperature_k,velocity_m_s,n2_mass_fraction,o2_mass_fraction"
    )?;
    for target in SNAPSHOTS {
        while network.time().value() < target {
            if steps >= 200_000 {
                return Err("accepted-step budget exhausted".into());
            }
            let report = network.step(
                Time::new(target - network.time().value()),
                StepControl::default(),
                || false,
            )?;
            for (i, value) in external_species.iter_mut().enumerate() {
                *value += report.boundary_change.species_masses[i];
            }
            external_energy += report.boundary_change.energy;
            maximum_cfl = maximum_cfl.max(report.maximum_cfl);
            fallbacks += report.hll_fallback_faces;
            halvings += report.halvings;
            steps += 1;
        }
        snapshot(&mut output, &network)?;
    }
    output.flush()?;
    let final_state = network.totals();
    let residual = std::array::from_fn::<_, 2, _>(|i| {
        final_state.species_masses[i] - initial.species_masses[i] - external_species[i]
    });
    let energy_residual = final_state.energy - initial.energy - external_energy;
    let physical_energy_residual = energy_residual
        + residual
            .iter()
            .zip(network.gas().energy_references())
            .map(|(dm, e)| dm * e)
            .sum::<f64>();
    let plenum = network.plenums()[0].state();
    let plenum_primitive = network.gas().primitive(plenum, || false)?;
    eprintln!(
        "N2/O2 NASA9 manifold: cells_per_branch={cells}, steps={steps}, time_s={:.9}, maximum_cfl={maximum_cfl:.6e}, hll_fallback_faces={fallbacks}, halvings={halvings}, n2_residual_kg={:.6e}, o2_residual_kg={:.6e}, referenced_energy_residual_j={energy_residual:.6e}, physical_energy_residual_j={physical_energy_residual:.6e}, final_plenum_n2_mass_fraction={:.9}, final_plenum_pressure_pa={:.6}",
        network.time().value(),
        residual[0],
        residual[1],
        plenum.partial_densities[0] / plenum.density(),
        plenum_primitive.pressure,
    );
    Ok(())
}

fn snapshot(output: &mut impl Write, network: &MixtureNetwork1d) -> Result<()> {
    for (branch_index, branch) in network.branches().iter().enumerate() {
        for i in 0..branch.cell_count() {
            let q = branch.cell(i).ok_or("missing admitted branch cell")?;
            let x = branch.mesh().positions()[i].midpoint(branch.mesh().positions()[i + 1]);
            write_cell(output, network, "branch", branch_index, x, q)?;
        }
    }
    for (i, plenum) in network.plenums().iter().enumerate() {
        write_cell(output, network, "plenum", i, 0.0, plenum.state())?;
    }
    Ok(())
}
fn write_cell(
    output: &mut impl Write,
    network: &MixtureNetwork1d,
    owner: &str,
    index: usize,
    x: f64,
    q: Cell<'_>,
) -> Result<()> {
    let p = network.gas().primitive(q, || false)?;
    writeln!(
        output,
        "{:.9},{owner},{index},{x:.9},{:.9},{:.9},{:.9},{:.12},{:.12}",
        network.time().value(),
        p.pressure,
        p.temperature,
        p.velocity,
        q.partial_densities[0] / p.density,
        q.partial_densities[1] / p.density
    )?;
    Ok(())
}
fn conical_mesh(cells: usize, length: f64, inlet: f64, outlet: f64) -> Result<Mesh1d> {
    let x: Vec<_> = (0..=cells)
        .map(|i| length * i as f64 / cells as f64)
        .collect();
    let radii: Vec<_> = x
        .iter()
        .map(|&z| inlet + (outlet - inlet) * (z / length))
        .collect();
    let areas = radii
        .iter()
        .map(|r| core::f64::consts::PI * r * r)
        .collect();
    let volumes = (0..cells)
        .map(|i| {
            core::f64::consts::PI
                * (x[i + 1] - x[i])
                * (radii[i] * radii[i] + radii[i] * radii[i + 1] + radii[i + 1] * radii[i + 1])
                / 3.0
        })
        .collect();
    Ok(Mesh1d::new(x, areas, volumes)?)
}
fn state(
    gas: &MixtureGas,
    mole_fractions: [f64; 2],
    temperature: f64,
    pressure: f64,
) -> Result<Conserved> {
    let concentration = pressure / (R * temperature);
    let partial = std::array::from_fn::<_, 2, _>(|i| {
        concentration * mole_fractions[i] * gas.source().models()[i].molar_mass().value()
    });
    Ok(gas.conserved(&partial, 0.0, pressure)?)
}
fn nitrogen_oxygen() -> Result<MixtureGas> {
    let models = vec![
        source(
            "N2",
            0.0280134,
            [
                2.210371497e4,
                -381.846182,
                6.08273836,
                -8.53091441e-3,
                1.384646189e-5,
                -9.62579362e-9,
                2.519705809e-12,
                710.846086,
                -10.76003744,
            ],
        )?,
        source(
            "O2",
            0.0319988,
            [
                -3.42556342e4,
                484.700097,
                1.119010961,
                4.29388924e-3,
                -6.83630052e-7,
                -2.0233727e-9,
                1.039040018e-12,
                -3391.45487,
                18.4969947,
            ],
        )?,
    ];
    let ids: Vec<_> = models.iter().map(|model| model.species().clone()).collect();
    let bank = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("N")?, ElementId::new("O")?],
            ids.clone(),
            vec![vec![2, 0], vec![0, 2]],
        )?,
        ChargeVector::new(ids.clone(), vec![0; ids.len()])?,
        models,
    )?;
    Ok(MixtureGas::new(
        bank,
        CaloricDomain::new(Temperature::new(200.0), Temperature::new(1000.0)),
        || false,
    )?)
}
fn source(name: &str, mass: f64, coefficients: [f64; 9]) -> Result<Nasa9StandardStateModelV1> {
    let mut parameters = std::collections::BTreeMap::new();
    for ((key, dims), value) in required_nasa9_parameter_dims()
        .into_iter()
        .take(9)
        .zip(coefficients)
    {
        parameters.insert(key.to_owned(), LawParameter { value, dims });
    }
    // Cantera's omitted species reference-pressure defaults to 1 atm:
    // https://cantera.org/3.2/yaml/species.html#species-thermo-models
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: AMBIENT_PRESSURE,
            dims: Pressure::DIMS,
        },
    );
    let card = ConstitutiveModelCard {
        law: LawId("nasa9-standard-state".into()), law_version: 1, parameters,
        state_schema_version: 0, initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 200.0, 1000.0), sources: vec![],
        provenance: Provenance {
            source: "NASA/TP-2002-211556; Cantera v3.2.0 data/airNASA9.yaml; blob 5975d1303819332c1bd77a783d9873b66f65c952".into(),
            license: "NASA technical report; Cantera BSD-3-Clause distribution".into(), artifact: None,
        },
    };
    Ok(Nasa9StandardStateModelV1::new(
        SpeciesId::new(name)?,
        MolarMass::new(mass),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(AMBIENT_PRESSURE),
            ElementalReferenceIdV1::new("nasa-glenn-2002-reference-elements")?,
        )?,
        vec![Nasa9RegionV1::from_card(card)?],
    )?)
}
