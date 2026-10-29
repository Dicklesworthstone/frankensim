//! Run `cargo run -p fs-gas --features experimental-euler --example mixture_manifold -- 64`.
//! Append `molecular` for cylindrical passages with species diffusion, heat
//! conduction and viscous transport across the finite-plenum connections.
//!
//! Pressurized nitrogen expands into an air-filled compact plenum and two daughter
//! passages. Compression waves, reflected waves and the composition front evolve
//! together; the closed system retains each species and physical total energy.
//! The 200--1000 K N2/O2 NASA-9 coefficients are retained from Cantera v3.2.0
//! data/airNASA9.yaml, blob 5975d1303819332c1bd77a783d9873b66f65c952:
//! https://github.com/Cantera/cantera/blob/v3.2.0/data/airNASA9.yaml
//!
//! Both modes use a rigid perfectly mixed 30 mL plenum and an N2/O2 air surrogate.
//! The default inviscid mode retains one conical daughter. Molecular mode uses
//! three cylinders and an explicit 10 mm plenum-side mixing distance at each
//! port; this compact closure does not resolve three-dimensional junction mixing.
//! Original N2/O2 collision parameters are retained from Cantera v3.2.0
//! data/gri30.yaml, blob 9a62004ec17bbec8ac11e97785932d426dd58610:
//! https://github.com/Cantera/cantera/blob/v3.2.0/data/gri30.yaml
//! Neither mode adds reaction, empirical junction/wall loss or experimental
//! validation. The molecular Eucken conductivity is an estimated closure.

use fs_evidence::ValidityDomain;
use fs_gas::multicomponent::dissipation::KineticTransport;
use fs_gas::multicomponent::network::dissipating::{DissipatingNetwork1d, DissipatingStepControl};
use fs_gas::multicomponent::network::{
    BranchEnd, MixtureNetwork1d, NetworkTotals, Plenum, PlenumConnection,
};
use fs_gas::multicomponent::{
    Boundaries, Boundary, CaloricDomain, Cell, Conserved, Mesh1d, MixtureEuler1d, MixtureGas,
    StepControl,
};
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{DynViscosity, Length, MolarMass, Pressure, Temperature, Time};
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::transport::kinetic::{
    KineticGasTransportV1, KineticTransportDomainV1, LennardJonesSpeciesV1,
};
use fs_thermochem::{
    ChargeVector, ElementId, ElementalMatrix, ElementalReferenceIdV1, Nasa9RegionV1,
    Nasa9StandardStateModelV1, ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1,
    StandardStatePhaseV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R, required_nasa9_parameter_dims,
};
use std::io::{BufWriter, Write};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
const AMBIENT_PRESSURE: f64 = 101_325.0;
const PLENUM_VOLUME: f64 = 30.0e-6;
const PLENUM_MIXING_LENGTH: f64 = 0.010;
const SNAPSHOTS: [f64; 5] = [0.0, 0.0005, 0.001, 0.002, 0.004];

fn main() -> Result<()> {
    let mut arguments = std::env::args().skip(1);
    let cells = arguments
        .next()
        .map(|value| value.parse::<usize>())
        .transpose()?
        .unwrap_or(64);
    let molecular = match arguments.next().as_deref() {
        None => false,
        Some("molecular") => true,
        Some(_) => return Err("usage: mixture_manifold [cell_count] [molecular]".into()),
    };
    if arguments.next().is_some() {
        return Err("usage: mixture_manifold [cell_count] [molecular]".into());
    }
    if !(16..=2048).contains(&cells) {
        return Err("cell count per branch must be in 16..=2048".into());
    }
    let gas = nitrogen_oxygen()?;
    let air = state(&gas, [0.79, 0.21], 300.0, AMBIENT_PRESSURE)?;
    let nitrogen = state(&gas, [1.0, 0.0], 330.0, 180_000.0)?;
    let profiles = [
        (0.6, 0.010, 0.010),
        (0.5, 0.008, 0.008),
        (0.7, 0.008, if molecular { 0.008 } else { 0.012 }),
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
    let network = MixtureNetwork1d::new(
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
    let mut network = if molecular {
        let law = nitrogen_oxygen_transport(gas.source().clone())?;
        eprintln!(
            "Molecular geometry: three stationary cylinders; lengths_m=[0.6,0.5,0.7], radii_m=[0.010,0.008,0.008], rigid_plenum_m3={PLENUM_VOLUME:.6e}, plenum_side_mixing_length_m={PLENUM_MIXING_LENGTH:.6e} at each port. The gradient distance also includes the endpoint half-cell width."
        );
        eprintln!(
            "Molecular source: Cantera v3.2.0 data/gri30.yaml, blob=9a62004ec17bbec8ac11e97785932d426dd58610; N2(sigma_m=3.621e-10,epsilon_over_k_k=97.53), O2(sigma_m=3.458e-10,epsilon_over_k_k=107.4). Calorics remain the original airNASA9 bank."
        );
        eprintln!(
            "Molecular closures: dilute nonpolar Chapman-Enskog/Lennard-Jones coefficients, Wilke viscosity, estimated Eucken conductivity, Maxwell-Stefan species/partial-enthalpy flux, Fourier conduction and longitudinal Newtonian stress; bulk_viscosity_pa_s=0 (Stokes). Declared domain: 200..1000 K, 1000..1e6 Pa. No Soret/Dufour, dense-gas or rotational-relaxation corrections; compact mixing distance is a model input, not resolved plenum geometry."
        );
        Manifold::Molecular(Box::new(DissipatingNetwork1d::new(
            network,
            law,
            vec![Length::new(PLENUM_MIXING_LENGTH); 3],
            || false,
        )?))
    } else {
        Manifold::Inviscid(network)
    };
    let initial = network.totals();
    let mut external_species = [0.0; 2];
    let mut external_energy = 0.0;
    let mut steps = 0_usize;
    let mut maximum_cfl = 0.0_f64;
    let mut fallbacks = 0_usize;
    let mut halvings = 0_usize;
    let mut molecular_steps = 0_usize;
    let mut molecular_rejected_steps = 0_usize;
    let mut coupling_retries = 0_usize;
    let mut maximum_parabolic_number = 0.0_f64;
    let mut molecular_n2_to_plenum = [0.0; 3];
    let mut molecular_energy_to_plenum = [0.0; 3];
    let mut daughter_snapshot_peak_overpressure = [0.0; 2];
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
            let duration = Time::new(target - network.time().value());
            let (report, boundary_change) = match &mut network {
                Manifold::Inviscid(owner) => {
                    let report = owner.step(duration, StepControl::default(), || false)?;
                    (report, report.boundary_change)
                }
                Manifold::Molecular(owner) => {
                    let report =
                        owner.step(duration, DissipatingStepControl::default(), || false)?;
                    coupling_retries += report.coupling_retries;
                    for half in [report.molecular.first, report.molecular.second] {
                        molecular_steps += half.accepted_steps;
                        molecular_rejected_steps += half.rejected_steps;
                        maximum_parabolic_number =
                            maximum_parabolic_number.max(half.maximum_parabolic_number);
                    }
                    for (i, exchange) in owner
                        .accepted_molecular_port_exchanges()
                        .ok_or("missing accepted molecular port report")?
                        .iter()
                        .enumerate()
                    {
                        molecular_n2_to_plenum[i] += exchange.species_masses_to_plenum[0];
                        molecular_energy_to_plenum[i] +=
                            exchange.absolute_energy_to_plenum(owner.gas())?;
                    }
                    (report.transport, report.boundary_change)
                }
            };
            for (i, value) in external_species.iter_mut().enumerate() {
                *value += boundary_change.species_masses[i];
            }
            external_energy += boundary_change.energy;
            maximum_cfl = maximum_cfl.max(report.maximum_cfl);
            fallbacks += report.hll_fallback_faces;
            halvings += report.halvings;
            steps += 1;
        }
        snapshot(
            &mut output,
            &network,
            &mut daughter_snapshot_peak_overpressure,
        )?;
    }
    output.flush()?;
    let final_state = network.totals();
    let residual = std::array::from_fn::<_, 2, _>(|i| {
        final_state.species_masses[i] - initial.species_masses[i] - external_species[i]
    });
    let energy_residual = final_state.energy - initial.energy - external_energy;
    let (gas, time, branches, plenums) = network.view();
    let physical_energy_residual = energy_residual
        + residual
            .iter()
            .zip(gas.energy_references())
            .map(|(dm, e)| dm * e)
            .sum::<f64>();
    let plenum = plenums[0].state();
    let plenum_primitive = gas.primitive(plenum, || false)?;
    let daughter_max_n2 = std::array::from_fn::<_, 2, _>(|i| {
        let branch = &branches[i + 1];
        (0..branch.cell_count())
            .map(|j| {
                let q = branch.cell(j).expect("admitted branch cell");
                q.partial_densities[0] / q.density()
            })
            .fold(0.0_f64, f64::max)
    });
    eprintln!(
        "N2/O2 NASA9 manifold: cells_per_branch={cells}, steps={steps}, time_s={:.9}, maximum_cfl={maximum_cfl:.6e}, hll_fallback_faces={fallbacks}, halvings={halvings}, n2_residual_kg={:.6e}, o2_residual_kg={:.6e}, referenced_energy_residual_j={energy_residual:.6e}, physical_energy_residual_j={physical_energy_residual:.6e}, final_plenum_n2_mass_fraction={:.9}, final_plenum_pressure_pa={:.6}",
        time.value(),
        residual[0],
        residual[1],
        plenum.partial_densities[0] / plenum.density(),
        plenum_primitive.pressure,
    );
    eprintln!(
        "Wave/composition response: daughter_snapshot_peak_overpressure_pa=[{:.6},{:.6}], final_daughter_max_n2_mass_fraction=[{:.9},{:.9}]",
        daughter_snapshot_peak_overpressure[0],
        daughter_snapshot_peak_overpressure[1],
        daughter_max_n2[0],
        daughter_max_n2[1],
    );
    if molecular {
        eprintln!(
            "Molecular network: accepted_substeps={molecular_steps}, rejected_substeps={molecular_rejected_steps}, coupling_retries={coupling_retries}, maximum_parabolic_number={maximum_parabolic_number:.6e}, port_n2_to_plenum_kg=[{:.6e},{:.6e},{:.6e}], port_physical_energy_to_plenum_j=[{:.6e},{:.6e},{:.6e}]",
            molecular_n2_to_plenum[0],
            molecular_n2_to_plenum[1],
            molecular_n2_to_plenum[2],
            molecular_energy_to_plenum[0],
            molecular_energy_to_plenum[1],
            molecular_energy_to_plenum[2],
        );
    }
    Ok(())
}

enum Manifold {
    Inviscid(MixtureNetwork1d),
    Molecular(Box<DissipatingNetwork1d<KineticTransport>>),
}
impl Manifold {
    fn view(&self) -> (&MixtureGas, Time, &[MixtureEuler1d], &[Plenum]) {
        match self {
            Self::Inviscid(owner) => (owner.gas(), owner.time(), owner.branches(), owner.plenums()),
            Self::Molecular(owner) => {
                (owner.gas(), owner.time(), owner.branches(), owner.plenums())
            }
        }
    }
    fn time(&self) -> Time {
        self.view().1
    }
    fn totals(&self) -> NetworkTotals {
        match self {
            Self::Inviscid(owner) => owner.totals(),
            Self::Molecular(owner) => owner.totals(),
        }
    }
}

fn snapshot(
    output: &mut impl Write,
    network: &Manifold,
    peak_overpressure: &mut [f64; 2],
) -> Result<()> {
    let (gas, time, branches, plenums) = network.view();
    for (branch_index, branch) in branches.iter().enumerate() {
        for i in 0..branch.cell_count() {
            let q = branch.cell(i).ok_or("missing admitted branch cell")?;
            let x = branch.mesh().positions()[i].midpoint(branch.mesh().positions()[i + 1]);
            let pressure = write_cell(output, gas, time, "branch", branch_index, x, q)?;
            if branch_index > 0 {
                peak_overpressure[branch_index - 1] =
                    peak_overpressure[branch_index - 1].max(pressure - AMBIENT_PRESSURE);
            }
        }
    }
    for (i, plenum) in plenums.iter().enumerate() {
        write_cell(output, gas, time, "plenum", i, 0.0, plenum.state())?;
    }
    Ok(())
}
fn write_cell(
    output: &mut impl Write,
    gas: &MixtureGas,
    time: Time,
    owner: &str,
    index: usize,
    x: f64,
    q: Cell<'_>,
) -> Result<f64> {
    let p = gas.primitive(q, || false)?;
    writeln!(
        output,
        "{:.9},{owner},{index},{x:.9},{:.9},{:.9},{:.9},{:.12},{:.12}",
        time.value(),
        p.pressure,
        p.temperature,
        p.velocity,
        q.partial_densities[0] / p.density,
        q.partial_densities[1] / p.density
    )?;
    Ok(p.pressure)
}

fn nitrogen_oxygen_transport(source: IdealGasEquilibriumV1) -> Result<KineticTransport> {
    // Original Cantera v3.2.0 GRI3 collision records, Git blob
    // 9a62004ec17bbec8ac11e97785932d426dd58610. Convert angstroms to metres;
    // epsilon/kB is already in kelvin. Molecular masses and Cp(T) remain owned
    // by the exact NASA9 source used by this network, not the GRI NASA7 fits.
    let model = KineticGasTransportV1::new(
        source,
        "gri30-cantera-3.2.0-n2-o2-lennard-jones",
        1,
        vec![
            LennardJonesSpeciesV1::new(
                SpeciesId::new("N2")?,
                Length::new(3.621e-10),
                Temperature::new(97.53),
            )?,
            LennardJonesSpeciesV1::new(
                SpeciesId::new("O2")?,
                Length::new(3.458e-10),
                Temperature::new(107.4),
            )?,
        ],
        KineticTransportDomainV1 {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 1000.0,
            minimum_pressure_pa: 1000.0,
            maximum_pressure_pa: 1e6,
        },
        || false,
    )?;
    // Explicit Stokes approximation; no bulk viscosity is inferred from labels.
    Ok(KineticTransport::new(model, DynViscosity::new(0.0))?)
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
