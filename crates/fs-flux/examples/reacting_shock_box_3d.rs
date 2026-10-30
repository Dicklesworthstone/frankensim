//! Run `cargo run --release -p fs-flux --no-default-features --features experimental-compressible --example reacting_shock_box_3d -- 4 100`.
//! Append `molecular` to include species diffusion, Fourier heat and full viscous
//! stress/work using the same source-bound three-dimensional field.
//! A warm lean H2/O2/N2 mixture contains a hotter pocket at equal pressure.
//! Original finite-rate heat release launches pressure waves in a rigid 3-D box.
//! The optional transport is a nonpolar Chapman--Enskog/Wilke/Eucken
//! approximation, including for polar H2O; polar corrections, rotational
//! relaxation and Soret/Dufour effects are absent. This coarse stationary
//! example is not a resolved flame/detonation calculation, an experimental
//! validation, or a physical shock-thickness claim.
//! The retained atomic-H Cp polynomial needs an explicitly admitted source-fit
//! deficit of 2e-8 J/(mol K) below the classical translational limit; the
//! original source Cp remains unchanged in calorics and conductivity.

use fs_flux::compressible::spatial::mixture::molecular::{
    DissipationControl, MolecularTransport3d, SpatialTransport,
};
use fs_flux::compressible::spatial::mixture::reactive::{
    CellChemistry, IsochoricKineticsSolve, ReactingEuler3d, ReactiveStepControl,
};
use fs_flux::compressible::spatial::mixture::{MixtureEuler3d, MixtureGas, Totals, state_into};
use fs_flux::compressible::spatial::{
    AxisBoundaries, Boundaries3d, BoundaryCondition, BoundarySide, CartesianGrid,
};
use fs_gas::chemistry::IsochoricKineticsMethod;
use fs_gas::multicomponent::CaloricDomain;
use fs_gas::multicomponent::dissipation::KineticTransport;
use fs_qty::{DynViscosity, Length, Temperature, Time};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use fs_thermochem::kinetics::hydrogen::{HydrogenCalorics, SOURCE_COMMIT, cantera_h2o2};
use fs_thermochem::transport::kinetic::{
    KineticGasTransportV1, KineticTransportDomainV1, LennardJonesSpeciesV1,
};
use fs_thermochem::{MolarThermalQuantityV1, SpeciesId};
use std::io::{BufWriter, Write};

const CP_FIT_DEFICIT_ALLOWANCE_J_MOL_K: f64 = 2e-8;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let n = args
        .next()
        .map(|s| s.parse::<usize>())
        .transpose()?
        .unwrap_or(4);
    let duration_us = args
        .next()
        .map(|s| s.parse::<f64>())
        .transpose()?
        .unwrap_or(100.0);
    let transport_mode = args.next();
    if !(2..=8).contains(&n)
        || !(1.0..=200.0).contains(&duration_us)
        || transport_mode
            .as_deref()
            .is_some_and(|mode| mode != "molecular")
        || args.next().is_some()
    {
        return Err(
            "usage: reacting_shock_box_3d [cells_per_axis in 2..=8] [duration_us in 1..=200] [molecular]"
                .into(),
        );
    }
    let use_molecular = transport_mode.is_some();
    let mechanism = cantera_h2o2(HydrogenCalorics::ContinuousAt300K)?;
    let gas = MixtureGas::new(
        mechanism.thermodynamics().clone(),
        CaloricDomain::new(Temperature::new(300.0), Temperature::new(3500.0)),
        || false,
    )?;
    let width = gas.species_count() + 4;
    let h = 0.1 / n as f64;
    let grid = CartesianGrid::new(n, n, n, h, h, h)?;
    let mut states = vec![0.0; n * n * n * width];
    let mut densities = vec![0.0; gas.species_count()];
    let mut hot_cells = 0;
    for (index, row) in states.chunks_exact_mut(width).enumerate() {
        let [x, y, z] = grid.cell_center(index).expect("valid cell index");
        let hot =
            (x - 0.035).powi(2) + (y - 0.045).powi(2) + (z - 0.055).powi(2) < 0.035_f64.powi(2);
        hot_cells += usize::from(hot);
        let temperature = if hot { 1200.0 } else { 900.0 };
        let concentration_unit = 101_325.0 / (R * temperature * 5.76);
        for (density, model) in densities.iter_mut().zip(gas.source().models()) {
            let amount_ratio = match model.species().as_str() {
                "H2" | "O2" => 1.0,
                "N2" => 3.76,
                _ => 0.0, // No radical, product or inert-species seed.
            };
            *density = concentration_unit * amount_ratio * model.molar_mass().value();
        }
        state_into(&gas, &densities, [0.0; 3], 101_325.0, row)?;
    }
    let walls = AxisBoundaries::Ends {
        lower: BoundaryCondition::Reflecting,
        upper: BoundaryCondition::Reflecting,
    };
    let bytes = MixtureEuler3d::required_workspace_bytes(&gas, &grid)?;
    let field = MixtureEuler3d::new(
        gas,
        grid,
        states,
        Boundaries3d {
            x: walls,
            y: walls,
            z: walls,
        },
        bytes,
        || false,
    )?;
    let volume_scale = grid.cell_volume() / 0.001;
    let chemistry = IsochoricKineticsSolve {
        method: IsochoricKineticsMethod::Sdirk2,
        minimum_temperature_k: 300.0,
        maximum_temperature_k: 3500.0,
        relative_tolerance: 1e-4,
        absolute_amount_tolerance_mol: 1e-17 * volume_scale,
        energy_tolerance_j: 1e-8 * volume_scale,
        entropy_tolerance_j_k: 1e-9 * volume_scale,
        initial_step_s: 1e-8,
        maximum_step_s: 5e-6,
        maximum_steps: 20_000,
        ..IsochoricKineticsSolve::default()
    };
    let reactor = CellChemistry::new(field.gas().clone(), mechanism, chemistry, || false)?;
    let bytes = ReactingEuler3d::required_workspace_bytes(&field, &reactor)?;
    let mut flow = ReactingEuler3d::with_chemistry(field, reactor, bytes, || false)?;
    let mut molecular = if use_molecular {
        let source = molecular_source(flow.flow().gas())?;
        let bytes = MolecularTransport3d::required_workspace_bytes(flow.flow(), &source)?;
        Some(MolecularTransport3d::new(
            flow.flow(),
            source,
            bytes,
            || false,
        )?)
    } else {
        None
    };
    let bytes = flow
        .workspace_bytes()
        .checked_add(
            molecular
                .as_ref()
                .map_or(0, MolecularTransport3d::workspace_bytes),
        )
        .ok_or("combined numerical workspace size overflow")?;
    let initial = flow.flow().totals();
    let initial_elements = element_amounts(flow.flow().gas(), initial);
    let mut support = [0.0; 3];
    let mut correction = [0.0; 3];
    let mut steps = 0;
    let mut retries = 0;
    let mut cell_solves = 0;
    let mut source_steps = 0;
    let mut maximum_cfl: f64 = 0.0;
    let mut molecular_steps = 0;
    let mut molecular_rejections = 0;
    let mut maximum_parabolic_number: f64 = 0.0;
    let mut maximum_gradient_residual_per_m: f64 = 0.0;
    let mut maximum_mass_frame_residual_kg_m2_s: f64 = 0.0;
    let duration = duration_us * 1e-6;
    while flow.flow().time().value() < duration {
        if steps >= 10_000 {
            return Err("accepted-step budget exhausted".into());
        }
        let remaining = Time::new(duration - flow.flow().time().value());
        let report = match molecular.as_mut() {
            Some(operator) => flow.step_with_dissipation(
                operator,
                remaining,
                ReactiveStepControl::default(),
                DissipationControl::default(),
                bytes,
                || false,
            )?,
            None => flow.step(remaining, ReactiveStepControl::default(), || false)?,
        };
        maximum_cfl = maximum_cfl.max(report.transport.maximum_cfl);
        retries += report.coupling_retries;
        cell_solves += report.chemistry.cell_solves;
        source_steps += report.chemistry.accepted_steps;
        if let Some(split) = report.dissipation {
            for half in [split.first, split.second] {
                molecular_steps += half.accepted_steps;
                molecular_rejections += half.rejected_steps;
                maximum_parabolic_number =
                    maximum_parabolic_number.max(half.maximum_parabolic_number);
                maximum_gradient_residual_per_m =
                    maximum_gradient_residual_per_m.max(half.maximum_gradient_residual_per_m);
                maximum_mass_frame_residual_kg_m2_s = maximum_mass_frame_residual_kg_m2_s
                    .max(half.maximum_mass_frame_residual_kg_m2_s);
            }
        }
        for side in BoundarySide::ALL {
            for impulse in flow.boundary_impulses(side) {
                for axis in 0..3 {
                    let value = impulse.momentum[axis] - correction[axis];
                    let next = support[axis] + value;
                    correction[axis] = (next - support[axis]) - value;
                    support[axis] = next;
                }
            }
        }
        steps += 1;
    }
    let gas = flow.flow().gas();
    let nitrogen = gas
        .source()
        .models()
        .iter()
        .position(|model| model.species().as_str() == "N2")
        .ok_or("pinned source is missing inert N2")?;
    let initial_nitrogen_fraction = initial.species_masses[nitrogen] / initial.mass();
    let after = flow.flow().totals();
    let elements = element_amounts(gas, after);
    let physical_energy_residual = after.energy - initial.energy
        + (0..gas.species_count())
            .map(|k| {
                (after.species_masses[k] - initial.species_masses[k]) * gas.energy_references()[k]
            })
            .sum::<f64>();
    let momentum_residual: [f64; 3] =
        std::array::from_fn(|axis| after.momentum[axis] - initial.momentum[axis] + support[axis]);
    eprintln!(
        "H2/O2/N2 pinned source {SOURCE_COMMIT}; explicit ContinuousAt300K calorics; SDIRK2 source rtol=1e-4; stationary adiabatic impermeable slip walls"
    );
    if use_molecular {
        eprintln!(
            "transport=nonpolar Chapman-Enskog/Wilke/Eucken approximation using pinned source Lennard-Jones data; explicit bulk_viscosity_pa_s=0; translational_cp_deficit_allowance_j_mol_k={CP_FIT_DEFICIT_ALLOWANCE_J_MOL_K:.6e}; polar H2O corrections, rotational relaxation, Soret/Dufour effects excluded"
        );
        eprintln!(
            "accepted_molecular_substeps={molecular_steps}, rejected_molecular_substeps={molecular_rejections}, maximum_parabolic_number={maximum_parabolic_number:.6e}, maximum_gradient_residual_per_m={maximum_gradient_residual_per_m:.6e}, maximum_mass_frame_residual_kg_m2_s={maximum_mass_frame_residual_kg_m2_s:.6e}"
        );
    } else {
        eprintln!("transport=inviscid; no molecular diffusion");
    }
    eprintln!(
        "grid={n}x{n}x{n}, hot_cells={hot_cells}, accepted_steps={steps}, coupling_retries={retries}, accepted_cell_solves={cell_solves}, accepted_source_steps={source_steps}, time_s={:.12e}, maximum_cfl={maximum_cfl:.6}, numerical_workspace_bytes={bytes}",
        flow.flow().time().value()
    );
    eprintln!(
        "mass_residual_kg={:.6e}, physical_energy_residual_j={physical_energy_residual:.6e}, momentum_residual_kg_m_s=[{:.6e},{:.6e},{:.6e}]",
        after.mass() - initial.mass(),
        momentum_residual[0],
        momentum_residual[1],
        momentum_residual[2]
    );
    for (k, element) in gas
        .source()
        .elemental_matrix()
        .elements()
        .iter()
        .enumerate()
    {
        eprintln!(
            "element_{}_residual_mol={:.6e}",
            element.as_str(),
            elements[k] - initial_elements[k]
        );
    }
    let mut out = BufWriter::new(std::io::stdout().lock());
    write!(out, "x_m,y_m,z_m,rho_kg_m3,u_m_s,v_m_s,w_m_s,p_pa,T_K")?;
    for model in gas.source().models() {
        write!(out, ",Y_{}", model.species().as_str())?;
    }
    writeln!(out)?;
    let mut nitrogen_range = [f64::INFINITY, f64::NEG_INFINITY];
    let mut maximum_nitrogen_fraction_change: f64 = 0.0;
    for index in 0..flow.flow().cell_count() {
        let [x, y, z] = grid.cell_center(index).expect("valid cell index");
        let q = flow.flow().cell(index).expect("valid cell index");
        let p = q.primitive(gas, || false)?;
        let nitrogen_fraction = q.partial_densities[nitrogen] / p.density;
        nitrogen_range[0] = nitrogen_range[0].min(nitrogen_fraction);
        nitrogen_range[1] = nitrogen_range[1].max(nitrogen_fraction);
        maximum_nitrogen_fraction_change = maximum_nitrogen_fraction_change
            .max((nitrogen_fraction - initial_nitrogen_fraction).abs());
        write!(
            out,
            "{x:.9},{y:.9},{z:.9},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}",
            p.density, p.velocity[0], p.velocity[1], p.velocity[2], p.pressure, p.temperature
        )?;
        for rho in q.partial_densities {
            write!(out, ",{:.12e}", rho / p.density)?;
        }
        writeln!(out)?;
    }
    out.flush()?;
    // N2 is chemically inert, but differential molecular diffusion can change
    // its local fraction. Only the inviscid mode preserves this uniform ratio.
    let nitrogen_measure = if use_molecular {
        "maximum_inert_N2_fraction_change"
    } else {
        "maximum_inert_N2_fraction_error"
    };
    eprintln!(
        "initial_inert_N2_mass_fraction={initial_nitrogen_fraction:.12e}, final_inert_N2_mass_fraction_range=[{:.12e},{:.12e}], {nitrogen_measure}={maximum_nitrogen_fraction_change:.6e}",
        nitrogen_range[0], nitrogen_range[1],
    );
    Ok(())
}

fn molecular_source(
    gas: &MixtureGas,
) -> Result<SpatialTransport<KineticTransport>, Box<dyn std::error::Error>> {
    // Original diameter [angstrom] and epsilon/kB [K] records, verified at
    // https://github.com/Cantera/cantera/blob/4a8358eb80cfeb50474386b5f9ec0b3a83519889/data/h2o2.yaml
    // Transport uses these customary units regardless of the YAML length unit:
    // https://cantera.org/stable/yaml/species.html#gas-transport
    // Species identities and molecular masses stay with the original gas bank.
    const COLLISION_SOURCE_COMMIT: &str = "4a8358eb80cfeb50474386b5f9ec0b3a83519889";
    if SOURCE_COMMIT != COLLISION_SOURCE_COMMIT {
        return Err("collision records require the retained h2o2 source revision".into());
    }
    let mut species = Vec::with_capacity(10);
    for (name, diameter_angstrom, well_depth_k) in [
        ("H2", 2.92, 38.0),
        ("H", 2.05, 145.0),
        ("O", 2.75, 80.0),
        ("O2", 3.458, 107.4),
        ("OH", 2.75, 80.0),
        ("H2O", 2.605, 572.4),
        ("HO2", 3.458, 107.4),
        ("H2O2", 3.458, 107.4),
        ("AR", 3.33, 136.5),
        ("N2", 3.621, 97.53),
    ] {
        species.push(LennardJonesSpeciesV1::new(
            SpeciesId::new(name)?,
            Length::new(diameter_angstrom * 1e-10),
            Temperature::new(well_depth_k),
        )?);
    }
    let model = KineticGasTransportV1::new(
        gas.source().clone(),
        format!("cantera-h2o2-{COLLISION_SOURCE_COMMIT}-nonpolar-eucken"),
        1,
        species,
        KineticTransportDomainV1 {
            minimum_temperature_k: 300.0,
            maximum_temperature_k: 3500.0,
            minimum_pressure_pa: 1e3,
            maximum_pressure_pa: 1e6,
        },
        || false,
    )?
    .with_translational_cp_deficit_allowance(MolarThermalQuantityV1::new(
        CP_FIT_DEFICIT_ALLOWANCE_J_MOL_K,
    ))?;
    Ok(SpatialTransport::new(
        gas.clone(),
        KineticTransport::new(model, DynViscosity::new(0.0))?,
        || false,
    )?)
}

fn element_amounts(gas: &MixtureGas, totals: Totals) -> Vec<f64> {
    let source = gas.source();
    let matrix = source.elemental_matrix();
    (0..matrix.row_count())
        .map(|element| {
            (0..gas.species_count())
                .map(|species| {
                    totals.species_masses[species] / source.models()[species].molar_mass().value()
                        * matrix
                            .get(element, species)
                            .expect("original canonical source axes")
                            as f64
                })
                .sum()
        })
        .collect()
}
