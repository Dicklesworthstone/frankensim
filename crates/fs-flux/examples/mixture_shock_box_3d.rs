//! Run `cargo run --release -p fs-flux --no-default-features --features experimental-compressible --example mixture_shock_box_3d -- 12`.
//! A hot nitrogen-rich pocket drives a 3-D shock into an oxygen-rich mixture.
//! CSV contains the evolved composition and thermodynamic field. This fixed-box
//! inviscid example transports species; it does not integrate chemical reactions.

use fs_flux::compressible::spatial::mixture::{MixtureEuler3d, MixtureGas, state_into};
use fs_flux::compressible::spatial::{
    AxisBoundaries, Boundaries3d, BoundaryCondition, BoundarySide, CartesianGrid, StepControl,
};
use fs_gas::multicomponent::CaloricDomain;
use fs_qty::{Temperature, Time};
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::kinetics::hydrogen::{HydrogenCalorics, SOURCE_COMMIT, cantera_h2o2};
use fs_thermochem::{
    ChargeVector, ElementId, ElementalMatrix, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
};
use std::io::{BufWriter, Write};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let n = args
        .next()
        .map(|s| s.parse::<usize>())
        .transpose()?
        .unwrap_or(12);
    if !(8..=32).contains(&n) || args.next().is_some() {
        return Err("usage: mixture_shock_box_3d [cells_per_axis in 8..=32]".into());
    }
    let gas = nitrogen_oxygen()?;
    let grid = CartesianGrid::new(n, n, n, 1.0 / n as f64, 1.0 / n as f64, 1.0 / n as f64)?;
    let molar_mass = [
        gas.source().models()[0].molar_mass().value(),
        gas.source().models()[1].molar_mass().value(),
    ];
    let pressure = |density: [f64; 2], temperature: f64| {
        R * temperature * (density[0] / molar_mass[0] + density[1] / molar_mass[1])
    };
    let ambient_density = [0.25, 0.75];
    let ambient_pressure = pressure(ambient_density, 450.0);
    let ambient = gas.primitive(
        gas.conserved(&ambient_density, 0.0, ambient_pressure)?
            .as_cell(),
        || false,
    )?;
    let duration = 0.7 / ambient.sound_speed;
    let mut states = vec![0.0; n * n * n * 6];
    for (index, row) in states.chunks_exact_mut(6).enumerate() {
        let [x, y, z] = grid.cell_center(index).expect("valid cell index");
        let hot = (x - 0.32).powi(2) + (y - 0.44).powi(2) + (z - 0.57).powi(2) < 0.14_f64.powi(2);
        let (density, temperature) = if hot {
            ([0.95, 0.05], 850.0)
        } else {
            (ambient_density, 450.0)
        };
        state_into(
            &gas,
            &density,
            [0.0; 3],
            pressure(density, temperature),
            row,
        )?;
    }
    let walls = AxisBoundaries::Ends {
        lower: BoundaryCondition::Reflecting,
        upper: BoundaryCondition::Reflecting,
    };
    let bytes = MixtureEuler3d::required_workspace_bytes(&gas, &grid)?;
    let mut flow = MixtureEuler3d::new(
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
    let initial = flow.totals();
    let mut support_impulse = [0.0; 3];
    let mut support_correction = [0.0; 3];
    let mut steps = 0;
    let mut maximum_cfl: f64 = 0.0;
    while flow.time().value() < duration {
        if steps >= 20_000 {
            return Err("accepted-step budget exhausted".into());
        }
        let report = flow.step(
            Time::new(duration - flow.time().value()),
            StepControl::default(),
            || false,
        )?;
        maximum_cfl = maximum_cfl.max(report.maximum_cfl);
        for side in BoundarySide::ALL {
            for impulse in flow.boundary_impulses(side) {
                for axis in 0..3 {
                    let value = impulse.momentum[axis] - support_correction[axis];
                    let next = support_impulse[axis] + value;
                    support_correction[axis] = (next - support_impulse[axis]) - value;
                    support_impulse[axis] = next;
                }
            }
        }
        steps += 1;
    }
    let final_state = flow.totals();
    let momentum_residual: [f64; 3] = std::array::from_fn(|axis| {
        final_state.momentum[axis] - initial.momentum[axis] + support_impulse[axis]
    });
    eprintln!(
        "N2/O2 thermodynamics from pinned source {SOURCE_COMMIT}; fixed-box inviscid mixture, no reactions or diffusion"
    );
    eprintln!(
        "grid={n}x{n}x{n}, accepted_steps={steps}, time_s={:.12e}, maximum_cfl={maximum_cfl:.6}, numerical_workspace_bytes={bytes}",
        flow.time().value(),
    );
    eprintln!(
        "species_mass_residual_kg=[{:.6e},{:.6e}], referenced_energy_residual_j={:.6e}, physical_energy_residual_j={:.6e}, momentum_residual_kg_m_s=[{:.6e},{:.6e},{:.6e}]",
        final_state.species_masses[0] - initial.species_masses[0],
        final_state.species_masses[1] - initial.species_masses[1],
        final_state.energy - initial.energy,
        final_state.absolute_energy(flow.gas())? - initial.absolute_energy(flow.gas())?,
        momentum_residual[0],
        momentum_residual[1],
        momentum_residual[2],
    );
    let mut out = BufWriter::new(std::io::stdout().lock());
    writeln!(
        out,
        "x_m,y_m,z_m,rho_kg_m3,u_m_s,v_m_s,w_m_s,p_pa,T_K,Y_N2,Y_O2"
    )?;
    for index in 0..flow.cell_count() {
        let [x, y, z] = grid.cell_center(index).expect("valid cell index");
        let q = flow.cell(index).expect("valid cell index");
        let p = q.primitive(flow.gas(), || false)?;
        writeln!(
            out,
            "{x:.9},{y:.9},{z:.9},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}",
            p.density,
            p.velocity[0],
            p.velocity[1],
            p.velocity[2],
            p.pressure,
            p.temperature,
            q.partial_densities[0] / p.density,
            q.partial_densities[1] / p.density,
        )?;
    }
    out.flush()?;
    Ok(())
}

fn nitrogen_oxygen() -> Result<MixtureGas, Box<dyn std::error::Error>> {
    // Reuse the retained source models without copying or refitting their NASA
    // coefficients. Selecting inert N2/O2 does not activate the source kinetics.
    let source = cantera_h2o2(HydrogenCalorics::Original)?;
    let models = ["N2", "O2"]
        .into_iter()
        .map(|name| {
            source
                .thermodynamics()
                .models()
                .iter()
                .find(|model| model.species().as_str() == name)
                .cloned()
                .ok_or("pinned source is missing the selected species")
        })
        .collect::<Result<Vec<_>, _>>()?;
    let ids: Vec<_> = models.iter().map(|model| model.species().clone()).collect();
    let bank = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("N")?, ElementId::new("O")?],
            ids.clone(),
            vec![vec![2, 0], vec![0, 2]],
        )?,
        ChargeVector::new(ids, vec![0, 0])?,
        models,
    )?;
    // Remain in the original low-temperature regions; do not implicitly bridge
    // independently fitted NASA region seams or extrapolate shocked states.
    Ok(MixtureGas::new(
        bank,
        CaloricDomain::new(Temperature::new(300.0), Temperature::new(990.0)),
        || false,
    )?)
}
