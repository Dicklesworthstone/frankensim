//! Source-loaded finite gas, existing slider-crank motion and actual shaft work.
//! See driven_gas_cycle/README.md for the explicit demonstration geometry and
//! source compilation. CSV is written only after the complete trajectory passes.

use fs_couple::prescribed_gas::{PrescribedGasDrive, PrescribedGasSolve, step_candidate};
use fs_matdb::{
    MaterialStateId, NormalizedMaterialCardPack, NormalizedPack, NormalizedSpeciesPack,
};
use fs_material::gas_species::input::ThermochemicalPhase;
use fs_material::gas_species::input::packs::{SourcePackBytes, read_source_parcel};
use fs_material::thermochemical_gas::energy::GasEnergySolve;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::GasNozzle;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::network::snapshot::SnapshotNozzleLink;
use fs_material::thermochemical_gas::energy::parcel::transfer::snapshot::GasSnapshotSolve;
use fs_mbd::otto::{OttoTopologyParams, step_otto_topology};
use std::{
    collections::BTreeMap,
    io::{Read, Write},
    path::Path,
};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

// Authored SI engineering geometry, not historical Otto dimensions. The short
// stroke keeps the supplied journey inside the existing dilute transport cards.
const CRANK_M: f64 = 0.001;
const ROD_M: f64 = 0.004;
const BORE_M: f64 = 0.05;
const CLEARANCE_M3: f64 = 0.00008;
const RECEIVER_M3: f64 = 0.002;
const RPM: f64 = 120.0;
const EXTERIOR_PA: f64 = 80000.0;
const NOZZLE_M2: f64 = 5e-8;

fn read(path: &Path, maximum: usize) -> Result<Vec<u8>> {
    let mut bytes = Vec::new();
    std::fs::File::open(path)?
        .take((maximum + 1) as u64)
        .read_to_end(&mut bytes)?;
    if bytes.len() > maximum {
        return Err(format!("{} exceeds {maximum} bytes", path.display()).into());
    }
    Ok(bytes)
}

fn volume(angle: f64) -> Result<f64> {
    let pose = step_otto_topology(OttoTopologyParams {
        crank_angle_rad: angle,
        crank_radius: CRANK_M,
        connecting_rod_length: ROD_M,
        engine_rpm: RPM,
    })
    .map_err(|error| format!("original slider-crank refused: {error:?}"))?;
    let piston_tdc = -CRANK_M - ROD_M;
    Ok(CLEARANCE_M3
        + core::f64::consts::PI * BORE_M * BORE_M / 4.0 * (pose.piston_pin_x - piston_tdc))
}

fn extensive_change(old: &ClosedGasParcel, new: &ClosedGasParcel) -> f64 {
    new.mass_kg()
        * (new.state().specific_internal_energy_j_kg()
            - old.state().specific_internal_energy_j_kg())
        + (new.mass_kg() - old.mass_kg()) * old.state().specific_internal_energy_j_kg()
}

struct Row {
    angle: f64,
    compression: bool,
    volume: f64,
    pressure: f64,
    temperature: f64,
    mass: f64,
    receiver_pressure: f64,
    receiver_temperature: f64,
    receiver_mass: f64,
    drive_work: f64,
    torque: f64,
    cumulative_work: f64,
    residual: f64,
}

fn simulate(mut gases: BTreeMap<u64, ClosedGasParcel>, count: usize) -> Result<Vec<Row>> {
    if !(8..=4096).contains(&count) {
        return Err("steps per stroke must lie in 8..=4096".into());
    }
    if gases[&1].state().state().pressure <= gases[&2].state().state().pressure {
        return Err("this discharge-stroke example requires initial cylinder pressure above receiver pressure".into());
    }
    let original = gases.clone();
    let gas_solve = GasEnergySolve {
        minimum_temperature_k: 200.0,
        maximum_temperature_k: 900.0,
        energy_tolerance_j_kg: 1e-7,
        maximum_iterations: 128,
    };
    let endpoints = BTreeMap::from([(1, gas_solve), (2, gas_solve)]);
    let controls = PrescribedGasSolve {
        exchange: GasSnapshotSolve {
            maximum_withdrawn_fraction: 0.05,
            mass_tolerance_kg: 1e-15,
            energy_tolerance_j: 1e-8,
        },
        total_energy_tolerance_j: 1e-7,
    };
    let mut rows = Vec::with_capacity(2 * count);
    let (mut work, mut exterior) = (0.0, 0.0);
    let mut previous = core::f64::consts::PI;
    for index in 1..=2 * count {
        let angle = core::f64::consts::PI + core::f64::consts::PI * (index as f64 / count as f64);
        let delta = angle - previous;
        let dt = delta / (RPM * core::f64::consts::TAU / 60.0);
        let compression = index <= count;
        let links = [SnapshotNozzleLink {
            id: 1,
            parcels: [1, 2],
            nozzle: GasNozzle {
                area_m2: if compression { NOZZLE_M2 } else { 0.0 },
                discharge_coefficient: 0.8,
            },
            nozzle_solve: gas_solve,
        }];
        let drives = BTreeMap::from([
            (
                1,
                PrescribedGasDrive {
                    final_volume_m3: volume(angle)?,
                    heat_j: 0.0,
                    exterior_pressure_pa: EXTERIOR_PA,
                },
            ),
            (
                2,
                PrescribedGasDrive {
                    final_volume_m3: RECEIVER_M3,
                    heat_j: 0.0,
                    exterior_pressure_pa: EXTERIOR_PA,
                },
            ),
        ]);
        let candidate=step_candidate(&gases,&links,dt,&endpoints,&drives,controls,||false)
            .map_err(|error|format!("crank interval {index}: {error}; refine time or select an admitted source/geometry"))?;
        let (next, step) = candidate.into_parts();
        work += step.net_drive_work_j;
        exterior += step.exterior_pressure_work_j;
        let residual = extensive_change(&original[&1], &next[&1])
            + extensive_change(&original[&2], &next[&2])
            + work
            + exterior;
        if !residual.is_finite() || residual.abs() > 1e-6 {
            return Err(
                format!("cumulative physical energy residual {residual} J exceeds 1e-6 J").into(),
            );
        }
        rows.push(Row {
            angle,
            compression,
            volume: next[&1].volume_m3(),
            pressure: next[&1].state().state().pressure,
            temperature: next[&1].state().state().temperature,
            mass: next[&1].mass_kg(),
            receiver_pressure: next[&2].state().state().pressure,
            receiver_temperature: next[&2].state().state().temperature,
            receiver_mass: next[&2].mass_kg(),
            drive_work: step.net_drive_work_j,
            torque: step.net_drive_work_j / delta,
            cumulative_work: work,
            residual,
        });
        gases = next;
        previous = angle;
    }
    eprintln!(
        "net work to prescribed drive: {work:.12e} J; mass delivered to finite receiver: {:.12e} kg; final storage balance: {:.12e} J",
        gases[&2].mass_kg() - original[&2].mass_kg(),
        rows.last().expect("nonempty trajectory").residual
    );
    Ok(rows)
}

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if !(5..=6).contains(&args.len()) {
        return Err("usage: driven_gas_cycle SPECIES_PACK NASA9_PACK TRANSPORT_PACK CYLINDER_RECORD RECEIVER_RECORD [STEPS_PER_STROKE]".into());
    }
    let count = if args.len() == 6 {
        args[5].parse::<usize>()?
    } else {
        256
    };
    let species = read(Path::new(&args[0]), 16 << 20)?;
    let nasa9 = read(Path::new(&args[1]), 16 << 20)?;
    let transport = read(Path::new(&args[2]), 16 << 20)?;
    // matdb-pack's pure-property sources produce a claim pack. Bind that
    // original pack to the explicit species, as the existing gas CLI does.
    // Already canonical material cards keep their exact original identity.
    let transport = if transport.starts_with(b"FSMCDPK\0") {
        transport
    } else {
        let identity = NormalizedSpeciesPack::from_bytes(&species)?;
        NormalizedMaterialCardPack::new(
            MaterialStateId {
                chemistry: identity.association().species().as_str().into(),
                phase: identity.association().standard_state_phase().into(),
                process: "caller-declared gas-component source association".into(),
                revision: 0,
            },
            NormalizedPack::from_bytes(&transport)?,
        )?
        .to_bytes()
    };
    let sources = [SourcePackBytes {
        species: &species,
        nasa9: &nasa9,
        transport: &transport,
    }];
    let load = |record: &str, volume| -> Result<ClosedGasParcel> {
        let bytes = read(Path::new(record), 8192)?;
        Ok(ClosedGasParcel::from_isochoric(read_source_parcel(
            std::str::from_utf8(&bytes)?,
            &sources,
            None,
            ThermochemicalPhase::FrozenGas,
            volume,
        )?))
    };
    let gases = BTreeMap::from([
        (1, load(&args[3], volume(core::f64::consts::PI)?)?),
        (2, load(&args[4], RECEIVER_M3)?),
    ]);
    let rows = simulate(gases, count)?;
    let mut out = std::io::BufWriter::new(std::io::stdout().lock());
    writeln!(
        out,
        "crank_angle_rad,stroke,volume_m3,cylinder_pressure_pa,cylinder_temperature_k,cylinder_mass_kg,receiver_pressure_pa,receiver_temperature_k,receiver_mass_kg,drive_work_j,average_drive_torque_nm,cumulative_drive_work_j,total_energy_residual_j"
    )?;
    for row in rows {
        writeln!(
            out,
            "{:.12e},{},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}",
            row.angle,
            if row.compression {
                "compression_discharge"
            } else {
                "closed_expansion"
            },
            row.volume,
            row.pressure,
            row.temperature,
            row.mass,
            row.receiver_pressure,
            row.receiver_temperature,
            row.receiver_mass,
            row.drive_work,
            row.torque,
            row.cumulative_work,
            row.residual
        )?;
    }
    out.flush()?;
    Ok(())
}
