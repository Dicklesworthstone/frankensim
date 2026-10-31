//! Synthetic material swaps through the original coupled circuit/thermal/gas solver.
//! The fixtures below are existing manufactured sources, not measured materials.
#[path = "../../tests/finite_gas_convection/source.rs"]
mod gas_source;
#[path = "../../../fs-couple/tests/circuit_enthalpy/common.rs"]
#[allow(dead_code)]
mod common;

use fs_airflow::graph::thermal::medium::convection::{ForcedConvectionSettings,
    finite_gas::{ForcedGasSurface, natural::NaturalGasSurface,
        network::{ConvectiveContact, exchange::{self, ExchangeContact, GasHeatLaw}}}};
#[cfg(test)]
use fs_airflow::graph::thermal::medium::convection::finite_gas::network;
use fs_convection::CorrelationId;
use fs_couple::cavity_phs::piston::{GasPiston, PistonInput, PistonSolve};
use fs_couple::cavity_phs::piston::thermal::{HeatExchangeSolve, source::ThermalMaterialPoint};
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::{Chamber, Contact, NetworkSolve};
use fs_couple::circuit_enthalpy::material::{CellResistor, ConductionLink, MaterialCell, MaterialNetwork};
use fs_couple::material_circuit::MaterialCircuit;
use fs_couple::circuit_enthalpy::contacts::material::ElectricalContactNetwork;
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::material::exchange::MaterialExchangeNetwork;
use fs_material::conductor::ElectricalPropertyBasis;
use fs_material::state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermochemical_gas::energy::{GasEnergySolve, parcel::ClosedGasParcel};
use fs_matdb::PropertyKey;
use fs_qty::Length;
use fs_thermal::enthalpy::{CellId, MaterialId, DepositionSimulation};
use std::collections::BTreeMap;
use std::error::Error;

type Result<T> = std::result::Result<T, Box<dyn Error>>;
const DT: f64 = 0.00005;
const BASES: [(usize, ElectricalPropertyBasis); 2] = [
    (1, ElectricalPropertyBasis::Resistivity), (2, ElectricalPropertyBasis::Resistivity),
];
fn fixture<R>(run: impl FnOnce(&MaterialNetwork<'_>, DepositionSimulation<'_>, MaterialCircuit) -> R) -> R {
    let (a, ca) = common::material_with_conductivity(1.0, 1.0, 0.005, 0.5);
    let (b, cb) = common::material_with_conductivity(3.0, 2.0, 0.002, 2.0);
    let cards = [a, b]; let curves = [ca, cb]; let context = common::point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = [0, 1].map(|i| MaterialCell {
        cell: common::DESTINATIONS[i].cell, material_id: MaterialId(i as u64), curve: curves[i].curve(),
        point: ThermalMaterialPoint { material: &cards[i], context: &context,
            temperature_axis: "T", liquid_fraction_axis: None },
        conductivity: Some((&key, MaterialPropertySelection::SingleClaimOnly)),
    });
    let resistors = [0, 1].map(|i| CellResistor { destination: common::DESTINATIONS[i],
        selection: MaterialPropertySelection::SingleClaimOnly, length_m: 1.0, area_m2: 1.0 });
    let links = [ConductionLink { link: 7, distance_a_m: 0.5, distance_b_m: 0.5 }];
    let thermal = common::initial(&curves, 999.0, None);
    let circuit = common::circuit(&thermal, &cards, &curves);
    run(&MaterialNetwork { cells: &cells, resistors: &resistors, links: &links }, thermal, circuit)
}
fn gases(fraction: f64) -> BTreeMap<u64, GasPiston> {
    [(10, fraction), (20, 0.5)].into_iter().map(|(id, x)| (id,
        GasPiston::new(ClosedGasParcel::new(gas_source::model(x), gas_source::point(400.0, 100000.0),
            0.001).unwrap(), 1.0, 0.01, 0.0).unwrap())).collect()
}
fn chambers(gases: &BTreeMap<u64, GasPiston>, dt: f64) -> Vec<Chamber<'_>> {
    gases.iter().map(|(&id, piston)| Chamber { id, piston,
        input: PistonInput { duration_s: dt, heat_j: 0.0, exterior_pressure_pa: 100000.0, external_force_n: 0.0 },
        piston_solve: PistonSolve { volume_bracket_m3: [0.0009, 0.0011],
            gas: GasEnergySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 900.0,
                energy_tolerance_j_kg: 1e-7, maximum_iterations: 128 },
            momentum_tolerance_n_s: 1e-11, energy_tolerance_j: 1e-7, maximum_iterations: 128 },
        exchange: HeatExchangeSolve { heat_bracket_j: [-0.1, 0.1], source_tolerance_j: 1e-10,
            energy_tolerance_j: 1e-7, maximum_iterations: 128 },
    }).collect()
}
fn contacts(speed: f64) -> [ConvectiveContact; 3] {
    [(101, 10, CellId(30)), (202, 10, CellId(10)), (303, 20, CellId(30))].map(|(id, chamber, cell)|
        ConvectiveContact { contact: Contact { id, chamber, cell }, surface: ForcedGasSurface {
            area_m2: 0.1, characteristic_length: Length::new(0.1), speed_m_s: speed,
            correlation: CorrelationId::FlatPlateLaminarAverage, settings: ForcedConvectionSettings::default(),
        } })
}
fn controls() -> NetworkSolve {
    NetworkSolve { network: common::config(), maximum_chambers: 2, maximum_contacts: 3,
        maximum_boundary_solves: 2000, total_energy_tolerance_j: 4e-6 }
}
fn storage(thermal: &DepositionSimulation<'_>, circuit: &MaterialCircuit, gas: &BTreeMap<u64, GasPiston>) -> f64 {
    thermal.energy_j() + 0.5*circuit.charges_c()[0].powi(2) + gas.values().map(|g|
        g.gas().mass_kg()*g.gas().state().specific_internal_energy_j_kg()
        + 0.5*g.mass_kg()*g.velocity_m_s().powi(2)).sum::<f64>()
}
#[derive(Debug)]
struct Outcome { cell_k: f64, charge_c: f64, velocity_m_s: f64, sound_m_s: f64, energy_defect_j: f64 }
fn with_runtime<R>(run: impl FnOnce(&fs_exec::Cx<'_>) -> R) -> R {
    let gate = fs_exec::CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| run(&fs_exec::Cx::new(&gate, arena,
        fs_exec::StreamKey { seed: 71, kernel_id: 3, tile: 0, iteration: 0 },
        fs_exec::Budget::INFINITE, fs_exec::ExecMode::Deterministic)))
}
fn heat_law(input: &str) -> Result<GasHeatLaw> {
    let geometry = contacts(1.0)[0].surface;
    if input == "natural" || input.starts_with("natural:") {
        let gravity: f64 = match input.strip_prefix("natural:") {
            Some(value) => value.parse()?, None => 9.80665,
        };
        if !gravity.is_finite() || gravity <= 0.0 {
            return Err("natural-convection gravity must be positive and finite".into());
        }
        Ok(GasHeatLaw::NaturalConvection(NaturalGasSurface {
            geometry: fs_convection::NaturalConvectionGeometry::VerticalPlate,
            area_m2: geometry.area_m2, characteristic_length: geometry.characteristic_length,
            gravity_m_s2: gravity,
        }))
    } else {
        let speed: f64 = input.parse()?;
        if !speed.is_finite() || speed <= 0.0 {
            return Err("forced-convection speed must be positive and finite".into());
        }
        Ok(GasHeatLaw::ForcedConvection(ForcedGasSurface { speed_m_s: speed, ..geometry }))
    }
}
#[cfg(test)]
fn run(fraction: f64, speed: f64, steps: usize, emit: bool) -> Result<Outcome> {
    run_with_law(fraction, GasHeatLaw::ForcedConvection(contacts(speed)[0].surface), steps, emit)
}
fn run_with_law(fraction: f64, law: GasHeatLaw, steps: usize, emit: bool) -> Result<Outcome> {
    with_runtime(|cx| fixture(|binding, mut thermal, mut circuit| {
        let electrical = ElectricalContactNetwork { bulk: binding, contacts: &[], electrical_bases: &BASES };
        let assembly = MaterialExchangeNetwork { electrical: &electrical, thermal_contacts: &[], radiation: None };
        let prepared = assembly.prepare(cx, &thermal)?;
        let mut gas = gases(fraction);
        let initial = storage(&thermal, &circuit, &gas);
        let patches = contacts(1.0).map(|patch| ExchangeContact { contact: patch.contact, law });
        let mut outside = 0.0;
        if emit { println!("time_s,cell30_k,cell10_k,charge_c,gas10_k,gas10_pa,wall10_m_s,sound10_m_s,energy_defect_j"); }
        for i in 0..steps {
            let target = (i+1) as f64*DT;
            let drives = chambers(&gas, target-thermal.time_s());
            let candidate = exchange::step_candidate(&prepared, cx, &thermal, target, &circuit, &[],
                &drives, &[], &patches, controls(), &[0.0; 2])?;
            let (next_thermal, next_circuit, report, _) = candidate.into_parts();
            let mut next_gas = BTreeMap::new();
            for entry in report.source.chambers {
                outside += entry.boundary.piston.exterior_pressure_work_j;
                next_gas.insert(entry.id, entry.piston);
            }
            thermal = next_thermal; circuit = next_circuit; gas = next_gas;
            if emit {
                let g = &gas[&10]; let medium = g.gas().state().state();
                println!("{target:.9},{:.12},{:.12},{:.12},{:.12},{:.12},{:.12},{:.12},{:.12}",
                    thermal.temperature_k(CellId(30))?, thermal.temperature_k(CellId(10))?,
                    circuit.charges_c()[0], medium.temperature, medium.pressure, g.velocity_m_s(),
                    medium.sound_speed, storage(&thermal, &circuit, &gas)+outside-initial);
            }
        }
        Ok(Outcome { cell_k: thermal.temperature_k(CellId(30))?, charge_c: circuit.charges_c()[0],
            velocity_m_s: gas[&10].velocity_m_s(), sound_m_s: gas[&10].gas().state().state().sound_speed,
            energy_defect_j: storage(&thermal, &circuit, &gas)+outside-initial })
    }))
}
fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.is_empty() || args.len() > 3 {
        return Err("usage: material_convective_network A_mole_fraction [speed_m_s|natural[:gravity_m_s2]] [steps=24]".into());
    }
    let fraction: f64 = args[0].parse()?;
    let law = heat_law(args.get(1).map_or("1", String::as_str))?;
    let steps: usize = args.get(2).map_or(Ok(24), |s| s.parse())?;
    if !fraction.is_finite() || !(0.01..=0.99).contains(&fraction) || !(1..=500).contains(&steps) {
        return Err("fraction must be 0.01..=0.99 and steps 1..=500".into());
    }
    match law {
        GasHeatLaw::ForcedConvection(surface) =>
            eprintln!("forced convection: relative speed {} m/s", surface.speed_m_s),
        GasHeatLaw::NaturalConvection(surface) =>
            eprintln!("natural convection: gravity {} m/s²", surface.gravity_m_s2),
        GasHeatLaw::Conduction(_) => eprintln!("conductive gas boundary"),
    }
    let out = run_with_law(fraction, law, steps, true)?;
    eprintln!("final cell={} K, charge={} C, velocity={} m/s, sound={} m/s, balance={} J",
        out.cell_k, out.charge_c, out.velocity_m_s, out.sound_m_s, out.energy_defect_j);
    Ok(())
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod exchange_tests;
