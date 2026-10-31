//! A finite real-gas reservoir, heated body and moving chamber in one trajectory.
//!
//! Run: cargo run -p fs-airflow --example real_gas_heated_chamber -- 5000000 500 200
//! Optional arguments: reservoir pressure [Pa], body temperature [K], steps.
//! The chamber starts at 4.8 MPa and 300 K; lower reservoir pressure reverses
//! the stream. A 2 kg piston exposes 0.01 m2, with 20 us physical timesteps.
//!
//! The native nozzle transfers mass/enthalpy once, then the native reciprocal
//! thermal/motion solve queries PR film properties at the actual mean gas/body
//! temperature and chamber pressure. All three endpoints advance together.
//! The split is first order and requires timestep refinement. Synthetic data,
//! homogeneous single-phase gases, an isothermal body and prescribed relative
//! flow are explicit assumptions; no jet thrust or fan-energy model is added.

use fs_airflow::graph::thermal::medium::convection::{
    ForcedConvectionSettings,
    finite_gas::{
        ForcedGasSurface,
        film::{FilmGasHeatLaw, FilmGasSurface, PengRobinsonFilmSource},
    },
};
use fs_convection::CorrelationId;
use fs_couple::cavity_phs::piston::{
    GasPiston, PistonInput, PistonSolve, thermal::HeatExchangeSolve,
};
use fs_matdb::{ClaimSet, MaterialCard, MaterialStateId};
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase, UniformEnthalpyStepInput,
};
use fs_material::thermochemical_gas::{
    energy::{GasEnergySolve, entropy::GasEntropySolve, parcel::transfer::nozzle::GasNozzle},
    nonideal::{
        PengRobinsonPhaseV1, PengRobinsonSpeciesBinding, ThermochemicalPengRobinsonModel,
        energy::parcel::{
            IsochoricRealGasParcel,
            transfer::{RealGasMassTransferSolve, nozzle::RealGasNozzleStepSolve},
        },
    },
};
use fs_qty::{Length, Pressure, Temperature};
use fs_thermochem::SpeciesId;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
#[path = "../../fs-material/examples/material_gas_spring/source.rs"]
mod source;

const DT: f64 = 2e-5;
const INITIAL_VOLUME: f64 = 0.001;
const EXTERIOR_PRESSURE: f64 = 4.8e6;
const BODY_MASS: f64 = 0.01;
const BODY_DENSITY: f64 = 4000.0;
const SURFACE_AREA: f64 = 0.02;

fn model() -> Result<ThermochemicalPengRobinsonModel> {
    let source = source::load(
        "A 0.75 0.028 3.5 0.001 1.7e-5 0.025\nB 0.25 0.004 2.5 0 2e-5 0.08\n",
        "real_gas_heated_chamber manufactured mixture",
    )?;
    Ok(ThermochemicalPengRobinsonModel::new(
        source,
        vec![
            PengRobinsonSpeciesBinding {
                species: SpeciesId::new("A")?,
                critical_temperature: Temperature::new(190.564),
                critical_pressure: Pressure::new(4.5992e6),
                acentric_factor: 0.01142,
            },
            PengRobinsonSpeciesBinding {
                species: SpeciesId::new("B")?,
                critical_temperature: Temperature::new(150.687),
                critical_pressure: Pressure::new(4.863e6),
                acentric_factor: -0.00219,
            },
        ],
        vec![0.0, 0.03, 0.03, 0.0],
    )?)
}

fn body_curve() -> Result<EquilibriumEnthalpyPhaseCurve> {
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic thermal body".into(),
            phase: "solid".into(),
            process: "manufactured example".into(),
            revision: 0,
        },
        ClaimSet::new(),
        vec![],
    )?;
    Ok(EquilibriumEnthalpyPhaseCurve::try_single_phase(
        card.content_hash(),
        SolidLiquidPhase::Solid,
        [200.0, 800.0]
            .map(|temperature_k| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 1000.0 * (temperature_k - 200.0),
                temperature_k,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: BODY_DENSITY,
            })
            .to_vec(),
    )?)
}

fn relative_energy(
    gas: &IsochoricRealGasParcel,
    reference: &IsochoricRealGasParcel,
) -> Result<f64> {
    // Both inventories retain the same composition and original law. Their
    // conserved total mass cancels the common formation-energy reference.
    Ok(gas.mass_kg()
        * gas
            .model()
            .specific_internal_energy_change_between_densities(
                reference.density_kg_m3(),
                reference.state().temperature_k(),
                gas.density_kg_m3(),
                gas.state().temperature_k(),
            )?)
}

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() > 3 {
        return Err(
            "usage: real_gas_heated_chamber [reservoir_pressure_pa] [body_temperature_k] [steps]"
                .into(),
        );
    }
    let reservoir_pressure: f64 = args.first().map_or(Ok(5e6), |v| v.parse())?;
    let temperature: f64 = args.get(1).map_or(Ok(500.0), |v| v.parse())?;
    let steps: usize = args.get(2).map_or(Ok(200), |v| v.parse())?;
    if !reservoir_pressure.is_finite()
        || reservoir_pressure <= 0.0
        || !temperature.is_finite()
        || !(225.0..=775.0).contains(&temperature)
        || !(1..=10_000).contains(&steps)
    {
        return Err(
            "positive finite pressure, body temperature 225..775 K and 1..10000 steps required"
                .into(),
        );
    }
    let model = model()?;
    let chamber = IsochoricRealGasParcel::new(
        model.clone(),
        source::point(300.0, EXTERIOR_PRESSURE)?,
        INITIAL_VOLUME,
        PengRobinsonPhaseV1::Vapor,
    )?;
    let reference = chamber.clone();
    let mut reservoir = IsochoricRealGasParcel::new(
        model.clone(),
        source::point(300.0, reservoir_pressure)?,
        0.002,
        PengRobinsonPhaseV1::Vapor,
    )?;
    let initial_pair_mass = chamber.mass_kg() + reservoir.mass_kg();
    let initial_relative_energy = relative_energy(&reservoir, &reference)?;
    let mut piston = GasPiston::new(chamber, 2.0, 0.01, 0.0)?;
    let curve = body_curve()?;
    let mut thermal = curve.state_at_specific_enthalpy(1000.0 * (temperature - 200.0))?;
    let initial_body_enthalpy = thermal.specific_enthalpy_j_kg();
    let energy = GasEnergySolve {
        minimum_temperature_k: 250.0,
        maximum_temperature_k: 500.0,
        energy_tolerance_j_kg: 1e-7,
        maximum_iterations: 128,
    };
    let flow = RealGasNozzleStepSolve {
        entropy: GasEntropySolve {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 500.0,
            entropy_tolerance_j_kg_k: 1e-9,
            maximum_iterations: 128,
        },
        transfer: RealGasMassTransferSolve {
            donor: energy,
            receiver: energy,
            mass_tolerance_kg: 1e-14,
            energy_tolerance_j: 1e-7,
        },
        minimum_density_fraction: 0.5,
        pressure_tolerance_pa: 1e-4,
        sonic_tolerance_m2_s2: 1e-4,
        maximum_iterations: 128,
        maximum_donor_fraction: 0.01,
    };
    let exchange = HeatExchangeSolve {
        heat_bracket_j: [-0.25, 0.25],
        source_tolerance_j: 1e-10,
        energy_tolerance_j: 2e-6,
        maximum_iterations: 128,
    };
    let property_source = PengRobinsonFilmSource {
        model: &model,
        phase: PengRobinsonPhaseV1::Vapor,
    };
    let surface = FilmGasSurface {
        model: &property_source,
        law: FilmGasHeatLaw::Forced(ForcedGasSurface {
            area_m2: SURFACE_AREA,
            characteristic_length: Length::new(0.1),
            speed_m_s: 0.5,
            correlation: CorrelationId::FlatPlateLaminarAverage,
            settings: ForcedConvectionSettings::default(),
        }),
    };
    let nozzle = GasNozzle {
        area_m2: 2e-7,
        discharge_coefficient: 0.8,
    };
    let mut wall_heat = 0.0;
    let mut gas_work = 0.0;
    println!(
        "step,time_s,chamber_temperature_k,chamber_pressure_pa,chamber_density_kg_m3,chamber_mass_kg,chamber_volume_m3,piston_velocity_m_s,reservoir_temperature_k,reservoir_pressure_pa,reservoir_density_kg_m3,reservoir_mass_kg,body_temperature_k,wall_heat_into_gas_j,gas_work_j,pair_mass_residual_kg,total_energy_residual_j"
    );
    for step in 0..=steps {
        if step > 0 {
            let volume = piston.gas().volume_m3();
            let candidate = piston.real_gas_nozzle_thermal_candidate(
                &reservoir,
                nozzle,
                PistonInput {
                    duration_s: DT,
                    heat_j: 0.0,
                    exterior_pressure_pa: EXTERIOR_PRESSURE,
                    external_force_n: 0.0,
                },
                UniformEnthalpyStepInput {
                    curve: &curve,
                    initial: thermal,
                    mass_kg: BODY_MASS,
                    volume_m3: BODY_MASS / BODY_DENSITY,
                    surface_area_m2: SURFACE_AREA,
                    duration_s: DT,
                    internal_heat_j: 0.0,
                },
                PistonSolve {
                    volume_bracket_m3: [0.95 * volume, 1.05 * volume],
                    gas: energy,
                    momentum_tolerance_n_s: 1e-9,
                    energy_tolerance_j: 1e-6,
                    maximum_iterations: 128,
                },
                flow,
                exchange,
                |gas, body| surface.resolve(gas, body),
                || false,
            )?;
            let (next, next_reservoir, next_body, report) = candidate.into_parts();
            wall_heat += report.thermal.exchanged_heat_j;
            gas_work += report.thermal.piston.gas.work_by_gas_j;
            piston = next;
            reservoir = next_reservoir;
            thermal = next_body.state;
        }
        let gas = piston.gas();
        let mass_residual = gas.mass_kg() + reservoir.mass_kg() - initial_pair_mass;
        let gas_energy = relative_energy(gas, &reference)?
            + relative_energy(&reservoir, &reference)?
            - initial_relative_energy;
        let body_energy = BODY_MASS * (thermal.specific_enthalpy_j_kg() - initial_body_enthalpy);
        let kinetic = 0.5 * piston.mass_kg() * piston.velocity_m_s().powi(2);
        let exterior_work = EXTERIOR_PRESSURE * (gas.volume_m3() - INITIAL_VOLUME);
        let total_residual = gas_energy + body_energy + kinetic + exterior_work;
        if step > 0
            && (total_residual.abs() > step as f64 * exchange.energy_tolerance_j
                || mass_residual.abs() > step as f64 * flow.transfer.mass_tolerance_kg)
        {
            return Err(format!("continuing three-owner energy/mass account failed at step {step}: {total_residual} J, {mass_residual} kg").into());
        }
        println!(
            "{step},{:.9},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{wall_heat:.12e},{gas_work:.12e},{mass_residual:.12e},{total_residual:.12e}",
            step as f64 * DT,
            gas.state().temperature_k(),
            gas.state().pressure_pa(),
            gas.density_kg_m3(),
            gas.mass_kg(),
            gas.volume_m3(),
            piston.velocity_m_s(),
            reservoir.state().temperature_k(),
            reservoir.state().pressure_pa(),
            reservoir.density_kg_m3(),
            reservoir.mass_kg(),
            thermal.temperature_k()
        );
    }
    Ok(())
}
