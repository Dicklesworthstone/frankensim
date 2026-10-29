//! Prescribed-HRR two-zone Steckler comparison with original gas/vent/wall owners.
//!
//! cargo run -p fs-airflow --features compartment-fire --example compartment_fire --
//!   data/reference/compartment_fire/steckler_1982_cases.csv 10 10 10 0.5 1800 4
//!
//! Required: cases CSV, case number, inner h [W/m2/K], outer h, maximum dt [s].
//! Optional: end time [s] (1800), wall cells (4). Coefficients are explicitly
//! chosen engineering inputs, not fitted material or validation parameters.
//! The CFAST deck's wall approximation is used without retuning. This prescribed
//! external heat-source case does not consume fuel/O2 or make products/soot;
//! species-limited burning is a separate source operator. Gas radiation and
//! weak-plume penetration remain outside this example. Radiation from the fire
//! is distributed uniformly by enclosure area; no gas emissivity is inferred.

#[path = "compartment_fire/air.rs"]
mod air;

use fs_airflow::compartment::fire::PrescribedFire;
use fs_airflow::compartment::wall::{
    GasContact, HeatTransferCoefficient, SpecificEnthalpy, ThermalConductivity, Wall1d,
    WallBoundary, WallGeometry, WallSolve,
};
use fs_airflow::compartment::{Ambient, Endpoint, NetworkControl, Room, WallVent, ZoneNetwork};
use fs_airflow::stratified::{DischargeCoefficient, WallOpening};
use fs_airflow::{SourceProvenance, ToleranceBasis};
use fs_gas::multicomponent::compartment::{Layer, TwoZone};
use fs_qty::{
    Acceleration, Area, Density, Dimensionless, Length, Power, Pressure, Temperature, Time,
};
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};
use std::error::Error;

struct WallLaw {
    initial_temperature: f64,
}
impl ThermalLaw for WallLaw {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = self.initial_temperature + h / 1000.0;
        if t.is_finite() && (200.0..=2000.0).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.001
    }
}

struct Case {
    id: String,
    width: f64,
    depth: f64,
    height: f64,
    door_width: f64,
    door_height: f64,
    heat: f64,
    ambient: f64,
    observed_upper: f64,
    observed_lower: f64,
    observed_interface: f64,
}
fn positive(s: &str) -> Result<f64, Box<dyn Error>> {
    let v: f64 = s.parse()?;
    if v.is_finite() && v > 0.0 {
        Ok(v)
    } else {
        Err("arguments and case quantities must be finite and positive".into())
    }
}
fn read_case(path: &str, number: &str) -> Result<Case, Box<dyn Error>> {
    let id = format!("steckler-1982-{number}");
    let contents = std::fs::read_to_string(path)?;
    let mut selected = None;
    for row in contents
        .lines()
        .filter(|r| !r.starts_with('#') && !r.trim().is_empty())
    {
        let fields: Vec<_> = row.split(',').map(str::trim).collect();
        if fields.first().copied() != Some(id.as_str()) {
            continue;
        }
        if fields.len() != 11 || row.contains('"') || selected.is_some() {
            return Err("expected one selected case with eleven unquoted fields".into());
        }
        let values = fields[1..]
            .iter()
            .map(|v| positive(v))
            .collect::<Result<Vec<_>, _>>()?;
        let c = Case {
            id: id.clone(),
            width: values[0],
            depth: values[1],
            height: values[2],
            door_width: values[3],
            door_height: values[4],
            heat: values[5],
            ambient: values[6],
            observed_upper: values[7],
            observed_lower: values[8],
            observed_interface: values[9],
        };
        if c.door_width > c.width
            || c.door_height > c.height
            || c.observed_upper <= c.ambient
            || c.observed_lower <= c.ambient
            || c.observed_interface >= c.height
        {
            return Err("inadmissible selected geometry or observed temperature rise".into());
        }
        selected = Some(c);
    }
    selected.ok_or_else(|| "requested case is absent from the retained CSV".into())
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if !(5..=7).contains(&args.len()) {
        return Err("usage: compartment_fire <cases.csv> <case-number> <inner-h-W/m2/K> <outer-h-W/m2/K> <max-dt-s> [end-s=1800] [wall-cells=4]".into());
    }
    let c = read_case(&args[0], &args[1])?;
    let h_inside = positive(&args[2])?;
    let h_outside = positive(&args[3])?;
    let maximum_dt = positive(&args[4])?;
    let end = args
        .get(5)
        .map(|s| positive(s))
        .transpose()?
        .unwrap_or(1800.0);
    let cells = args
        .get(6)
        .map(|s| s.parse::<usize>())
        .transpose()?
        .unwrap_or(4);
    if end > 1800.0 || !(1..=32).contains(&cells) {
        return Err("case supports at most 1800 s; wall cells must be 1..=32".into());
    }
    let (gas, fractions) = air::source()?;
    let floor_area = c.width * c.depth;
    let room = TwoZone::new(
        gas,
        Area::new(floor_area),
        Length::new(c.height),
        Length::new(c.height),
        Pressure::new(101325.0),
        [Temperature::new(c.ambient); 2],
        [fractions.clone(), fractions.clone()],
        || false,
    )?;
    let discharge = DischargeCoefficient::new(
        Dimensionless::new(0.7),
        SourceProvenance::new(
            "NIST TN 1889v1, section 4.1, nominal CFAST wall-opening Cd=0.7",
            "NIST-TN-1889v1:4.1",
        ),
        0.0,
        ToleranceBasis::EngineeringAllowance,
    )?;
    let mut network = ZoneNetwork::new(
        vec![Room {
            gas: room,
            floor_elevation: Length::new(0.0),
        }],
        Ambient {
            pressure: Pressure::new(101325.0),
            reference_height: Length::new(0.0),
            temperature: Temperature::new(c.ambient),
            mass_fractions: fractions,
        },
        vec![WallVent {
            left: Endpoint::Room(0),
            right: Endpoint::Ambient,
            opening: WallOpening::new(
                Length::new(0.0),
                Length::new(c.door_height),
                Length::new(c.door_width),
                discharge,
            )?,
        }],
        vec![],
        Acceleration::new(9.81),
        || false,
    )?;
    let fire = PrescribedFire {
        total_power: Power::new(c.heat),
        convective_power: Power::new(0.8 * c.heat),
        // Equivalent diameter of the authored CFAST 0.09 m2 burner; the
        // original experimental diameter is 0.30 m, retained separately in data.
        diameter: Length::new((4.0 * 0.09 / core::f64::consts::PI).sqrt()),
        base_height: Length::new(0.02),
    };
    let perimeter = 2.0 * (c.width + c.depth);
    let vertical_area = perimeter * c.height - c.door_width * c.door_height;
    let areas = [vertical_area, floor_area, floor_area]; // walls, ceiling, floor
    let enclosure_area: f64 = areas.iter().sum();
    let law = WallLaw {
        initial_temperature: c.ambient,
    };
    let mut walls = areas
        .iter()
        .map(|a| {
            Wall1d::new(
                &law,
                WallGeometry {
                    area: Area::new(*a),
                    thickness: Length::new(0.013),
                    density: Density::new(200.0),
                    conductivity: ThermalConductivity::new(1.0),
                    cells,
                },
                SpecificEnthalpy::new(0.0),
                || false,
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    let mut work = walls
        .iter()
        .map(Wall1d::workspace)
        .collect::<Result<Vec<_>, _>>()?;
    let initial_energy = network.rooms()[0].gas.physical_internal_energy().value()
        + walls.iter().map(|w| w.simulation().energy_j()).sum::<f64>();
    let mut expected_energy_change = 0.0;
    let mut max_mass_residual = 0.0_f64;
    let mut max_energy_residual = 0.0_f64;
    let mut max_pressure_residual = 0.0_f64;
    let mut total_halvings = 0;
    let mut steps = 0;
    let mut average_upper = 0.0;
    let mut average_lower = 0.0;
    let mut average_interface = 0.0;
    let mut observation_duration = 0.0;
    let mut next_output = 0.0;
    eprintln!(
        "Case={}, original NASA9 N2/O2 200..1000 K, prescribed HRR={} W, convective fraction=0.8",
        c.id, c.heat
    );
    eprintln!(
        "CFAST authored wall approximation: k=1 W/m/K, rho=200 kg/m3, cp=1000 J/kg/K, thickness=0.013 m; cells={cells}."
    );
    eprintln!(
        "Explicit h_inside={h_inside}, h_outside={h_outside} W/m2/K; nominal Cd=0.7 with zero analysis allowance and unquantified physical discrepancy."
    );
    eprintln!(
        "Uniform area-weighted incident fire radiation; no participating-gas radiation, pyrolysis, burning chemistry or fitted coefficient."
    );
    println!(
        "time_s,lower_temperature_k,upper_temperature_k,interface_height_m,pressure_pa,wall_inner_cell_temperature_k,ceiling_inner_cell_temperature_k,floor_inner_cell_temperature_k"
    );
    while network.rooms()[0].gas.time().value() < end {
        if steps >= 200000 {
            return Err("bounded example step budget exhausted".into());
        }
        let before_time = network.rooms()[0].gas.time().value();
        let duration = Time::new(maximum_dt.min(end - before_time));
        let (report, candidates) = network.step_with_sources(
            duration,
            NetworkControl::default(),
            |rooms, dt, rates, cancelled| {
                let room = &rooms[0].gas;
                let source = fire.evaluate(room, Acceleration::new(9.81), &mut *cancelled)?;
                rates[0] = source.rates;
                let z = room.interface_height().value();
                let lower_side = perimeter * z - c.door_width * z.min(c.door_height);
                let lower_ceiling = if room.layer(Layer::Upper).mass_kg() == 0.0 {
                    floor_area
                } else {
                    0.0
                };
                let lower_floor = if room.layer(Layer::Lower).mass_kg() > 0.0 {
                    floor_area
                } else {
                    0.0
                };
                let lower_areas = [lower_side, lower_ceiling, lower_floor];
                let mut candidates = Vec::with_capacity(3);
                for i in 0..3 {
                    let contact_areas = [lower_areas[i], areas[i] - lower_areas[i]];
                    let contact = |k, layer| GasContact {
                        area: Area::new(contact_areas[k]),
                        temperature: room.layer(layer).temperature(),
                        htc: HeatTransferCoefficient::new(h_inside),
                    };
                    let candidate = walls[i].candidate(
                        dt,
                        WallBoundary {
                            gas: [contact(0, Layer::Lower), contact(1, Layer::Upper)],
                            ambient_temperature: Temperature::new(c.ambient),
                            ambient_htc: HeatTransferCoefficient::new(h_outside),
                            external_radiant_power: Power::new(
                                source.radiant_power.value() * areas[i] / enclosure_area,
                            ),
                        },
                        WallSolve::default(),
                        &mut work[i],
                        &mut *cancelled,
                    )?;
                    for k in 0..2 {
                        rates[0][k].add_heat(candidate.report.gas_heat_power[k])?;
                    }
                    candidates.push(candidate);
                }
                Ok(candidates)
            },
            || false,
        )?;
        expected_energy_change += report.external_physical_energy_j;
        max_mass_residual = max_mass_residual.max(report.mass_residual_kg.abs());
        max_pressure_residual = max_pressure_residual.max(report.maximum_pressure_residual_pa);
        let mut step_residual = report.physical_energy_residual_j;
        for (wall, candidate) in walls.iter_mut().zip(candidates) {
            let r = candidate.report;
            expected_energy_change += r.ambient_heat_energy.value()
                + r.external_radiant_energy.value()
                - r.gas_heat_energy.iter().map(|q| q.value()).sum::<f64>();
            step_residual += r.audit.residual_j;
            // Infallible publication only after the complete network accepted.
            *wall = candidate.wall;
        }
        max_energy_residual = max_energy_residual.max(step_residual.abs());
        total_halvings += report.halvings;
        steps += 1;
        let room = &network.rooms()[0].gas;
        let time = room.time().value();
        let overlap = (time.min(1800.0) - before_time.max(1500.0)).max(0.0);
        average_upper += overlap * room.layer(Layer::Upper).temperature().value();
        average_lower += overlap * room.layer(Layer::Lower).temperature().value();
        average_interface += overlap * room.interface_height().value();
        observation_duration += overlap;
        if time >= next_output || time == end {
            println!(
                "{time:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9}",
                room.layer(Layer::Lower).temperature().value(),
                room.layer(Layer::Upper).temperature().value(),
                room.interface_height().value(),
                room.pressure().value(),
                walls[0].cell_temperature(0)?.value(),
                walls[1].cell_temperature(0)?.value(),
                walls[2].cell_temperature(0)?.value()
            );
            next_output = (time / 10.0).floor() * 10.0 + 10.0;
        }
    }
    let final_energy = network.rooms()[0].gas.physical_internal_energy().value()
        + walls.iter().map(|w| w.simulation().energy_j()).sum::<f64>();
    eprintln!(
        "steps={steps},whole_step_halvings={total_halvings},max_mass_residual_kg={max_mass_residual:.6e},max_combined_step_energy_residual_j={max_energy_residual:.6e},max_transport_pressure_residual_pa={max_pressure_residual:.6e},combined_energy_residual_j={:.6e}",
        final_energy - initial_energy - expected_energy_change
    );
    if end == 1800.0 && observation_duration > 0.0 {
        let predicted_upper = average_upper / observation_duration;
        let predicted_rise = predicted_upper - c.ambient;
        let observed_rise = c.observed_upper - c.ambient;
        eprintln!(
            "Steckler reduced-layer comparison, mean over 1500..1800 s: upper_k={predicted_upper:.9},interface_m={:.9},predicted_upper_rise_k={predicted_rise:.9},observed_upper_rise_k={observed_rise:.9},rise_ratio={:.9}",
            average_interface / observation_duration,
            predicted_rise / observed_rise
        );
        eprintln!(
            "Other retained zone quantities: lower_k={:.9},observed_lower_k={:.9},observed_interface_m={:.9},predicted_hot_layer_depth_m={:.9},observed_hot_layer_depth_m={:.9}",
            average_lower / observation_duration,
            c.observed_lower,
            c.observed_interface,
            c.height - average_interface / observation_duration,
            c.height - c.observed_interface
        );
        eprintln!(
            "Reduced-QoI standard uncertainty is unknown; this comparison is not a calibrated model or validation pass."
        );
    } else {
        eprintln!("Short run: the public 1500..1800 s comparison interval was not completed.");
    }
    Ok(())
}
