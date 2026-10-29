//! Heat-driven decomposition on the original one-dimensional thermal graph.
//!
//! cargo run -p fs-thermal --example heated_pyrolysis_slab --
//!   [noncharring|authored-char] [end_s=60] [max_dt_s=0.01] [cells=48] [incident_W_m2=50000]
//!
//! The non-charring scalar inputs are the retained NIST 2021 MaCFP PMMA set.
//! This is an authored heating scenario, not a reproduction of an experiment:
//! 0.01 m² by 6 mm, initially 300.15 K, one irradiated/convecting face, insulated
//! back and sides. The char case changes the residue yield to 0.2 and char Cp
//! to 1200 J/(kg K) to exercise retained residue; it is NOT a wood property set.
//!
//! Effective escaping volatile enthalpy closes the supplied constant heat of
//! pyrolysis exactly. It is not an MMA gas caloric model and must not be passed
//! to a gas solver with a different formation/enthalpy reference. The separate
//! source-bound gas adapter is needed for that coupling.
//!
//! A bounded surface balance retains front half-cell conduction, convection,
//! absorbed irradiation and emitted radiation. Its net heat is held over each
//! requested interval. The original thermal graph supplies interior conduction;
//! the reactive owner solves local temperature and finite-rate mass loss. This
//! is first order and fixed geometry: no shrinkage, moving surface, flame,
//! porosity-dependent conductivity or experimental validation is implied.

use fs_thermal::enthalpy::pyrolysis::{
    Arrhenius, CondensedCaloric, PyrolysisControls, PyrolysisMaterial, PyrolysisSimulation,
};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, ThermalLink, TimedEvent,
    series_conductance_w_k,
};
use std::error::Error;
use std::io::{self, BufWriter, Write};

const INITIAL_T: f64 = 300.15;
const AREA: f64 = 0.01;
const THICKNESS: f64 = 0.006;
const DENSITY: f64 = 1100.0;
const CONDUCTIVITY: f64 = 0.20;
const CP_VIRGIN: f64 = 2200.0;
const PYROLYSIS_ENTHALPY: f64 = 870_000.0;
const EMISSIVITY: f64 = 0.9;
const HTC: f64 = 10.0;
const SIGMA: f64 = 5.670_374_419e-8;
const DOMAIN: [f64; 2] = [250.0, 1000.0];

struct Caloric {
    cp: f64,
}
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = INITIAL_T + h / self.cp;
        if t.is_finite() && (DOMAIN[0]..=DOMAIN[1]).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / self.cp
    }
}
impl CondensedCaloric for Caloric {
    fn specific_enthalpy_j_kg(&self, t: f64) -> Result<f64, DepositionError> {
        if !t.is_finite() || !(DOMAIN[0]..=DOMAIN[1]).contains(&t) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(self.cp * (t - INITIAL_T))
    }
}

fn insulated_boundary() -> Boundary {
    Boundary {
        ambient_k: INITIAL_T,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w: 0.0,
    }
}

// The monotonically increasing surface residual includes the half-cell
// resistance. The identical accepted watts enter the original thermal owner.
fn front_heat(center_t: f64, cells: usize, incident: f64) -> Result<(f64, f64), Box<dyn Error>> {
    let dx = THICKNESS / cells as f64;
    let conductance = series_conductance_w_k(AREA, [0.25 * dx; 2], [CONDUCTIVITY; 2])?;
    let absorbed = EMISSIVITY * AREA * incident;
    let residual = |t: f64| {
        conductance * (t - center_t)
            + EMISSIVITY
                * SIGMA
                * AREA
                * (t - INITIAL_T)
                * (t + INITIAL_T)
                * (t * t + INITIAL_T * INITIAL_T)
            + HTC * AREA * (t - INITIAL_T)
            - absorbed
    };
    let mut lo = center_t.min(INITIAL_T);
    let mut hi = center_t.max(INITIAL_T) + absorbed / conductance;
    if residual(lo) > 0.0 || residual(hi) < 0.0 || !hi.is_finite() {
        return Err("front surface heat balance is not bracketed".into());
    }
    for _ in 0..80 {
        let mid = lo + 0.5 * (hi - lo);
        let r = residual(mid);
        if r.abs() <= 1e-10 || mid == lo || mid == hi {
            return Ok((mid, conductance * (mid - center_t)));
        }
        if r < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Err("front surface heat balance exhausted its iteration budget".into())
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() > 5 {
        return Err("usage: heated_pyrolysis_slab [noncharring|authored-char] [end_s] [max_dt_s] [cells] [incident_W_m2]".into());
    }
    let name = args.first().map(String::as_str).unwrap_or("noncharring");
    let char_yield = match name {
        "noncharring" => 0.0,
        "authored-char" => 0.2,
        _ => return Err("material scenario must be noncharring or authored-char".into()),
    };
    let number = |index: usize, default: f64| -> Result<f64, Box<dyn Error>> {
        let value = args
            .get(index)
            .map(|s| s.parse())
            .transpose()?
            .unwrap_or(default);
        if !f64::is_finite(value) || value <= 0.0 {
            return Err("numerical arguments must be finite and positive".into());
        }
        Ok(value)
    };
    let end = number(1, 60.0)?;
    let maximum_dt = number(2, 0.01)?;
    let cells = args
        .get(3)
        .map(|s| s.parse::<usize>())
        .transpose()?
        .unwrap_or(48);
    let incident = number(4, 50_000.0)?;
    if end > 600.0 || maximum_dt > 0.1 || !(1..=128).contains(&cells) || incident > 100_000.0 {
        return Err("example bounds: end<=600 s, max_dt<=0.1 s, cells=1..128, incident<=100000 W/m2; depletion can require a moving-geometry model".into());
    }
    let intervals = (end / maximum_dt).ceil();
    if !intervals.is_finite() || !(1.0..=1_000_000.0).contains(&intervals) {
        return Err("example requires one through one million boundary intervals".into());
    }
    let intervals = intervals as usize;
    let virgin = Caloric { cp: CP_VIRGIN };
    let char_caloric = Caloric { cp: 1200.0 };
    let volatile = |t: f64| -> Result<f64, DepositionError> {
        let hv = virgin.specific_enthalpy_j_kg(t)?;
        let hc = char_caloric.specific_enthalpy_j_kg(t)?;
        Ok((hv + PYROLYSIS_ENTHALPY - char_yield * hc) / (1.0 - char_yield))
    };
    let dx = THICKNESS / cells as f64;
    let mass = DENSITY * AREA * dx;
    let conductance = series_conductance_w_k(AREA, [0.5 * dx; 2], [CONDUCTIVITY; 2])?;
    let mut events = Vec::with_capacity(2 * cells - 1);
    for i in 0..cells {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                segment_id: 0,
                body_id: 0,
                material_id: MaterialId(1),
                mass_kg: mass,
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [i as f64 * dx, 0.0, 0.0],
                    end_m: [(i + 1) as f64 * dx, 0.0, 0.0],
                    volume_m3: AREA * dx,
                    surface_area_m2: 2.0 * AREA,
                },
                boundary: insulated_boundary(),
            }),
        });
        if i > 0 {
            events.push(TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: i as u64,
                    a: CellId(i as u64 - 1),
                    b: CellId(i as u64),
                    area_m2: AREA,
                    conductance_w_k: conductance,
                }),
            });
        }
    }
    let thermal = DepositionSimulation::new(
        events,
        &[(MaterialId(1), &virgin)],
        DepositionBudget {
            max_cells: 128,
            max_links: 127,
            max_events: 255,
            max_substeps: 10_000,
        },
    )?
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: maximum_dt,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )?
    .0;
    let mut solid = PyrolysisSimulation::new(
        thermal,
        &[PyrolysisMaterial {
            material_id: MaterialId(1),
            virgin: &virgin,
            char_caloric: if char_yield > 0.0 {
                Some(&char_caloric)
            } else {
                None
            },
            volatile: &volatile,
            char_mass_yield: char_yield,
            kinetics: Arrhenius {
                prefactor_s_inv: 2.85e13,
                activation_energy_j_mol: 191_000.0,
            },
            temperature_range_k: DOMAIN,
            source: "MaCFP NIST PMMA2021 scalar inputs at9ae21bebd5dbc7478b8637fbc4df1c81291ec6ea; authored geometry/boundary/effective volatile enthalpy; optional authored char sensitivity",
        }],
        || Ok(()),
    )?;
    let initial_mass = solid.mass_kg();
    let initial_energy = solid.energy_j();
    let controls = PyrolysisControls {
        maximum_dt_s: maximum_dt,
        ..PyrolysisControls::default()
    };
    let mut input_j = 0.0;
    let mut volatile_j = 0.0;
    let mut volatile_kg = 0.0;
    let mut steps = 0usize;
    let mut substeps = 0usize;
    let mut maximum_cell_residual = 0.0_f64;
    let mut next_output = 0.0;
    let mut out = BufWriter::new(io::stdout().lock());
    writeln!(
        out,
        "time_s,front_face_k,front_cell_k,back_cell_k,solid_mass_kg,char_mass_kg,emitted_mass_kg,volatile_rate_kg_s,net_front_power_w,solid_energy_j,emitted_enthalpy_j,external_heat_j"
    )?;
    for interval in 1..=intervals {
        let time = solid.time_s();
        // Partition the requested duration directly; repeated clock addition
        // must not manufacture a near-zero final reaction interval.
        let target = end * (interval as f64 / intervals as f64);
        let dt = target - time;
        let (_, power) = front_heat(solid.temperature_k(CellId(0))?, cells, incident)?;
        let mut boundary = insulated_boundary();
        boundary.source_w = power;
        solid.set_boundary(CellId(0), boundary, || Ok(()))?;
        let (next, audit) = solid.advance_to(
            target,
            PyrolysisControls {
                maximum_dt_s: dt,
                ..controls
            },
            || Ok(()),
        )?;
        input_j += audit.thermal.external_heat_j;
        volatile_j += audit.volatile_enthalpy_j;
        volatile_kg += audit.volatile_mass_kg;
        maximum_cell_residual = maximum_cell_residual.max(audit.maximum_cell_energy_residual_j);
        substeps += audit.thermal.substeps;
        steps += 1;
        solid = next;
        if solid.time_s() >= next_output || solid.time_s() == end {
            let front_t = solid.temperature_k(CellId(0))?;
            let (face_t, _) = front_heat(front_t, cells, incident)?;
            let char_mass: f64 = solid.inventories().values().map(|s| s.char_mass_kg).sum();
            writeln!(
                out,
                "{:.9},{face_t:.9},{front_t:.9},{:.9},{:.12e},{char_mass:.12e},{volatile_kg:.12e},{:.12e},{power:.9},{:.12e},{volatile_j:.12e},{input_j:.12e}",
                solid.time_s(),
                solid.temperature_k(CellId(cells as u64 - 1))?,
                solid.mass_kg(),
                audit.volatile_mass_kg / dt,
                solid.energy_j()
            )?;
            next_output = solid.time_s().floor() + 1.0;
        }
    }
    out.flush()?;
    let mass_residual = solid.mass_kg() + volatile_kg - initial_mass;
    let energy_residual = solid.energy_j() - initial_energy + volatile_j - input_j;
    eprintln!(
        "scenario={name}; status=estimated authored slab; no experimental validation; geometry=fixed"
    );
    eprintln!(
        "NIST input A=2.85e13/s E=191000J/mol cp=2200J/kg/K rho=1100kg/m3 k=0.2W/m/K pyrolysis_enthalpy=870000J/kg emissivity=0.9"
    );
    eprintln!(
        "area_m2={AREA}; thickness_m={THICKNESS}; cells={cells}; incident_w_m2={incident}; htc_w_m2_k={HTC}; char_yield={char_yield}; max_dt_s={maximum_dt}"
    );
    eprintln!(
        "steps={steps}; thermal_substeps={substeps}; final_solid_kg={:.12e}; emitted_volatile_kg={volatile_kg:.12e}; mass_residual_kg={mass_residual:.12e}; energy_residual_j={energy_residual:.12e}; maximum_cell_energy_residual_j={maximum_cell_residual:.12e}",
        solid.mass_kg()
    );
    if mass_residual.abs() > 1e-10 || energy_residual.abs() > 1e-3 {
        return Err("whole-slab mass/energy balance exceeded the example gate".into());
    }
    Ok(())
}
