//! Three resolved passes, reheating, constrained cooling and prescribed release.
//!
//! Synthetic coefficients, not calibrated filament/adhesion data. Axially uniform
//! layers are perfectly bonded by declaration. One-way thermal-to-mechanical
//! coupling: mechanical dissipation is reported, NOT returned as thermal heat.
//! The mechanics cadence resolves every birth and release; refine it separately
//! from thermal subcycling and section subdivision. Existing geometry does not
//! feed back into thermal contact. The bed release is prescribed, not predicted.
//!
//! This development example compiles the ORIGINAL std-only conduction producer
//! by path, avoiding a new production fs-solid -> fs-conduction dependency. Normal
//! applications can import the same producer from fs_conduction::transient.
//! Run: cargo run -p fs-solid --example deposited_thermovisco_strip -- 0.25
//! Test: cargo test -p fs-solid --example deposited_thermovisco_strip

#[path = "../../fs-conduction/src/transient/deposition.rs"]
mod thermal;

use std::collections::BTreeMap;
use std::error::Error;
use fs_material::thermovisco::{Limits, ThermalKnot, ThermoMaxwell};
use fs_solid::deposition::{AxialControl, BendingControl, Budget, Config, Control, LayerGeometry, LayeredStrip, Report};

const LENGTH_M: f64 = 0.01;
const WIDTH_M: f64 = 0.002;
const HEIGHT_M: f64 = 0.0002;
const DENSITY: f64 = 1200.0;
const CP: f64 = 1800.0;
const RELEASE_S: f64 = 30.0;
const RELEASE_DURATION_S: f64 = 2.0;
const END_S: f64 = 60.0;

struct ConstantCp;
impl thermal::ThermalLaw for ConstantCp {
    fn temperature_k(&self, specific_h_j_kg: f64) -> Result<f64, thermal::DepositionError> {
        let t = 300.0+specific_h_j_kg/CP;
        if t.is_finite() && (300.0..=550.0).contains(&t) { Ok(t) }
        else { Err(thermal::DepositionError::MaterialDomain) }
    }
    fn maximum_temperature_slope(&self) -> f64 { 1.0/CP }
}

fn material() -> Result<ThermoMaxwell, Box<dyn Error>> {
    // e_inf and branch pairs are compatible with visco::GeneralizedMaxwell.
    Ok(ThermoMaxwell::new(2e8, vec![(6e8, 8.0), (2e8, 80.0)], vec![
        ThermalKnot { temperature_k: 300.0, modulus_scale: 1.0, relaxation_scale: 1.0, free_strain: 0.0 },
        ThermalKnot { temperature_k: 400.0, modulus_scale: 0.25, relaxation_scale: 0.05, free_strain: 0.001 },
        ThermalKnot { temperature_k: 550.0, modulus_scale: 0.02, relaxation_scale: 0.005, free_strain: 0.0025 },
    ], Limits { max_strain: 0.03, max_strain_rate_s: 1.0 })?)
}
fn geometry(id: u64) -> LayerGeometry {
    LayerGeometry { area_m2: WIDTH_M*HEIGHT_M, centroid_z_m: (id as f64-2.0)*HEIGHT_M, thickness_m: HEIGHT_M }
}
fn boundary(on_bed: bool) -> thermal::Boundary {
    thermal::Boundary { ambient_k: 300.0, htc_w_m2_k: 20.0, source_w: 0.0,
        reservoir: on_bed.then_some(thermal::Reservoir {
            area_m2: LENGTH_M*WIDTH_M, conductance_w_k: 0.0004, temperature_k: 330.0,
        }) }
}
fn schedule() -> Vec<thermal::TimedEvent> {
    let mut events = Vec::new();
    for id in 1..=3 {
        let time_s = 10.0*(id-1) as f64; let g = geometry(id);
        events.push(thermal::TimedEvent { time_s, event: thermal::Event::Birth(thermal::DepositedCell {
            id: thermal::CellId(id), material_id: thermal::MaterialId(7), body_id: 1, segment_id: id,
            mass_kg: DENSITY*g.area_m2*LENGTH_M, specific_h_j_kg: CP*200.0,
            geometry: thermal::CellGeometry { start_m: [0.0, 0.0, g.centroid_z_m],
                end_m: [LENGTH_M, 0.0, g.centroid_z_m], volume_m3: g.area_m2*LENGTH_M,
                surface_area_m2: 2.0*(LENGTH_M*WIDTH_M+LENGTH_M*HEIGHT_M+WIDTH_M*HEIGHT_M) },
            boundary: boundary(id == 1),
        }) });
        if id > 1 {
            events.push(thermal::TimedEvent { time_s, event: thermal::Event::Contact(thermal::ThermalLink {
                id: id-1, a: thermal::CellId(id-1), b: thermal::CellId(id),
                area_m2: LENGTH_M*WIDTH_M, conductance_w_k: 0.001,
            }) });
        }
    }
    events.push(thermal::TimedEvent { time_s: RELEASE_S,
        event: thermal::Event::Boundary(thermal::CellId(1), boundary(false)) });
    events
}
fn config() -> Config {
    Config { force_tolerance_n: 1e-8, moment_tolerance_nm: 1e-10, energy_tolerance_j: 1e-9,
        max_strain: 0.03, max_rotation_rad: 0.1 }
}
struct Run<'a> {
    mechanics: LayeredStrip<'a>,
    thermal: thermal::DepositionSimulation<'a>,
    report: Report,
    before_release: [f64; 2],
    first_released: [f64; 2],
    first_pass_before_second_k: f64,
    first_pass_reheated_k: f64,
    thermal_defect_j: f64,
    mechanical_defect_l1_j: f64,
    physical_dissipation_j: f64,
    numerical_dissipation_j: f64,
}

fn simulate<'a>(law: &'a ThermoMaxwell, cp: &'a ConstantCp, dt: f64) -> Result<Run<'a>, Box<dyn Error>> {
    if !dt.is_finite() || !(0.01..=0.25).contains(&dt) {
        return Err("mechanics step must lie in 0.01..0.25 seconds".into());
    }
    let thermal = thermal::DepositionSimulation::new(schedule(), &[(thermal::MaterialId(7), cp)],
        thermal::DepositionBudget { max_cells: 3, max_links: 2, max_events: 6, max_substeps: 1000 })?;
    let (mut thermal, initial_audit) = thermal.advance_to(0.0, thermal::StepConfig {
        max_dt_s: dt, energy_tolerance_j: 1e-9 }, || Ok(()))?;
    let mut mechanics = LayeredStrip::new(LENGTH_M, Budget { max_layers: 3, max_branch_states: 12 })?
        .deposit(1, geometry(1), law, thermal.temperature_k(thermal::CellId(1))?)?;
    let mut incoming = initial_audit.incoming_energy_j; let mut external = initial_audit.external_heat_j;
    let mut report = Report::default(); let mut before_release = [0.0; 2]; let mut first_released = [0.0; 2];
    let mut before_second = 0.0; let mut reheated = 0.0;
    let mut release_load = [0.0; 2];
    let mut mechanical_defect_l1_j = 0.0; let mut physical_dissipation_j = 0.0; let mut numerical_dissipation_j = 0.0;
    while thermal.time_s() < END_S {
        let start = thermal.time_s();
        let stop = [10.0, 20.0, RELEASE_S, RELEASE_S+RELEASE_DURATION_S, END_S].into_iter().filter(|t| *t > start)
            .fold((start+dt).min(END_S), f64::min);
        // Both candidates use the same interval. Publish neither if either fails.
        let (candidate_thermal, audit) = thermal.advance_to(stop, thermal::StepConfig {
            max_dt_s: dt, energy_tolerance_j: 1e-9 }, || Ok(()))?;
        let temperatures: BTreeMap<_, _> = mechanics.layers().keys().map(|id|
            candidate_thermal.temperature_k(thermal::CellId(*id)).map(|t| (*id, t))).collect::<Result<_, _>>()?;
        // Resolve a two-second unloading ramp, rather than interpreting an
        // instantaneous stiffness/reaction removal as physical viscous loss.
        let control = if start < RELEASE_S { Control::HELD } else {
            let remaining = (1.0-(stop-RELEASE_S)/RELEASE_DURATION_S).max(0.0);
            Control { axial: AxialControl::ForceN(remaining*release_load[0]),
                bending: BendingControl::MomentNm(remaining*release_load[1]) }
        };
        let (mut candidate_mechanics, candidate_report) = mechanics.advance_to(stop, &temperatures, control, config(), || Ok(()))?;
        // Material born at the END has no preceding mechanics/thermal aging.
        for id in 1..=3 {
            if candidate_thermal.cells().contains_key(&thermal::CellId(id)) && !candidate_mechanics.layers().contains_key(&id) {
                candidate_mechanics = candidate_mechanics.deposit(id, geometry(id), law,
                    candidate_thermal.temperature_k(thermal::CellId(id))?)?;
            }
        }
        if stop == 10.0 { before_second = candidate_thermal.temperature_k(thermal::CellId(1))?; }
        if start == 10.0 { reheated = candidate_thermal.temperature_k(thermal::CellId(1))?; }
        if stop == RELEASE_S {
            before_release = candidate_mechanics.deformed_position_m(LENGTH_M, 0.0)?;
            release_load = [candidate_report.axial_force_n, candidate_report.bending_moment_nm];
        }
        if start == RELEASE_S { first_released = candidate_mechanics.deformed_position_m(LENGTH_M, 0.0)?; }
        incoming += audit.incoming_energy_j; external += audit.external_heat_j;
        mechanical_defect_l1_j += candidate_report.energy_j.residual.abs();
        physical_dissipation_j += candidate_report.energy_j.dissipated;
        numerical_dissipation_j += candidate_report.energy_j.numerical_dissipation;
        thermal = candidate_thermal; mechanics = candidate_mechanics; report = candidate_report;
    }
    let thermal_defect_j = thermal.energy_j()-incoming-external;
    Ok(Run { mechanics, thermal, report, before_release, first_released,
        first_pass_before_second_k: before_second, first_pass_reheated_k: reheated,
        thermal_defect_j, mechanical_defect_l1_j, physical_dissipation_j, numerical_dissipation_j })
}

fn main() -> Result<(), Box<dyn Error>> {
    let dt: f64 = std::env::args().nth(1).map_or(Ok(0.25), |s| s.parse())?;
    let law = material()?; let cp = ConstantCp; let run = simulate(&law, &cp, dt)?;
    let tip = run.mechanics.deformed_position_m(LENGTH_M, 0.0)?;
    println!("Synthetic one-way thermo-Maxwell strip; all lengths in metres, temperatures K.");
    println!("time_s={}, layers={}, mass_kg={:.12e}", run.thermal.time_s(), run.mechanics.layers().len(), run.thermal.mass_kg());
    println!("first_pass_before_second_K={:.9}, reheated_K={:.9}", run.first_pass_before_second_k, run.first_pass_reheated_k);
    println!("held_tip_m={:?}, first_released_tip_m={:?}, cooled_tip_m={tip:?}", run.before_release, run.first_released);
    println!("axial_strain={:.12e}, curvature_per_m={:.12e}, force_residual_N={:.12e}, moment_residual_Nm={:.12e}",
        run.mechanics.axial_strain(), run.mechanics.curvature_per_m(), run.report.force_residual_n, run.report.moment_residual_nm);
    println!("thermal_defect_J={:.12e}, mechanical_defect_L1_J={:.12e}, dashpot_dissipation_J={:.12e}, numerical_dissipation_J={:.12e}",
        run.thermal_defect_j, run.mechanical_defect_l1_j, run.physical_dissipation_j, run.numerical_dissipation_j);
    Ok(())
}

#[cfg(test)]
mod coupled_tests {
    use super::*;
    #[test]
    fn actual_deposition_cools_reheats_and_deforms_only_after_prescribed_release() {
        let law = material().unwrap(); let cp = ConstantCp; let run = simulate(&law, &cp, 0.25).unwrap();
        assert_eq!(run.thermal.time_s(), run.mechanics.time_s());
        assert_eq!(run.mechanics.layers().len(), 3);
        assert!(run.first_pass_reheated_k > run.first_pass_before_second_k);
        assert_eq!(run.before_release[1], 0.0); assert!(run.first_released[1].abs() > 1e-9);
        assert!(run.thermal_defect_j.abs() < 1e-9 && run.mechanical_defect_l1_j < 1e-8);
        assert!(run.physical_dissipation_j > 0.0 && run.numerical_dissipation_j > 0.0);
        assert!(run.report.force_residual_n.abs() < 1e-8 && run.report.moment_residual_nm.abs() < 1e-10);
        assert!((run.thermal.mass_kg()-3.0*DENSITY*WIDTH_M*HEIGHT_M*LENGTH_M).abs() < 1e-14);
    }
    #[test]
    fn independent_coupling_cadence_refinement_reduces_geometry_change() {
        let law = material().unwrap(); let cp = ConstantCp;
        let a = simulate(&law, &cp, 0.25).unwrap();
        let b = simulate(&law, &cp, 0.125).unwrap();
        let c = simulate(&law, &cp, 0.0625).unwrap();
        let coarse = (a.mechanics.curvature_per_m()-b.mechanics.curvature_per_m()).abs();
        let fine = (b.mechanics.curvature_per_m()-c.mechanics.curvature_per_m()).abs();
        assert!(fine < coarse, "coarse {coarse}, fine {fine}");
    }
}
