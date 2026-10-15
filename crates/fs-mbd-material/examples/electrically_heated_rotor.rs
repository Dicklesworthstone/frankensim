//! A data-file-driven current -> heat -> free expansion -> joint-response path.
//! See docs/MATERIAL_ELECTROTHERMAL_MECHANISM.md for physical scope and usage.
#[path = "electrically_heated_rotor/source.rs"]
mod source;

use std::io::{Read, Write};
use fs_exec::CancelGate;
use fs_ga::{Se3, Vec3, Wrench};
use fs_lbm::conjugate::{ElectricSetup, Electrode, ElectrodeDrive, VoxelDomain};
use fs_lbm::conjugate::enthalpy::JouleEnthalpyInput;
use fs_material::state_point::MaterialPropertySelection;
use fs_mbd::articulated::{BaseState, JointModel, forward_dynamics};
use fs_mbd_material::{RegionBinding, VolumeMoments};
use fs_mbd_material::articulated::{MaterialLink, ResolvedArticulated, resolve_articulated};
use fs_mbd_material::thermal::heating::{HeatedArticulated, LinkHeating, ThermalRegion};
use fs_thermal::enthalpy::*;
use source::Source;

type Error = Box<dyn std::error::Error>;
fn invalid(s: &str) -> Error { std::io::Error::new(std::io::ErrorKind::InvalidInput, s.to_string()).into() }
const KEYS: [&str; 3] = ["left-electrode", "rotor", "right-electrode"];
const CELLS: [CellId; 3] = [CellId(0), CellId(1), CellId(2)];
const LINKS: [u64; 2] = [70, 71];

fn reference(source: &Source) -> Result<ResolvedArticulated, Error> {
    let p = source::point(source.reference_k)?;
    let cube = VolumeMoments::cuboid(Vec3::new(1., 1., 1.))?;
    let mut regions = Vec::new();
    for (slot, key) in KEYS.into_iter().enumerate() {
        regions.push(RegionBinding { key, material: &source.card, point: &p,
            geometry: VolumeMoments::new(1., Vec3::new(slot as f64 + 0.5, 0.5, 0.5), cube.inertia_m5())?,
            selection: MaterialPropertySelection::SingleClaimOnly });
    }
    let mount = [regions[0].clone(), regions[2].clone()];
    Ok(resolve_articulated(&[
        MaterialLink { name: "mount", parent: None, parent_from_child_zero: Se3::identity(),
            joint: JointModel::FIXED, regions: &mount },
        MaterialLink { name: "rotor", parent: Some(0), parent_from_child_zero: Se3::identity(),
            joint: JointModel::revolute(Vec3::new(0., 0., 1.), None)?, regions: &regions[1..2] },
    ], 3, || false)?)
}
fn heated<'a>(source: &'a Source, reference: &'a ResolvedArticulated) -> Result<HeatedArticulated<'a>, Error> {
    let region = |slot: usize| ThermalRegion { key: KEYS[slot], card: &source.card, caloric: &source.caloric,
        cell_id: CELLS[slot], material_id: MaterialId(slot as u64), segment_id: slot as u64,
        body_id: if slot == 1 { 1 } else { 0 },
        geometry: CellGeometry { start_m: [slot as f64, 0.5, 0.5], end_m: [slot as f64 + 1., 0.5, 0.5],
            volume_m3: 1., surface_area_m2: 6. },
        initial_specific_enthalpy_j_kg: 0., boundary: Boundary { ambient_k: source.reference_k,
            htc_w_m2_k: 0., reservoir: None, source_w: 0. },
        expansion_selection: MaterialPropertySelection::SingleClaimOnly, fixed_point_m: Vec3::default(),
    };
    let events = (0..2).map(|i| TimedEvent { time_s: 0., event: Event::Contact(ThermalLink {
        id: LINKS[i], a: CELLS[i], b: CELLS[i+1], area_m2: 1., conductance_w_k: 0.,
    }) }).collect();
    Ok(HeatedArticulated::new(reference, vec![
        LinkHeating { link: "mount", regions: vec![region(0), region(2)] },
        LinkHeating { link: "rotor", regions: vec![region(1)] },
    ], events, DepositionBudget { max_cells: 3, max_links: 2, max_events: 5, max_substeps: 100_000 },
        StepConfig { max_dt_s: 0.25, energy_tolerance_j: 1e-8 }, || Ok(()))?)
}
fn electrodes(voltage: f64) -> ElectricSetup {
    ElectricSetup { conductors: vec![], electrodes: vec![
        Electrode { lo: [0,0,0], hi: [1,1,1], drive: ElectrodeDrive::Voltage(voltage) },
        Electrode { lo: [2,0,0], hi: [3,1,1], drive: ElectrodeDrive::Voltage(0.) },
    ] }
}

#[derive(Clone, Copy, Debug)]
struct Row {
    time: f64, voltage: f64, temperature: [f64; 3], current: f64,
    thermal_energy: f64, source_work: f64, volume: f64, mass: f64, acceleration: f64,
}

fn run(source: &Source, voltage: f64, off: f64, end: f64, max_dt: f64,
    mut output: impl FnMut(Row) -> Result<(), Error>,
) -> Result<(), Error> {
    if !voltage.is_finite() || !off.is_finite() || !end.is_finite()
        || !max_dt.is_finite() || max_dt <= 0. || off <= 0. || end <= off
        || end / max_dt > 10_000. {
        return Err(invalid("finite voltage, 0 < shutoff < end, positive step, at most 10000 requested steps required"));
    }
    let reference = reference(source)?;
    let mut body = heated(source, &reference)?;
    let domain = VoxelDomain::new(3,1,1,1.)?;
    let on = electrodes(voltage);
    let shorted = electrodes(0.);
    let gate = CancelGate::new();
    let config = StepConfig { max_dt_s: max_dt, energy_tolerance_j: 1e-8 };
    let mut work = 0.;
    let mut total_substeps = 0usize;
    while body.thermal().time_s() < end {
        let time = body.thermal().time_s();
        let driven = time < off;
        let target = (time + max_dt).min(if driven { off } else { end });
        if target <= time { return Err(invalid("time step cannot advance the physical clock")); }
        let remaining = 100_000usize.checked_sub(total_substeps)
            .filter(|n| *n > 0).ok_or_else(|| invalid("total example thermal work budget exhausted"))?;
        let input = JouleEnthalpyInput { domain: &domain, setup: if driven { &on } else { &shorted },
            contacts: &[], cell_ids: &CELLS, max_cells: 3, max_current_solves: remaining + 1,
            power_tolerance_w: 1e-8, energy_tolerance_j: 1e-8 };
        let (next, report) = body.advance_joule(target, config, input, &LINKS, &gate,
            |_,_,t| source.transport(t), |patch, sides| {
                // Actual adjacent coefficients, explicit reference half-cell
                // series geometry, and a perfect thermal interface. No new law.
                let ka = sides[0].response.conductivity_w_m_k();
                let kb = sides[1].response.conductivity_w_m_k();
                Ok((patch.link.area_m2 / (0.5/ka + 0.5/kb), (ka,kb)))
            })?;
        total_substeps += report.joule.thermal.substeps;
        work += report.joule.electrode_energy_j;
        if !work.is_finite() { return Err(invalid("accumulated source work is unrepresentable")); }
        let mechanical = &next.configuration().bodies()[1];
        let dynamics = forward_dynamics(next.configuration().model(), BaseState::stationary(Se3::identity()),
            &[0.], &[0.], &[1.], Vec3::default(), &[Wrench::default();2])?;
        let row = Row { time: next.thermal().time_s(), voltage: if driven { voltage } else { 0. },
            temperature: [next.thermal().temperature_k(CELLS[0])?, next.thermal().temperature_k(CELLS[1])?,
                next.thermal().temperature_k(CELLS[2])?],
            current: report.joule.endpoint.electric.electrodes[0].current_a,
            thermal_energy: next.thermal().energy_j(), source_work: work,
            volume: mechanical.regions()[0].geometry().volume_m3(), mass: mechanical.inertia().mass(),
            acceleration: dynamics.generalized_acceleration[0] };
        // Publish a complete accepted interval. Earlier CSV rows are retained
        // on a later failure; no incomplete row or partial mechanism is emitted.
        output(row)?;
        body = next;
    }
    Ok(())
}

fn main() -> Result<(), Error> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() != 5 {
        return Err(invalid("usage: electrically_heated_rotor MATERIAL.txt VOLTAGE_V SHUTOFF_S END_S MAX_DT_S"));
    }
    let file = std::fs::File::open(&args[0])?;
    let mut text = String::new();
    file.take(source::MAX_BYTES as u64 + 1).read_to_string(&mut text)?;
    let source = Source::parse(&text, &args[0])?;
    let [voltage, off, end, dt] = [args[1].parse()?, args[2].parse()?, args[3].parse()?, args[4].parse()?];
    let stdout = std::io::stdout();
    let mut out = stdout.lock();
    writeln!(out, "time_s,applied_voltage_v,left_k,rotor_k,right_k,current_a,thermal_j,electrode_work_j,rotor_volume_m3,rotor_mass_kg,test_acceleration_rad_s2")?;
    run(&source, voltage, off, end, dt, |r| {
        writeln!(out, "{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.12},{:.12},{:.12},{:.12},{:.12}",
            r.time, r.voltage, r.temperature[0], r.temperature[1], r.temperature[2], r.current,
            r.thermal_energy, r.source_work, r.volume, r.mass, r.acceleration)?;
        Ok(())
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    const CP10: &str = include_str!("electrically_heated_rotor/cp10.txt");
    const CP20: &str = include_str!("electrically_heated_rotor/cp20.txt");
    #[test]
    fn file_only_material_swap_reaches_current_heat_and_joint_response() {
        let mut traces = Vec::new();
        for text in [CP10, CP20] {
            let source = Source::parse(text, "same fixture label").unwrap();
            let mut rows = Vec::new();
            run(&source,20.,2.,4.,0.25,|row| {rows.push(row);Ok(())}).unwrap();
            assert_eq!(rows.len(),16);
            let hot = rows[7]; let cool = rows[15];
            assert_eq!(hot.time,2.); assert_eq!(cool.time,4.);
            assert!(hot.current>0.); assert_eq!(cool.current,0.);
            assert!(cool.temperature[1]<hot.temperature[1]);
            assert!(cool.temperature[0]>hot.temperature[0]);
            assert!(cool.acceleration>hot.acceleration);
            assert_eq!(hot.mass,cool.mass);
            assert_eq!(hot.source_work,cool.source_work);
            assert!((hot.thermal_energy-cool.thermal_energy).abs()<1e-8);
            traces.push(rows);
        }
        assert!(traces[0][7].temperature[1]>traces[1][7].temperature[1]);
        assert!(traces[0][7].current<traces[1][7].current);
        assert!(traces[0][7].acceleration<traces[1][7].acceleration);
    }
    #[test]
    fn malformed_data_and_unresolvable_clocks_refuse() {
        for text in ["", source::HEADER, "wrong header\n300 2 10 .002 5 2\n500 2 10 .002 15 6",
            "frankensim-electrothermal-material-v1\n300 2 NaN .002 5 2\n500 2 10 .002 15 6",
            "frankensim-electrothermal-material-v1\n300 2 10 .002 5 2\n300 2 10 .002 15 6"] {
            assert!(Source::parse(text,"bad").is_err());
        }
        let source=Source::parse(CP10,"fixture").unwrap();
        assert!(run(&source,20.,2.,4.,0.,|_|Ok(())).is_err());
        assert!(run(&source,20.,2.,4.,1e-9,|_|Ok(())).is_err());
        assert!(run(&source,f64::NAN,2.,4.,0.25,|_|Ok(())).is_err());
    }
}
