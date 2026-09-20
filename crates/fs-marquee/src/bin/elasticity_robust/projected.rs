//! Authored-load, same-material mode of the existing robust elasticity binary.
use super::{load_cases, numbers, parse, writer};
use fs_cutfem::DesignBoxEdge;
use fs_topols::robust_descent::{
    MultiLoadProjectedAttempt, MultiLoadProjectedOptimizer, MultiLoadProjectedProgress,
    MultiLoadProjectedSettings, MultiLoadProjectedState,
};
use fs_topols::volume::VolumeProjectionSettings;
use fs_topols::{GridSdf, OptimizeSettings, RobustAggregate, RobustLoadCase};
use std::error::Error;
use std::fs::File;
use std::io::{Read, Write};
use std::path::Path;

fn quoted(value: &str) -> String {
    let mut out = String::new();
    out.push('"');
    for ch in value.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            ch if ch.is_control() => {
                use std::fmt::Write as _;
                let _ = write!(out, "\\u{:04x}", ch as u32);
            }
            ch => out.push(ch),
        }
    }
    out.push('"');
    out
}

fn state_json(state: &MultiLoadProjectedState) -> String {
    format!(
        "{{\"objective\":{:.17e},\"case_compliances\":{},\"active_case\":{},\"volume\":{:.17e},\"snapshot\":\"{:#018x}\"}}",
        state.objective, numbers(&state.case_compliances),
        state.active_case.map_or_else(|| "null".into(), |index| index.to_string()),
        state.volume, state.snapshot,
    )
}

fn write_field(path: &Path, field: &GridSdf) -> Result<(), Box<dyn Error>> {
    let mut file = writer(path)?;
    writeln!(file, "x_normalized,y_normalized,phi_normalized")?;
    for j in 0..=field.n() {
        for i in 0..=field.n() {
            let [x, y] = field.pos(i, j);
            writeln!(file, "{x:.17e},{y:.17e},{:.17e}", field.node(i, j))?;
        }
    }
    file.flush()?;
    Ok(())
}

// Admit exactly the exporter layout: no omitted, duplicated, reordered or
// non-finite samples. Parse coordinates instead of silently ignoring them.
fn read_field(path: &Path, n: usize) -> Result<GridSdf, Box<dyn Error>> {
    const MAX_BYTES: u64 = 8 * 1024 * 1024;
    let mut text = String::new();
    File::open(path)?.take(MAX_BYTES + 1).read_to_string(&mut text)?;
    if text.len() as u64 > MAX_BYTES {
        return Err("initial level-set CSV exceeds 8 MiB".into());
    }
    let mut lines = text.lines().filter(|line| !line.trim().is_empty());
    if lines.next().map(str::trim) != Some("x_normalized,y_normalized,phi_normalized") {
        return Err("initial field requires x_normalized,y_normalized,phi_normalized header".into());
    }
    let mut field = GridSdf::from_fn(n, &|_, _| 0.0);
    for j in 0..=n {
        for i in 0..=n {
            let row = lines.next().ok_or("initial field has missing nodal rows")?;
            let values = row.split(',').map(str::trim).collect::<Vec<_>>();
            if values.len() != 3 { return Err("initial field rows require x,y,phi".into()); }
            let x: f64 = values[0].parse()?;
            let y: f64 = values[1].parse()?;
            let phi: f64 = values[2].parse()?;
            if [x, y] != field.pos(i, j) || !phi.is_finite() {
                return Err(format!("initial field has wrong coordinates or non-finite phi at node ({i},{j})").into());
            }
            *field.node_mut(i, j) = phi;
        }
    }
    if lines.next().is_some() { return Err("initial field has extra nodal rows".into()); }
    Ok(field)
}

fn write_loads(path: &Path, cases: &[RobustLoadCase]) -> Result<(), Box<dyn Error>> {
    let mut file = writer(path)?;
    writeln!(file, "# edge,start,end,fx,fy,weight")?;
    for case in cases {
        let edge = match case.edge() {
            DesignBoxEdge::Left => "left",
            DesignBoxEdge::Right => "right",
            DesignBoxEdge::Top => "top",
            DesignBoxEdge::Bottom => "bottom",
        };
        let [start, end] = case.interval();
        let [fx, fy] = case.traction();
        writeln!(file, "{edge},{start:.17e},{end:.17e},{fx:.17e},{fy:.17e},{:.17e}", case.weight())?;
    }
    file.flush()?;
    Ok(())
}

fn write_attempts(
    file: &mut impl Write, iteration: usize, attempts: &[MultiLoadProjectedAttempt],
) -> std::io::Result<()> {
    for attempt in attempts {
        let state = attempt.state.as_ref().map_or_else(|| "null".into(), state_json);
        let refusal = attempt.refusal.as_deref().map_or_else(|| "null".into(), quoted);
        let projection = attempt.projection.map_or_else(|| "null".into(), |report| format!(
            "{{\"volume\":{:.17e},\"shift\":{:.17e},\"evaluations\":{}}}",
            report.volume, report.shift, report.evaluations,
        ));
        writeln!(file,
            "{{\"iteration\":{iteration},\"candidate\":{},\"move_cells\":{:.17e},\"projection\":{projection},\"state\":{state},\"refusal\":{refusal}}}",
            attempt.index, attempt.move_cells,
        )?;
    }
    file.flush()
}

pub(super) fn run(args: &[String]) -> Result<u8, Box<dyn Error>> {
    if !(2..=9).contains(&args.len()) {
        return Err("usage: fs-marquee-elasticity-robust --projected OUTPUT_DIR LOAD_CASES.csv [LEVEL=4] [UPDATES=30] [AREA=0.45] [CANDIDATES=6] [AGGREGATE=worst] [MAX_SOLVES] [INITIAL_FIELD.csv]".into());
    }
    let output = Path::new(&args[0]);
    if output.try_exists()? { return Err("output directory already exists; refusing to overwrite it".into()); }
    let level: u32 = parse(args, 2, 4)?;
    let iterations: usize = parse(args, 3, 30)?;
    let volfrac: f64 = parse(args, 4, 0.45)?;
    let max_candidates: usize = parse(args, 5, 6)?;
    let aggregate = match args.get(6).map(String::as_str).unwrap_or("worst") {
        "sum" | "weighted-sum" => RobustAggregate::WeightedSum,
        "worst" | "worst-weighted" => RobustAggregate::WorstWeightedCase,
        _ => return Err("AGGREGATE must be sum or worst".into()),
    };
    if !(2..=7).contains(&level) || !(1..=200).contains(&iterations)
        || !(1..=16).contains(&max_candidates)
        || !(volfrac.is_finite() && volfrac > 0.0 && volfrac < 1.0)
    {
        return Err("projected mode requires LEVEL 2..=7, UPDATES 1..=200, CANDIDATES 1..=16 and AREA (0,1)".into());
    }
    let cases = load_cases(Path::new(&args[1]))?;
    let default_solves = iterations.checked_mul(max_candidates)
        .and_then(|attempts| attempts.checked_add(1))
        .and_then(|families| families.checked_mul(cases.len()))
        .ok_or("complete study solve budget overflow")?;
    let max_solves: usize = parse(args, 7, default_solves)?;
    let n = 1usize << level;
    let field = match args.get(8) {
        Some(path) => read_field(Path::new(path), n)?,
        None => GridSdf::from_fn(n, &|_, y| (y - 0.5).abs() - 0.42),
    };
    let fixed = field.nodes().iter().copied().enumerate().filter(|(index, _)| {
        let i = index % (n + 1);
        let j = index / (n + 1);
        i == 0 || i == n || j == 0 || j == n
    }).collect();
    let projection = VolumeProjectionSettings {
        target: volfrac, tolerance: 1e-4, max_shift: 2.0, max_evaluations: 64,
    };
    let settings = OptimizeSettings {
        level, iterations, volfrac, move_cells: 0.35,
        nucleation_period: 4, hole_radius_cells: 1.5, ..OptimizeSettings::default()
    };
    let mut optimizer = MultiLoadProjectedOptimizer::new(
        field.clone(), &cases, settings, aggregate, fixed, projection,
        MultiLoadProjectedSettings { max_candidates, max_solves, ..MultiLoadProjectedSettings::default() },
    )?;
    std::fs::create_dir(output)?;
    write_field(&output.join("input-level-set.csv"), &field)?;
    write_field(&output.join("baseline-level-set.csv"), optimizer.geometry())?;
    write_loads(&output.join("load-cases.csv"), &cases)?;
    let mut trace = writer(&output.join("trajectory.jsonl"))?;
    let mut attempts = writer(&output.join("attempts.jsonl"))?;
    let (status, exit, failure) = loop {
        let iteration = optimizer.next_iteration();
        match optimizer.advance_one() {
            Ok(MultiLoadProjectedProgress::Accepted(step)) => {
                write_attempts(&mut attempts, iteration, &step.attempts)?;
                writeln!(trace,
                    "{{\"iteration\":{},\"previous\":{},\"state\":{},\"projection_shift\":{:.17e},\"proposal_drift_h\":{:.17e},\"proposal_nucleation_count\":{},\"proposal_load_pad_nodes\":{},\"solves_started\":{}}}",
                    step.iteration, state_json(&step.previous), state_json(&step.state),
                    step.projection.shift, step.proposal_audit.interface_drift_h,
                    step.proposal_events.len(), step.proposal_load_pad_nodes, optimizer.solves_started(),
                )?;
                trace.flush()?;
            }
            Ok(MultiLoadProjectedProgress::IterationLimit) => break ("iteration_limit", 0, None),
            Ok(MultiLoadProjectedProgress::NoDescent(rows)) => {
                write_attempts(&mut attempts, iteration, &rows)?;
                break ("no_descent", 11, None);
            }
            Ok(MultiLoadProjectedProgress::SolveBudget(rows)) => {
                write_attempts(&mut attempts, iteration, &rows)?;
                break ("solve_budget", 13, None);
            }
            Err(error) => break ("refused", 12, Some(error.to_string())),
        }
    };
    // Keep the last fully accepted state even after numerical or budget stops.
    // Summary is the LAST export: a failed field/trace write cannot print success.
    trace.flush()?;
    attempts.flush()?;
    write_field(&output.join("level-set.csv"), optimizer.geometry())?;
    let aggregate_name = match aggregate {
        RobustAggregate::WeightedSum => "weighted_sum",
        RobustAggregate::WorstWeightedCase => "worst_weighted_case",
    };
    let failure_json = failure.as_deref().map_or_else(|| "null".into(), quoted);
    let summary = format!(
        "{{\"schema\":\"projected-multiload-v1\",\"model\":\"normalized_unit_square_plane_strain\",\"authority\":\"estimated\",\"status\":\"{status}\",\"aggregate\":\"{aggregate_name}\",\"level\":{level},\"requested_updates\":{iterations},\"accepted_updates\":{},\"load_cases\":{},\"area_target\":{volfrac:.17e},\"area_tolerance\":{:.17e},\"candidate_budget\":{max_candidates},\"max_solves\":{max_solves},\"solves_started\":{},\"baseline\":{},\"final\":{},\"refusal\":{failure_json},\"claims\":{{\"physical_validation\":false,\"kkt_convergence\":false,\"global_optimum\":false,\"continuum_volume_certificate\":false,\"three_dimensional\":false}}}}",
        optimizer.next_iteration(), cases.len(), projection.tolerance, optimizer.solves_started(),
        state_json(optimizer.baseline()), state_json(&optimizer.current()),
    );
    let mut file = writer(&output.join("summary.json"))?;
    writeln!(file, "{summary}")?;
    file.flush()?;
    println!("{summary}");
    Ok(exit)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn path(name: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!("frankensim-projected-{name}-{}-{}.csv",
            std::process::id(), std::thread::current().name().unwrap_or("test")))
    }

    #[test]
    fn exported_field_round_trips_exact_nodal_bits() {
        let path = path("field-round-trip");
        let field = GridSdf::from_fn(8, &|x, y| x - y + 0.0123);
        write_field(&path, &field).unwrap();
        let replay = read_field(&path, 8).unwrap();
        assert_eq!(field.nodes().iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
            replay.nodes().iter().map(|v| v.to_bits()).collect::<Vec<_>>());
    }

    #[test]
    fn wrong_coordinates_and_truncated_fields_refuse() {
        let path = path("wrong-coordinates");
        std::fs::write(&path, "x_normalized,y_normalized,phi_normalized\n0.125,0,-0.1\n").unwrap();
        assert!(read_field(&path, 8).is_err());
        let path = path.with_extension("truncated.csv");
        std::fs::write(&path, "x_normalized,y_normalized,phi_normalized\n0,0,-0.1\n").unwrap();
        assert!(read_field(&path, 8).is_err());
    }
}
