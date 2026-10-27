//! The accepted local heat sequence drives one retained native solid trajectory.
//! Source states and problems are retained before motion starts so their borrowed
//! identities outlive the run. No equilibrium reset or alternate integrator is used.
use super::*;
use fs_couple::solid_enthalpy::SolidEnthalpyState;
use fs_solid::loads3::transient::{TetTransientRun3, TransientConfig3};

#[derive(Clone, Copy, Debug)]
pub(super) struct Options {
    pub maximum_step_s: f64,
    pub mass_damping_per_s: f64,
    /// Explicit Rayleigh stiffness coefficient [s], not an acoustic loss factor.
    pub stiffness_damping_s: f64,
}

impl Options {
    fn validate(self) -> Result<()> {
        if !self.maximum_step_s.is_finite()
            || self.maximum_step_s <= 0.0
            || !self.mass_damping_per_s.is_finite()
            || self.mass_damping_per_s < 0.0
            || !self.stiffness_damping_s.is_finite()
            || self.stiffness_damping_s < 0.0
        {
            return Err("inertial history needs a finite positive maximum step and finite nonnegative Rayleigh coefficients".into());
        }
        Ok(())
    }

    pub fn extract(args: &[String]) -> Result<(Vec<String>, Option<Self>)> {
        let mut remaining = Vec::new();
        let mut selected = None;
        let mut stiffness = None;
        let mut words = args.iter();
        while let Some(word) = words.next() {
            if word == "--inertial-stiffness-damping" {
                if stiffness.is_some() {
                    return Err("--inertial-stiffness-damping may be supplied only once".into());
                }
                stiffness = Some(words.next()
                    .ok_or("--inertial-stiffness-damping needs a coefficient in seconds")?
                    .parse::<f64>()?);
                continue;
            }
            if word != "--inertial-history" {
                remaining.push(word.clone());
                continue;
            }
            if selected.is_some() {
                return Err("--inertial-history may be supplied only once".into());
            }
            let options = Self {
                maximum_step_s: words
                    .next()
                    .ok_or("--inertial-history needs MAX_STEP_S MASS_DAMPING_PER_S")?
                    .parse()?,
                mass_damping_per_s: words
                    .next()
                    .ok_or("--inertial-history needs MAX_STEP_S MASS_DAMPING_PER_S")?
                    .parse()?,
                stiffness_damping_s: 0.0,
            };
            options.validate()?;
            selected = Some(options);
        }
        if let Some(value) = stiffness {
            let options = selected.as_mut()
                .ok_or("--inertial-stiffness-damping requires --inertial-history")?;
            options.stiffness_damping_s = value;
            options.validate()?;
        }
        Ok((remaining, selected))
    }

    /// This is a declared timestep ceiling, not a bound on integration error.
    pub fn admit_interval(self, interval_s: f64) -> Result<()> {
        self.validate()?;
        // Subtracting up to 1024 authored endpoints amplifies clock roundoff.
        if !interval_s.is_finite()
            || interval_s <= 0.0
            || interval_s > self.maximum_step_s * (1.0 + 4096.0 * f64::EPSILON)
        {
            return Err("heat trajectory exceeds the inertial maximum step; increase finite-gas STEPS or shorten HEAT_SECONDS".into());
        }
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(super) struct Frame {
    pub time_s: f64,
    pub step: usize,
    /// Full node-major displacement [m], velocity [m/s], acceleration [m/s²].
    pub motion: [Vec<[f64; 3]>; 3],
    pub kinetic_j: f64,
    pub elastic_j: f64,
    /// Cumulative external, thermal-expansion, stiffness work, damping loss [J].
    pub work_j: [f64; 4],
    pub residual_j: f64,
}

pub(super) struct History {
    pub frames: Vec<Frame>,
    nodes_m: Vec<[f64; 3]>,
    options: Options,
}

fn frame(
    run: &TetTransientRun3<'_, '_, '_, '_>,
    work_j: [f64; 4],
    initial_energy_j: f64,
    cx: &Cx<'_>,
) -> Result<Frame> {
    let response = run.recover_response(FRAME, 0.02, cx)?;
    let kinetic_j = run.state().kinetic_energy_j;
    let elastic_j = response.stress.elastic_energy_j();
    let residual_j = response.kinetic_plus_elastic_energy_j - initial_energy_j + work_j[3]
        - work_j[0]
        - work_j[1]
        - work_j[2];
    if work_j
        .iter()
        .chain([&kinetic_j, &elastic_j, &residual_j])
        .any(|v| !v.is_finite())
    {
        return Err("inertial history contains a nonfinite energy account".into());
    }
    Ok(Frame {
        time_s: run.state().time_s,
        step: run.state().step,
        motion: run.nodal_motion(cx)?,
        kinetic_j,
        elastic_j,
        work_j,
        residual_j,
    })
}

pub(super) fn evolve(
    source: &SolidEnthalpy<'_>,
    mesh: &Mesh,
    states: &[SolidEnthalpyState<'_>],
    options: Options,
    cx: &Cx<'_>,
) -> Result<History> {
    options.validate()?;
    cx.checkpoint().map_err(|_| "heated valve cancelled")?;
    if !(2..=1025).contains(&states.len()) || states[0].time_s() != 0.0 {
        return Err(
            "inertial history requires the initial state and 1..=1024 accepted heat endpoints"
                .into(),
        );
    }
    for pair in states.windows(2) {
        options.admit_interval(pair[1].time_s() - pair[0].time_s())?;
    }
    let problems: Vec<_> = states
        .iter()
        .map(|state| state.problem(&mesh.fixed))
        .collect();
    let loads = source
        .geometry()
        .assemble(&[], VolumeLoad3::None, &[], cx)?;
    let zero = vec![[0.0; 3]; mesh.nodes.len()];
    let controls = TransientConfig3 {
        maximum_steps: states.len() - 1,
        maximum_displacement_m: 0.0002,
        absolute_energy_tolerance_j: 1e-10,
        relative_energy_tolerance: 1e-8,
        mass_damping_per_s: options.mass_damping_per_s,
        // Both coefficients are explicit history inputs, never copied from acoustic loss.
        stiffness_damping_s: options.stiffness_damping_s,
        ..TransientConfig3::default()
    };
    let mut run = loads.begin_transient(
        &problems[0],
        Some(states[0].thermal_strains()),
        states[0].time_s(),
        &zero,
        &zero,
        controls,
        cx,
    )?;
    let initial_energy = run
        .recover_response(FRAME, 0.02, cx)?
        .kinetic_plus_elastic_energy_j;
    let mut work = [0.0; 4];
    let mut frames = vec![frame(&run, work, initial_energy, cx)?];
    for index in 1..states.len() {
        let (candidate, report) = states[index].advance_transient_candidate(
            &states[index - 1],
            &run,
            &problems[index],
            &loads,
            cx,
        )?;
        work[0] += report.native.mechanical_work_j + report.native.support_work_j;
        work[1] += report.thermal_material_work_j;
        work[2] += report.stiffness_material_work_j;
        work[3] += report.native.damping_loss_j;
        let accepted = frame(&candidate, work, initial_energy, cx)?;
        let scale = (accepted.kinetic_j + accepted.elastic_j)
            .abs()
            .max(work.iter().map(|v| v.abs()).sum::<f64>());
        let allowance = index as f64 * controls.absolute_energy_tolerance_j
            + controls.relative_energy_tolerance * scale;
        if accepted.time_s.to_bits() != states[index].time_s().to_bits()
            || accepted.step != index
            || accepted.residual_j.abs() > allowance
        {
            return Err("retained motion must share the accepted heat clock and close cumulative mechanical work".into());
        }
        frames.push(accepted);
        run = candidate;
    }
    cx.checkpoint().map_err(|_| "heated valve cancelled")?;
    Ok(History {
        frames,
        nodes_m: mesh.nodes.clone(),
        options,
    })
}

impl History {
    pub fn write(&self, path: &str) -> Result<()> {
        let mut out = std::io::BufWriter::new(
            std::fs::OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(path)?,
        );
        writeln!(
            out,
            "# retained native Newmark motion from stress-free rest; maximum_step_s={}; mass_damping_per_s={}; stiffness_damping_s={}; fixed reference mesh/mass and zero supports; prescribed one-way thermal/material work; global cumulative energies are repeated on each node row; acoustic output remains frozen-state experiments",
            self.options.maximum_step_s, self.options.mass_damping_per_s,
            self.options.stiffness_damping_s
        )?;
        writeln!(
            out,
            "time_s,step,node,reference_x_m,reference_y_m,reference_z_m,ux_m,uy_m,uz_m,vx_m_s,vy_m_s,vz_m_s,ax_m_s2,ay_m_s2,az_m_s2,kinetic_J,elastic_J,external_work_J,thermal_expansion_work_J,stiffness_work_J,damping_loss_J,mechanical_residual_J"
        )?;
        for frame in &self.frames {
            for (node, reference) in self.nodes_m.iter().enumerate() {
                write!(out, "{:.17e},{},{}", frame.time_s, frame.step, node)?;
                for value in reference
                    .iter()
                    .chain(frame.motion.iter().flat_map(|field| field[node].iter()))
                    .chain([&frame.kinetic_j, &frame.elastic_j])
                    .chain(&frame.work_j)
                    .chain([&frame.residual_j])
                {
                    write!(out, ",{value:.17e}")?;
                }
                writeln!(out)?;
            }
        }
        out.flush()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn extract(words: &[&str]) -> Result<(Vec<String>, Option<Options>)> {
        Options::extract(&words.iter().map(|word| (*word).to_owned()).collect::<Vec<_>>())
    }

    #[test]
    fn rayleigh_input_preserves_legacy_arguments_and_accepts_either_option_order() {
        let (remaining, legacy) = extract(&["--demo", "1", "--inertial-history", "0.0003", "20"]).unwrap();
        assert_eq!(remaining, vec!["--demo".to_owned(), "1".to_owned()]);
        assert_eq!(legacy.unwrap().stiffness_damping_s, 0.0);
        for words in [
            vec!["--inertial-history", "0.0003", "20", "--inertial-stiffness-damping", "0.00001"],
            vec!["--inertial-stiffness-damping", "0.00001", "--inertial-history", "0.0003", "20"],
        ] {
            let (remaining, options) = extract(&words).unwrap();
            assert!(remaining.is_empty());
            let options = options.unwrap();
            assert_eq!(options.mass_damping_per_s, 20.0);
            assert_eq!(options.stiffness_damping_s, 0.00001);
            options.admit_interval(0.00025).unwrap();
        }
    }

    #[test]
    fn stiffness_damping_requires_explicit_history_and_finite_nonnegative_data() {
        assert!(extract(&["--inertial-stiffness-damping", "0.001"]).is_err());
        assert!(extract(&["--inertial-stiffness-damping"]).is_err());
        for value in ["-1", "NaN", "inf"] {
            assert!(extract(&["--inertial-history", "0.0003", "20",
                "--inertial-stiffness-damping", value]).is_err());
        }
        assert!(extract(&["--inertial-history", "0.0003", "20",
            "--inertial-stiffness-damping", "0.001",
            "--inertial-stiffness-damping", "0.001"]).is_err());
    }
}
