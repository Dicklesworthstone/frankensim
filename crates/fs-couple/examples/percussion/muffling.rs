//! Fixed spatial mufflers in actual shell/head coordinates, upstream of sound.
//! These ideal bilateral viscous supports are explicitly supplied SI inputs,
//! not measured finger/gel constitutives, switched chokes or output envelopes.
use super::{Error, ModePair, ShellReduction, TensionedDisk};
use fs_couple::render::plate::impact::damping::{ViscousDamper, MAX_VISCOUS_DAMPERS};
use fs_couple::render::plate::impact::linear::wire::film_shapes;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Surface { Batter, Resonant, Shell }
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Muffler {
    pub surface: Surface,
    pub position_m: [f64; 2],
    pub resistance_n_s_m: f64,
}

/// Repeat --muffler SURFACE X_m Y_m C_Ns/m; preserve all other arguments.
pub fn options(args: &mut Vec<String>) -> Result<Vec<Muffler>, Error> {
    let mut result = Vec::new(); let mut i = 0;
    while i < args.len() {
        if args[i] != "--muffler" { i += 1; continue; }
        if result.len() == MAX_VISCOUS_DAMPERS || args.len()-i < 5 {
            return Err("--muffler requires surface, x_m, y_m, resistance_Ns/m; at most 16 ports".into());
        }
        let surface = match args[i+1].as_str() {
            "batter" => Surface::Batter, "resonant" => Surface::Resonant, "shell" => Surface::Shell,
            _ => return Err("muffler surface must be batter, resonant or shell".into()),
        };
        let spec = Muffler { surface, position_m: [args[i+2].parse()?, args[i+3].parse()?],
            resistance_n_s_m: args[i+4].parse()? };
        admit(&[spec], surface == Surface::Shell)?;
        result.push(spec); args.drain(i..i+5);
    }
    Ok(result)
}
fn admit(specs: &[Muffler], shell: bool) -> Result<(), Error> {
    if specs.len() > MAX_VISCOUS_DAMPERS {
        return Err("muffler count exceeds the physical port budget".into());
    }
    for spec in specs {
        if !spec.position_m.iter().all(|x| x.is_finite())
            || !spec.resistance_n_s_m.is_finite() || spec.resistance_n_s_m < 0.0 {
            return Err("mufflers require finite XY positions and nonnegative SI resistance".into());
        }
        if (spec.surface == Surface::Shell) != shell {
            return Err("shell mufflers require a splash command; batter/resonant mufflers require a drum or snare".into());
        }
    }
    Ok(())
}
pub fn admit_command(specs: &[Muffler], command: &str) -> Result<(), Error> {
    admit(specs, matches!(command, "splash"|"splash-wav"|"splash-mic"))
}

/// Actual vertical shell displacement. A mounting hole/outside point refuses;
/// no nearest-triangle snap and no normalization of the mechanical mode row.
pub fn shell_ports(specs: &[Muffler], reduction: &ShellReduction,
    nodes: &[[f64;3]], triangles: &[[usize;3]]) -> Result<Vec<ViscousDamper>, Error>
{
    admit(specs, true)?;
    specs.iter().map(|spec| {
        let (triangle, bary) = super::playing::shell_location(nodes, triangles, spec.position_m)?;
        let mut weights = vec![0.0]; // No direct damping on the separate stick.
        weights.extend(reduction.point_port(triangle, bary, [0.0,0.0,-1.0])?);
        Ok(ViscousDamper { weights, damping_n_s_m: spec.resistance_n_s_m })
    }).collect()
}

/// Original head addresses stay fixed. Second striker and all wire coordinates
/// occupy the zero suffix; the cavity compiler later appends exact gas zeros.
pub fn head_ports(specs: &[Muffler], films: &[TensionedDisk], modes: &[Vec<ModePair>],
    total: usize) -> Result<Vec<ViscousDamper>, Error>
{
    admit(specs, false)?;
    if films.len() != 2 || modes.len() != 2
        || modes.iter().any(Vec::is_empty)
        || total < 1+modes.iter().map(Vec::len).sum::<usize>() {
        return Err("head mufflers need the unchanged two-head structural layout".into());
    }
    specs.iter().map(|spec| {
        let head = usize::from(spec.surface == Surface::Resonant);
        let shapes = film_shapes(&films[head], &modes[head], &[spec.position_m])?.remove(0);
        if shapes.iter().all(|b| *b == 0.0) {
            return Err("muffler lies on a nonmoving head station; rim attachment is not simulated".into());
        }
        let start = if head == 0 { 1 } else { 1+modes[0].len() };
        let mut weights = vec![0.0; total];
        weights[start..start+shapes.len()].copy_from_slice(&shapes);
        Ok(ViscousDamper { weights, damping_n_s_m: spec.resistance_n_s_m })
    }).collect()
}

#[cfg(test)]
mod tests;
