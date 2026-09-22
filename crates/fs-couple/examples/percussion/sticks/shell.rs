//! A second work-conjugate vertical tip on the original curved shell.
use super::{Error, ImpactBody, Obstacle, Port, Stroke};
use fs_plate::shell::reduction::ShellReduction;

/// Retain the original [first stick, shell modes] prefix, then append an
/// independent inertial stick. Both contacts react on the SAME nonlinear shell.
/// Positions are exact reference-XY surface locations; no edge/hole snapping.
pub fn build(stroke: Stroke, reduction: &ShellReduction, nodes: &[[f64;3]],
    triangles: &[[usize;3]], coordinate: usize, total: usize,
) -> Result<(ImpactBody, Obstacle, Port), Error> {
    if reduction.mode_count() == 0 || coordinate != 1+reduction.mode_count()
        || total != coordinate+1 {
        return Err("second cymbal stick requires the unchanged first-stick/shell prefix and one appended inertia".into());
    }
    let position = stroke.position_m.ok_or("second cymbal stick requires an explicit physical XY station")?;
    let (triangle,bary) = crate::playing::shell_location(nodes,triangles,position)?;
    // The shell reduction owns rotations, curvature and mass normalization.
    // Do not normalize this row: b.v is the actual physical tip-site velocity.
    let shapes = reduction.point_port(triangle,bary,[0.0,0.0,-1.0])?;
    if shapes.iter().all(|b| *b == 0.0) {
        return Err("second cymbal stick has no retained vertical surface participation".into());
    }
    let (body,weight) = crate::stick_with_speed(stroke.speed_m_s)?;
    let mut row = vec![0.0;total];row[coordinate] = weight;
    for (i,b) in shapes.iter().enumerate() {row[1+i] = -b;}
    let contact = crate::elastic_contact(row)?;
    eprintln!("second physical cymbal stick: xy_m={position:?}, launch_m_s={}, coordinate={coordinate}; independent inertia, shared curved shell and stand felt; no direct radiation or stick-stick collision",stroke.speed_m_s);
    Ok((body,contact,Port {coordinate,weight}))
}

#[cfg(test)]
mod tests;
