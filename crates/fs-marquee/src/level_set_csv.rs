//! The normalized nodal CSV shared by the elasticity study exporters/readers.
use fs_topols::GridSdf;
use std::error::Error;
use std::fs::File;
use std::io::Read;
use std::path::Path;

/// Read one complete x-fastest dyadic lattice, without resampling its geometry.
///
/// Coordinates are checked, not discarded. Finite phi bits round-trip through
/// the existing 17-digit exporter. Import starts a NEW study; it does not restore
/// an optimizer's multiplier, iteration ordinal or search history.
///
/// # Errors
/// Refuses invalid lattice sizes, reads above 8 MiB, incompatible headers,
/// non-finite values, and missing, extra, duplicated or reordered nodes.
pub fn read_field(path: &Path, n: usize) -> Result<GridSdf, Box<dyn Error>> {
    if !n.is_power_of_two() || !(2..=256).contains(&n) {
        return Err("initial field requires a dyadic lattice with 2..=256 cells per side".into());
    }
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
