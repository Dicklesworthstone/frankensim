//! Conservative initial-stress work on the original tapered Rayleigh field.
use super::{
    PlateError, RoundInertia, RoundStation, bad, hermite, round_element_with_reference_material,
};
use fs_sparse::{Coo, Csr};

/// Full physical-coordinate tangent and inertia of a straight round shaft.
/// Coordinates are `[w, theta]` at every node, with `theta = dw/dx`.
/// No supports, rigid-mode projection or spectral truncation are applied.
/// Apply the same physical constraints to both matrices before using the
/// existing modal/stability owners. Compression can make `k` indefinite.
#[derive(Debug, Clone)]
pub struct RoundBeamPencil {
    /// Elastic plus initial-stress tangent.
    pub k: Csr,
    /// Original consistent translational and rotary inertia.
    pub m: Csr,
    /// Axial coordinates [m], including every original material interface.
    pub nodes_m: Vec<f64>,
}

/// Tapered, graded Rayleigh element with signed axial section force [N].
/// The force varies linearly between the supplied endpoints, positive in
/// tension and negative in compression. It adds the conservative potential
/// `1/2 * integral(N(x) * w'(x)^2 dx)` to the ORIGINAL Hermite displacement.
/// Three Gauss points integrate this degree-five polynomial exactly up to
/// rounding. Elastic moduli, density and inertia are unchanged; zero force
/// preserves [`round_element_with_reference_material`] bit for bit.
///
/// The caller supplies an equilibrated straight reference and its balancing
/// dead loads. Nonuniform forces require corresponding distributed/nodal axial
/// loads. This function neither solves that equilibrium nor infers supports,
/// follower forces, finite rotations or post-buckling motion. In particular,
/// dead axial force contributes work to an infinitesimal rigid rotation; it
/// must not enter `RoundBeamModes`' unloaded rigid-nullspace projection.
///
/// # Errors
/// Invalid material/geometry, nonfinite force or nonrepresentable tangent.
pub fn round_element_with_axial_force(
    length_m: f64,
    radii_m: [f64; 2],
    young_pa: [f64; 2],
    reference_density_kg_m3: [f64; 2],
    axial_force_n: [f64; 2],
) -> Result<([f64; 16], [f64; 16]), PlateError> {
    if axial_force_n.iter().any(|n| !n.is_finite()) {
        return Err(bad("round beam axial force endpoints must be finite"));
    }
    let (mut k, m) = round_element_with_reference_material(
        length_m,
        radii_m,
        young_pa,
        reference_density_kg_m3,
    )?;
    if axial_force_n == [0.0; 2] {
        return Ok((k, m));
    }
    for (t, weight) in [
        (0.1127016653792583, 5.0 / 18.0),
        (0.5, 4.0 / 9.0),
        (0.8872983346207417, 5.0 / 18.0),
    ] {
        let force = interpolate(axial_force_n, t);
        let (_, slope, _) = hermite(t, length_m);
        for i in 0..4 {
            for j in i..4 {
                let value = (length_m * weight * force) * slope[i] * slope[j];
                k[4 * i + j] += value;
                if i != j {
                    k[4 * j + i] += value;
                }
            }
        }
    }
    if k.iter().any(|v| !v.is_finite()) {
        return Err(bad("round beam initial-stress tangent overflow"));
    }
    Ok((k, m))
}

fn interpolate(pair: [f64; 2], t: f64) -> f64 {
    if pair[0] == pair[1] {
        pair[0]
    } else {
        (1.0 - t) * pair[0] + t * pair[1]
    }
}

/// Assemble the original tapered shaft with spatial elasticity, inertia and
/// signed equilibrated axial force. Each field needs one endpoint pair per
/// original geometry interval. Interfaces may jump; subdivision interpolates
/// within each original interval and never averages across an interface.
/// Attached bodies use the original reciprocal Hermite displacement/slope
/// maps, including bodies between nodes. Mass is independent of prestress.
///
/// The result retains ALL physical coordinates so supports and the signed
/// tangent can enter the existing modal and stability owners. Geometry
/// stations can be placed at arbitrary support locations. Unlike the unloaded
/// single-pin `RoundBeamModes`, no rigid rotation is removed or forced to zero.
/// The cold-work bound is 2..=33 original stations, 1..=16 subdivisions per
/// interval, at most 66 resulting nodes and eight attached bodies. These are
/// resolution/work bounds, not a slenderness or mesh-convergence certificate.
/// See [`round_element_with_axial_force`] for the physical applicability.
///
/// # Errors
/// Incomplete/nonphysical fields, invalid geometry/attachments, exceeded mesh
/// bounds or nonfinite element/assembled entries. No partial pencil escapes.
pub fn assemble_prestressed_round_beam(
    stations: &[RoundStation],
    subdivisions: usize,
    young_pa: &[[f64; 2]],
    reference_density_kg_m3: &[[f64; 2]],
    axial_force_n: &[[f64; 2]],
    attached: &[RoundInertia],
) -> Result<RoundBeamPencil, PlateError> {
    if !(2..=33).contains(&stations.len())
        || !(1..=16).contains(&subdivisions)
        || stations
            .iter()
            .any(|s| !s.x_m.is_finite() || !s.radius_m.is_finite() || s.radius_m < 0.0)
        || stations
            .windows(2)
            .any(|s| s[0].x_m >= s[1].x_m || s[0].radius_m == 0.0 && s[1].radius_m == 0.0)
    {
        return Err(bad("invalid bounded prestressed round beam geometry"));
    }
    let intervals = stations.len() - 1;
    if young_pa.len() != intervals
        || reference_density_kg_m3.len() != intervals
        || axial_force_n.len() != intervals
        || young_pa
            .iter()
            .flatten()
            .chain(reference_density_kg_m3.iter().flatten())
            .any(|v| !v.is_finite() || *v <= 0.0)
        || axial_force_n.iter().flatten().any(|v| !v.is_finite())
    {
        return Err(bad(
            "round beam needs complete positive E/rho and finite signed force fields",
        ));
    }
    let start = stations[0].x_m;
    let end = stations[intervals].x_m;
    if !(end - start).is_finite()
        || attached.len() > 8
        || attached.iter().any(|a| {
            !a.x_m.is_finite()
                || !(start..=end).contains(&a.x_m)
                || !a.mass_kg.is_finite()
                || a.mass_kg < 0.0
                || !a.rotary_kg_m2.is_finite()
                || a.rotary_kg_m2 < 0.0
                || a.mass_kg == 0.0 && a.rotary_kg_m2 == 0.0
        })
    {
        return Err(bad(
            "invalid prestressed round beam length or attached inertia",
        ));
    }
    let count = intervals * subdivisions + 1;
    if count > 66 {
        return Err(bad("prestressed round beam cold mesh exceeds 66 nodes"));
    }
    let mut nodes = Vec::with_capacity(count);
    for pair in stations.windows(2) {
        for i in 0..subdivisions {
            nodes.push(pair[0].x_m + (pair[1].x_m - pair[0].x_m) * i as f64 / subdivisions as f64);
        }
    }
    nodes.push(end);
    if nodes.windows(2).any(|x| !x[0].is_finite() || x[0] >= x[1]) {
        return Err(bad("prestressed round beam mesh lost axial resolution"));
    }
    let nd = 2 * nodes.len();
    let mut k = Coo::new(nd, nd);
    let mut m = Coo::new(nd, nd);
    for (cell, x) in nodes.windows(2).enumerate() {
        let interval = cell / subdivisions;
        let a = stations[interval];
        let b = stations[interval + 1];
        let field = |pair| [x[0], x[1]].map(|at| interpolate(pair, (at - a.x_m) / (b.x_m - a.x_m)));
        let (ke, me) = round_element_with_axial_force(
            x[1] - x[0],
            field([a.radius_m, b.radius_m]),
            field(young_pa[interval]),
            field(reference_density_kg_m3[interval]),
            field(axial_force_n[interval]),
        )?;
        for i in 0..4 {
            for j in 0..4 {
                k.push(2 * cell + i, 2 * cell + j, ke[4 * i + j]);
                m.push(2 * cell + i, 2 * cell + j, me[4 * i + j]);
            }
        }
    }
    for body in attached {
        let cell = nodes
            .partition_point(|x| *x <= body.x_m)
            .saturating_sub(1)
            .min(nodes.len() - 2);
        let length = nodes[cell + 1] - nodes[cell];
        let (n, d, _) = hermite((body.x_m - nodes[cell]) / length, length);
        for i in 0..4 {
            for j in 0..4 {
                m.push(
                    2 * cell + i,
                    2 * cell + j,
                    body.mass_kg * n[i] * n[j] + body.rotary_kg_m2 * d[i] * d[j],
                );
            }
        }
    }
    let k = k.assemble();
    let m = m.assemble();
    if (0..nd).any(|row| {
        k.row(row)
            .1
            .iter()
            .chain(m.row(row).1)
            .any(|v| !v.is_finite())
    }) {
        return Err(bad(
            "prestressed round beam tangent/inertia accumulation overflow",
        ));
    }
    Ok(RoundBeamPencil {
        k,
        m,
        nodes_m: nodes,
    })
}

#[cfg(test)]
#[path = "prestress_tests.rs"]
mod tests;
