//! Burton-Miller receiver transfer of every band mode, each at its own frequency (the narrow-band
//! modal observer of the bar_s1 test). The bar's boundary mesh is the structural grid's own boundary
//! (or a nested coarser one, `bem_stride`); the vertex displacements of each mode give its normal
//! velocity on the panels. Rigid coordinates get no transfer here: a zero frequency has no narrow-band
//! hold. The observer's `observer_rigid_rows` row is the only place they are heard.
use crate::bar::{Mesh, TAU, hz, mesh as build_mesh};
use crate::card::{Card, Tube};
use crate::tube::{IsolatedMode, MouthField, mouth_centre, mouth_points};
use fs_bem::helmholtz::{
    Formulation, HelmholtzError, Medium, exterior_pressure_at_points, solve_radiation_batch,
};
use fs_bem::panel3d::SpherePanels;
use fs_math::{c64::C64, det};
use fs_modal::ModePair;
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

/// Outward-oriented boundary triangles of a structured grid, as grid-node index triples. Each grid face
/// is split along its shorter diagonal: under the cut a side face is sheared, and fs-bem resolves a
/// wavelength along a panel's longest edge, so the longer diagonal would only lower the floor. A tie
/// (every rectangular face, so every face of an uncut bar) keeps the a-c split.
pub(crate) fn boundary(mesh: &Mesh) -> Vec<[usize; 3]> {
    split_boundary(mesh, true)
}

/// The split before the shorter-diagonal rule, every face along a-c, for the tests that compare.
#[cfg(test)]
pub(crate) fn boundary_along_ac(mesh: &Mesh) -> Vec<[usize; 3]> {
    split_boundary(mesh, false)
}

fn split_boundary(mesh: &Mesh, shorter: bool) -> Vec<[usize; 3]> {
    let squared = |u: usize, v: usize| {
        let (p, q) = (mesh.nodes[u], mesh.nodes[v]);
        (0..3).map(|k| (p[k] - q[k]) * (p[k] - q[k])).sum::<f64>()
    };
    let counts = mesh.counts;
    let index = |q: [usize; 3]| (q[0] * (counts[1] + 1) + q[1]) * (counts[2] + 1) + q[2];
    let mut faces = Vec::new();
    for axis in 0..3 {
        let other: Vec<_> = (0..3).filter(|&a| a != axis).collect();
        for high in [false, true] {
            for i in 0..counts[other[0]] {
                for j in 0..counts[other[1]] {
                    let mut q = [0; 3];
                    q[axis] = if high { counts[axis] } else { 0 };
                    q[other[0]] = i;
                    q[other[1]] = j;
                    let a = index(q);
                    q[other[0]] += 1;
                    let b = index(q);
                    q[other[1]] += 1;
                    let c = index(q);
                    q[other[0]] -= 1;
                    let d = index(q);
                    let pair = if shorter && squared(b, d) < squared(a, c) {
                        [[a, b, d], [b, c, d]]
                    } else {
                        [[a, b, c], [a, c, d]]
                    };
                    for mut face in pair {
                        let p = face.map(|v| mesh.nodes[v]);
                        let ab: [f64; 3] = core::array::from_fn(|k| p[1][k] - p[0][k]);
                        let ac: [f64; 3] = core::array::from_fn(|k| p[2][k] - p[0][k]);
                        let normal_component = ab[(axis + 1) % 3] * ac[(axis + 2) % 3]
                            - ab[(axis + 2) % 3] * ac[(axis + 1) % 3];
                        if (normal_component > 0.0) != high {
                            face.swap(1, 2);
                        }
                        faces.push(face);
                    }
                }
            }
        }
    }
    faces
}

/// The BEM surface and, per mode, the mode's normal velocity on each panel (unit generalized modal
/// velocity). The surface nodes are a nested subset of the finest structural grid, so no nearest-point
/// search touches the mode shapes.
pub fn surface_and_shapes(
    card: &Card,
    finest: &Mesh,
    modes: &[ModePair],
) -> Result<(SpherePanels, Vec<Vec<f64>>), String> {
    surface_and_shapes_at(card, finest, modes, [card.bem_stride; 3])
}

/// The same with its own stride along each axis.
pub fn surface_and_shapes_at(
    card: &Card,
    finest: &Mesh,
    modes: &[ModePair],
    stride: [usize; 3],
) -> Result<(SpherePanels, Vec<Vec<f64>>), String> {
    let fields: Vec<&[f64]> = modes.iter().map(|m| m.phi.as_slice()).collect();
    surface_and_fields_at(card, finest, &fields, stride)
}

/// The same for any fields laid out as a mode's `phi`.
pub fn surface_and_fields_at(
    card: &Card,
    finest: &Mesh,
    fields: &[&[f64]],
    stride: [usize; 3],
) -> Result<(SpherePanels, Vec<Vec<f64>>), String> {
    let counts: [usize; 3] = core::array::from_fn(|a| finest.counts[a] / stride[a]);
    let coarse = build_mesh(card, counts);
    let faces = boundary(&coarse);
    let surface =
        SpherePanels::from_triangles(faces.iter().map(|f| f.map(|v| coarse.nodes[v])).collect())
            .map_err(|e| format!("BEM surface: {e:?}"))?;
    let fine = finest.counts;
    let map: Vec<usize> = (0..coarse.nodes.len())
        .map(|node| {
            let k = node % (counts[2] + 1);
            let j = node / (counts[2] + 1) % (counts[1] + 1);
            let i = node / ((counts[1] + 1) * (counts[2] + 1));
            ((i * stride[0]) * (fine[1] + 1) + j * stride[1]) * (fine[2] + 1) + k * stride[2]
        })
        .collect();
    let shapes = fields
        .iter()
        .map(|phi| {
            faces
                .iter()
                .zip(surface.normals())
                .map(|(face, normal)| {
                    face.iter()
                        .map(|&v| (0..3).map(|c| phi[3 * map[v] + c] * normal[c]).sum::<f64>())
                        .sum::<f64>()
                        / 3.0
                })
                .collect()
        })
        .collect();
    Ok((surface, shapes))
}

pub struct Transfers {
    /// Receiver pressure per unit generalized modal velocity, exp(-iwt), Pa s/(m sqrt(kg)).
    pub pressure: Vec<C64>,
    /// Each mode's own radiated power per unit generalized modal velocity, the same solve's
    /// `radiated_power` (1/2 Re SUM p conj(v) A over the bar's panels). With mass-normalized shapes it is the
    /// mode's amplitude decay rate from radiation, in 1/s. Below the BEM's resolution a weak radiator comes
    /// out negative (about -5e-3/s on a 323 mm aluminium bar's fundamental, whose true value is near 1e-6/s),
    /// so `render::bar_modes` adds max(W, 0) and only under the card's `bar_self_radiation on`.
    pub self_radiation: Vec<f64>,
    pub panels: usize,
    /// Smallest panels-per-wavelength over the modes (the BEM refuses below its floor).
    pub min_panels_per_wavelength: f64,
}

/// Solve every mode at its own frequency, in parallel; the first refusal is returned with the mode's frequency.
pub fn modal_transfers(
    card: &Card,
    medium: Medium,
    surface: &SpherePanels,
    shapes: &[Vec<f64>],
    modes: &[ModePair],
) -> Result<Transfers, String> {
    let next = AtomicUsize::new(0);
    let results: Vec<Mutex<Option<Result<(C64, f64, f64), String>>>> =
        modes.iter().map(|_| Mutex::new(None)).collect();
    let threads = std::thread::available_parallelism()
        .map_or(4, |n| n.get())
        .min(16);
    std::thread::scope(|scope| {
        for _ in 0..threads {
            scope.spawn(|| {
                loop {
                    let i = next.fetch_add(1, Ordering::Relaxed);
                    if i >= modes.len() {
                        break;
                    }
                    let f = hz(&modes[i]);
                    let velocity: Vec<C64> = shapes[i].iter().map(|&v| C64::from_re(v)).collect();
                    let solved = solve_radiation_batch(
                        surface,
                        TAU * f / medium.sound_speed,
                        medium,
                        &[velocity.as_slice()],
                        Formulation::BurtonMiller,
                    );
                    let out = match solved {
                        Ok(s) => exterior_pressure_at_points(surface, &s[0], medium, &[card.receiver_m])
                            .map(|p| (p[0], s[0].panels_per_wavelength, s[0].radiated_power))
                            .map_err(|e| format!("receiver pressure at {f:.1} Hz: {e:?}")),
                        Err(HelmholtzError::TooCoarse { .. }) => Err(format!(
                            "BEM surface is too coarse at {f:.1} Hz: lower band_hz or refine the mesh"
                        )),
                        Err(e) => Err(format!("BEM solve at {f:.1} Hz: {e:?}")),
                    };
                    if let Ok(mut slot) = results[i].lock() {
                        *slot = Some(out);
                    }
                }
            });
        }
    });
    let mut pressure = Vec::with_capacity(modes.len());
    let mut self_radiation = Vec::with_capacity(modes.len());
    let mut min_ppw = f64::INFINITY;
    for slot in results {
        let solved = slot
            .into_inner()
            .map_err(|_| "poisoned BEM result".to_string())?
            .ok_or_else(|| "a mode was never solved".to_string())??;
        if !(solved.0.re.is_finite() && solved.0.im.is_finite()) {
            return Err("non-finite receiver pressure".into());
        }
        min_ppw = min_ppw.min(solved.1);
        pressure.push(solved.0);
        self_radiation.push(solved.2);
    }
    Ok(Transfers {
        pressure,
        self_radiation,
        panels: surface.centroids().len(),
        min_panels_per_wavelength: min_ppw,
    })
}

/// Normal velocity on each panel that a rigid bar must answer with, per unit outward volume velocity spread over
/// `sources`: the negative of the free field's velocity from a source p = -i w rho exp(ikr) / (4 pi r), exp(-iwt).
pub fn reflected_velocity(surface: &SpherePanels, sources: &[[f64; 3]], k: f64) -> Vec<C64> {
    let share = 1.0 / sources.len() as f64;
    surface
        .centroids()
        .iter()
        .zip(surface.normals())
        .map(|(x, n)| {
            let mut sum = C64::ZERO;
            for s in sources {
                let d: [f64; 3] = core::array::from_fn(|c| x[c] - s[c]);
                let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
                let cos_n = (d[0] * n[0] + d[1] * n[1] + d[2] * n[2]) / r;
                let travel = C64::new(det::cos(k * r), det::sin(k * r));
                let spread = C64::new(-1.0 / r, k);
                sum = sum + travel * spread * C64::from_re(cos_n * share / (2.0 * TAU * r));
            }
            sum
        })
        .collect()
}

/// At each tube mode's frequency: the mouth-mean pressure per unit modal velocity of every bar mode, and the
/// mouth-mean pressure the rigid bar reflects per unit mouth flow. One factorization per tube mode, in parallel.
pub fn mouth_fields(
    card: &Card,
    tube: &Tube,
    medium: Medium,
    surface: &SpherePanels,
    shapes: &[Vec<f64>],
    isolated: &[IsolatedMode],
) -> Result<Vec<MouthField>, String> {
    let points = mouth_points(mouth_centre(card, tube), tube.radius_m);
    let share = 1.0 / points.len() as f64;
    let velocities: Vec<Vec<C64>> = shapes
        .iter()
        .map(|shape| shape.iter().map(|&v| C64::from_re(v)).collect())
        .collect();
    let next = AtomicUsize::new(0);
    let results: Vec<Mutex<Option<Result<MouthField, String>>>> =
        isolated.iter().map(|_| Mutex::new(None)).collect();
    let threads = std::thread::available_parallelism()
        .map_or(4, std::num::NonZero::get)
        .min(16);
    std::thread::scope(|scope| {
        for _ in 0..threads {
            scope.spawn(|| {
                loop {
                    let m = next.fetch_add(1, Ordering::Relaxed);
                    if m >= isolated.len() {
                        break;
                    }
                    let f = isolated[m].omega_rad_s / TAU;
                    let k = isolated[m].omega_rad_s / medium.sound_speed;
                    let reflected = reflected_velocity(surface, &points, k);
                    let mut fields: Vec<&[C64]> = velocities.iter().map(Vec::as_slice).collect();
                    fields.push(&reflected);
                    let out = solve_radiation_batch(
                        surface,
                        k,
                        medium,
                        &fields,
                        Formulation::BurtonMiller,
                    )
                    .map_err(|e| format!("tube mode at {f:.1} Hz: BEM solve: {e:?}"))
                    .and_then(|solved| {
                        let mut mean = Vec::with_capacity(solved.len());
                        for s in &solved {
                            let p = exterior_pressure_at_points(surface, s, medium, &points)
                                .map_err(|e| format!("mouth pressure at {f:.1} Hz: {e:?}"))?;
                            let sum = p.iter().fold(C64::ZERO, |a, b| a + *b);
                            mean.push(sum * C64::from_re(share));
                        }
                        let reflection = mean.pop().ok_or("no reflection field")?;
                        Ok(MouthField {
                            per_bar_mode: mean,
                            reflection,
                        })
                    });
                    if let Ok(mut slot) = results[m].lock() {
                        *slot = Some(out);
                    }
                }
            });
        }
    });
    results
        .into_iter()
        .map(|slot| {
            slot.into_inner()
                .map_err(|_| "poisoned BEM result".to_string())?
                .ok_or_else(|| "a tube mode was never solved".to_string())?
        })
        .collect()
}
