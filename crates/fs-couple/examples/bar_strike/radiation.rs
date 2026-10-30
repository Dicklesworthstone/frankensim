//! Burton-Miller receiver transfer of every band mode, each at its own frequency (the narrow-band
//! modal observer of the bar_s1 test). The bar's boundary mesh is the structural grid's own boundary
//! (or a nested coarser one, `bem_stride`); the vertex displacements of each mode give its normal
//! velocity on the panels. Rigid coordinates get no transfer (direct rigid-body radiation is not modelled).
use crate::bar::{Mesh, TAU, hz, mesh as build_mesh};
use crate::card::Card;
use fs_bem::helmholtz::{
    Formulation, HelmholtzError, Medium, exterior_pressure_at_points, solve_radiation_batch,
};
use fs_bem::panel3d::SpherePanels;
use fs_math::c64::C64;
use fs_modal::ModePair;
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

/// Outward-oriented boundary triangles of a structured grid, as grid-node index triples.
fn boundary(mesh: &Mesh) -> Vec<[usize; 3]> {
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
                    for mut face in [[a, b, c], [a, c, d]] {
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
    let stride = card.bem_stride;
    let counts = finest.counts.map(|c| c / stride);
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
            ((i * stride) * (fine[1] + 1) + j * stride) * (fine[2] + 1) + k * stride
        })
        .collect();
    let shapes = modes
        .iter()
        .map(|mode| {
            faces
                .iter()
                .zip(surface.normals())
                .map(|(face, normal)| {
                    face.iter()
                        .map(|&v| {
                            (0..3)
                                .map(|c| mode.phi[3 * map[v] + c] * normal[c])
                                .sum::<f64>()
                        })
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
    let results: Vec<Mutex<Option<Result<(C64, f64), String>>>> =
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
                            .map(|p| (p[0], s[0].panels_per_wavelength))
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
    }
    Ok(Transfers {
        pressure,
        panels: surface.centroids().len(),
        min_panels_per_wavelength: min_ppw,
    })
}
