//! The bar as a TET10 solid: mesh from the card's thickness profile, orthotropic wood, band modes
//! of the free body, and the six rigid coordinates the strike needs. Composition of existing owners
//! (fs-conduction structured grid, fs-solid TET10, fs-modal certified slicing); no new numerics.
use crate::card::Card;
use fs_exec::Cx;
use fs_material::OrthotropicElastic;
use fs_modal::{ModePair, SliceOptions, WindowPlan, slice_windows};
use fs_solid::linear3::tet10::{
    BoundaryProjection, QuadraticTetProblem, TetQuadrature, promote_p1_to_p2,
};
use fs_solid::linear3::{
    TetAssemblyBudget, TetElasticAssembly, TetElasticError, TetElasticMaterial, TetMaterialField,
};

pub const TAU: f64 = core::f64::consts::TAU;
/// Lower edge of the elastic window: above the six rigid modes, below the first elastic one.
const ELASTIC_LOW_HZ: f64 = 20.0;

pub struct Mesh {
    pub nodes: Vec<[f64; 3]>,
    pub tets: Vec<[usize; 4]>,
    pub counts: [usize; 3],
}

/// Mesh point of the reference-cube point `u`: x along the length, y across the width (centred), z
/// through the thickness (centred). The cut removes wood from the lower face only, so the top face
/// stays at z = +T/2 and the lower face sits at z = T/2 - h(x).
pub fn point(card: &Card, u: [f64; 3]) -> [f64; 3] {
    let x = card.length_m * u[0];
    let h = card.profile.thickness(card.length_m, card.thickness_m, x);
    [
        x,
        card.width_m * (u[1] - 0.5),
        0.5 * card.thickness_m - h * (1.0 - u[2]),
    ]
}

/// Inverse of [`point`] on the bar's own geometry.
fn reference(card: &Card, p: [f64; 3]) -> [f64; 3] {
    let h = card
        .profile
        .thickness(card.length_m, card.thickness_m, p[0]);
    [
        p[0] / card.length_m,
        p[1] / card.width_m + 0.5,
        1.0 - (0.5 * card.thickness_m - p[2]) / h,
    ]
}

pub fn mesh(card: &Card, counts: [usize; 3]) -> Mesh {
    let (complex, nodes) = fs_conduction::fixtures::structured_tets(counts, &|i, j, k| {
        point(
            card,
            [
                i as f64 / counts[0] as f64,
                j as f64 / counts[1] as f64,
                k as f64 / counts[2] as f64,
            ],
        )
    });
    Mesh {
        nodes,
        tets: complex.tets.iter().map(|t| t.map(|i| i as usize)).collect(),
        counts,
    }
}

pub fn rung_counts(card: &Card, rung: usize) -> [usize; 3] {
    card.mesh_counts.map(|c| c * rung)
}

/// Boundary mid-edge nodes are snapped back onto the bar's surface: on the curved lower face a straight
/// midpoint would leave a P1 geometry error under the quadratic field. On a prism this is the identity.
struct Surface<'a>(&'a Card);

impl BoundaryProjection for Surface<'_> {
    fn project(&self, p: [f64; 3]) -> Result<[f64; 3], TetElasticError> {
        let mut u = reference(self.0, p);
        let axis = (0..3)
            .min_by(|&a, &b| u[a].min(1.0 - u[a]).total_cmp(&u[b].min(1.0 - u[b])))
            .unwrap_or(0);
        u[axis] = if u[axis] < 0.5 { 0.0 } else { 1.0 };
        Ok(point(self.0, u))
    }
}

/// The card's wood as an oriented orthotropic law. Mesh axes: x = length (L), y = width (T), z = thickness (R).
pub fn material(card: &Card, density: f64) -> Result<TetElasticMaterial, String> {
    let law = OrthotropicElastic::new(card.youngs_pa, card.poisson, card.shear_pa, 1e-3)
        .map_err(|e| format!("bar card wood is not an admissible orthotropic solid: {e}"))?;
    let identity = fs_blake3::hash_domain(
        "org.frankensim.bar-strike.material.v1",
        format!("{}\ndensity={density}", card.text).as_bytes(),
    );
    TetElasticMaterial::try_new_oriented_orthotropic(
        density,
        &law,
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        identity.0,
        fs_blake3::hash_domain("org.frankensim.bar-strike.axes.v1", b"x=L,y=T,z=R").0,
    )
    .map_err(|e| format!("oriented orthotropic material: {e}"))
}

pub fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = fs_exec::CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 270,
                kernel_id: 270,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        f(&cx)
    })
}

/// The P1 mesh promoted to TET10. Promotion keeps the P1 node ids first, so the strike node, the
/// classifier and the BEM vertex map read the same vertex DOFs as on P1. Returns all nodes, the element
/// count and the free-free assembly.
pub fn assemble(
    card: &Card,
    mesh: &Mesh,
    mat: &TetElasticMaterial,
    cx: &Cx<'_>,
) -> Result<(Vec<[f64; 3]>, usize, TetElasticAssembly), String> {
    let (nodes, elements) = promote_p1_to_p2(&mesh.nodes, &mesh.tets, Some(&Surface(card)), cx)
        .map_err(|e| format!("P1 to TET10 promotion: {e:?}"))?;
    let a = QuadraticTetProblem {
        nodes_m: &nodes,
        elements: &elements,
        materials: TetMaterialField::Uniform(mat),
        fixed_dofs: &[],
        budget: TetAssemblyBudget::standard(),
        // Affine minima: stiffness degree 2, consistent mass degree 4.
        quadrature: TetQuadrature {
            stiffness_points_per_dir: 3,
            mass_points_per_dir: 4,
        },
    }
    .assemble(cx)
    .map_err(|e| format!("free-free TET10 assembly: {e:?}"))?;
    Ok((nodes, elements.len(), a))
}

fn slice_options() -> SliceOptions {
    let mut o = SliceOptions {
        max_lanczos: 2048,
        max_restarts: 48,
        ritz_tol: 1e-13,
        mass_diagonal_equilibration: true,
        ..SliceOptions::default()
    };
    // A window holds a few modes: a bounded Krylov budget per window, and a 1 GiB factor cap each.
    o.max_lanczos = 512;
    o.ldlt.max_factor_bytes = Some(1 << 30);
    o
}

/// Every elastic mode of the free body from 20 Hz to the card's band, in ordered parallel windows.
///
/// A window edge that falls within the certificate's half-width of an eigenvalue leaves that window
/// unresolved (`WindowUnresolved`, seen on the A5 check cards). A different window count moves every edge,
/// so the card's count is tried first and then the next three. Each attempt is a complete certified window
/// plan, so a success is as good as any other; if all four refuse, the last refusal is returned.
pub fn band_modes(
    card: &Card,
    a: &TetElasticAssembly,
    cx: &Cx<'_>,
) -> Result<Vec<ModePair>, String> {
    let mut last = String::new();
    for n in card.windows..card.windows + 4 {
        match band_modes_in(card, a, cx, n) {
            Ok(modes) => return Ok(modes),
            Err(e) => last = e,
        }
    }
    Err(last)
}

fn band_modes_in(
    card: &Card,
    a: &TetElasticAssembly,
    cx: &Cx<'_>,
    n: usize,
) -> Result<Vec<ModePair>, String> {
    let edges: Vec<f64> = (0..=n)
        .map(|i| ELASTIC_LOW_HZ + (card.band_hz - ELASTIC_LOW_HZ) * i as f64 / n as f64)
        .collect();
    let plan = WindowPlan {
        windows: edges
            .windows(2)
            .map(|w| ((TAU * w[0]).powi(2), (TAU * w[1]).powi(2)))
            .collect(),
        threads: n,
        maximum_windows: n,
        maximum_bytes: 16 << 30,
        maximum_relative_halfwidth: 1e-8,
        coordinates: None,
    };
    let report = slice_windows(&a.stiffness, &a.mass, &plan, &slice_options(), cx)
        .map_err(|e| format!("certified band windows: {e:?}"))?;
    let mut modes: Vec<ModePair> = report.reports.into_iter().flat_map(|w| w.modes).collect();
    canonical_phase(&mut modes);
    Ok(modes)
}

pub fn hz(m: &ModePair) -> f64 {
    m.lambda.max(0.0).sqrt() / TAU
}

pub fn hz_error(m: &ModePair) -> f64 {
    let f = hz(m);
    ((m.interval.1.max(0.0).sqrt() / TAU - f).abs())
        .max((f - m.interval.0.max(0.0).sqrt() / TAU).abs())
}

/// Canonical phase for reproducible renders.
fn canonical_phase(modes: &mut [ModePair]) {
    for m in modes {
        let at = m
            .phi
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.abs().total_cmp(&b.1.abs()))
            .map_or(0, |x| x.0);
        if m.phi[at] < 0.0 {
            for v in &mut m.phi {
                *v = -*v;
            }
        }
    }
}

/// The grid node on the top face at the card's strike station, mid-width.
pub fn strike_node(card: &Card, mesh: &Mesh) -> usize {
    let [nx, ny, nz] = mesh.counts;
    let i = (card.strike_fraction * nx as f64).round() as usize;
    (i * (ny + 1) + ny / 2) * (nz + 1) + nz
}

/// Station of [`strike_node`] along the bar, in metres.
pub fn strike_x(card: &Card, mesh: &Mesh) -> f64 {
    mesh.nodes[strike_node(card, mesh)][0]
}

/// Vertical (through-thickness) bending modes, picked by shape: the cross-section mean displacement
/// along the thickness axis dominates. Frequencies alone never label a mode.
pub fn vertical_modes(mesh: &Mesh, modes: &[ModePair]) -> Vec<usize> {
    let [_, ny, nz] = mesh.counts;
    let end = (ny + 1) * (nz + 1);
    let mut vertical = Vec::new();
    for (i, m) in modes.iter().enumerate() {
        let mut mean = [0.0; 3];
        let mut rms = 0.0;
        let mut weight = 0.0;
        for node in 0..end {
            let j = node / (nz + 1);
            let k = node % (nz + 1);
            let w = (if j == 0 || j == ny { 0.5 } else { 1.0 })
                * (if k == 0 || k == nz { 0.5 } else { 1.0 });
            let u: [f64; 3] = core::array::from_fn(|c| m.phi[3 * node + c]);
            for c in 0..3 {
                mean[c] += w * u[c];
            }
            rms += w * u.iter().map(|v| v * v).sum::<f64>();
            weight += w;
        }
        mean = mean.map(|v| v / weight);
        let amplitude = (rms / weight).sqrt();
        if mean[2].abs() > 0.3 * amplitude && mean[2].abs() > 2.0 * mean[1].abs() {
            vertical.push(i);
        }
    }
    vertical
}

/// Preserve all six free-body coordinates in contact: dropping them would spuriously anchor a nominally
/// free-free bar. The analytic affine fields are M-orthonormalized on the assembled body, cut included.
/// `nodes` are the assembly's nodes (the TET10 mid-edge nodes included). Returns, per rigid coordinate,
/// its (negated) displacement at the strike node along the thickness axis.
pub fn rigid_projectors(
    card: &Card,
    mesh: &Mesh,
    nodes: &[[f64; 3]],
    a: &TetElasticAssembly,
) -> Result<Vec<f64>, String> {
    let mut basis: Vec<Vec<f64>> = Vec::new();
    let mut mass_basis: Vec<Vec<f64>> = Vec::new();
    let mut scratch = vec![0.0; a.free_dofs.len()];
    let mut projections = Vec::new();
    let node = strike_node(card, mesh);
    for axis in 0..6 {
        let mut field: Vec<f64> = nodes
            .iter()
            .flat_map(|&[x, y, z]| match axis {
                0 => [1.0, 0.0, 0.0],
                1 => [0.0, 1.0, 0.0],
                2 => [0.0, 0.0, 1.0],
                3 => [0.0, -z, y],
                4 => [z, 0.0, -x],
                _ => [-y, x, 0.0],
            })
            .collect();
        for _ in 0..2 {
            for (b, mb) in basis.iter().zip(&mass_basis) {
                let dot: f64 = field.iter().zip(mb).map(|(x, y)| x * y).sum();
                for (x, y) in field.iter_mut().zip(b) {
                    *x -= dot * y;
                }
            }
        }
        a.mass.spmv(&field, &mut scratch);
        let norm = field
            .iter()
            .zip(&scratch)
            .map(|(x, y)| x * y)
            .sum::<f64>()
            .sqrt();
        if !(norm.is_finite() && norm > 0.0) {
            return Err(format!("rigid coordinate {axis} has no mass norm"));
        }
        for x in &mut field {
            *x /= norm;
        }
        a.mass.spmv(&field, &mut scratch);
        let check: f64 = field.iter().zip(&scratch).map(|(x, y)| x * y).sum();
        if (check - 1.0).abs() >= 1e-10 {
            return Err(format!(
                "rigid coordinate {axis} is not M-orthonormal: {check}"
            ));
        }
        projections.push(-field[3 * node + 2]);
        basis.push(field);
        mass_basis.push(scratch.clone());
    }
    Ok(projections)
}
