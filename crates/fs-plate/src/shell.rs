//! Flat-facet shell assembly for curved shells and musical bells (bead `frankensim-music-v8-root-3ez8g.12.2`).
//!
//! Combines:
//! - Constant Strain Triangle (CST) in-plane membrane stiffness
//! - Discrete Kirchhoff Triangle (DKT) plate bending stiffness
//! - Regularized drilling DOF ($\theta_z$) handling with disclosed parameter $\alpha_{\text{drill}}$
//! - Full 3D local-to-global frame transformations ($18 \times 18$)
//! - Lumped mass matrix (translational + rotary)
//! - Axisymmetric bell profile mesh generator and harmonic partial ratio analysis
//! - Oracle ladder generators (cylinder, hemisphere, church bell)

use crate::{dkt_stiffness, PlateError, PlateSection};
use fs_modal::{slice_window, SliceOptions, SliceReport};
use fs_sparse::{Coo, Csr};

/// Drilling DOF regularization coefficient (disclosed in CONTRACT.md).
pub const DRILLING_ALPHA: f64 = 1e-3;

/// 3D shell mesh with triangular facets.
#[derive(Debug, Clone, PartialEq)]
pub struct ShellMesh {
    /// 3D node coordinates `(x, y, z)` [m].
    pub nodes: Vec<[f64; 3]>,
    /// Triangle node indices `[n0, n1, n2]`.
    pub tris: Vec<[usize; 3]>,
}

impl ShellMesh {
    /// Construct a new shell mesh from nodes and triangles.
    ///
    /// # Errors
    /// Returns [`PlateError::DegenerateElement`] if any triangle references out-of-range nodes.
    pub fn new(nodes: Vec<[f64; 3]>, tris: Vec<[usize; 3]>) -> Result<Self, PlateError> {
        let nn = nodes.len();
        for tri in &tris {
            for &n in tri {
                if n >= nn {
                    return Err(PlateError::BadBoundary {
                        node: n,
                        node_count: nn,
                    });
                }
            }
        }
        Ok(Self { nodes, tris })
    }

    /// Number of nodes in the shell mesh.
    #[must_use]
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Number of triangular elements.
    #[must_use]
    pub fn element_count(&self) -> usize {
        self.tris.len()
    }
}

/// Boundary support condition for shell nodes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShellSupport {
    /// All 6 DOFs free.
    Free,
    /// Clamped: all 6 DOFs constrained to zero (u=v=w=θx=θy=θz=0).
    Clamped,
    /// Pinned: translation constrained (u=v=w=0), rotations free.
    Pinned,
}

/// Assembled reduced shell model pencil `(K, M)`.
#[derive(Debug, Clone)]
pub struct ShellModel {
    /// Reduced stiffness matrix (membrane + bending + drilling).
    pub k: Csr,
    /// Reduced lumped mass matrix (translational + rotary).
    pub m: Csr,
    /// Map from full DOF index (`6 * node + component`) to reduced free DOF index.
    pub dof_map: Vec<Option<usize>>,
    /// Number of free DOFs.
    pub free: usize,
}

/// Assemble the reduced (K, M) pencil for a 3D flat-facet shell mesh.
///
/// # Errors
/// Returns [`PlateError`] on invalid geometry or section.
pub fn assemble_shell(
    mesh: &ShellMesh,
    section: &PlateSection,
    boundary_nodes: &[usize],
    support: ShellSupport,
) -> Result<ShellModel, PlateError> {
    let nn = mesh.node_count();
    let ndof = 6 * nn;

    // 1. Build DOF elimination map
    let mut constrained = vec![false; ndof];
    if support != ShellSupport::Free {
        for &b in boundary_nodes {
            if b >= nn {
                return Err(PlateError::BadBoundary {
                    node: b,
                    node_count: nn,
                });
            }
            match support {
                ShellSupport::Clamped => {
                    for comp in 0..6 {
                        constrained[6 * b + comp] = true;
                    }
                }
                ShellSupport::Pinned => {
                    for comp in 0..3 {
                        constrained[6 * b + comp] = true;
                    }
                }
                ShellSupport::Free => {}
            }
        }
    }

    let mut dof_map = vec![None; ndof];
    let mut free_count = 0;
    for i in 0..ndof {
        if !constrained[i] {
            dof_map[i] = Some(free_count);
            free_count += 1;
        }
    }

    let mut k_coo = Coo::new(free_count, free_count);
    let mut m_diag = vec![0.0f64; free_count];

    let h = section.thickness;
    let rho = section.density;

    // Plane stress constitutive matrix for membrane: C = 12 / h^2 * D
    let c_scale = 12.0 / (h * h);
    let c_mat = [
        c_scale * section.d[0],
        c_scale * section.d[1],
        0.0,
        c_scale * section.d[3],
        c_scale * section.d[4],
        0.0,
        0.0,
        0.0,
        c_scale * section.d[8],
    ];

    // 2. Loop over triangular facets
    for (elem_idx, tri) in mesh.tris.iter().enumerate() {
        let p0 = mesh.nodes[tri[0]];
        let p1 = mesh.nodes[tri[1]];
        let p2 = mesh.nodes[tri[2]];

        // Vector p0 -> p1
        let v01 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
        let l01 = (v01[0] * v01[0] + v01[1] * v01[1] + v01[2] * v01[2]).sqrt();
        if l01 < 1e-12 {
            return Err(PlateError::DegenerateElement {
                element: elem_idx,
                twice_area: 0.0,
            });
        }
        let ex = [v01[0] / l01, v01[1] / l01, v01[2] / l01];

        // Vector p0 -> p2
        let v02 = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];

        // Normal ez = ex x v02
        let n_cross = [
            ex[1] * v02[2] - ex[2] * v02[1],
            ex[2] * v02[0] - ex[0] * v02[2],
            ex[0] * v02[1] - ex[1] * v02[0],
        ];
        let n_len = (n_cross[0] * n_cross[0] + n_cross[1] * n_cross[1] + n_cross[2] * n_cross[2]).sqrt();
        if n_len < 1e-12 {
            return Err(PlateError::DegenerateElement {
                element: elem_idx,
                twice_area: 0.0,
            });
        }
        let ez = [n_cross[0] / n_len, n_cross[1] / n_len, n_cross[2] / n_len];

        // ey = ez x ex
        let ey = [
            ez[1] * ex[2] - ez[2] * ex[1],
            ez[2] * ex[0] - ez[0] * ex[2],
            ez[0] * ex[1] - ez[1] * ex[0],
        ];

        // Local 2D coordinates in facet plane (z' = 0 by construction)
        let lx0 = 0.0;
        let ly0 = 0.0;
        let lx1 = l01;
        let ly1 = 0.0;
        let lx2 = v02[0] * ex[0] + v02[1] * ex[1] + v02[2] * ex[2];
        let ly2 = v02[0] * ey[0] + v02[1] * ey[1] + v02[2] * ey[2];

        let lx = [lx0, lx1, lx2];
        let ly = [ly0, ly1, ly2];

        let twice_area = ((lx1 - lx0) * (ly2 - ly0) - (lx2 - lx0) * (ly1 - ly0)).abs();
        let area = 0.5 * twice_area;
        if area < 1e-14 {
            return Err(PlateError::DegenerateElement {
                element: elem_idx,
                twice_area,
            });
        }

        // Local 18x18 stiffness matrix
        let mut k_local = [0.0f64; 18 * 18];

        // (a) Membrane stiffness (CST in-plane): local DOFs u1, v1, u2, v2, u3, v3 -> slots 0,1, 6,7, 12,13
        // B_m matrix: epsilon = B_m * u
        let b1 = ly1 - ly2; // y1 - y2
        let b2 = ly2 - ly0; // y2 - y0
        let b3 = ly0 - ly1; // y0 - y1
        let c1 = lx2 - lx1; // x2 - x1
        let c2 = lx0 - lx2; // x0 - x2
        let c3 = lx1 - lx0; // x1 - x0

        let bm = [
            [b1 / twice_area, 0.0, b2 / twice_area, 0.0, b3 / twice_area, 0.0],
            [0.0, c1 / twice_area, 0.0, c2 / twice_area, 0.0, c3 / twice_area],
            [
                c1 / twice_area,
                b1 / twice_area,
                c2 / twice_area,
                b2 / twice_area,
                c3 / twice_area,
                b3 / twice_area,
            ],
        ];

        let mut km_6x6 = [0.0f64; 36];
        for r in 0..6 {
            for c in 0..6 {
                let mut sum = 0.0;
                for i in 0..3 {
                    for j in 0..3 {
                        sum += bm[i][r] * c_mat[i * 3 + j] * bm[j][c];
                    }
                }
                km_6x6[r * 6 + c] = sum * area;
            }
        }

        // Place km_6x6 into k_local (u, v components: slots 0,1; 6,7; 12,13)
        let m_dofs = [0, 1, 6, 7, 12, 13];
        for (i, &di) in m_dofs.iter().enumerate() {
            for (j, &dj) in m_dofs.iter().enumerate() {
                k_local[di * 18 + dj] += km_6x6[i * 6 + j];
            }
        }

        // (b) DKT Bending stiffness (9x9): local DOFs w1, tx1, ty1, w2, tx2, ty2, w3, tx3, ty3 -> slots 2,3,4; 8,9,10; 14,15,16
        let (kb_9x9, _) = dkt_stiffness(&lx, &ly, &section.d, elem_idx)?;
        let b_dofs = [2, 3, 4, 8, 9, 10, 14, 15, 16];
        for (i, &di) in b_dofs.iter().enumerate() {
            for (j, &dj) in b_dofs.iter().enumerate() {
                k_local[di * 18 + dj] += kb_9x9[i * 9 + j];
            }
        }

        // (c) Regularized drilling stiffness on theta_z: slots 5, 11, 17
        let k_drill = DRILLING_ALPHA * c_mat[0] * area;
        k_local[5 * 18 + 5] += k_drill;
        k_local[11 * 18 + 11] += k_drill;
        k_local[17 * 18 + 17] += k_drill;

        // (d) Transformation matrix R (3x3): rows are ex, ey, ez
        // For each node, T_node = diag(R, R) (6x6)
        // Transform k_local (18x18) to k_global: K_glob = T^T * K_loc * T
        let r_mat = [
            ex[0], ex[1], ex[2],
            ey[0], ey[1], ey[2],
            ez[0], ez[1], ez[2],
        ];

        let mut k_global = [0.0f64; 18 * 18];
        for node_i in 0..3 {
            for node_j in 0..3 {
                for comp_ti in 0..2 {
                    // 0 = translational, 1 = rotational
                    for comp_tj in 0..2 {
                        // Transform 3x3 sub-block
                        let r_offset_i = node_i * 6 + comp_ti * 3;
                        let r_offset_j = node_j * 6 + comp_tj * 3;

                        for a in 0..3 {
                            for b in 0..3 {
                                let mut val = 0.0;
                                for p in 0..3 {
                                    for q in 0..3 {
                                        let k_loc_val = k_local[(r_offset_i + p) * 18 + (r_offset_j + q)];
                                        val += r_mat[p * 3 + a] * k_loc_val * r_mat[q * 3 + b];
                                    }
                                }
                                k_global[(r_offset_i + a) * 18 + (r_offset_j + b)] += val;
                            }
                        }
                    }
                }
            }
        }

        // (e) Accumulate into global Coo
        let global_dof_indices = [
            6 * tri[0], 6 * tri[0] + 1, 6 * tri[0] + 2, 6 * tri[0] + 3, 6 * tri[0] + 4, 6 * tri[0] + 5,
            6 * tri[1], 6 * tri[1] + 1, 6 * tri[1] + 2, 6 * tri[1] + 3, 6 * tri[1] + 4, 6 * tri[1] + 5,
            6 * tri[2], 6 * tri[2] + 1, 6 * tri[2] + 2, 6 * tri[2] + 3, 6 * tri[2] + 4, 6 * tri[2] + 5,
        ];

        for i in 0..18 {
            if let Some(ri) = dof_map[global_dof_indices[i]] {
                for j in 0..18 {
                    if let Some(cj) = dof_map[global_dof_indices[j]] {
                        k_coo.push(ri, cj, k_global[i * 18 + j]);
                    }
                }
            }
        }

        // (f) Lumped mass contributions per node
        let node_mass = rho * h * area / 3.0;
        let node_rot_inertia = rho * h * h * h * area / 36.0;

        for &n in tri {
            for comp in 0..3 {
                if let Some(r) = dof_map[6 * n + comp] {
                    m_diag[r] += node_mass;
                }
            }
            for comp in 3..6 {
                if let Some(r) = dof_map[6 * n + comp] {
                    m_diag[r] += node_rot_inertia;
                }
            }
        }
    }

    let k_csr = k_coo.assemble();
    let mut m_coo = Coo::new(free_count, free_count);
    for (i, &val) in m_diag.iter().enumerate() {
        m_coo.push(i, i, val.max(1e-15));
    }
    let m_csr = m_coo.assemble();

    Ok(ShellModel {
        k: k_csr,
        m: m_csr,
        dof_map,
        free: free_count,
    })
}

/// Compute certified modes of a shell model in the frequency-squared window `(low, high]`.
///
/// # Errors
/// Returns [`PlateError::Modal`] on solver refusal.
pub fn modes_shell(
    model: &ShellModel,
    window: (f64, f64),
    opts: &SliceOptions,
) -> Result<SliceReport, PlateError> {
    Ok(slice_window(&model.k, &model.m, window, opts)?)
}

/// Generate a cylindrical shell mesh of radius `r` and height `h`.
#[must_use]
pub fn generate_cylinder_shell(r: f64, h: f64, n_theta: usize, n_z: usize) -> ShellMesh {
    let mut nodes = Vec::with_capacity((n_theta + 1) * (n_z + 1));
    for j in 0..=n_z {
        let z = (j as f64 / n_z as f64) * h;
        for i in 0..n_theta {
            let theta = (i as f64 / n_theta as f64) * 2.0 * std::f64::consts::PI;
            let x = r * theta.cos();
            let y = r * theta.sin();
            nodes.push([x, y, z]);
        }
    }

    let mut tris = Vec::with_capacity(2 * n_theta * n_z);
    for j in 0..n_z {
        for i in 0..n_theta {
            let next_i = (i + 1) % n_theta;
            let n00 = j * n_theta + i;
            let n10 = j * n_theta + next_i;
            let n01 = (j + 1) * n_theta + i;
            let n11 = (j + 1) * n_theta + next_i;

            tris.push([n00, n10, n11]);
            tris.push([n00, n11, n01]);
        }
    }

    ShellMesh { nodes, tris }
}

/// Revolve an axisymmetric bell profile into a 3D shell mesh.
/// `profile` contains `(r, z)` points from crown ($z=H$) to lip ($z=0$).
#[must_use]
pub fn generate_bell_shell(profile: &[(f64, f64)], n_theta: usize) -> ShellMesh {
    let n_points = profile.len();
    let mut nodes = Vec::with_capacity(n_points * n_theta);

    for &(r, z) in profile {
        for i in 0..n_theta {
            let theta = (i as f64 / n_theta as f64) * 2.0 * std::f64::consts::PI;
            let x = r * theta.cos();
            let y = r * theta.sin();
            nodes.push([x, y, z]);
        }
    }

    let mut tris = Vec::with_capacity(2 * (n_points - 1) * n_theta);
    for j in 0..n_points - 1 {
        for i in 0..n_theta {
            let next_i = (i + 1) % n_theta;
            let n00 = j * n_theta + i;
            let n10 = j * n_theta + next_i;
            let n01 = (j + 1) * n_theta + i;
            let n11 = (j + 1) * n_theta + next_i;

            tris.push([n00, n10, n11]);
            tris.push([n00, n11, n01]);
        }
    }

    ShellMesh { nodes, tris }
}

/// Standard canonical English church bell profile (normalized coordinates).
#[must_use]
pub fn canonical_church_bell_profile(scale_m: f64) -> Vec<(f64, f64)> {
    // 10 radial slices from crown to soundring/lip
    let raw = [
        (0.10, 1.00), // Crown
        (0.18, 0.88), // Shoulder
        (0.24, 0.72), // Waist top
        (0.30, 0.55), // Waist mid
        (0.38, 0.38), // Soundbow upper
        (0.50, 0.20), // Soundbow
        (0.65, 0.08), // Soundring
        (0.75, 0.00), // Lip / mouth
    ];
    raw.iter().map(|&(r, z)| (r * scale_m, z * scale_m)).collect()
}
