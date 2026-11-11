//! The bending branch of the bore wall: the modal wall of a brass bell, coupled to the nonlinear bore cell by cell.
//!
//! The breathing ring of [`super::bore_wall`] carries the hoop stretch of the wall, which is stiff. Measured bells move
//! where the metal bends (Gorman, Rokni, Moore, Kausel and Chatziioannou, ISMA 2014: "it is much easier for the metal of
//! the bell to bend than for the radius of the metal tubing to increase"; the broad axial resonances that carry the
//! audible effect sit between 500 and 800 Hz). This module meshes the flaring part of the bore as a shell of
//! revolution, takes its modes from the existing MITC3+ shell owner (`fs_plate::shell`), and gives the bore a
//! [`WallDriver`]: every cell's pressure drives every mode, and every mode sweeps volume from every cell.
//!
//! The coupling (derivation note satwall/BELLWALL-DERIVATION-NOTE.md, section 2, generalised from one ring per cell):
//!
//! - mass-normalised modes `phi_j` (`omega_j`), the unit pressure on cell `i`'s band of facets gives the row
//!   `a_ji = int n . phi_j dA` [m^2 / sqrt(kg)] (`ShellPressurePort`, `n` the normal out of the gas);
//! - `q_j'' + c_j q_j' + omega_j^2 q_j = F_j`, `F_j = sum_i a_ji (p_mid,i - p_room)`, `c_j = eta omega_j`;
//! - the volume cell `i` sweeps is `dV_i = sum_j a_ji q_j`, and the bore rescales the cell for it (work `p_mid dV`);
//! - the force is the MIDPOINT one, `p_mid,i = p_i - k_i dV_i / 2` with `k_i = rho c^2 / V_i` of the room gas, which
//!   makes `(I + S diag(Gamma) / 2) F = b - S d / 2` one linear system of the mode count, `S = A K A^T`, re-factored
//!   whenever a cell's `k_i` changes (the gas moves it by a few percent over a loud note);
//! - each mode's step is exact for its held force, so `W_j' - W_j + D_j = F_j dq_j` and the account closes exactly:
//!   gas work = stored modal energy + dissipated work + work against the room pressure.
//!
//! The mesh is a polygon of flat facets (the card's `azimuths` points round, every quad cut by the same diagonal). A flat facet that carries a
//! radial breathing varying along the axis has an in-plane shear a circle does not, so the polygon is STIFFER than the circle,
//! by about 0.4 to 0.65 of `alpha^2 lambda^2` in the squared frequency (`alpha = pi / azimuths`, `lambda` the axial wave number
//! times the radius). Measured on a cylinder of R/t 100 and L/R 2 against Leissa's axisymmetric roots: the m = 1 hoop root is
//! 5.26 %, 3.00 %, 1.37 %, 0.80 % and 0.39 % too high at 12, 16, 24, 32 and 48 points, and the m = 2 root 4.00 %, 2.36 % and
//! 1.19 % too high at 24, 32 and 48: second order in the azimuth count, and growing with the axial wave number (at 24 points the
//! m = 3 hoop mode is 8.1 % high, the m = 5 mode 20 %). The azimuth count is therefore a convergence quantity of the card, like
//! the cut-off.
//!
//! Every quad is cut by the same diagonal, so the mesh keeps the full rotation symmetry of its azimuth count: the circumferential
//! harmonics stay separate, and every non-axisymmetric mode has a band participation at roundoff (at most 2.1e-11 of the largest
//! on the cylinder). The one-sided diagonals leave the mesh without a mirror, and the torsion (swirl) modes, which sweep no
//! volume in a shell of revolution, then show a weak participation (2.3e-2 of the largest at 24 points): they are dropped by
//! their swirl share. Cutting the diagonals alternately was tried and is worse: the azimuth-Nyquist harmonic then shares the
//! symmetry class of harmonic 0, mixes with the hoop modes (the ring-mean share of the m = 1 hoop mode fell to 0.36, of another to
//! 0.21, on 16 points and 20 rings) and has a spurious band participation of its own.
//!
//! The polygon is inscribed (its nodes lie on the gas radius): its swept area per unit radial displacement is
//! `(az / 2 pi) sin(2 pi / az)` of the circle's, 0.9745 at 16 points and 0.9936 at 32, and its midsurface sits at the gas
//! radius, not at that radius plus half the wall thickness. Both biases fall with the square of the azimuth count and are part
//! of what the card test's azimuth doubling moves. The clamp is at the first bore face at or after the throat (up to one cell
//! downstream of it); the receipt names the first covered cell. A radius step between two trailing cones is continuous in the
//! shell where the gas bore steps.
//!
//! A homogeneous brass wall has one loss factor `eta`, so its modal damping is diagonal, `c_j = eta omega_j` (the loss
//! stiffness is proportional to the stiffness). Modes that are not axisymmetric are dropped. Not carried: braces, the rim bead, the player's hand, radiation loss into the outside air, a non-uniform
//! thickness. The mount is one assumed case: the throat ring clamped, the rest free.

use super::bore_wall::{BELL_MAX_AZIMUTHS, BELL_MIN_AZIMUTHS, BellShellRecord as BellShellSpec};
use fs_duct::nonlinear::NonlinearDuctError;
use fs_duct::nonlinear::bore::{WallDriver, WallStepInput};
use fs_duct::{Duct, Segment};
use fs_math::det;
use fs_modal::SliceOptions;
use fs_phs::linear_exact::exponential;
use fs_plate::shell::mitc::{
    ConstraintFrame, MitcBudget, NodalDirectors, ShearShellSection, ShellConstraint, ShellMass, assemble_mitc3p,
};
use fs_plate::shell::pressure::{OrientedShellFacet, ShellFace, ShellPressurePort};
use fs_plate::shell::profile::{ProfileBudget, ProfileStation, revolve};
use fs_plate::shell::{ShellMesh, modes_shell};

/// The revolved mesh is turned this far about its axis [rad]: a mesh whose points sit at multiples of a quarter turn has pencil
/// entries that are exactly zero on one side of the diagonal and at roundoff on the other, and the sparse direct factor of the
/// modal solve refuses a structurally asymmetric pattern. The turn changes no physics (the bell is axisymmetric).
pub const AZIMUTH_OFFSET_RAD: f64 = 0.1;
/// Transverse-shear correction of the MITC3+ section.
const SHEAR_CORRECTION: f64 = 5.0 / 6.0;
/// A mode participates when its largest cell row is above this fraction of the largest row of any mode.
const PARTICIPATION: f64 = 1.0e-8;
/// The most modes a cut-off may admit: the slice returns all of them, and each is held as six numbers per mesh node.
const MAX_WINDOW_MODES: usize = 2048;
/// A mode is a torsion (swirl) mode when more than this share of its displacement power is tangential: in a shell of revolution of
/// an isotropic material the torsion decouples from the meridional motion and sweeps no volume (the hoop modes measured 0.000 to
/// 0.001, the torsion mode 1.000).
const SWIRL_SHARE: f64 = 0.5;

/// Brass constants of the bell shell: density, Young's modulus and Poisson ratio.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShellMaterial {
    /// Density [kg/m^3].
    pub density_kg_m3: f64,
    /// Young's modulus [Pa].
    pub young_pa: f64,
    /// Poisson ratio.
    pub poisson: f64,
}

/// The modal wall: a [`WallDriver`] for a bore whose cells it covers.
pub struct ModalBellWall {
    cells: usize,
    /// `rows[j][i]`: mode `j`'s participation in cell `i` [m^2 / sqrt(kg)].
    rows: Vec<Vec<f64>>,
    omega_rad_s: Vec<f64>,
    /// Per mode `[Phi_qq, Phi_qv, Phi_vq, Phi_vv, Gamma_q, Gamma_v]` of the exact held-force step.
    step: Vec<[f64; 6]>,
    state: Vec<[f64; 2]>,
    staged: Vec<[f64; 2]>,
    staged_lost: f64,
    dissipated: f64,
    /// A step is staged and not yet committed.
    pending: bool,
    /// Modes the cut-off admitted (every circumferential harmonic), before the axisymmetric ones were kept.
    modes_below_cutoff: usize,
    // Built from the gas of the step: `k_i`, `S` and the inverse of `I + S diag(Gamma) / 2` (row-major, m x m); rebuilt when `k_i` changes.
    factored: Option<Factored>,
}

struct Factored {
    stiffness_per_volume: Vec<f64>,
    s: Vec<f64>,
    inverse: Vec<f64>,
}

fn refuse(what: impl Into<String>) -> String {
    what.into()
}

impl ModalBellWall {
    /// The brass bell of `duct`: its trailing run of flaring cones, meshed at the bore's own faces and clamped at the
    /// throat. `faces` are the bore mesh's face positions [m] and `cell_volumes` its cell volumes (the bell cells are the
    /// ones from the first face at or after the throat).
    ///
    /// # Errors
    /// A bore without a flaring bell, a mesh or modal solve the shell owner refuses, or a spec outside its range.
    pub fn from_duct(
        duct: &Duct,
        faces: &[f64],
        spec: BellShellSpec,
        material: ShellMaterial,
        time_step_s: f64,
    ) -> Result<Self, String> {
        let cells = faces.len().checked_sub(1).ok_or_else(|| refuse("the bore has no cells"))?;
        let azimuths = spec.azimuths;
        if azimuths % 2 != 0 || !(BELL_MIN_AZIMUTHS..=BELL_MAX_AZIMUTHS).contains(&azimuths) {
            return Err(refuse(format!(
                "the bell mesh needs an even azimuth count from {BELL_MIN_AZIMUTHS} to {BELL_MAX_AZIMUTHS}, not {azimuths}"
            )));
        }
        // The trailing run of cones that widen.
        let mut start_index = duct.segments.len();
        while start_index > 0 {
            match &duct.segments[start_index - 1] {
                Segment::Cone {
                    inlet_radius,
                    outlet_radius,
                    ..
                } if outlet_radius > inlet_radius => start_index -= 1,
                _ => break,
            }
        }
        if start_index == duct.segments.len() {
            return Err(refuse(
                "bore_wall shell-brass needs a flaring bell: a trailing run of widening cones in the bore",
            ));
        }
        let mut origin = 0.0;
        for segment in &duct.segments[..start_index] {
            origin += segment.axial_length();
        }
        let bell = &duct.segments[start_index..];
        let radius_at = |x: f64| -> f64 {
            let mut at = origin;
            for segment in bell {
                if let Segment::Cone {
                    inlet_radius,
                    outlet_radius,
                    length,
                } = segment
                {
                    if x <= at + length + 1.0e-12 {
                        let t = ((x - at) / length).clamp(0.0, 1.0);
                        return (1.0 - t) * inlet_radius + t * outlet_radius;
                    }
                    at += length;
                }
            }
            0.0
        };
        let first_face = faces
            .iter()
            .position(|&x| x >= origin - 1.0e-12)
            .ok_or_else(|| refuse("the bore has no face at or after its bell"))?;
        let stations_count = faces.len() - first_face;
        if stations_count < 3 {
            return Err(refuse("the bell holds fewer than two bore cells"));
        }
        let stations: Vec<ProfileStation> = faces[first_face..]
            .iter()
            .map(|&x| ProfileStation {
                radius_m: radius_at(x),
                height_m: x - faces[first_face],
                thickness_m: spec.thickness_m,
            })
            .collect();
        let budget = ProfileBudget {
            max_nodes: stations_count * azimuths,
            max_triangles: 2 * (stations_count - 1) * azimuths,
            max_feature_evaluations: 0,
        };
        let shell = revolve(
            &stations,
            azimuths,
            material.young_pa,
            material.poisson,
            material.density_kg_m3,
            &[],
            &[],
            budget,
        )
        .map_err(|e| refuse(format!("bell mesh refused: {e}")))?;
        let (turn_cos, turn_sin) = (det::cos(AZIMUTH_OFFSET_RAD), det::sin(AZIMUTH_OFFSET_RAD));
        let turned: Vec<[f64; 3]> = shell
            .mesh
            .nodes
            .iter()
            .map(|n| [turn_cos * n[0] - turn_sin * n[1], turn_sin * n[0] + turn_cos * n[1], n[2]])
            .collect();
        let mesh = ShellMesh::new(turned, shell.mesh.tris.clone()).map_err(|e| refuse(format!("bell mesh refused: {e}")))?;
        // Directors: the unit normal of the winding, `(-dz r_hat + dr z_hat) / ds`, from the stations' finite differences.
        let mut vectors = Vec::with_capacity(shell.mesh.nodes.len());
        for j in 0..stations_count {
            let (lo, hi) = (j.saturating_sub(1), (j + 1).min(stations_count - 1));
            let (dr, dz) = (
                stations[hi].radius_m - stations[lo].radius_m,
                stations[hi].height_m - stations[lo].height_m,
            );
            let ds = det::sqrt(dr * dr + dz * dz);
            let (nr, nz) = (-dz / ds, dr / ds);
            for i in 0..azimuths {
                let theta = AZIMUTH_OFFSET_RAD + core::f64::consts::TAU * i as f64 / azimuths as f64;
                vectors.push([nr * det::cos(theta), nr * det::sin(theta), nz]);
            }
        }
        let directors = NodalDirectors::new(&mesh, vectors, 0.9)
            .map_err(|e| refuse(format!("bell directors refused: {e}")))?;
        let constraints: Vec<ShellConstraint> = (0..azimuths)
            .map(|node| ShellConstraint {
                node,
                mask: [true; 6],
                frame: ConstraintFrame::Global,
            })
            .collect();
        // Band `j` of the mesh is the ring between stations `j` and `j + 1`: triangles `2 j az .. 2 (j + 1) az`.
        let mut bands: Vec<Option<Vec<usize>>> = vec![None; cells];
        for (j, band) in (0..stations_count - 1).enumerate() {
            let cell = first_face + j;
            if cell < cells {
                bands[cell] = Some(((2 * band * azimuths)..(2 * (band + 1) * azimuths)).collect());
            }
        }
        Self::from_shell(&mesh, &directors, &constraints, &bands, spec, material, time_step_s)
    }

    /// Any shell mesh with explicit supports and, per bore cell, the facets of its band (`None`: a rigid cell). The
    /// normals of the facets' winding must point INTO the gas; the swept volume is then positive for an outward motion.
    ///
    /// # Errors
    /// A mesh, support or modal solve the shell owner refuses, no coupling mode below the cut-off, or a spec out of range.
    pub fn from_shell(
        mesh: &ShellMesh,
        directors: &NodalDirectors,
        constraints: &[ShellConstraint],
        bands: &[Option<Vec<usize>>],
        spec: BellShellSpec,
        material: ShellMaterial,
        time_step_s: f64,
    ) -> Result<Self, String> {
        if !(spec.thickness_m > 0.0
            && spec.thickness_m.is_finite()
            && spec.loss_factor >= 0.0
            && spec.loss_factor.is_finite()
            && spec.max_frequency_hz > 1.0
            && spec.max_frequency_hz.is_finite()
            && time_step_s > 0.0
            && time_step_s.is_finite())
        {
            return Err(refuse(
                "the bell shell needs a positive thickness, a nonnegative loss factor, a cut-off above 1 Hz and a positive step",
            ));
        }
        let g = material.young_pa / (2.0 * (1.0 + material.poisson));
        let section = ShearShellSection::from_orthotropic_solid(
            material.young_pa,
            material.young_pa,
            material.poisson,
            g,
            g,
            g,
            spec.thickness_m,
            material.density_kg_m3,
            0.0,
            SHEAR_CORRECTION,
        )
        .map_err(|e| refuse(format!("bell section refused: {e}")))?;
        let sections = vec![section; mesh.tris.len()];
        let assembly = MitcBudget {
            max_nodes: mesh.nodes.len(),
            max_triangles: mesh.tris.len(),
        };
        let model = assemble_mitc3p(
            mesh,
            directors,
            &sections,
            constraints,
            ShellMass::Consistent,
            assembly,
            1.0e-3,
            || false,
        )
        .map_err(|e| refuse(format!("bell shell assembly refused: {e}")))?;
        let high = (core::f64::consts::TAU * spec.max_frequency_hz).powi(2);
        let spectral = SliceOptions {
            mass_diagonal_equilibration: true,
            ..Default::default()
        };
        let report = modes_shell(model.pencil(), (1.0, high), &spectral)
            .map_err(|e| refuse(format!("bell modal solve refused: {e}")))?;
        let modes_below_cutoff = report.modes.len();
        if modes_below_cutoff > MAX_WINDOW_MODES {
            return Err(refuse(format!(
                "the cut-off admits {modes_below_cutoff} modes of the bell mesh (more than {MAX_WINDOW_MODES}): lower the cut-off"
            )));
        }
        let azimuths = spec.azimuths;
        if azimuths == 0 || mesh.nodes.len() % azimuths != 0 {
            return Err(refuse("the bell mesh is not whole rings of the card's azimuth count"));
        }
        // The 1-D bore drives and is driven by the axisymmetric modes only. The mesh keeps the rotation symmetry of its azimuth
        // count, so every non-axisymmetric mode has a band participation at roundoff and falls to the participation filter below;
        // the torsion modes, which sweep no volume, keep a weak participation from the polygon's one-sided diagonals and are
        // dropped here by their swirl share.
        let mut physical = Vec::new();
        let mut omega = Vec::new();
        for mode in &report.modes {
            let full = model.expand(&mode.phi).map_err(|e| refuse(format!("bell mode expansion refused: {e}")))?;
            if !is_swirl_mode(mesh, &full) {
                physical.push(full);
                omega.push(det::sqrt(mode.lambda));
            }
        }
        if physical.is_empty() {
            return Err(refuse("no mode of the bell shell other than a torsion mode lies below the cut-off"));
        }
        // Participation rows: one pressure port per cell band.
        let port_budget = ProfileBudget {
            max_nodes: mesh.nodes.len(),
            max_triangles: mesh.tris.len(),
            max_feature_evaluations: 0,
        };
        let cells = bands.len();
        let mut rows = vec![vec![0.0; cells]; physical.len()];
        for (i, band) in bands.iter().enumerate() {
            let Some(triangles) = band else { continue };
            let facets: Vec<OrientedShellFacet> = triangles
                .iter()
                .map(|&triangle| OrientedShellFacet {
                    triangle,
                    side: ShellFace::Negative,
                })
                .collect();
            let port = ShellPressurePort::new(mesh, facets, port_budget, || false)
                .map_err(|e| refuse(format!("bell pressure port refused: {e}")))?;
            let shapes = port
                .project_modes(&physical, usize::MAX, || false)
                .map_err(|e| refuse(format!("bell pressure projection refused: {e}")))?;
            for (j, value) in shapes.into_iter().enumerate() {
                rows[j][i] = value;
            }
        }
        let largest = rows.iter().flatten().fold(0.0_f64, |m, x| m.max(x.abs()));
        if !(largest > 0.0) {
            return Err(refuse("no mode of the bell shell sweeps volume from the bore (below the cut-off)"));
        }
        let keep: Vec<usize> = (0..rows.len())
            .filter(|&j| rows[j].iter().fold(0.0_f64, |m, x| m.max(x.abs())) > PARTICIPATION * largest)
            .collect();
        let rows: Vec<Vec<f64>> = keep.iter().map(|&j| rows[j].clone()).collect();
        let omega_rad_s: Vec<f64> = keep.iter().map(|&j| omega[j]).collect();
        // The exact held-force step of each mode: exp([[0, 1, 0], [-w^2, -c, 1], [0, 0, 0]] dt) on (q, v, F).
        let mut step = Vec::with_capacity(omega_rad_s.len());
        for &w in &omega_rad_s {
            let c = spec.loss_factor * w;
            let dt = time_step_s;
            let generator = [0.0, dt, 0.0, -dt * w * w, -dt * c, dt, 0.0, 0.0, 0.0];
            let e = exponential(&generator, 3).map_err(|_| refuse("a modal step exponential is not representable"))?;
            if e.iter().any(|x| !x.is_finite()) {
                return Err(refuse("a modal step exponential is not representable"));
            }
            step.push([e[0], e[1], e[3], e[4], e[2], e[5]]);
        }
        let modes = omega_rad_s.len();
        Ok(Self {
            cells,
            rows,
            omega_rad_s,
            step,
            state: vec![[0.0, 0.0]; modes],
            staged: vec![[0.0, 0.0]; modes],
            staged_lost: 0.0,
            dissipated: 0.0,
            pending: false,
            modes_below_cutoff,
            factored: None,
        })
    }

    /// Bore cells with a wall band (the bell's).
    #[must_use]
    pub fn covered_cells(&self) -> usize {
        (0..self.cells)
            .filter(|&i| self.rows.iter().any(|row| row[i] != 0.0))
            .count()
    }

    /// The first bore cell with a wall band (`None` when no cell is covered).
    #[must_use]
    pub fn first_covered_cell(&self) -> Option<usize> {
        (0..self.cells).find(|&i| self.rows.iter().any(|row| row[i] != 0.0))
    }

    /// Modes the cut-off admitted, of every circumferential harmonic; the coupling modes are the axisymmetric ones among them.
    #[must_use]
    pub const fn modes_below_cutoff(&self) -> usize {
        self.modes_below_cutoff
    }

    /// The coupling modes' frequencies [Hz], ascending.
    #[must_use]
    pub fn frequencies_hz(&self) -> Vec<f64> {
        self.omega_rad_s.iter().map(|w| w / core::f64::consts::TAU).collect()
    }

    /// The coupling modes' participation rows (`[mode][cell]`, m^2 / sqrt(kg)).
    #[must_use]
    pub fn participation(&self) -> &[Vec<f64>] {
        &self.rows
    }

    fn factor(&mut self, stiffness: Vec<f64>) -> Result<(), NonlinearDuctError> {
        let modes = self.rows.len();
        // S = A K A^T.
        let mut s = vec![0.0; modes * modes];
        for j in 0..modes {
            for l in 0..modes {
                s[j * modes + l] = (0..self.cells)
                    .map(|i| self.rows[j][i] * stiffness[i] * self.rows[l][i])
                    .sum();
            }
        }
        // M = I + S diag(Gamma_q) / 2, inverted by Gauss-Jordan with partial pivoting.
        let mut m = vec![0.0; modes * modes];
        for j in 0..modes {
            for l in 0..modes {
                m[j * modes + l] = f64::from(u8::from(j == l)) + 0.5 * s[j * modes + l] * self.step[l][4];
            }
        }
        let mut inverse = vec![0.0; modes * modes];
        for j in 0..modes {
            inverse[j * modes + j] = 1.0;
        }
        for col in 0..modes {
            let pivot = (col..modes)
                .max_by(|&a, &b| m[a * modes + col].abs().total_cmp(&m[b * modes + col].abs()))
                .unwrap_or(col);
            let scale = m[pivot * modes + col];
            if !(scale.abs() > 0.0 && scale.is_finite()) {
                return Err(NonlinearDuctError::UnsupportedGeometry(
                    "the modal wall's midpoint system is singular",
                ));
            }
            if pivot != col {
                for k in 0..modes {
                    m.swap(pivot * modes + k, col * modes + k);
                    inverse.swap(pivot * modes + k, col * modes + k);
                }
            }
            for k in 0..modes {
                m[col * modes + k] /= scale;
                inverse[col * modes + k] /= scale;
            }
            for row in 0..modes {
                if row != col {
                    let factor = m[row * modes + col];
                    if factor != 0.0 {
                        for k in 0..modes {
                            m[row * modes + k] -= factor * m[col * modes + k];
                            inverse[row * modes + k] -= factor * inverse[col * modes + k];
                        }
                    }
                }
            }
        }
        self.factored = Some(Factored {
            stiffness_per_volume: stiffness,
            s,
            inverse,
        });
        Ok(())
    }
}

impl ModalBellWall {
    /// A modal wall from explicit rows (`rows[mode][cell]`), frequencies and a loss factor: the tests' lumped cases.
    #[cfg(test)]
    fn from_parts(rows: Vec<Vec<f64>>, omega_rad_s: Vec<f64>, loss_factor: f64, time_step_s: f64) -> Self {
        let modes = rows.len();
        let cells = rows[0].len();
        let step = omega_rad_s
            .iter()
            .map(|&w| {
                let dt = time_step_s;
                let c = loss_factor * w;
                let generator = [0.0, dt, 0.0, -dt * w * w, -dt * c, dt, 0.0, 0.0, 0.0];
                let e = exponential(&generator, 3).unwrap();
                [e[0], e[1], e[3], e[4], e[2], e[5]]
            })
            .collect();
        Self {
            cells,
            rows,
            omega_rad_s,
            step,
            state: vec![[0.0, 0.0]; modes],
            staged: vec![[0.0, 0.0]; modes],
            staged_lost: 0.0,
            dissipated: 0.0,
            pending: false,
            modes_below_cutoff: modes,
            factored: None,
        }
    }
}

/// Whether a mode (six numbers per node, global axes) is a torsion mode about the z axis: more than [`SWIRL_SHARE`] of its
/// displacement power is tangential (perpendicular to the axis and to the radius of its node).
fn is_swirl_mode(mesh: &ShellMesh, mode: &[f64]) -> bool {
    let (mut tangential_power, mut total_power) = (0.0_f64, 0.0_f64);
    for node in 0..mesh.nodes.len().min(mode.len() / 6) {
        let [x, y, _] = mesh.nodes[node];
        let rho = det::sqrt(x * x + y * y);
        let (ux, uy, uz) = (mode[6 * node], mode[6 * node + 1], mode[6 * node + 2]);
        let tangential = if rho > 0.0 { (uy * x - ux * y) / rho } else { 0.0 };
        tangential_power += tangential * tangential;
        total_power += ux * ux + uy * uy + uz * uz;
    }
    tangential_power > SWIRL_SHARE * total_power
}

impl WallDriver for ModalBellWall {
    fn cells(&self) -> usize {
        self.cells
    }

    fn stage(
        &mut self,
        input: &WallStepInput<'_>,
        swept: &mut [f64],
        mean_pressure: &mut [f64],
    ) -> Result<(), NonlinearDuctError> {
        if input.pressure.len() != self.cells
            || input.gas_stiffness.len() != self.cells
            || input.volume.len() != self.cells
            || swept.len() != self.cells
            || mean_pressure.len() != self.cells
        {
            return Err(NonlinearDuctError::UnsupportedGeometry(
                "the modal wall covers a different number of cells than the bore",
            ));
        }
        self.pending = false;
        let stiffness: Vec<f64> = input.gas_stiffness.iter().zip(input.volume).map(|(k, v)| k / v).collect();
        if !stiffness.iter().all(|k| k.is_finite() && *k >= 0.0) {
            return Err(NonlinearDuctError::UnsupportedGeometry(
                "the modal wall needs a finite gas stiffness and a positive cell volume",
            ));
        }
        if self.factored.as_ref().is_none_or(|f| f.stiffness_per_volume != stiffness) {
            self.factor(stiffness)?;
        }
        let Some(factored) = self.factored.as_ref() else {
            return Err(NonlinearDuctError::UnsupportedGeometry("the modal wall is not factored"));
        };
        let modes = self.rows.len();
        // b_j = sum_i a_ji (p_i - p_room); d_j: the displacement with no force.
        let b: Vec<f64> = self
            .rows
            .iter()
            .map(|row| {
                row.iter()
                    .zip(input.pressure)
                    .map(|(a, p)| a * (p - input.room_pressure))
                    .sum()
            })
            .collect();
        let d: Vec<f64> = (0..modes)
            .map(|j| (self.step[j][0] - 1.0) * self.state[j][0] + self.step[j][1] * self.state[j][1])
            .collect();
        let rhs: Vec<f64> = (0..modes)
            .map(|j| b[j] - 0.5 * (0..modes).map(|l| factored.s[j * modes + l] * d[l]).sum::<f64>())
            .collect();
        let force: Vec<f64> = (0..modes)
            .map(|j| (0..modes).map(|l| factored.inverse[j * modes + l] * rhs[l]).sum())
            .collect();
        let mut lost = 0.0;
        let mut dq = vec![0.0; modes];
        for j in 0..modes {
            let [q, v] = self.state[j];
            let [pqq, pqv, pvq, pvv, gq, gv] = self.step[j];
            let q1 = pqq * q + pqv * v + gq * force[j];
            let v1 = pvq * q + pvv * v + gv * force[j];
            dq[j] = q1 - q;
            let w2 = self.omega_rad_s[j] * self.omega_rad_s[j];
            let stored = |q: f64, v: f64| 0.5 * (v * v + w2 * q * q);
            lost += force[j] * dq[j] - (stored(q1, v1) - stored(q, v));
            self.staged[j] = [q1, v1];
        }
        if !self.staged.iter().all(|[q, v]| q.is_finite() && v.is_finite()) {
            return Err(NonlinearDuctError::UnsupportedGeometry("the modal wall left the finite range"));
        }
        self.staged_lost = lost;
        self.pending = true;
        for i in 0..self.cells {
            let dv: f64 = (0..modes).map(|j| self.rows[j][i] * dq[j]).sum();
            swept[i] = dv;
            mean_pressure[i] = input.pressure[i] - 0.5 * factored.stiffness_per_volume[i] * dv;
        }
        Ok(())
    }

    fn commit(&mut self) {
        if !self.pending {
            return;
        }
        self.pending = false;
        core::mem::swap(&mut self.state, &mut self.staged);
        self.dissipated += self.staged_lost;
    }

    fn stored_energy_j(&self) -> f64 {
        self.state
            .iter()
            .zip(&self.omega_rad_s)
            .map(|([q, v], w)| 0.5 * (v * v + w * w * q * q))
            .sum()
    }

    fn dissipated_j(&self) -> f64 {
        self.dissipated
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TAU: f64 = core::f64::consts::TAU;

    // A full cylinder shell of revolution (radius r, length l), `az` points round, `nz` rings of facets, shear diaphragm
    // (radial and tangential displacement fixed, axial free, rotations free) at both ends; the winding normals point
    // INTO the cylinder. One band per ring of facets.
    fn cylinder(
        r: f64,
        l: f64,
        az: usize,
        nz: usize,
    ) -> (ShellMesh, NodalDirectors, Vec<ShellConstraint>, Vec<Option<Vec<usize>>>) {
        let mut nodes = Vec::new();
        let mut vectors = Vec::new();
        for k in 0..=nz {
            for i in 0..az {
                let theta = AZIMUTH_OFFSET_RAD + TAU * i as f64 / az as f64;
                nodes.push([r * det::cos(theta), r * det::sin(theta), l * k as f64 / nz as f64]);
                vectors.push([-det::cos(theta), -det::sin(theta), 0.0]);
            }
        }
        let mut tris = Vec::new();
        for j in 0..nz {
            for i in 0..az {
                let (a, b) = (j * az + i, j * az + (i + 1) % az);
                let (c, d) = ((j + 1) * az + i, (j + 1) * az + (i + 1) % az);
                tris.push([a, c, d]);
                tris.push([a, d, b]);
            }
        }
        let mesh = ShellMesh::new(nodes, tris).unwrap();
        // The admission tolerance: the dot product of the normals of two facets at a node, cos(2 pi / az), has to reach it.
        let directors = NodalDirectors::new(&mesh, vectors, 0.8).unwrap();
        let mut constraints = Vec::new();
        for k in [0, nz] {
            // Radial and tangential translation fixed = both in-plane translations (x and y) fixed, in any frame: the global frame.
            for i in 0..az {
                constraints.push(ShellConstraint {
                    node: k * az + i,
                    mask: [true, true, false, false, false, false],
                    frame: ConstraintFrame::Global,
                });
            }
        }
        let bands = (0..nz)
            .map(|j| Some(((2 * j * az)..(2 * (j + 1) * az)).collect()))
            .collect();
        (mesh, directors, constraints, bands)
    }

    // The steel-like cylinder of bl1: R 1 m, L 2 m, R/t 100.
    const CYLINDER: (f64, f64, f64) = (1.0, 2.0, 0.01);

    fn cylinder_material() -> ShellMaterial {
        ShellMaterial {
            density_kg_m3: 7800.0,
            young_pa: 200.0e9,
            poisson: 0.3,
        }
    }

    // The unit of Omega: the ring frequency sqrt(E / (rho (1 - nu^2))) / R [rad/s] of the cylinder.
    fn ring_frequency() -> f64 {
        let material = cylinder_material();
        det::sqrt(material.young_pa / (material.density_kg_m3 * (1.0 - material.poisson * material.poisson))) / CYLINDER.0
    }

    // One axisymmetric mode of the cylinder: Omega, its axial half-wave count (sign changes of the ring-mean radial displacement
    // plus one), the swirl (tangential) share of its displacement power and its largest band participation relative to the
    // largest of the listing.
    struct AxisymmetricMode {
        omega: f64,
        half_waves: usize,
        swirl: f64,
        participation: f64,
    }

    // The modes of the shear-diaphragm cylinder (`az` points round, `nz` rings of facets) with Omega in (`low`, `high`] whose
    // radial displacement is more than half circumferential harmonic 0, ascending in Omega, with the band participation of each.
    fn axisymmetric_modes(az: usize, nz: usize, low: f64, high: f64) -> Vec<AxisymmetricMode> {
        let (r, l, h) = CYLINDER;
        let material = cylinder_material();
        let scale = ring_frequency();
        let (mesh, directors, constraints, bands) = cylinder(r, l, az, nz);
        let g = material.young_pa / (2.0 * (1.0 + material.poisson));
        let section = ShearShellSection::from_orthotropic_solid(
            material.young_pa,
            material.young_pa,
            material.poisson,
            g,
            g,
            g,
            h,
            material.density_kg_m3,
            0.0,
            SHEAR_CORRECTION,
        )
        .unwrap();
        let sections = vec![section; mesh.tris.len()];
        let model = assemble_mitc3p(
            &mesh,
            &directors,
            &sections,
            &constraints,
            ShellMass::Consistent,
            MitcBudget {
                max_nodes: mesh.nodes.len(),
                max_triangles: mesh.tris.len(),
            },
            1.0e-3,
            || false,
        )
        .unwrap();
        let spectral = SliceOptions {
            mass_diagonal_equilibration: true,
            ..Default::default()
        };
        let report = modes_shell(model.pencil(), ((low * scale).powi(2), (high * scale).powi(2)), &spectral).unwrap();
        let mut kept: Vec<(f64, usize, f64, Vec<f64>)> = Vec::new();
        for mode in &report.modes {
            let full = model.expand(&mode.phi).unwrap();
            let (mut total, mut axisymmetric, mut swirl, mut power) = (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
            let mut means = vec![0.0_f64; nz + 1];
            for k in 0..=nz {
                let mut radial = Vec::with_capacity(az);
                for i in 0..az {
                    let node = k * az + i;
                    let theta = AZIMUTH_OFFSET_RAD + TAU * i as f64 / az as f64;
                    let (x, y, z) = (full[6 * node], full[6 * node + 1], full[6 * node + 2]);
                    radial.push(x * det::cos(theta) + y * det::sin(theta));
                    let tangential = -x * det::sin(theta) + y * det::cos(theta);
                    swirl += tangential * tangential;
                    power += x * x + y * y + z * z;
                }
                means[k] = radial.iter().sum::<f64>() / az as f64;
                axisymmetric += az as f64 * means[k] * means[k];
                total += radial.iter().map(|u| u * u).sum::<f64>();
            }
            if axisymmetric > 0.5 * total {
                let peak = means.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
                let (mut crossings, mut last) = (0, 0.0);
                for x in &means {
                    if x.abs() > 1.0e-3 * peak {
                        if last != 0.0 && (*x > 0.0) != (last > 0.0) {
                            crossings += 1;
                        }
                        last = *x;
                    }
                }
                kept.push((det::sqrt(mode.lambda) / scale, crossings + 1, swirl / power, full));
            }
        }
        let physical: Vec<Vec<f64>> = kept.iter().map(|k| k.3.clone()).collect();
        let port_budget = ProfileBudget {
            max_nodes: mesh.nodes.len(),
            max_triangles: mesh.tris.len(),
            max_feature_evaluations: 0,
        };
        let mut peak_row = vec![0.0_f64; physical.len()];
        for band in bands.iter().flatten() {
            let facets: Vec<OrientedShellFacet> = band
                .iter()
                .map(|&triangle| OrientedShellFacet {
                    triangle,
                    side: ShellFace::Negative,
                })
                .collect();
            let port = ShellPressurePort::new(&mesh, facets, port_budget, || false).unwrap();
            for (j, value) in port.project_modes(&physical, usize::MAX, || false).unwrap().into_iter().enumerate() {
                peak_row[j] = peak_row[j].max(value.abs());
            }
        }
        let largest = peak_row.iter().fold(0.0_f64, |m, x| m.max(*x));
        kept.iter()
            .zip(&peak_row)
            .map(|(k, row)| AxisymmetricMode {
                omega: k.0,
                half_waves: k.1,
                swirl: k.2,
                participation: row / largest,
            })
            .collect()
    }

    // Leissa (1973), NASA SP-288 ch. 2, Table 2.1, the Flugge-Lur'e-Byrne operator at circumferential wave n = 0 (the same
    // entries as fs-plate/tests/shell_vibration_locking.rs): the torsional component decouples, and the axial/radial pair
    // [[lam^2, nu lam - k lam^3], [. , 1 + k (lam^4 + 1)]] has the two roots Omega^2 returned (lowest first), Omega = omega R
    // sqrt(rho (1 - nu^2) / E), lam = m pi R / L, k = h^2 / (12 R^2).
    fn flugge_axisymmetric_omega_squared(m: usize, r: f64, l: f64, h: f64, nu: f64) -> [f64; 2] {
        let lam = m as f64 * core::f64::consts::PI * r / l;
        let k = h * h / (12.0 * r * r);
        let (a, b, c) = (lam * lam, nu * lam - k * lam * lam * lam, 1.0 + k * (lam.powi(4) + 1.0));
        let mean = 0.5 * (a + c);
        let spread = det::sqrt(0.25 * (a - c) * (a - c) + b * b);
        [mean - spread, mean + spread]
    }

    // bl1, re-posed (STAMP-bell-tests.txt, "bl1 re-posed", bounds written before this run). The first form of bl1 (every coupling
    // mode of the 24-point cylinder below Omega 1.02 equals a closed-form root within 4 %, in equal number) FAILED: the polygon of
    // flat facets is stiffer than the circle for a radial breathing that varies along the axis, by an error of order alpha^2
    // lambda^2 (alpha = pi / az, lambda = m pi R / L), which pushes the higher hoop roots out of the window. What the shell owes
    // is therefore a rate: on the cylinder (R/t 100, L/R 2, 40 rings, hoop modes of Omega between 0.90 and
    // 1.01) the hoop m = 1 root exceeds Leissa's from above at every azimuth count, the excess falls by a factor between 3 and 5
    // at each doubling of the azimuth count (12, 24, 48; second order is 4), the same for m = 2 between 24 and 48, the excess is at
    // most 0.5 % at 48 points, the hoop modes found at 48 points carry the half-wave counts 1, 2, 3 .. without gap or duplicate
    // (at least three), every hoop mode (participation above half of the largest) has a swirl share below 1e-2, and every torsion
    // mode (swirl share above one half) has a participation below 0.05 (the one-sided diagonals' artefact: 2.26e-2 at 24 points).
    #[test]
    fn bl1_the_hoop_roots_of_a_cylinder_converge_at_second_order_in_the_azimuth_count() {
        let (r, l, h) = CYLINDER;
        let nu = cylinder_material().poisson;
        let root = |m: usize| det::sqrt(flugge_axisymmetric_omega_squared(m, r, l, h, nu)[0]);
        let hoop = |az: usize| -> Vec<AxisymmetricMode> {
            axisymmetric_modes(az, 40, 0.90, 1.01)
                .into_iter()
                .filter(|mode| mode.participation > 0.5)
                .collect()
        };
        let (coarse, middle, fine) = (hoop(12), hoop(24), hoop(48));
        let excess = |modes: &[AxisymmetricMode], m: usize| -> f64 {
            let mode = modes.iter().find(|x| x.half_waves == m).unwrap_or_else(|| panic!("no hoop mode of {m} half-waves"));
            mode.omega / root(m) - 1.0
        };
        let (e1_coarse, e1_middle, e1_fine) = (excess(&coarse, 1), excess(&middle, 1), excess(&fine, 1));
        let (e2_middle, e2_fine) = (excess(&middle, 2), excess(&fine, 2));
        let list = |modes: &[AxisymmetricMode]| -> Vec<String> { modes.iter().map(|x| format!("{:.5}(m {})", x.omega, x.half_waves)).collect() };
        println!(
            "bl1: hoop modes 12 az {:?}, 24 az {:?}, 48 az {:?}; m = 1 excess {e1_coarse:.5} {e1_middle:.5} {e1_fine:.5} (ratios {:.3} {:.3}); m = 2 excess {e2_middle:.5} {e2_fine:.5} (ratio {:.3})",
            list(&coarse),
            list(&middle),
            list(&fine),
            e1_coarse / e1_middle,
            e1_middle / e1_fine,
            e2_middle / e2_fine
        );
        for (name, e) in [("12", e1_coarse), ("24", e1_middle), ("48", e1_fine), ("m2 24", e2_middle), ("m2 48", e2_fine)] {
            assert!(e > 0.0, "the polygon is stiffer than the circle: excess {e} at {name} points");
        }
        for (name, ratio) in [("12 to 24", e1_coarse / e1_middle), ("24 to 48", e1_middle / e1_fine), ("m = 2, 24 to 48", e2_middle / e2_fine)] {
            assert!((3.0..=5.0).contains(&ratio), "the excess falls by {ratio} ({name}), not by a factor of 3 to 5");
        }
        assert!(e1_fine <= 0.005, "the m = 1 excess at 48 points is {e1_fine}");
        let counts: Vec<usize> = fine.iter().map(|x| x.half_waves).collect();
        assert!(counts.len() >= 3 && counts == (1..=counts.len()).collect::<Vec<_>>(), "half-wave counts {counts:?}");
        let all = axisymmetric_modes(24, 40, 0.90, 1.01);
        assert!(all.iter().any(|mode| mode.swirl > SWIRL_SHARE), "the torsion mode is in the window");
        for mode in &all {
            if mode.participation > 0.5 {
                assert!(mode.swirl < 1.0e-2, "a hoop mode with swirl share {}", mode.swirl);
            }
            if mode.swirl > SWIRL_SHARE {
                assert!(mode.participation < 0.05, "a torsion mode with participation {}", mode.participation);
            }
        }
    }

    // The lumped cell, the ring written as one mode (a = A / sqrt(m), the numbers of the o1 case of the fs-duct ring tests:
    // r 15 mm, w 5 mm, E 1e8, rho 1200, t 0.3 mm, Q 20, a 100 Pa pulse): the scheme against the continuous coupled system
    // dp' = -k a v, q' = v, v' = a dp - omega^2 q - c v (k = rho c^2 / V). Bounds fixed before the run: the gap is below 2e-3
    // of the coupling's own scale dp eps / (1 + eps) (the model satwall/o1-model2.py gives gaps 1e2 times smaller); every
    // step's account closes within 1e-12 of its terms; the energy H = stored + dp^2 / (2 k) never exceeds its start by 1e-9.
    #[test]
    fn bl2_one_mode_one_cell_is_the_lumped_ring_and_its_account_closes() {
        let (r, w, t, young, rho_wall) = (15.0e-3, 5.0e-3, 0.3e-3, 1.0e8, 1200.0);
        let area = TAU * r * w;
        let volume = core::f64::consts::PI * r * r * w;
        let mass = rho_wall * t * area;
        let stiffness = young * t / (r * r) * area;
        let omega = det::sqrt(stiffness / mass);
        let a = area / det::sqrt(mass);
        let gas = 1.4 * 101_325.0; // rho c^2 [Pa]
        let k = gas / volume; // [Pa / m^3]
        let loss = 1.0 / 20.0; // Q = 20
        let dt = 1.0 / 512_000.0;
        let c = loss * omega;
        let mut wall = ModalBellWall::from_parts(vec![vec![a]], vec![omega], loss, dt);
        let dp0 = 100.0;
        let room = 101_325.0;
        let mut dp = dp0;
        let epsilon = k * a * a / (omega * omega);
        let omega_c = omega * (1.0 + epsilon).sqrt();
        let steps = (4.0 * TAU / omega_c / dt).round() as usize;
        let energy = |wall: &ModalBellWall, dp: f64| wall.stored_energy_j() + dp * dp / (2.0 * k);
        let h0 = energy(&wall, dp);
        let (mut swept, mut mean) = ([0.0], [0.0]);
        let (mut worst_account, mut highest, mut smallest_loss) = (0.0_f64, 0.0_f64, f64::INFINITY);
        let mut total_loss = 0.0;
        for _ in 0..steps {
            let stored_before = wall.stored_energy_j();
            let lost_before = wall.dissipated_j();
            let input = WallStepInput {
                pressure: &[room + dp],
                gas_stiffness: &[gas],
                volume: &[volume],
                room_pressure: room,
            };
            wall.stage(&input, &mut swept, &mut mean).unwrap();
            wall.commit();
            // The step's account: gas work p_mid dV = stored change + dissipated + work against the room pressure.
            let lost = wall.dissipated_j() - lost_before;
            let account = mean[0] * swept[0]
                - ((wall.stored_energy_j() - stored_before) + lost + room * swept[0]);
            worst_account = worst_account.max(account.abs() / ((mean[0] * swept[0]).abs() + (room * swept[0]).abs()));
            smallest_loss = smallest_loss.min(lost);
            total_loss += lost;
            // The gas: the volume swept lowers the pressure by k dV (first order, as the bore's rescale does).
            dp -= k * swept[0];
            highest = highest.max(energy(&wall, dp) / h0);
        }
        let exact = fs_phs::linear_exact::exponential(
            &[0.0, 0.0, -k * a, 0.0, 0.0, 1.0, a, -omega * omega, -c].map(|x| x * steps as f64 * dt),
            3,
        )
        .unwrap();
        let continuous = exact[0] * dp0;
        let gap = (dp - continuous).abs();
        let bound = 2.0e-3 * dp0 * epsilon / (1.0 + epsilon);
        println!(
            "bl2 eps {epsilon:.4}: scheme {dp:+.6} continuous {continuous:+.6} gap {gap:.3e} bound {bound:.3e}; worst account {worst_account:.2e}, highest H/H0 {highest:.12}, smallest step loss {smallest_loss:.3e} of a mean {:.3e}",
            total_loss / steps as f64
        );
        assert!(gap <= bound + 1e-9 * dp0, "{gap} {bound}");
        assert!(worst_account <= 1.0e-12, "{worst_account}");
        assert!(highest <= 1.0 + 1.0e-9, "{highest}");
        assert!(smallest_loss >= -1.0e-9 * total_loss / steps as f64, "{smallest_loss}");
        assert!(total_loss > 0.0);
    }

    // bl3, three cells of different volume and two modes with rows of mixed sign, against the exact exponential of the continuous
    // coupled system dp_i' = -k_i sum_j a_ji v_j, q_j' = v_j, v_j' = -omega_j^2 q_j - c_j v_j + sum_i a_ji dp_i. It pins what bl2
    // cannot: the matrix S, its inverse, the orientation of Gamma, the per-cell stiffness. Bounds fixed before the run
    // (STAMP-bell-tests.txt): every cell's pressure at the end is within 2e-3 of the coupling's own scale (the largest initial dp
    // times eps / (1 + eps), eps the largest k_i a_ji^2 / omega_j^2); every step's account closes within 1e-12 of its terms; the
    // coupling moves the pressures by more than 10 times the allowed gap (the test sees the coupling).
    #[test]
    fn bl3_three_cells_and_two_modes_follow_the_exact_coupled_system() {
        let gas = 1.4 * 101_325.0;
        let volumes = [2.0e-6, 4.0e-6, 8.0e-6];
        let k: Vec<f64> = volumes.iter().map(|v| gas / v).collect();
        let omegas = [TAU * 600.0, TAU * 1500.0];
        let a0 = 1.4e-3;
        let rows = vec![vec![a0, 0.7 * a0, 0.3 * a0], vec![0.5 * a0, -0.4 * a0, 0.9 * a0]];
        let loss = 0.02;
        let dt = 1.0 / 512_000.0;
        let room = 101_325.0;
        let dp0 = [200.0, -100.0, 50.0];
        let mut wall = ModalBellWall::from_parts(rows.clone(), omegas.to_vec(), loss, dt);
        let mut dp = dp0;
        let (mut swept, mut mean) = ([0.0; 3], [0.0; 3]);
        let steps = (4.0 * TAU / omegas[0] / dt).round() as usize;
        let mut worst_account = 0.0_f64;
        for _ in 0..steps {
            let (stored_before, lost_before) = (wall.stored_energy_j(), wall.dissipated_j());
            let pressure = [room + dp[0], room + dp[1], room + dp[2]];
            let input = WallStepInput {
                pressure: &pressure,
                gas_stiffness: &[gas; 3],
                volume: &volumes,
                room_pressure: room,
            };
            wall.stage(&input, &mut swept, &mut mean).unwrap();
            wall.commit();
            let lost = wall.dissipated_j() - lost_before;
            let work: f64 = (0..3).map(|i| mean[i] * swept[i]).sum();
            let against_room: f64 = (0..3).map(|i| room * swept[i]).sum();
            let account = work - ((wall.stored_energy_j() - stored_before) + lost + against_room);
            worst_account = worst_account.max(account.abs() / (work.abs() + against_room.abs()));
            for i in 0..3 {
                dp[i] -= k[i] * swept[i];
            }
        }
        // The exact system on the state (dp_0..2, q_0..1, v_0..1).
        let n = 7;
        let mut g = vec![0.0_f64; n * n];
        for i in 0..3 {
            for j in 0..2 {
                g[i * n + 5 + j] = -k[i] * rows[j][i];
                g[(5 + j) * n + i] = rows[j][i];
            }
        }
        for j in 0..2 {
            g[(3 + j) * n + 5 + j] = 1.0;
            g[(5 + j) * n + 3 + j] = -omegas[j] * omegas[j];
            g[(5 + j) * n + 5 + j] = -loss * omegas[j];
        }
        let t = steps as f64 * dt;
        let e = fs_phs::linear_exact::exponential(&g.iter().map(|x| x * t).collect::<Vec<_>>(), n).unwrap();
        let exact: Vec<f64> = (0..3).map(|i| (0..3).map(|c| e[i * n + c] * dp0[c]).sum()).collect();
        let epsilon = (0..2)
            .flat_map(|j| (0..3).map(move |i| (j, i)))
            .map(|(j, i)| k[i] * rows[j][i] * rows[j][i] / (omegas[j] * omegas[j]))
            .fold(0.0_f64, f64::max);
        let scale = dp0.iter().fold(0.0_f64, |m, x| m.max(x.abs())) * epsilon / (1.0 + epsilon);
        let gap = (0..3).map(|i| (dp[i] - exact[i]).abs()).fold(0.0_f64, f64::max);
        let moved = (0..3).map(|i| (exact[i] - dp0[i]).abs()).fold(0.0_f64, f64::max);
        println!(
            "bl3 eps {epsilon:.4}: scheme {dp:?} exact {exact:?}; gap {gap:.3e} of bound {:.3e}; the coupling moved the pressures by {moved:.3e}; worst account {worst_account:.2e}",
            2.0e-3 * scale
        );
        assert!(gap <= 2.0e-3 * scale, "{gap} against {}", 2.0e-3 * scale);
        assert!(moved > 10.0 * 2.0e-3 * scale, "the test does not see the coupling: {moved}");
        assert!(worst_account <= 1.0e-12, "{worst_account}");
    }

    // bl4: on the cylinder of 16 points and 20 rings the wall keeps exactly the hoop modes the oracle helper finds (participation above
    // half of the largest): the torsion mode (swirl share above one half, weak artefact participation) and every mode of another
    // harmonic are dropped; the weakest kept row is above half of the largest row; the cut-off admitted more modes than were kept;
    // the first covered cell is the first bore cell.
    #[test]
    fn bl4_the_wall_keeps_the_hoop_modes_and_drops_the_torsion_mode() {
        let (r, l, h) = CYLINDER;
        let material = cylinder_material();
        let scale = ring_frequency();
        let (mesh, directors, constraints, bands) = cylinder(r, l, 16, 20);
        let spec = BellShellSpec {
            thickness_m: h,
            loss_factor: 1.0e-4,
            max_frequency_hz: 1.15 * scale / TAU,
            azimuths: 16,
        };
        let wall = ModalBellWall::from_shell(&mesh, &directors, &constraints, &bands, spec, material, 1.0e-6).unwrap();
        let expected = axisymmetric_modes(16, 20, 0.05, 1.15)
            .iter()
            .filter(|mode| mode.participation > 0.5)
            .count();
        let rows = wall.participation();
        let largest = rows.iter().flatten().fold(0.0_f64, |m, x| m.max(x.abs()));
        let weakest = rows
            .iter()
            .map(|row| row.iter().fold(0.0_f64, |m, x| m.max(x.abs())))
            .fold(f64::INFINITY, f64::min);
        println!(
            "bl4: {} modes kept (the oracle helper finds {expected} hoop modes) of {} below the cut-off; weakest row {:.3} of the largest",
            rows.len(),
            wall.modes_below_cutoff(),
            weakest / largest
        );
        assert!(expected >= 3 && rows.len() == expected, "{} kept, {expected} expected", rows.len());
        assert!(weakest > 0.5 * largest, "a weak (mesh artefact) row was kept: {}", weakest / largest);
        assert!(wall.modes_below_cutoff() > rows.len());
        assert_eq!(wall.first_covered_cell(), Some(0));
    }
}
