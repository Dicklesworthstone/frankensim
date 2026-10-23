//! Hysteretic (constant loss factor) material loss of an assembled plate,
//! projected ONCE on the complete certified source modes.
//!
//! Wood's loss factor is roughly independent of frequency in the measured
//! audio bands, so a material share c with loss factor `eta_c` adds the loss
//! stiffness `eta_c K_c` to the elastic stiffness `K = sum_c K_c`. This module
//! splits the plate's elastic stiffness into named shares (orthotropic panel
//! bending components, and stiffener groups as differences of two assemblies),
//! forms each share's modal strain energy `E_c = Phi^T K_c Phi` on the same
//! mass-normalized source modes, and combines them with per-mode loss factors
//! into `H = sum_c sqrt(eta_c,a eta_c,b) E_c` (diagonal `sum_c eta_c,a E_c,aa`).
//!
//! The reduction consumes the diagonal: mode a's loss factor
//! `eta_a = H_aa / omega_a^2` and source damping `c_a = H_aa / omega_a` [1/s]
//! at its own frequency (the modal strain-energy approximation). The neglected
//! off-diagonal terms are measured, never dropped silently: per mode the
//! largest normalized coupling `|H_ab| / sqrt(H_aa H_bb)` and, for its nearest
//! lower and upper neighbours, `|H_ab| / |omega_a^2 - omega_b^2|`, the ratio that
//! decides whether the diagonal holds.
//!
//! This is an offline projection over the source slice (at most
//! [`MAX_SHARE_MODES`] modes). It creates no runtime dashpot channels; the
//! Kelvin-Voigt tensor image and its 64-connection limit are unchanged.
use super::{AcousticRealizeError, invalid, owner, poll};
use fs_exec::CancelGate;
use fs_plate::{PlateMesh, PlateModel};

/// The complete certified source slice a board reduction admits; not a
/// runtime connection count.
pub const MAX_SHARE_MODES: usize = 512;
/// Upper bound on the named stiffness shares of one projection.
pub const MAX_SHARES: usize = 64;

/// One named share of an assembled plate stiffness.
#[derive(Clone, Copy, Debug)]
pub enum StiffnessShare<'a> {
    /// Panel bending tensors [N m], one row-major 3x3 tensor per mesh triangle
    /// in mesh coordinates (DKT curvature order `xx, yy, 2xy`). Zero tensors
    /// contribute nothing; a tensor need not be definite.
    Panel(&'a [[f64; 9]]),
    /// The stiffness of `plus` minus that of `minus`: two assemblies of the
    /// same mesh, supports and DOF map that differ only in this share (for
    /// example a stiffener group present and absent).
    Difference { plus: &'a PlateModel, minus: &'a PlateModel },
}

/// One orthotropic panel section split into its four material-axis bending
/// components, each rotated to mesh coordinates exactly as fs-plate rotates the
/// whole section (`PlateSection::orthotropic_plane_stress_at_angle`): the L
/// (`D11`), R (`D22`), in-plane shear LR (`D66`) and Poisson (`D12`) parts.
/// Their sum is that section's bending tensor.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OrthotropicBendingParts {
    pub longitudinal: [f64; 9],
    pub radial: [f64; 9],
    pub shear: [f64; 9],
    pub poisson: [f64; 9],
}

impl OrthotropicBendingParts {
    /// The same constants fs-plate takes: `E_L, E_R, nu_LR, G_LR` [Pa],
    /// thickness [m] and material axis 1 at `angle_rad` from mesh x toward y.
    ///
    /// # Errors
    /// Nonpositive or nonfinite moduli or thickness, a nonfinite angle or a
    /// nonpositive `1 - nu_LR nu_RL` refuse.
    pub fn new(e_l: f64, e_r: f64, nu_lr: f64, g_lr: f64, thickness: f64, angle_rad: f64)
        -> Result<Self, AcousticRealizeError> {
        if [e_l, e_r, g_lr, thickness].iter().any(|v| !v.is_finite() || *v <= 0.0)
            || !nu_lr.is_finite() || !angle_rad.is_finite() {
            return Err(invalid("orthotropic bending parts need positive finite moduli and thickness"));
        }
        let denom = 1.0 - nu_lr * (nu_lr * e_r / e_l);
        if !denom.is_finite() || denom <= 0.0 {
            return Err(invalid("orthotropic bending parts need a positive 1 - nu_LR nu_RL"));
        }
        let h3 = thickness * thickness * thickness / 12.0;
        let (d11, d22, d12, d33) = (e_l * h3 / denom, e_r * h3 / denom, nu_lr * e_r * h3 / denom, g_lr * h3);
        // fs-plate keeps the aligned and the isotropic paths unrotated.
        let aligned = angle_rad == 0.0 || (e_l == e_r && g_lr == e_l / (2.0 * (1.0 + nu_lr)));
        let rotate = |d: [f64; 9]| -> [f64; 9] {
            if aligned { return d; }
            let (s, c) = angle_rad.sin_cos();
            let (cc, ss, cs) = (c * c, s * s, c * s);
            let t = [[cc, ss, cs], [ss, cc, -cs], [-2.0 * cs, 2.0 * cs, cc - ss]];
            let mut out = [0.0; 9];
            for i in 0..3 {
                for j in i..3 {
                    let mut value = 0.0;
                    for a in 0..3 { for b in 0..3 { value += t[a][i] * d[3 * a + b] * t[b][j]; } }
                    out[3 * i + j] = value;
                    out[3 * j + i] = value;
                }
            }
            out
        };
        Ok(Self {
            longitudinal: rotate([d11, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
            radial: rotate([0.0, 0.0, 0.0, 0.0, d22, 0.0, 0.0, 0.0, 0.0]),
            shear: rotate([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, d33]),
            poisson: rotate([0.0, d12, 0.0, d12, 0.0, 0.0, 0.0, 0.0, 0.0]),
        })
    }
}

/// Modal strain energy of every named share on the same source modes:
/// `E_c[a][b] = phi_a^T K_c phi_b` [1/s^2] for mass-normalized `phi`.
#[derive(Clone, Debug, PartialEq)]
pub struct ModalEnergyShares {
    names: Vec<String>,
    omega: Vec<f64>,
    energy: Vec<Vec<f64>>,
}

/// One neighbour's coupling `|H_ab| / |omega_a^2 - omega_b^2|`; infinite for a
/// repeated frequency with nonzero coupling.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NeighbourCoupling {
    pub mode: usize,
    pub ratio: f64,
}

/// The diagonal hysteretic image consumed by the reduction, with the measured
/// size of the off-diagonal terms it neglects.
#[derive(Clone, Debug, PartialEq)]
pub struct HystereticModalLoss {
    /// `eta_a = H_aa / omega_a^2`.
    pub loss_factor: Vec<f64>,
    /// `c_a = H_aa / omega_a` [1/s]: the diagonal viscous term at the mode's
    /// own frequency.
    pub source_damping: Vec<f64>,
    /// `max_{b != a} |H_ab| / sqrt(H_aa H_bb)` (zero for a single mode).
    pub max_coupling: Vec<f64>,
    /// The nearest lower and upper modes in frequency, when they exist.
    pub neighbours: Vec<[Option<NeighbourCoupling>; 2]>,
}

/// Project every share's modal strain energy once on the supplied source modes.
///
/// # Errors
/// An empty or oversized basis, misaligned or nonfinite modes or frequencies,
/// unnamed or repeated shares, a panel share with the wrong triangle count or
/// nonfinite entries, a difference whose assemblies do not share the model's
/// DOF map, a degenerate element and cancellation refuse.
pub fn project_energy_shares(mesh: &PlateMesh, model: &PlateModel, shares: &[(&str, StiffnessShare<'_>)],
    phi: &[Vec<f64>], omega: &[f64], gate: &CancelGate) -> Result<ModalEnergyShares, AcousticRealizeError> {
    let n = phi.len();
    let free = model.free;
    if n == 0 || n > MAX_SHARE_MODES || omega.len() != n {
        return Err(invalid("hysteretic projection needs 1..=512 modes with one frequency each"));
    }
    if phi.iter().any(|v| v.len() != free || v.iter().any(|x| !x.is_finite()))
        || omega.iter().any(|w| !w.is_finite() || *w <= 0.0) {
        return Err(invalid("hysteretic projection needs finite reduced-DOF modes and positive frequencies"));
    }
    if shares.is_empty() || shares.len() > MAX_SHARES {
        return Err(invalid("hysteretic projection needs 1..=64 named stiffness shares"));
    }
    let mut names: Vec<String> = Vec::with_capacity(shares.len());
    for (name, _) in shares {
        if name.trim().is_empty() || names.iter().any(|seen| seen.as_str() == *name) {
            return Err(invalid("every stiffness share needs a unique nonempty name"));
        }
        names.push((*name).to_string());
    }
    let mut energy = Vec::with_capacity(shares.len());
    for (_, share) in shares {
        poll(gate)?;
        // y_b = K_c phi_b on the reduced DOFs.
        let mut y = vec![vec![0.0; free]; n];
        match share {
            StiffnessShare::Panel(tensors) => {
                if tensors.len() != mesh.tris.len() || tensors.iter().flatten().any(|v| !v.is_finite()) {
                    return Err(invalid("a panel share needs one finite tensor per mesh triangle"));
                }
                for (e, tri) in mesh.tris.iter().enumerate() {
                    poll(gate)?;
                    let d = &tensors[e];
                    if d.iter().all(|v| *v == 0.0) { continue; }
                    let (c, _) = fs_plate::dkt_stiffness(&tri.map(|v| mesh.nodes[v].0),
                        &tri.map(|v| mesh.nodes[v].1), d, e).map_err(owner)?;
                    let dofs: [Option<usize>; 9] = std::array::from_fn(|i| model.dof_map[3 * tri[i / 3] + i % 3]);
                    for (b, out) in y.iter_mut().enumerate() {
                        let local = dofs.map(|v| v.map_or(0.0, |i| phi[b][i]));
                        for (a, row) in dofs.iter().enumerate() {
                            let Some(row) = *row else { continue };
                            out[row] += (0..9).map(|k| c[9 * a + k] * local[k]).sum::<f64>();
                        }
                    }
                }
            }
            StiffnessShare::Difference { plus, minus } => {
                if plus.free != free || minus.free != free || plus.dof_map != model.dof_map
                    || minus.dof_map != model.dof_map {
                    return Err(invalid("a stiffness difference needs two assemblies on the model's own DOF map"));
                }
                let mut less = vec![0.0; free];
                for (b, out) in y.iter_mut().enumerate() {
                    poll(gate)?;
                    plus.k.spmv(&phi[b], out);
                    minus.k.spmv(&phi[b], &mut less);
                    for (o, l) in out.iter_mut().zip(&less) { *o -= l; }
                }
            }
        }
        let mut matrix = vec![0.0; n * n];
        for a in 0..n {
            poll(gate)?;
            for b in a..n {
                let value: f64 = phi[a].iter().zip(&y[b]).map(|(p, q)| p * q).sum();
                if !value.is_finite() { return Err(invalid("hysteretic share energy overflows")); }
                matrix[n * a + b] = value;
                matrix[n * b + a] = value;
            }
        }
        energy.push(matrix);
    }
    Ok(ModalEnergyShares { names, omega: omega.to_vec(), energy })
}

/// [`project_energy_shares`] for an offline caller with no cancellation source
/// (board generation): the projection is never interrupted.
///
/// # Errors
/// Exactly those of [`project_energy_shares`] other than cancellation.
pub fn project_energy_shares_offline(mesh: &PlateMesh, model: &PlateModel, shares: &[(&str, StiffnessShare<'_>)],
    phi: &[Vec<f64>], omega: &[f64]) -> Result<ModalEnergyShares, AcousticRealizeError> {
    project_energy_shares(mesh, model, shares, phi, omega, &CancelGate::new_clock_free())
}

impl ModalEnergyShares {
    /// Share names in projection order.
    #[must_use] pub fn names(&self) -> &[String] { &self.names }
    /// Source angular frequencies [rad/s].
    #[must_use] pub fn omega(&self) -> &[f64] { &self.omega }
    /// Number of source modes.
    #[must_use] pub fn modes(&self) -> usize { self.omega.len() }
    /// `phi_a^T K_c phi_b` [1/s^2] of share `share`.
    #[must_use] pub fn energy(&self, share: usize, a: usize, b: usize) -> f64 {
        self.energy[share][self.omega.len() * a + b]
    }
    /// Share `share`'s fraction of mode a's strain energy, `E_c,aa / omega_a^2`.
    #[must_use] pub fn fraction(&self, share: usize, a: usize) -> f64 {
        self.energy(share, a, a) / (self.omega[a] * self.omega[a])
    }
    /// The sum of every share's fraction of mode a: one when the shares
    /// partition the elastic stiffness that produced the modes.
    #[must_use] pub fn fraction_sum(&self, a: usize) -> f64 {
        (0..self.names.len()).map(|c| self.fraction(c, a)).sum()
    }

    /// Combine the shares with per-mode loss factors: `eta[c][a]` is share c's
    /// loss factor at mode a (zero for a share whose loss is omitted).
    ///
    /// # Errors
    /// A loss table that does not cover every share and mode, negative or
    /// nonfinite loss factors, and a non-dissipative or nonfinite diagonal refuse.
    pub fn hysteretic(&self, eta: &[Vec<f64>]) -> Result<HystereticModalLoss, AcousticRealizeError> {
        let n = self.omega.len();
        if eta.len() != self.names.len() || eta.iter().any(|row| row.len() != n
            || row.iter().any(|v| !v.is_finite() || *v < 0.0)) {
            return Err(invalid("hysteretic loss needs one nonnegative finite loss factor per share and mode"));
        }
        let h = |a: usize, b: usize| -> f64 {
            (0..self.names.len()).map(|c| (eta[c][a] * eta[c][b]).sqrt() * self.energy(c, a, b)).sum()
        };
        let diagonal: Vec<f64> = (0..n).map(|a| h(a, a)).collect();
        if diagonal.iter().any(|v| !v.is_finite() || *v < 0.0) {
            return Err(invalid("hysteretic material loss must be finite and dissipative on every mode"));
        }
        let mut max_coupling = vec![0.0_f64; n];
        for a in 0..n {
            for b in 0..n {
                if a == b { continue; }
                let off = h(a, b).abs();
                let scale = (diagonal[a] * diagonal[b]).sqrt();
                let ratio = if scale > 0.0 { off / scale } else if off > 0.0 { f64::INFINITY } else { 0.0 };
                max_coupling[a] = max_coupling[a].max(ratio);
            }
        }
        let neighbours = (0..n).map(|a| {
            let coupling = |b: usize| {
                let gap = (self.omega[a] * self.omega[a] - self.omega[b] * self.omega[b]).abs();
                let off = h(a, b).abs();
                NeighbourCoupling { mode: b,
                    ratio: if gap > 0.0 { off / gap } else if off > 0.0 { f64::INFINITY } else { 0.0 } }
            };
            let lower = (0..n).filter(|&b| b != a && self.omega[b] <= self.omega[a])
                .max_by(|&x, &y| self.omega[x].total_cmp(&self.omega[y]).then(x.cmp(&y)));
            let upper = (0..n).filter(|&b| b != a && self.omega[b] > self.omega[a])
                .min_by(|&x, &y| self.omega[x].total_cmp(&self.omega[y]).then(x.cmp(&y)));
            [lower.map(coupling), upper.map(coupling)]
        }).collect();
        Ok(HystereticModalLoss {
            loss_factor: (0..n).map(|a| diagonal[a] / (self.omega[a] * self.omega[a])).collect(),
            source_damping: (0..n).map(|a| diagonal[a] / self.omega[a]).collect(),
            max_coupling, neighbours,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_plate::{AssemblyOptions, EdgeSupport, PlateSection, SliceOptions, Stiffener};

    const E_L: f64 = 11.9e9;
    const E_R: f64 = 0.93e9;
    const NU: f64 = 0.37;
    const G_LR: f64 = 0.76e9;
    const H: f64 = 0.008;
    const RHO: f64 = 450.0;
    const ANGLE: f64 = 0.61;

    /// A clamped orthotropic panel with one eccentric stiffener, its source
    /// modes, and the same assembly without the stiffener and with its E doubled.
    struct Board { mesh: PlateMesh, model: PlateModel, bare: PlateModel, stiff: PlateModel,
        phi: Vec<Vec<f64>>, omega: Vec<f64> }
    fn board() -> Board {
        let (nx, ny) = (8, 6);
        let mesh = PlateMesh::rectangle(0.6, 0.45, nx, ny);
        let boundary = PlateMesh::rectangle_boundary(nx, ny);
        let section = PlateSection::orthotropic_plane_stress_at_angle(E_L, E_R, NU, G_LR, H, RHO, ANGLE).unwrap();
        let rib = |e: f64| Stiffener { nodes: (0..=nx).map(|i| 3 * (nx + 1) + i).collect(), e, g: 0.6e9,
            area: 0.025 * 0.02, inertia: 0.025 * 0.02_f64.powi(3) / 12.0, torsion: 3.0e-8,
            eccentricity: -0.014, density: 400.0 };
        let options = AssemblyOptions { support: EdgeSupport::Clamped, pretension: 0.0 };
        let assemble = |ribs: &[Stiffener]| fs_plate::assemble(&mesh, &section, &boundary, ribs, &options).unwrap();
        let model = assemble(&[rib(10e9)]);
        let report = fs_plate::modes(&model, (0.0, (std::f64::consts::TAU * 2500.0).powi(2)),
            &SliceOptions::default()).unwrap();
        assert!(report.modes.len() >= 6, "{} modes", report.modes.len());
        let phi = report.modes.iter().map(|pair| pair.phi.clone()).collect();
        let omega = report.modes.iter().map(|pair| pair.lambda.sqrt()).collect();
        Board { bare: assemble(&[]), stiff: assemble(&[rib(20e9)]), mesh, model, phi, omega }
    }
    fn parts(mesh: &PlateMesh) -> Vec<OrthotropicBendingParts> {
        vec![OrthotropicBendingParts::new(E_L, E_R, NU, G_LR, H, ANGLE).unwrap(); mesh.tris.len()]
    }
    /// Panel L, R, LR, Poisson, rib bending and rib torsion: a partition of K.
    fn shares(b: &Board) -> ModalEnergyShares {
        let p = parts(&b.mesh);
        let pick = |f: fn(&OrthotropicBendingParts) -> [f64; 9]| p.iter().map(f).collect::<Vec<_>>();
        let (l, r, lr, nu) = (pick(|x| x.longitudinal), pick(|x| x.radial), pick(|x| x.shear), pick(|x| x.poisson));
        let gate = CancelGate::new_clock_free();
        // Bending is linear in E: (E doubled) - (E) is the rib's bending share;
        // the rest of (with) - (without) is its torsion.
        project_energy_shares(&b.mesh, &b.model, &[("panel-l", StiffnessShare::Panel(&l)),
            ("panel-r", StiffnessShare::Panel(&r)), ("panel-lr", StiffnessShare::Panel(&lr)),
            ("panel-12", StiffnessShare::Panel(&nu)),
            ("rib", StiffnessShare::Difference { plus: &b.model, minus: &b.bare }),
            ("rib-bending", StiffnessShare::Difference { plus: &b.stiff, minus: &b.model })],
            &b.phi, &b.omega, &gate).unwrap()
    }

    // G0: the four material-axis parts sum to fs-plate's own rotated section.
    #[test]
    fn g0_bending_parts_sum_to_the_fs_plate_section() {
        for angle in [0.0, ANGLE, -1.2] {
            let p = OrthotropicBendingParts::new(E_L, E_R, NU, G_LR, H, angle).unwrap();
            let d = PlateSection::orthotropic_plane_stress_at_angle(E_L, E_R, NU, G_LR, H, RHO, angle).unwrap().d;
            let scale = d.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
            for i in 0..9 {
                let sum = p.longitudinal[i] + p.radial[i] + p.shear[i] + p.poisson[i];
                assert!((sum - d[i]).abs() <= 1e-13 * scale, "angle {angle} entry {i}: {sum:e} vs {:e}", d[i]);
            }
        }
        assert!(OrthotropicBendingParts::new(E_L, E_R, 4.0, G_LR, H, 0.0).is_err());
        assert!(OrthotropicBendingParts::new(E_L, -1.0, NU, G_LR, H, 0.0).is_err());
        assert!(OrthotropicBendingParts::new(E_L, E_R, NU, G_LR, H, f64::NAN).is_err());
    }

    // G0: the shares partition the stiffness that produced the modes, so every
    // mode's fractions sum to one; a uniform loss factor gives every mode that
    // loss factor; and a source damping is its loss factor times omega.
    #[test]
    fn g0_shares_partition_each_mode_and_a_uniform_loss_factor_is_returned() {
        let b = board();
        let e = shares(&b);
        let n = e.modes();
        // rib = rib-bending + rib-torsion; the partition uses rib-torsion = rib - rib-bending.
        for a in 0..n {
            let sum = e.fraction(0, a) + e.fraction(1, a) + e.fraction(2, a) + e.fraction(3, a) + e.fraction(4, a);
            assert!((sum - 1.0).abs() < 1e-9, "mode {a}: {sum}");
            assert!(e.fraction(5, a) >= 0.0 && e.fraction(5, a) <= e.fraction(4, a) + 1e-12, "mode {a}");
        }
        let eta = 0.013;
        let uniform: Vec<Vec<f64>> = (0..6).map(|c| vec![if c == 5 { 0.0 } else { eta }; n]).collect();
        let loss = e.hysteretic(&uniform).unwrap();
        for a in 0..n {
            assert!((loss.loss_factor[a] / eta - 1.0).abs() < 1e-9, "mode {a}: {}", loss.loss_factor[a]);
            assert!((loss.source_damping[a] / (loss.loss_factor[a] * b.omega[a]) - 1.0).abs() < 1e-12);
            // Uniform loss is proportional to K, which the modes diagonalize.
            assert!(loss.max_coupling[a] < 1e-6, "mode {a}: {}", loss.max_coupling[a]);
        }
    }

    // G1: changing one share's loss factor moves each mode's loss factor by
    // exactly that change times the mode's energy fraction in the share.
    #[test]
    fn g1_one_share_change_moves_each_mode_by_its_energy_fraction() {
        let b = board();
        let e = shares(&b);
        let n = e.modes();
        let base: Vec<Vec<f64>> = vec![vec![0.0124; n], vec![0.0188; n], vec![0.0155; n], vec![0.0; n],
            vec![0.0; n], vec![0.0; n]];
        let mut swapped = base.clone();
        swapped[0] = vec![0.0070; n];
        let (before, after) = (e.hysteretic(&base).unwrap(), e.hysteretic(&swapped).unwrap());
        for a in 0..n {
            let predicted = (0.0070 - 0.0124) * e.fraction(0, a);
            assert!(((after.loss_factor[a] - before.loss_factor[a]) - predicted).abs() < 1e-15 + 1e-12 * predicted.abs(),
                "mode {a}");
            assert!(after.loss_factor[a] < before.loss_factor[a], "mode {a}");
        }
        // Nonuniform loss is not proportional to K: coupling appears and is measured.
        assert!(before.max_coupling.iter().any(|c| *c > 1e-6));
        for (a, pair) in before.neighbours.iter().enumerate() {
            assert!(pair.iter().flatten().all(|c| c.mode != a && c.ratio.is_finite() && c.ratio >= 0.0));
        }
        assert!(before.neighbours[0][0].is_none() && before.neighbours[n - 1][1].is_none());
    }

    // G0: malformed inputs refuse before any projection.
    #[test]
    fn g0_malformed_projection_inputs_refuse() {
        let b = board();
        let gate = CancelGate::new_clock_free();
        let p = parts(&b.mesh);
        let l: Vec<[f64; 9]> = p.iter().map(|x| x.longitudinal).collect();
        let panel = [("panel-l", StiffnessShare::Panel(&l))];
        assert!(project_energy_shares(&b.mesh, &b.model, &panel, &[], &[], &gate).is_err());
        assert!(project_energy_shares(&b.mesh, &b.model, &panel, &b.phi, &b.omega[1..], &gate).is_err());
        assert!(project_energy_shares(&b.mesh, &b.model, &[], &b.phi, &b.omega, &gate).is_err());
        assert!(project_energy_shares(&b.mesh, &b.model, &[("a", StiffnessShare::Panel(&l)),
            ("a", StiffnessShare::Panel(&l))], &b.phi, &b.omega, &gate).is_err());
        assert!(project_energy_shares(&b.mesh, &b.model, &[("short", StiffnessShare::Panel(&l[1..]))],
            &b.phi, &b.omega, &gate).is_err());
        let mut bad = b.omega.clone();
        bad[0] = -1.0;
        assert!(project_energy_shares(&b.mesh, &b.model, &panel, &b.phi, &bad, &gate).is_err());
        let e = project_energy_shares(&b.mesh, &b.model, &panel, &b.phi, &b.omega, &gate).unwrap();
        let n = e.modes();
        assert!(e.hysteretic(&[vec![0.01; n - 1]]).is_err());
        assert!(e.hysteretic(&[vec![-0.01; n]]).is_err());
        assert!(e.hysteretic(&[vec![f64::NAN; n]]).is_err());
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(project_energy_shares(&b.mesh, &b.model, &panel, &b.phi, &b.omega, &cancelled).is_err());
    }
}
