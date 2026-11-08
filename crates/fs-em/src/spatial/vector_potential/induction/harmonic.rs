//! Three-dimensional time-harmonic induction on the original Whitney A–phi space.
//!
//! Peak phasors use `exp(-i omega t)`, as in `fs_phs::harmonic`. With the
//! integrated electric potential `q = phi / (-i omega)`, the original coupled
//! conductivity form M gives `(K - V - i omega M) [A, q] = [j_impressed, 0]`.
//! V is the same sigma*(v cross curl) operator used by transient induction.
//! Floating charge, relative tree gauge, harmonic circulation and anisotropic
//! material coefficients are retained; no dense whole-field reduction occurs.
//!
//! This solves the AC component in fixed Eulerian geometry and linear media,
//! with prescribed steady velocity and homogeneous AC tangential-A traces.
//! Any original DC field/heat is separate. Permanent-magnet operating envelopes,
//! saturation, hysteresis, terminal circuits, transported material occupancy,
//! displacement current and full-wave radiation need their own owners.
//! A skin-depth admission limit screens resolution; it is not an error bound.
//!
//! Formulation reference: Ansys Maxwell, "AC Magnetic (Eddy Current) A-Phi
//! Solver" (the MQS reduction omits displacement and wave-propagation terms):
//! <https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v261/en/Subsystems/Maxwell/Content/Maxwell/ACMagnetic_APhiSolver.htm>.

/// Original electrical networks coupled to the resolved 3D harmonic field.
pub mod circuit;
/// Consistent mean and twice-frequency mechanical loads on original bodies.
pub mod loads;
/// Closed source-bound winding terminals with current and voltage drive.
pub mod winding;

use super::super::super::{cross, norm, sub};
use super::*;
use fs_math::c64::C64;

/// Work limits and physical residual allowances for one AC frequency.
#[derive(Clone, Copy, Debug)]
pub struct SpatialHarmonicInductionControls {
    pub relative_tolerance: f64,
    pub absolute_residual_a: f64,
    pub restart: usize,
    pub maximum_cycles: usize,
    /// Includes BOTH real blocks, restart vectors, Hessenberg data and history.
    pub maximum_krylov_scalars: usize,
    /// Additional field/recovery/motion payload, excluding retained preparation.
    pub maximum_field_scalars: usize,
    pub current_tolerance_a: f64,
    pub flux_tolerance_wb: f64,
    /// Applies to mean real and reactive power closure [W / var].
    pub power_tolerance_w: f64,
    pub maximum_cell_magnetic_peclet: f64,
    /// max(diameter * sqrt(omega * mu_bound * sigma_bound / 2)).
    /// Principal-coefficient bounds also screen orthotropic materials.
    pub maximum_cell_skin_ratio: f64,
}

/// Original local-vertex fields; phasors are peaks, forces and powers are means.
#[derive(Clone, Debug)]
pub struct SpatialHarmonicInductionCell {
    pub cell: usize,
    pub component: usize,
    pub volume_m3: f64,
    pub potential_v: [C64; 4],
    pub laboratory_electric_field_v_m: [[C64; 3]; 4],
    pub conductor_electric_field_v_m: [[C64; 3]; 4],
    pub current_density_a_m2: [[C64; 3]; 4],
    pub mean_lorentz_force_density_n_m3: [[f64; 3]; 4],
    pub mean_ohmic_power_w: f64,
    /// Positive means field-to-prescribed-motion work.
    pub mean_mechanical_power_w: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct SpatialHarmonicInductionBalance {
    /// AC perturbation storage, V*(Br.nu.Br + Bi.nu.Bi)/4.
    pub mean_magnetic_energy_j: f64,
    pub mean_ohmic_power_w: f64,
    pub mean_mechanical_power_w: f64,
    /// Im(integral((conj(J) cross B).v)/2), needed by moving AC balance.
    pub reactive_mechanical_power_var: f64,
    pub mean_supplied_power_w: f64,
    /// Im(sum(conj(j_impressed) * (-i omega A)) / 2).
    pub reactive_supplied_power_var: f64,
    pub power_residual_w: f64,
    pub reactive_power_residual_var: f64,
}

/// An accepted AC response retaining its original geometry and material sources.
/// Solving does not mutate the transient history or advance its clock.
pub struct SpatialHarmonicInductionResponse<
    's,
    'm,
    C: SpatialConductorMaterial = ResolvedConductorState,
> {
    source: &'s SpatialInduction<'m, C>,
    omega_rad_s: f64,
    // Original reduced A/q coordinates, real block followed by imaginary.
    // Retained for exact linear superposition followed by original residual
    // and physical recovery checks; never a detached field authority.
    solution: Vec<f64>,
    potential_wb: Vec<C64>,
    flux_density_t: Vec<[C64; 3]>,
    cells: Vec<SpatialHarmonicInductionCell>,
    balance: SpatialHarmonicInductionBalance,
    report: SolveReport,
    maximum_ampere_residual_a: f64,
    maximum_charge_residual_a: f64,
    maximum_cell_magnetic_peclet: f64,
    maximum_cell_skin_ratio: f64,
}
impl<C: SpatialConductorMaterial> SpatialHarmonicInductionResponse<'_, '_, C> {
    pub fn problem(&self) -> &SpatialVectorMagnetostatics {
        self.source.problem()
    }
    pub fn conductors(&self) -> &[Option<C>] {
        self.source.conductors()
    }
    pub fn omega_rad_s(&self) -> f64 {
        self.omega_rad_s
    }
    pub fn potential_wb(&self) -> &[C64] {
        &self.potential_wb
    }
    pub fn flux_density_t(&self) -> &[[C64; 3]] {
        &self.flux_density_t
    }
    pub fn cells(&self) -> &[SpatialHarmonicInductionCell] {
        &self.cells
    }
    pub fn balance(&self) -> &SpatialHarmonicInductionBalance {
        &self.balance
    }
    pub fn report(&self) -> &SolveReport {
        &self.report
    }
    pub fn maximum_ampere_residual_a(&self) -> f64 {
        self.maximum_ampere_residual_a
    }
    pub fn maximum_charge_residual_a(&self) -> f64 {
        self.maximum_charge_residual_a
    }
    pub fn maximum_cell_magnetic_peclet(&self) -> f64 {
        self.maximum_cell_magnetic_peclet
    }
    pub fn maximum_cell_skin_ratio(&self) -> f64 {
        self.maximum_cell_skin_ratio
    }
}

impl<'m, C: SpatialConductorMaterial> SpatialInduction<'m, C> {
    /// Resolve one sinusoidal impressed-current field, including induced
    /// currents and their charge potentials. Frequency is angular [rad/s],
    /// sources are peak phasors [A/m²], velocity is held [m/s].
    pub fn solve_harmonic(
        &self,
        impressed_current_a_m2: &[[C64; 3]],
        velocity_m_s: &[[f64; 3]],
        omega_rad_s: f64,
        controls: SpatialHarmonicInductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHarmonicInductionResponse<'_, 'm, C>> {
        self.solve_harmonic_seeded(
            impressed_current_a_m2,
            velocity_m_s,
            omega_rad_s,
            controls,
            None,
            cancelled,
        )
    }

    // Cached unit winding solutions are only initial guesses. The same
    // operator, original-equation gates and field/power recovery below remain
    // authoritative at the actual combined source amplitude.
    #[allow(clippy::too_many_arguments)]
    fn solve_harmonic_seeded(
        &self,
        impressed_current_a_m2: &[[C64; 3]],
        velocity_m_s: &[[f64; 3]],
        omega_rad_s: f64,
        controls: SpatialHarmonicInductionControls,
        seed: Option<&[f64]>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHarmonicInductionResponse<'_, 'm, C>> {
        poll(cancelled)?;
        let p = &self.prepared;
        let problem = self.problem();
        let g = geometry(problem);
        let mesh = g.mesh();
        let n = p.mass.nrows();
        let restart = admit(
            self,
            impressed_current_a_m2,
            velocity_m_s,
            omega_rad_s,
            controls,
        )?;
        if seed.is_some_and(|values| {
            values.len() != 2 * n || values.iter().any(|value| !value.is_finite())
        }) {
            return Err(invalid(
                "finite complete original harmonic coordinates required",
            ));
        }
        let omega = omega_rad_s;
        let re: Vec<_> = impressed_current_a_m2
            .iter()
            .map(|j| j.map(|v| v.re))
            .collect();
        let im: Vec<_> = impressed_current_a_m2
            .iter()
            .map(|j| j.map(|v| v.im))
            .collect();
        let load_re = p.current_load(g, &re, controls.current_tolerance_a, cancelled)?;
        let load_im = p.current_load(g, &im, controls.current_tolerance_a, cancelled)?;
        let mut rhs = vec![0.0; 2 * n];
        for (edge, index) in p.edge_unknown.iter().enumerate() {
            if let Some(i) = *index {
                rhs[i] = load_re[edge];
                rhs[n + i] = load_im[edge];
            }
        }
        let mut maximum_cell_skin_ratio = 0.0_f64;
        for cell in &p.cells {
            poll(cancelled)?;
            let vertices = mesh.tets[cell.cell].map(|v| problem.positions_m()[v as usize]);
            let mut diameter = 0.0_f64;
            for a in 0..4 {
                for b in a + 1..4 {
                    diameter = diameter.max(norm(sub(vertices[a], vertices[b])));
                }
            }
            // Factorized square roots avoid overflowing omega*mu*sigma when
            // the final dimensionless ratio is representable.
            let ratio = finite(
                diameter
                    * omega.sqrt()
                    * problem.permeability_bound[cell.cell].sqrt()
                    * (0.5 * cell.sigma.maximum_conductivity_s_m()).sqrt(),
            )?;
            maximum_cell_skin_ratio = maximum_cell_skin_ratio.max(ratio);
            if ratio > controls.maximum_cell_skin_ratio {
                return Err(invalid(
                    "harmonic induction skin-depth resolution allowance exceeded",
                ));
            }
        }
        // Stationary frequency sweeps reuse the prepared K/M directly: no
        // cell-quadrature assembly of a zero motion matrix or zero field.
        let (motion, maximum_cell_magnetic_peclet) =
            if velocity_m_s.iter().flatten().all(|v| *v == 0.0) {
                (Csr::from_parts(n, n, vec![0; n + 1], vec![], vec![]), 0.0)
            } else {
                // Only the homogeneous motion action is needed for AC.
                let (motion, _, peclet) = p.motion(
                    g,
                    velocity_m_s,
                    &vec![[0.0; 3]; mesh.tets.len()],
                    &problem.permeability_bound,
                    controls.maximum_cell_magnetic_peclet,
                    cancelled,
                )?;
                (motion, peclet)
            };
        let transpose = ops::transpose(&motion);
        let preconditioner = HarmonicJacobi::new(&p.magnetic, &p.mass, &motion, omega)?;
        let rhs_norm = finite(rhs.iter().fold(0.0_f64, |s, v| s.hypot(*v)))?;
        let target = controls
            .absolute_residual_a
            .max(finite(controls.relative_tolerance * rhs_norm)?);
        let mut x = seed.map_or_else(|| vec![0.0; 2 * n], <[f64]>::to_vec);
        let report;
        {
            let cancel = RefCell::new(&mut *cancelled);
            let stopped = Cell::new(false);
            let check = || {
                if stopped.get() {
                    true
                } else {
                    let yes = (*cancel.borrow_mut())();
                    stopped.set(yes);
                    yes
                }
            };
            let operator = HarmonicOp {
                magnetic: &p.magnetic,
                mass: &p.mass,
                motion: &motion,
                transpose: &transpose,
                omega,
                check: &check,
            };
            let residual_rhs: std::borrow::Cow<'_, [f64]> = if seed.is_some() {
                let mut image = vec![0.0; 2 * n];
                operator.apply(&x, &mut image);
                if check() {
                    return Err(SpatialMagneticError::Cancelled);
                }
                for (value, b) in image.iter_mut().zip(&rhs) {
                    *value = finite(b - *value)?;
                }
                std::borrow::Cow::Owned(image)
            } else {
                std::borrow::Cow::Borrowed(&rhs)
            };
            let initial_residual = finite(
                residual_rhs
                    .iter()
                    .fold(0.0_f64, |r, value| r.hypot(*value)),
            )?;
            let mut solved = None;
            if initial_residual > target {
                let mut solver = FgmresState::new(&residual_rhs, restart);
                let tolerance = finite((0.25 * target / initial_residual).min(0.5))?;
                if !positive(tolerance) {
                    return Err(invalid(
                        "unrepresentable harmonic induction residual target",
                    ));
                }
                let outcome = solver.run_controlled(
                    &operator,
                    &preconditioner,
                    &residual_rhs,
                    tolerance,
                    controls.maximum_cycles,
                    &mut || check(),
                );
                if outcome.progress != SolverRunProgress::Complete || check() {
                    return Err(SpatialMagneticError::Cancelled);
                }
                for (value, correction) in x.iter_mut().zip(solver.x) {
                    *value = finite(*value + correction)?;
                }
                solved = Some(outcome.report);
            }
            let mut ax = vec![0.0; 2 * n];
            operator.apply(&x, &mut ax);
            if check() {
                return Err(SpatialMagneticError::Cancelled);
            }
            let residual = finite(rhs.iter().zip(ax).fold(0.0_f64, |r, (b, a)| r.hypot(b - a)))?;
            let scale = if rhs_norm == 0.0 { 1.0 } else { rhs_norm };
            let mut history = solved.as_ref().map_or_else(Vec::new, |r| {
                r.history
                    .iter()
                    .map(|value| value * (initial_residual / scale))
                    .collect()
            });
            history.push(residual / scale);
            report = SolveReport::from_claim(
                solved.as_ref().map_or(0, |r| r.iters),
                ResidualClaim::TrueEuclidean(residual / scale),
                target / scale,
                history,
            );
            if !report.converged_euclidean() {
                return Err(SpatialMagneticError::Nonconvergence(report));
            }
        }
        let mut a_re = vec![0.0; mesh.edges.len()];
        let mut a_im = a_re.clone();
        for (edge, index) in p.edge_unknown.iter().enumerate() {
            if let Some(i) = *index {
                a_re[edge] = finite(x[i])?;
                a_im[edge] = finite(x[n + i])?;
            }
        }
        let b_re = g.flux(&a_re, cancelled)?;
        let b_im = g.flux(&a_im, cancelled)?;
        // Reuse the exact original current/force/quadrature owner. Its unit
        // interval receives RATES, so returned heat/work numbers are powers.
        let rates_re: Vec<_> = x[n..]
            .iter()
            .map(|v| finite(omega * v))
            .collect::<Result<_>>()?;
        let rates_im: Vec<_> = x[..n]
            .iter()
            .map(|v| finite(-omega * v))
            .collect::<Result<_>>()?;
        let derivative_re: Vec<_> = a_im
            .iter()
            .map(|v| finite(omega * v))
            .collect::<Result<_>>()?;
        let derivative_im: Vec<_> = a_re
            .iter()
            .map(|v| finite(-omega * v))
            .collect::<Result<_>>()?;
        let current_re = p.recover(
            g,
            &derivative_re,
            &rates_re,
            &b_re,
            velocity_m_s,
            1.0,
            cancelled,
        )?;
        let current_im = p.recover(
            g,
            &derivative_im,
            &rates_im,
            &b_im,
            velocity_m_s,
            1.0,
            cancelled,
        )?;
        let effort_re = apply(&problem.stiffness, &a_re, cancelled)?;
        let effort_im = apply(&problem.stiffness, &a_im, cancelled)?;
        let mut maximum_ampere_residual_a = 0.0_f64;
        for edge in 0..mesh.edges.len() {
            poll(cancelled)?;
            if !p.fixed_edges[edge] {
                let real =
                    finite(effort_re[edge] - load_re[edge] - current_re.edge_current_a[edge])?;
                let imaginary =
                    finite(effort_im[edge] - load_im[edge] - current_im.edge_current_a[edge])?;
                maximum_ampere_residual_a =
                    maximum_ampere_residual_a.max(finite(real.hypot(imaginary))?);
            }
        }
        let mut maximum_charge_residual_a = 0.0_f64;
        for (real, imaginary) in current_re
            .charge_current_a
            .iter()
            .zip(&current_im.charge_current_a)
        {
            maximum_charge_residual_a =
                maximum_charge_residual_a.max(finite(real.hypot(*imaginary))?);
        }
        if maximum_ampere_residual_a > controls.current_tolerance_a
            || maximum_charge_residual_a > controls.current_tolerance_a
        {
            return Err(invalid(
                "harmonic full Ampere or floating charge equation failed",
            ));
        }
        let flux_re = apply(&problem.curl, &a_re, cancelled)?;
        let flux_im = apply(&problem.curl, &a_im, cancelled)?;
        for row in mesh.d2().rows {
            poll(cancelled)?;
            let real = finite(row.iter().map(|&(f, s)| f64::from(s) * flux_re[f]).sum())?;
            let imaginary = finite(row.iter().map(|&(f, s)| f64::from(s) * flux_im[f]).sum())?;
            if finite(real.hypot(imaginary))? > controls.flux_tolerance_wb {
                return Err(invalid(
                    "harmonic induction violates discrete div-curl allowance",
                ));
            }
        }
        let mut mean_magnetic_energy_j = 0.0;
        for cell in 0..mesh.tets.len() {
            poll(cancelled)?;
            let mu = problem.permeability[cell];
            let energy = 0.25
                * g.geometry.geometry.vol_signed[cell].abs()
                * (dot(b_re[cell], mu.field_strength_a_m(b_re[cell])?)
                    + dot(b_im[cell], mu.field_strength_a_m(b_im[cell])?));
            mean_magnetic_energy_j = finite(mean_magnetic_energy_j + energy)?;
        }
        let mean_ohmic_power_w = finite(0.5 * current_re.heat_j + 0.5 * current_im.heat_j)?;
        let mean_mechanical_power_w =
            finite(0.5 * current_re.mechanical_work_j + 0.5 * current_im.mechanical_work_j)?;
        let mut reactive_mechanical_power_var = 0.0;
        let mut cells = Vec::with_capacity(p.cells.len());
        for (real, imaginary) in current_re.cells.into_iter().zip(current_im.cells) {
            poll(cancelled)?;
            let velocity = mesh.tets[real.cell].map(|v| velocity_m_s[v as usize]);
            // Im(conj(J) cross B) = Jr cross Bi - Ji cross Br.
            let quadrature_force = std::array::from_fn(|k| {
                sub(
                    cross(real.current_density_a_m2[k], b_im[real.cell]),
                    cross(imaginary.current_density_a_m2[k], b_re[real.cell]),
                )
            });
            reactive_mechanical_power_var = finite(
                reactive_mechanical_power_var
                    + 0.5 * product(velocity, quadrature_force, real.volume_m3),
            )?;
            cells.push(SpatialHarmonicInductionCell {
                cell: real.cell,
                component: real.component,
                volume_m3: real.volume_m3,
                potential_v: std::array::from_fn(|k| {
                    C64::new(real.potential_v[k], imaginary.potential_v[k])
                }),
                laboratory_electric_field_v_m: pair(
                    real.laboratory_electric_field_v_m,
                    imaginary.laboratory_electric_field_v_m,
                ),
                conductor_electric_field_v_m: pair(
                    real.conductor_electric_field_v_m,
                    imaginary.conductor_electric_field_v_m,
                ),
                current_density_a_m2: pair(
                    real.current_density_a_m2,
                    imaginary.current_density_a_m2,
                ),
                mean_lorentz_force_density_n_m3: std::array::from_fn(|k| {
                    std::array::from_fn(|j| {
                        f64::midpoint(
                            real.lorentz_force_density_n_m3[k][j],
                            imaginary.lorentz_force_density_n_m3[k][j],
                        )
                    })
                }),
                mean_ohmic_power_w: finite(0.5 * real.heat_j + 0.5 * imaginary.heat_j)?,
                mean_mechanical_power_w: finite(
                    0.5 * real.mechanical_work_j + 0.5 * imaginary.mechanical_work_j,
                )?,
            });
        }
        let mut mean_supplied_power_w = 0.0;
        let mut reactive_supplied_power_var = 0.0;
        for edge in 0..mesh.edges.len() {
            poll(cancelled)?;
            mean_supplied_power_w = finite(
                mean_supplied_power_w
                    + 0.5
                        * (load_re[edge] * derivative_re[edge]
                            + load_im[edge] * derivative_im[edge]),
            )?;
            reactive_supplied_power_var = finite(
                reactive_supplied_power_var
                    + 0.5
                        * (load_re[edge] * derivative_im[edge]
                            - load_im[edge] * derivative_re[edge]),
            )?;
        }
        let power_residual_w =
            finite(mean_ohmic_power_w + mean_mechanical_power_w - mean_supplied_power_w)?;
        let reactive_power_residual_var = finite(
            reactive_supplied_power_var - reactive_mechanical_power_var
                + 2.0 * omega * mean_magnetic_energy_j,
        )?;
        if power_residual_w.abs() > controls.power_tolerance_w
            || reactive_power_residual_var.abs() > controls.power_tolerance_w
        {
            return Err(invalid(
                "harmonic magnetic, Ohmic, Lorentz and impressed power do not close",
            ));
        }
        poll(cancelled)?;
        Ok(SpatialHarmonicInductionResponse {
            source: self,
            omega_rad_s: omega,
            solution: x,
            potential_wb: a_re
                .into_iter()
                .zip(a_im)
                .map(|(a, b)| C64::new(a, b))
                .collect(),
            flux_density_t: b_re
                .into_iter()
                .zip(b_im)
                .map(|(a, b)| std::array::from_fn(|j| C64::new(a[j], b[j])))
                .collect(),
            cells,
            balance: SpatialHarmonicInductionBalance {
                mean_magnetic_energy_j,
                mean_ohmic_power_w,
                mean_mechanical_power_w,
                reactive_mechanical_power_var,
                mean_supplied_power_w,
                reactive_supplied_power_var,
                power_residual_w,
                reactive_power_residual_var,
            },
            report,
            maximum_ampere_residual_a,
            maximum_charge_residual_a,
            maximum_cell_magnetic_peclet,
            maximum_cell_skin_ratio,
        })
    }
}

fn pair(real: [[f64; 3]; 4], imaginary: [[f64; 3]; 4]) -> [[C64; 3]; 4] {
    std::array::from_fn(|k| std::array::from_fn(|j| C64::new(real[k][j], imaginary[k][j])))
}

fn admit<C: SpatialConductorMaterial>(
    state: &SpatialInduction<'_, C>,
    current: &[[C64; 3]],
    velocity: &[[f64; 3]],
    omega: f64,
    c: SpatialHarmonicInductionControls,
) -> Result<usize> {
    let mesh = state.problem().mesh();
    if !positive(omega)
        || !positive(c.relative_tolerance)
        || c.relative_tolerance >= 1.0
        || [
            c.absolute_residual_a,
            c.current_tolerance_a,
            c.flux_tolerance_wb,
            c.power_tolerance_w,
            c.maximum_cell_magnetic_peclet,
            c.maximum_cell_skin_ratio,
        ]
        .into_iter()
        .any(|v| !positive(v))
        || c.restart == 0
        || c.maximum_cycles == 0
        || current.len() != mesh.tets.len()
        || current
            .iter()
            .flatten()
            .any(|v| !v.re.is_finite() || !v.im.is_finite())
        || velocity.len() != mesh.vertex_count
        || velocity.iter().flatten().any(|v| !v.is_finite())
    {
        return Err(invalid(
            "finite peak sources, motion, positive frequency and harmonic controls required",
        ));
    }
    if state
        .problem()
        .materials
        .iter()
        .any(|m| !matches!(m, SpatialMagneticMaterial::Linear(_)))
    {
        return Err(invalid(
            "harmonic induction requires linear magnetic sources without recoil envelopes",
        ));
    }
    let n = state
        .prepared
        .mass
        .nrows()
        .checked_mul(2)
        .ok_or_else(|| invalid("harmonic induction dimension overflow"))?;
    let restart = c.restart.min(n.max(1));
    let memory = restart
        .checked_mul(2)
        .and_then(|v| v.checked_add(12))
        .and_then(|v| v.checked_mul(n))
        .and_then(|v| {
            restart
                .checked_add(8)
                .and_then(|r| r.checked_mul(restart))
                .and_then(|r| v.checked_add(r))
        })
        .and_then(|v| {
            c.maximum_cycles
                .checked_mul(8)
                .and_then(|r| v.checked_add(r))
        });
    let fields = n
        .checked_mul(24)
        .and_then(|v| {
            mesh.edges
                .len()
                .checked_mul(16)
                .and_then(|r| v.checked_add(r))
        })
        .and_then(|v| {
            mesh.tets
                .len()
                .checked_mul(512)
                .and_then(|r| v.checked_add(r))
        });
    if memory.is_none_or(|v| v > c.maximum_krylov_scalars)
        || fields.is_none_or(|v| v > c.maximum_field_scalars)
        || restart.checked_mul(c.maximum_cycles).is_none()
    {
        return Err(invalid(
            "harmonic induction field/Krylov work or storage budget exceeded",
        ));
    }
    Ok(restart)
}

/// Fused [[K-V, omega M], [-omega M, K-V]] application. No whole-field
/// temporary vectors or allocations occur in a Krylov application.
struct HarmonicOp<'a> {
    magnetic: &'a Csr,
    mass: &'a Csr,
    motion: &'a Csr,
    transpose: &'a Csr,
    omega: f64,
    check: &'a dyn Fn() -> bool,
}
impl HarmonicOp<'_> {
    fn action(&self, x: &[f64], y: &mut [f64], motion: &Csr, omega: f64) {
        let n = self.mass.nrows();
        for i in 0..n {
            if (self.check)() {
                y.fill(f64::NAN);
                return;
            }
            let row = |matrix: &Csr| {
                let (columns, values) = matrix.row(i);
                columns
                    .iter()
                    .zip(values)
                    .fold((0.0, 0.0), |(real, imaginary), (&j, v)| {
                        (real + v * x[j], imaginary + v * x[n + j])
                    })
            };
            let (kr, ki) = row(self.magnetic);
            let (vr, vi) = row(motion);
            let (mr, mi) = row(self.mass);
            y[i] = kr - vr + omega * mi;
            y[n + i] = ki - vi - omega * mr;
        }
    }
}
impl LinearOp for HarmonicOp<'_> {
    fn n(&self) -> usize {
        2 * self.mass.nrows()
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        self.action(x, y, self.motion, self.omega);
    }
    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        self.action(x, y, self.transpose, -self.omega);
    }
}

/// Exact inverse of each complex diagonal, including electric-potential rows.
struct HarmonicJacobi(Vec<C64>);
impl HarmonicJacobi {
    fn new(k: &Csr, m: &Csr, v: &Csr, omega: f64) -> Result<Self> {
        let entries = (0..m.nrows())
            .map(|i| {
                let real = finite(diagonal(k, i) - diagonal(v, i))?;
                let imaginary = finite(omega * diagonal(m, i))?;
                let scale = real.abs().max(imaginary.abs());
                if !positive(scale) {
                    return Err(invalid(
                        "nonzero representable harmonic Jacobi scale required",
                    ));
                }
                let (a, b) = (real / scale, imaginary / scale);
                let square = a * a + b * b;
                Ok(C64::new(
                    finite((a / square) / scale)?,
                    finite((b / square) / scale)?,
                ))
            })
            .collect::<Result<_>>()?;
        Ok(Self(entries))
    }
}
impl FlexiblePreconditioner for HarmonicJacobi {
    fn apply(&self, _: usize, residual: &[f64], output: &mut [f64]) {
        let n = self.0.len();
        for (i, z) in self.0.iter().enumerate() {
            output[i] = z.re * residual[i] - z.im * residual[n + i];
            output[n + i] = z.im * residual[i] + z.re * residual[n + i];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn harmonic_motion_operator_has_the_actual_signed_real_block_transpose() {
        let k = Csr::from_parts(
            2,
            2,
            vec![0, 2, 4],
            vec![0, 1, 0, 1],
            vec![4.0, -1.0, -1.0, 3.0],
        );
        let m = Csr::from_parts(
            2,
            2,
            vec![0, 2, 4],
            vec![0, 1, 0, 1],
            vec![2.0, 0.5, 0.5, 1.0],
        );
        let v = Csr::from_parts(
            2,
            2,
            vec![0, 2, 4],
            vec![0, 1, 0, 1],
            vec![0.2, 0.7, -0.4, -0.3],
        );
        let transpose = ops::transpose(&v);
        let op = HarmonicOp {
            magnetic: &k,
            mass: &m,
            motion: &v,
            transpose: &transpose,
            omega: 2.0,
            check: &|| false,
        };
        let dense = [
            [3.8, -1.7, 4.0, 1.0],
            [-0.6, 3.3, 1.0, 2.0],
            [-4.0, -1.0, 3.8, -1.7],
            [-1.0, -2.0, -0.6, 3.3],
        ];
        let x = [0.2, -0.5, 1.3, -0.7];
        let mut applied = [0.0; 4];
        let mut transposed = [0.0; 4];
        op.apply(&x, &mut applied);
        op.apply_transpose(&x, &mut transposed);
        for i in 0..4 {
            let a = (0..4).map(|j| dense[i][j] * x[j]).sum::<f64>();
            let at = (0..4).map(|j| dense[j][i] * x[j]).sum::<f64>();
            assert!((applied[i] - a).abs() < 2e-14);
            assert!((transposed[i] - at).abs() < 2e-14);
        }
        let jacobi = HarmonicJacobi::new(&k, &m, &v, 2.0).unwrap();
        let mut inverse = [0.0; 4];
        jacobi.apply(0, &x, &mut inverse);
        for i in 0..2 {
            let diagonal = C64::new(dense[i][i], dense[i + 2][i]);
            let restored = diagonal * C64::new(inverse[i], inverse[i + 2]);
            assert!((restored.re - x[i]).abs() < 2e-14);
            assert!((restored.im - x[i + 2]).abs() < 2e-14);
        }
    }
}
