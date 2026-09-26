/// Reusable discretization and material. The background grid is built once,
/// including when a guarded iteration backtracks over several geometries.
use super::*;

struct ComplianceKernel {
    grid: Quadtree,
    material: IsotropicElastic,
    lambda: f64,
    mu: f64,
    support: EdgeBand,
    fixture: Cantilever,
    settings: OptimizeSettings,
    mass: fs_sparse::Csr,
    stiffness: fs_sparse::Csr,
}

struct EvaluatedDesign {
    phi: GridSdf,
    solution: CutElasticitySolution,
    compliance: f64,
    volume: f64,
}

struct ShapeDirection {
    smooth: Vec<f64>,
    mean_energy: f64,
    topological: Option<Vec<f64>>,
}

struct GeometryMove {
    normal_multiplier: f64,
    hole_multiplier: f64,
    scale: f64,
    with_holes: bool,
}

struct GeometryTrial {
    phi: GridSdf,
    audit: RedistanceAudit,
    events: Vec<NucleationEvent>,
    load_pad_nodes: usize,
}

/// Validate evolution controls before allocating a background grid or moving φ.
fn validate_evolution(phi: &GridSdf, settings: OptimizeSettings) -> Result<(), CutFemError> {
    let n = 1u32.checked_shl(settings.level).ok_or_else(|| {
        invalid_input("optimizer grid level exceeds the u32 lattice address space")
    })?;
    if usize::try_from(n).ok() != Some(phi.n()) {
        return Err(invalid_input("SDF lattice must match the CutFEM grid"));
    }
    if phi.nodes().iter().any(|v| !v.is_finite()) {
        return Err(invalid_input("optimizer level set must contain only finite nodal values"));
    }
    if !(settings.volfrac.is_finite() && settings.volfrac > 0.0 && settings.volfrac <= 1.0) {
        return Err(invalid_input("optimizer volume fraction must lie in (0, 1]"));
    }
    if !(settings.band_cells.is_finite() && settings.band_cells > 0.0
        && settings.move_cells.is_finite() && settings.move_cells >= 0.0
        && settings.move_cells <= 0.5 * settings.band_cells)
    {
        return Err(invalid_input("optimizer move must be finite, nonnegative, and at most half the positive band width"));
    }
    if !(settings.ell0.is_finite() && settings.ell0 >= 0.0
        && settings.mu_al.is_finite() && settings.mu_al > 0.0
        && settings.sobolev_alpha.is_finite() && settings.sobolev_alpha >= 0.0)
    {
        return Err(invalid_input("optimizer multiplier, penalty, and smoothing controls are invalid"));
    }
    if settings.nucleation_period > 0
        && !(settings.hole_radius_cells.is_finite() && settings.hole_radius_cells > 0.0)
    {
        return Err(invalid_input("enabled nucleation requires a finite positive hole radius"));
    }
    Ok(())
}

impl ComplianceKernel {
    fn new(phi: &GridSdf, fixture: Cantilever, settings: OptimizeSettings) -> Result<Self, CutFemError> {
        let support = cantilever_support(fixture)?;
        let (material, lambda, mu) = validated_plane_strain_material(settings)?;
        validate_evolution(phi, settings)?;
        let grid = Quadtree::uniform(settings.level);
        let (mass, stiffness) = mass_stiffness(phi.n());
        Ok(Self { grid, material, lambda, mu, support, fixture, settings, mass, stiffness })
    }

    fn evaluate(&self, phi: GridSdf) -> Result<EvaluatedDesign, CutFemError> {
        if phi.nodes().iter().any(|v| !v.is_finite()) {
            return Err(invalid_input("evolution produced a non-finite level set"));
        }
        let clamp = |x: f64, _y: f64| x < 1e-9;
        let traction = |_: f64, _: f64| [0.0, -self.fixture.load];
        let solver = CutElasticity {
            grid: &self.grid,
            sdf: &phi,
            material: &self.material,
            nitsche_beta: 20.0,
            ghost_gamma: 0.5,
            stabilization_scaling: fs_cutfem::CutStabilizationScaling::LongitudinalModulus,
            quad_depth: 2,
            clamp: Some(&clamp),
            boundary_traction: None,
            traction_free_interface: true,
            solver_tol: SOLVER_TOL,
            solver_max_iters: SOLVER_MAX_ITERS,
        };
        let solution = solver.solve_with_boundary_traction(
            &|_, _| [0.0, 0.0],
            &|_, _| [0.0, 0.0],
            BoundaryTraction::EdgeBand { support: self.support, value: &traction },
        )?;
        let compliance = solution.compliance();
        let volume = material_volume(&self.grid, &phi);
        if !(compliance.is_finite() && compliance >= 0.0 && volume.is_finite() && volume > 0.0) {
            return Err(invalid_input("cantilever solve produced invalid compliance or material area"));
        }
        Ok(EvaluatedDesign { phi, solution, compliance, volume })
    }

    fn direction(&self, state: &EvaluatedDesign, iteration: usize) -> Result<ShapeDirection, CutFemError> {
        let grid = &self.grid;
        let phi = &state.phi;
        let sol = &state.solution;
        let (lambda, mu) = (self.lambda, self.mu);
        let n = phi.n();
        let h = phi.h();
        let stride = n + 1;
        let mut energy = vec![0.0f64; stride * stride];
        let mut seeded = vec![false; stride * stride];
        for j in 0..=n {
            for i in 0..=n {
                let k = i + j * stride;
                let p = phi.pos(i, j);
                // Sample slightly inside material along −∇φ.
                let g = phi.gradient_at(p);
                let norm = g[0].hypot(g[1]);
                if !norm.is_finite() {
                    return Err(invalid_input("level-set gradient overflowed at a sensitivity probe"));
                }
                let gn = norm.max(1e-12);
                let q = [
                    (p[0] - 0.75 * h * g[0] / gn).clamp(0.0, 1.0),
                    (p[1] - 0.75 * h * g[1] / gn).clamp(0.0, 1.0),
                ];
                if phi.value_at(q) > 0.0 {
                    continue;
                }
                let (eps, ok) = strain_at(grid, phi, sol, q);
                if !ok {
                    continue;
                }
                let sxx = (lambda + 2.0 * mu) * eps[0] + lambda * eps[1];
                let syy = lambda * eps[0] + (lambda + 2.0 * mu) * eps[1];
                let sxy = 2.0 * mu * eps[2];
                energy[k] = 0.5 * (sxx * eps[0] + syy * eps[1] + 2.0 * sxy * eps[2]);
                if !energy[k].is_finite() {
                    return Err(invalid_input("strain-energy sensitivity overflowed"));
                }
                seeded[k] = phi.node(i, j).abs() <= 2.0 * h;
            }
        }
        // 4. Extend off the interface, smooth, form v_n = w − ℓ.
        extend_velocity(phi, &mut energy, &seeded);
        if energy.iter().any(|value| !value.is_finite()) {
            return Err(invalid_input("normal velocity extension produced non-finite sensitivities"));
        }
        let (smooth, _iters) = fs_adjoint::sobolev::sobolev_smooth(
            &self.mass,
            &self.stiffness,
            self.settings.sobolev_alpha * h * h,
            &energy,
            1e-10,
        );
        if smooth.iter().any(|v| !v.is_finite()) {
            return Err(invalid_input("Sobolev smoothing produced non-finite shape sensitivities"));
        }
        #[allow(clippy::cast_precision_loss)]
        let mean_energy = smooth.iter().sum::<f64>() / smooth.len() as f64;
        if !mean_energy.is_finite() {
            return Err(invalid_input("mean shape energy is not finite"));
        }
        // Both shape and hole sensitivities belong to this exact solved φ.
        // Never probe the old displacement field on an already-advected design.
        let topological = if self.settings.nucleation_period > 0 && iteration > 0
            && iteration % self.settings.nucleation_period == 0
        {
            let mut dt_field = vec![f64::INFINITY; stride * stride];
            for j in 0..=n {
                for i in 0..=n {
                    let k = i + j * stride;
                    let p = phi.pos(i, j);
                    if phi.value_at(p) > -2.0 * h {
                        continue;
                    }
                    let (eps, ok) = strain_at(grid, phi, sol, p);
                    if !ok {
                        continue;
                    }
                    let sxx = (lambda + 2.0 * mu) * eps[0] + lambda * eps[1];
                    let syy = lambda * eps[0] + (lambda + 2.0 * mu) * eps[1];
                    let sxy = 2.0 * mu * eps[2];
                    let derivative = topological_derivative(lambda, mu, [sxx, syy, sxy], eps);
                    if !derivative.is_finite() {
                        return Err(invalid_input("hole sensitivity overflowed on an active material node"));
                    }
                    dt_field[k] = derivative;
                }
            }
            if dt_field.iter().any(|v| v.is_nan() || *v == f64::NEG_INFINITY) {
                return Err(invalid_input("non-finite hole sensitivity on an active material node"));
            }
            Some(dt_field)
        } else {
            None
        };
        Ok(ShapeDirection { smooth, mean_energy, topological })
    }

    fn propose(
        &self,
        state: &EvaluatedDesign,
        direction: &ShapeDirection,
        movement: GeometryMove,
    ) -> Result<GeometryTrial, CutFemError> {
        let settings = self.settings;
        let h = state.phi.h();
        let mut phi = state.phi.clone();
        let vn: Vec<f64> = direction.smooth.iter().map(|w| w - movement.normal_multiplier).collect();
        if vn.iter().any(|v| !v.is_finite()) {
            return Err(invalid_input("shape velocity overflowed"));
        }
        let vmax = vn.iter().fold(0.0f64, |m, v| m.max(v.abs())).max(1e-12);
        let band = build_band(&phi, settings.band_cells);
        let duration = movement.scale * settings.move_cells * h / vmax;
        if !duration.is_finite() {
            return Err(invalid_input("shape advection duration overflowed"));
        }
        advect(&mut phi, &band, &Velocity::Normal(&vn), duration, 0.45);
        let mut load_pad_nodes = retain_cantilever_load_pad(&mut phi, self.support);
        let mut audit = redistance(&mut phi, settings.band_cells);
        load_pad_nodes += retain_cantilever_load_pad(&mut phi, self.support);
        let mut events = Vec::new();
        if let Some(dt_field) = direction.topological.as_ref().filter(|_| movement.with_holes) {
            events = nucleate(
                &mut phi, dt_field, movement.hole_multiplier,
                settings.hole_radius_cells * h,
                6.0 * settings.hole_radius_cells * h, 2,
            );
            if !events.is_empty() {
                load_pad_nodes += retain_cantilever_load_pad(&mut phi, self.support);
                let first_drift = audit.interface_drift_h;
                audit = redistance(&mut phi, settings.band_cells);
                // Retain both redistancing diagnostics; their sum is not a
                // bound on the complete update, which also inserts holes.
                audit.interface_drift_h += first_drift;
                load_pad_nodes += retain_cantilever_load_pad(&mut phi, self.support);
            }
        }
        if phi.nodes().iter().any(|v| !v.is_finite()) || !audit.interface_drift_h.is_finite() {
            return Err(invalid_input("level-set evolution or redistancing audit is non-finite"));
        }
        Ok(GeometryTrial { phi, audit, events, load_pad_nodes })
    }
}

fn append_iteration(
    report: &mut OptimizeReport,
    state: &EvaluatedDesign,
    iteration: usize,
    ell: f64,
    audit: RedistanceAudit,
    events: Vec<NucleationEvent>,
    load_pad_nodes: usize,
) {
    let snap = fnv(&state.phi);
    let mut row = String::new();
    let compliance = state.compliance;
    let volume = state.volume;
    let _ = write!(
        row,
        "{{\"iter\":{iteration},\"compliance\":{compliance:.17e},\"volume\":{volume:.17e},\
         \"ell\":{ell:.17e},\"drift_h\":{:.17e},\"load_pad_nodes\":{load_pad_nodes},\
         \"snapshot\":\"{snap:#018x}\"}}",
        audit.interface_drift_h
    );
    report.rows.push(row);
    report.compliance.push(compliance);
    report.volume.push(volume);
    report.ell.push(ell);
    report.audits.push(audit);
    report.events.extend(events);
    report.snapshots.push(snap);
    report.load_pad_nodes.push(load_pad_nodes);
}

/// Run the level-set compliance descent. Every returned row describes the
/// evaluated geometry whose snapshot appears in that row, including the final
/// geometry after nucleation. Evolution occurs on a trial copy; a refused solve
/// leaves `phi` at the last successfully evaluated design.
///
/// The initial state costs one solve, followed by one solve per iteration; the
/// accepted solution is reused for the next iteration's sensitivities. This
/// function preserves the legacy fixed-step policy; it does not enforce descent.
///
/// # Errors
/// Invalid load/material/evolution settings refuse before mutation. Physics
/// failures propagate as [`CutFemError`]; no uncomputed candidate replaces `phi`.
pub fn optimize_compliance(
    phi: &mut GridSdf,
    fixture: Cantilever,
    settings: OptimizeSettings,
) -> Result<OptimizeReport, CutFemError> {
    // Empty runs still admit inputs but do not allocate or solve a PDE.
    cantilever_support(fixture)?;
    validated_plane_strain_material(settings)?;
    validate_evolution(phi, settings)?;
    let mut report = OptimizeReport::default();
    if settings.iterations == 0 {
        return Ok(report);
    }
    let kernel = ComplianceKernel::new(phi, fixture, settings)?;
    let mut current = kernel.evaluate(phi.clone())?;
    let mut ell = settings.ell0;
    for iteration in 0..settings.iterations {
        let direction = kernel.direction(&current, iteration)?;
        let GeometryTrial { phi: trial_phi, audit, events, load_pad_nodes } =
            kernel.propose(&current, &direction, GeometryMove {
                normal_multiplier: ell, hole_multiplier: ell, scale: 1.0, with_holes: true,
            })?;
        let candidate = kernel.evaluate(trial_phi)?;
        let next_ell = ell + settings.mu_al * direction.mean_energy.abs().max(1e-30)
            * (candidate.volume - settings.volfrac) / settings.volfrac;
        if !next_ell.is_finite() {
            return Err(invalid_input("volume multiplier update overflowed"));
        }
        ell = next_ell.max(0.0);
        append_iteration(&mut report, &candidate, iteration, ell, audit, events, load_pad_nodes);
        *phi = candidate.phi.clone();
        current = candidate;
    }
    Ok(report)
}

/// Strain at a point from the CutFEM solution (bilinear gradient on
/// the containing cell); `ok = false` outside the active mesh.
fn strain_at(
    grid: &Quadtree,
    phi: &GridSdf,
    sol: &CutElasticitySolution,
    p: [f64; 2],
) -> ([f64; 3], bool) {
    let level = grid.max_level();
    let nf = f64::from(1u32 << level);
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let ci = ((p[0] * nf).floor().clamp(0.0, nf - 1.0)) as u32;
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let cj = ((p[1] * nf).floor().clamp(0.0, nf - 1.0)) as u32;
    let cell = (level, ci, cj);
    let (lo, hi) = grid.rect(cell);
    let corners = grid.corner_nodes(cell);
    let nodal = sol.nodal();
    let mut vals = [[0.0f64; 2]; 4];
    for (a, c) in corners.iter().enumerate() {
        match nodal.get(c) {
            Some(u) => vals[a] = *u,
            None => return ([0.0; 3], false),
        }
    }
    let _ = phi;
    let hx = hi[0] - lo[0];
    let hy = hi[1] - lo[1];
    let xi = ((p[0] - lo[0]) / hx).clamp(0.0, 1.0);
    let et = ((p[1] - lo[1]) / hy).clamp(0.0, 1.0);
    let g = [
        [-(1.0 - et) / hx, -(1.0 - xi) / hy],
        [(1.0 - et) / hx, -xi / hy],
        [et / hx, xi / hy],
        [-et / hx, (1.0 - xi) / hy],
    ];
    let mut gu = [[0.0f64; 2]; 2];
    for a in 0..4 {
        for c in 0..2 {
            gu[c][0] += g[a][0] * vals[a][c];
            gu[c][1] += g[a][1] * vals[a][c];
        }
    }
    (
        [gu[0][0], gu[1][1], f64::midpoint(gu[0][1], gu[1][0])],
        true,
    )
}

#[cfg(test)]
mod shape_derivative_gate {
    //! Per-iteration check of the optimizer's shape derivative against central
    //! differences of the discrete CutFEM compliance (bead q61wp.75 item 3).
    use super::*;

    /// Continuum prediction dJ/de for phi -> phi - e*s: the boundary moves out
    /// of the material by e*s/|grad phi|, so dJ/de = -int_G 2 w s/|grad phi| ds,
    /// with w = sigma:eps/2 sampled 0.75h inside, as the optimizer samples it.
    fn predicted(kernel: &ComplianceKernel, state: &EvaluatedDesign, s: &[f64], sub: usize) -> f64 {
        let phi = &state.phi;
        let (n, h) = (phi.n(), phi.h());
        let stride = n + 1;
        let (lambda, mu) = (kernel.lambda, kernel.mu);
        let nodal = |p: [f64; 2]| {
            // Bilinear interpolation of the nodal velocity field.
            let x = (p[0] / h).clamp(0.0, n as f64 - 1e-12);
            let y = (p[1] / h).clamp(0.0, n as f64 - 1e-12);
            let (i, j) = (x.floor() as usize, y.floor() as usize);
            let (a, b) = (x - i as f64, y - j as f64);
            let v = |i: usize, j: usize| s[i + j * stride];
            v(i, j) * (1.0 - a) * (1.0 - b) + v(i + 1, j) * a * (1.0 - b) + v(i, j + 1) * (1.0 - a) * b + v(i + 1, j + 1) * a * b
        };
        let energy = |q: [f64; 2]| {
            let g = phi.gradient_at(q);
            let gn = g[0].hypot(g[1]).max(1e-12);
            let r = [(q[0] - 0.75 * h * g[0] / gn).clamp(0.0, 1.0), (q[1] - 0.75 * h * g[1] / gn).clamp(0.0, 1.0)];
            let (eps, ok) = strain_at(&kernel.grid, phi, &state.solution, r);
            if !ok { return None; }
            let sxx = (lambda + 2.0 * mu) * eps[0] + lambda * eps[1];
            let syy = lambda * eps[0] + (lambda + 2.0 * mu) * eps[1];
            let sxy = 2.0 * mu * eps[2];
            Some((0.5 * (sxx * eps[0] + syy * eps[1] + 2.0 * sxy * eps[2]), gn))
        };
        let m = n * sub;
        let d = 1.0 / m as f64;
        let mut total = 0.0;
        for j in 0..m {
            for i in 0..m {
                let c = [[i, j], [i + 1, j], [i + 1, j + 1], [i, j + 1]].map(|[a, b]| [a as f64 * d, b as f64 * d]);
                let v = c.map(|p| phi.value_at(p));
                let mut crossings = Vec::new();
                for e in 0..4 {
                    let (p, q, fp, fq) = (c[e], c[(e + 1) % 4], v[e], v[(e + 1) % 4]);
                    if (fp < 0.0) != (fq < 0.0) {
                        let t = fp / (fp - fq);
                        crossings.push([p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])]);
                    }
                }
                for pair in crossings.chunks_exact(2) {
                    let (a, b) = (pair[0], pair[1]);
                    let mid = [0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1])];
                    // The loaded right edge and clamped left edge are not free boundary.
                    if mid[0] < 1e-9 || mid[0] > 1.0 - 1e-9 { continue; }
                    if let Some((w, gn)) = energy(mid) {
                        total += -2.0 * w * nodal(mid) / gn * (b[0] - a[0]).hypot(b[1] - a[1]);
                    }
                }
            }
        }
        total
    }

    fn measured(kernel: &ComplianceKernel, phi: &GridSdf, s: &[f64], e: f64) -> f64 {
        let shifted = |sign: f64| {
            let mut p = phi.clone();
            for (v, sk) in p.nodes_mut().iter_mut().zip(s) { *v -= sign * e * sk; }
            kernel.evaluate(p).expect("perturbed design solves").compliance
        };
        (shifted(1.0) - shifted(-1.0)) / (2.0 * e)
    }

    fn start(level: u32, iterations: usize) -> (ComplianceKernel, EvaluatedDesign, OptimizeSettings) {
        let n = 1usize << level;
        let settings = OptimizeSettings { level, iterations, nucleation_period: 0, ..OptimizeSettings::default() };
        let phi = GridSdf::from_fn(n, &|_, y| (y - 0.5).abs() - 0.35);
        let kernel = ComplianceKernel::new(&phi, Cantilever { load: 1.0, band: 0.125 }, settings).unwrap();
        let current = kernel.evaluate(phi).unwrap();
        (kernel, current, settings)
    }

    /// (central difference, relative FD noise, continuum prediction) along the
    /// optimizer's own smoothed velocity at this exact design.
    fn gate(kernel: &ComplianceKernel, state: &EvaluatedDesign, s: &[f64]) -> (f64, f64, f64) {
        let smax = s.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        let e = 0.1 * state.phi.h() / smax;
        let (coarse, fine) = (measured(kernel, &state.phi, s, e), measured(kernel, &state.phi, s, 0.5 * e));
        (fine, ((coarse - fine) / fine).abs(), predicted(kernel, state, s, 4))
    }

    /// Consistency band of the discrete derivative against the continuum shape
    /// derivative at level 5. Measured 2026-09-25 by the refinement test below:
    /// ratio-1 = 0.468, 0.262, 0.164 at levels 4, 5, 6 (converging at about
    /// h^0.75). The band is the level-5 value with 1.5x headroom.
    const LEVEL5_CONSISTENCY: f64 = 0.4;
    /// Measured FD noise between e and e/2 peaked at 4e-3; 1% leaves headroom.
    const FD_NOISE: f64 = 0.01;

    #[test]
    fn shape_derivative_consistency_error_shrinks_under_refinement() {
        let mut errors = Vec::new();
        for level in [4u32, 5, 6] {
            let (kernel, current, _) = start(level, 1);
            let s = kernel.direction(&current, 0).unwrap().smooth;
            let (fd, noise, pr) = gate(&kernel, &current, &s);
            assert!(noise <= FD_NOISE, "level {level}: FD noise {noise:e}");
            errors.push(fd / pr - 1.0);
        }
        println!("{{\"ratio_minus_one\":{errors:?}}}");
        assert!(errors.iter().all(|&e| e > 0.0), "{errors:?}");
        assert!(errors[1] < 0.75 * errors[0] && errors[2] < 0.75 * errors[1], "{errors:?}");
    }

    #[test]
    fn every_sampled_iteration_moves_compliance_as_its_shape_derivative_predicts() {
        let (kernel, mut current, settings) = start(5, 7);
        let mut ell = settings.ell0;
        let mut checked = 0;
        for iteration in 0..settings.iterations {
            let direction = kernel.direction(&current, iteration).unwrap();
            if iteration % 3 == 0 {
                let (fd, noise, pr) = gate(&kernel, &current, &direction.smooth);
                println!("{{\"iteration\":{iteration},\"fd\":{fd:e},\"predicted\":{pr:e},\"noise\":{noise:e}}}");
                assert!(fd < 0.0 && pr < 0.0, "iteration {iteration}: the velocity is not a descent direction");
                assert!(noise <= FD_NOISE, "iteration {iteration}: FD noise {noise:e}");
                assert!((fd / pr - 1.0).abs() <= LEVEL5_CONSISTENCY, "iteration {iteration}: ratio {}", fd / pr);
                // Falsifier: the classic factor-of-two slip (w instead of
                // sigma:eps) must fall outside the band.
                assert!((fd / (0.5 * pr) - 1.0).abs() > LEVEL5_CONSISTENCY, "the band cannot see a factor-2 error");
                checked += 1;
            }
            let trial = kernel.propose(&current, &direction, GeometryMove {
                normal_multiplier: ell, hole_multiplier: ell, scale: 1.0, with_holes: true,
            }).unwrap();
            let candidate = kernel.evaluate(trial.phi).unwrap();
            ell = (ell + settings.mu_al * direction.mean_energy.abs().max(1e-30)
                * (candidate.volume - settings.volfrac) / settings.volfrac).max(0.0);
            current = candidate;
        }
        assert_eq!(checked, 3);
    }
}
