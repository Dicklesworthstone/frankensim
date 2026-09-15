use super::*;
use std::fmt::Write as _;

struct MultiState {
    phi: GridSdf,
    solutions: Vec<CutElasticitySolution>,
    compliances: Vec<f64>,
    active: Option<usize>,
    objective: f64,
    volume: f64,
}

struct Direction {
    smooth: Vec<f64>,
    mean_energy: f64,
    topological: Option<Vec<f64>>,
}

struct Trial {
    phi: GridSdf,
    audit: RedistanceAudit,
    events: Vec<NucleationEvent>,
    load_pad_nodes: usize,
}

struct Kernel {
    grid: Quadtree,
    material: IsotropicElastic,
    lambda: f64,
    mu: f64,
    supports: Vec<EdgeBand>,
    load_cases: Vec<RobustLoadCase>,
    aggregate: RobustAggregate,
    settings: OptimizeSettings,
    mass: fs_sparse::Csr,
    stiffness: fs_sparse::Csr,
}

impl Kernel {
    fn new(
        phi: &GridSdf,
        load_cases: &[RobustLoadCase],
        settings: OptimizeSettings,
        aggregate: RobustAggregate,
    ) -> Result<Self, CutFemError> {
        let supports = validate(phi, load_cases, settings)?;
        let (material, lambda, mu) = material(settings)?;
        let grid = Quadtree::uniform(settings.level);
        let (mass, stiffness) = mass_stiffness(phi.n());
        Ok(Self {
            grid, material, lambda, mu, supports, load_cases: load_cases.to_vec(),
            aggregate, settings, mass, stiffness,
        })
    }

    fn evaluate(&self, phi: GridSdf) -> Result<MultiState, CutFemError> {
        if phi.nodes().iter().any(|value| !value.is_finite()) {
            return Err(invalid("multi-load evolution produced non-finite level-set nodes"));
        }
        let clamp = |x: f64, _y: f64| x < 1e-9;
        let solver = CutElasticity {
            grid: &self.grid,
            sdf: &phi,
            material: &self.material,
            nitsche_beta: 20.0,
            ghost_gamma: 0.5,
            stabilization_scaling: CutStabilizationScaling::LongitudinalModulus,
            quad_depth: 2,
            clamp: Some(&clamp),
            boundary_traction: None,
            traction_free_interface: true,
            solver_tol: SOLVER_TOL,
            solver_max_iters: SOLVER_MAX_ITERS,
        };
        let mut solutions = Vec::with_capacity(self.load_cases.len());
        let mut compliances = Vec::with_capacity(self.load_cases.len());
        for (case, support) in self.load_cases.iter().zip(&self.supports) {
            let value = case.traction();
            let traction = move |_: f64, _: f64| value;
            let solution = solver.solve_with_boundary_traction(
                &|_, _| [0.0, 0.0],
                &|_, _| [0.0, 0.0],
                BoundaryTraction::EdgeBand { support: *support, value: &traction },
            )?;
            let compliance = solution.compliance();
            if !(compliance.is_finite() && compliance >= 0.0) {
                return Err(invalid("multi-load case produced invalid compliance"));
            }
            compliances.push(compliance);
            solutions.push(solution);
        }
        let active = match self.aggregate {
            RobustAggregate::WeightedSum => None,
            RobustAggregate::WorstWeightedCase => Some(active_case(&compliances, &self.load_cases)),
        };
        let objective = match self.aggregate {
            RobustAggregate::WeightedSum => compliances.iter().zip(&self.load_cases)
                .map(|(&c, load)| load.weight() * c).sum::<f64>(),
            RobustAggregate::WorstWeightedCase => {
                let index = active.expect("worst-weighted mode has an active case");
                self.load_cases[index].weight() * compliances[index]
            }
        };
        let volume = material_volume(&self.grid, &phi);
        if !(objective.is_finite() && objective >= 0.0 && volume.is_finite() && volume > 0.0) {
            return Err(invalid("multi-load evaluation produced invalid aggregate or area"));
        }
        Ok(MultiState { phi, solutions, compliances, active, objective, volume })
    }

    fn direction(&self, state: &MultiState, iteration: usize) -> Result<Direction, CutFemError> {
        let phi = &state.phi;
        let n = phi.n();
        let h = phi.h();
        let stride = n + 1;
        let mut energy = vec![0.0f64; stride * stride];
        let mut seeded = vec![false; stride * stride];
        for j in 0..=n {
            for i in 0..=n {
                let k = i + j * stride;
                let p = phi.pos(i, j);
                let g = phi.gradient_at(p);
                let norm = g[0].hypot(g[1]);
                if !norm.is_finite() {
                    return Err(invalid("multi-load level-set gradient overflowed"));
                }
                let gn = norm.max(1e-12);
                let q = [
                    (p[0] - 0.75 * h * g[0] / gn).clamp(0.0, 1.0),
                    (p[1] - 0.75 * h * g[1] / gn).clamp(0.0, 1.0),
                ];
                if phi.value_at(q) > 0.0 { continue; }
                let value = match state.active {
                    Some(index) => {
                        let (eps, ok) = strain_at(&self.grid, &state.solutions[index], q);
                        if ok { self.load_cases[index].weight() * stress_energy(self.lambda, self.mu, eps).1 } else { 0.0 }
                    }
                    None => {
                        let mut total = 0.0;
                        for (solution, load) in state.solutions.iter().zip(&self.load_cases) {
                            let (eps, ok) = strain_at(&self.grid, solution, q);
                            if ok { total += load.weight() * stress_energy(self.lambda, self.mu, eps).1; }
                        }
                        total
                    }
                };
                if !value.is_finite() {
                    return Err(invalid("multi-load strain-energy sensitivity overflowed"));
                }
                energy[k] = value;
                seeded[k] = phi.node(i, j).abs() <= 2.0 * h;
            }
        }
        extend_velocity(phi, &mut energy, &seeded);
        if energy.iter().any(|value| !value.is_finite()) {
            return Err(invalid("multi-load velocity extension produced non-finite values"));
        }
        let (smooth, _) = fs_adjoint::sobolev::sobolev_smooth(
            &self.mass,
            &self.stiffness,
            self.settings.sobolev_alpha * h * h,
            &energy,
            1e-10,
        );
        if smooth.iter().any(|value| !value.is_finite()) {
            return Err(invalid("multi-load Sobolev smoothing produced non-finite values"));
        }
        #[allow(clippy::cast_precision_loss)]
        let mean_energy = smooth.iter().sum::<f64>() / smooth.len() as f64;
        if !mean_energy.is_finite() {
            return Err(invalid("multi-load mean shape energy is non-finite"));
        }

        let topological = if self.settings.nucleation_period > 0 && iteration > 0
            && iteration % self.settings.nucleation_period == 0
        {
            let mut dt_field = vec![f64::INFINITY; stride * stride];
            for j in 0..=n {
                for i in 0..=n {
                    let k = i + j * stride;
                    let p = phi.pos(i, j);
                    if phi.value_at(p) > -2.0 * h { continue; }
                    let mut value = 0.0;
                    let mut any = false;
                    match state.active {
                        Some(index) => {
                            let (eps, ok) = strain_at(&self.grid, &state.solutions[index], p);
                            if ok {
                                let (stress, _) = stress_energy(self.lambda, self.mu, eps);
                                value = self.load_cases[index].weight()
                                    * topological_derivative(self.lambda, self.mu, stress, eps);
                                any = true;
                            }
                        }
                        None => {
                            for (solution, load) in state.solutions.iter().zip(&self.load_cases) {
                                let (eps, ok) = strain_at(&self.grid, solution, p);
                                if ok {
                                    let (stress, _) = stress_energy(self.lambda, self.mu, eps);
                                    value += load.weight()
                                        * topological_derivative(self.lambda, self.mu, stress, eps);
                                    any = true;
                                }
                            }
                        }
                    }
                    if any {
                        if !value.is_finite() {
                            return Err(invalid("multi-load topological sensitivity overflowed"));
                        }
                        dt_field[k] = value;
                    }
                }
            }
            Some(dt_field)
        } else {
            None
        };
        Ok(Direction { smooth, mean_energy, topological })
    }

    fn propose(&self, state: &MultiState, direction: &Direction, ell: f64) -> Result<Trial, CutFemError> {
        let h = state.phi.h();
        let mut phi = state.phi.clone();
        let vn: Vec<f64> = direction.smooth.iter().map(|w| w - ell).collect();
        if vn.iter().any(|value| !value.is_finite()) {
            return Err(invalid("multi-load shape velocity overflowed"));
        }
        let vmax = vn.iter().fold(0.0f64, |m, value| m.max(value.abs())).max(1e-12);
        let band = build_band(&phi, self.settings.band_cells);
        let duration = self.settings.move_cells * h / vmax;
        if !duration.is_finite() {
            return Err(invalid("multi-load advection duration overflowed"));
        }
        advect(&mut phi, &band, &Velocity::Normal(&vn), duration, 0.45);
        let mut load_pad_nodes = retain_load_pads(&mut phi, &self.supports);
        let mut audit = redistance(&mut phi, self.settings.band_cells);
        load_pad_nodes += retain_load_pads(&mut phi, &self.supports);
        let mut events = Vec::new();
        if let Some(dt_field) = &direction.topological {
            events = nucleate(
                &mut phi,
                dt_field,
                ell,
                self.settings.hole_radius_cells * h,
                6.0 * self.settings.hole_radius_cells * h,
                2,
            );
            if !events.is_empty() {
                load_pad_nodes += retain_load_pads(&mut phi, &self.supports);
                let first_drift = audit.interface_drift_h;
                audit = redistance(&mut phi, self.settings.band_cells);
                audit.interface_drift_h += first_drift;
                load_pad_nodes += retain_load_pads(&mut phi, &self.supports);
            }
        }
        if phi.nodes().iter().any(|value| !value.is_finite()) || !audit.interface_drift_h.is_finite() {
            return Err(invalid("multi-load evolution produced non-finite geometry or audit"));
        }
        Ok(Trial { phi, audit, events, load_pad_nodes })
    }
}

fn append_iteration(
    trajectory: &mut OptimizeReport,
    case_history: &mut Vec<Vec<f64>>,
    active_history: &mut Vec<Option<usize>>,
    state: &MultiState,
    iteration: usize,
    ell: f64,
    audit: RedistanceAudit,
    events: Vec<NucleationEvent>,
    load_pad_nodes: usize,
) {
    let snap = fnv(&state.phi);
    let cases = state.compliances.iter().map(|value| format!("{value:.17e}"))
        .collect::<Vec<_>>().join(",");
    let active_json = state.active.map_or_else(|| "null".to_string(), |index| index.to_string());
    let objective = state.objective;
    let volume = state.volume;
    let mut row = String::new();
    let _ = write!(row,
        "{{\"iter\":{iteration},\"compliance\":{objective:.17e},\"case_compliances\":[{cases}],\"active_case\":{active_json},\"volume\":{volume:.17e},\"ell\":{ell:.17e},\"drift_h\":{:.17e},\"load_pad_nodes\":{load_pad_nodes},\"snapshot\":\"{snap:#018x}\"}}",
        audit.interface_drift_h);
    trajectory.rows.push(row);
    trajectory.compliance.push(objective);
    trajectory.volume.push(volume);
    trajectory.ell.push(ell);
    trajectory.audits.push(audit);
    trajectory.events.extend(events);
    trajectory.snapshots.push(snap);
    trajectory.load_pad_nodes.push(load_pad_nodes);
    case_history.push(state.compliances.clone());
    active_history.push(state.active);
}

/// Run simultaneous independent-load level-set descent with evaluated-state
/// publication. The initial geometry costs one set of per-case solves; every
/// iteration then costs one additional set after the complete geometry update.
/// A refused candidate leaves `phi` at the last successfully evaluated design.
///
/// # Errors
/// Propagates typed material/load/support/settings and canonical CutFEM refusals.
pub fn optimize_compliance_multi_load(
    phi: &mut GridSdf,
    load_cases: &[RobustLoadCase],
    settings: OptimizeSettings,
    aggregate: RobustAggregate,
) -> Result<RobustDescentReport, CutFemError> {
    validate(phi, load_cases, settings)?;
    material(settings)?;
    let mut trajectory = OptimizeReport::default();
    let mut case_history = Vec::with_capacity(settings.iterations);
    let mut active_history = Vec::with_capacity(settings.iterations);
    if settings.iterations == 0 {
        return Ok(RobustDescentReport {
            trajectory, case_compliances: case_history, active_case: active_history, aggregate,
        });
    }
    let kernel = Kernel::new(phi, load_cases, settings, aggregate)?;
    let mut current = kernel.evaluate(phi.clone())?;
    let mut ell = settings.ell0;
    for iteration in 0..settings.iterations {
        let direction = kernel.direction(&current, iteration)?;
        let Trial { phi: trial_phi, audit, events, load_pad_nodes } =
            kernel.propose(&current, &direction, ell)?;
        let candidate = kernel.evaluate(trial_phi)?;
        let next_ell = ell + settings.mu_al * direction.mean_energy.abs().max(1e-30)
            * (candidate.volume - settings.volfrac) / settings.volfrac;
        if !next_ell.is_finite() {
            return Err(invalid("multi-load volume multiplier update overflowed"));
        }
        ell = next_ell.max(0.0);
        append_iteration(
            &mut trajectory,
            &mut case_history,
            &mut active_history,
            &candidate,
            iteration,
            ell,
            audit,
            events,
            load_pad_nodes,
        );
        *phi = candidate.phi.clone();
        current = candidate;
    }
    Ok(RobustDescentReport {
        trajectory,
        case_compliances: case_history,
        active_case: active_history,
        aggregate,
    })
}
