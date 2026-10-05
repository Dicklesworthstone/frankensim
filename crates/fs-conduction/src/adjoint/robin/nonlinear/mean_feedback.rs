//! Area-mean-dependent Robin laws on a checked, retained FEM field.
//!
//! For h_j(mean_j(T)) and r_j(mean_j(T)), the state correction is u_j v_j^T,
//! u_j = (dh_j/dmean_j) M_j (T - r_j) - h_j (dr_j/dmean_j) M_j 1.
//! v_j contains the area-mean weights. The factors are NOT generally
//! parallel. Apply the true transpose; do not treat this operator as SPD.

use fs_solver::{FgmresState, LinearOp, norm2};
use fs_sparse::Csr;

use super::super::{
    ConductionError, ConductionProblem, Cx, DofMap, LinearConfig, RobinGradient,
    RobinResponse, ThermalInterfaces, add, admit_linear,
    assemble_operator_scaled_with_interfaces, bind_ports, checked, failed,
    invalid, poll, reduce, true_residual, vector,
};

struct Update { left: Vec<f64>, right: Vec<f64> }

struct FeedbackOp<'a> {
    forward: &'a Csr,
    reverse: &'a Csr,
    updates: &'a [Update],
    transposed: bool,
}
impl FeedbackOp<'_> {
    fn oriented_apply(&self, x: &[f64], y: &mut [f64], transposed: bool) {
        if transposed { self.reverse.spmv(x, y); }
        else { self.forward.spmv(x, y); }
        for update in self.updates {
            let (left, right) = if transposed { (&update.right, &update.left) }
                else { (&update.left, &update.right) };
            let factor: f64 = right.iter().zip(x).map(|(a, b)| a * b).sum();
            for (value, a) in y.iter_mut().zip(left) { *value += a * factor; }
        }
    }
}
impl LinearOp for FeedbackOp<'_> {
    fn n(&self) -> usize { self.forward.nrows() }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        self.oriented_apply(x, y, self.transposed);
    }
    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        self.oriented_apply(x, y, !self.transposed);
    }
}

fn zeros(n: usize) -> Result<Vec<f64>, ConductionError> {
    let mut values = Vec::new();
    values.try_reserve_exact(n).map_err(|_| invalid("mean-HTC feedback allocation refused"))?;
    values.resize(n, 0.0);
    Ok(values)
}

fn feedback_residual(op: &impl LinearOp, x: &[f64], rhs: &[f64])
    -> Result<f64, ConductionError>
{
    let mut ax = zeros(rhs.len())?;
    op.apply(x, &mut ax);
    let mut scale = 0.0_f64;
    for &v in ax.iter().chain(rhs) { scale = scale.max(checked(v)?.abs()); }
    if scale == 0.0 { return Ok(0.0); }
    let residual: Vec<_> = rhs.iter().zip(&ax).map(|(b, a)| b / scale - a / scale).collect();
    let normalized: Vec<_> = rhs.iter().map(|b| b / scale).collect();
    checked(norm2(&residual) / norm2(&normalized).max(f64::MIN_POSITIVE))
}

impl RobinResponse {
    /// Differentiate a nodal-temperature goal through area-mean-dependent h.
    ///
    /// References are state-independent in this convenience wrapper. For a
    /// combined convection/radiation reference that depends on mean wall T,
    /// use [`Self::pullback_mean_robin`] with BOTH sets of slopes instead.
    /// The same residual, material, contact, work and cancellation gates apply.
    #[allow(clippy::too_many_arguments)]
    pub fn pullback_mean_htc(
        cx: &Cx<'_>, problem: ConductionProblem<'_>, interfaces: Option<&ThermalInterfaces>,
        linear: LinearConfig, temperature: &[f64], regions: &[&str],
        slopes_w_m2_k2: &[f64], nodal_weights: &[f64], max_feedback_entries: usize,
    ) -> Result<RobinGradient, ConductionError> {
        poll(cx, 0)?;
        if regions.len() > 64 { return Err(invalid("mean-HTC feedback exceeds 64 regions")); }
        Self::pullback_mean_robin(cx, problem, interfaces, linear, temperature, regions,
            slopes_w_m2_k2, &[0.0; 64][..regions.len()], nodal_weights, max_feedback_entries)
    }

    /// Differentiate a nodal goal through h(mean wall T) AND reference(mean wall T).
    ///
    /// The supplied field is checked against the production residual with its
    /// actual material assignment, matching contacts and Dirichlet lift. No
    /// new primal solve is performed. Each selected uniform Robin region gets
    /// dh/dmean in W/(m^2 K^2) and dreference/dmean in K/K. The caller owns
    /// the law and binds the coefficients and BOTH slopes to the same accepted
    /// state. Unselected Robin laws stay constant. Smooth k(T) is differentiated
    /// too; a material kink refuses. The consistent face mass is retained:
    /// differentiating only an integrated patch watt law would be incorrect
    /// on a non-isothermal trace.
    ///
    /// The returned reference and log(h) entries are independent-control
    /// partials through the complete state Jacobian. A physical control p
    /// contributes references * (partial reference/partial p at fixed mean)
    /// plus log_htc * (partial log h/partial p at fixed mean). These are NOT
    /// total state derivatives; feeding them back again double-counts feedback.
    /// Nodal-load derivatives vanish at prescribed nodes.
    ///
    /// This is a discrete local derivative, not an inverse/stability certificate
    /// or a continuum/error/uncertainty bound. `max_feedback_entries` bounds
    /// retained factors (two vectors per region), not total allocator/RSS.
    /// At most 64 regions are admitted. Inner-iteration budgets are exact;
    /// restart is capped at 32, with checkpoints between cycles and in traces.
    #[allow(clippy::too_many_arguments)]
    pub fn pullback_mean_robin(
        cx: &Cx<'_>, problem: ConductionProblem<'_>, interfaces: Option<&ThermalInterfaces>,
        linear: LinearConfig, temperature: &[f64], regions: &[&str],
        slopes_w_m2_k2: &[f64], reference_slopes: &[f64],
        nodal_weights: &[f64], max_feedback_entries: usize,
    ) -> Result<RobinGradient, ConductionError> {
        poll(cx, 0)?;
        admit_linear(linear)?;
        if linear.restart == 0 { return Err(invalid("mean-HTC feedback needs a positive FGMRES restart")); }
        let n = problem.mesh.vertex_count();
        if regions.len() > 64 || n.checked_mul(2).and_then(|v| v.checked_mul(regions.len()))
            .is_none_or(|v| v > max_feedback_entries)
        { return Err(invalid("mean-HTC feedback exceeds the declared vector-entry budget or 64 regions")); }
        vector(cx, temperature, n)?;
        vector(cx, nodal_weights, n)?;
        vector(cx, slopes_w_m2_k2, regions.len())?;
        vector(cx, reference_slopes, regions.len())?;
        let dofs = DofMap::new(problem.boundary, n)?;
        if dofs.fixed().is_empty() && !problem.boundary.has_robin() {
            return Err(ConductionError::SingularPureNeumann);
        }
        for &v in dofs.fixed() {
            poll(cx, v)?;
            if temperature[v] != dofs.prescribed()[v] {
                return Err(invalid("mean-HTC feedback requires the retained prescribed temperatures"));
            }
        }
        let ports = bind_ports(cx, problem, regions)?;
        let system = assemble_operator_scaled_with_interfaces(cx, problem.mesh,
            problem.boundary, problem.material, problem.source, temperature,
            None, interfaces, problem.element_materials)?;
        let (matrix, rhs) = reduce(&system, &dofs);
        let primal_residual = true_residual(&matrix, &dofs.gather(temperature), &rhs)?;
        if primal_residual >= linear.tolerance { return Err(failed(0, primal_residual, linear)); }
        let (jacobian, tangent) = super::prepare(cx, problem, interfaces, temperature, &dofs)?;
        if !tangent.smooth { return Err(invalid("mean-HTC feedback cannot choose a derivative at a material kink")); }
        let mut updates = Vec::with_capacity(ports.len());
        for ((port, &slope), &reference_slope) in ports.iter().zip(slopes_w_m2_k2).zip(reference_slopes) {
            let mut left = zeros(n)?;
            let mut right = zeros(n)?;
            for (vertices, area) in &port.faces {
                poll(cx, 0)?;
                for (a, &v) in vertices.iter().enumerate() {
                    add(&mut right[v], area / port.area_m2 / 3.0)?;
                    // Zero preserves the original mean-HTC arithmetic path.
                    if reference_slope != 0.0 {
                        add(&mut left[v], -port.htc_w_m2_k * reference_slope * (area / 3.0))?;
                    }
                    for (b, &w) in vertices.iter().enumerate() {
                        let mass = (area / 12.0) * if a == b { 2.0 } else { 1.0 };
                        add(&mut left[v], slope * mass * (temperature[w] - port.reference_k))?;
                    }
                }
            }
            updates.push(Update { left: dofs.gather(&left), right: dofs.gather(&right) });
        }
        let op = FeedbackOp { forward: &jacobian, reverse: &tangent.transpose,
            updates: &updates, transposed: true };
        let rhs = dofs.gather(nodal_weights);
        let scale = rhs.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        let mut lambda = zeros(n)?;
        let (relative_residual, iterations) = if scale == 0.0 { (0.0, 0) } else {
            let rhs: Vec<_> = rhs.iter().map(|v| v / scale).collect();
            let pre = crate::solve::spd_preconditioner(&tangent.transpose);
            let restart = linear.restart.min(dofs.n()).min(32);
            let mut state = FgmresState::new(&rhs, restart);
            while state.rel_residual() >= linear.tolerance && state.iters < linear.max_iterations {
                poll(cx, state.iters)?;
                let before = state.iters;
                state.restart = restart.min(linear.max_iterations - before);
                state.run(&op, &pre, &rhs, linear.tolerance, 1);
                if state.iters == before { break; }
            }
            poll(cx, state.iters)?;
            let residual = feedback_residual(&op, &state.x, &rhs)?;
            if residual >= linear.tolerance { return Err(failed(state.iters, residual, linear)); }
            for (i, &v) in dofs.free().iter().enumerate() {
                if i % 512 == 0 { poll(cx, i)?; }
                lambda[v] = checked(state.x[i] * scale)?;
            }
            (residual, state.iters)
        };
        let mut references = zeros(ports.len())?;
        let mut log_htc = zeros(ports.len())?;
        for (i, port) in ports.iter().enumerate() {
            for (vertices, area) in &port.faces {
                poll(cx, i)?;
                for (a, &v) in vertices.iter().enumerate() {
                    add(&mut references[i], lambda[v] * port.htc_w_m2_k * (area / 3.0))?;
                    for (b, &w) in vertices.iter().enumerate() {
                        let mass = port.htc_w_m2_k * (area / 12.0) * if a == b { 2.0 } else { 1.0 };
                        add(&mut log_htc[i], lambda[v] * mass * (port.reference_k - temperature[w]))?;
                    }
                }
            }
        }
        poll(cx, iterations)?;
        Ok(RobinGradient { references, log_htc, nodal_load: lambda, relative_residual, iterations })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_sparse::{Coo, ops};

    #[test]
    fn nonuniform_trace_correction_uses_the_actual_transpose() {
        let mut coo = Coo::new(3, 3);
        for (i, d) in [2.0, 3.0, 4.0].into_iter().enumerate() { coo.push(i, i, d); }
        coo.push(0, 1, 0.7);
        let a = coo.assemble();
        let at = ops::transpose(&a);
        let updates = [Update { left: vec![1.0, 2.0, -0.5], right: vec![0.2, 0.3, 0.5] },
            Update { left: vec![0.0, -0.3, 0.2], right: vec![0.0, 0.5, 0.5] }];
        let op = FeedbackOp { forward: &a, reverse: &at, updates: &updates, transposed: false };
        let x = [0.7, -1.0, 2.0];
        let y = [1.2, 0.4, -0.3];
        let mut ax = [0.0; 3]; let mut aty = [0.0; 3];
        op.apply(&x, &mut ax); op.apply_transpose(&y, &mut aty);
        let left: f64 = ax.iter().zip(y).map(|(a,b)| a*b).sum();
        let right: f64 = aty.iter().zip(x).map(|(a,b)| a*b).sum();
        assert!((left-right).abs() < 1e-14);
        let mut wrong = [0.0; 3]; op.apply(&y, &mut wrong);
        assert!(wrong.iter().zip(aty).any(|(a,b)| (a-b).abs() > 0.1));
        let dual = FeedbackOp { transposed: true, ..op };
        let mut actual = [0.0; 3]; dual.apply(&y, &mut actual);
        assert_eq!(actual, aty);
        assert_eq!(feedback_residual(&dual, &y, &actual).unwrap(), 0.0);
        assert!(feedback_residual(&dual, &[f64::NAN; 3], &actual).is_err());
    }

    use crate::{ConductionMesh, ConductivityModel, ConductivityTable, ScalarField,
        ThermalBoundary, ThermalBoundaryBuilder, ThermalBc, ConductionSolution,
        SolveConfig, InitialGuess};
    use crate::fixtures::{box_grid, on_box_face};
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

    fn with_gate<T>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> T) -> T {
        ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(gate, arena,
            StreamKey { seed: 51, kernel_id: 819, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic)))
    }
    fn config() -> SolveConfig {
        let mut config = SolveConfig::default();
        config.initial = InitialGuess::Uniform(320.0);
        config.linear.tolerance = 1e-10;
        config.stop.residual_rtol = 1e-12;
        config.stop.step_atol = 0.0;
        config
    }

    // Actual three-dimensional FEM fixed point; the spatially varying source
    // makes the wall trace non-isothermal, so replacing u v^T by v u^T fails.
    fn physical(cx: &Cx<'_>, source_shift: f64, ambient: f64)
        -> (ConductionMesh, ConductivityModel, ScalarField, ThermalBoundary, ConductionSolution)
    {
        physical_reference(cx, source_shift, ambient, 0.0)
    }

    fn physical_reference(cx: &Cx<'_>, source_shift: f64, ambient: f64, beta: f64)
        -> (ConductionMesh, ConductivityModel, ScalarField, ThermalBoundary, ConductionSolution)
    {
        let (complex, positions) = box_grid([3, 2, 2], [0.1, 0.04, 0.03]);
        let mesh = ConductionMesh::new(complex, positions).unwrap();
        let material = ConductivityModel::isotropic(ConductivityTable::declared_curve(
            vec![(250.0, 3.0), (450.0, 13.0)]).unwrap());
        let source = ScalarField::Nodal(mesh.positions().iter()
            .map(|p| 1000.0 + source_shift + 50000.0*p[1]).collect());
        let mut h = 40.0;
        let mut reference = ambient;
        for _ in 0..100 {
            let boundary = ThermalBoundaryBuilder::new(&mesh)
                .region("hot", |f| on_box_face(f.centroid[0], 0.0), ThermalBc::dirichlet(340.0).unwrap()).unwrap()
                .region("wall", |f| on_box_face(f.centroid[0], 0.1), ThermalBc::robin(h, reference).unwrap()).unwrap()
                .adiabatic_remainder().finish().unwrap();
            let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
                material: &material, element_materials: None, source: &source };
            let solution = crate::solve(cx, problem, config()).unwrap();
            let mut integral = 0.0; let mut area = 0.0;
            for face in mesh.boundary().iter().filter(|f| on_box_face(f.centroid[0], 0.1)) {
                area += face.area;
                integral += face.area/3.0 * face.vertices.iter()
                    .map(|&v| solution.temperature[v as usize]).sum::<f64>();
            }
            let next = 40.0 + 0.8*(integral/area - ambient - 20.0);
            let next_reference = ambient + beta*(integral/area - ambient);
            assert!(next > 0.0);
            if (next-h).abs() <= 1e-12*h && (next_reference-reference).abs() <= 1e-12*reference {
                return (mesh, material, source, boundary, solution);
            }
            h = 0.5*(h+next);
            if beta != 0.0 { reference = 0.5*(reference+next_reference); }
        }
        panic!("test's physical fixed point did not converge");
    }

    #[test]
    fn mean_feedback_matches_fresh_nonlinear_fem_power_and_ambient_resolves() {
        with_gate(&CancelGate::new_clock_free(), |cx| {
            let (mesh, material, source, boundary, solution) = physical(cx, 0.0, 293.0);
            let n = mesh.vertex_count();
            let mut weights = vec![0.0; n]; weights[n-1] = 1.0;
            let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
                material: &material, element_materials: None, source: &source };
            let gradient = RobinResponse::pullback_mean_htc(cx, problem, None, config().linear,
                &solution.temperature, &["wall"], &[0.8], &weights, 2*n).unwrap();
            let density_bar = |g: &RobinGradient| mesh.complex().tets.iter().enumerate()
                .map(|(e,tet)| mesh.element_volume(e)/4.0 * tet.iter()
                    .map(|&v| g.nodal_load[v as usize]).sum::<f64>()).sum::<f64>();
            let power = density_bar(&gradient);
            let low = physical(cx, -1.0, 293.0).4.temperature[n-1];
            let high = physical(cx, 1.0, 293.0).4.temperature[n-1];
            let expected = (high-low)/2.0;
            assert!((power-expected).abs() < 1e-4*expected.abs().max(1e-5), "{power:e} vs {expected:e}");
            let port = bind_ports(cx, problem, &["wall"]).unwrap();
            let ambient = gradient.references[0] - gradient.log_htc[0]*0.8/port[0].htc_w_m2_k;
            let low = physical(cx, 0.0, 292.99).4.temperature[n-1];
            let high = physical(cx, 0.0, 293.01).4.temperature[n-1];
            let expected = (high-low)/0.02;
            assert!((ambient-expected).abs() < 2e-5, "{ambient:e} vs {expected:e}");
            let frozen = RobinResponse::pullback_mean_htc(cx, problem, None, config().linear,
                &solution.temperature, &["wall"], &[0.0], &weights, 2*n).unwrap();
            assert!((power-density_bar(&frozen)).abs() > 1e-3*power.abs(), "fixture must distinguish frozen h");
            for &v in DofMap::new(&boundary, n).unwrap().fixed() { assert_eq!(gradient.nodal_load[v], 0.0); }
            assert!(gradient.relative_residual < config().linear.tolerance);
            assert!(RobinResponse::pullback_mean_htc(cx, problem, None, config().linear,
                &solution.temperature, &["wall"], &[0.8], &weights, 2*n-1).is_err());
            let zero = RobinResponse::pullback_mean_htc(cx, problem, None, config().linear,
                &solution.temperature, &["wall"], &[0.8], &vec![0.0;n], 2*n).unwrap();
            assert_eq!(zero.iterations, 0);
            assert!(zero.nodal_load.iter().all(|&x| x == 0.0));
            let short = LinearConfig { max_iterations: 1, restart: 1, ..config().linear };
            assert!(matches!(RobinResponse::pullback_mean_htc(cx, problem, None, short,
                &solution.temperature, &["wall"], &[0.8], &weights, 2*n),
                Err(ConductionError::LinearSolveFailed { krylov_iterations: 1, .. })));
            let gate = CancelGate::new_clock_free(); gate.request();
            with_gate(&gate, |cancelled| assert!(matches!(RobinResponse::pullback_mean_htc(
                cancelled, problem, None, config().linear, &solution.temperature,
                &["wall"], &[0.8], &weights, 2*n), Err(ConductionError::Cancelled { .. }))));
        });
    }

    #[test]
    fn mean_reference_feedback_matches_source_and_ambient_fem_resolves() {
        with_gate(&CancelGate::new_clock_free(), |cx| {
            let beta = 0.35;
            let (mesh, material, source, boundary, solution) = physical_reference(cx, 0.0, 293.0, beta);
            let n = mesh.vertex_count();
            let mut weights = vec![0.0; n]; weights[n-1] = 1.0;
            let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
                material: &material, element_materials: None, source: &source };
            let gradient = RobinResponse::pullback_mean_robin(cx, problem, None, config().linear,
                &solution.temperature, &["wall"], &[0.8], &[beta], &weights, 2*n).unwrap();
            let density_bar = |g: &RobinGradient| mesh.complex().tets.iter().enumerate()
                .map(|(e,tet)| mesh.element_volume(e)/4.0 * tet.iter()
                    .map(|&v| g.nodal_load[v as usize]).sum::<f64>()).sum::<f64>();
            let power = density_bar(&gradient);
            let low = physical_reference(cx, -1.0, 293.0, beta).4.temperature[n-1];
            let high = physical_reference(cx, 1.0, 293.0, beta).4.temperature[n-1];
            let expected = (high-low)/2.0;
            assert!((power-expected).abs() < 1e-4*expected.abs().max(1e-5), "{power:e} vs {expected:e}");
            let h = bind_ports(cx, problem, &["wall"]).unwrap()[0].htc_w_m2_k;
            let ambient = gradient.references[0]*(1.0-beta) - gradient.log_htc[0]*0.8/h;
            let low = physical_reference(cx, 0.0, 292.99, beta).4.temperature[n-1];
            let high = physical_reference(cx, 0.0, 293.01, beta).4.temperature[n-1];
            assert!((ambient-(high-low)/0.02).abs() < 2e-5);
            let frozen = RobinResponse::pullback_mean_htc(cx, problem, None, config().linear,
                &solution.temperature, &["wall"], &[0.8], &weights, 2*n).unwrap();
            assert!((power-density_bar(&frozen)).abs() > 1e-3*power.abs(), "fixture must distinguish a frozen reference");
            for &v in DofMap::new(&boundary, n).unwrap().fixed() { assert_eq!(gradient.nodal_load[v], 0.0); }
            for slopes in [&[][..], &[f64::NAN][..]] {
                assert!(RobinResponse::pullback_mean_robin(cx, problem, None, config().linear,
                    &solution.temperature, &["wall"], &[0.8], slopes, &weights, 2*n).is_err());
            }
            let zero_reference = RobinResponse::pullback_mean_robin(cx, problem, None, config().linear,
                &solution.temperature, &["wall"], &[0.8], &[0.0], &weights, 2*n).unwrap();
            assert_eq!(zero_reference, frozen, "compatibility wrapper must be bit-identical");
        });
    }
}
