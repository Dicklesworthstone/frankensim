//! Deterministic Fermat connector for two planar dielectric interfaces.
//!
//! This is the geometry kernel used by the tracer's finite smooth-dielectric
//! next-event proposal.  It solves for the stationary optical path from an
//! exterior source, through two prescribed interface planes, to an exterior
//! target.  Scene admission and visibility remain the caller's responsibility:
//! a solution on the infinite support planes is not by itself evidence that it
//! lies on the finite faces of a closed solid.

use fs_geom::{Point3, Vec3};

const MAX_NEWTON_STEPS: usize = 32;
const GRADIENT_TOLERANCE: f64 = 2.0e-11;
const STEP_TOLERANCE: f64 = 2.0e-12;
const MIN_SEGMENT_M: f64 = 1.0e-10;
const MIN_PIVOT: f64 = 1.0e-14;

#[derive(Clone, Copy, Debug)]
pub(super) struct PlaneFrame {
    pub(super) origin: Point3,
    pub(super) tangent: Vec3,
    pub(super) bitangent: Vec3,
    /// Unit normal pointing out of the dielectric.
    pub(super) outward_normal: Vec3,
}

impl PlaneFrame {
    pub(super) fn point(self, coordinates: [f64; 2]) -> Point3 {
        self.origin
            .offset(self.tangent.scale(coordinates[0]))
            .offset(self.bitangent.scale(coordinates[1]))
    }

    pub(super) fn coordinates(self, point: Point3) -> [f64; 2] {
        let displacement = point.delta_from(self.origin);
        [
            displacement.dot(self.tangent),
            displacement.dot(self.bitangent),
        ]
    }
}

#[derive(Clone, Copy, Debug)]
pub(super) struct PlanarDielectricConnection {
    pub(super) entry_point: Point3,
    pub(super) exit_point: Point3,
    pub(super) incident_direction: Vec3,
    pub(super) internal_direction: Vec3,
    pub(super) outgoing_direction: Vec3,
    pub(super) optical_path_length_m: f64,
}

#[derive(Clone, Copy, Debug)]
pub(super) struct PlanarDielectricReflectionConnection {
    pub(super) entry_point: Point3,
    pub(super) reflection_point: Point3,
    pub(super) exit_point: Point3,
    pub(super) incident_direction: Vec3,
    pub(super) internal_incident_direction: Vec3,
    pub(super) internal_reflected_direction: Vec3,
    pub(super) outgoing_direction: Vec3,
    pub(super) optical_path_length_m: f64,
}

/// Exact target-area to source-solid-angle Jacobian for a stationary
/// transmission/reflection/transmission connection.
pub(super) fn reflected_target_area_per_source_solid_angle(
    source: Point3,
    target: Point3,
    entry: PlaneFrame,
    reflection: PlaneFrame,
    exit: PlaneFrame,
    eta_glass: f64,
    connection: PlanarDielectricReflectionConnection,
    target_tangent: Vec3,
    target_bitangent: Vec3,
) -> Option<f64> {
    if !finite_vec(target_tangent)
        || !finite_vec(target_bitangent)
        || (target_tangent.norm() - 1.0).abs() > 2.0e-10
        || (target_bitangent.norm() - 1.0).abs() > 2.0e-10
        || target_tangent.dot(target_bitangent).abs() > 2.0e-10
    {
        return None;
    }
    let p = entry.coordinates(connection.entry_point);
    let r = reflection.coordinates(connection.reflection_point);
    let q = exit.coordinates(connection.exit_point);
    let parameters = [p[0], p[1], r[0], r[1], q[0], q[1]];
    let system = reflected_fermat_system(
        source, target, entry, reflection, exit, eta_glass, parameters,
    )?;
    let (_, source_direction_hessian) =
        direction_hessian(connection.entry_point.delta_from(source))?;
    let (_, target_direction_hessian) =
        direction_hessian(target.delta_from(connection.exit_point))?;
    let exit_basis = [exit.tangent, exit.bitangent];
    let target_basis = [target_tangent, target_bitangent];
    let mut direction_derivatives = [Vec3::new(0.0, 0.0, 0.0); 2];
    for target_axis in 0..2 {
        let mut rhs = [0.0; 6];
        for row in 0..2 {
            rhs[row + 4] = bilinear(
                exit_basis[row],
                target_direction_hessian,
                target_basis[target_axis],
            );
        }
        let derivative = solve_square(system.hessian, rhs)?;
        let entry_derivative = add(
            entry.tangent.scale(derivative[0]),
            entry.bitangent.scale(derivative[1]),
        );
        direction_derivatives[target_axis] =
            matrix_vector(source_direction_hessian, entry_derivative);
    }
    let solid_angle_per_area = connection
        .incident_direction
        .dot(cross(direction_derivatives[0], direction_derivatives[1]))
        .abs();
    (solid_angle_per_area.is_finite() && solid_angle_per_area > 0.0)
        .then(|| 1.0 / solid_angle_per_area)
}

/// Find the stationary path that transmits into a homogeneous dielectric,
/// reflects once from a prescribed internal face, and transmits back out.
pub(super) fn solve_three_planar_interfaces_one_reflection(
    source: Point3,
    target: Point3,
    entry: PlaneFrame,
    reflection: PlaneFrame,
    exit: PlaneFrame,
    eta_glass: f64,
    entry_seed: Point3,
    _reflection_seed: Point3,
    exit_seed: Point3,
) -> Option<PlanarDielectricReflectionConnection> {
    if !finite_point(source)
        || !finite_point(target)
        || !valid_frame(entry)
        || !valid_frame(reflection)
        || !valid_frame(exit)
        || !eta_glass.is_finite()
        || eta_glass <= 1.0
        || source.delta_from(entry.origin).dot(entry.outward_normal) <= 0.0
        || target.delta_from(exit.origin).dot(exit.outward_normal) <= 0.0
    {
        return None;
    }
    // Unfold the specular reflection. Reflection is an isometry, so the
    // stationary T-R-T path is exactly the stationary two-interface path to
    // the mirrored target through the mirrored exit plane. This replaces a
    // six-variable Newton solve with the existing four-variable solve without
    // changing the optical path or estimator.
    let unfolded_target = reflect_point_across_plane(target, reflection);
    let unfolded_exit = reflect_frame_across_plane(exit, reflection);
    let unfolded = solve_two_planar_interfaces(
        source,
        unfolded_target,
        entry,
        unfolded_exit,
        eta_glass,
        entry_seed,
        reflect_point_across_plane(exit_seed, reflection),
    )?;

    let unfolded_internal = unfolded.exit_point.delta_from(unfolded.entry_point);
    let denominator = unfolded_internal.dot(reflection.outward_normal);
    if !denominator.is_finite() || denominator.abs() <= MIN_PIVOT {
        return None;
    }
    let reflection_fraction = reflection
        .origin
        .delta_from(unfolded.entry_point)
        .dot(reflection.outward_normal)
        / denominator;
    if !reflection_fraction.is_finite() || !(0.0..1.0).contains(&reflection_fraction) {
        return None;
    }
    let reflection_point = unfolded
        .entry_point
        .offset(unfolded_internal.scale(reflection_fraction));
    let exit_point = reflect_point_across_plane(unfolded.exit_point, reflection);
    let p = entry.coordinates(unfolded.entry_point);
    let r = reflection.coordinates(reflection_point);
    let q = exit.coordinates(exit_point);
    finish_reflected_connection(
        source,
        target,
        entry,
        reflection,
        exit,
        eta_glass,
        [p[0], p[1], r[0], r[1], q[0], q[1]],
    )
}

/// Exact local change of variables from a target-area proposal to source
/// solid angle for one admitted stationary connection.
///
/// The derivative follows from the implicit-function theorem applied to the
/// same Fermat gradient and analytic Hessian used by the solve.  This avoids a
/// finite-difference PDF: an approximate density would bias next-event MIS
/// even when the geometric path itself had converged.
pub(super) fn target_area_per_source_solid_angle(
    source: Point3,
    target: Point3,
    entry: PlaneFrame,
    exit: PlaneFrame,
    eta_glass: f64,
    connection: PlanarDielectricConnection,
    target_tangent: Vec3,
    target_bitangent: Vec3,
) -> Option<f64> {
    if !finite_vec(target_tangent)
        || !finite_vec(target_bitangent)
        || (target_tangent.norm() - 1.0).abs() > 2.0e-10
        || (target_bitangent.norm() - 1.0).abs() > 2.0e-10
        || target_tangent.dot(target_bitangent).abs() > 2.0e-10
    {
        return None;
    }
    let entry_coordinates = entry.coordinates(connection.entry_point);
    let exit_coordinates = exit.coordinates(connection.exit_point);
    let parameters = [
        entry_coordinates[0],
        entry_coordinates[1],
        exit_coordinates[0],
        exit_coordinates[1],
    ];
    let system = fermat_system(source, target, entry, exit, eta_glass, parameters)?;
    let (_, source_direction_hessian) =
        direction_hessian(connection.entry_point.delta_from(source))?;
    let (_, target_direction_hessian) =
        direction_hessian(target.delta_from(connection.exit_point))?;
    let exit_basis = [exit.tangent, exit.bitangent];
    let target_basis = [target_tangent, target_bitangent];
    let mut direction_derivatives = [Vec3::new(0.0, 0.0, 0.0); 2];
    for target_axis in 0..2 {
        // H dx/dy = -dg/dy. Only the exit-interface rows depend on
        // target position and dg_q/dT = -B_exit^T H_target.
        let mut rhs = [0.0; 4];
        for row in 0..2 {
            rhs[row + 2] = bilinear(
                exit_basis[row],
                target_direction_hessian,
                target_basis[target_axis],
            );
        }
        let derivative = solve_4x4(system.hessian, rhs)?;
        let entry_point_derivative = add(
            entry.tangent.scale(derivative[0]),
            entry.bitangent.scale(derivative[1]),
        );
        direction_derivatives[target_axis] =
            matrix_vector(source_direction_hessian, entry_point_derivative);
    }
    let solid_angle_per_area = connection
        .incident_direction
        .dot(cross(direction_derivatives[0], direction_derivatives[1]))
        .abs();
    (solid_angle_per_area.is_finite() && solid_angle_per_area > 0.0)
        .then(|| 1.0 / solid_angle_per_area)
}

/// Find the stationary transmission path through two planar interfaces.
///
/// `entry_seed` and `exit_seed` need only be reasonable points on their
/// respective support planes.  The Newton system is the projected gradient of
/// `|S-P| + eta |P-Q| + |Q-T|`; its analytic Hessian makes the result
/// deterministic and avoids a finite-difference Snell residual.  A damped
/// line search refuses non-convergence rather than returning a plausible but
/// non-stationary path.
pub(super) fn solve_two_planar_interfaces(
    source: Point3,
    target: Point3,
    entry: PlaneFrame,
    exit: PlaneFrame,
    eta_glass: f64,
    entry_seed: Point3,
    exit_seed: Point3,
) -> Option<PlanarDielectricConnection> {
    if !finite_point(source)
        || !finite_point(target)
        || !valid_frame(entry)
        || !valid_frame(exit)
        || !eta_glass.is_finite()
        || eta_glass <= 1.0
    {
        return None;
    }
    let entry_side = source.delta_from(entry.origin).dot(entry.outward_normal);
    let exit_side = target.delta_from(exit.origin).dot(exit.outward_normal);
    if entry_side <= 0.0 || exit_side <= 0.0 {
        return None;
    }

    let mut parameters = [
        entry.coordinates(entry_seed)[0],
        entry.coordinates(entry_seed)[1],
        exit.coordinates(exit_seed)[0],
        exit.coordinates(exit_seed)[1],
    ];
    let mut previous_objective =
        optical_length(source, target, entry, exit, eta_glass, parameters)?;

    for _ in 0..MAX_NEWTON_STEPS {
        let system = fermat_system(source, target, entry, exit, eta_glass, parameters)?;
        let gradient_norm = system
            .gradient
            .into_iter()
            .map(f64::abs)
            .fold(0.0, f64::max);
        if gradient_norm <= GRADIENT_TOLERANCE {
            return finish_connection(source, target, entry, exit, eta_glass, parameters);
        }
        let rhs = system.gradient.map(|value| -value);
        let step = solve_4x4(system.hessian, rhs)?;
        let step_norm = step.into_iter().map(f64::abs).fold(0.0, f64::max);
        if !step_norm.is_finite() {
            return None;
        }

        let directional_derivative = dot4(system.gradient, step);
        if !directional_derivative.is_finite() || directional_derivative >= 0.0 {
            return None;
        }
        let mut scale = 1.0;
        let mut accepted = None;
        for _ in 0..24 {
            let candidate = std::array::from_fn(|index| parameters[index] + scale * step[index]);
            if let Some(objective) =
                optical_length(source, target, entry, exit, eta_glass, candidate)
                && objective <= previous_objective + 1.0e-4 * scale * directional_derivative
            {
                accepted = Some((candidate, objective));
                break;
            }
            scale *= 0.5;
        }
        let (candidate, objective) = accepted?;
        parameters = candidate;
        previous_objective = objective;
        if scale * step_norm <= STEP_TOLERANCE {
            let system = fermat_system(source, target, entry, exit, eta_glass, parameters)?;
            if system
                .gradient
                .into_iter()
                .map(f64::abs)
                .fold(0.0, f64::max)
                <= 8.0 * GRADIENT_TOLERANCE
            {
                return finish_connection(source, target, entry, exit, eta_glass, parameters);
            }
            return None;
        }
    }
    None
}

#[derive(Clone, Copy)]
struct FermatSystem {
    gradient: [f64; 4],
    hessian: [[f64; 4]; 4],
}

#[derive(Clone, Copy)]
struct ReflectedFermatSystem {
    gradient: [f64; 6],
    hessian: [[f64; 6]; 6],
}

fn reflected_fermat_system(
    source: Point3,
    target: Point3,
    entry: PlaneFrame,
    reflection: PlaneFrame,
    exit: PlaneFrame,
    eta_glass: f64,
    parameters: [f64; 6],
) -> Option<ReflectedFermatSystem> {
    let p = entry.point([parameters[0], parameters[1]]);
    let r = reflection.point([parameters[2], parameters[3]]);
    let q = exit.point([parameters[4], parameters[5]]);
    let (u_source, h_source) = direction_hessian(p.delta_from(source))?;
    let (u_first, h_first) = direction_hessian(r.delta_from(p))?;
    let (u_second, h_second) = direction_hessian(q.delta_from(r))?;
    let (u_target, h_target) = direction_hessian(target.delta_from(q))?;
    let bases = [
        [entry.tangent, entry.bitangent],
        [reflection.tangent, reflection.bitangent],
        [exit.tangent, exit.bitangent],
    ];
    let vector_gradient = [
        sub(u_source, u_first.scale(eta_glass)),
        sub(u_first.scale(eta_glass), u_second.scale(eta_glass)),
        sub(u_second.scale(eta_glass), u_target),
    ];
    let diagonal = [
        matrix_add(h_source, matrix_scale(h_first, eta_glass)),
        matrix_scale(matrix_add(h_first, h_second), eta_glass),
        matrix_add(matrix_scale(h_second, eta_glass), h_target),
    ];
    let adjacent = [
        matrix_scale(h_first, -eta_glass),
        matrix_scale(h_second, -eta_glass),
    ];
    let mut gradient = [0.0; 6];
    let mut hessian = [[0.0; 6]; 6];
    for vertex in 0..3 {
        for row in 0..2 {
            let matrix_row = 2 * vertex + row;
            gradient[matrix_row] = bases[vertex][row].dot(vector_gradient[vertex]);
            for column in 0..2 {
                hessian[matrix_row][2 * vertex + column] =
                    bilinear(bases[vertex][row], diagonal[vertex], bases[vertex][column]);
            }
        }
    }
    for edge in 0..2 {
        for row in 0..2 {
            for column in 0..2 {
                let value = bilinear(bases[edge][row], adjacent[edge], bases[edge + 1][column]);
                hessian[2 * edge + row][2 * (edge + 1) + column] = value;
                hessian[2 * (edge + 1) + column][2 * edge + row] = value;
            }
        }
    }
    gradient
        .iter()
        .chain(hessian.iter().flatten())
        .all(|value| value.is_finite())
        .then_some(ReflectedFermatSystem { gradient, hessian })
}

fn finish_reflected_connection(
    source: Point3,
    target: Point3,
    entry: PlaneFrame,
    reflection: PlaneFrame,
    exit: PlaneFrame,
    eta_glass: f64,
    parameters: [f64; 6],
) -> Option<PlanarDielectricReflectionConnection> {
    let entry_point = entry.point([parameters[0], parameters[1]]);
    let reflection_point = reflection.point([parameters[2], parameters[3]]);
    let exit_point = exit.point([parameters[4], parameters[5]]);
    let incident_direction = unit(entry_point.delta_from(source))?;
    let internal_incident_direction = unit(reflection_point.delta_from(entry_point))?;
    let internal_reflected_direction = unit(exit_point.delta_from(reflection_point))?;
    let outgoing_direction = unit(target.delta_from(exit_point))?;
    if incident_direction.dot(entry.outward_normal) >= 0.0
        || internal_incident_direction.dot(entry.outward_normal) >= 0.0
        || internal_incident_direction.dot(reflection.outward_normal) <= 0.0
        || internal_reflected_direction.dot(reflection.outward_normal) >= 0.0
        || internal_reflected_direction.dot(exit.outward_normal) <= 0.0
        || outgoing_direction.dot(exit.outward_normal) <= 0.0
    {
        return None;
    }
    let system = reflected_fermat_system(
        source, target, entry, reflection, exit, eta_glass, parameters,
    )?;
    if system
        .gradient
        .into_iter()
        .map(f64::abs)
        .fold(0.0, f64::max)
        > 8.0 * GRADIENT_TOLERANCE
    {
        return None;
    }
    Some(PlanarDielectricReflectionConnection {
        entry_point,
        reflection_point,
        exit_point,
        incident_direction,
        internal_incident_direction,
        internal_reflected_direction,
        outgoing_direction,
        optical_path_length_m: reflected_optical_length(
            source, target, entry, reflection, exit, eta_glass, parameters,
        )?,
    })
}

fn reflected_optical_length(
    source: Point3,
    target: Point3,
    entry: PlaneFrame,
    reflection: PlaneFrame,
    exit: PlaneFrame,
    eta_glass: f64,
    parameters: [f64; 6],
) -> Option<f64> {
    let p = entry.point([parameters[0], parameters[1]]);
    let r = reflection.point([parameters[2], parameters[3]]);
    let q = exit.point([parameters[4], parameters[5]]);
    let first = p.delta_from(source).norm();
    let internal_first = r.delta_from(p).norm();
    let internal_second = q.delta_from(r).norm();
    let last = target.delta_from(q).norm();
    let value = first + eta_glass * (internal_first + internal_second) + last;
    (first > MIN_SEGMENT_M
        && internal_first > MIN_SEGMENT_M
        && internal_second > MIN_SEGMENT_M
        && last > MIN_SEGMENT_M
        && value.is_finite())
    .then_some(value)
}

fn fermat_system(
    source: Point3,
    target: Point3,
    entry: PlaneFrame,
    exit: PlaneFrame,
    eta_glass: f64,
    parameters: [f64; 4],
) -> Option<FermatSystem> {
    let p = entry.point([parameters[0], parameters[1]]);
    let q = exit.point([parameters[2], parameters[3]]);
    let (u_source, h_source) = direction_hessian(p.delta_from(source))?;
    let (u_internal, h_internal) = direction_hessian(q.delta_from(p))?;
    let (u_target, h_target) = direction_hessian(target.delta_from(q))?;
    let entry_basis = [entry.tangent, entry.bitangent];
    let exit_basis = [exit.tangent, exit.bitangent];
    let gradient_p = sub(u_source, u_internal.scale(eta_glass));
    let gradient_q = sub(u_internal.scale(eta_glass), u_target);
    let mut gradient = [0.0; 4];
    for row in 0..2 {
        gradient[row] = entry_basis[row].dot(gradient_p);
        gradient[row + 2] = exit_basis[row].dot(gradient_q);
    }

    let h_pp = matrix_add(h_source, matrix_scale(h_internal, eta_glass));
    let h_pq = matrix_scale(h_internal, -eta_glass);
    let h_qq = matrix_add(matrix_scale(h_internal, eta_glass), h_target);
    let mut hessian = [[0.0; 4]; 4];
    for row in 0..2 {
        for column in 0..2 {
            hessian[row][column] = bilinear(entry_basis[row], h_pp, entry_basis[column]);
            hessian[row][column + 2] = bilinear(entry_basis[row], h_pq, exit_basis[column]);
            hessian[row + 2][column] = bilinear(exit_basis[row], h_pq, entry_basis[column]);
            hessian[row + 2][column + 2] = bilinear(exit_basis[row], h_qq, exit_basis[column]);
        }
    }
    gradient
        .iter()
        .chain(hessian.iter().flatten())
        .all(|value| value.is_finite())
        .then_some(FermatSystem { gradient, hessian })
}

fn finish_connection(
    source: Point3,
    target: Point3,
    entry: PlaneFrame,
    exit: PlaneFrame,
    eta_glass: f64,
    parameters: [f64; 4],
) -> Option<PlanarDielectricConnection> {
    let entry_point = entry.point([parameters[0], parameters[1]]);
    let exit_point = exit.point([parameters[2], parameters[3]]);
    let incident_direction = unit(entry_point.delta_from(source))?;
    let internal_direction = unit(exit_point.delta_from(entry_point))?;
    let outgoing_direction = unit(target.delta_from(exit_point))?;
    // These oriented-side checks reject stationary extensions that cross a
    // support plane in the wrong direction.
    if incident_direction.dot(entry.outward_normal) >= 0.0
        || internal_direction.dot(entry.outward_normal) >= 0.0
        || internal_direction.dot(exit.outward_normal) <= 0.0
        || outgoing_direction.dot(exit.outward_normal) <= 0.0
    {
        return None;
    }
    let system = fermat_system(source, target, entry, exit, eta_glass, parameters)?;
    if system
        .gradient
        .into_iter()
        .map(f64::abs)
        .fold(0.0, f64::max)
        > 8.0 * GRADIENT_TOLERANCE
    {
        return None;
    }
    Some(PlanarDielectricConnection {
        entry_point,
        exit_point,
        incident_direction,
        internal_direction,
        outgoing_direction,
        optical_path_length_m: optical_length(source, target, entry, exit, eta_glass, parameters)?,
    })
}

fn optical_length(
    source: Point3,
    target: Point3,
    entry: PlaneFrame,
    exit: PlaneFrame,
    eta_glass: f64,
    parameters: [f64; 4],
) -> Option<f64> {
    let p = entry.point([parameters[0], parameters[1]]);
    let q = exit.point([parameters[2], parameters[3]]);
    let first = p.delta_from(source).norm();
    let internal = q.delta_from(p).norm();
    let last = target.delta_from(q).norm();
    let value = first + eta_glass * internal + last;
    (first > MIN_SEGMENT_M && internal > MIN_SEGMENT_M && last > MIN_SEGMENT_M && value.is_finite())
        .then_some(value)
}

fn direction_hessian(displacement: Vec3) -> Option<(Vec3, [[f64; 3]; 3])> {
    let length = displacement.norm();
    if !length.is_finite() || length <= MIN_SEGMENT_M {
        return None;
    }
    let direction = displacement.scale(1.0 / length);
    let values = [direction.x, direction.y, direction.z];
    let hessian = std::array::from_fn(|row| {
        std::array::from_fn(|column| {
            (if row == column { 1.0 } else { 0.0 } - values[row] * values[column]) / length
        })
    });
    Some((direction, hessian))
}

fn solve_4x4(mut matrix: [[f64; 4]; 4], mut rhs: [f64; 4]) -> Option<[f64; 4]> {
    for column in 0..4 {
        let pivot = (column..4).max_by(|left, right| {
            matrix[*left][column]
                .abs()
                .total_cmp(&matrix[*right][column].abs())
        })?;
        if !matrix[pivot][column].is_finite() || matrix[pivot][column].abs() <= MIN_PIVOT {
            return None;
        }
        matrix.swap(column, pivot);
        rhs.swap(column, pivot);
        for row in column + 1..4 {
            let factor = matrix[row][column] / matrix[column][column];
            for inner in column..4 {
                matrix[row][inner] -= factor * matrix[column][inner];
            }
            rhs[row] -= factor * rhs[column];
        }
    }
    let mut solution = [0.0; 4];
    for row in (0..4).rev() {
        let tail = (row + 1..4)
            .map(|column| matrix[row][column] * solution[column])
            .sum::<f64>();
        solution[row] = (rhs[row] - tail) / matrix[row][row];
    }
    solution
        .iter()
        .all(|value| value.is_finite())
        .then_some(solution)
}

fn solve_square<const N: usize>(mut matrix: [[f64; N]; N], mut rhs: [f64; N]) -> Option<[f64; N]> {
    for column in 0..N {
        let pivot = (column..N).max_by(|left, right| {
            matrix[*left][column]
                .abs()
                .total_cmp(&matrix[*right][column].abs())
        })?;
        if !matrix[pivot][column].is_finite() || matrix[pivot][column].abs() <= MIN_PIVOT {
            return None;
        }
        matrix.swap(column, pivot);
        rhs.swap(column, pivot);
        for row in column + 1..N {
            let factor = matrix[row][column] / matrix[column][column];
            for inner in column..N {
                matrix[row][inner] -= factor * matrix[column][inner];
            }
            rhs[row] -= factor * rhs[column];
        }
    }
    let mut solution = [0.0; N];
    for row in (0..N).rev() {
        let tail = (row + 1..N)
            .map(|column| matrix[row][column] * solution[column])
            .sum::<f64>();
        solution[row] = (rhs[row] - tail) / matrix[row][row];
    }
    solution
        .iter()
        .all(|value| value.is_finite())
        .then_some(solution)
}

fn valid_frame(frame: PlaneFrame) -> bool {
    if !finite_point(frame.origin)
        || !finite_vec(frame.tangent)
        || !finite_vec(frame.bitangent)
        || !finite_vec(frame.outward_normal)
    {
        return false;
    }
    let tolerance = 2.0e-10;
    (frame.tangent.norm() - 1.0).abs() <= tolerance
        && (frame.bitangent.norm() - 1.0).abs() <= tolerance
        && (frame.outward_normal.norm() - 1.0).abs() <= tolerance
        && frame.tangent.dot(frame.bitangent).abs() <= tolerance
        && frame.tangent.dot(frame.outward_normal).abs() <= tolerance
        && frame.bitangent.dot(frame.outward_normal).abs() <= tolerance
}

fn unit(value: Vec3) -> Option<Vec3> {
    let norm = value.norm();
    (norm.is_finite() && norm > MIN_SEGMENT_M).then(|| value.scale(1.0 / norm))
}

fn finite_point(value: Point3) -> bool {
    value.x.is_finite() && value.y.is_finite() && value.z.is_finite()
}

fn finite_vec(value: Vec3) -> bool {
    value.x.is_finite() && value.y.is_finite() && value.z.is_finite()
}

fn sub(left: Vec3, right: Vec3) -> Vec3 {
    Vec3::new(left.x - right.x, left.y - right.y, left.z - right.z)
}

fn add(left: Vec3, right: Vec3) -> Vec3 {
    Vec3::new(left.x + right.x, left.y + right.y, left.z + right.z)
}

fn matrix_add(left: [[f64; 3]; 3], right: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    std::array::from_fn(|row| std::array::from_fn(|column| left[row][column] + right[row][column]))
}

fn matrix_scale(matrix: [[f64; 3]; 3], scale: f64) -> [[f64; 3]; 3] {
    matrix.map(|row| row.map(|value| value * scale))
}

fn bilinear(left: Vec3, matrix: [[f64; 3]; 3], right: Vec3) -> f64 {
    let right = [right.x, right.y, right.z];
    let product = matrix.map(|row| {
        row.into_iter()
            .zip(right)
            .map(|(value, coordinate)| value * coordinate)
            .sum::<f64>()
    });
    [left.x, left.y, left.z]
        .into_iter()
        .zip(product)
        .map(|(coordinate, value)| coordinate * value)
        .sum()
}

fn matrix_vector(matrix: [[f64; 3]; 3], right: Vec3) -> Vec3 {
    let right = [right.x, right.y, right.z];
    let product = matrix.map(|row| {
        row.into_iter()
            .zip(right)
            .map(|(value, coordinate)| value * coordinate)
            .sum::<f64>()
    });
    Vec3::new(product[0], product[1], product[2])
}

fn reflect_point_across_plane(point: Point3, plane: PlaneFrame) -> Point3 {
    let signed_distance = point.delta_from(plane.origin).dot(plane.outward_normal);
    point.offset(plane.outward_normal.scale(-2.0 * signed_distance))
}

fn reflect_vector_across_plane(vector: Vec3, plane: PlaneFrame) -> Vec3 {
    add(
        vector,
        plane
            .outward_normal
            .scale(-2.0 * vector.dot(plane.outward_normal)),
    )
}

fn reflect_frame_across_plane(frame: PlaneFrame, mirror: PlaneFrame) -> PlaneFrame {
    PlaneFrame {
        origin: reflect_point_across_plane(frame.origin, mirror),
        tangent: reflect_vector_across_plane(frame.tangent, mirror),
        bitangent: reflect_vector_across_plane(frame.bitangent, mirror),
        outward_normal: reflect_vector_across_plane(frame.outward_normal, mirror),
    }
}

fn cross(left: Vec3, right: Vec3) -> Vec3 {
    Vec3::new(
        left.y * right.z - left.z * right.y,
        left.z * right.x - left.x * right.z,
        left.x * right.y - left.y * right.x,
    )
}

fn dot4(left: [f64; 4], right: [f64; 4]) -> f64 {
    left.into_iter().zip(right).map(|(a, b)| a * b).sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn frame(origin: Point3, outward_normal: Vec3) -> PlaneFrame {
        let reference = if outward_normal.x.abs() <= outward_normal.y.abs()
            && outward_normal.x.abs() <= outward_normal.z.abs()
        {
            Vec3::new(1.0, 0.0, 0.0)
        } else if outward_normal.y.abs() <= outward_normal.z.abs() {
            Vec3::new(0.0, 1.0, 0.0)
        } else {
            Vec3::new(0.0, 0.0, 1.0)
        };
        let raw_tangent = super::cross(reference, outward_normal);
        let tangent = raw_tangent.scale(1.0 / raw_tangent.norm());
        let bitangent = super::cross(outward_normal, tangent);
        PlaneFrame {
            origin,
            tangent,
            bitangent,
            outward_normal,
        }
    }

    fn assert_near(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "actual={actual:.17e} expected={expected:.17e} tolerance={tolerance:.3e}"
        );
    }

    #[test]
    fn normal_parallel_connection_is_exactly_axial() {
        let source = Point3::new(0.0, 0.0, -1.0);
        let target = Point3::new(0.0, 0.0, 2.0);
        let entry = frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, -1.0));
        let exit = frame(Point3::new(0.0, 0.0, 0.2), Vec3::new(0.0, 0.0, 1.0));
        let connection = solve_two_planar_interfaces(
            source,
            target,
            entry,
            exit,
            1.5,
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 0.2),
        )
        .unwrap();
        assert_eq!(connection.incident_direction, Vec3::new(0.0, 0.0, 1.0));
        assert_eq!(connection.internal_direction, Vec3::new(0.0, 0.0, 1.0));
        assert_eq!(connection.outgoing_direction, Vec3::new(0.0, 0.0, 1.0));
        assert_near(connection.optical_path_length_m, 3.1, 8.0 * f64::EPSILON);
    }

    #[test]
    fn parallel_off_axis_connection_obeys_snell_and_reciprocity() {
        let source = Point3::new(-0.3, 0.1, -0.8);
        let target = Point3::new(0.9, -0.2, 1.4);
        let entry = frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, -1.0));
        let exit = frame(Point3::new(0.0, 0.0, 0.25), Vec3::new(0.0, 0.0, 1.0));
        let forward = solve_two_planar_interfaces(
            source,
            target,
            entry,
            exit,
            1.52,
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.1, 0.0, 0.25),
        )
        .unwrap();
        let sin_external = super::cross(forward.incident_direction, entry.outward_normal).norm();
        let sin_internal = super::cross(forward.internal_direction, entry.outward_normal).norm();
        assert_near(sin_external, 1.52 * sin_internal, 2.0e-10);

        let reverse = solve_two_planar_interfaces(
            target,
            source,
            PlaneFrame {
                origin: exit.origin,
                tangent: exit.tangent,
                bitangent: exit.bitangent,
                outward_normal: exit.outward_normal,
            },
            PlaneFrame {
                origin: entry.origin,
                tangent: entry.tangent,
                bitangent: entry.bitangent,
                outward_normal: entry.outward_normal,
            },
            1.52,
            forward.exit_point,
            forward.entry_point,
        )
        .unwrap();
        assert!(reverse.entry_point.delta_from(forward.exit_point).norm() <= 2.0e-10);
        assert!(reverse.exit_point.delta_from(forward.entry_point).norm() <= 2.0e-10);
    }

    #[test]
    fn nonparallel_exit_plane_has_stationary_fermat_path() {
        let source = Point3::new(-0.4, 0.1, -0.9);
        let target = Point3::new(0.8, -0.15, 1.2);
        let entry = frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, -1.0));
        let exit_normal = {
            let raw = Vec3::new(0.25, 0.0, 1.0);
            raw.scale(1.0 / raw.norm())
        };
        let exit = frame(Point3::new(0.0, 0.0, 0.25), exit_normal);
        let connection = solve_two_planar_interfaces(
            source,
            target,
            entry,
            exit,
            1.5,
            Point3::new(-0.05, 0.0, 0.0),
            Point3::new(0.1, 0.0, 0.225),
        )
        .unwrap();
        let entry_external_sine =
            super::cross(connection.incident_direction, entry.outward_normal).norm();
        let entry_internal_sine =
            super::cross(connection.internal_direction, entry.outward_normal).norm();
        let exit_internal_sine =
            super::cross(connection.internal_direction, exit.outward_normal).norm();
        let exit_external_sine =
            super::cross(connection.outgoing_direction, exit.outward_normal).norm();
        assert_near(entry_external_sine, 1.5 * entry_internal_sine, 3.0e-10);
        assert_near(1.5 * exit_internal_sine, exit_external_sine, 3.0e-10);
    }

    #[test]
    fn implicit_area_jacobian_matches_central_difference() {
        let source = Point3::new(-0.3, 0.1, -0.8);
        let target = Point3::new(0.9, -0.2, 1.4);
        let entry = frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, -1.0));
        let exit = frame(Point3::new(0.0, 0.0, 0.25), Vec3::new(0.0, 0.0, 1.0));
        let connection = solve_two_planar_interfaces(
            source,
            target,
            entry,
            exit,
            1.52,
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.1, 0.0, 0.25),
        )
        .unwrap();
        let area_per_solid_angle = target_area_per_source_solid_angle(
            source,
            target,
            entry,
            exit,
            1.52,
            connection,
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        )
        .unwrap();
        let epsilon = 1.0e-5;
        let shifted = |axis: Vec3, sign: f64| {
            let shifted_target = target.offset(axis.scale(sign * epsilon));
            solve_two_planar_interfaces(
                source,
                shifted_target,
                entry,
                exit,
                1.52,
                connection.entry_point,
                connection.exit_point,
            )
            .unwrap()
            .incident_direction
        };
        let derivative_x = sub(
            shifted(Vec3::new(1.0, 0.0, 0.0), 1.0),
            shifted(Vec3::new(1.0, 0.0, 0.0), -1.0),
        )
        .scale(0.5 / epsilon);
        let derivative_y = sub(
            shifted(Vec3::new(0.0, 1.0, 0.0), 1.0),
            shifted(Vec3::new(0.0, 1.0, 0.0), -1.0),
        )
        .scale(0.5 / epsilon);
        let finite_difference = 1.0
            / connection
                .incident_direction
                .dot(super::cross(derivative_x, derivative_y))
                .abs();
        assert_near(
            area_per_solid_angle,
            finite_difference,
            2.0e-7 * finite_difference,
        );
    }

    #[test]
    fn one_internal_reflection_is_reciprocal_and_has_exact_area_jacobian() {
        let source = Point3::new(-0.35, 0.08, -0.9);
        let target = Point3::new(0.55, -0.12, -0.7);
        let lower = frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, -1.0));
        let upper = frame(Point3::new(0.0, 0.0, 0.25), Vec3::new(0.0, 0.0, 1.0));
        let connection = solve_three_planar_interfaces_one_reflection(
            source,
            target,
            lower,
            upper,
            lower,
            1.52,
            Point3::new(-0.2, 0.04, 0.0),
            Point3::new(0.0, 0.0, 0.25),
            Point3::new(0.25, -0.05, 0.0),
        )
        .unwrap();
        let reverse = solve_three_planar_interfaces_one_reflection(
            target,
            source,
            lower,
            upper,
            lower,
            1.52,
            connection.exit_point,
            connection.reflection_point,
            connection.entry_point,
        )
        .unwrap();
        assert!(reverse.entry_point.delta_from(connection.exit_point).norm() <= 2.0e-10);
        assert!(
            reverse
                .reflection_point
                .delta_from(connection.reflection_point)
                .norm()
                <= 2.0e-10
        );
        assert!(reverse.exit_point.delta_from(connection.entry_point).norm() <= 2.0e-10);

        let target_tangent = Vec3::new(1.0, 0.0, 0.0);
        let target_bitangent = Vec3::new(0.0, 1.0, 0.0);
        let exact = reflected_target_area_per_source_solid_angle(
            source,
            target,
            lower,
            upper,
            lower,
            1.52,
            connection,
            target_tangent,
            target_bitangent,
        )
        .unwrap();
        let epsilon = 1.0e-5;
        let shifted = |axis: Vec3, sign: f64| {
            solve_three_planar_interfaces_one_reflection(
                source,
                target.offset(axis.scale(sign * epsilon)),
                lower,
                upper,
                lower,
                1.52,
                connection.entry_point,
                connection.reflection_point,
                connection.exit_point,
            )
            .unwrap()
            .incident_direction
        };
        let derivative_x =
            sub(shifted(target_tangent, 1.0), shifted(target_tangent, -1.0)).scale(0.5 / epsilon);
        let derivative_y = sub(
            shifted(target_bitangent, 1.0),
            shifted(target_bitangent, -1.0),
        )
        .scale(0.5 / epsilon);
        let finite_difference = 1.0
            / connection
                .incident_direction
                .dot(super::cross(derivative_x, derivative_y))
                .abs();
        assert_near(exact, finite_difference, 4.0e-7 * finite_difference);
    }
}
