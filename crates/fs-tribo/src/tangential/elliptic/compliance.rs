//! Direct two-axis body response for the original associated ellipse.
//!
//! In principal traction coordinates, `F = D y` and the sliding equations are
//! `(D C D + lambda I) y = D du_free`, `|y| = 1`. A bounded scalar secular
//! solve replaces the generic matrix adapter's repeated constitutive probes.
//! No compliance is averaged, regularized, or inferred from material names.

use super::*;

/// One original symmetric positive-semidefinite body compliance [m/N].
/// Coordinates are the same tangent frame as the supplied released motion.
/// Admission is floating point: symmetry is exact and a negative normalized
/// determinant refuses. Rank-deficient static forces need not be unique.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SymmetricTangentialCompliance {
    response: SymmetricResponse,
}

impl SymmetricTangentialCompliance {
    pub fn new(matrix_m_per_n: [[f64; 2]; 2]) -> Result<Self, TriboError> {
        Ok(Self {
            response: SymmetricResponse::new(matrix_m_per_n)?,
        })
    }

    #[must_use]
    pub const fn matrix_m_per_n(self) -> [[f64; 2]; 2] {
        self.response.matrix
    }

    pub(super) fn material(self, axis: [f64; 2]) -> Result<([f64; 3], f64), TriboError> {
        self.response.material(axis)
    }
}

// Dimension-neutral symmetric response for the shared projection kernel.
// Public wrappers attach physical units and require PSD compliance or SPD
// impulse mobility without inventing a time interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct SymmetricResponse {
    pub(super) matrix: [[f64; 2]; 2],
    pub(super) scale: f64,
    pub(super) normalized: [f64; 3],
    pub(super) determinant: f64,
}

impl SymmetricResponse {
    pub(super) fn new(m: [[f64; 2]; 2]) -> Result<Self, TriboError> {
        if m.iter().flatten().any(|v| !v.is_finite())
            || m[0][1] != m[1][0]
            || m[0][0] < 0.0
            || m[1][1] < 0.0
        {
            return Err(invalid("elliptic_symmetric_compliance"));
        }
        if m[0][1] != 0.0 && (m[0][0] == 0.0 || m[1][1] == 0.0) {
            return Err(invalid("elliptic_compliance_not_psd"));
        }
        let scale = m[0][0].max(m[1][1]);
        let normalized = if scale == 0.0 {
            if m[0][1] != 0.0 {
                return Err(invalid("elliptic_compliance_not_psd"));
            }
            [0.0; 3]
        } else {
            [m[0][0] / scale, m[0][1] / scale, m[1][1] / scale]
        };
        let [a, b, d] = normalized;
        let determinant = a * d - b * b;
        if !determinant.is_finite() || determinant < 0.0 {
            return Err(invalid("elliptic_compliance_not_psd"));
        }
        if scale > 0.0
            && ((m[0][0] > 0.0 && a == 0.0)
                || (m[1][1] > 0.0 && d == 0.0)
                || (m[0][1] != 0.0 && b == 0.0))
        {
            return Err(invalid("elliptic_compliance_scale"));
        }
        Ok(Self {
            matrix: m,
            scale,
            normalized,
            determinant,
        })
    }

    pub(super) fn material(self, axis: [f64; 2]) -> Result<([f64; 3], f64), TriboError> {
        // Pivot on the unit normalized diagonal. Rotating this Gram factor
        // avoids subtracting nearly equal terms in a rank-one diagonal entry.
        let [a, b, d] = self.normalized;
        let r = det::sqrt(self.determinant);
        let rows = if a >= d {
            [[1.0, b], [0.0, r]]
        } else {
            [[r, 0.0], [b, 1.0]]
        };
        let x = to_principal(axis, rows[0])?;
        let y = to_principal(axis, rows[1])?;
        let matrix = [
            x[0] * x[0] + y[0] * y[0],
            x[0] * x[1] + y[0] * y[1],
            x[1] * x[1] + y[1] * y[1],
        ];
        let length2 = axis[0] * axis[0] + axis[1] * axis[1];
        Ok((matrix, self.determinant * length2 * length2))
    }
}

impl EllipticCoulomb {
    /// Resolve the original associated law against a full body compliance.
    /// `input.compliance_m_per_n` must be zero: the explicit matrix is the sole
    /// compliance owner. Exactly isotropic matrices use the original scalar
    /// path bit for bit. General SPD matrices have a unique static reaction.
    ///
    /// Rank-one compliance admits unique sliding and reactions determined by
    /// a one-axis friction law. Nonzero feasible static demands with two active
    /// friction axes conservatively refuse, including a boundary solution;
    /// the method does not choose among possible null-space forces.
    /// Zero released motion retains the unloaded zero-force convention. This
    /// local law supplies no body trajectory, normal impulse or arrest decision.
    pub fn advance_with_compliance(
        &self,
        state: &EllipticCoulombState,
        input: &TangentialStepInput,
        compliance: SymmetricTangentialCompliance,
    ) -> Result<TangentialStep<EllipticCoulombState>, TriboError> {
        self.check(state)?;
        if input.compliance_m_per_n != 0.0 {
            return Err(invalid("elliptic_competing_compliance"));
        }
        let m = compliance.response.matrix;
        if m[0][1] == 0.0 && m[0][0] == m[1][1] {
            return self.advance(
                state,
                &TangentialStepInput {
                    compliance_m_per_n: m[0][0],
                    ..*input
                },
            );
        }
        input.solve.validate()?;
        nonnegative_finite(input.normal_n, "normal_n")?;
        positive_finite(input.dt_s, "dt_s")?;
        if let Some(t) = input.contact_temperature_k {
            positive_finite(t, "contact_temperature_k")?;
        }
        if let Some(p) = input.nominal_pressure_pa {
            nonnegative_finite(p, "nominal_pressure_pa")?;
        }
        let free = to_principal(state.axis, input.slip_increment_m)?;
        let static_caps = capacities(self.static_mu, input.normal_n)?;
        let kinetic_caps = capacities(self.kinetic_mu, input.normal_n)?;
        let (c, determinant) = compliance.material(state.axis)?;
        let (force, slip, regime, iterations) = resolve_matrix(
            free,
            static_caps,
            kinetic_caps,
            c,
            determinant,
            compliance.response,
            state.axis,
            input.slip_increment_m,
            input.normal_n > 0.0,
            ProjectionControls::from(input.solve),
        )?;
        let traction_n = from_principal(state.axis, force)?;
        let slip_increment_m = if force == [0.0; 2] {
            input.slip_increment_m
        } else {
            from_principal(state.axis, slip)?
        };
        let force_norm = det::hypot(traction_n[0], traction_n[1]);
        if !force_norm.is_finite() || force_norm > input.solve.maximum_force_n {
            return Err(invalid("elliptic_maximum_force_n"));
        }
        let mut next = state.clone();
        next.contact_temperature_k = input.contact_temperature_k;
        let mut compatibility_residual_m = [0.0; 2];
        for i in 0..2 {
            next.slip_m[i] += slip_increment_m[i];
            finite(next.slip_m[i], "elliptic_accumulated_slip_m")?;
            compatibility_residual_m[i] = slip_increment_m[i]
                - (input.slip_increment_m[i] - (m[i][0] * traction_n[0] + m[i][1] * traction_n[1]));
            finite(
                compatibility_residual_m[i],
                "elliptic_compatibility_residual_m",
            )?;
        }
        let dissipated_j = force[0] * slip[0] + force[1] * slip[1];
        nonnegative_finite(dissipated_j, "elliptic_dissipated_j")?;
        let actual_work_j =
            traction_n[0] * slip_increment_m[0] + traction_n[1] * slip_increment_m[1];
        finite(actual_work_j, "elliptic_actual_work_j")?;
        let closure_residual_j = dissipated_j - actual_work_j;
        finite(closure_residual_j, "elliptic_closure_residual_j")?;
        Ok(TangentialStep {
            traction_n,
            slip_increment_m,
            state: next,
            regime,
            stored_change_j: 0.0,
            dissipated_j,
            closure_residual_j,
            compatibility_residual_m,
            iterations,
        })
    }

    /// Kinetic continuation with full body recoil, retaining the original static
    /// coefficients in the returned checkpoint. The mechanics owner admits rest.
    pub fn advance_kinetic_with_compliance(
        &self,
        state: &EllipticCoulombState,
        input: &TangentialStepInput,
        compliance: SymmetricTangentialCompliance,
    ) -> Result<TangentialStep<EllipticCoulombState>, TriboError> {
        self.check(state)?;
        let kinetic = Self {
            static_mu: self.kinetic_mu,
            ..*self
        };
        let mut history = state.clone();
        history.owner = kinetic;
        let mut step = kinetic.advance_with_compliance(&history, input, compliance)?;
        step.state.owner = *self;
        Ok(step)
    }
}

fn matrix_slip(
    free: [f64; 2],
    force: [f64; 2],
    c: [f64; 3],
    scale: f64,
) -> Result<[f64; 2], TriboError> {
    let slip = [
        free[0] - scale * (c[0] * force[0] + c[1] * force[1]),
        free[1] - scale * (c[1] * force[0] + c[2] * force[1]),
    ];
    for x in slip {
        finite(x, "elliptic_resolved_slip_m")?;
    }
    Ok(slip)
}

// Positive eigenvalues and a unit maximum-eigenvalue axis. The determinant is
// retained from the original Gram factor, avoiding cancellation of a*d-b*b.
pub(super) fn spectrum(c: [f64; 3], determinant: f64) -> Result<([f64; 2], [f64; 2]), TriboError> {
    let [a, b, d] = c;
    if b == 0.0 {
        return Ok(if a >= d {
            ([a, d], [1.0, 0.0])
        } else {
            ([d, a], [0.0, 1.0])
        });
    }
    let delta = 0.5 * (a - d);
    let gap = det::hypot(delta, b);
    let largest = 0.5 * a + 0.5 * d + gap;
    positive_finite(largest, "elliptic_compliance_eigenvalue")?;
    let axis = if delta >= 0.0 {
        [gap + delta, b]
    } else {
        [b, gap - delta]
    };
    let length = det::hypot(axis[0], axis[1]);
    positive_finite(length, "elliptic_compliance_eigenaxis")?;
    Ok(([largest, determinant / largest], axis.map(|v| v / length)))
}

#[allow(clippy::too_many_arguments)]
pub(super) fn resolve_matrix(
    free: [f64; 2],
    static_caps: [f64; 2],
    caps: [f64; 2],
    c: [f64; 3],
    determinant: f64,
    matrix: SymmetricResponse,
    material_axis: [f64; 2],
    free_original: [f64; 2],
    normal_active: bool,
    solve: ProjectionControls,
) -> Result<Resolution, TriboError> {
    if !normal_active || static_caps == [0.0; 2] {
        return Ok(([0.0; 2], free, FrictionRegime::Sliding, 0));
    }
    if free == [0.0; 2] {
        return Ok(([0.0; 2], free, FrictionRegime::Sticking, 0));
    }
    if matrix.determinant > 0.0 {
        let [a, b, d] = matrix.normalized;
        let g = free_original.map(|v| v / matrix.scale);
        let required = [
            (d * g[0] - b * g[1]) / matrix.determinant,
            (a * g[1] - b * g[0]) / matrix.determinant,
        ];
        let required = to_principal(material_axis, required)?;
        if inside(required, static_caps) {
            return Ok((required, [0.0; 2], FrictionRegime::Sticking, 0));
        }
    } else if static_caps[0] == 0.0 || static_caps[1] == 0.0 {
        let active = usize::from(static_caps[0] == 0.0);
        let direction = if active == 0 {
            material_axis
        } else {
            [-material_axis[1], material_axis[0]]
        };
        let m = matrix.matrix;
        let column = [
            m[0][0] * direction[0] + m[0][1] * direction[1],
            m[1][0] * direction[0] + m[1][1] * direction[1],
        ];
        let pivot = usize::from(column[1].abs() > column[0].abs());
        if column[pivot] != 0.0 {
            let required = free_original[pivot] / column[pivot];
            let consistent = (0..2).all(|i| {
                (free_original[i] - column[i] * required).abs()
                    <= 64.0 * f64::EPSILON * (free_original[i].abs() + (column[i] * required).abs())
            });
            if required.abs() <= static_caps[active] && consistent {
                let mut force = [0.0; 2];
                force[active] = required;
                return Ok((force, [0.0; 2], FrictionRegime::Sticking, 0));
            }
        }
    } else {
        let [a, b, d] = matrix.normalized;
        let pivot = usize::from(d > a);
        let other = 1 - pivot;
        let g = free_original;
        let null_residual = g[other] - b * g[pivot];
        if null_residual.abs() <= 64.0 * f64::EPSILON * (g[other].abs() + (b * g[pivot]).abs()) {
            let direction = if pivot == 0 { [1.0, b] } else { [b, 1.0] };
            let direction = to_principal(material_axis, direction)?;
            let capacity = static_caps[0].max(static_caps[1]);
            let support = capacity
                * det::hypot(
                    (static_caps[0] / capacity) * direction[0],
                    (static_caps[1] / capacity) * direction[1],
                );
            if (g[pivot] / matrix.scale).abs() <= support {
                return Err(invalid("elliptic_rank_deficient_static_reaction"));
            }
        }
    }
    if caps == [0.0; 2] {
        return Ok(([0.0; 2], free, FrictionRegime::Sliding, 0));
    }
    if caps[0] == 0.0 || caps[1] == 0.0 {
        let active = usize::from(caps[0] == 0.0);
        let response = matrix.scale * c[2 * active];
        let mut force = [0.0; 2];
        let mut slip;
        if response == 0.0 {
            if free[active] == 0.0 {
                return Err(invalid("elliptic_nonunique_kinetic_reaction"));
            }
            force[active] = caps[active].copysign(free[active]);
            slip = matrix_slip(free, force, c, matrix.scale)?;
        } else {
            let required = free[active] / response;
            force[active] = required.clamp(-caps[active], caps[active]);
            slip = matrix_slip(free, force, c, matrix.scale)?;
            if required.abs() <= caps[active] {
                slip[active] = 0.0;
            }
        }
        return Ok((force, slip, FrictionRegime::Sliding, 0));
    }
    let capacity = caps[0].max(caps[1]);
    let distance = free[0].abs().max(free[1].abs());
    let r = caps.map(|v| v / capacity);
    for x in r {
        if x * x == 0.0 {
            return Err(invalid("elliptic_unresolved_axis_scale"));
        }
    }
    let h = matrix.scale * (capacity / distance);
    positive_finite(h, "elliptic_scaled_compliance")?;
    let weighted = [c[0] * r[0] * r[0], c[1] * r[0] * r[1], c[2] * r[1] * r[1]];
    let weighted_det = determinant * r[0] * r[0] * r[1] * r[1];
    if determinant > 0.0 && weighted_det == 0.0 {
        return Err(invalid("elliptic_unresolved_axis_scale"));
    }
    let (eigenvalues, axis) = spectrum(weighted, weighted_det)?;
    let rhs = [r[0] * (free[0] / distance), r[1] * (free[1] / distance)];
    let rhs_eigen = to_principal(axis, rhs)?;
    let mut lower = 0.0;
    let mut upper = det::hypot(rhs[0], rhs[1]);
    for iteration in 1..=solve.maximum_iterations {
        let q = f64::midpoint(lower, upper);
        let mut y = [0.0; 2];
        for i in 0..2 {
            let denominator = h * eigenvalues[i] + q;
            positive_finite(denominator, "elliptic_projection_denominator")?;
            y[i] = rhs_eigen[i] / denominator;
        }
        let y = from_principal(axis, y)?;
        let force = [caps[0] * y[0], caps[1] * y[1]];
        let slip = matrix_slip(free, force, c, matrix.scale)?;
        let expected = support(caps, slip)?;
        let residual = det::hypot(force[0] - expected[0], force[1] - expected[1]);
        let tolerance = solve.tolerance(
            det::hypot(force[0], force[1]),
            det::hypot(expected[0], expected[1]),
        )?;
        if residual <= tolerance && (0..2).all(|i| force[i] * slip[i] >= 0.0) && slip != [0.0; 2] {
            return Ok((force, slip, FrictionRegime::Sliding, iteration));
        }
        if q == lower || q == upper {
            break;
        }
        if det::hypot(y[0], y[1]) > 1.0 {
            lower = q;
        } else {
            upper = q;
        }
    }
    Err(invalid("elliptic_sliding_root_unresolved"))
}
