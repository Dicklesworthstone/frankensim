//! Bounded Robin coordinates over the original FMM layers and Krylov owner.
//! The unknown x satisfies (A diag(a) + ik B diag(b)) x = rhs, with
//! a=(1+r)/2, b=(1-r)/2. Its residual is exactly A p - ik B (rho*c*v)
//! after the prescribed particular trace is restored. No admittance is formed.

use super::*;
use crate::ComplexSum;
use crate::helmholtz::bc::{
    Admittance, BoundarySolution, IncidentTraces, PreparedIterativeBoundary, ScatteringSolution,
    validate_closed,
};

enum ReflectionFactor<'a> {
    Rigid(PreparedIterativeBoundary<'a>),
    Local(IterativeRadiation<'a>),
}

/// Reusable matrix-free local boundary, including exact pressure release.
///
/// Reflection uses exp(-i omega t). Passive coefficients satisfy |r|<=1
/// within roundoff; coefficients and wall resistance are never repaired.
/// The original Burton-Miller FMM, panel rule, work budgets and GMRES owner
/// remain authoritative. Near-field right-preconditioner blocks are those
/// of the transformed equation, including at a=0; no dense fallback exists.
/// The reported residual measures the boundary equation under the selected
/// FMM approximation. It is not a bound on FMM or spatial discretization error.
pub struct PreparedIterativeReflectionBoundary<'a> {
    surface: &'a SpherePanels,
    k: f64,
    medium: Medium,
    pressure_coefficient: Vec<C64>,
    flow_coefficient: Vec<C64>,
    inverse_norm: Vec<f64>,
    declared_active: bool,
    factor: ReflectionFactor<'a>,
}

impl<'a> PreparedIterativeReflectionBoundary<'a> {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        surface: &'a SpherePanels,
        k: f64,
        medium: Medium,
        reflection: &[C64],
        declared_active: bool,
        formulation: Formulation,
        options: &IterativeOptions,
        gate: &CancelGate,
    ) -> Result<Self, IterativeError> {
        let started = Instant::now();
        poll(gate)?;
        admit(surface, k, medium, options)?;
        validate_closed(surface)?;
        if formulation == Formulation::PlainCbie {
            return Err(invalid(
                "reflection boundaries keep the Burton-Miller iterative owner",
            ));
        }
        let rho_c = medium.density * medium.sound_speed;
        if !rho_c.is_finite() || rho_c <= 0.0 || reflection.len() != surface.areas().len() {
            return Err(invalid(
                "finite characteristic impedance and one reflection per panel required",
            ));
        }
        let mut pressure_coefficient = Vec::with_capacity(reflection.len());
        let mut flow_coefficient = Vec::with_capacity(reflection.len());
        let mut inverse_norm = Vec::with_capacity(reflection.len());
        for (index, &r) in reflection.iter().enumerate() {
            if index % 256 == 0 {
                poll(gate)?;
            }
            if !finite(r) || (!declared_active && r.abs() > 1.0 + 64.0 * f64::EPSILON) {
                return Err(invalid(
                    "one finite, declared-passive reflection per panel required",
                ));
            }
            let a = (C64::ONE + r).scale(0.5);
            let b = (C64::ONE - r).scale(0.5);
            let inverse = 1.0 / (a.norm_sq() + b.norm_sq());
            if !inverse.is_finite() || inverse <= 0.0 {
                return Err(invalid("unrepresentable reflection normalization"));
            }
            pressure_coefficient.push(a);
            flow_coefficient.push(b);
            inverse_norm.push(inverse);
        }
        let factor = if reflection.iter().all(|&r| r == C64::ONE) {
            ReflectionFactor::Rigid(PreparedIterativeBoundary::new(
                surface,
                k,
                medium,
                Admittance::Rigid,
                formulation,
                options,
                gate,
            )?)
        } else {
            // Build the layer owner once, without requiring an A-only block
            // to be invertible before the actual Robin block is considered.
            let mut setup_options = options.clone();
            setup_options.precond = Preconditioner::None;
            let mut operator =
                IterativeRadiation::prepare(surface, k, medium, formulation, &setup_options, gate)?;
            prepare_blocks(
                &mut operator,
                &pressure_coefficient,
                &flow_coefficient,
                options.precond,
                gate,
            )?;
            operator.options.precond = options.precond;
            for word in [0x5245_464c, u64::from(declared_active)] {
                operator.stamp = (operator.stamp ^ word).wrapping_mul(0x100000001b3);
            }
            for r in reflection {
                for word in [r.re.to_bits(), r.im.to_bits()] {
                    operator.stamp = (operator.stamp ^ word).wrapping_mul(0x100000001b3);
                }
            }
            operator.setup_s = started.elapsed().as_secs_f64();
            ReflectionFactor::Local(operator)
        };
        poll(gate)?;
        Ok(Self {
            surface,
            k,
            medium,
            pressure_coefficient,
            flow_coefficient,
            inverse_norm,
            declared_active,
            factor,
        })
    }

    /// Prescribed outward drive velocity. At r=-1 its pressure work vanishes;
    /// the finite-impedance drive convention is retained on every other patch.
    pub fn solve(
        &self,
        drive_velocity: &[C64],
        gate: &CancelGate,
    ) -> Result<(BoundarySolution, SolveReport, ResourceReceipt), IterativeError> {
        poll(gate)?;
        if drive_velocity.len() != self.surface.areas().len()
            || drive_velocity.iter().any(|&v| !finite(v))
        {
            return Err(invalid(
                "reflection boundary needs one finite drive per panel",
            ));
        }
        if let ReflectionFactor::Rigid(owner) = &self.factor {
            let (mut response, report, receipt) = owner.solve(drive_velocity, gate)?;
            response.declared_active = self.declared_active;
            return Ok((response, report, receipt));
        }
        let rho_c = self.medium.density * self.medium.sound_speed;
        let h: Vec<_> = drive_velocity
            .iter()
            .zip(&self.pressure_coefficient)
            .map(|(&v, &a)| a * v.scale(rho_c))
            .collect();
        let (radiation, mut amplitude, report, receipt) = self.trace(&h, gate)?;
        let mut driving_work = ComplexSum::default();
        for i in 0..amplitude.len() {
            if i % 256 == 0 {
                poll(gate)?;
            }
            amplitude[i] = amplitude[i]
                + self.flow_coefficient[i].conj()
                    * drive_velocity[i].scale(rho_c * self.inverse_norm[i]);
            driving_work.add(
                (radiation.pressure[i] * drive_velocity[i].conj())
                    .scale(0.5 * self.surface.areas()[i]),
            );
        }
        let driving_power_w = driving_work.finish().re;
        let absorbed_power_w = self.wall_work(&amplitude, gate)?;
        if !driving_power_w.is_finite() {
            return Err(invalid("nonfinite reflection boundary drive work"));
        }
        Ok((
            BoundarySolution {
                radiation,
                absorbed_power_w,
                driving_power_w,
                declared_active: self.declared_active,
            },
            report,
            receipt,
        ))
    }

    /// Both incident pressure and normal velocity enter the original equation.
    /// Wall work uses the total coherent field, including interference.
    pub fn scatter(
        &self,
        incident: &IncidentTraces,
        gate: &CancelGate,
    ) -> Result<(ScatteringSolution, SolveReport, ResourceReceipt), IterativeError> {
        poll(gate)?;
        incident.validate_binding(self.surface, self.k, self.medium)?;
        if let ReflectionFactor::Rigid(owner) = &self.factor {
            let (mut response, report, receipt) = owner.scatter(incident, gate)?;
            response.declared_active = self.declared_active;
            return Ok((response, report, receipt));
        }
        let rho_c = self.medium.density * self.medium.sound_speed;
        let h: Vec<_> = (0..self.surface.areas().len())
            .map(|i| {
                -(self.flow_coefficient[i] * incident.pressure_pa()[i]
                    + self.pressure_coefficient[i] * incident.normal_velocity_m_s()[i].scale(rho_c))
            })
            .collect();
        let (scattered, mut amplitude, report, receipt) = self.trace(&h, gate)?;
        for i in 0..amplitude.len() {
            if i % 256 == 0 {
                poll(gate)?;
            }
            amplitude[i] = amplitude[i]
                + (self.pressure_coefficient[i].conj() * incident.pressure_pa()[i]
                    - self.flow_coefficient[i].conj()
                        * incident.normal_velocity_m_s()[i].scale(rho_c))
                .scale(self.inverse_norm[i]);
        }
        let absorbed_power_w = self.wall_work(&amplitude, gate)?;
        Ok((
            ScatteringSolution {
                scattered,
                incident_pressure_pa: incident.pressure_pa().to_vec(),
                incident_normal_velocity_m_s: incident.normal_velocity_m_s().to_vec(),
                absorbed_power_w,
                declared_active: self.declared_active,
            },
            report,
            receipt,
        ))
    }

    fn trace(
        &self,
        h: &[C64],
        gate: &CancelGate,
    ) -> Result<(RadiationSolution, Vec<C64>, SolveReport, ResourceReceipt), IterativeError> {
        let started = Instant::now();
        let ReflectionFactor::Local(operator) = &self.factor else {
            return Err(invalid(
                "rigid reflection must use the original Neumann owner",
            ));
        };
        let n = h.len();
        let mut pressure = Vec::with_capacity(n);
        let mut flow = Vec::with_capacity(n);
        for (i, &value) in h.iter().enumerate() {
            if i % 256 == 0 {
                poll(gate)?;
            }
            let particular = value.scale(self.inverse_norm[i]);
            pressure.push(-self.flow_coefficient[i].conj() * particular);
            flow.push(self.pressure_coefficient[i].conj() * particular);
        }
        let rhs = self.image(operator, &pressure, &flow, gate)?;
        let mut state = operator.start_from_rhs(&rhs, &[], started, gate)?;
        operator.advance_with_operator(
            &mut state,
            operator.options.max_iters.saturating_add(1),
            gate,
            &mut |_| false,
            &|x, gate| {
                let p: Vec<_> = x
                    .iter()
                    .zip(&self.pressure_coefficient)
                    .map(|(&x, &a)| a * x)
                    .collect();
                let w: Vec<_> = x
                    .iter()
                    .zip(&self.flow_coefficient)
                    .map(|(&x, &b)| b * x)
                    .collect();
                self.image(operator, &p, &w, gate)
            },
        )?;
        let (x, report, receipt) = match operator.finish_coordinates(&state) {
            Ok(result) => result,
            Err(IterativeError::NotConverged {
                pressure: x,
                report,
                receipt,
            }) => {
                // The public diagnostic remains physical pressure, not the
                // internal mixed pressure/flow coordinate.
                return Err(IterativeError::NotConverged {
                    pressure: x
                        .iter()
                        .enumerate()
                        .map(|(i, &x)| {
                            self.pressure_coefficient[i] * x
                                + self.flow_coefficient[i].conj() * h[i].scale(self.inverse_norm[i])
                        })
                        .collect(),
                    report,
                    receipt,
                });
            }
            Err(error) => return Err(error),
        };
        let inverse_rho_c = 1.0 / (self.medium.density * self.medium.sound_speed);
        pressure.clear();
        flow.clear();
        for i in 0..n {
            if i % 256 == 0 {
                poll(gate)?;
            }
            let particular = h[i].scale(self.inverse_norm[i]);
            pressure.push(
                self.pressure_coefficient[i] * x[i] + self.flow_coefficient[i].conj() * particular,
            );
            flow.push(
                (-self.flow_coefficient[i] * x[i]
                    + self.pressure_coefficient[i].conj() * particular)
                    .scale(inverse_rho_c),
            );
        }
        if pressure.iter().chain(&flow).any(|&z| !finite(z)) {
            return Err(invalid("nonfinite reflection boundary trace"));
        }
        let (radiated_power, radiated_power_roundoff_interval, power_diagnostics) =
            evaluate_radiated_power(&pressure, &flow, self.surface.areas(), self.medium);
        poll(gate)?;
        Ok((
            RadiationSolution {
                pressure,
                velocity: flow,
                k: self.k,
                medium: self.medium,
                surface_fingerprint: surface_fingerprint(self.surface),
                panels_per_wavelength: operator.ppw,
                radiated_power,
                radiated_power_roundoff_interval,
                power_diagnostics,
                condition_lower_bound: 1.0,
                dense_cap_utilization: n as f64 / MAX_DENSE_PANELS as f64,
            },
            x,
            report,
            receipt,
        ))
    }

    // A*p + ik*B*w through the same prepared kernels and exact near rules.
    // Zero coordinate arms skip their FMM work at either endpoint.
    fn image(
        &self,
        operator: &IterativeRadiation<'_>,
        pressure: &[C64],
        flow: &[C64],
        gate: &CancelGate,
    ) -> Result<Vec<C64>, IterativeError> {
        poll(gate)?;
        let mut out = if pressure.iter().all(|&p| p == C64::ZERO) {
            vec![C64::ZERO; pressure.len()]
        } else {
            operator.apply(pressure, gate)?
        };
        if flow.iter().any(|&w| w != C64::ZERO) {
            let q: Vec<_> = flow.iter().map(|&w| w * C64::new(0.0, self.k)).collect();
            for (value, b) in out
                .iter_mut()
                .zip(operator.boundary_rhs_from_flux(&q, gate)?)
            {
                *value = *value + b;
            }
        }
        if out.iter().any(|&z| !finite(z)) {
            return Err(invalid(
                "nonfinite reflection boundary-integral application",
            ));
        }
        Ok(out)
    }

    fn wall_work(&self, amplitude: &[C64], gate: &CancelGate) -> Result<f64, IterativeError> {
        let inverse_rho_c = 1.0 / (self.medium.density * self.medium.sound_speed);
        let mut work = ComplexSum::default();
        let mut apparent = 0.0;
        for (i, (&x, &area)) in amplitude.iter().zip(self.surface.areas()).enumerate() {
            if i % 256 == 0 {
                poll(gate)?;
            }
            let p = self.pressure_coefficient[i] * x;
            let v = (self.flow_coefficient[i] * x).scale(inverse_rho_c);
            work.add((p * v.conj()).scale(0.5 * area));
            apparent += 0.5 * p.abs() * v.abs() * area;
        }
        let work = work.finish().re;
        if !work.is_finite()
            || !apparent.is_finite()
            || (!self.declared_active && work < -1e-10 * apparent.max(f64::MIN_POSITIVE))
        {
            return Err(invalid("reflection boundary failed passive-work check"));
        }
        Ok(work)
    }
}

fn invalid(what: &'static str) -> IterativeError {
    IterativeError::Invalid { what }
}

fn prepare_blocks(
    operator: &mut IterativeRadiation<'_>,
    a: &[C64],
    b: &[C64],
    preconditioner: Preconditioner,
    gate: &CancelGate,
) -> Result<(), IterativeError> {
    let Preconditioner::NearFieldBlockJacobi { max_block_panels } = preconditioner else {
        return Ok(());
    };
    let n = operator.surface.areas().len();
    let groups: Vec<Vec<usize>> = match &operator.layers {
        Layers::BurtonMiller { double, .. } => double
            .near_field_blocks()
            .flat_map(|leaf| leaf.chunks(max_block_panels).map(<[usize]>::to_vec))
            .collect(),
        Layers::BurtonMillerCloud(cloud) => {
            CloudPartition::new(&cloud.double, n, cloud.weights.len())
                .leaf_blocks(n, max_block_panels)
        }
        Layers::PlainCbie(_) => {
            return Err(invalid("reflection preconditioner requires Burton-Miller"));
        }
    };
    let near_entries = groups
        .iter()
        .try_fold(operator.near_entries, |count, group| {
            group
                .len()
                .checked_mul(group.len())
                .and_then(|size| count.checked_add(size))
        })
        .ok_or_else(|| invalid("reflection near-field extent overflow"))?;
    if near_entries > operator.options.max_near_entries {
        return Err(IterativeError::NearFieldBudget {
            required_entries: near_entries,
            maximum_entries: operator.options.max_near_entries,
        });
    }
    for indices in groups {
        poll(gate)?;
        let m = indices.len();
        let mut matrix = vec![C64::ZERO; m * m];
        for (row, &i) in indices.iter().enumerate() {
            poll(gate)?;
            for (column, &j) in indices.iter().enumerate() {
                let (left, right) = near_pair(operator, i, j)?;
                matrix[row * m + column] = left * a[j] + C64::new(0.0, operator.k) * right * b[j];
            }
        }
        if matrix.iter().any(|&z| !finite(z)) {
            return Err(invalid("nonfinite transformed near-field block"));
        }
        let inverse = lu_complex(&matrix, m).map_err(|_| IterativeError::SingularNearBlock {
            block: operator.blocks.len(),
        })?;
        operator.blocks.push(NearBlock { indices, inverse });
    }
    operator.near_entries = near_entries;
    operator.stamp = (operator.stamp ^ max_block_panels as u64).wrapping_mul(0x100000001b3);
    poll(gate)
}

// Entries are needed only within the same bounded FMM leaf groups. These
// are the original operator's self closures and direct-interaction rows.
fn near_pair(
    operator: &IterativeRadiation<'_>,
    i: usize,
    j: usize,
) -> Result<(C64, C64), IterativeError> {
    let surface = operator.surface;
    match &operator.layers {
        Layers::BurtonMiller {
            diagonal_a,
            diagonal_b,
            near,
            ..
        } => {
            let (mut a, mut b) = if i == j {
                (diagonal_a[i], diagonal_b[i])
            } else {
                let (s, d, dp, h) = kernels(
                    operator.k,
                    surface.centroids()[i],
                    surface.normals()[i],
                    surface.centroids()[j],
                    surface.normals()[j],
                );
                (
                    (d + operator.alpha * h).scale(surface.areas()[j]),
                    (s + operator.alpha * dp).scale(surface.areas()[j]),
                )
            };
            if let Some(table) = near {
                a = a
                    + table.get(i, j, NearOperator::Double).unwrap_or(C64::ZERO)
                    + operator.alpha
                        * table
                            .get(i, j, NearOperator::Hypersingular)
                            .unwrap_or(C64::ZERO);
                b = b + table.get(i, j, NearOperator::Single).unwrap_or(C64::ZERO);
            }
            Ok((a, b))
        }
        Layers::BurtonMillerCloud(cloud) => {
            if i == j {
                return Ok((cloud.diagonal_a[i], cloud.diagonal_b[i]));
            }
            let entry = |kind| {
                cloud
                    .table
                    .get(i, j, kind)
                    .ok_or_else(|| invalid("reflection near table lacks a same-leaf panel"))
            };
            let (_, _, dp, _) = kernels(
                operator.k,
                surface.centroids()[i],
                surface.normals()[i],
                surface.centroids()[j],
                surface.normals()[j],
            );
            Ok((
                entry(NearOperator::Double)? + operator.alpha * entry(NearOperator::Hypersingular)?,
                entry(NearOperator::Single)? + operator.alpha * dp.scale(surface.areas()[j]),
            ))
        }
        Layers::PlainCbie(_) => Err(invalid("reflection near entries require Burton-Miller")),
    }
}

#[cfg(test)]
#[path = "helmholtz_bc_reflection_iterative_tests.rs"]
mod tests;
