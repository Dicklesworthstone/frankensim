//! Original nodal multi-tone loading and its work-conjugate observations.
use super::*;

/// Accepted continuous-tone solid interval. All moments retain the original
/// full nodal layout; fixed support rows have zero velocity and do no work.
#[derive(Clone, Debug)]
pub struct PreparedTensorSolidToneStep {
    pub motion: PreparedTensorSolidStep,
    /// Integral_0^1 cos(omega*h*s) v_node(h*s) ds [m/s], in tone order.
    pub cosine_nodal_velocity_moment_mps: Vec<Vec<[f64; 3]>>,
    /// Integral_0^1 sin(omega*h*s) v_node(h*s) ds [m/s], in tone order.
    pub sine_nodal_velocity_moment_mps: Vec<Vec<[f64; 3]>>,
}

impl<'s> PreparedTensorSolidBody<'s> {
    /// Prepare continuous multi-tone nodal forces on the same immutable source,
    /// complete coordinate basis, tensor histories and per-element heat owner.
    /// Each frequency must lie inside every original material law's band.
    /// The force scale conditions source coordinates only and has units N.
    pub fn new_tones(
        projection: ProjectedViscoSolid<'s>,
        dt_s: f64,
        angular_frequencies_rad_s: &[f64],
        force_scale_n: f64,
        controls: PreparedTensorControls,
        gate: &CancelGate,
    ) -> Result<Self> {
        Self::prepare(
            projection,
            dt_s,
            Some(force_scale_n),
            Some(angular_frequencies_rad_s),
            controls,
            gate,
        )
    }

    pub fn angular_frequencies_rad_s(&self) -> Option<&[f64]> {
        self.body.angular_frequencies_rad_s()
    }

    /// Move unchanged nodal force amplitudes to the next interval origin.
    /// Returned phases are independent source data; mechanical history is unchanged.
    pub fn advance_force_phase(
        &self,
        cosine_nodal_force_n: &[Vec<[f64; 3]>],
        sine_nodal_force_n: &[Vec<[f64; 3]>],
        gate: &CancelGate,
    ) -> Result<(Vec<Vec<[f64; 3]>>, Vec<Vec<[f64; 3]>>)> {
        poll(gate)?;
        let nodes = self.projection.source().reference_nodes_m().len();
        let frequencies = self
            .angular_frequencies_rad_s()
            .ok_or(PreparedSolidError::Invalid(
                "continuous-tone preparation required",
            ))?;
        if cosine_nodal_force_n.len() != frequencies.len()
            || sine_nodal_force_n.len() != frequencies.len()
            || cosine_nodal_force_n
                .iter()
                .chain(sine_nodal_force_n)
                .any(|v| v.len() != nodes)
        {
            return Err(PreparedSolidError::Invalid(
                "complete nodal coefficient fields required",
            ));
        }
        let flatten = |forces: &[Vec<[f64; 3]>]| {
            forces
                .iter()
                .map(|force| force.iter().flatten().copied().collect::<Vec<_>>())
                .collect::<Vec<_>>()
        };
        let (cosine, sine) = self.body.rotate_force_phase(
            &flatten(cosine_nodal_force_n),
            &flatten(sine_nodal_force_n),
            3 * nodes,
            gate,
        )?;
        let expand = |forces: Vec<Vec<f64>>| {
            forces
                .into_iter()
                .map(|force| {
                    force
                        .chunks_exact(3)
                        .map(|v| [v[0], v[1], v[2]])
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>()
        };
        poll(gate)?;
        Ok((expand(cosine), expand(sine)))
    }
}

impl PreparedTensorSolidHistory<'_, '_> {
    /// Advance all supplied original nodal harmonics simultaneously:
    /// F(t0+t)=sum(A_j*cos(omega_j*t)+B_j*sin(omega_j*t)).
    ///
    /// Moment expansion uses the same displacement basis as force projection.
    /// Independent continuous source work must equal the original nodal
    /// coefficient/moment contraction before any accepted history is replaced.
    pub fn advance_tones(
        &mut self,
        cosine_nodal_force_n: &[Vec<[f64; 3]>],
        sine_nodal_force_n: &[Vec<[f64; 3]>],
        gate: &CancelGate,
    ) -> Result<PreparedTensorSolidToneStep> {
        poll(gate)?;
        let frequencies =
            self.solid
                .angular_frequencies_rad_s()
                .ok_or(PreparedSolidError::Invalid(
                    "continuous-tone preparation required",
                ))?;
        if cosine_nodal_force_n.len() != frequencies.len()
            || sine_nodal_force_n.len() != frequencies.len()
        {
            return Err(PreparedSolidError::Invalid(
                "one nodal coefficient pair per tone required",
            ));
        }
        let projection = &self.solid.projection;
        let mut cosine = Vec::with_capacity(frequencies.len());
        let mut sine = Vec::with_capacity(frequencies.len());
        for (a, b) in cosine_nodal_force_n.iter().zip(sine_nodal_force_n) {
            poll(gate)?;
            cosine.push(projection.project_nodal_force(a, gate)?);
            sine.push(projection.project_nodal_force(b, gate)?);
        }
        let mut candidate = self.inner.clone();
        let step = candidate.advance_tones(&cosine, &sine, gate)?;
        let mean_nodal_velocity_mps =
            projection.expand_nodal(&step.motion.mean_velocity_mps, gate)?;
        let mut cosine_moments = Vec::with_capacity(frequencies.len());
        let mut sine_moments = Vec::with_capacity(frequencies.len());
        let (mut work, mut correction) = (0.0, 0.0);
        let dt = self.solid.body.duration_s();
        for tone in 0..frequencies.len() {
            poll(gate)?;
            let a = projection.expand_nodal(&step.cosine_velocity_moment_mps[tone], gate)?;
            let b = projection.expand_nodal(&step.sine_velocity_moment_mps[tone], gate)?;
            for (((fa, fb), va), vb) in cosine_nodal_force_n[tone]
                .iter()
                .zip(&sine_nodal_force_n[tone])
                .zip(&a)
                .zip(&b)
            {
                poll(gate)?;
                for axis in 0..3 {
                    for term in [dt * fa[axis] * va[axis], dt * fb[axis] * vb[axis]] {
                        let increment = term - correction;
                        let next = work + increment;
                        correction = (next - work) - increment;
                        work = next;
                    }
                }
            }
            cosine_moments.push(a);
            sine_moments.push(b);
        }
        let defect_j = work - step.motion.energy.supplied_work_j;
        if !defect_j.is_finite() || defect_j.abs() > self.solid.work_allowance_j {
            return Err(PreparedSolidError::Work {
                defect_j,
                allowance_j: self.solid.work_allowance_j,
            });
        }
        poll(gate)?;
        self.inner = candidate;
        Ok(PreparedTensorSolidToneStep {
            motion: PreparedTensorSolidStep {
                material: step.motion,
                mean_nodal_velocity_mps,
                nodal_work_error_j: defect_j,
            },
            cosine_nodal_velocity_moment_mps: cosine_moments,
            sine_nodal_velocity_moment_mps: sine_moments,
        })
    }
}
