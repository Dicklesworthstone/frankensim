//! Prescribed world-component motion on the original fixed reference mesh.
//!
//! Partitioning the full equation gives
//! M_ff a_f + C_ff v_f + K_ff u_f = f_f + fT_f
//!     - M_fc a_c - C_fc v_c - K_fc u_c.
//! The support acceleration follows the same Newmark predictors as the free
//! state. Thus the full trapezoidal work identity includes reaction work, with
//! no fictitious force or omission of consistent-mass coupling.

use super::*;

#[derive(Clone)]
pub(super) struct PrescribedMotion3 {
    /// Full node-major fields, zero on free components.
    pub(super) motion: [Vec<f64>; 3],
    thermal: Arc<[f64]>,
    /// Support force ON the solid at the accepted endpoint; zero on free DOFs.
    reactions: Vec<f64>,
}

impl PrescribedMotion3 {
    pub(super) fn initial(
        problem: &TetLinearElasticProblem<'_>,
        free: &[usize],
        u: &[[f64; 3]],
        v: &[[f64; 3]],
        a: &[[f64; 3]],
        thermal: Vec<f64>,
        cx: &Cx<'_>,
    ) -> Result<Self, TransientError3> {
        check(cx)?;
        if a.len() != problem.nodes_m.len() {
            return Err(TransientError3::Invalid(
                "prescribed initial acceleration needs one xyz vector per node",
            ));
        }
        for (node, row) in a.iter().enumerate() {
            check(cx)?;
            for c in 0..3 {
                number(row[c])?;
                if free.binary_search(&(3 * node + c)).is_ok() && row[c] != 0.0 {
                    return Err(TransientError3::Invalid(
                        "free initial acceleration is solved, not prescribed",
                    ));
                }
            }
        }
        let full = problem
            .nodes_m
            .len()
            .checked_mul(3)
            .ok_or(TransientError3::Budget)?;
        let mut motion = std::array::from_fn(|_| vec![0.0; full]);
        for &d in problem.fixed_dofs {
            check(cx)?;
            motion[0][d] = u[d / 3][d % 3];
            motion[1][d] = v[d / 3][d % 3];
            motion[2][d] = a[d / 3][d % 3];
        }
        let thermal = if thermal.is_empty() {
            vec![0.0; full]
        } else {
            thermal
        };
        check(cx)?;
        Ok(Self {
            motion,
            thermal: Arc::from(thermal),
            reactions: vec![0.0; full],
        })
    }
}

impl<'g, 'm> TetLoadCase3<'g, 'm> {
    /// Initialize a sparse inertial trajectory with prescribed support motion.
    /// Initial displacement and velocity are FULL node-major xyz fields [m,
    /// m/s]. `prescribed_acceleration_m_s2` supplies fixed-component acceleration
    /// and must be zero on free components; their initial acceleration is solved
    /// from the complete partitioned equation, including M_fc, C_fc and K_fc.
    ///
    /// Every later endpoint must use `advance_with_supports` (or its history
    /// adapter). Fixed components and material/thermal state remain frozen.
    /// Prescribed motion is absolute in the reference world's Cartesian axes;
    /// no accelerating-frame or velocity-integration convention is inferred.
    /// Existing displacement, element, solve and cancellation bounds apply.
    #[allow(clippy::too_many_arguments)]
    pub fn begin_transient_with_supports<'p, 's>(
        &self,
        problem: &'p TetLinearElasticProblem<'s>,
        thermal: Option<TetThermalStrainField<'_>>,
        time_s: f64,
        displacement_m: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        prescribed_acceleration_m_s2: &[[f64; 3]],
        config: TransientConfig3,
        cx: &Cx<'_>,
    ) -> Result<TetTransientRun3<'g, 'm, 'p, 's>, TransientError3> {
        self.begin_transient_impl(
            problem,
            thermal,
            time_s,
            displacement_m,
            velocity_m_s,
            Some(prescribed_acceleration_m_s2),
            config,
            cx,
        )
    }
}

struct Endpoint {
    motion: [Vec<f64>; 3],
    kinetic: f64,
    strain: f64,
    reactions: Vec<f64>,
}

impl TetTransientRun3<'_, '_, '_, '_> {
    /// Advance to explicitly prescribed absolute displacement components [m].
    /// DOFs must be sorted, unique and fixed in the original problem. An omitted
    /// fixed component is prescribed ZERO at this endpoint, not held at its old
    /// value. Support velocity and acceleration use beta=1/4, gamma=1/2 Newmark:
    /// a_c1=(u_c1-u_c0-h*v_c0-h²*a_c0/4)/(h²/4), v_c1=v_c0+h*(a_c0+a_c1)/2.
    /// Recovering acceleration amplifies displacement-sample roundoff by h^-2;
    /// arbitrarily small steps do not guarantee accurate prescribed derivatives.
    ///
    /// The same effective sparse M_ff/C_ff/K_ff solve advances free components.
    /// Full element actions supply inertial coupling, reactions, kinetic/strain
    /// energy, prescribed-load work and support work. No global dense matrix is
    /// allocated. Cancellation/refusal commits neither support nor free motion.
    /// A discontinuous displacement is not a finite force jump: resolve physical
    /// motion smoothly with an adequate step ladder. No impulse/event or temporal
    /// accuracy certificate is implied. Support-control adjoints are not supplied.
    pub fn advance_with_supports(
        &mut self,
        load_next: &TetLoadCase3<'_, '_>,
        prescribed_dofs_m: &[(usize, f64)],
        h_s: f64,
        cx: &Cx<'_>,
    ) -> Result<TransientStep3, TransientError3> {
        check(cx)?;
        let Some(previous) = &self.prescribed else {
            return Err(TransientError3::Invalid(
                "initialize prescribed-support motion before supplying endpoints",
            ));
        };
        if self.state.step >= self.config.maximum_steps {
            return Err(TransientError3::Budget);
        }
        let next_time = self.state.time_s + h_s;
        let b = 0.25 * h_s * h_s;
        if !h_s.is_finite()
            || h_s <= 0.0
            || !next_time.is_finite()
            || next_time <= self.state.time_s
            || !b.is_finite()
            || b == 0.0
        {
            return Err(TransientError3::Invalid(
                "prescribed step must reach a distinct finite time with representable squared step",
            ));
        }
        if prescribed_dofs_m.len() > self.problem.fixed_dofs.len() {
            return Err(TransientError3::Invalid(
                "too many prescribed support components",
            ));
        }
        let full = previous.motion[0].len();
        let mut prior = None;
        for &(d, u) in prescribed_dofs_m {
            check(cx)?;
            if d >= full
                || prior.is_some_and(|p| d <= p)
                || self.free_dofs().binary_search(&d).is_ok()
                || !u.is_finite()
                || u.abs() > self.config.maximum_displacement_m
            {
                return Err(TransientError3::Invalid(
                    "prescribed components must be finite, bounded, sorted, unique and fixed",
                ));
            }
            prior = Some(d);
        }
        let mut candidate = previous.clone();
        candidate.motion[0].fill(0.0);
        for &(d, u) in prescribed_dofs_m {
            check(cx)?;
            candidate.motion[0][d] = u;
        }
        for &d in self.problem.fixed_dofs {
            check(cx)?;
            let predictor = number(
                previous.motion[0][d] + h_s * previous.motion[1][d] + b * previous.motion[2][d],
            )?;
            candidate.motion[2][d] = number((candidate.motion[0][d] - predictor) / b)?;
            candidate.motion[1][d] = number(
                previous.motion[1][d]
                    + 0.5 * h_s * (previous.motion[2][d] + candidate.motion[2][d]),
            )?;
        }
        self.advance_impl(load_next, Some(candidate), h_s, cx)
    }

    pub(super) fn subtract_support_force(
        &self,
        supports: &PrescribedMotion3,
        rhs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), TransientError3> {
        let actions = self.problem.apply_full_kinematics(&supports.motion, cx)?;
        for (i, &d) in self.free_dofs().iter().enumerate() {
            if i % 256 == 0 {
                check(cx)?;
            }
            let forcing = number(
                actions.ma[d]
                    + self.config.mass_damping_per_s * actions.mv[d]
                    + self.config.stiffness_damping_s * actions.kv[d]
                    + actions.ku[d],
            )?;
            rhs[i] = number(rhs[i] - forcing)?;
        }
        check(cx)
    }

    fn full_motion(
        &self,
        supports: &PrescribedMotion3,
        reduced: [&[f64]; 3],
        cx: &Cx<'_>,
    ) -> Result<[Vec<f64>; 3], TransientError3> {
        let mut motion = supports.motion.clone();
        for (i, &d) in self.free_dofs().iter().enumerate() {
            if i % 256 == 0 {
                check(cx)?;
            }
            for field in 0..3 {
                motion[field][d] = reduced[field][i];
            }
        }
        for node in motion[0].chunks_exact(3) {
            check(cx)?;
            if node[0].hypot(node[1]).hypot(node[2]) > self.config.maximum_displacement_m {
                return Err(TransientError3::Invalid(
                    "full nodal motion exceeds displacement allowance",
                ));
            }
        }
        Ok(motion)
    }

    fn endpoint(
        &self,
        supports: &PrescribedMotion3,
        reduced: [&[f64]; 3],
        full_force: &[(usize, f64)],
        cx: &Cx<'_>,
    ) -> Result<Endpoint, TransientError3> {
        let motion = self.full_motion(supports, reduced, cx)?;
        let actions = self.problem.apply_full_kinematics(&motion, cx)?;
        let kinetic = number(0.5 * inner(&motion[1], &actions.mv))?;
        let strain = number(0.5 * inner(&motion[0], &actions.ku))?;
        let mut reactions = vec![0.0; motion[0].len()];
        for &d in self.problem.fixed_dofs {
            check(cx)?;
            reactions[d] = number(
                actions.ma[d]
                    + self.config.mass_damping_per_s * actions.mv[d]
                    + self.config.stiffness_damping_s * actions.kv[d]
                    + actions.ku[d]
                    - supports.thermal[d],
            )?;
        }
        for &(d, f) in full_force {
            check(cx)?;
            if self.free_dofs().binary_search(&d).is_err() {
                reactions[d] = number(reactions[d] - f)?;
            }
        }
        check(cx)?;
        Ok(Endpoint {
            motion,
            kinetic,
            strain,
            reactions,
        })
    }

    pub(super) fn initialize_support_energy(&mut self, cx: &Cx<'_>) -> Result<(), TransientError3> {
        let Some(supports) = &self.prescribed else {
            return Ok(());
        };
        let endpoint = self.endpoint(
            supports,
            [
                &self.state.displacement,
                &self.state.velocity,
                &self.state.acceleration,
            ],
            &self.full_force,
            cx,
        )?;
        let energy = number(endpoint.kinetic + endpoint.strain)?;
        check(cx)?;
        self.state.kinetic_energy_j = endpoint.kinetic;
        self.state.strain_energy_j = endpoint.strain;
        self.state.total_energy_j = energy;
        self.prescribed
            .as_mut()
            .expect("retained prescribed state")
            .reactions = endpoint.reactions;
        Ok(())
    }

    /// Full work identity. The original numerical state is untouched until the
    /// caller also accepts its effective-equation residual and energy allowance.
    pub(super) fn supported_accounting(
        &self,
        next: &mut PrescribedMotion3,
        reduced: [&[f64]; 3],
        load: &TetLoadCase3<'_, '_>,
        h: f64,
        cx: &Cx<'_>,
    ) -> Result<(f64, f64, f64, f64, f64, f64), TransientError3> {
        let previous = self.prescribed.as_ref().ok_or(TransientError3::Invalid(
            "missing previous prescribed support state",
        ))?;
        let old = self.full_motion(
            previous,
            [
                &self.state.displacement,
                &self.state.velocity,
                &self.state.acceleration,
            ],
            cx,
        )?;
        let endpoint = self.endpoint(next, reduced, load.nodal_forces_n(), cx)?;
        let full = old[0].len();
        let mut dq = vec![0.0; full];
        let mut midpoint = std::array::from_fn(|_| vec![0.0; full]);
        for d in 0..full {
            if d % 256 == 0 {
                check(cx)?;
            }
            dq[d] = number(endpoint.motion[0][d] - old[0][d])?;
            midpoint[1][d] = number(0.5 * old[1][d] + 0.5 * endpoint.motion[1][d])?;
        }
        let actions = self.problem.apply_full_kinematics(&midpoint, cx)?;
        let mut damping = 0.0;
        for d in 0..full {
            if d % 256 == 0 {
                check(cx)?;
            }
            let cv = number(
                self.config.mass_damping_per_s * actions.mv[d]
                    + self.config.stiffness_damping_s * actions.kv[d],
            )?;
            damping = number(midpoint[1][d].mul_add(cv, damping))?;
        }
        let damping = number(h * damping)?;
        let mut mechanical = 0.0;
        for forces in [&self.full_force[..], load.nodal_forces_n()] {
            for &(d, f) in forces {
                check(cx)?;
                mechanical = number((0.5 * f).mul_add(dq[d], mechanical))?;
            }
        }
        let thermal = number(inner(&next.thermal, &dq))?;
        let mut support = 0.0;
        for &d in self.problem.fixed_dofs {
            check(cx)?;
            let average = number(0.5 * previous.reactions[d] + 0.5 * endpoint.reactions[d])?;
            support = number(average.mul_add(dq[d], support))?;
        }
        check(cx)?;
        next.reactions = endpoint.reactions;
        Ok((
            endpoint.kinetic,
            endpoint.strain,
            mechanical,
            support,
            thermal,
            damping,
        ))
    }
}
