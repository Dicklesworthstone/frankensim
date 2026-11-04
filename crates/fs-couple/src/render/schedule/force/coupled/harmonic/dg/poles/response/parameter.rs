//! Native source-parameter directions on the complete loaded pole pencil.
//! Only the physical coefficients vary; numerical row scales and the mode's
//! pivot-one chart stay fixed. Every original acoustic/algebraic state survives.

use super::*;

/// One signed physical direction per caller-declared unit of theta. Source
/// indices retain construction order. These are derivatives, not replacements;
/// the caller owns the parameter units and the interval where its law is valid.
/// No dense pencil derivative or change of source topology is requested.
#[derive(Clone, Copy, Debug)]
pub enum DgPoleParameterDirection<'a> {
    /// Original spring and resistance derivatives at fixed work-conjugate
    /// columns. Rest extension does not enter this linear free-mode pencil.
    Connection {
        connection: usize,
        stiffness_n_m: f64,
        damping_n_s_m: f64,
    },
    /// Original Maxwell strength and relaxation-time derivatives. The complete
    /// sqrt(strength) chain in the energy-normalized memory coordinate is kept.
    /// Nonzero strength derivative at zero strength has no differentiable chart.
    Maxwell {
        arm: usize,
        stiffness: f64,
        relaxation_time_s: f64,
    },
    /// Complete d(flow_column)/dtheta in original mechanical mode order. Both
    /// modal pressure reaction and DG volume inflow change reciprocally.
    InletFlowColumn(&'a [f64]),
}

/// Local simple-pole and complete shape derivatives per unit theta. The SAME
/// prepared component is held at one, so its state derivative is exactly zero.
/// All other coordinates retain their original units; this is neither a moving
/// energy normalization nor a finite-amplitude power/storage interpretation.
#[derive(Clone, Debug)]
pub struct DgLoadedModalPoleDerivative {
    /// ds/dtheta for exp(s*t); d(omega)/dtheta = i*ds/dtheta.
    pub laplace_pole_per_s: C64,
    pub state: Vec<C64>,
    /// d(s*q)/dtheta, including pole motion as well as displacement motion.
    pub modal_velocity_m_sqrt_kg_per_s: Vec<C64>,
    pub radiation_velocity: Vec<C64>,
    /// Includes the explicit inlet-column derivative when it is the parameter.
    pub inlet_flow_m3_s: C64,
    pub inlet_pressure_pa: C64,
    /// Fixed complex bilinear weights^T*d(state)/dtheta, without conjugation.
    pub observations: Vec<C64>,
    pub linear_residual: f64,
    /// Independent original T*x_theta + (T_s*s_theta + T_theta)*x audit.
    /// A rowwise backward check, not a forward sensitivity/conditioning bound.
    pub original_relative_residual: f64,
}

// Original real columns share the existing double-double observation owner.
fn project(
    column: impl Iterator<Item = f64>,
    state: &[C64],
    gate: &CancelGate,
) -> Result<(C64, f64), DgPoleResponseError> {
    let (mut re, mut im, mut bound) = (Dd::ZERO, Dd::ZERO, Dd::ZERO);
    for (i, (w, &x)) in column.zip(state).enumerate() {
        if i % 64 == 0 {
            check_cancel(gate)?;
        }
        let weight = Dd::from_f64(w);
        re = re + weight * Dd::from_f64(x.re);
        im = im + weight * Dd::from_f64(x.im);
        bound = bound + weight.abs() * Dd::from_f64(x.abs());
    }
    let value = C64::new(re.to_f64(), im.to_f64());
    let bound = bound.to_f64();
    if !complex_finite(value) || !bound.is_finite() {
        return Err(DgPoleResponseError::Invalid(
            "native parameter contraction overflows",
        ));
    }
    Ok((value, bound))
}

impl PreparedDgPoleResponse<'_, '_> {
    /// Differentiate native connection, Maxwell or reciprocal inlet parameters
    /// through the SAME complete loaded system and cached bordered factor.
    /// The actual material/terminal/DG source bands remain in force. No fitted
    /// damping, diagonal loss projection, isolated impedance inverse, omitted
    /// field or finite-difference derivative is used.
    ///
    /// A nonzero connection-stiffness direction refuses if ANY prepared
    /// constant-loss connection is present: its g(K)b columns also depend on K,
    /// and their full matrix-function derivative is not supplied by this API.
    /// Resistance-only, Maxwell and inlet directions retain those columns
    /// exactly. Directions do not alter their fixed source projection vectors.
    ///
    /// Observations are complete, fixed complex raw-state weights as in
    /// forcing_residues. All outputs use the same fixed pivot-one gauge as
    /// right_mode(); no derivative of numerical row scaling is taken. A zero
    /// direction returns a zero tangent. Repeated/defective poles may refuse.
    ///
    /// In addition to the response admission, physical ceilings bound
    /// 16*(2*d^2+k*(8*d+m)+20*d)+512*(k+m) bytes and
    /// 256*(k+1)*d^2+64*k*(m+4)*d+512*d work units, with declared vector/sensor
    /// ceilings k,m. Sources/caller inputs and the numerical factor are excluded.
    /// A refused/cancelled batch publishes no partial result or state change.
    pub fn parameter_derivatives(
        &self,
        directions: &[DgPoleParameterDirection<'_>],
        observations: &[&[C64]],
        gate: &CancelGate,
    ) -> Result<Vec<DgLoadedModalPoleDerivative>, DgPoleResponseError> {
        check_cancel(gate)?;
        let d = self.source.dimension();
        let (n, k, m) = (
            d as u128,
            self.options.numerical.max_vectors as u128,
            self.options.max_observations as u128,
        );
        if 16 * (2 * n * n + k * (8 * n + m) + 20 * n) + 512 * (k + m)
            > self.options.max_physical_bytes as u128
            || 256 * (k + 1) * n * n + 64 * k * (m + 4) * n + 512 * n
                > self.options.max_physical_work as u128
            || directions.is_empty()
            || directions.len() > self.options.numerical.max_vectors
            || observations.len() > self.options.max_observations
            || observations
                .iter()
                .any(|w| w.len() != d || w.iter().any(|&v| !complex_finite(v)))
        {
            return Err(DgPoleResponseError::Invalid(
                "complete native parameter/observation batch exceeds its physical admission",
            ));
        }
        let mut original = Vec::with_capacity(directions.len());
        let mut mapped = Vec::with_capacity(directions.len());
        for &direction in directions {
            let action = self.parameter_action(direction, gate)?;
            let mut rhs = action.0.clone();
            self.source
                .pencil
                .scale_rhs(&mut rhs, &mut || gate.is_requested())
                .map_err(DgLoadedPoleError::from)?;
            original.push(action);
            mapped.push(rhs);
        }
        let refs: Vec<_> = mapped.iter().map(Vec::as_slice).collect();
        let numerical = self
            .numerical
            .parameter_derivatives(&refs, || gate.is_requested())?;
        let tolerance = self
            .options
            .numerical
            .linear_residual_tolerance
            .min(self.source.source.config.solve_relative_tolerance);
        let mut output = Vec::with_capacity(directions.len());
        for ((derivative, &direction), (action, envelope)) in
            numerical.into_iter().zip(directions).zip(original)
        {
            check_cancel(gate)?;
            let s = self.mode.laplace_pole_per_s;
            let ds = derivative.value;
            let relative = self
                .source
                .pencil
                .parameter_residual(
                    s,
                    &self.mode.state,
                    ds,
                    &derivative.vector,
                    &action,
                    &envelope,
                    &mut || gate.is_requested(),
                )
                .map_err(DgLoadedPoleError::from)?;
            if relative > tolerance {
                return Err(DgLoadedPoleError::Residual {
                    relative,
                    tolerance,
                }
                .into());
            }
            let velocity = |span: Range<usize>| {
                span.map(|i| s * derivative.vector[i] + ds * self.mode.state[i])
                    .collect::<Vec<_>>()
            };
            let modal = velocity(self.source.coordinates.modal_displacement.clone());
            let radiation = velocity(self.source.coordinates.radiation_displacement.clone());
            let mut flow = project(self.source.flow_column.iter().copied(), &modal, gate)?.0;
            if let DgPoleParameterDirection::InletFlowColumn(column) = direction {
                flow = flow
                    + project(
                        column.iter().copied(),
                        &self.mode.modal_velocity_m_sqrt_kg_per_s,
                        gate,
                    )?
                    .0;
            }
            let pressure = project(
                self.source.pressure_column.iter().copied(),
                &derivative.vector[self.source.coordinates.duct_energy.clone()],
                gate,
            )?
            .0;
            if modal.iter().chain(&radiation).any(|&v| !complex_finite(v)) || !complex_finite(flow)
            {
                return Err(DgPoleResponseError::Invalid(
                    "native parameter velocity or flow overflows",
                ));
            }
            let mut observed = Vec::with_capacity(observations.len());
            for &weights in observations {
                observed.push(contract(weights, &derivative.vector, gate)?);
            }
            output.push(DgLoadedModalPoleDerivative {
                laplace_pole_per_s: ds,
                state: derivative.vector,
                modal_velocity_m_sqrt_kg_per_s: modal,
                radiation_velocity: radiation,
                inlet_flow_m3_s: flow,
                inlet_pressure_pa: pressure,
                observations: observed,
                linear_residual: derivative.linear_residual,
                original_relative_residual: relative,
            });
        }
        check_cancel(gate)?;
        Ok(output)
    }

    // Physical T_theta*x at FIXED s, before the prepared row equilibration.
    // The envelope retains uncancelled source terms for the independent audit.
    fn parameter_action(
        &self,
        direction: DgPoleParameterDirection<'_>,
        gate: &CancelGate,
    ) -> Result<(Vec<C64>, Vec<f64>), DgPoleResponseError> {
        check_cancel(gate)?;
        let source = self.source.source;
        let n = source.mode_count();
        let q = &self.mode.state[self.source.coordinates.modal_displacement.clone()];
        let s = self.mode.laplace_pole_per_s;
        let mut rows = vec![(Dd::ZERO, Dd::ZERO, Dd::ZERO); self.source.dimension()];
        let mut add = |i: usize, scale: f64, value: C64, bound: f64| {
            let a = Dd::from_f64(scale);
            rows[i].0 = rows[i].0 + a * Dd::from_f64(value.re);
            rows[i].1 = rows[i].1 + a * Dd::from_f64(value.im);
            rows[i].2 = rows[i].2 + a.abs() * Dd::from_f64(bound);
        };
        match direction {
            DgPoleParameterDirection::Connection {
                connection,
                stiffness_n_m: dk,
                damping_n_s_m: dc,
            } => {
                if connection >= source.connections.len()
                    || !dk.is_finite()
                    || !dc.is_finite()
                    || (dk != 0.0
                        && source
                            .connections
                            .iter()
                            .any(|law| law.constant_loss_reference_omega_rad_s.is_some()))
                {
                    return Err(DgPoleResponseError::Invalid(
                        "finite original connection direction required; stiffness-dependent constant-loss columns need their full derivative",
                    ));
                }
                let b = &source.columns[connection];
                let c = source
                    .harmonic_damping_column(connection)
                    .map_err(DgLoadedPoleError::from)?;
                let (extension, be) = project(b.iter().copied(), q, gate)?;
                let (velocity, bc) = project(c.iter().copied(), q, gate)?;
                for i in 0..n {
                    if i % 64 == 0 {
                        check_cancel(gate)?;
                    }
                    add(i, dk * b[i], extension, be);
                    add(i, dc * c[i], s * velocity, s.abs() * bc);
                }
            }
            DgPoleParameterDirection::Maxwell {
                arm,
                stiffness: dk,
                relaxation_time_s: dt,
            } => {
                let arm_source =
                    source
                        .relaxation_branches
                        .get(arm)
                        .ok_or(DgPoleResponseError::Invalid(
                            "native Maxwell direction names an absent arm",
                        ))?;
                if !dk.is_finite() || !dt.is_finite() || (arm_source.stiffness == 0.0 && dk != 0.0)
                {
                    return Err(DgPoleResponseError::Invalid(
                        "finite Maxwell direction with a differentiable sqrt-strength chart required",
                    ));
                }
                let root = det::sqrt(arm_source.stiffness);
                let dr = if dk == 0.0 { 0.0 } else { 0.5 * (dk / root) };
                let rate = arm_source.relaxation_time_s.recip();
                let rate_derivative = -rate * (dt / arm_source.relaxation_time_s);
                let memory = self.source.coordinates.material_memory.start + arm;
                let z = self.mode.state[memory];
                let (extension, bound) =
                    project(arm_source.projection.iter().step_by(2).copied(), q, gate)?;
                for i in 0..n {
                    if i % 64 == 0 {
                        check_cancel(gate)?;
                    }
                    let b = arm_source.projection[2 * i];
                    add(i, dk * b, extension, bound);
                    add(i, -dr * b, z, z.abs());
                }
                add(memory, rate_derivative, z, z.abs());
                add(memory, -rate_derivative * root, extension, bound);
                add(memory, -rate * dr, extension, bound);
            }
            DgPoleParameterDirection::InletFlowColumn(column) => {
                if column.len() != n || column.iter().any(|v| !v.is_finite()) {
                    return Err(DgPoleResponseError::Invalid(
                        "native inlet direction requires the complete finite modal column",
                    ));
                }
                let (pressure, bound) = project(
                    self.source.pressure_column.iter().copied(),
                    &self.mode.state[self.source.coordinates.duct_energy.clone()],
                    gate,
                )?;
                for (i, &b) in column.iter().enumerate() {
                    add(i, b, pressure, bound);
                }
                let (volume, bound) = project(column.iter().copied(), q, gate)?;
                for (i, &p) in self.source.pressure_column.iter().enumerate() {
                    if i % 64 == 0 {
                        check_cancel(gate)?;
                    }
                    add(
                        self.source.coordinates.duct_energy.start + i,
                        -p,
                        s * volume,
                        s.abs() * bound,
                    );
                }
            }
        }
        let mut vector = Vec::with_capacity(rows.len());
        let mut envelope = Vec::with_capacity(rows.len());
        for (re, im, bound) in rows {
            let value = C64::new(re.to_f64(), im.to_f64());
            let bound = bound.to_f64();
            if !complex_finite(value) || !bound.is_finite() || bound < 0.0 {
                return Err(DgPoleResponseError::Invalid(
                    "native parameter action or envelope overflows",
                ));
            }
            vector.push(value);
            envelope.push(bound);
        }
        Ok((vector, envelope))
    }
}

#[cfg(test)]
mod tests;
