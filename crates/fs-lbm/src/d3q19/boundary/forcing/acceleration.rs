//! Density-dependent acceleration in the original Guo source history.

use super::*;

impl BoundaryGrid3 {
    /// Select a uniform held body acceleration, in lattice units `dx/dt²`.
    ///
    /// Each fluid cell supplies its actual density in `F = rho * a`. This
    /// replaces both previous force components and preserves physical density,
    /// momentum and force-corrected stress without advancing time or applying
    /// an impulse. Initialize density/velocity before this selection and rebind
    /// after any later manual density edit. Static solid cells carry no force.
    ///
    /// Advance with [`Self::step_with_endpoint_body_acceleration`] so the force
    /// follows transported density. Supply the same acceleration for a held
    /// load or the next time's acceleration for a smooth prescribed history.
    /// An ordinary `step` instead retains the selected force density.
    ///
    /// Closed walls and periodic boundaries are supported. Open-face density
    /// reconstruction and finite-inertia immersed-body feedback need their own
    /// coupled density-dependent source closure. Low-Mach lattice requirements
    /// remain unchanged; pressure, viscosity and wall reactions are not added
    /// again as body forces.
    pub fn rebind_body_acceleration(
        &mut self,
        acceleration: [f64; 3],
        gate: &CancelGate,
    ) -> Result<(), BoundaryForcingError3> {
        // The zero-motion frame already owns the atomic rho*a preparation and
        // physical Guo reframe, including reused population/force allocations.
        RigidFrame3::default().prepare_step(self, acceleration, gate)
    }

    /// Advance with prescribed acceleration at the next physical endpoint.
    ///
    /// Collision uses the accepted force `F_n`. After streaming, every fluid
    /// cell supplies its NEW density to `F_(n+1) = rho_(n+1) * acceleration`.
    /// The populations receive no source reframe at this elapsed endpoint:
    /// physical momentum `sum(f_q*c_q) + F/2` receives the trapezoidal impulse
    /// `(F_n + F_(n+1))/2`. Bind the initial acceleration with
    /// [`Self::rebind_body_acceleration`] to preserve its physical velocity.
    ///
    /// This is prescribed acceleration, not a held force density or an
    /// implicit velocity-dependent force. Closed/periodic boundaries and the
    /// original forcing-compatible collision laws are supported. Endpoint
    /// admission is fused into streaming; the accepted force allocation is
    /// reused after first use. Refusal/cancellation changes neither accepted
    /// populations nor forces. No additional grid clone or second collision
    /// is introduced.
    pub fn step_with_endpoint_body_acceleration(
        &mut self,
        acceleration: [f64; 3],
        gate: &CancelGate,
    ) -> Result<(), BoundaryStepError3> {
        let endpoint = self.endpoint_body_acceleration(acceleration, gate)?;
        self.step_with_endpoint_forcing_inner(endpoint, None, gate)
            .map(|_| ())
    }

    /// Acceleration step with the original selected-wall impulse, torque and
    /// work measured from the same collision that supplies the accepted flow.
    pub fn step_with_endpoint_body_acceleration_and_wall_momentum(
        &mut self,
        acceleration: [f64; 3],
        links: &[BoundaryLink3],
        moment_origin: [f64; 3],
        gate: &CancelGate,
    ) -> Result<WallMomentumExchange3, BoundaryStepError3> {
        self.validate_wall_measurement(links, moment_origin)?;
        let endpoint = self.endpoint_body_acceleration(acceleration, gate)?;
        self.step_with_endpoint_forcing_inner(endpoint, Some((links, moment_origin)), gate)
    }

    /// Pooled counterpart using identical native cell arithmetic and endpoint
    /// density-dependent force publication.
    pub fn step_pooled_with_endpoint_body_acceleration<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        acceleration: [f64; 3],
    ) -> Result<(), BoundaryStepError3> {
        let endpoint = self.endpoint_body_acceleration(acceleration, gate)?;
        self.step_pooled_coupling_inner(
            pool,
            gate,
            Some((&[], [0.0; 3])),
            None,
            None,
            Some(endpoint),
        )
        .map(|_| ())
    }

    /// Pooled acceleration step with the same selected-wall loads as the
    /// serial counterpart.
    pub fn step_pooled_with_endpoint_body_acceleration_and_wall_momentum<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        acceleration: [f64; 3],
        links: &[BoundaryLink3],
        moment_origin: [f64; 3],
    ) -> Result<WallMomentumExchange3, BoundaryStepError3> {
        self.validate_wall_measurement(links, moment_origin)?;
        let endpoint = self.endpoint_body_acceleration(acceleration, gate)?;
        self.step_pooled_coupling_inner(
            pool,
            gate,
            Some((links, moment_origin)),
            None,
            None,
            Some(endpoint),
        )
        .map(|(exchange, _)| exchange)
    }

    fn endpoint_body_acceleration(
        &self,
        value: [f64; 3],
        gate: &CancelGate,
    ) -> Result<EndpointForcing3, BoundaryStepError3> {
        endpoint_checkpoint(gate)?;
        let invalid = |reason| BoundaryStepError3::EndpointForce { reason };
        if !value.into_iter().all(f64::is_finite) {
            return Err(invalid("endpoint body acceleration must be finite"));
        }
        if !self.collision_model.supports_body_force() {
            return Err(invalid(
                "body acceleration requires a forcing-compatible collision model",
            ));
        }
        if Face3::ALL
            .into_iter()
            .any(|face| self.boundaries.face(face).is_open())
        {
            return Err(invalid(
                "body acceleration requires closed or periodic boundaries",
            ));
        }
        let initial_field = self
            .force_field
            .is_none()
            .then(|| vec![[0.0; 3]; self.f[0].len() * TILE_CELLS]);
        Ok(EndpointForcing3::Acceleration {
            value,
            initial_field,
        })
    }

    // Called only after all endpoint cells and the final checkpoint admit.
    // Recompute rho from the now-accepted population order to reuse the force
    // allocation without a second full-size proposal or a retained rho array.
    pub(in super::super) fn publish_endpoint_body_acceleration(
        &mut self,
        acceleration: [f64; 3],
        mut initial_field: Option<Vec<[f64; 3]>>,
    ) {
        let forces = self
            .force_field
            .as_mut()
            .or(initial_field.as_mut())
            .expect("accepted or first-use acceleration force storage exists");
        for (tile, forces) in forces.chunks_exact_mut(TILE_CELLS).enumerate() {
            for (lane, force) in forces.iter_mut().enumerate() {
                *force = if self.solid[tile] & (1u64 << lane) == 0 {
                    let density: f64 = self.f.iter().map(|field| field[tile].0[lane]).sum();
                    acceleration.map(|value| density * value)
                } else {
                    [0.0; 3]
                };
            }
        }
        if let Some(field) = initial_field {
            self.force_field = Some(field);
        }
        self.force = [0.0; 3];
        self.force_field_shape = None;
    }
}
