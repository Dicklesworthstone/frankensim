//! Density-dependent acceleration in the original Guo source history.

use super::*;

mod affine_density;

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
    /// Closed, periodic and native velocity/pressure boundaries are supported.
    /// Open-face endpoints solve their density-dependent source together with
    /// the existing regularized reconstruction. Finite-inertia immersed-body
    /// feedback needs its own coupled source closure. Low-Mach requirements
    /// remain unchanged; pressure, viscosity and wall reactions are not added
    /// again as body forces.
    pub fn rebind_body_acceleration(
        &mut self,
        acceleration: [f64; 3],
        gate: &CancelGate,
    ) -> Result<(), BoundaryForcingError3> {
        // Reuse the atomic zero-motion preparation, without broadening the
        // public rotating-frame owner's boundary-admission contract.
        RigidFrame3::prepare_body_acceleration(self, acceleration, gate)
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
    /// implicit velocity-dependent force. Closed, periodic and native open
    /// boundaries use the original forcing-compatible collision laws. At a
    /// wall/open seam, only open-owned populations change, so the final density
    /// and its force are solved together. This preserves the existing partial
    /// reconstruction convention: full-face targets need not hold at seams.
    /// Reservoir transport remains separate from the body-force impulse.
    ///
    /// Endpoint admission is fused into streaming and open reconstruction;
    /// the accepted force allocation is reused after first use. Refusal or
    /// cancellation changes neither accepted populations nor forces. No
    /// additional grid clone or second collision is introduced.
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
        let initial_field = self
            .force_field
            .is_none()
            .then(|| vec![[0.0; 3]; self.f[0].len() * TILE_CELLS]);
        Ok(EndpointForcing3::Acceleration {
            value,
            initial_field,
        })
    }

    // The native source frame is linear in force at fixed reconstruction
    // velocity. A full face preserves its template density. At a wall seam,
    // replacing only open-owned rows gives rho = B + C*rho for F = rho*a.
    pub(in super::super) fn open_body_acceleration_populations(
        &self,
        fields: &[Vec<Tile>; Q3],
        cell: [usize; 3],
        rho: f64,
        velocity: [f64; 3],
        stress: [[f64; 3]; 3],
        acceleration: [f64; 3],
    ) -> Result<[f64; Q3], BoundaryStepError3> {
        let (tile, lane) = self.addr(cell[0], cell[1], cell[2]);
        let wall_mask = self.link_masks[tile].0[lane];
        let open_mask = self.open_link_masks[tile].0[lane];
        let before = core::array::from_fn(|q| fields[q][tile].0[lane]);
        let template = regularized_populations(rho, velocity, stress);
        let density = if wall_mask == 0 {
            // The native Guo reframe retains the full template's zeroth moment.
            template.iter().sum::<f64>()
        } else {
            let base = reconstruct_open_populations(
                template, velocity, [0.0; 3], before, wall_mask, open_mask,
            );
            let mut response = [0.0; Q3];
            crate::thermal::forcing::reframe3(&mut response, velocity, [0.0; 3], acceleration);
            for (q, value) in response.iter_mut().enumerate() {
                if open_mask & (1u32 << q) == 0 {
                    *value = 0.0;
                }
            }
            affine_density::solve(&base, &response)
                .ok_or(BoundaryStepError3::Unphysical { cell })?
        };
        let force = acceleration.map(|value| density * value);
        let reconstructed =
            reconstruct_open_populations(template, velocity, force, before, wall_mask, open_mask);
        // Admit against the ACTUAL reconstructed mass used at publication.
        // A bad closure, overflowing source or physical velocity stays private.
        let actual_force = EndpointSource3::Acceleration(acceleration)
            .force(&reconstructed, tile, lane)
            .ok_or(BoundaryStepError3::Unphysical { cell })?;
        if !affine_density::matches(&reconstructed, density)
            || !finite_forced_populations(&reconstructed, actual_force)
        {
            return Err(BoundaryStepError3::Unphysical { cell });
        }
        Ok(reconstructed)
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
