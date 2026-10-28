//! Multiple moving surfaces and one neck share ONE compression coordinate.
//! Material/geometry owners supply signed volume shapes; the original Helmholtz
//! owner supplies storage and loss, and CoupledModalSystem owns every step.
//! Pressure diagnostics are read-only; they neither drive a second integrator
//! nor add another energy account. All motion is incremental about a fixed
//! reference equilibrium. This is not nonlinear gas dynamics or moving geometry.
use super::*;
use crate::modal_acoustic_time::ModalAcousticState;
use crate::render::plate::impact::radiation::Model as RadiationModel;
use crate::render::schedule::force::coupled::{
    CoupledModalFrame, CoupledModalSystem, ModalColumnConnection, ModalCouplingConfig,
};
use fs_exec::CancelGate;

/// One material/geometry-derived modal surface, with its ORIGINAL signed basis.
/// Source receipts stay with the upstream owner (e.g. MembraneVibration), not
/// inside the modal time integrator. Positive shapes point out of this cavity.
pub struct CavitySurface {
    /// Original clock, states, frequencies, loss and observation transfers.
    pub model: ModalAcousticTimeModel,
    /// Integrated outward normal-displacement shapes [m2/sqrt(kg)]. These are
    /// not point displacements or magnitudes; opposite contributions cancel.
    pub outward_volume_shapes: Vec<f64>,
}

/// Physical observations at an accepted endpoint, not midpoint reaction forces.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SharedCavityObservation {
    /// Increase in effective cavity volume [m3], including outward neck travel.
    pub volume_change_m3: f64,
    /// Signed small-signal change divided by fixed reference cavity volume.
    pub fractional_volume_change: f64,
    /// Uniform INTERNAL pressure minus the equilibrium pressure [Pa].
    pub pressure_increment_pa: f64,
    /// Sum of moving-surface outward volume rates [m3/s].
    pub surface_volume_rate_m3_s: f64,
    /// Outward neck flow [m3/s], suitable for the original acoustic receiver.
    pub neck_outflow_m3_s: f64,
    /// Compression storage already included in network connection energy [J].
    pub compression_energy_j: f64,
    /// Neck kinetic energy already included in network modal energy [J].
    pub neck_kinetic_energy_j: f64,
}

/// Original modal components coupled through a single material-gas volume.
/// Surface order is retained; the free neck is appended last. No frequency
/// tuning, per-surface gas springs or one-way pressure playback occurs here.
pub struct SharedCavity {
    source: HelmholtzModalCavity,
    system: CoupledModalSystem,
    column: Vec<f64>,
    surface_modes: usize,
    external: Vec<f64>,
    ports: crate::port_adoption::ModulePorts,
    old_q: Vec<f64>,
    radiation_neck_port: Option<usize>,
}
impl SharedCavity {
    /// Keep every surface's states and connect their SUMMED volume displacement
    /// to one neck. All components must have the supplied mechanical clock.
    /// The original pair lowering supplies each shape and constitutive number;
    /// the existing whole-column solver retains inter-surface cross terms.
    pub fn new(
        source: HelmholtzModalCavity,
        surfaces: Vec<CavitySurface>,
        sample_rate_hz: u32,
        neck_initial: ModalAcousticState,
        neck_budget: ModalAcousticTimeBudget,
        config: ModalCouplingConfig,
        gate: &CancelGate,
    ) -> Result<Self, ModalCouplingError> {
        if gate.is_requested() { return Err(ModalCouplingError::Cancelled); }
        if surfaces.is_empty() || surfaces.len()>64 {
            return Err(network_invalid("shared cavity requires 1..=64 moving surfaces"));
        }
        let surface_modes = surfaces.iter().try_fold(0usize, |n,s| n.checked_add(s.model.modes().len()))
            .ok_or_else(||network_invalid("shared cavity mode count overflow"))?;
        let count = surface_modes.checked_add(1).ok_or_else(||network_invalid("shared cavity mode count overflow"))?;
        if count>4096 || count>config.max_modes {
            return Err(network_invalid("shared cavity exceeds the original mode budget"));
        }
        let neck_index=surfaces.len();
        let mut charts = Vec::with_capacity(count);
        let mut models=Vec::with_capacity(neck_index+1);
        let mut column=vec![0.0;count];
        let mut offset=0;
        let mut loss=None;
        for (index,surface) in surfaces.into_iter().enumerate() {
            if gate.is_requested() { return Err(ModalCouplingError::Cancelled); }
            if surface.outward_volume_shapes.len()!=surface.model.modes().len() {
                return Err(network_invalid("surface volume projection must match its original modes"));
            }
            let pair=source.connections(index,neck_index,&surface.outward_volume_shapes)?;
            charts.extend(surface.outward_volume_shapes.iter().map(|&projection|
                fs_ports::PortCoordinateChart::MassNormalized {
                    reference_mass_kg: 1.0, volume_projection_m2_per_sqrt_kg: Some(projection),
                }));
            let spring=&pair[0];
            let end=offset+surface.model.modes().len();
            column[offset..end].copy_from_slice(&spring.left.shapes);
            if index==0 {
                // A connection is left MINUS right; the neck's negative shape
                // makes outward surface and neck displacements ADD to volume.
                column[surface_modes]=-spring.right.shapes[0];
                if let Some(damper)=pair.get(1) {
                    let mut c=vec![0.0;count];
                    c[surface_modes]=damper.left.shapes[0]-damper.right.shapes[0];
                    loss=Some(ModalColumnConnection {column:c,stiffness_n_m:damper.stiffness_n_m,
                        damping_n_s_m:damper.damping_n_s_m,rest_extension_m:damper.rest_extension_m});
                }
            }
            offset=end;models.push(surface.model);
        }
        let mut neck=source.neck_component(sample_rate_hz,neck_budget)?;
        neck.restore_states(&[neck_initial]).map_err(|source|ModalCouplingError::Component {component:neck_index,source})?;
        models.push(neck);
        let mut links=vec![ModalColumnConnection {column:column.clone(),
            stiffness_n_m:source.compression_stiffness_n_m(),damping_n_s_m:0.0,rest_extension_m:0.0}];
        if let Some(loss)=loss {links.push(loss);}
        let system=CoupledModalSystem::new_with_columns(models,links,config,gate)?;
        let names: Vec<_> = (0..surface_modes)
            .map(|i| format!("surface/generalized/{i}"))
            .collect();
        let mut kinds: Vec<_> = names
            .iter()
            .map(|s| (s.as_str(), fs_ports::PortKind::MechanicalForceVelocity))
            .collect();
        kinds.push(("neck/exterior", fs_ports::PortKind::FluidPressureFlux));
        let mut ports = crate::port_adoption::ModulePorts::new("shared-cavity", kinds);
        charts.push(fs_ports::PortCoordinateChart::CanonicalSi);
        ports.charts = charts;
        // The retained signed-volume shapes pull each surface into pressure /
        // volume-flow coordinates. One extra endpoint owns compression storage.
        ports.junctions.push(crate::port_adoption::junction(
            "shared-cavity",
            "compression-volume",
            neck_index + 2,
            fs_ports::PortKind::FluidPressureFlux,
            false,
        ));
        let result = Self {
            source,
            system,
            column,
            surface_modes,
            external: vec![0.0; count],
            ports,
            old_q: vec![0.0; count],
            radiation_neck_port: None,
        };
        result.observe()?;
        if gate.is_requested() {return Err(ModalCouplingError::Cancelled);}
        Ok(result)
    }
    /// Original coherent gas, geometry, neck inertia, compliance and loss.
    pub fn source(&self)->&HelmholtzModalCavity {&self.source}
    /// Original accepted bank; no mutable access that could bypass coupling.
    pub fn network(&self)->&CoupledModalSystem {&self.system}
    /// Number of external surface forces, excluding the appended neck.
    pub fn surface_mode_count(&self)->usize {self.surface_modes}

    /// Attach a reciprocal exterior load to the SAME cavity/surface time solve.
    ///
    /// `source_modes` selects original flattened surface coordinates and must
    /// contain the neck coordinate, `surface_mode_count()`, exactly once. Body
    /// columns use their existing mass-normalized velocities. The neck column
    /// uses physical OUTWARD volume flow [m3/s], so its conjugate effort is
    /// pressure [Pa]. Its coupling is pulled back by Q = (S/sqrt(M)) q_dot;
    /// this preserves both cross-source reactions and instantaneous power.
    ///
    /// The source must use `ExteriorMouth`: the attached impedance owns the
    /// exterior added mass and radiation loss, with no duplicated end correction.
    /// This is a cold attachment. The existing exact propagator, cancellation,
    /// pole/work limits and energy accounting own all acoustic states. A zero
    /// impedance leaves the original trajectory unchanged. Geometry, bandwidth
    /// and fit accuracy remain obligations of the supplied model's producer.
    /// The pressure supplied to `step` is then additional INCIDENT pressure,
    /// never the already-computed reaction returned by `neck_exterior_pressure_pa`.
    pub fn with_exterior_radiation(mut self, model: &RadiationModel,
        source_modes: &[usize], maximum_poles: usize, gate: &CancelGate)
        -> Result<Self, ModalCouplingError>
    {
        if gate.is_requested() { return Err(ModalCouplingError::Cancelled); }
        model.validate().map_err(|_| network_invalid("invalid complete exterior cavity load"))?;
        if self.source.loss() != CavityNeckLoss::ExteriorMouth
            || self.radiation_neck_port.is_some()
            || source_modes.len() != model.ports
            || source_modes.iter().filter(|&&i| i == self.surface_modes).count() != 1 {
            return Err(network_invalid("exterior cavity load needs an internal-only neck and exactly one physical neck port"));
        }
        let neck_port = source_modes.iter().position(|&i| i == self.surface_modes)
            .expect("one admitted neck port");
        let projection = self.source.neck_volume_shape_m2_per_sqrt_kg();
        let mut pulled_back = model.clone();
        for pole in &mut pulled_back.poles {
            if gate.is_requested() { return Err(ModalCouplingError::Cancelled); }
            let original = pole.coupling[neck_port];
            let scaled = original * projection;
            if !scaled.is_finite() || (original != 0.0 && scaled == 0.0) {
                return Err(network_invalid("exterior neck pressure/flow projection is not representable"));
            }
            pole.coupling[neck_port] = scaled;
        }
        self.system = self.system.with_radiation_load(&pulled_back, source_modes, maximum_poles, gate)?;
        self.radiation_neck_port = Some(neck_port);
        Ok(self)
    }

    /// Accepted exterior reaction pressure [Pa] opposing OUTWARD neck flow.
    /// Includes the signed contribution of every retained surface source. This
    /// is a pure observation of the joint acoustic state, not a microphone or
    /// another force to feed back into `step`; its work is already internal.
    pub fn neck_exterior_pressure_pa(&self) -> Result<f64, ModalCouplingError> {
        let Some(port) = self.radiation_neck_port else { return Ok(0.0); };
        let Some(model) = self.system.radiation_model() else { return Ok(0.0); };
        let projection = self.source.neck_volume_shape_m2_per_sqrt_kg();
        let mut pressure = 0.0;
        for (pole, state) in model.poles.iter().zip(self.system.radiation_state().chunks_exact(2)) {
            pressure += (pole.coupling[port] / projection) * state[1];
        }
        if !pressure.is_finite() { return Err(network_invalid("exterior neck reaction pressure overflow")); }
        Ok(pressure)
    }

    /// One simultaneous accepted step. Surface forces are flattened in input
    /// surface/mode order [N/sqrt(kg)]. Positive EXTERIOR neck pressure [Pa]
    /// acts inward. It is not another internal-pressure drive: that reaction
    /// already comes from the collective compression connection.
    /// No allocation occurs here. Rejection leaves all accepted states/clock
    /// unchanged; scratch forces may be overwritten and are recomputed on retry.
    pub fn step(
        &mut self,
        surface_force: &[f64],
        exterior_neck_pressure_pa: f64,
        gate: &CancelGate,
    ) -> Result<&CoupledModalFrame, ModalCouplingError> {
        if gate.is_requested() {return Err(ModalCouplingError::Cancelled);}
        if surface_force.len()!=self.surface_modes || surface_force.iter().any(|v|!v.is_finite())
            || !exterior_neck_pressure_pa.is_finite() {
            return Err(network_invalid("finite surface force per original mode and exterior neck pressure required"));
        }
        let neck_force=-self.source.neck_volume_shape_m2_per_sqrt_kg()*exterior_neck_pressure_pa;
        if !neck_force.is_finite() {return Err(network_invalid("exterior neck pressure force overflow"));}
        self.external[..self.surface_modes].copy_from_slice(surface_force);
        self.external[self.surface_modes]=neck_force;
        for (old, state) in self
            .old_q
            .iter_mut()
            .zip(self.system.components().iter().flat_map(|m| m.states()))
        {
            *old = state.displacement_m_sqrt_kg;
        }
        let before = self.system.total_energy_j()?;
        let (sample, after, loss, tolerance) = {
            let frame = self.system.step_under_gate(&self.external, gate)?;
            (
                frame.sample,
                frame.total_energy_j(),
                frame.total_dissipation_j(),
                frame.energy_tolerance_j,
            )
        };
        let old_q = &self.old_q;
        let external = &self.external;
        self.ports.accept(
            crate::port_adoption::StepAccount {
                tick: sample,
                duration_s: self.system.sample_period_s(),
                storage_change_j: Some(after - before),
                loss_j: Some(loss),
                included_work_j: 0.0,
                entropy_change_j_k: 0.0,
                tolerance_j: Some(tolerance),
                unavailable: None,
            },
            self.system
                .components()
                .iter()
                .flat_map(|m| m.states())
                .enumerate()
                .map(|(i, state)| {
                    crate::port_adoption::PortTransfer::work(
                        external[i] * (state.displacement_m_sqrt_kg - old_q[i]),
                        true,
                    )
                }),
        );
        Ok(self.system.last_frame().expect("step accepted"))
    }

    /// Pure observation of the already accepted network. An observation error
    /// does not undo a prior step; it reports unrepresentable diagnostic units.
    /// It never advances, deposits heat, claims radiation or counts energy twice.
    pub fn observe(&self)->Result<SharedCavityObservation,ModalCouplingError> {
        let mut extension=0.0;let mut surface_rate=0.0;let mut index=0;
        for model in self.system.components() {for state in model.states() {
            extension+=self.column[index]*state.displacement_m_sqrt_kg;
            if index<self.surface_modes {surface_rate+=self.column[index]*state.velocity_m_sqrt_kg_per_s;}
            index+=1;
        }}
        let neck=self.system.components().last().expect("admitted neck").states()[0];
        let volume=extension*self.source.neck_area_m2();
        let pressure=-(self.source.compression_stiffness_n_m()*extension)/self.source.neck_area_m2();
        let observation=SharedCavityObservation {
            volume_change_m3:volume,fractional_volume_change:volume/self.source.geometry().volume_m3,
            pressure_increment_pa:pressure,surface_volume_rate_m3_s:surface_rate*self.source.neck_area_m2(),
            neck_outflow_m3_s:self.source.neck_volume_shape_m2_per_sqrt_kg()*neck.velocity_m_sqrt_kg_per_s,
            compression_energy_j:0.5*self.source.compression_stiffness_n_m()*extension*extension,
            neck_kinetic_energy_j:0.5*neck.velocity_m_sqrt_kg_per_s*neck.velocity_m_sqrt_kg_per_s,
        };
        if [observation.volume_change_m3,observation.fractional_volume_change,observation.pressure_increment_pa,
            observation.surface_volume_rate_m3_s,observation.neck_outflow_m3_s,observation.compression_energy_j,
            observation.neck_kinetic_energy_j].iter().any(|v|!v.is_finite()) {
            return Err(network_invalid("shared cavity physical observation overflows"));
        }
        Ok(observation)
    }
}

#[cfg(test)]
mod tests;

crate::port_adoption::impl_ported_module!(SharedCavity, ports);
