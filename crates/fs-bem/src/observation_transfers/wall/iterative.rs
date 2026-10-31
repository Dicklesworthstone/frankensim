//! Material diffraction observers through the existing matrix-free boundary.
//! FMM error, GMRES residual, receiver quadrature and FIR convergence remain
//! distinct. A failed source solve cannot produce a partial frequency or bake.
use super::*;
use crate::helmholtz::bc::{IncidentField, IncidentTraces, Monopole, PreparedIterativeBoundary};
use crate::helmholtz::{IterativeError, IterativeOptions, Preconditioner, ResourceReceipt};
use fs_exec::CancelGate;
use fs_solver::krylov::SolveReport;

mod response;
pub use response::{
    BoundaryPortImpedance, IterativeBoundaryResponse, IterativeCoherentDirectivity,
    IterativeDirectivityOptions, IterativeFarFieldSample,
};

/// Explicit spatial/numerical controls; the material and physical sources are
/// separate inputs. The same solver settings apply throughout the FFT ladder.
#[derive(Clone, Debug)]
pub struct IterativeObservationOptions {
    pub highest_omega: f64,
    pub minimum_clearance_m: f64,
    pub receiver_quadrature: Options,
    pub solver: IterativeOptions,
    pub lower_window: Option<LowerBandWindow>,
}

#[derive(Clone, Debug)]
pub enum IterativeObservationError {
    Transfer(TransferError),
    /// Original convergence, FMM-budget or cancellation diagnostics retained.
    Solver(Box<IterativeError>),
}
impl From<TransferError> for IterativeObservationError {
    fn from(error: TransferError) -> Self {
        Self::Transfer(error)
    }
}
impl From<HelmholtzError> for IterativeObservationError {
    fn from(error: HelmholtzError) -> Self {
        Self::Transfer(match error {
            HelmholtzError::Cancelled => TransferError::Cancelled,
            other => TransferError::Helmholtz(other),
        })
    }
}
impl From<IterativeError> for IterativeObservationError {
    fn from(error: IterativeError) -> Self {
        Self::Solver(Box::new(error))
    }
}
impl core::fmt::Display for IterativeObservationError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(formatter, "iterative observation: {self:?}")
    }
}
impl std::error::Error for IterativeObservationError {}

#[derive(Clone, Debug)]
pub struct IterativeSourceReport {
    pub source: SourceId,
    /// Recomputed original-equation residual; not an FMM error bound.
    pub solve: SolveReport,
    /// Original receipt includes shared setup time; do not sum that time as
    /// though each source rebuilt a separate FMM/preconditioner.
    pub resources: ResourceReceipt,
    /// Total-wall passive work for this unit source coordinate [W].
    pub absorbed_power_w: f64,
}

#[derive(Clone, Debug)]
pub struct IterativeObservationSample {
    pub omega_rad_s: f64,
    pub operator_hash: ContentHash,
    /// Source-major, receiver-major p, dp/dx, dp/dy, dp/dz.
    pub transfer: Vec<C64>,
    /// One actual solve per source. Empty for DC or an exactly zero window;
    /// neither case manufactures solver convergence reports.
    pub sources: Vec<IterativeSourceReport>,
}

/// Matrix-free source-to-receiver diffraction. It binds the same material,
/// source and window owners as the dense adapter, with explicit FMM/GMRES
/// settings. Receiver integration retains its original finite work ceiling.
pub struct IterativeBoundaryObservation<'a> {
    observation: LocalBoundaryObservation<'a>,
    solver: IterativeOptions,
}

#[derive(Clone, Debug)]
pub struct IterativeBoundaryObservationBake {
    observation: BoundaryObservationBake,
    solver: IterativeOptions,
}
impl IterativeBoundaryObservationBake {
    pub fn observation(&self) -> &BoundaryObservationBake {
        &self.observation
    }
    pub fn solver_options(&self) -> &IterativeOptions {
        &self.solver
    }
    pub fn observer(&self) -> Result<ObservationFir, TransferError> {
        self.observation.observer()
    }
}

impl<'a> IterativeBoundaryObservation<'a> {
    /// Sources and receivers are admitted once against the immutable closed
    /// surface. No dense factorization or FMM setup occurs in this constructor.
    pub fn for_scattering(
        surface: &'a SpherePanels,
        medium: Medium,
        sources: &[PointSourceSpec],
        receivers: &[ReceiverSpec],
        boundary: LocalBoundary,
        options: IterativeObservationOptions,
        gate: &CancelGate,
    ) -> Result<Self, IterativeObservationError> {
        admit(surface, medium, &boundary, &options, gate)?;
        let owner = ScatteringObservation::new_admitted(
            surface,
            medium,
            sources,
            receivers,
            options.highest_omega,
            options.minimum_clearance_m,
            options.receiver_quadrature,
            Some(gate),
            || gate.is_requested(),
        )?;
        Ok(Self {
            observation: LocalBoundaryObservation {
                excitation: Excitation::Point(owner),
                boundary,
                highest_omega: options.highest_omega,
                lower_window: options.lower_window,
            },
            solver: options.solver,
        })
    }

    /// Original unit-flow/signed structural traces, with exactly the same
    /// disjoint prescribed-drive/absorbing-patch requirement as the dense path.
    #[allow(clippy::too_many_arguments)]
    pub fn for_surface(
        surface: &'a SpherePanels,
        medium: Medium,
        sources: &[SourceId],
        traces: &[Vec<f64>],
        receivers: &[ReceiverSpec],
        boundary: LocalBoundary,
        options: IterativeObservationOptions,
        gate: &CancelGate,
    ) -> Result<Self, IterativeObservationError> {
        admit(surface, medium, &boundary, &options, gate)?;
        let source_entries =
            sources
                .len()
                .checked_mul(surface.areas().len())
                .ok_or(TransferError::Invalid(
                    "iterative source basis shape overflow",
                ))?;
        cap(source_entries, MAX_OBSERVATION_ENTRIES)?;
        let owner = BemObservation::new_admitted_controlled(
            surface,
            medium,
            sources,
            traces,
            receivers,
            options.minimum_clearance_m,
            false,
            Some(gate),
        )?;
        for field in &owner.fields {
            poll(gate)?;
            if field
                .iter()
                .enumerate()
                .any(|(index, &value)| value != C64::ZERO && !boundary.rigid_at(index))
            {
                return Err(TransferError::Invalid(
                    "absorbing/compliant wall overlaps a prescribed source trace",
                )
                .into());
            }
        }
        let geometry =
            Geometry::new_iterative(surface, &owner.points, options.minimum_clearance_m, gate)?;
        Ok(Self {
            observation: LocalBoundaryObservation {
                excitation: Excitation::Surface {
                    owner,
                    geometry,
                    quadrature: options.receiver_quadrature,
                },
                boundary,
                highest_omega: options.highest_omega,
                lower_window: options.lower_window,
            },
            solver: options.solver,
        })
    }

    pub fn sources(&self) -> &[SourceId] {
        self.observation.sources()
    }
    pub fn receivers(&self) -> &[ReceiverSpec] {
        self.observation.receivers()
    }
    pub fn geometric_delays_s(&self) -> Vec<f64> {
        self.observation.geometric_delays_s()
    }
    pub fn lower_band_window(&self) -> Option<LowerBandWindow> {
        self.observation.lower_band_window()
    }
    pub fn solver_options(&self) -> &IterativeOptions {
        &self.solver
    }

    pub fn operator_hash(&self) -> ContentHash {
        let mut hash = DomainHasher::new("fs-bem:iterative-local-boundary-observation:v1");
        hash.update(self.observation.operator_hash().as_bytes());
        hash_options(&mut hash, &self.solver);
        hash.finalize()
    }

    /// Prepare the complete harmonic p/u source basis through this matrix-free
    /// owner. The returned original-equation solve reports remain available
    /// alongside the field bank; algebraic statistics do not certify FMM error.
    /// All subsequent coherence and microphone-array evaluations reuse these
    /// transfers, without rerunning GMRES or constructing a dense BEM matrix.
    pub fn prepare_harmonic(
        &self,
        omega: f64,
        gate: &CancelGate,
    ) -> Result<(PreparedHarmonicObservation, Vec<IterativeSourceReport>), IterativeObservationError>
    {
        poll(gate)?;
        let medium = match &self.observation.excitation {
            Excitation::Point(owner) => owner.medium,
            Excitation::Surface { owner, .. } => owner.medium,
        };
        harmonic::admit_frequency(omega, medium)?;
        let sample = self.physical_frequency_sample(omega, gate)?;
        let field = PreparedHarmonicObservation::from_samples(
            omega,
            medium,
            self.sources(),
            self.receivers(),
            sample.operator_hash,
            sample.transfer,
        )?;
        Ok((field, sample.sources))
    }

    /// Native exp(-i omega t) physical transfer. One prepared Robin FMM and
    /// preconditioner plus one receiver quadrature table serve every source.
    /// No dense fallback, field truncation or unconverged pressure is returned.
    pub fn physical_frequency_sample(
        &self,
        omega: f64,
        gate: &CancelGate,
    ) -> Result<IterativeObservationSample, IterativeObservationError> {
        poll(gate)?;
        self.observation.admit_frequency(omega)?;
        if omega == 0.0 {
            if self.observation.boundary.supported_band_hz().0 > 0.0 {
                return Err(
                    TransferError::Invalid("DC is outside local-wall source support").into(),
                );
            }
            return Ok(self.zero_sample(omega));
        }
        let (_, medium, geometry, quadrature) = self.field_geometry();
        let (_, output) = self.visit_solutions(
            omega,
            gate,
            || {
                Ok((
                    geometry.prepare_velocity_cancellable(
                        omega / medium.sound_speed,
                        medium,
                        quadrature,
                        gate,
                    )?,
                    IterativeObservationSample {
                        omega_rad_s: omega,
                        operator_hash: self.operator_hash(),
                        transfer: Vec::with_capacity(self.observation.channel_count()),
                        sources: Vec::with_capacity(self.sources().len()),
                    },
                ))
            },
            |(rows, output), index, radiation, report| {
                self.append_receiver_field(omega, index, rows, &radiation, &mut output.transfer)?;
                output.sources.push(report);
                Ok(())
            },
        )?;
        self.admit_sample(&output)?;
        Ok(output)
    }

    fn field_geometry(&self) -> (&'a SpherePanels, Medium, &Geometry<'a>, Options) {
        match &self.observation.excitation {
            Excitation::Point(owner) => (
                owner.surface,
                owner.medium,
                &owner.receiver_geometry,
                owner.quadrature,
            ),
            Excitation::Surface {
                owner,
                geometry,
                quadrature,
            } => (owner.surface, owner.medium, geometry, *quadrature),
        }
    }

    // One original solve traversal for streaming observation and retained
    // boundary responses. The consumer is prepared after the boundary owner,
    // preserving material refusal and receiver-admission order.
    fn visit_solutions<T>(
        &self,
        omega: f64,
        gate: &CancelGate,
        prepare: impl FnOnce() -> Result<T, IterativeObservationError>,
        mut receive: impl FnMut(
            &mut T,
            usize,
            helmholtz::RadiationSolution,
            IterativeSourceReport,
        ) -> Result<(), IterativeObservationError>,
    ) -> Result<T, IterativeObservationError> {
        poll(gate)?;
        self.observation.admit_frequency(omega)?;
        if omega == 0.0 {
            return Err(TransferError::Invalid(
                "retained boundary fields require positive frequency",
            )
            .into());
        }
        let beta = self
            .observation
            .boundary
            .normalized_admittance(omega, || gate.is_requested())?;
        let law = Admittance::Local {
            beta: &beta,
            declared_active: false,
        };
        let (surface, medium, _, _) = self.field_geometry();
        let k = omega / medium.sound_speed;
        let prepared = PreparedIterativeBoundary::new(
            surface,
            k,
            medium,
            law,
            helmholtz::Formulation::BurtonMiller,
            &self.solver,
            gate,
        )?;
        let mut output = prepare()?;
        for (index, source) in self.sources().iter().enumerate() {
            poll(gate)?;
            let (radiation, solve, resources, absorbed_power_w) = match &self.observation.excitation
            {
                Excitation::Point(owner) => {
                    let monopoles = [Monopole {
                        position_m: owner.point_sources[index].position_m,
                        volume_velocity_m3_s: C64::ONE,
                    }];
                    let incident = IncidentTraces::from_admitted_field(
                        surface,
                        k,
                        medium,
                        IncidentField::Monopoles(&monopoles),
                        Some(gate),
                    )?;
                    let (response, solve, resources) = prepared.scatter(&incident, gate)?;
                    (
                        response.scattered,
                        solve,
                        resources,
                        response.absorbed_power_w,
                    )
                }
                Excitation::Surface { owner, .. } => {
                    let (response, solve, resources) =
                        prepared.solve(&owner.fields[index], gate)?;
                    (
                        response.radiation,
                        solve,
                        resources,
                        response.absorbed_power_w,
                    )
                }
            };
            receive(
                &mut output,
                index,
                radiation,
                IterativeSourceReport {
                    source: source.clone(),
                    solve,
                    resources,
                    absorbed_power_w,
                },
            )?;
        }
        poll(gate)?;
        Ok(output)
    }

    fn append_receiver_field(
        &self,
        omega: f64,
        index: usize,
        rows: &crate::near_field::Prepared<'a>,
        radiation: &helmholtz::RadiationSolution,
        output: &mut Vec<C64>,
    ) -> Result<(), IterativeObservationError> {
        let medium = self.field_geometry().1;
        let k = omega / medium.sound_speed;
        let received = rows.evaluate_velocity(radiation)?;
        let momentum = C64::new(0.0, omega * medium.density);
        for (receiver, (&pressure, velocity)) in self.receivers().iter().zip(
            received
                .scalar
                .pressure
                .iter()
                .zip(&received.particle_velocity_m_s),
        ) {
            let (pressure, velocity) = match &self.observation.excitation {
                Excitation::Point(owner) => {
                    let monopoles = [Monopole {
                        position_m: owner.point_sources[index].position_m,
                        volume_velocity_m3_s: C64::ONE,
                    }];
                    let direct = IncidentField::Monopoles(&monopoles).pressure_velocity_at(
                        k,
                        medium,
                        receiver.position_m,
                    )?;
                    (
                        pressure + direct.pressure_pa,
                        std::array::from_fn(|axis| velocity[axis] + direct.velocity_m_s[axis]),
                    )
                }
                Excitation::Surface { .. } => (pressure, *velocity),
            };
            output.push(pressure);
            output.extend(velocity.iter().map(|&component| momentum * component));
        }
        Ok(())
    }

    fn admit_sample(
        &self,
        output: &IterativeObservationSample,
    ) -> Result<(), IterativeObservationError> {
        if output.transfer.len() != self.observation.channel_count()
            || output.transfer.iter().any(|&value| !finite(value))
        {
            return Err(
                TransferError::Invalid("nonfinite or incomplete iterative observation").into(),
            );
        }
        Ok(())
    }

    pub fn band_limited_frequency_sample(
        &self,
        omega: f64,
        gate: &CancelGate,
    ) -> Result<IterativeObservationSample, IterativeObservationError> {
        poll(gate)?;
        self.observation.admit_frequency(omega)?;
        let weight = self
            .lower_band_window()
            .map_or(1.0, |window| window.weight(omega / core::f64::consts::TAU));
        if weight == 0.0 {
            return Ok(self.zero_sample(omega));
        }
        let mut sample = self.physical_frequency_sample(omega, gate)?;
        if weight != 1.0 {
            for value in &mut sample.transfer {
                *value = value.scale(weight);
            }
        }
        Ok(sample)
    }

    /// Every nonzero sampled source must pass the original GMRES convergence
    /// gate before the common FIR owner can inspect temporal convergence.
    /// FMM approximation and spatial error remain producer uncertainty.
    pub fn bake(
        &self,
        declaration: &TransferDeclaration,
        gate: &CancelGate,
    ) -> Result<IterativeBoundaryObservationBake, IterativeObservationError> {
        poll(gate)?;
        self.observation
            .admit_declaration(declaration, self.operator_hash())?;
        let transfer = bake_observation_transfers_with_errors(
            |omega| {
                self.band_limited_frequency_sample(omega, gate)
                    .map(|sample| sample.transfer)
            },
            declaration,
            self.lower_band_window()
                .map_or(0.0, |window| window.pass_hz),
            || gate.is_requested(),
        )?;
        Ok(IterativeBoundaryObservationBake {
            observation: BoundaryObservationBake {
                transfer,
                lower_window: self.lower_band_window(),
                boundary_hash: self.observation.boundary.content_hash(),
            },
            solver: self.solver.clone(),
        })
    }

    fn zero_sample(&self, omega: f64) -> IterativeObservationSample {
        IterativeObservationSample {
            omega_rad_s: omega,
            operator_hash: self.operator_hash(),
            transfer: vec![C64::ZERO; self.observation.channel_count()],
            sources: Vec::new(),
        }
    }
}

fn poll(gate: &CancelGate) -> Result<(), IterativeObservationError> {
    if gate.is_requested() {
        Err(TransferError::Cancelled.into())
    } else {
        Ok(())
    }
}

fn admit(
    surface: &SpherePanels,
    medium: Medium,
    boundary: &LocalBoundary,
    options: &IterativeObservationOptions,
    gate: &CancelGate,
) -> Result<(), IterativeObservationError> {
    poll(gate)?;
    options.receiver_quadrature.validate()?;
    admit_binding_for_route(
        surface,
        medium,
        options.highest_omega,
        boundary,
        options.lower_window,
        true,
    )?;
    crate::helmholtz::iterative::admit(
        surface,
        options.highest_omega / medium.sound_speed,
        medium,
        &options.solver,
    )?;
    if let Some(directional) = &options.solver.directional {
        directional.validate().map_err(IterativeError::from)?;
    }
    poll(gate)
}

fn hash_options(hash: &mut DomainHasher, options: &IterativeOptions) {
    hash.update(&options.tol.to_bits().to_le_bytes());
    for value in [
        options.restart,
        options.max_iters,
        options.fmm_order,
        options.leaf_capacity,
        options.max_near_entries,
        options.near_field_cache_bytes,
        options.max_translation_work,
        options.max_coefficient_slots,
        options.workers,
        options.m2l_cache_bytes,
    ] {
        hash.update(&(value as u64).to_le_bytes());
    }
    match options.precond {
        Preconditioner::None => hash.update(&[0]),
        Preconditioner::NearFieldBlockJacobi { max_block_panels } => {
            hash.update(&[1]);
            hash.update(&(max_block_panels as u64).to_le_bytes());
        }
    }
    match options.panel_rule {
        helmholtz::PanelRule::Centroid => hash.update(&[0]),
        helmholtz::PanelRule::TriangleNear { max_pairs } => {
            hash.update(&[1]);
            hash.update(&(max_pairs as u64).to_le_bytes());
        }
        helmholtz::PanelRule::Triangle { max_pairs } => {
            hash.update(&[2]);
            hash.update(&(max_pairs as u64).to_le_bytes());
        }
    }
    match &options.directional {
        None => hash.update(&[0]),
        Some(direction) => {
            hash.update(&[1]);
            for value in [
                direction.tolerance,
                direction.low_freq_threshold,
                direction.compression_tolerance,
            ] {
                hash.update(&value.to_bits().to_le_bytes());
            }
            for value in [
                direction.max_cones_per_level,
                direction.order,
                direction.cone_resolution,
                direction.leaf_capacity,
                direction.max_rank,
                direction.max_coefficient_slots,
                direction.max_setup_work,
                direction.w5a_budget.max_translation_work,
                direction.w5a_budget.max_coefficient_slots,
                direction.max_near_pairs,
                direction.setup_workers,
                direction.max_translation_work,
            ] {
                hash.update(&(value as u64).to_le_bytes());
            }
        }
    }
}

#[cfg(test)]
mod tests;
