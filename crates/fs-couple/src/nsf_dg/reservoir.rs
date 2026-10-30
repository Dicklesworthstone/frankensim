//! Same-interval reaction of the resolved escape aperture and all bore outlets.
//! One original joint radiation owner receives independent bore and mouth flows.
//! Spatial/temporal acoustic projection errors remain explicit work defects;
//! streamed mass and conservative total energy are never reclassified.

use super::*;
use crate::bernoulli_aperture::radiation::environment::multiport::bake::mixed::DrivenMouthSurface;

/// Explicit placement, radiation quadrature and pressure-control search domain.
/// The original reservoir tiles must form one axis-aligned rectangle. Its size
/// is derived from those cells, never independently authored or rescaled.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NsfDgReservoirRadiation {
    /// Grid tangential edge zero in the existing common rigid baffle [m].
    /// Coordinates follow the reservoir observer's cyclic axes and low-face sign.
    pub origin_m: [f64; 2],
    /// Original rectangular radiation quadrature; independent of fluid meshing.
    pub cells: [usize; 2],
    pub terminal_radial_rings: usize,
    pub terminal_angular_points: usize,
    /// Ghost pressure offsets from the original source reference [Pa].
    /// This is an admissible solve domain, not an amplitude clamp.
    pub pressure_bracket_pa: [f64; 2],
}

pub(super) struct ReservoirRadiation {
    aperture: usize,
    input: NsfDgReservoirRadiation,
    temperature_k: f64,
}

/// Accepted acoustic exchange at the actual resolved transverse opening.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NsfDgReservoirFrame {
    /// Complete accepted fluid interval, including all original NSF substeps.
    /// Work uses the cached conservative FV predictor flux, not endpoint traces.
    pub resolved_exchange: fs_lbm::nsf3::NsfReservoirExchange,
    pub ghost_pressure_offset_pa: f64,
    pub mean_acoustic_pressure_pa: f64,
    pub outward_flow_m3_s: f64,
    /// Original joint exterior's mean acoustic pressure [Pa].
    pub radiation_pressure_pa: f64,
    /// Work delivered to its original mouth coordinate [J].
    pub radiation_work_j: f64,
    pub pressure_residual_pa: f64,
    /// Original resolved work minus mean pressure times held flow [J].
    pub projection_covariance_j: f64,
    /// Original resolved acoustic work minus exterior mouth work [J].
    pub acoustic_work_residual_j: f64,
    /// Streamed mass minus reference-density swept continuum volume [kg].
    pub mass_projection_difference_kg: f64,
}

impl<L: Lattice3<27>> NsfDg<L> {
    /// Bind the original transverse opening to the SAME computed self/mutual
    /// radiation matrix as the bore terminals. Source and bore flow remain
    /// independent. Both efforts are solved before any physical state advances.
    ///
    /// This is a uniform acoustic-pressure projection on a prescribed-ghost
    /// boundary, with explicitly checked work error. It does not make that
    /// original boundary transparent or add mean gas transport into the bore.
    /// Arbitrary mouth shapes/rotations and a free-field body are not inferred.
    pub fn with_reservoir_radiation(
        self, input: NsfDgReservoirRadiation, gate: &CancelGate,
    ) -> Result<Self, NsfDgError> {
        if self.zone.reservoir_aperture_count() != 1 {
            return Err(NsfDgError::Input("select the passive reservoir explicitly when multiple apertures exist"));
        }
        self.with_reservoir_radiation_at(0, input, gate)
    }

    /// Select the original escape owner independently of a prescribed breath inlet.
    pub fn with_reservoir_radiation_at(
        mut self, aperture: usize, input: NsfDgReservoirRadiation, gate: &CancelGate,
    ) -> Result<Self, NsfDgError> {
        poll(gate)?;
        if self.steps != 0 || self.reservoir_radiation.is_some()
            || self.finite_reservoir_aperture() == Some(aperture)
            || !input.origin_m.iter().chain(input.pressure_bracket_pa.iter()).all(|x| x.is_finite())
            || input.pressure_bracket_pa[0] >= input.pressure_bracket_pa[1]
            || input.cells.iter().any(|n| !(1..=512).contains(n))
        {
            return Err(NsfDgError::Input("resolved reservoir radiation needs one cold finite placement, quadrature and pressure bracket"));
        }
        let opening = self.zone.reservoir_aperture_at(aperture).ok_or(NsfDgError::Input(
            "resolved radiation requires the original transverse reservoir"))?;
        let record = self.load.radiation_record().ok_or(NsfDgError::Input(
            "resolved mouth reaction requires the bore's original computed exterior"))?;
        let gas = record.gas_record().state();
        let near = |a: f64, b: f64| (a-b).abs() <= 32.0*f64::EPSILON*a.abs().max(b.abs());
        if gas != &self.reference
            || self.zone.reservoir_reference_pressure_pa().to_bits() != self.reference.pressure.to_bits()
            || !near(opening.state.density_kg_m3, gas.density)
            || !near(opening.state.temperature_k, gas.temperature)
            || opening.state.velocity_m_s != [0.0; 3]
        {
            return Err(NsfDgError::Input("joint reservoir reaction needs one quiescent source/exterior reference gas"));
        }
        let temperature_k = opening.state.temperature_k;
        let a = (opening.axis + 1) % 3;
        let b = (opening.axis + 2) % 3;
        let mut low = [usize::MAX; 2];
        let mut high = [0usize; 2];
        for cell in &opening.cells {
            for (d, axis) in [a, b].into_iter().enumerate() {
                low[d] = low[d].min(cell[axis]);
                high[d] = high[d].max(cell[axis]);
            }
        }
        let width = high[0].checked_sub(low[0]).and_then(|n| n.checked_add(1));
        let height = high[1].checked_sub(low[1]).and_then(|n| n.checked_add(1));
        if width.zip(height).and_then(|(x,y)| x.checked_mul(y)) != Some(opening.cells.len()) {
            return Err(NsfDgError::Input("joint radiation requires the complete original rectangular reservoir tile set"));
        }
        // The original aperture constructor already enforces unique cells on
        // one face. Count equals its bounding rectangle only for a complete set.
        let dx = self.zone.law().dx_m;
        let sign = if opening.high { 1.0 } else { -1.0 };
        let surface = DrivenMouthSurface {
            center_m: [
                input.origin_m[0] + 0.5*(low[0] as f64 + high[0] as f64 + 1.0)*dx,
                input.origin_m[1] + sign*0.5*(low[1] as f64 + high[1] as f64 + 1.0)*dx,
            ],
            size_m: [width.expect("checked rectangle") as f64*dx,
                height.expect("checked rectangle") as f64*dx],
            cells: input.cells,
            terminal_radial_rings: input.terminal_radial_rings,
            terminal_angular_points: input.terminal_angular_points,
        };
        for offset in input.pressure_bracket_pa {
            self.zone.reservoir_aperture_state_for(aperture, self.reference.pressure + offset, temperature_k)
                .map_err(NsfDgError::Flow)?;
        }
        // Reuse the original source's integer-Hz bake ceiling. Its runtime
        // retains the exact physical dt, including a fractional reciprocal.
        let rate = 2.0*record.band_hz()[1];
        if !rate.is_finite() || rate.fract() != 0.0 || !(1.0..=f64::from(u32::MAX)).contains(&rate) {
            return Err(NsfDgError::Input("joint radiation lost its original source bake clock"));
        }
        self.load = self.load.with_independent_mouth(surface, rate as u32, gate)
            .map_err(|error| if gate.is_requested() { NsfDgError::Cancelled }
                else { NsfDgError::Acoustic(AcousticRealizeError::Nonlinear(error.to_string())) })?;
        poll(gate)?;
        self.reservoir_radiation = Some(ReservoirRadiation { aperture, input, temperature_k });
        Ok(self)
    }

    pub fn reservoir_radiation(&self) -> Option<&NsfDgReservoirRadiation> {
        self.reservoir_radiation.as_ref().map(|owner| &owner.input)
    }

    /// Stable original owner of the passive escape opening.
    pub fn reservoir_radiation_aperture(&self) -> Option<usize> {
        self.reservoir_radiation.as_ref().map(|owner| owner.aperture)
    }

    pub(super) fn step_with_reservoir_drive(
        &mut self, drive: NsfDgDrive, gate: &CancelGate,
    ) -> Result<NsfDgFrame, NsfDgError> {
        poll(gate)?;
        self.admit_finite_reservoir_drive(drive)?;
        let (external, rise, breath) = match drive {
            NsfDgDrive::Boundary(value) => (value, 0.0, None),
            NsfDgDrive::UpstreamPressureRise(value) =>
                (NsfBoundaryDrive::IncomingCharacteristic(0.0), value, None),
            NsfDgDrive::Reservoir { .. } => return Err(NsfDgError::Input(
                "one reservoir opening cannot own prescribed breath pressure and passive exterior pressure simultaneously")),
            NsfDgDrive::ReservoirAperture { aperture, pressure_pa, temperature_k } => {
                let state = self.zone.reservoir_aperture_state_for(aperture, pressure_pa, temperature_k)
                    .map_err(NsfDgError::Flow)?;
                (NsfBoundaryDrive::IncomingCharacteristic(0.0), 0.0, Some((aperture, state)))
            }
            NsfDgDrive::FiniteReservoir { aperture }
            | NsfDgDrive::FiniteReservoirActuation { aperture, .. } => {
                let state = self.finite_reservoir_state(aperture)?;
                (NsfBoundaryDrive::IncomingCharacteristic(0.0), 0.0, Some((aperture, state)))
            }
        };
        let misplaced = match drive {
            NsfDgDrive::Boundary(_) => self.upstream.is_some(),
            NsfDgDrive::UpstreamPressureRise(_) => self.upstream.is_none(),
            NsfDgDrive::Reservoir { .. } | NsfDgDrive::ReservoirAperture { .. }
            | NsfDgDrive::FiniteReservoir { .. }
            | NsfDgDrive::FiniteReservoirActuation { .. } => false,
        };
        if !drive.is_finite() || self.steps >= self.limits.maximum_steps || misplaced
        {
            return Err(NsfDgError::Input("resolved joint exterior requires a finite correctly located source within its horizon"));
        }
        let owner = self.reservoir_radiation.as_ref().ok_or(NsfDgError::Input(
            "resolved joint exterior is not attached"))?;
        let aperture = owner.aperture;
        if breath.is_some_and(|(source, _)| source == aperture) {
            return Err(NsfDgError::Input("prescribed breath and passive exterior require distinct original apertures"));
        }
        let pressure_bracket = owner.input.pressure_bracket_pa;
        let temperature = owner.temperature_k;
        let relation = self.load.independent_mouth_response().map_err(NsfDgError::Acoustic)?;
        let upstream_relation = self.upstream.as_ref()
            .map(|owner| owner.load.port_response().map_err(NsfDgError::Acoustic)).transpose()?;
        let upstream_reference = self.upstream.as_ref().map(|owner| owner.reference);
        let upstream_inlet = self.upstream.as_ref().map(|owner| owner.load.inlet());
        let limits = self.limits;
        let side = usize::from(self.high);
        let other = 1-side;
        let patch = self.patch_boundaries();
        let reference = self.reference;
        let mut evaluations = 0usize;
        let mut last_residual = 0.0;
        let zone = &mut self.zone;
        let mut evaluate = |pressure: f64, incoming: [f64; 2]| -> Result<([f64; 3], [f64; 3]), NsfDgError> {
            poll(gate)?;
            // Reserve the final original interval that will actually publish.
            if evaluations + 1 >= limits.maximum_evaluations {
                return Err(NsfDgError::Convergence { residual_pa: last_residual, evaluations });
            }
            evaluations += 1;
            let mut drives = incoming.map(NsfBoundaryDrive::IncomingCharacteristic);
            if upstream_relation.is_none() { drives[other] = external; }
            let trial = preview(zone, drives, side, patch, aperture, breath, pressure, temperature, limits, gate)?;
            let bore = coupled_exchange(&trial, side)?;
            let mouth = exchange(&trial, aperture)?;
            let flow = [bore.normal_volume_flow_m3_s, mouth.normal_volume_flow_m3_s];
            let pressure = [bore.mean_acoustic_pressure_pa, mouth.mean_acoustic_pressure_pa];
            let mut residuals = [0.0; 3];
            let mut scales = [0.0; 3];
            for i in 0..2 {
                let terms = [relation.impedance_pa_s_m3[i][0]*flow[0],
                    relation.impedance_pa_s_m3[i][1]*flow[1]];
                residuals[i] = pressure[i] - (relation.history_pressure_pa[i] + terms[0] + terms[1]);
                scales[i] = pressure[i].abs() + relation.history_pressure_pa[i].abs()
                    + terms[0].abs() + terms[1].abs();
            }
            if let Some(upstream) = upstream_relation {
                let face = coupled_exchange(&trial, other)?;
                let dynamic = upstream.impedance_pa_s_m3*face.normal_volume_flow_m3_s;
                residuals[2] = face.mean_acoustic_pressure_pa - (upstream.history_pressure_pa + dynamic + rise);
                scales[2] = face.mean_acoustic_pressure_pa.abs() + upstream.history_pressure_pa.abs()
                    + dynamic.abs() + rise.abs();
            }
            if !residuals.into_iter().chain(scales).all(f64::is_finite) {
                return Err(NsfDgError::Input("resolved joint exterior pressure equation overflowed"));
            }
            last_residual = residuals.into_iter().map(f64::abs).fold(0.0, f64::max);
            Ok((residuals, scales))
        };
        let (pressure_offset, (downstream, upstream_incoming)) =
            upstream::solve_pressure(limits, pressure_bracket, |pressure| {
                if upstream_relation.is_some() {
                    let (incoming_upstream, (incoming_main, value)) =
                        upstream::solve_scalar(limits, None, |incoming_upstream| {
                            let mut incoming = [incoming_upstream; 2];
                            let (incoming_main, value) =
                                upstream::solve_scalar(limits, None, |incoming_main| {
                                    incoming[side] = incoming_main;
                                    let value = evaluate(pressure, incoming)?;
                                    Ok((value.0[0], value.1[0], value))
                                })?;
                            Ok((value.0[2], value.1[2], (incoming_main, value)))
                        })?;
                    Ok((value.0[1], value.1[1], (incoming_main, Some(incoming_upstream))))
                } else {
                    let (incoming_main, value) =
                        upstream::solve_scalar(limits, None, |incoming_main| {
                            let mut incoming = [0.0; 2];
                            incoming[side] = incoming_main;
                            let value = evaluate(pressure, incoming)?;
                            Ok((value.0[0], value.1[0], value))
                        })?;
                    Ok((value.0[1], value.1[1], (incoming_main, None)))
                }
            })?;
        drop(evaluate);
        evaluations += 1;
        let mut drives = [external; 2];
        drives[side] = NsfBoundaryDrive::IncomingCharacteristic(downstream);
        if let Some(value) = upstream_incoming { drives[other] = NsfBoundaryDrive::IncomingCharacteristic(value); }
        let fluid = preview(zone, drives, side, patch, aperture, breath, pressure_offset, temperature, limits, gate)?;
        let faces = *fluid.boundary_exchange();
        let coupled_patch = fluid.patch_exchange_on(side != 0).copied();
        let source_patch = fluid.patch_exchange_on(other != 0).copied();
        let bore = coupled_exchange(&fluid, side)?;
        let source_face = coupled_exchange(&fluid, other)?;
        let mouth = exchange(&fluid, aperture)?;
        let resolved_exchange = *fluid.reservoir_exchange_at(aperture).ok_or(NsfDgError::Input(
            "joint exterior lost its conservative reservoir exchange"))?;
        let dt = self.load.propagator().dt();
        let inlet = self.load.inlet();
        let acoustic = self.load.preview_step_with_mouth_flow(
            bore.normal_volume_flow_m3_s, mouth.normal_volume_flow_m3_s).map_err(NsfDgError::Acoustic)?;
        let bore_pressure = acoustic.output().port_outputs.iter().find(|port| port.id == inlet)
            .ok_or(NsfDgError::Input("joint DG inlet observation disappeared"))?.conjugate_rate;
        let radiation = *acoustic.independent_mouth_frame().ok_or(NsfDgError::Input(
            "joint exterior lost its independent mouth coordinate"))?;
        let bore_relation = DgPortResponse {
            history_pressure_pa: relation.history_pressure_pa[0]
                + relation.impedance_pa_s_m3[0][1]*mouth.normal_volume_flow_m3_s,
            impedance_pa_s_m3: relation.impedance_pa_s_m3[0][0],
        };
        let mouth_relation = DgPortResponse {
            history_pressure_pa: relation.history_pressure_pa[1]
                + relation.impedance_pa_s_m3[1][0]*bore.normal_volume_flow_m3_s,
            impedance_pa_s_m3: relation.impedance_pa_s_m3[1][1],
        };
        let bore_audit = upstream::audit(bore, bore_pressure, acoustic.inlet_work_j(),
            0.0, bore_relation, reference, dt, limits)?;
        let mouth_audit = upstream::audit(mouth, radiation.mean_pressure_pa, radiation.load_work_j,
            0.0, mouth_relation, reference, dt, limits)?;
        let source = self.upstream.as_mut().map(|owner|
            owner.load.preview_step(source_face.normal_volume_flow_m3_s)
                .map_err(NsfDgError::Acoustic)).transpose()?;
        let upstream_values = if let Some(source) = &source {
            let pressure = source.output().port_outputs.iter()
                .find(|port| Some(port.id) == upstream_inlet)
                .ok_or(NsfDgError::Input("joint upstream pressure observation disappeared"))?.conjugate_rate;
            let audit = upstream::audit(source_face, pressure, source.inlet_work_j(), rise,
                upstream_relation.expect("original upstream response"),
                upstream_reference.expect("original upstream reference"), dt, limits)?;
            Some((pressure, audit, source.terminal_frames().to_vec()))
        } else { None };
        let flow = *fluid.energy_ledger();
        let reservoir_mass_out_kg = fluid.reservoir_mass_out_kg();
        let reservoir_flows = fluid.reservoir_flows().to_vec();
        let reservoir_exchanges = fluid.reservoir_exchanges().to_vec();
        let mass_change = fluid.grid().total_mass_kg() - fluid.previous_mass_kg();
        let energy_change = fluid.grid().total_energy_j() - fluid.previous_energy_j();
        let terminals = acoustic.terminal_frames().to_vec();
        if ![mass_change, energy_change, reservoir_mass_out_kg].iter().all(|x| x.is_finite()) {
            return Err(NsfDgError::Input("resolved joint exterior extensive accounts overflowed"));
        }
        let finite_candidate = finite_reservoir::prepare(
            self.finite_reservoir.as_ref(), &fluid, drive, dt, gate,
        )?;
        poll(gate)?;
        fluid.accept();
        let duct = acoustic.accept();
        let upstream = source.zip(upstream_values).map(|(source, (pressure, audit, terminals))|
            NsfDgUpstreamFrame {
                face: other, incoming_characteristic_m_s: upstream_incoming.expect("solved upstream"),
                duct: source.accept(), terminals,
                inlet_flow_m3_s: source_face.normal_volume_flow_m3_s, inlet_pressure_pa: pressure,
                source_work_j: audit.source_work, pressure_residual_pa: audit.pressure,
                projection_covariance_j: audit.covariance, acoustic_work_residual_j: audit.work,
                port_work_residual_j: audit.port, mass_projection_difference_kg: audit.mass,
            });
        let finite_reservoir = finite_reservoir::publish(&mut self.finite_reservoir, finite_candidate);
        self.steps += 1;
        Ok(NsfDgFrame {
            step: self.steps, evaluations, incoming_characteristic_m_s: downstream, drive,
            upstream, faces, coupled_patch, source_patch, reservoir_mass_out_kg, reservoir_flows,
            reservoir_exchanges, finite_reservoir,
            reservoir_radiation: Some(NsfDgReservoirFrame {
                resolved_exchange,
                ghost_pressure_offset_pa: pressure_offset,
                mean_acoustic_pressure_pa: mouth.mean_acoustic_pressure_pa,
                outward_flow_m3_s: mouth.normal_volume_flow_m3_s,
                radiation_pressure_pa: radiation.mean_pressure_pa,
                radiation_work_j: radiation.load_work_j,
                pressure_residual_pa: mouth_audit.pressure,
                projection_covariance_j: mouth_audit.covariance,
                acoustic_work_residual_j: mouth_audit.work,
                mass_projection_difference_kg: mouth_audit.mass,
            }),
            flow, flow_mass_change_kg: mass_change, flow_energy_change_j: energy_change,
            duct, terminals, inlet_flow_m3_s: bore.normal_volume_flow_m3_s,
            inlet_pressure_pa: bore_pressure, pressure_residual_pa: bore_audit.pressure,
            projection_covariance_j: bore_audit.covariance, acoustic_work_residual_j: bore_audit.work,
            port_work_residual_j: bore_audit.port, mass_projection_difference_kg: bore_audit.mass,
        })
    }
}

#[allow(clippy::too_many_arguments)]
fn preview<'a, L: Lattice3<27>>(
    zone: &'a mut NsfGrid3<L>, drives: [NsfBoundaryDrive; 2], side: usize,
    patch: [Option<NsfBoundaryDrive>; 2], aperture: usize,
    breath: Option<(usize, NsfOpenState)>, pressure_offset: f64, temperature: f64,
    limits: NsfDgLimits, gate: &CancelGate,
) -> Result<NsfCharacteristicTrial<'a, L>, NsfDgError> {
    let reference = zone.characteristic_reference_gas(side == 1).map_err(NsfDgError::Flow)?;
    let state = zone.reservoir_aperture_state_for(aperture, reference.pressure + pressure_offset, temperature)
        .map_err(NsfDgError::Flow)?;
    // Stage both original ghost controls together; no candidate sees a lagged
    // escape pressure or publishes the breath state before joint acceptance.
    let controls = [breath.unwrap_or((aperture, state)), (aperture, state)];
    let controls = if breath.is_some() { &controls[..] } else { &controls[1..] };
    preview_interval(zone, drives, patch, controls, limits, gate)
}

fn exchange<L: Lattice3<27>>(
    trial: &NsfCharacteristicTrial<'_, L>, aperture: usize,
) -> Result<NsfBoundaryExchange, NsfDgError> {
    let source = trial.reservoir_exchange_at(aperture).ok_or(NsfDgError::Input(
        "resolved exterior has no original conservative reservoir work quadrature"))?;
    // Match pressure and audit work from the SAME cached FV predictors that
    // update total energy. Endpoint midpoint-pressure times midpoint-velocity
    // is a separate observation quadrature and cannot replace this transfer.
    // Streamed mass remains independent; no enthalpy, conduction or viscous
    // work is folded into the acoustic projection.
    Ok(NsfBoundaryExchange {
        mass_out_kg: source.mass_out_kg,
        mean_pressure_pa: source.mean_pressure_pa,
        normal_volume_flow_m3_s: source.normal_volume_flow_m3_s,
        pressure_work_out_j: source.pressure_work_out_j,
        mean_acoustic_pressure_pa: source.mean_acoustic_pressure_pa,
        acoustic_pressure_work_out_j: source.acoustic_pressure_work_out_j,
        ..NsfBoundaryExchange::default()
    })
}
