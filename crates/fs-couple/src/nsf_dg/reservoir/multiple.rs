//! Simultaneous original NSF, DG and shared exterior pressure equations.
//! All rectangular escape channels belong to one original mutual radiation
//! owner. Newton--Krylov changes trial controls, never physical histories.
use super::*;
use fs_lbm::nsf3::MAX_NSF_RESERVOIR_APERTURES;
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{NewtonKrylovConfig, NewtonKrylovState, SolverRunProgress};
use std::cell::{Cell, RefCell};

/// Explicit numerical policy for the simultaneous pressure equations.
/// Coordinates are normalized by each authored bracket's midpoint/half-width;
/// residuals remain Pa. Every source evaluation, including numerical Jacobian
/// actions, consumes the existing NSF interval/cell-work budget. A trial or
/// difference stencil outside ANY authored bracket refuses without clipping.
/// Original source errors remain refusals, including during globalization.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NsfDgRadiationSolve {
    /// Maximum normalized-coordinate perturbation of the original numerical
    /// Jacobian action; finite and strictly between zero and one.
    pub difference_step: f64,
    /// Maximum original Newton outer attempts, in 1..=256.
    pub maximum_iterations: usize,
    /// Original explicit iteration, tolerance and globalization controls.
    /// Its convergence flag does not replace the final per-port physical audit.
    pub newton: NewtonKrylovConfig,
}

impl NsfDgRadiationSolve {
    /// Canonical admission shared by native input and programmatic attachment.
    pub fn validate(&self) -> Result<(), NsfDgError> {
        self.newton.validate().map_err(NsfDgError::Nonlinear)?;
        if !self.difference_step.is_finite()
            || !(0.0 < self.difference_step && self.difference_step < 1.0)
            || !(1..=256).contains(&self.maximum_iterations)
            || !(1..=MAX_NSF_RESERVOIR_APERTURES + 2).contains(&self.newton.linear_restart)
            || !(1..=256).contains(&self.newton.max_linear_cycles)
        {
            return Err(NsfDgError::Input(
                "joint reservoir solve needs bounded original Newton work and a normalized difference step in (0,1)",
            ));
        }
        Ok(())
    }
}

pub(in crate::nsf_dg) struct ReservoirRadiations {
    inputs: Vec<(usize, NsfDgReservoirRadiation)>,
    temperatures_k: Vec<f64>,
    solve: NsfDgRadiationSolve,
    // Warm-start controls publish only with all physical participants. A cold
    // later upstream attachment can change the dimension before the first step.
    previous_coordinates: Vec<f64>,
}

fn numerical(error: FiniteDifferenceError<NsfDgError>) -> NsfDgError {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(message) => NsfDgError::Input(message),
    }
}

fn coordinates(normalized: &[f64], domains: &[[f64; 2]]) -> Result<Vec<f64>, NsfDgError> {
    if normalized.len() != domains.len() {
        return Err(NsfDgError::Input("joint reservoir pressure coordinate count changed"));
    }
    normalized.iter().zip(domains).map(|(&value, &[low, high])| {
        let midpoint = f64::midpoint(low, high);
        let half = 0.5 * high - 0.5 * low;
        let physical = midpoint + half * value;
        if !value.is_finite() || !(-1.0..=1.0).contains(&value)
            || !midpoint.is_finite() || !half.is_finite() || half <= 0.0
            || !physical.is_finite() || physical < low || physical > high
        {
            return Err(NsfDgError::Input(
                "joint reservoir Newton trial or difference stencil leaves its authored pressure/characteristic bracket",
            ));
        }
        Ok(physical)
    }).collect()
}

impl<L: Lattice3<27>> NsfDg<L> {
    /// Attach two or more original escape apertures to ONE computed rectangular
    /// and circular self/mutual exterior. Their distinct flows remain independent
    /// of the bore inlet, finite source and optional upstream duct. All original
    /// geometry, gas, pressure-domain and radiation work admissions still apply.
    ///
    /// The authored order fixes the shared exterior channel order. This cold
    /// attachment is exclusive with the original singleton attachment. Each
    /// rectangle retains one uniform pressure mode, while per-tile microphones
    /// remain observations of the original resolved flow.
    pub fn with_reservoir_radiations(
        mut self, inputs: Vec<(usize, NsfDgReservoirRadiation)>,
        solve: NsfDgRadiationSolve, gate: &CancelGate,
    ) -> Result<Self, NsfDgError> {
        poll(gate)?;
        solve.validate()?;
        if !(2..=MAX_NSF_RESERVOIR_APERTURES).contains(&inputs.len())
            || self.steps != 0 || self.reservoir_radiation.is_some()
            || self.reservoir_radiations.is_some()
            || inputs.iter().enumerate().any(|(i, (aperture, _))|
                inputs[..i].iter().any(|(prior, _)| prior == aperture))
        {
            return Err(NsfDgError::Input(
                "joint reservoirs require distinct cold apertures and one exclusive exterior attachment",
            ));
        }
        let mut surfaces = Vec::with_capacity(inputs.len());
        let mut temperatures_k = Vec::with_capacity(inputs.len());
        let mut rate = None;
        for &(aperture, input) in &inputs {
            let (surface, temperature, source_rate) =
                self.reservoir_radiation_surface(aperture, input, gate)?;
            if rate.is_some_and(|previous| previous != source_rate) {
                return Err(NsfDgError::Input("joint reservoirs changed the original exterior bake clock"));
            }
            rate = Some(source_rate);
            surfaces.push(surface);
            temperatures_k.push(temperature);
        }
        self.load = self.load.with_independent_mouths(
            &surfaces, rate.expect("nonempty admitted apertures"), gate,
        ).map_err(|error| if gate.is_requested() { NsfDgError::Cancelled }
            else { NsfDgError::Acoustic(AcousticRealizeError::Nonlinear(error.to_string())) })?;
        poll(gate)?;
        self.reservoir_radiations = Some(ReservoirRadiations {
            inputs, temperatures_k, solve, previous_coordinates: Vec::new(),
        });
        Ok(self)
    }

    /// Ordered declarations of the plural owner; empty on the legacy singleton.
    pub fn reservoir_radiations(&self) -> &[(usize, NsfDgReservoirRadiation)] {
        self.reservoir_radiations.as_ref().map_or(&[], |owner| owner.inputs.as_slice())
    }

    /// Original passive ownership at one aperture, for either attachment mode.
    pub fn reservoir_radiation_at(&self, aperture: usize) -> Option<&NsfDgReservoirRadiation> {
        self.reservoir_radiation.as_ref().filter(|owner| owner.aperture == aperture)
            .map(|owner| &owner.input).or_else(|| self.reservoir_radiations()
                .iter().find(|(index, _)| *index == aperture).map(|(_, input)| input))
    }

    /// Explicit original nonlinear policy; absent on the legacy singleton path.
    pub fn radiation_solve(&self) -> Option<&NsfDgRadiationSolve> {
        self.reservoir_radiations.as_ref().map(|owner| &owner.solve)
    }

    pub(in crate::nsf_dg) fn step_with_multiple_reservoir_drive(
        &mut self, drive: NsfDgDrive, gate: &CancelGate,
    ) -> Result<NsfDgFrame, NsfDgError> {
        poll(gate)?;
        self.admit_finite_reservoir_drive(drive)?;
        let (external, rise, breath) = match drive {
            NsfDgDrive::Boundary(value) if self.upstream.is_none() => (value, 0.0, None),
            NsfDgDrive::UpstreamPressureRise(value) if self.upstream.is_some() =>
                (NsfBoundaryDrive::IncomingCharacteristic(0.0), value, None),
            NsfDgDrive::ReservoirAperture { aperture, pressure_pa, temperature_k } => {
                let state = self.zone.reservoir_aperture_state_for(aperture, pressure_pa, temperature_k)
                    .map_err(NsfDgError::Flow)?;
                (NsfBoundaryDrive::IncomingCharacteristic(0.0), 0.0, Some((aperture, state)))
            }
            NsfDgDrive::FiniteReservoir { aperture }
            | NsfDgDrive::FiniteReservoirActuation { aperture, .. } =>
                (NsfBoundaryDrive::IncomingCharacteristic(0.0), 0.0,
                    Some((aperture, self.finite_reservoir_state(aperture)?))),
            _ => return Err(NsfDgError::Input(
                "joint reservoirs require a correctly located distinct physical source",
            )),
        };
        if !drive.is_finite() || self.steps >= self.limits.maximum_steps {
            return Err(NsfDgError::Input("joint reservoir drive is nonfinite or its horizon is exhausted"));
        }
        let owner = self.reservoir_radiations.as_ref()
            .ok_or(NsfDgError::Input("joint reservoirs are not attached"))?;
        if breath.is_some_and(|(aperture, _)|
            owner.inputs.iter().any(|(escape, _)| *escape == aperture))
        {
            return Err(NsfDgError::Input("breath and every passive escape must have distinct original owners"));
        }
        let inputs = owner.inputs.clone();
        let temperatures = owner.temperatures_k.clone();
        let solve = owner.solve;
        let relation = self.load.independent_mouths_response().map_err(NsfDgError::Acoustic)?;
        let width = inputs.len() + 1;
        if relation.history_pressure_pa.len() != width
            || relation.impedance_pa_s_m3.len() != width * width
            || !relation.history_pressure_pa.iter().chain(&relation.impedance_pa_s_m3)
                .all(|value| value.is_finite())
        {
            return Err(NsfDgError::Input("joint exterior lost its complete original pressure basis"));
        }
        let upstream_relation = self.upstream.as_ref()
            .map(|owner| owner.load.port_response().map_err(NsfDgError::Acoustic)).transpose()?;
        let upstream_reference = self.upstream.as_ref().map(|owner| owner.reference);
        let upstream_inlet = self.upstream.as_ref().map(|owner| owner.load.inlet());
        let limits = self.limits;
        let side = usize::from(self.high);
        let other = 1 - side;
        let patch = self.patch_boundaries();
        let reference = self.reference;
        let mut domains = vec![limits.incoming_bracket_m_s];
        domains.extend(inputs.iter().map(|(_, input)| input.pressure_bracket_pa));
        if upstream_relation.is_some() { domains.push(limits.incoming_bracket_m_s); }
        let seed = if owner.previous_coordinates.len() == domains.len() {
            owner.previous_coordinates.clone()
        } else {
            domains.iter().map(|&[low, high]| {
                if low <= 0.0 && high >= 0.0 {
                    -f64::midpoint(low, high) / (0.5 * high - 0.5 * low)
                } else { 0.0 }
            }).collect()
        };
        let evaluations = Cell::new(0usize);
        let last_residual = Cell::new(0.0);
        let normalized = {
            let zone = RefCell::new(&mut self.zone);
            let problem = FallibleFiniteDifference::new(domains.len(), solve.difference_step,
                |normalized: &[f64]| {
                    poll(gate)?;
                    let physical = coordinates(normalized, &domains)?;
                    if evaluations.get() + 1 >= limits.maximum_evaluations {
                        return Err(NsfDgError::Convergence {
                            residual_pa: last_residual.get(), evaluations: evaluations.get(),
                        });
                    }
                    evaluations.set(evaluations.get() + 1);
                    let mut zone = zone.borrow_mut();
                    let trial = preview_multiple(&mut **zone, &physical, &inputs, &temperatures,
                        external, breath, side, patch, upstream_relation.is_some(),
                        reference.pressure, limits, gate)?;
                    let bore = coupled_exchange(&trial, side)?;
                    let mut flows = vec![bore.normal_volume_flow_m3_s];
                    let mut pressures = vec![bore.mean_acoustic_pressure_pa];
                    for &(aperture, _) in &inputs {
                        let mouth = exchange(&trial, aperture)?;
                        flows.push(mouth.normal_volume_flow_m3_s);
                        pressures.push(mouth.mean_acoustic_pressure_pa);
                    }
                    let mut residuals = Vec::with_capacity(domains.len());
                    for row in 0..width {
                        let mut expected = relation.history_pressure_pa[row];
                        for (column, flow) in flows.iter().enumerate() {
                            expected += relation.impedance_pa_s_m3[row * width + column] * flow;
                        }
                        residuals.push(pressures[row] - expected);
                    }
                    if let Some(upstream) = upstream_relation {
                        let face = coupled_exchange(&trial, other)?;
                        residuals.push(face.mean_acoustic_pressure_pa - (upstream.history_pressure_pa
                            + upstream.impedance_pa_s_m3 * face.normal_volume_flow_m3_s + rise));
                    }
                    if residuals.iter().any(|value| !value.is_finite()) {
                        return Err(NsfDgError::Input("joint reservoir pressure equation overflowed"));
                    }
                    last_residual.set(residuals.iter().copied().map(f64::abs).fold(0.0, f64::max));
                    Ok(residuals)
                }).map_err(numerical)?;
            // The original constructor evaluates once; the source callback
            // itself polls cancellation and charges that physical preview.
            let state = NewtonKrylovState::new(&problem, seed, solve.newton);
            if let Some(error) = problem.take_failure() { return Err(numerical(error)); }
            let mut state = state.map_err(NsfDgError::Nonlinear)?;
            let report = state.run_controlled(&problem, solve.maximum_iterations,
                &mut || gate.is_requested() || problem.has_failed());
            if let Some(error) = problem.take_failure() { return Err(numerical(error)); }
            poll(gate)?;
            if report.progress != SolverRunProgress::Complete {
                return Err(NsfDgError::Input("joint reservoir nonlinear callbacks did not complete"));
            }
            if !report.report.converged {
                return Err(NsfDgError::NonlinearConvergence {
                    evaluations: evaluations.get(), report: report.report,
                });
            }
            state.x
        };
        let physical = coordinates(&normalized, &domains)?;
        evaluations.set(evaluations.get() + 1);
        let fluid = preview_multiple(&mut self.zone, &physical, &inputs, &temperatures,
            external, breath, side, patch, upstream_relation.is_some(),
            reference.pressure, limits, gate)?;
        let faces = *fluid.boundary_exchange();
        let coupled_patch = fluid.patch_exchange_on(side != 0).copied();
        let source_patch = fluid.patch_exchange_on(other != 0).copied();
        let bore = coupled_exchange(&fluid, side)?;
        let source_face = coupled_exchange(&fluid, other)?;
        let mouths = inputs.iter().map(|(aperture, _)| exchange(&fluid, *aperture))
            .collect::<Result<Vec<_>, _>>()?;
        let mouth_flows: Vec<_> = mouths.iter().map(|mouth| mouth.normal_volume_flow_m3_s).collect();
        let mut flows = vec![bore.normal_volume_flow_m3_s];
        flows.extend_from_slice(&mouth_flows);
        let scalar_relation = |row: usize| {
            let mut history = relation.history_pressure_pa[row];
            for (column, flow) in flows.iter().enumerate() {
                if column != row { history += relation.impedance_pa_s_m3[row * width + column] * flow; }
            }
            DgPortResponse { history_pressure_pa: history,
                impedance_pa_s_m3: relation.impedance_pa_s_m3[row * width + row] }
        };
        let dt = self.load.propagator().dt();
        let inlet = self.load.inlet();
        let acoustic = self.load.preview_step_with_mouth_flows(
            bore.normal_volume_flow_m3_s, &mouth_flows).map_err(NsfDgError::Acoustic)?;
        let bore_pressure = acoustic.output().port_outputs.iter().find(|port| port.id == inlet)
            .ok_or(NsfDgError::Input("joint DG inlet observation disappeared"))?.conjugate_rate;
        let radiation = acoustic.independent_mouths_frame()
            .ok_or(NsfDgError::Input("joint exterior lost its independent mouth receipts"))?;
        if radiation.mouths.len() != mouths.len() {
            return Err(NsfDgError::Input("joint exterior receipt count changed"));
        }
        let bore_audit = upstream::audit(bore, bore_pressure, acoustic.inlet_work_j(),
            0.0, scalar_relation(0), reference, dt, limits)?;
        let mut reservoir_radiations = Vec::with_capacity(inputs.len());
        for (index, (&mouth, receipt)) in mouths.iter().zip(&radiation.mouths).enumerate() {
            let aperture = inputs[index].0;
            let audit = upstream::audit(mouth, receipt.mean_pressure_pa, receipt.load_work_j,
                0.0, scalar_relation(index + 1), reference, dt, limits)?;
            let resolved_exchange = *fluid.reservoir_exchange_at(aperture)
                .ok_or(NsfDgError::Input("joint exterior lost an original conservative reservoir exchange"))?;
            reservoir_radiations.push(NsfDgReservoirFrame {
                resolved_exchange, ghost_pressure_offset_pa: physical[index + 1],
                mean_acoustic_pressure_pa: mouth.mean_acoustic_pressure_pa,
                outward_flow_m3_s: mouth.normal_volume_flow_m3_s,
                radiation_pressure_pa: receipt.mean_pressure_pa,
                radiation_work_j: receipt.load_work_j,
                pressure_residual_pa: audit.pressure, projection_covariance_j: audit.covariance,
                acoustic_work_residual_j: audit.work, mass_projection_difference_kg: audit.mass,
            });
        }
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
        if ![mass_change, energy_change, reservoir_mass_out_kg].iter().all(|value| value.is_finite()) {
            return Err(NsfDgError::Input("joint reservoir extensive accounts overflowed"));
        }
        let finite_candidate = finite_reservoir::prepare(
            self.finite_reservoir.as_ref(), &fluid, drive, dt, gate,
        )?;
        poll(gate)?;
        // Nothing below can refuse: one fluid candidate, both original ducts,
        // one shared exterior memory, and optional finite source publish once.
        fluid.accept();
        let duct = acoustic.accept();
        let upstream = source.zip(upstream_values).map(|(source, (pressure, audit, terminals))|
            NsfDgUpstreamFrame {
                face: other, incoming_characteristic_m_s: physical[width],
                duct: source.accept(), terminals,
                inlet_flow_m3_s: source_face.normal_volume_flow_m3_s, inlet_pressure_pa: pressure,
                source_work_j: audit.source_work, pressure_residual_pa: audit.pressure,
                projection_covariance_j: audit.covariance, acoustic_work_residual_j: audit.work,
                port_work_residual_j: audit.port, mass_projection_difference_kg: audit.mass,
            });
        let finite_reservoir = finite_reservoir::publish(&mut self.finite_reservoir, finite_candidate);
        self.reservoir_radiations.as_mut().expect("original joint owner")
            .previous_coordinates = normalized;
        self.steps += 1;
        Ok(NsfDgFrame {
            step: self.steps, evaluations: evaluations.get(), incoming_characteristic_m_s: physical[0],
            drive, upstream, faces, coupled_patch, source_patch, bore_exchange: bore,
            reservoir_mass_out_kg,
            reservoir_flows, reservoir_exchanges, reservoir_radiation: None,
            reservoir_radiations, finite_reservoir, flow, flow_mass_change_kg: mass_change,
            flow_energy_change_j: energy_change, duct, terminals,
            inlet_flow_m3_s: bore.normal_volume_flow_m3_s, inlet_pressure_pa: bore_pressure,
            pressure_residual_pa: bore_audit.pressure, projection_covariance_j: bore_audit.covariance,
            acoustic_work_residual_j: bore_audit.work, port_work_residual_j: bore_audit.port,
            mass_projection_difference_kg: bore_audit.mass,
        })
    }
}

#[allow(clippy::too_many_arguments)]
fn preview_multiple<'a, L: Lattice3<27>>(
    zone: &'a mut NsfGrid3<L>, physical: &[f64],
    inputs: &[(usize, NsfDgReservoirRadiation)], temperatures: &[f64],
    external: NsfBoundaryDrive, breath: Option<(usize, NsfOpenState)>, side: usize,
    patch: [Option<NsfBoundaryDrive>; 2], upstream: bool, reference_pressure: f64,
    limits: NsfDgLimits, gate: &CancelGate,
) -> Result<NsfCharacteristicTrial<'a, L>, NsfDgError> {
    let mut drives = [external; 2];
    drives[side] = NsfBoundaryDrive::IncomingCharacteristic(physical[0]);
    if upstream {
        drives[1 - side] = NsfBoundaryDrive::IncomingCharacteristic(physical[inputs.len() + 1]);
    }
    let mut controls = Vec::with_capacity(inputs.len() + usize::from(breath.is_some()));
    if let Some(breath) = breath { controls.push(breath); }
    for (index, &(aperture, _)) in inputs.iter().enumerate() {
        let state = zone.reservoir_aperture_state_for(
            aperture, reference_pressure + physical[index + 1], temperatures[index],
        ).map_err(NsfDgError::Flow)?;
        controls.push((aperture, state));
    }
    preview_interval(zone, drives, patch, &controls, limits, gate)
}
