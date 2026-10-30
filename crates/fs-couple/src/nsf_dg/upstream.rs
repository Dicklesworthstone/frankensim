//! Two original acoustic loads around one original resolved gas volume.
//!
//! Both face pressures are solved on the same NSF interval. The upstream ideal
//! series source supplies p_zone - p_load = pressure_rise; its independent work
//! is -pressure_rise * Q_out * dt. Every duct's inlet points AWAY from the zone.
//! Duct mass/entropy advection is not represented by this acoustic interface.

use super::*;

pub(super) struct UpstreamLoad {
    pub(super) load: DgPortLoad,
    pub(super) reference: GasState,
}

/// Finite upstream acoustic state and its independently declared source work.
#[derive(Debug, Clone)]
pub struct NsfDgUpstreamFrame {
    /// Original low/high face index; the main bore occupies the other face.
    pub face: usize,
    /// Accepted original incoming characteristic [m/s].
    pub incoming_characteristic_m_s: f64,
    /// Original aggregate DG/branch accounting.
    pub duct: DgStepOutput,
    /// Original upstream terminal ordering and pressure/flow accounts.
    pub terminals: Vec<DgTerminalFrame>,
    /// Held flow OUT of the zone and into the upstream duct [m³/s].
    pub inlet_flow_m3_s: f64,
    /// Actual upstream DG interval-mean inlet pressure [Pa].
    pub inlet_pressure_pa: f64,
    /// Positive work supplied to the complete coupled system [J].
    pub source_work_j: f64,
    /// Resolved pressure minus upstream pressure and declared source rise [Pa].
    pub pressure_residual_pa: f64,
    /// Original resolved covariance, never reclassified as heat [J].
    pub projection_covariance_j: f64,
    /// Resolved outward work + source work - upstream DG inlet work [J].
    pub acoustic_work_residual_j: f64,
    /// Original integrated DG inlet work minus its pressure/flow product [J].
    pub port_work_residual_j: f64,
    /// Original lattice mass flux minus reference-density swept volume [kg].
    pub mass_projection_difference_kg: f64,
}

impl<L: Lattice3<27>> NsfDg<L> {
    /// Attach an original finite acoustic duct to the other resolved face or
    /// its explicitly selected source aperture. Surrounding opening flux is
    /// never assigned to this duct or to the series pressure-rise source.
    /// Its propagation, material wall, molecular and terminal histories remain
    /// intact. Construction requires cold matching gas, area and physical clock.
    /// No steady mass/enthalpy flow into that acoustic duct is inferred.
    pub fn with_upstream_load(
        mut self,
        load: DgPortLoad,
        gate: &CancelGate,
    ) -> Result<Self, NsfDgError> {
        poll(gate)?;
        let high = !self.high;
        let reference = self
            .zone
            .characteristic_reference_gas(high)
            .map_err(NsfDgError::Flow)?;
        let traces = if self.source_patch().is_some() { self.zone.open_patch_traces_on(high) }
            else { self.zone.open_boundary_traces(high) }
            .ok_or(NsfDgError::Input(
                "upstream coupling requires original selected source quadrature",
            ))?;
        let area = traces
            .first()
            .ok_or(NsfDgError::Input("upstream face is empty"))?
            .area_m2
            * traces.len() as f64;
        let other_area = load.duct().inlet_area_m2();
        if self.steps != 0
            || self.upstream.is_some()
            || load.state().step_count() != 0
            || load.has_driven_mouth()
            || load.duct().inlet_gas() != &reference
            || !area.is_finite()
            || area <= 0.0
            || !other_area.is_finite()
            || (area - other_area).abs() > 32.0 * f64::EPSILON * area.max(other_area)
            || !same_clock(self.load.propagator().dt(), load.propagator().dt())
        {
            return Err(NsfDgError::Input(
                "upstream DG attachment requires one cold matching original gas, area and clock",
            ));
        }
        poll(gate)?;
        self.upstream = Some(UpstreamLoad { load, reference });
        Ok(self)
    }

    /// The finite original upstream acoustic and material owner, if present.
    pub fn upstream_load(&self) -> Option<&DgPortLoad> {
        self.upstream.as_ref().map(|owner| &owner.load)
    }

    /// Jointly solve BOTH original characteristic faces and acoustic histories.
    /// The declared pressure rise is placed between upstream duct and zone.
    /// It is neither an absolute boundary pressure nor a finite player model.
    /// Nested safeguarded scalar roots reuse each original trial from the same
    /// accepted state. For the explicit source aperture, accepted incoming
    /// controls seed the roots and successive outer trials reuse the solved
    /// main-face incoming, with a full bracket fallback at the original
    /// residual bound. Existing complete-face and downstream-aperture owners
    /// retain their original cold-start trial sequence.
    /// The existing total evaluation/cell-step budget covers
    /// every inner, outer and final-restaging evaluation together.
    pub fn step_with_upstream_pressure_rise(
        &mut self,
        pressure_rise_pa: f64,
        gate: &CancelGate,
    ) -> Result<NsfDgFrame, NsfDgError> {
        self.step_with_upstream_conditions(pressure_rise_pa, &[],
            NsfDgDrive::UpstreamPressureRise(pressure_rise_pa), gate)
    }

    pub(super) fn step_with_upstream_conditions(
        &mut self, pressure_rise_pa: f64, reservoir_states: &[(usize, NsfOpenState)],
        control: NsfDgDrive, gate: &CancelGate,
    ) -> Result<NsfDgFrame, NsfDgError> {
        poll(gate)?;
        self.admit_finite_reservoir_drive(control)?;
        if self.reservoir_radiation.is_some() {
            return self.step_with_reservoir_drive(control, gate);
        }
        if !pressure_rise_pa.is_finite() || self.steps >= self.limits.maximum_steps {
            return Err(NsfDgError::Input(
                "upstream pressure rise is nonfinite or the horizon is exhausted",
            ));
        }
        let patch = self.patch_boundaries();
        let side = usize::from(self.high);
        let other = 1 - side;
        let seeds = if self.source_patch().is_some() {
            let mut accepted = self.zone.boundary_drives().ok_or(NsfDgError::Input(
                "two-sided coupling requires its original accepted boundary controls",
            ))?;
            for face in 0..2 {
                if patch[face].is_some() {
                    accepted[face] = self.zone.characteristic_patch_drive_on(face != 0)
                        .ok_or(NsfDgError::Input(
                            "two-sided coupling lost an original aperture control"))?;
                }
            }
            accepted.map(|drive| Some(match drive {
                NsfBoundaryDrive::IncomingCharacteristic(value) => value,
                NsfBoundaryDrive::PressureOffset(_)
                | NsfBoundaryDrive::PressureTemperature { .. } => f64::midpoint(
                    self.limits.incoming_bracket_m_s[0], self.limits.incoming_bracket_m_s[1]),
            }))
        } else {
            // Preserve the exact preexisting numerical trajectory when no
            // source aperture requested this new conditional-root strategy.
            [None; 2]
        };
        let upstream = self.upstream.as_mut().ok_or(NsfDgError::Input(
            "upstream pressure rise requires its original finite acoustic load",
        ))?;
        let main_relation = self.load.port_response().map_err(NsfDgError::Acoustic)?;
        let upstream_relation = upstream
            .load
            .port_response()
            .map_err(NsfDgError::Acoustic)?;
        let mut relations = [main_relation; 2];
        relations[other] = upstream_relation;
        let mut rises = [0.0; 2];
        rises[other] = pressure_rise_pa;
        let limits = self.limits;
        let mut evaluations = 0;
        let mut last_residual = 0.0;
        let zone = &mut self.zone;
        let mut evaluate = |incoming: [f64; 2]| -> Result<([f64; 2], [f64; 2]), NsfDgError> {
            poll(gate)?;
            // Reserve one actual original interval for final coupled restaging.
            if evaluations + 1 >= limits.maximum_evaluations {
                return Err(NsfDgError::Convergence {
                    residual_pa: last_residual,
                    evaluations,
                });
            }
            evaluations += 1;
            let trial = preview_interval(zone,
                incoming.map(NsfBoundaryDrive::IncomingCharacteristic), patch,
                reservoir_states, limits, gate)?;
            let mut residuals = [0.0; 2];
            let mut scales = [0.0; 2];
            for face in 0..2 {
                let exchange = coupled_exchange(&trial, face)?;
                let dynamic = relations[face].impedance_pa_s_m3 * exchange.normal_volume_flow_m3_s;
                residuals[face] = exchange.mean_acoustic_pressure_pa
                    - (relations[face].history_pressure_pa + dynamic + rises[face]);
                scales[face] = exchange.mean_acoustic_pressure_pa.abs()
                    + relations[face].history_pressure_pa.abs()
                    + dynamic.abs()
                    + rises[face].abs();
            }
            if !residuals.into_iter().chain(scales).all(f64::is_finite) {
                return Err(NsfDgError::Input(
                    "two-sided original pressure relation overflowed",
                ));
            }
            last_residual = residuals.into_iter().map(f64::abs).fold(0.0, f64::max);
            Ok((residuals, scales))
        };
        let mut main_seed = seeds[side];
        let (upstream_incoming, downstream_incoming) = solve_scalar(limits, seeds[other], |incoming_upstream| {
            let mut incoming = [incoming_upstream; 2];
            let (downstream, (residuals, scales)) = solve_scalar(limits, main_seed, |incoming_main| {
                incoming[side] = incoming_main;
                let value = evaluate(incoming)?;
                Ok((value.0[side], value.1[side], value))
            })?;
            // All outer trials share the original accepted fluid snapshot.
            // Rechecking this conditional root avoids repeating its entire
            // wide-bracket search when the remote face has little influence.
            if main_seed.is_some() { main_seed = Some(downstream); }
            Ok((residuals[other], scales[other], downstream))
        })?;
        // No numerical trial escapes with accepted state. Reevaluate exactly
        // the joint root once and keep all three candidates unpublished.
        let mut incoming = [upstream_incoming; 2];
        incoming[side] = downstream_incoming;
        evaluations += 1;
        let fluid = preview_interval(zone,
            incoming.map(NsfBoundaryDrive::IncomingCharacteristic), patch,
            reservoir_states, limits, gate)?;
        let dt = self.load.propagator().dt();
        let main_inlet = self.load.inlet();
        let upstream_inlet = upstream.load.inlet();
        let faces = *fluid.boundary_exchange();
        let coupled_patch = fluid.patch_exchange_on(side != 0).copied();
        let source_patch = fluid.patch_exchange_on(other != 0).copied();
        let main_face = coupled_exchange(&fluid, side)?;
        let source_face = coupled_exchange(&fluid, other)?;
        let main = self
            .load
            .preview_step(main_face.normal_volume_flow_m3_s)
            .map_err(NsfDgError::Acoustic)?;
        let source = upstream
            .load
            .preview_step(source_face.normal_volume_flow_m3_s)
            .map_err(NsfDgError::Acoustic)?;
        let main_p = main
            .output()
            .port_outputs
            .iter()
            .find(|p| p.id == main_inlet)
            .ok_or(NsfDgError::Input(
                "original main DG pressure observation disappeared",
            ))?
            .conjugate_rate;
        let source_p = source
            .output()
            .port_outputs
            .iter()
            .find(|p| p.id == upstream_inlet)
            .ok_or(NsfDgError::Input(
                "original upstream DG pressure observation disappeared",
            ))?
            .conjugate_rate;
        let main_audit = audit(
            main_face,
            main_p,
            main.inlet_work_j(),
            0.0,
            main_relation,
            self.reference,
            dt,
            limits,
        )?;
        let source_audit = audit(
            source_face,
            source_p,
            source.inlet_work_j(),
            pressure_rise_pa,
            upstream_relation,
            upstream.reference,
            dt,
            limits,
        )?;
        let flow = *fluid.energy_ledger();
        let reservoir_mass_out_kg = fluid.reservoir_mass_out_kg();
        let reservoir_flows = fluid.reservoir_flows().to_vec();
        let reservoir_exchanges = fluid.reservoir_exchanges().to_vec();
        let mass_change = fluid.grid().total_mass_kg() - fluid.previous_mass_kg();
        let energy_change = fluid.grid().total_energy_j() - fluid.previous_energy_j();
        if ![mass_change, energy_change, reservoir_mass_out_kg].iter().all(|x| x.is_finite()) {
            return Err(NsfDgError::Input(
                "two-sided resolved extensive accounts overflowed",
            ));
        }
        let terminals = main.terminal_frames().to_vec();
        let upstream_terminals = source.terminal_frames().to_vec();
        let finite_candidate = finite_reservoir::prepare(
            self.finite_reservoir.as_ref(), &fluid, control, dt, gate,
        )?;
        poll(gate)?;
        fluid.accept();
        let duct = main.accept();
        let source_duct = source.accept();
        let finite_reservoir = finite_reservoir::publish(&mut self.finite_reservoir, finite_candidate);
        self.steps += 1;
        Ok(NsfDgFrame {
            step: self.steps,
            evaluations,
            incoming_characteristic_m_s: downstream_incoming,
            drive: control,
            upstream: Some(NsfDgUpstreamFrame {
                face: other,
                incoming_characteristic_m_s: upstream_incoming,
                duct: source_duct,
                terminals: upstream_terminals,
                inlet_flow_m3_s: source_face.normal_volume_flow_m3_s,
                inlet_pressure_pa: source_p,
                source_work_j: source_audit.source_work,
                pressure_residual_pa: source_audit.pressure,
                projection_covariance_j: source_audit.covariance,
                acoustic_work_residual_j: source_audit.work,
                port_work_residual_j: source_audit.port,
                mass_projection_difference_kg: source_audit.mass,
            }),
            faces,
            coupled_patch,
            source_patch,
            reservoir_mass_out_kg,
            reservoir_flows,
            reservoir_exchanges,
            reservoir_radiation: None,
            finite_reservoir,
            flow,
            flow_mass_change_kg: mass_change,
            flow_energy_change_j: energy_change,
            duct,
            terminals,
            inlet_flow_m3_s: main_face.normal_volume_flow_m3_s,
            inlet_pressure_pa: main_p,
            pressure_residual_pa: main_audit.pressure,
            projection_covariance_j: main_audit.covariance,
            acoustic_work_residual_j: main_audit.work,
            port_work_residual_j: main_audit.port,
            mass_projection_difference_kg: main_audit.mass,
        })
    }
}

// Numerical root only: constitutive evaluation, pressure, time and work remain
// with the original owners. The bracket is a caller-declared domain, not a clamp.
pub(super) fn solve_scalar<T>(
    limits: NsfDgLimits,
    initial: Option<f64>,
    evaluate: impl FnMut(f64) -> Result<(f64, f64, T), NsfDgError>,
) -> Result<(f64, T), NsfDgError> {
    solve_in_domain(limits, limits.incoming_bracket_m_s, initial, false, evaluate)
}

pub(super) fn solve_pressure<T>(
    limits: NsfDgLimits,
    pressure_bracket_pa: [f64; 2],
    evaluate: impl FnMut(f64) -> Result<(f64, f64, T), NsfDgError>,
) -> Result<(f64, T), NsfDgError> {
    solve_in_domain(limits, pressure_bracket_pa, None, true, evaluate)
}

fn solve_in_domain<T>(
    limits: NsfDgLimits,
    mut bracket: [f64; 2],
    initial: Option<f64>,
    pressure_coordinate: bool,
    mut evaluate: impl FnMut(f64) -> Result<(f64, f64, T), NsfDgError>,
) -> Result<(f64, T), NsfDgError> {
    let seed = if let Some(initial) = initial {
        let incoming = if initial.is_finite() && initial > bracket[0] && initial < bracket[1] {
            initial
        } else { f64::midpoint(bracket[0], bracket[1]) };
        let (residual, scale, value) = evaluate(incoming)?;
        if residual.abs()
            <= limits.pressure_absolute_tolerance_pa + limits.pressure_relative_tolerance * scale
        {
            return Ok((incoming, value));
        }
        Some((incoming, residual))
    } else { None };
    let mut residuals = [0.0; 2];
    let mut weights = [0.0; 2];
    let mut previous_replaced_side = None;
    let mut incoming = bracket[0];
    let remaining = limits.maximum_evaluations - usize::from(seed.is_some());
    for evaluation in 1..=remaining {
        let (residual, scale, value) = evaluate(incoming)?;
        if residual.abs()
            <= limits.pressure_absolute_tolerance_pa + limits.pressure_relative_tolerance * scale
        {
            return Ok((incoming, value));
        }
        if evaluation == 1 {
            residuals[0] = residual;
            weights[0] = residual;
            incoming = bracket[1];
            continue;
        }
        if evaluation == 2 {
            residuals[1] = residual;
            weights[1] = residual;
            if residuals[0].is_sign_positive() == residuals[1].is_sign_positive() {
                return Err(NsfDgError::NoBracket {
                    residual_pa: residuals,
                });
            }
            // The checked interior seed narrows the same signed physical
            // bracket. Its evaluation is already charged to the shared budget.
            if let Some((seed, seed_residual)) = seed {
                let seed_side = usize::from(seed_residual.is_sign_positive() != residuals[0].is_sign_positive());
                bracket[seed_side] = seed;
                residuals[seed_side] = seed_residual;
                weights[seed_side] = seed_residual;
            }
        } else {
            let side = usize::from(residual.is_sign_positive() != residuals[0].is_sign_positive());
            bracket[side] = incoming;
            residuals[side] = residual;
            weights[side] = residual;
            // Illinois false position retains the physical signed bracket.
            // Repeatedly replacing one endpoint halves the OTHER interpolation
            // weight so an asymptotic endpoint cannot stall the scalar solve.
            // The actual residuals still determine convergence and admission.
            weights[1 - side] = if previous_replaced_side == Some(side) {
                0.5 * weights[1 - side]
            } else {
                residuals[1 - side]
            };
            previous_replaced_side = Some(side);
        }
        let proposed =
            bracket[0] - weights[0] * ((bracket[1] - bracket[0]) / (weights[1] - weights[0]));
        incoming = if proposed.is_finite() && proposed > bracket[0] && proposed < bracket[1] {
            proposed
        } else {
            f64::midpoint(bracket[0], bracket[1])
        };
        if incoming <= bracket[0] || incoming >= bracket[1] {
            return Err(if pressure_coordinate {
                NsfDgError::PressureResolution {
                    pressure_bracket_pa: bracket, residual_pa: residuals,
                }
            } else {
                NsfDgError::Resolution {
                    incoming_bracket_m_s: bracket, residual_pa: residuals,
                }
            });
        }
    }
    Err(NsfDgError::Input(
        "coupled scalar solve exhausted its local evaluation allowance",
    ))
}

pub(super) struct Audit {
    pub(super) pressure: f64,
    pub(super) covariance: f64,
    pub(super) work: f64,
    pub(super) port: f64,
    pub(super) mass: f64,
    pub(super) source_work: f64,
}
#[allow(clippy::too_many_arguments)]
pub(super) fn audit(
    face: NsfBoundaryExchange,
    pressure: f64,
    inlet_work: f64,
    rise: f64,
    relation: DgPortResponse,
    reference: GasState,
    dt: f64,
    limits: NsfDgLimits,
) -> Result<Audit, NsfDgError> {
    let q = face.normal_volume_flow_m3_s;
    let source_work = -rise * (q * dt);
    let residual = face.mean_acoustic_pressure_pa - pressure - rise;
    let scale = face.mean_acoustic_pressure_pa.abs()
        + pressure.abs()
        + rise.abs()
        + relation.history_pressure_pa.abs()
        + (relation.impedance_pa_s_m3 * q).abs();
    let port = inlet_work - pressure * (q * dt);
    let covariance = face.acoustic_pressure_work_out_j - face.mean_acoustic_pressure_pa * (q * dt);
    let work = face.acoustic_pressure_work_out_j + source_work - inlet_work;
    let work_scale = face.acoustic_pressure_work_out_j.abs() + source_work.abs() + inlet_work.abs();
    let allowance =
        limits.projection_absolute_tolerance_j + limits.projection_relative_tolerance * work_scale;
    let mass = face.mass_out_kg - reference.density * (q * dt);
    if ![
        pressure,
        inlet_work,
        source_work,
        residual,
        scale,
        port,
        covariance,
        work,
        work_scale,
        allowance,
        mass,
    ]
    .iter()
    .all(|x| x.is_finite())
        || residual.abs()
            > limits.pressure_absolute_tolerance_pa + limits.pressure_relative_tolerance * scale
        || port.abs() > 1e-10 * (inlet_work.abs() + scale * q.abs() * dt).max(1e-300)
    {
        return Err(NsfDgError::Input(
            "original two-sided DG pressure or integrated work disagrees with its source relation",
        ));
    }
    if work.abs() > allowance {
        return Err(NsfDgError::Projection {
            residual_j: work,
            allowance_j: allowance,
        });
    }
    Ok(Audit {
        pressure: residual,
        covariance,
        work,
        port,
        mass,
        source_work,
    })
}
