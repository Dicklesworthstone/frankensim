//! Original coupled intervals and one shared controller for hydraulic storage.

use super::*;
use crate::laminar_enthalpy::adaptive::{
    LaminarAdaptiveCandidate, LaminarAdaptiveControls, LaminarAdaptiveError,
    LaminarAdaptiveFailure, LaminarAdaptiveReport, LaminarAdaptiveWork,
};
use fs_math::dd::Dd;

/// Physical circuit accuracy, independent of the original algebraic tolerances.
#[derive(Clone, Copy, Debug)]
pub struct LaminarWindkesselAdaptiveControls {
    pub coupled: LaminarAdaptiveControls,
    pub pressure_law_tolerance_pa: f64,
    pub flow_law_tolerance_m2_s: f64,
    /// Compliance and added port-pressure local error [Pa]. Reservoir pressure
    /// offsets never contribute to the relative error allowance.
    pub pressure_absolute_tolerance_pa: f64,
    pub flow_absolute_tolerance_m2_s: f64,
    /// Signed and resolved-magnitude integral of native Q over physical depth [m³].
    pub transport_absolute_tolerance_m3: f64,
    /// Each source-owned C and L storage error, independently, after depth [J].
    pub storage_absolute_tolerance_j: f64,
    /// Each proximal/distal/reservoir/port transfer, signed and magnitude [J].
    pub transfer_absolute_tolerance_j: f64,
    /// Independent per-port AND aggregate external energy closure [J].
    pub window_tolerance_j: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct LaminarWindkesselTimeEstimate {
    pub maximum_error_ratio: f64,
    pub maximum_pressure_difference_pa: f64,
    pub maximum_flow_difference_m2_s: f64,
    pub maximum_transport_difference_m3: f64,
    pub maximum_storage_difference_j: f64,
    pub maximum_transfer_difference_j: f64,
}

/// Accepted-only physical-depth work. Numerical BDF defect is separate and is
/// included only in the external circuit identity, never deposited as heat.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LaminarWindkesselWindowBalance {
    /// Sum Q*dt*depth, positive from fluid to this external port [m³].
    pub transported_volume_m3: f64,
    /// Sum |Q|*dt*depth over accepted fine intervals [m³].
    pub absolute_transported_volume_m3: f64,
    pub compliance_energy_change_j: f64,
    pub inertance_energy_change_j: f64,
    pub proximal_dissipation_j: f64,
    pub distal_dissipation_j: f64,
    pub reservoir_work_j: f64,
    pub port_work_j: f64,
    pub bdf_energy_defect_j: f64,
    /// Port - reservoir - both losses - both storage changes - BDF defect.
    pub balance_residual_j: f64,
}

pub type LaminarWindkesselAdaptiveCandidate<'m, G, T> =
    LaminarAdaptiveCandidate<'m, G, LaminarWindkesselExchange<T>>;

impl<'s, 'mesh, 'source, F: FluidTransport + ?Sized> LaminarEnthalpy<'s, 'mesh, 'source, F> {
    /// Advance the original fluid, caloric inventory and held external circuit
    /// through one bounded, atomic adaptive window. The existing BDF1 coarse /
    /// two-fine controller owns all work ceilings and accept/reject decisions.
    ///
    /// Coarse and first-half solves borrow the same accepted port history. The
    /// second half borrows the first half's actual compliance and velocity.
    /// A rejection discards all three candidates; no mutable circuit state is
    /// advanced by a thermal probe or discarded comparison. Only accepted fine
    /// receipts enter the final external energy account. Continuation uses the
    /// last retained step's compliance coordinates with the returned flow/time.
    ///
    /// Hydraulic pressure/flow/storage, integrated port throughput and each
    /// signed/magnitude work channel
    /// join the original velocity/caloric/transfer norm. BDF defects remain
    /// numerical. The final original C/L energy change is measured from paired
    /// endpoint states, independently of summed per-step receipts. No external
    /// resistor heat, new solver, extrapolated field or passivity claim is added.
    /// Reservoir references and original mechanical/circuit data remain held;
    /// this entry adds no scheduled drive events or Boussinesq/load combination.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_windkessel_adaptive<'m, G: ThermalFluidTransport, T, E>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        input: LaminarWindkesselInput<'_>,
        controls: LaminarWindkesselAdaptiveControls,
        resolve_material: impl Fn(MaterialFluxRheologyPoint, f64, &F) -> Result<G, MaterialError>,
        exchange: impl Fn(
            &Cx<'_>,
            &DepositionSimulation<'m>,
            &MaterialStokesWindkesselSolution<
                MaterialNavierStokesRheologyStep<LaminarThermalMaterial<G>>,
            >,
        ) -> Result<LaminarHeatExchange<T>, E>,
    ) -> Result<LaminarWindkesselAdaptiveCandidate<'m, G, T>, LaminarAdaptiveError<E>> {
        let empty = LaminarAdaptiveWork::default();
        checkpoint(cx).map_err(|e| failure(e, empty))?;
        if ![
            controls.pressure_absolute_tolerance_pa,
            controls.flow_absolute_tolerance_m2_s,
            controls.transport_absolute_tolerance_m3,
            controls.storage_absolute_tolerance_j,
            controls.transfer_absolute_tolerance_j,
            controls.window_tolerance_j,
        ]
        .iter()
        .all(|x| x.is_finite() && *x > 0.0)
            || !nonnegative(controls.pressure_law_tolerance_pa)
            || !nonnegative(controls.flow_law_tolerance_m2_s)
            || input.ports.is_empty()
            || input.ports.len() != input.loads.len()
            || input.ports.len() != input.previous_compliance_pressures_pa.len()
            || input
                .ports
                .iter()
                .try_fold(0usize, |n, p| n.checked_add(p.len()))
                .is_none_or(|n| n > self.system.mesh.edges.len())
        {
            return Err(LaminarAdaptiveError {
                failure: LaminarAdaptiveFailure::Invalid(
                    "finite circuit tolerances and complete bounded paired port history required",
                ),
                work: empty,
            });
        }
        // The original mean flux is canonical; never infer old Q from a future
        // energy defect or replace it with an independent caller flow history.
        let rates = self
            .system
            .edge_volume_transport(input.fluid.history.velocity_m_s)
            .map_err(|e| failure(LaminarEnthalpyError::FaceTransport(e), empty))?;
        let mut initial_q = Vec::with_capacity(input.ports.len());
        for port in input.ports {
            let mut q = Dd::from_f64(0.0);
            for &edge in *port {
                checkpoint(cx).map_err(|e| failure(e, empty))?;
                let rate = rates
                    .get(edge)
                    .filter(|_| self.system.mesh.is_boundary(edge))
                    .ok_or_else(|| {
                        failure(
                            LaminarEnthalpyError::Invalid("original boundary port edge required"),
                            empty,
                        )
                    })?;
                q = q + Dd::from_f64(rate.net_m2_s());
            }
            initial_q.push(q.to_f64());
        }
        drop(rates);
        let mut solved = self.advance_adaptive_with_port_history(
            cx,
            thermal,
            target_s,
            input.fluid.history,
            &[],
            controls.coupled,
            |old,
             query,
             history,
             previous: Option<&LaminarEnthalpyReport<LaminarWindkesselExchange<T>>>,
             step,
             guess| {
                let pressure =
                    previous.map(|report| report.exchange.circuit.compliance_pressures_pa());
                self.advance_windkessel(
                    cx,
                    old,
                    query.time_s,
                    LaminarWindkesselInput {
                        fluid: LaminarEnthalpyInput {
                            history,
                            ..input.fluid
                        },
                        previous_compliance_pressures_pa: pressure
                            .as_deref()
                            .unwrap_or(input.previous_compliance_pressures_pa),
                        ..input
                    },
                    LaminarWindkesselControls {
                        coupled: step,
                        absolute_pressure_tolerance_pa: controls.pressure_law_tolerance_pa,
                        absolute_flow_tolerance_m2_s: controls.flow_law_tolerance_m2_s,
                    },
                    guess,
                    &resolve_material,
                    &exchange,
                )
            },
            |previous, coarse, first, fine, correction| {
                compare(
                    cx,
                    previous.map(|r| &r.exchange.circuit),
                    &coarse.report().exchange.circuit,
                    &first.report().exchange.circuit,
                    &fine.report().exchange.circuit,
                    input,
                    &initial_q,
                    controls,
                    correction,
                )
                .map(Some)
            },
        )?;
        let balances = window_account(
            cx,
            &solved.report,
            input,
            &initial_q,
            self.binding.depth_m(),
            controls.window_tolerance_j,
        )
        .map_err(|e| failure(e, solved.report.work))?;
        checkpoint(cx).map_err(|e| failure(e, solved.report.work))?;
        solved.report.hydraulic = Some(balances);
        Ok(solved)
    }
}

#[allow(clippy::too_many_arguments)]
fn compare<E>(
    cx: &Cx<'_>,
    previous: Option<&LaminarWindkesselReport>,
    coarse: &LaminarWindkesselReport,
    first: &LaminarWindkesselReport,
    fine: &LaminarWindkesselReport,
    input: LaminarWindkesselInput<'_>,
    initial_q: &[f64],
    controls: LaminarWindkesselAdaptiveControls,
    correction: f64,
) -> Result<LaminarWindkesselTimeEstimate, LaminarEnthalpyError<E>> {
    let mut result = LaminarWindkesselTimeEstimate::default();
    let ratio = |difference, change, absolute| {
        crate::laminar_enthalpy::adaptive::estimate::ratio(
            difference,
            change,
            absolute,
            controls.coupled.relative_tolerance,
            correction,
        )
    };
    for (i, load) in input.loads.iter().enumerate() {
        checkpoint(cx)?;
        let x0 = previous.map_or(input.previous_compliance_pressures_pa[i], |p| {
            p.balances[i].compliance_pressure_pa
        });
        let q0 = previous.map_or(initial_q[i], |p| p.outward_volume_rates_m2_s[i]);
        let xc = coarse.balances[i].compliance_pressure_pa;
        let xf = fine.balances[i].compliance_pressure_pa;
        let qc = coarse.outward_volume_rates_m2_s[i];
        let qf = fine.outward_volume_rates_m2_s[i];
        // Subtract the held reservoir and start's resistive/compliance state.
        // Inertance is the actual endpoint response; no initial acceleration
        // or pressure derivative is invented for the relative scale.
        let pressure_origin = binding::checked_sum(&[
            load.reservoir_pressure_pa,
            x0,
            load.proximal_resistance_pa_s_per_m2 * q0,
        ])?;
        for (a, b, origin) in [
            (xc, xf, x0),
            (
                coarse.pressure_offsets_pa[i],
                fine.pressure_offsets_pa[i],
                pressure_origin,
            ),
        ] {
            let difference = (Dd::from_f64(b) - Dd::from_f64(a)).to_f64().abs();
            let change = (a - origin).abs().max((b - origin).abs());
            result.maximum_pressure_difference_pa =
                result.maximum_pressure_difference_pa.max(difference);
            result.maximum_error_ratio = result.maximum_error_ratio.max(ratio(
                difference,
                change,
                controls.pressure_absolute_tolerance_pa,
            )?);
        }
        let difference = (Dd::from_f64(qf) - Dd::from_f64(qc)).to_f64().abs();
        result.maximum_flow_difference_m2_s = result.maximum_flow_difference_m2_s.max(difference);
        result.maximum_error_ratio = result.maximum_error_ratio.max(ratio(
            difference,
            (qc - q0).abs().max((qf - q0).abs()),
            controls.flow_absolute_tolerance_m2_s,
        )?);
        for (coefficient, a, b, origin) in [
            (load.compliance_m2_per_pa, xc, xf, x0),
            (load.inertance_pa_s2_per_m2, qc, qf, q0),
        ] {
            let difference = energy_change(coefficient, b, a, fine.depth_m)?.abs();
            let change = energy_change(coefficient, a, origin, fine.depth_m)?
                .abs()
                .max(energy_change(coefficient, b, origin, fine.depth_m)?.abs());
            result.maximum_storage_difference_j =
                result.maximum_storage_difference_j.max(difference);
            result.maximum_error_ratio = result.maximum_error_ratio.max(ratio(
                difference,
                change,
                controls.storage_absolute_tolerance_j,
            )?);
        }
        let a = transfers(coarse, i)?;
        let b = transfers(first, i)?;
        let c = transfers(fine, i)?;
        for channel in 0..5 {
            let (a, b, c, allowance) = if channel < 4 {
                (
                    a[channel],
                    b[channel],
                    c[channel],
                    controls.transfer_absolute_tolerance_j,
                )
            } else {
                (
                    integrated(coarse.outward_volume_rates_m2_s[i], coarse)?,
                    integrated(first.outward_volume_rates_m2_s[i], first)?,
                    integrated(fine.outward_volume_rates_m2_s[i], fine)?,
                    controls.transport_absolute_tolerance_m3,
                )
            };
            for (coarse, fine) in [
                (a, b + c),
                (
                    Dd::from_f64(a.to_f64().abs()),
                    Dd::from_f64(b.to_f64().abs()) + Dd::from_f64(c.to_f64().abs()),
                ),
            ] {
                let difference = (fine - coarse).to_f64().abs();
                if channel < 4 {
                    result.maximum_transfer_difference_j =
                        result.maximum_transfer_difference_j.max(difference);
                } else {
                    result.maximum_transport_difference_m3 =
                        result.maximum_transport_difference_m3.max(difference);
                }
                result.maximum_error_ratio = result.maximum_error_ratio.max(ratio(
                    difference,
                    coarse.to_f64().abs().max(fine.to_f64().abs()),
                    allowance,
                )?);
            }
        }
    }
    Ok(result)
}

fn transfers<E>(
    report: &LaminarWindkesselReport,
    port: usize,
) -> Result<[Dd; 5], LaminarEnthalpyError<E>> {
    let b = report.balances[port];
    let mut result = [Dd::from_f64(0.0); 5];
    for (out, power) in result.iter_mut().zip([
        b.proximal_dissipation_w_per_m,
        b.distal_dissipation_w_per_m,
        b.reservoir_power_w_per_m,
        b.port_power_w_per_m,
        b.bdf_energy_defect_w_per_m,
    ]) {
        *out = integrated(power, report)?;
    }
    Ok(result)
}

fn integrated<E>(
    rate: f64,
    report: &LaminarWindkesselReport,
) -> Result<Dd, LaminarEnthalpyError<E>> {
    Ok(Dd::from_f64(
        binding::product(&[rate.abs(), report.duration_s, report.depth_m])?.copysign(rate),
    ))
}

fn energy_change<E>(
    coefficient: f64,
    next: f64,
    previous: f64,
    depth: f64,
) -> Result<f64, LaminarEnthalpyError<E>> {
    // Difference of squares BEFORE dimensional multiplication retains small
    // changes on a large original C/L storage; no endpoint energy subtraction.
    let difference = (Dd::from_f64(next) - Dd::from_f64(previous)).to_f64();
    let sum = (Dd::from_f64(next) + Dd::from_f64(previous)).to_f64();
    let magnitude = binding::product(&[0.5, coefficient, difference.abs(), sum.abs(), depth])?;
    Ok(magnitude * difference.signum() * sum.signum())
}

fn window_account<T, E>(
    cx: &Cx<'_>,
    report: &LaminarAdaptiveReport<LaminarWindkesselExchange<T>>,
    input: LaminarWindkesselInput<'_>,
    initial_q: &[f64],
    depth: f64,
    tolerance: f64,
) -> Result<Vec<LaminarWindkesselWindowBalance>, LaminarEnthalpyError<E>> {
    let final_circuit = &report
        .steps
        .last()
        .ok_or(LaminarEnthalpyError::Invalid("empty hydraulic window"))?
        .report
        .exchange
        .circuit;
    let mut result = Vec::with_capacity(input.loads.len());
    let mut global_residual = Dd::from_f64(0.0);
    let mut global_scale = Dd::from_f64(0.0);
    for (i, load) in input.loads.iter().enumerate() {
        checkpoint(cx)?;
        let compliance = energy_change(
            load.compliance_m2_per_pa,
            final_circuit.balances[i].compliance_pressure_pa,
            input.previous_compliance_pressures_pa[i],
            depth,
        )?;
        let inertance = energy_change(
            load.inertance_pa_s2_per_m2,
            final_circuit.outward_volume_rates_m2_s[i],
            initial_q[i],
            depth,
        )?;
        let mut total = [Dd::from_f64(0.0); 5];
        let mut volume = Dd::from_f64(0.0);
        let mut absolute_volume = Dd::from_f64(0.0);
        for step in &report.steps {
            checkpoint(cx)?;
            let circuit = &step.report.exchange.circuit;
            let amount = integrated(circuit.outward_volume_rates_m2_s[i], circuit)?;
            volume = volume + amount;
            absolute_volume = absolute_volume + Dd::from_f64(amount.to_f64().abs());
            for (sum, amount) in total
                .iter_mut()
                .zip(transfers(&step.report.exchange.circuit, i)?)
            {
                *sum = *sum + amount;
            }
        }
        if !volume.to_f64().is_finite() || !absolute_volume.to_f64().is_finite() {
            return Err(DepositionError::Arithmetic.into());
        }
        let residual = total[3]
            - total[2]
            - total[0]
            - total[1]
            - Dd::from_f64(compliance)
            - Dd::from_f64(inertance)
            - total[4];
        let scale = binding::checked_sum(&[
            compliance.abs(),
            inertance.abs(),
            total[0].to_f64().abs(),
            total[1].to_f64().abs(),
            total[2].to_f64().abs(),
            total[3].to_f64().abs(),
            total[4].to_f64().abs(),
        ])?;
        check_balance(residual.to_f64(), scale, tolerance)?;
        global_residual = global_residual + residual;
        global_scale = global_scale + Dd::from_f64(scale);
        result.push(LaminarWindkesselWindowBalance {
            transported_volume_m3: volume.to_f64(),
            absolute_transported_volume_m3: absolute_volume.to_f64(),
            compliance_energy_change_j: compliance,
            inertance_energy_change_j: inertance,
            proximal_dissipation_j: total[0].to_f64(),
            distal_dissipation_j: total[1].to_f64(),
            reservoir_work_j: total[2].to_f64(),
            port_work_j: total[3].to_f64(),
            bdf_energy_defect_j: total[4].to_f64(),
            balance_residual_j: residual.to_f64(),
        });
    }
    check_balance(global_residual.to_f64(), global_scale.to_f64(), tolerance)?;
    Ok(result)
}

fn check_balance<E>(
    residual: f64,
    scale: f64,
    absolute: f64,
) -> Result<(), LaminarEnthalpyError<E>> {
    let allowance = absolute + 64.0 * f64::EPSILON * scale;
    if !residual.is_finite() || !allowance.is_finite() {
        return Err(DepositionError::Arithmetic.into());
    }
    if residual.abs() > allowance {
        return Err(DepositionError::EnergyBalance(residual).into());
    }
    Ok(())
}

fn failure<E>(
    error: LaminarEnthalpyError<E>,
    work: LaminarAdaptiveWork,
) -> LaminarAdaptiveError<E> {
    LaminarAdaptiveError {
        failure: LaminarAdaptiveFailure::Native {
            error,
            endpoint_allowance: 0,
        },
        work,
    }
}
