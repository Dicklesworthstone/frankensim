//! Physical boundary drives on the original continuous thermal impedance.
//!
//! At a fixed operating point, inlet enthalpy contributes mdot * delta_h.
//! A COMMON fractional change of every prescribed mass rate contributes the
//! original, reference-reduced net advective watts times that fraction. This
//! preserves fixed cell mass and any admitted stationary species balance.
//! Independent valve rates would require a new continuity/composition solve.
//!
//! Ordinary ambient/reservoir temperatures act through their original
//! conductance at fixed cell energy. A complete nonlinear boundary law with
//! no bath-temperature derivative is explicitly unavailable as a bath drive;
//! its original state tangent still participates in heat and flow responses.
//! All inputs are peak phasors using exp(-i omega t). Amplitudes must remain
//! in the admitted linear branch. No physical clock or energy credit advances.

use super::*;
use crate::enthalpy::advection::MassBoundary;

/// The two original thermal boundary channels on a cell.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum ThermalBathChannel {
    Ambient,
    Reservoir,
}

/// Frozen original bath; None requires the boundary law's own input derivative.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalBathPort {
    pub cell: CellId,
    pub channel: ThermalBathChannel,
    pub temperature_k: f64,
    /// Inward watts per kelvin, with the original cell energy held fixed.
    pub inward_power_per_kelvin_w_k: Option<f64>,
    row: usize,
}

/// A real inlet at the same physical clock as the prepared equilibrium tangent.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalInletPort {
    pub boundary_id: u64,
    pub cell: CellId,
    pub mass_rate_kg_s: f64,
    pub specific_enthalpy_j_kg: f64,
    row: usize,
}

/// Original tangent plus immutable physical input columns, with O(cells+inlets)
/// added storage. No state or frequency perturbation is used to build them.
pub struct DrivenThermalImpedanceModel<'m> {
    original: ThermalImpedanceModel,
    baths: Vec<ThermalBathPort>,
    inlets: Vec<ThermalInletPort>,
    inlet_laws: Vec<&'m dyn ThermalLaw>,
    advective_power_w: Vec<f64>,
    has_advection: bool,
}

/// One native LU, shared by all combinations of physical input phasors.
pub struct DrivenThermalImpedanceFactor<'s, 'm> {
    original: ThermalImpedanceFactor<'s>,
    model: &'s DrivenThermalImpedanceModel<'m>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct DrivenThermalImpedanceResponse {
    /// Original cell order, caloric response and full original-row residuals.
    pub field: ThermalImpedanceResponse,
    /// Total linearized inward drive in watts, before state feedback.
    /// This is the original equation RHS, not an additional heat credit.
    pub forcing_w: Vec<C64>,
}

impl ImplicitConductionChain<'_> {
    /// Compile physical bath, inlet-enthalpy and common-pump inputs on the
    /// SAME original operating tangent. The numerical payload bound is
    /// 64*n*n + 512*n + 128*k bytes for n cells and k inlets, excluding the
    /// original graph/source maps and allocator overhead. Input assembly adds
    /// a bounded 64*(n+k) work allowance to the original dense-work bound.
    ///
    /// Inlet laws are retained to admit the actual peak enthalpy extrema.
    /// Inlet composition, pressure convention and enthalpy datum stay fixed.
    /// The same fractional flow change applies to all original connections;
    /// a stopped pump remains stopped. Pump reversal and independent valves
    /// are outside this local fixed-mass linear response.
    pub fn prepare_driven_thermal_impedance<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        budget: ThermalImpedanceBudget,
        maximum_inlets: usize,
        mut checkpoint: impl FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<DrivenThermalImpedanceModel<'m>> {
        checkpoint()?;
        let mut inlet_count = 0_usize;
        if let Some(advection) = &self.advection {
            for boundary in advection.flow.boundaries() {
                checkpoint()?;
                if matches!(boundary, MassBoundary::Inlet { .. }) {
                    inlet_count = inlet_count.checked_add(1)
                        .ok_or_else(|| invalid("thermal inlet count overflow"))?;
                }
            }
        }
        let n = self.ids.len() as u128;
        if inlet_count > maximum_inlets
            || 64*n*n + 512*n + 128*inlet_count as u128
                > budget.maximum_dense_bytes as u128
            || n*n*(n+16) + 64*(n+inlet_count as u128)
                > budget.maximum_dense_work as u128
        {
            return Err(invalid("driven thermal response exceeds inlet/storage budget"));
        }
        let (original, endpoint) = self.prepare_thermal_impedance_and_endpoint(
            thermal, budget, &mut checkpoint,
        )?;
        let mut covered = vec![0.0; self.ids.len()];
        for (i, id) in self.ids.iter().enumerate() {
            checkpoint()?;
            covered[i] = thermal.cells[id].boundary.reservoir.map_or(0.0, |r| r.area_m2);
        }
        for link in thermal.links.values() {
            checkpoint()?;
            let areas = thermal.contact_areas(link);
            for (id, area) in [(link.a, areas[0]), (link.b, areas[1])] {
                let row = self.index[&id];
                covered[row] = finite(covered[row] + area)?;
            }
        }
        let mut baths = Vec::with_capacity(2*self.ids.len());
        for (row, &id) in self.ids.iter().enumerate() {
            checkpoint()?;
            let cell = &thermal.cells[&id];
            let exposed_area = finite(cell.geometry.surface_area_m2 - covered[row])?;
            if exposed_area < 0.0 {
                return Err(DepositionError::Arithmetic.into());
            }
            let ambient_gain = if self.ambient_laws.contains_key(&id)
                || self.ambient_exchanges.contains_key(&id)
            {
                None
            } else {
                Some(endpoint.boundary_values.ambient(
                    id, exposed_area, cell.boundary.htc_w_m2_k,
                )?)
            };
            baths.push(ThermalBathPort {
                cell: id,
                channel: ThermalBathChannel::Ambient,
                temperature_k: boundary::temperature_at(
                    &self.ambient_histories, id, cell.boundary.ambient_k, thermal.time_s,
                )?,
                inward_power_per_kelvin_w_k: ambient_gain,
                row,
            });
            if let Some(reservoir) = cell.boundary.reservoir {
                let gain = if self.reservoir_laws.contains_key(&id)
                    || self.reservoir_exchanges.contains_key(&id)
                {
                    None
                } else {
                    Some(*endpoint.boundary_values.reservoir_w_k.get(&id)
                        .ok_or(DepositionError::Arithmetic)?)
                };
                baths.push(ThermalBathPort {
                    cell: id,
                    channel: ThermalBathChannel::Reservoir,
                    temperature_k: boundary::temperature_at(
                        &self.reservoir_histories, id, reservoir.temperature_k, thermal.time_s,
                    )?,
                    inward_power_per_kelvin_w_k: gain,
                    row,
                });
            }
        }
        let mut inlets = Vec::with_capacity(inlet_count);
        let mut inlet_laws = Vec::with_capacity(inlet_count);
        if let Some(advection) = &self.advection {
            let stage = advection.stage(thermal.time_s, &mut checkpoint)?;
            for boundary in stage.flow().boundaries() {
                checkpoint()?;
                if let MassBoundary::Inlet {
                    id, cell, kg_s, material_id, specific_h_j_kg,
                } = *boundary {
                    inlets.push(ThermalInletPort {
                        boundary_id: id, cell, mass_rate_kg_s: kg_s,
                        specific_enthalpy_j_kg: specific_h_j_kg, row: self.index[&cell],
                    });
                    inlet_laws.push(*thermal.materials.get(&material_id)
                        .ok_or(DepositionError::Arithmetic)?);
                }
            }
        }
        // The original boundary list is in stable connection-ID order. The
        // scheduled sampler changes only rates/h; it never reorders its IDs.
        let mut powers: Vec<Sum> = (0..self.ids.len()).map(|_| Sum::default()).collect();
        if let Some(flow) = &endpoint.stationary_advection {
            // The stream powers come beside the transport (flow.sources), and
            // equilibrium.rs reads them the same way; sources_w is empty here.
            flow.sources.visit(&mut |id, watts| {
                checkpoint()?;
                powers[self.index[&id]].add(watts)
            })?;
            for (&id, &power) in &flow.transport.heat.bulk.sources_w {
                checkpoint()?;
                powers[self.index[&id]].add(power)?;
            }
            for transfer in &flow.transport.transfers {
                checkpoint()?;
                powers[self.index[&transfer.a]].add(-transfer.power_w)?;
                powers[self.index[&transfer.b]].add(transfer.power_w)?;
            }
        }
        let mut advective_power_w = Vec::with_capacity(self.ids.len());
        for power in powers {
            checkpoint()?;
            advective_power_w.push(power.total()?);
        }
        checkpoint()?;
        Ok(DrivenThermalImpedanceModel {
            original, baths, inlets, inlet_laws, advective_power_w,
            has_advection: self.advection.is_some(),
        })
    }
}

impl<'m> DrivenThermalImpedanceModel<'m> {
    pub fn original(&self) -> &ThermalImpedanceModel { &self.original }
    pub fn baths(&self) -> &[ThermalBathPort] { &self.baths }
    pub fn inlets(&self) -> &[ThermalInletPort] { &self.inlets }

    /// Original per-cell advective watts: the exact common fractional-rate
    /// input column. Uses the native datum-reduced fluxes, including their
    /// admitted continuity-roundoff convention, never gross h*mdot subtraction.
    pub fn common_flow_input_w(&self) -> &[f64] { &self.advective_power_w }

    /// Assemble physical inputs deterministically in original ID order.
    ///
    /// All original heat/source feedback is already in the tangent. heat_w
    /// contains only ADDITIONAL prescribed heat. Bath maps use (cell, channel);
    /// inlet maps use original connection IDs. Missing entries mean zero.
    /// Source domains are checked at inlet h +/- |delta_h| and bath T remains
    /// positive. These checks do not certify finite-amplitude linear accuracy.
    pub fn forcing(
        &self,
        heat_w: &[C64],
        bath_temperature_k: &BTreeMap<(CellId, ThermalBathChannel), C64>,
        inlet_enthalpy_j_kg: &BTreeMap<u64, C64>,
        common_flow_fraction: C64,
        mut checkpoint: impl FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<Vec<C64>> {
        checkpoint()?;
        let n = self.original.ids.len();
        if heat_w.len() != n || bath_temperature_k.len() > self.baths.len()
            || inlet_enthalpy_j_kg.len() > self.inlets.len()
        {
            return Err(invalid("driven thermal input dimensions exceed original ports"));
        }
        complex(common_flow_fraction)?;
        if !self.has_advection && common_flow_fraction != C64::ZERO {
            return Err(invalid("common flow drive requires original prescribed advection"));
        }
        // |delta rate/rate| < 1 keeps every directed edge on the same branch.
        if common_flow_fraction.abs() >= 1.0 {
            return Err(invalid("thermal flow phasor reaches a zero-flow/reversal corner"));
        }
        let mut real: Vec<Sum> = (0..n).map(|_| Sum::default()).collect();
        let mut imag: Vec<Sum> = (0..n).map(|_| Sum::default()).collect();
        for (row, &heat) in heat_w.iter().enumerate() {
            checkpoint()?;
            complex(heat)?;
            let flow = drive_power(common_flow_fraction, self.advective_power_w[row])?;
            real[row].add(heat.re)?;
            real[row].add(flow.re)?;
            imag[row].add(heat.im)?;
            imag[row].add(flow.im)?;
        }
        for (&key, &delta) in bath_temperature_k {
            checkpoint()?;
            complex(delta)?;
            let index = self.baths.binary_search_by_key(&key, |port| (port.cell, port.channel))
                .map_err(|_| invalid("unknown original thermal bath drive"))?;
            let port = self.baths[index];
            if delta == C64::ZERO { continue; }
            if !positive(port.temperature_k - delta.abs()) {
                return Err(DepositionError::MaterialDomain.into());
            }
            let gain = port.inward_power_per_kelvin_w_k
                .ok_or_else(|| invalid("original boundary law has no bath-input derivative"))?;
            let power = drive_power(delta, gain)?;
            real[port.row].add(power.re)?;
            imag[port.row].add(power.im)?;
        }
        for (&id, &delta) in inlet_enthalpy_j_kg {
            checkpoint()?;
            complex(delta)?;
            let index = self.inlets.binary_search_by_key(&id, |port| port.boundary_id)
                .map_err(|_| invalid("unknown original thermal inlet drive"))?;
            let port = self.inlets[index];
            if delta == C64::ZERO { continue; }
            let amplitude = finite(delta.abs())?;
            let law = self.inlet_laws[index];
            for h in [
                finite(port.specific_enthalpy_j_kg - amplitude)?,
                finite(port.specific_enthalpy_j_kg + amplitude)?,
            ] {
                checkpoint()?;
                if !positive(law.temperature_k(h)?) {
                    return Err(DepositionError::MaterialDomain.into());
                }
            }
            let power = drive_power(delta, port.mass_rate_kg_s)?;
            real[port.row].add(power.re)?;
            imag[port.row].add(power.im)?;
        }
        let mut result = Vec::with_capacity(n);
        for (re, im) in real.into_iter().zip(imag) {
            checkpoint()?;
            result.push(C64::new(re.total()?, im.total()?));
        }
        checkpoint()?;
        Ok(result)
    }

    pub fn factor(
        &self,
        omega_rad_s: f64,
        checkpoint: impl FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<DrivenThermalImpedanceFactor<'_, 'm>> {
        Ok(DrivenThermalImpedanceFactor {
            original: self.original.factor(omega_rad_s, checkpoint)?,
            model: self,
        })
    }
}

impl DrivenThermalImpedanceFactor<'_, '_> {
    pub fn omega_rad_s(&self) -> f64 { self.original.omega_rad_s() }

    pub fn response(
        &self,
        heat_w: &[C64],
        bath_temperature_k: &BTreeMap<(CellId, ThermalBathChannel), C64>,
        inlet_enthalpy_j_kg: &BTreeMap<u64, C64>,
        common_flow_fraction: C64,
        mut checkpoint: impl FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<DrivenThermalImpedanceResponse> {
        let forcing_w = self.model.forcing(
            heat_w, bath_temperature_k, inlet_enthalpy_j_kg, common_flow_fraction,
            &mut checkpoint,
        )?;
        let field = self.original.response(&forcing_w, &mut checkpoint)?;
        checkpoint()?;
        Ok(DrivenThermalImpedanceResponse { field, forcing_w })
    }
}

fn drive_power(value: C64, gain: f64) -> Result<C64> {
    let result = complex(value.scale(gain))?;
    if gain != 0.0 && ((value.re != 0.0 && result.re == 0.0)
        || (value.im != 0.0 && result.im == 0.0))
    {
        return Err(invalid("physical thermal input power is not representable"));
    }
    Ok(result)
}
