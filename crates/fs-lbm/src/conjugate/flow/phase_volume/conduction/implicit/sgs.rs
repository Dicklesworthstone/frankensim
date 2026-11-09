//! SGS heat diffusion on the existing implicit fluid--fluid thermal links.
//!
//! Each side retains its original k(h) and caloric Cp(h). Actual original M/V
//! converts accepted native nu_t into k_t=rho Cp nu_t/Pr_t before the two
//! half-cell resistances are combined. Newton differentiates the complete heat
//! flux, including original k and Cp tangents, at fixed native nu_t and geometry.
//! The outer thermal/mechanical owner separately closes that frozen closure.
//! Baths, finite-solid contacts and walls keep their molecular transport law;
//! this module supplies no thermal wall function or modeled dissipation heat.

use super::*;
use fs_conduction::enthalpy::ContactPatch;
use fs_conduction::enthalpy::conduction::implicit::link::{
    ImplicitLinkLaw, ImplicitLinkState,
};

#[derive(Clone, Copy)]
pub(super) struct SgsThermalField<'a> {
    pub(super) eddy_viscosity_m2_s: &'a [f64],
    pub(super) turbulent_prandtl_number: f64,
}

pub(super) fn admit_prandtl_number(value: f64) -> Result<(), ChtError> {
    if !value.is_finite() || value <= 0.0 {
        return Err(refuse("SGS heat requires a positive finite turbulent Prandtl number"));
    }
    Ok(())
}

impl SgsThermalField<'_> {
    fn admit(self, domain: &VoxelDomain, mass: &[f64], gate: &CancelGate) -> Result<(), ChtError> {
        poll(gate)?;
        admit_prandtl_number(self.turbulent_prandtl_number)?;
        let n = domain.cell_count();
        if mass.len() != n || self.eddy_viscosity_m2_s.len() != n {
            return Err(refuse("SGS heat requires the complete original mass and native closure fields"));
        }
        for (cell, &eddy) in self.eddy_viscosity_m2_s.iter().enumerate() {
            poll_every(gate, cell)?;
            if !eddy.is_finite() || eddy < 0.0 || (!domain.is_fluid(cell) && eddy != 0.0) {
                return Err(refuse("SGS heat needs finite nonnegative, zero-solid native eddy viscosity"));
            }
            if domain.is_fluid(cell) {
                finite_positive("flow.sgs_heat_mass", mass[cell])?;
                finite_positive("flow.sgs_heat_density", mass[cell] / domain.volume(cell))?;
                self.scale(domain, mass, cell)?;
            }
        }
        poll(gate)
    }

    fn scale(self, domain: &VoxelDomain, mass: &[f64], cell: usize) -> Result<f64, ChtError> {
        let eddy = self.eddy_viscosity_m2_s[cell];
        if eddy == 0.0 {
            return Ok(0.0);
        }
        let scale = (mass[cell] / domain.volume(cell)) * (eddy / self.turbulent_prandtl_number);
        finite_positive("flow.sgs_heat_density_diffusivity", scale)?;
        Ok(scale)
    }
}

pub(super) struct SgsThermalLink<'a, 'm> {
    source: &'a PreparedPhaseConductivity<'m>,
    density_diffusivity: [f64; 2],
}

/// Source Cp and its analytic enthalpy tangent on the SAME one-sided interval
/// as the compiled conductivity. The prepared upper bound selects left; all
/// other exact caloric knots select right. No chord replaces affine Cp(T).
fn sensible_capacity(
    source: &PreparedPhaseConductivity<'_>,
    sample: PhaseConductivitySample,
) -> Result<[f64; 2], DepositionError> {
    let slope = sample.temperature_enthalpy_derivative_kg_k_j();
    if sample.liquid_mass_fraction() != 1.0 || !slope.is_finite() || slope <= 0.0 {
        return Err(DepositionError::MaterialDomain);
    }
    let capacity = slope.recip();
    if !capacity.is_finite() || capacity <= 0.0 {
        return Err(DepositionError::Arithmetic);
    }
    let h = sample.specific_enthalpy_j_kg();
    let curve = source.caloric_curve();
    let knots = curve.knots();
    let left = h == source.enthalpy_bounds_j_kg()[1];
    let index = knots
        .partition_point(|knot| if left {
            knot.specific_enthalpy_j_kg < h
        } else {
            knot.specific_enthalpy_j_kg <= h
        })
        .max(1)
        .min(knots.len() - 1) - 1;
    let derivative = match curve.sensible_heat_capacities_on_interval(index) {
        Some([a, b]) => {
            ((b - a) / (knots[index + 1].temperature_k - knots[index].temperature_k)) * slope
        }
        None => 0.0,
    };
    if !derivative.is_finite() {
        return Err(DepositionError::Arithmetic);
    }
    Ok([capacity, derivative])
}

impl SgsThermalLink<'_, '_> {
    fn total_conductivity(
        &self,
        samples: [PhaseConductivitySample; 2],
    ) -> Result<[[f64; 2]; 2], DepositionError> {
        let mut result = [[0.0; 2]; 2];
        for side in 0..2 {
            let molecular = samples[side].conductivity_w_m_k();
            let derivative = samples[side].conductivity_enthalpy_derivative_kg_w_m_k_j();
            result[side] = [molecular, derivative];
            let scale = self.density_diffusivity[side];
            if scale == 0.0 {
                continue;
            }
            let [capacity, capacity_derivative] = sensible_capacity(self.source, samples[side])?;
            let modeled = scale * capacity;
            let total = molecular + modeled;
            let tangent = derivative + scale * capacity_derivative;
            if !modeled.is_finite() || modeled <= 0.0 || !total.is_finite() || !tangent.is_finite() {
                return Err(DepositionError::Arithmetic);
            }
            result[side] = [total, tangent];
        }
        Ok(result)
    }

    fn response(
        &self,
        samples: [PhaseConductivitySample; 2],
        geometry: ConductiveLink,
        area_m2: f64,
    ) -> Result<ImplicitLinkState, DepositionError> {
        if !area_m2.is_finite() || area_m2 <= 0.0
            || geometry.half_distances_m.iter().any(|d| !d.is_finite() || *d <= 0.0)
            || !geometry.contact_resistance_m2_k_w.is_finite()
            || geometry.contact_resistance_m2_k_w < 0.0
        {
            return Err(DepositionError::MaterialDomain);
        }
        let total = self.total_conductivity(samples)?;
        if total.iter().any(|k| k[0] == 0.0) {
            return Ok(ImplicitLinkState {
                conductance_w_k: 0.0,
                power_enthalpy_derivatives_w_kg_j: [0.0; 2],
            });
        }
        let resistance = [
            geometry.half_distances_m[0] / total[0][0],
            geometry.half_distances_m[1] / total[1][0],
        ];
        if resistance.iter().any(|r| !r.is_finite() || *r <= 0.0) {
            return Err(DepositionError::Arithmetic);
        }
        let rc = geometry.contact_resistance_m2_k_w;
        let largest = resistance[0].max(resistance[1]).max(rc);
        let scaled = [resistance[0] / largest, resistance[1] / largest];
        let denominator = scaled[0] + scaled[1] + rc / largest;
        let conductance = (area_m2 / largest) / denominator;
        let delta = samples[0].temperature_k() - samples[1].temperature_k();
        let mut derivatives = [0.0; 2];
        for side in 0..2 {
            let temperature_slope = samples[side].temperature_enthalpy_derivative_kg_k_j();
            let resistance_slope = (scaled[side] / denominator) * (total[side][1] / total[side][0]);
            derivatives[side] = conductance * resistance_slope.mul_add(
                delta,
                if side == 0 { temperature_slope } else { -temperature_slope },
            );
        }
        if !conductance.is_finite() || conductance <= 0.0 || !delta.is_finite()
            || derivatives.iter().any(|value| !value.is_finite())
        {
            return Err(DepositionError::Arithmetic);
        }
        Ok(ImplicitLinkState {
            conductance_w_k: conductance,
            power_enthalpy_derivatives_w_kg_j: derivatives,
        })
    }

    pub(super) fn conductance_from_samples_w_k(
        &self,
        samples: [PhaseConductivitySample; 2],
        geometry: ConductiveLink,
        area_m2: f64,
        gate: &CancelGate,
    ) -> Result<f64, ChtError> {
        poll(gate)?;
        let state = self.response(samples, geometry, area_m2).map_err(|e| native_error(gate, e))?;
        poll(gate)?;
        Ok(state.conductance_w_k)
    }

    #[cfg(test)]
    pub(super) fn conductance_w_k(
        &self,
        h: [f64; 2],
        geometry: ConductiveLink,
        area_m2: f64,
        gate: &CancelGate,
    ) -> Result<f64, ChtError> {
        poll(gate)?;
        let samples = [
            self.source.state_at_specific_enthalpy(h[0]).map_err(|e| refuse(e.to_string()))?,
            self.source.state_at_specific_enthalpy(h[1]).map_err(|e| refuse(e.to_string()))?,
        ];
        self.conductance_from_samples_w_k(samples, geometry, area_m2, gate)
    }
}

impl ImplicitLinkLaw for SgsThermalLink<'_, '_> {
    fn caloric_laws(&self) -> [&dyn ThermalLaw; 2] {
        [self.source.caloric_curve(), self.source.caloric_curve()]
    }

    fn state_at_specific_enthalpies(
        &self,
        h: [f64; 2],
        geometry: ConductiveLink,
        patch: ContactPatch,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitLinkState, DepositionError> {
        checkpoint()?;
        if geometry.link_id != patch.link.id {
            return Err(DepositionError::MaterialDomain);
        }
        let samples = [
            self.source.state_at_specific_enthalpy(h[0]).map_err(|_| DepositionError::MaterialDomain)?,
            self.source.state_at_specific_enthalpy(h[1]).map_err(|_| DepositionError::MaterialDomain)?,
        ];
        let response = self.response(samples, geometry, patch.link.area_m2)?;
        checkpoint()?;
        Ok(response)
    }
}

/// Bind only original fluid--fluid links with a nonzero accepted closure.
/// The same law objects serve native Newton and independent endpoint readback.
#[allow(clippy::too_many_arguments)]
pub(super) fn prepare_links<'a, 'm>(
    domain: &VoxelDomain,
    conductivity: &'a PreparedPhaseConductivity<'m>,
    initial: &DepositionSimulation<'_>,
    mass_kg: &[f64],
    field: SgsThermalField<'_>,
    max_faces: usize,
    gate: &CancelGate,
) -> Result<Vec<(u64, SgsThermalLink<'a, 'm>)>, ChtError> {
    field.admit(domain, mass_kg, gate)?;
    if initial.links().len() > max_faces {
        return Err(refuse("SGS thermal link work exceeds the original conductive-face budget"));
    }
    let mut result = Vec::new();
    for (index, (&id, link)) in initial.links().iter().enumerate() {
        poll_every(gate, index)?;
        let cells = [
            usize::try_from(link.a.0).map_err(|_| refuse("invalid SGS heat cell identity"))?,
            usize::try_from(link.b.0).map_err(|_| refuse("invalid SGS heat cell identity"))?,
        ];
        if cells.iter().any(|&cell| cell >= domain.cell_count()) {
            return Err(refuse("SGS thermal link is outside the original grid"));
        }
        if cells.iter().any(|&cell| !domain.is_fluid(cell)) {
            continue;
        }
        let scales = [
            field.scale(domain, mass_kg, cells[0])?,
            field.scale(domain, mass_kg, cells[1])?,
        ];
        if scales == [0.0; 2] {
            continue;
        }
        result.try_reserve(1).map_err(|_| refuse("SGS thermal link allocation failed"))?;
        result.push((id, SgsThermalLink { source: conductivity, density_diffusivity: scales }));
    }
    poll(gate)?;
    Ok(result)
}

/// Non-cancelling extensive error from the closure actually used by the heat
/// graph versus the accepted momentum closure, at the SAME thermal endpoint.
/// Each link contributes both cell magnitudes, 2*|delta Q|. Double-double edge
/// differences and outward-rounded positive accumulation retain small defects;
/// An unrepresentable nonzero transfer refuses instead of becoming zero.
#[allow(clippy::too_many_arguments)]
pub(super) fn heat_flux_mismatch_j(
    domain: &VoxelDomain,
    conductivity: &PreparedPhaseConductivity<'_>,
    mass_kg: &[f64],
    used: SgsThermalField<'_>,
    actual_eddy: &[f64],
    samples: &[Option<PhaseConductivitySample>],
    duration_s: f64,
    max_faces: usize,
    gate: &CancelGate,
) -> Result<f64, ChtError> {
    used.admit(domain, mass_kg, gate)?;
    let actual = SgsThermalField { eddy_viscosity_m2_s: actual_eddy, ..used };
    actual.admit(domain, mass_kg, gate)?;
    finite_positive("flow.sgs_heat_duration", duration_s)?;
    if samples.len() != domain.cell_count() {
        return Err(refuse("SGS heat feedback needs the original thermal endpoint samples"));
    }
    let mut identical = true;
    for (cell, (&used, &actual)) in used.eddy_viscosity_m2_s.iter().zip(actual_eddy).enumerate() {
        poll_every(gate, cell)?;
        if used != actual {
            identical = false;
            break;
        }
    }
    if identical {
        poll(gate)?;
        return Ok(0.0);
    }
    let mut count = 0_usize;
    let mut total = 0.0;
    for a in 0..domain.cell_count() {
        poll_every(gate, a)?;
        if !domain.is_fluid(a) {
            continue;
        }
        for face in [1, 3, 5] {
            let Some(b) = domain.neighbor(a, face) else { continue };
            if !domain.is_fluid(b) {
                continue;
            }
            count = count.checked_add(1).ok_or_else(|| refuse("SGS heat face count overflow"))?;
            if count > max_faces {
                return Err(refuse("SGS heat feedback exceeds the original conductive-face budget"));
            }
            let pair = [
                samples[a].ok_or_else(|| refuse("missing fluid SGS thermal endpoint"))?,
                samples[b].ok_or_else(|| refuse("missing fluid SGS thermal endpoint"))?,
            ];
            let geometry = ConductiveLink {
                link_id: 0,
                half_distances_m: [
                    0.5 * domain.widths(a)[face / 2],
                    0.5 * domain.widths(b)[face / 2],
                ],
                contact_resistance_m2_k_w: 0.0,
            };
            let mut conductance = [0.0; 2];
            for (slot, field) in [used, actual].into_iter().enumerate() {
                let law = SgsThermalLink {
                    source: conductivity,
                    density_diffusivity: [
                        field.scale(domain, mass_kg, a)?,
                        field.scale(domain, mass_kg, b)?,
                    ],
                };
                conductance[slot] = law.conductance_from_samples_w_k(
                    pair, geometry, domain.face_area(a, face / 2), gate,
                )?;
            }
            let difference = Dd::from_f64(conductance[1]) - Dd::from_f64(conductance[0]);
            let delta = pair[1].temperature_k() - pair[0].temperature_k();
            finite("flow.sgs_heat_temperature_difference", delta)?;
            if difference == Dd::ZERO || delta == 0.0 {
                continue;
            }
            let magnitude = (Dd::from_f64(duration_s) * difference * Dd::from_f64(delta)).abs().to_f64();
            finite_positive("flow.sgs_heat_flux_mismatch", magnitude)?;
            let outward = 2.0 * magnitude.next_up();
            total = (total + outward).next_up();
            finite_positive("flow.sgs_heat_coupling_bound", total)?;
        }
    }
    poll(gate)?;
    Ok(total)
}

#[cfg(test)]
mod tests;
