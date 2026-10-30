//! Molecular and thermal diffusion memory in ONE thermodynamic pressure law.
//!
//! Let cp(s)=cp_inf+sum Cj/(1+s*tau_j), using the original material owner.
//! Each original circular thermal quadrature mode has weight w_n and
//! d_n=cp_eq*lambda_n. Eliminating temperature gives the per-volume compliance
//! kappa(s)=1/(rho*R*T) * [1-R*sum w_n*s/(s*cp(s)+d_n)].
//! With no wall, one w=1,d=0 mode recovers the material's exact relaxing
//! wavenumber. With no reservoirs this is the existing finite Bessel law.
//! Adding independent molecular and equilibrium thermal shunts would instead
//! double-count compressibility and give the wrong isothermal limit.
//!
//! Pressure uses frozen compliance. For each radial mode let
//! t_n=y_n+sqrt(w_n*R/cv_inf)*x_p. Its physical loss is
//! (d_n/cp_inf)*t_n^2 + sum (z_nj-sqrt(Cj/cp_inf)*t_n)^2/tau_j.
//! Stamp the NEGATIVE outer products of those vectors. Thus the original DG
//! exponential advances pressure, flow, thermal and molecular states together,
//! retaining positive storage and nonnegative independent loss squares.
//! This follows the linear ideal-gas energy equation with relaxing internal
//! temperatures; it does not add a second bulk-viscosity or SK absorption term.
//! Axial conduction, finite wall heat capacity, mean flow and evolving species
//! are absent. The existing Bessel band admission concerns its equilibrium
//! spectral measure; it is not a new error certificate for the coupled law.
//!
//! A spatial wall uses the SAME positive axial quadrature functional in place
//! of x_p. Its radius-dependent Bessel modes therefore couple to molecular
//! reservoirs at their original physical point. Uniform gas makes the frozen
//! pressure mass map a scalar rescaling, so the normalized functional is
//! unchanged: sqrt(w*C_inf) M_inf^(-1/2) = sqrt(w*C_eq) M_eq^(-1/2).
//! The axial projectors sum to identity only to the selected quadrature's
//! area-weighted mass accuracy. Isothermal/equilibrium storage on a curved
//! bore therefore needs axial refinement; the original mass is never lumped.
//! This is a convergent local-cylinder Webster wall discretization on tapered
//! geometry. It does not introduce spherical-cone or three-dimensional physics.
//! The pointwise field attachment additionally resolves the same law on the
//! original primitive field and integrates its frozen mass. A local matrix
//! congruence replaces the scalar rescaling; see [`DgDuct::with_molecular_relaxation_field`].

use super::viscothermal::ViscothermalMemory;
use super::*;
use fs_material::relaxation::{RelaxationAbsorption, RelaxationState};

mod field;
mod update;
pub(super) mod heat;

#[derive(Clone, Debug)]
pub(super) struct MolecularMemory {
    state: Option<RelaxationState>,
    field_law: Option<RelaxationAbsorption>,
    bulk: f64,
    modes: Vec<ThermalMode>,
}

#[derive(Clone, Debug)]
struct ThermalMode {
    start: usize,
    pressure_scale: f64,
    heat_rate: f64,
    projection: PressureProjection,
    reservoir_scales: Vec<f64>,
    rates: Vec<f64>,
}

#[derive(Clone, Debug)]
enum PressureProjection {
    /// One identical local material law for each exact-mass coordinate.
    Uniform,
    /// One original positive spatial quadrature point in an element.
    Quadrature { first: usize, weights: Vec<f64> },
}

impl PressureProjection {
    fn width(&self, field_dofs: usize) -> usize {
        match self {
            Self::Uniform => field_dofs,
            Self::Quadrature { .. } => 1,
        }
    }

    fn value(&self, field: &[f64], i: usize) -> f64 {
        match self {
            Self::Uniform => field[i],
            Self::Quadrature { first, weights } => weights
                .iter()
                .enumerate()
                .map(|(j, weight)| weight * field[*first + j])
                .sum(),
        }
    }

    fn entries(&self, i: usize, scale: f64) -> Vec<(usize, f64)> {
        match self {
            Self::Uniform => vec![(i, scale)],
            Self::Quadrature { first, weights } => weights
                .iter()
                .enumerate()
                .map(|(j, weight)| (*first + j, scale * weight))
                .collect(),
        }
    }
}

fn invalid(what: &'static str) -> DuctError {
    DuctError::BadParameter { what }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), DuctError> {
    if cancelled() {
        Err(DuctError::Cancelled { completed_steps: 0 })
    } else {
        Ok(())
    }
}

fn stamp_square(operator: &mut [f64], n: usize, terms: &[(usize, f64)], rate: f64) {
    for &(i, a) in terms {
        for &(j, b) in terms {
            operator[i * n + j] -= rate * a * b;
        }
    }
}

fn restore_pressure_boundaries(duct: &DgDuct, operator: &mut [f64], n: usize) {
    let m = duct.field_dofs;
    for i in 0..m {
        operator[i * n..i * n + m].fill(0.0);
    }
    let local = duct.mesh.order + 1;
    let last = duct.elements.len() - 1;
    for (boundary, element, node, gas) in [
        (duct.inlet, 0, 0, &duct.gas),
        (duct.outlet, last, local - 1, &duct.outlet_gas),
    ] {
        if boundary != DgBoundary::Matched {
            continue;
        }
        let e = &duct.elements[element];
        let resistance = gas.density * gas.sound_speed / e.areas[node];
        let factor = det::sqrt(0.5 / resistance);
        let entries: Vec<_> = (0..local)
            .map(|i| {
                (
                    element * local + i,
                    factor * e.mp_inv_sqrt[i * local + node],
                )
            })
            .collect();
        stamp_square(operator, n, &entries, 1.0);
    }
}

fn stamp_mode(operator: &mut [f64], n: usize, mode: &ThermalMode, field_dofs: usize) {
    let width = mode.projection.width(field_dofs);
    for i in 0..width {
        let mut temperature = mode.projection.entries(i, mode.pressure_scale);
        temperature.push((mode.start + i, 1.0));
        if mode.heat_rate != 0.0 {
            stamp_square(operator, n, &temperature, mode.heat_rate);
        }
        for (j, (&b, &rate)) in mode.reservoir_scales.iter().zip(&mode.rates).enumerate() {
            let mut slip: Vec<_> = temperature
                .iter()
                .map(|&(index, value)| (index, -b * value))
                .collect();
            slip.push((mode.start + (j + 1) * width + i, 1.0));
            stamp_square(operator, n, &slip, rate);
        }
    }
}

impl DgDuct {
    /// Couple the exact evaluated molecular reservoirs to the original pressure
    /// and the original circular thermal diffusion histories.
    ///
    /// Attach AFTER viscothermal memory and BEFORE compliant walls. Uniform
    /// exact-mass and spatial wall quadrature are both supported: tapered and
    /// profiled bores retain their original local radius and wall measure.
    /// Bare variable-area ducts support the adiabatic molecular law. A
    /// nonuniform gas field refuses instead of using a mean material state.
    /// The original equilibrium GasState and source field remain unchanged;
    /// only the internal pressure energy map uses frozen compressibility.
    /// Matched boundaries retain their original equilibrium resistance, which
    /// is not a broadband nonreflecting termination in a dispersive gas.
    ///
    /// This is cold construction from an evaluated law. Live uniform changes
    /// use `rebind_uniform_molecular_gas_field` after source-backed field
    /// attachment; this evaluated-only route retains no original declaration
    /// to re-evaluate and cannot refresh. Empty reservoirs
    /// preserve the existing operator and all original coefficients bitwise.
    pub fn with_molecular_relaxation(self, state: RelaxationState) -> Result<Self, DuctError> {
        self.with_molecular_relaxation_controlled(state, &mut || false)
    }

    /// Cancellable counterpart; no candidate is published on any refusal.
    pub fn with_molecular_relaxation_controlled(
        mut self,
        state: RelaxationState,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, DuctError> {
        poll(cancelled)?;
        if self.molecular.is_some() || self.compliant_wall.is_some() {
            return Err(invalid(
                "attach molecular gas once, before compliant wall mechanics",
            ));
        }
        if state.gas() != &self.gas
            || self
                .gas_field
                .as_ref()
                .is_some_and(|field| field.uniform_state().is_none())
        {
            return Err(invalid(
                "molecular gas must exactly match the complete uniform DG material",
            ));
        }
        let equilibrium_bulk = self.gas.density * self.gas.sound_speed * self.gas.sound_speed;
        if state
            .processes()
            .iter()
            .all(|p| p.heat_capacity_j_kg_k == 0.0)
        {
            self.molecular = Some(MolecularMemory {
                state: Some(state),
                field_law: None,
                bulk: equilibrium_bulk,
                modes: Vec::new(),
            });
            return Ok(self);
        }
        let cv = state.frozen_cv_j_kg_k();
        let cp = cv + self.gas.specific_gas_constant;
        let speed = state.frozen_sound_speed_m_s();
        let bulk = self.gas.density * speed * speed;
        let scale = speed / self.gas.sound_speed;
        let active: Vec<_> = state
            .processes()
            .iter()
            .filter(|p| p.heat_capacity_j_kg_k > 0.0)
            .collect();
        let molecular = active.len();
        let reservoir_scales: Vec<_> = active
            .iter()
            .map(|p| det::sqrt(p.heat_capacity_j_kg_k / cp))
            .collect();
        let rates: Vec<_> = active.iter().map(|p| 1.0 / p.relaxation_time_s).collect();
        if [cv, cp, bulk, scale]
            .iter()
            .chain(&rates)
            .chain(&reservoir_scales)
            .any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(invalid(
                "molecular thermodynamic coordinates are unrepresentable",
            ));
        }
        let m = self.field_dofs;
        let old_total = self.total_dofs();
        let gas = self.gas;
        let mut modes = Vec::new();
        let mut append_modes =
            |terms: &[(f64, f64)], projection: PressureProjection| -> Result<(), DuctError> {
                let weight_sum: f64 = terms.iter().map(|r| r.0).sum();
                if terms.is_empty()
                    || (weight_sum - 1.0).abs() > 4e-12
                    || terms
                        .iter()
                        .any(|(w, r)| !w.is_finite() || *w <= 0.0 || !r.is_finite() || *r < 0.0)
                {
                    return Err(invalid(
                        "original thermal measure lost its positive storage moment",
                    ));
                }
                for &(weight, heat_rate) in terms {
                    modes.push(ThermalMode {
                        start: 0,
                        pressure_scale: det::sqrt(weight * gas.specific_gas_constant / cv),
                        heat_rate,
                        projection: projection.clone(),
                        reservoir_scales: reservoir_scales.clone(),
                        rates: rates.clone(),
                    });
                }
                Ok(())
            };
        let radial_terms = |terms: &[fs_phs::ZkDiffusiveTerm]| -> Vec<(f64, f64)> {
            terms
                .iter()
                .map(|t| {
                    (
                        t.gain_per_s / (t.rate_per_s * (gas.gamma - 1.0)),
                        t.rate_per_s * gas.specific_heat_cp / cp,
                    )
                })
                .collect()
        };
        let kept: Vec<usize> = match &mut self.viscothermal {
            Some(ViscothermalMemory::Uniform(model)) => {
                append_modes(
                    &radial_terms(model.shunt_terms()),
                    PressureProjection::Uniform,
                )?;
                (0..2 * m + model.series_terms().len() * m).collect()
            }
            None => {
                append_modes(&[(1.0, 0.0)], PressureProjection::Uniform)?;
                (0..2 * m).collect()
            }
            Some(ViscothermalMemory::Spatial(wall)) => {
                let (kept, channels) = wall.split_for_molecular(m);
                for channel in channels {
                    poll(cancelled)?;
                    append_modes(
                        &radial_terms(&channel.terms),
                        PressureProjection::Quadrature {
                            first: channel.field_start,
                            weights: channel.projection,
                        },
                    )?;
                }
                kept
            }
        };
        let prefix = kept.len();
        let mut total = prefix;
        for mode in &mut modes {
            mode.start = total;
            total = mode
                .projection
                .width(m)
                .checked_mul(molecular + 1)
                .and_then(|extra| total.checked_add(extra))
                .filter(|v| *v <= MAX_DG_PROPAGATOR_DOFS)
                .ok_or_else(|| invalid("molecular plus thermal gas exceeds the DG state budget"))?;
        }
        let mut operator = vec![0.0; total * total];
        for i in 0..prefix {
            poll(cancelled)?;
            for j in 0..prefix {
                operator[i * total + j] = self.l_tilde[kept[i] * old_total + kept[j]];
            }
        }
        // The original pressure-pressure block is only matched-boundary loss
        // plus the thermal arms being replaced. Re-form its original boundary
        // squares; subtracting large thermal rates would leave roundoff of
        // either sign in a physically lossless pressure direction.
        restore_pressure_boundaries(&self, &mut operator, total);
        if operator.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "original pressure boundary loss is unrepresentable",
            ));
        }
        for i in 0..prefix {
            poll(cancelled)?;
            for j in 0..prefix {
                if i < m {
                    operator[i * total + j] *= scale;
                }
                if j < m {
                    operator[i * total + j] *= scale;
                }
            }
        }
        for mode in &modes {
            poll(cancelled)?;
            stamp_mode(&mut operator, total, mode, m);
        }
        if operator.iter().any(|v| !v.is_finite()) {
            return Err(invalid("coupled molecular/thermal operator is not finite"));
        }
        for element in &mut self.elements {
            for value in &mut element.mp_sqrt {
                *value /= scale;
            }
            for value in &mut element.mp_inv_sqrt {
                *value *= scale;
            }
        }
        for (_, port) in &mut self.ports {
            port.resize(total, 0.0);
            for value in &mut port[..m] {
                *value *= scale;
            }
        }
        self.l_tilde = operator;
        self.memory_dofs = total - 2 * m;
        self.g_block = None;
        self.molecular = Some(MolecularMemory {
            state: Some(state),
            field_law: None,
            bulk,
            modes,
        });
        self.operator_identity = hash_operator(
            self.mesh,
            self.inlet,
            self.outlet,
            &self.elements,
            &self.l_tilde,
            &self.ports,
        );
        poll(cancelled)?;
        Ok(self)
    }

}

impl<M: DgMedium> DgDuct<M> {
    /// Original material evaluation; no fitted relaxation parameters.
    #[must_use]
    pub fn molecular_relaxation(&self) -> Option<&RelaxationState> {
        self.molecular.as_ref().and_then(|m| m.state.as_ref())
    }

    pub(super) fn pressure_bulk_modulus(&self) -> f64 {
        self.molecular.as_ref().map_or(
            self.gas.density_kg_m3() * self.gas.sound_speed_m_s() * self.gas.sound_speed_m_s(),
            |m| m.bulk,
        )
    }
    pub(super) fn molecular_replaces_thermal(&self) -> bool {
        self.molecular.as_ref().is_some_and(|m| !m.modes.is_empty())
    }
    pub(super) fn molecular_thermal_loss(&self, state: &DgState) -> f64 {
        let Some(memory) = &self.molecular else {
            return 0.0;
        };
        let mut power = 0.0;
        for mode in &memory.modes {
            for i in 0..mode.projection.width(self.field_dofs) {
                let temperature = state.tilde[mode.start + i]
                    + mode.pressure_scale * mode.projection.value(&state.tilde, i);
                power += mode.heat_rate * temperature * temperature;
            }
        }
        power
    }

    /// Instantaneous nonnegative internal molecular loss [W], separate from
    /// wall heat conduction, wall shear, solid viscosity and outlet absorption.
    pub fn molecular_dissipated_power_w(&self, state: &DgState) -> Result<f64, DuctError> {
        if state.operator_identity != self.operator_identity
            || state.tilde.len() != self.total_dofs()
            || state.tilde.iter().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "molecular loss needs a finite state of its exact DG operator",
            ));
        }
        let Some(memory) = &self.molecular else {
            return Ok(0.0);
        };
        let mut power = 0.0;
        let m = self.field_dofs;
        for mode in &memory.modes {
            let width = mode.projection.width(m);
            for i in 0..width {
                let temperature = state.tilde[mode.start + i]
                    + mode.pressure_scale * mode.projection.value(&state.tilde, i);
                for (j, (&b, &rate)) in mode.reservoir_scales.iter().zip(&mode.rates).enumerate() {
                    let slip = state.tilde[mode.start + (j + 1) * width + i] - b * temperature;
                    power += rate * slip * slip;
                }
            }
        }
        if !power.is_finite() {
            return Err(invalid("molecular loss is unrepresentable"));
        }
        Ok(power)
    }
}

#[cfg(test)]
mod tests;
