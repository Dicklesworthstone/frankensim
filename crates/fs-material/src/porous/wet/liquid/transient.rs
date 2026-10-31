//! Conservative one-dimensional evolution of prescribed wetting-liquid storage.
//!
//! This owner closes the gap between the instantaneous van Genuchten--Mualem
//! state and callers that need a changing saturation field. Each finite volume
//! owns porosity, thickness and intrinsic normal permeability. Gas pressure and
//! temperature are prescribed over one step; liquid pressure follows
//! p_l = p_g - Pc(S_l). Darcy face flux uses harmonic hydraulic mobility
//! k*kr/mu and an explicit signed body-force pressure gradient.
//!
//! The scheme is deliberately first-order explicit finite volume. It preserves
//! liquid volume exactly up to floating-point reduction, never clips saturation,
//! never opens a dry seal with an epsilon, and refuses a step that leaves any
//! constitutive/source domain. Adaptive stepping and a fully implicit Richards
//! solve are separate numerical owners.

use super::{LiquidCapillaryPressure, WetPoreLiquidLaw, WetPoreLiquidState, fraction, invalid, positive};
use crate::MaterialError;

/// One rigid finite volume of porous material, per unit transverse area.
#[derive(Clone, Debug, PartialEq)]
pub struct LiquidSaturationCell {
    pub thickness_m: f64,
    pub porosity: f64,
    pub intrinsic_permeability_m2: f64,
    pub law: WetPoreLiquidLaw,
}
impl LiquidSaturationCell {
    fn validate(&self) -> Result<(), MaterialError> {
        if !positive(self.thickness_m)
            || !self.porosity.is_finite()
            || self.porosity <= 0.0
            || self.porosity > 1.0
            || !self.intrinsic_permeability_m2.is_finite()
            || self.intrinsic_permeability_m2 < 0.0
        {
            return Err(invalid(
                "transient liquid cell needs positive thickness, porosity in (0,1], and nonnegative intrinsic permeability",
            ));
        }
        Ok(())
    }
    fn capacity_m(&self) -> f64 {
        self.porosity * self.thickness_m
    }
}

/// Prescribed thermohydraulic coordinates held over one explicit step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LiquidSaturationStepInput {
    pub duration_s: f64,
    pub temperature_k: f64,
    pub gas_pressure_pa: f64,
    pub dynamic_viscosity_pa_s: f64,
    /// Signed body-force contribution rho*g·x_hat [Pa/m].
    pub body_force_pressure_gradient_pa_m: f64,
    /// Absolute liquid pressures immediately outside the left/right faces.
    pub boundary_liquid_pressure_pa: [f64; 2],
}
impl LiquidSaturationStepInput {
    fn validate(self) -> Result<(), MaterialError> {
        if !positive(self.duration_s)
            || !positive(self.temperature_k)
            || !positive(self.gas_pressure_pa)
            || !positive(self.dynamic_viscosity_pa_s)
            || !self.body_force_pressure_gradient_pa_m.is_finite()
            || self
                .boundary_liquid_pressure_pa
                .iter()
                .any(|&p| !positive(p))
        {
            return Err(invalid(
                "transient liquid step requires positive dt/T/gas pressure/viscosity/boundary liquid pressure and finite body-force gradient",
            ));
        }
        Ok(())
    }
}

/// Accepted conservative step report. Flux is superficial volume flux [m/s],
/// positive left-to-right. Face count is cell_count+1.
#[derive(Clone, Debug, PartialEq)]
pub struct LiquidSaturationStep {
    pub duration_s: f64,
    pub face_flux_m_s: Vec<f64>,
    pub liquid_pressure_pa_before: Vec<Option<f64>>,
    pub liquid_pressure_pa_after: Vec<Option<f64>>,
    pub liquid_volume_per_area_before_m: f64,
    pub liquid_volume_per_area_after_m: f64,
    pub boundary_volume_per_area_m: f64,
    pub balance_residual_m: f64,
}

/// Immutable candidate so a rejected solve never mutates accepted saturation.
#[derive(Clone, Debug, PartialEq)]
pub struct LiquidSaturationCandidate {
    state: LiquidSaturationColumn,
    step: LiquidSaturationStep,
}
impl LiquidSaturationCandidate {
    pub const fn state(&self) -> &LiquidSaturationColumn { &self.state }
    pub const fn step(&self) -> &LiquidSaturationStep { &self.step }
    pub fn into_parts(self) -> (LiquidSaturationColumn, LiquidSaturationStep) {
        (self.state, self.step)
    }
}

/// Rigid one-dimensional porous column with conservative liquid storage.
#[derive(Clone, Debug, PartialEq)]
pub struct LiquidSaturationColumn {
    cells: Vec<LiquidSaturationCell>,
    saturations: Vec<f64>,
    states: Vec<WetPoreLiquidState>,
    elapsed_s: f64,
}
impl LiquidSaturationColumn {
    pub fn new(
        cells: Vec<LiquidSaturationCell>,
        saturations: Vec<f64>,
    ) -> Result<Self, MaterialError> {
        if cells.is_empty() || cells.len() != saturations.len() || cells.len() > 4096 {
            return Err(invalid(
                "transient liquid column needs 1..4096 cells with one saturation per cell",
            ));
        }
        let mut states = Vec::with_capacity(cells.len());
        for (cell, &s) in cells.iter().zip(&saturations) {
            cell.validate()?;
            if !fraction(s) {
                return Err(invalid(
                    "initial liquid saturation must be finite and inside each cell calibration",
                ));
            }
            states.push(cell.law.sample(s)?);
        }
        Ok(Self { cells, saturations, states, elapsed_s: 0.0 })
    }

    pub fn cells(&self) -> &[LiquidSaturationCell] { &self.cells }
    pub fn saturations(&self) -> &[f64] { &self.saturations }
    /// Cached constitutive states for allocation-free hydraulic rebinding.
    pub fn states(&self) -> &[WetPoreLiquidState] { &self.states }
    pub const fn elapsed_s(&self) -> f64 { self.elapsed_s }

    /// Compatibility copy of the cached material states.
    pub fn liquid_states(&self) -> Result<Vec<WetPoreLiquidState>, MaterialError> {
        Ok(self.states.clone())
    }

    pub fn liquid_volume_per_area_m(&self) -> Result<f64, MaterialError> {
        compensated_sum(
            self.cells
                .iter()
                .zip(&self.saturations)
                .map(|(cell, &s)| cell.capacity_m() * s),
        )
    }

    /// Advance one explicit conservative interval.
    ///
    /// A face whose harmonic mobility is exactly zero has exactly zero flux;
    /// its undefined/infinite dry-limit capillary pressure is never regularized.
    /// Every nonzero-mobility cell must have a finite liquid pressure and all
    /// actual T/p/speed coordinates must fit that cell's calibration.
    pub fn step_candidate(
        &self,
        input: LiquidSaturationStepInput,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<LiquidSaturationCandidate, MaterialError> {
        input.validate()?;
        if cancelled() {
            return Err(invalid("transient liquid saturation step cancelled"));
        }
        let n = self.cells.len();
        let before_states = &self.states;
        let before_pressures = pressures(before_states, input.gas_pressure_pa)?;
        let mut mobility = Vec::with_capacity(n);
        for (i, ((cell, state), pressure)) in self
            .cells
            .iter()
            .zip(before_states)
            .zip(&before_pressures)
            .enumerate()
        {
            if i % 128 == 0 && cancelled() {
                return Err(invalid("transient liquid saturation step cancelled"));
            }
            let kr = state.liquid_relative_permeability();
            let value = (cell.intrinsic_permeability_m2 * kr) / input.dynamic_viscosity_pa_s;
            if !value.is_finite() || value < 0.0 {
                return Err(invalid("transient liquid hydraulic mobility is unrepresentable"));
            }
            if value > 0.0 {
                let p = pressure.ok_or_else(|| invalid(
                    "mobile liquid state requires finite capillary/liquid pressure",
                ))?;
                admit_state(state, input.temperature_k, p, 0.0)?;
            }
            mobility.push(value);
        }

        let mut faces = vec![0.0; n + 1];
        if mobility[0] > 0.0 {
            let p0 = before_pressures[0].ok_or_else(|| invalid("missing mobile left-cell liquid pressure"))?;
            faces[0] = face_flux(
                mobility[0],
                input.boundary_liquid_pressure_pa[0],
                p0,
                0.5 * self.cells[0].thickness_m,
                input.body_force_pressure_gradient_pa_m,
            )?;
        }
        for i in 0..n - 1 {
            let m = harmonic(mobility[i], mobility[i + 1])?;
            if m > 0.0 {
                let left = before_pressures[i].ok_or_else(|| invalid("missing mobile left liquid pressure"))?;
                let right = before_pressures[i + 1].ok_or_else(|| invalid("missing mobile right liquid pressure"))?;
                faces[i + 1] = face_flux(
                    m,
                    left,
                    right,
                    0.5 * (self.cells[i].thickness_m + self.cells[i + 1].thickness_m),
                    input.body_force_pressure_gradient_pa_m,
                )?;
            }
        }
        if mobility[n - 1] > 0.0 {
            let pn = before_pressures[n - 1].ok_or_else(|| invalid("missing mobile right-cell liquid pressure"))?;
            faces[n] = face_flux(
                mobility[n - 1],
                pn,
                input.boundary_liquid_pressure_pa[1],
                0.5 * self.cells[n - 1].thickness_m,
                input.body_force_pressure_gradient_pa_m,
            )?;
        }

        let mut next_saturation = Vec::with_capacity(n);
        for i in 0..n {
            let delta = input.duration_s * (faces[i] - faces[i + 1]) / self.cells[i].capacity_m();
            let next = self.saturations[i] + delta;
            if !next.is_finite() || self.cells[i].law.sample(next).is_err() {
                return Err(invalid(
                    "explicit liquid saturation step leaves constitutive support; reduce timestep or use an implicit owner",
                ));
            }
            next_saturation.push(next);
        }
        let after_states = self
            .cells
            .iter()
            .zip(&next_saturation)
            .map(|(cell, &s)| cell.law.sample(s))
            .collect::<Result<Vec<_>, _>>()?;
        let next = Self {
            cells: self.cells.clone(),
            saturations: next_saturation,
            states: after_states,
            elapsed_s: self.elapsed_s + input.duration_s,
        };
        if !next.elapsed_s.is_finite() {
            return Err(invalid("transient liquid time is unrepresentable"));
        }
        let after_pressures = pressures(&next.states, input.gas_pressure_pa)?;
        for i in 0..n {
            let state = &next.states[i];
            let speed = faces[i].abs().max(faces[i + 1].abs());
            if state.liquid_relative_permeability() > 0.0 {
                let p = after_pressures[i].ok_or_else(|| invalid("accepted mobile liquid pressure is missing"))?;
                admit_state(state, input.temperature_k, p, speed)?;
            } else if speed != 0.0 {
                return Err(invalid("nonzero liquid flux crosses an accepted dry seal"));
            }
        }
        let before = self.liquid_volume_per_area_m()?;
        let after = next.liquid_volume_per_area_m()?;
        let boundary = input.duration_s * (faces[0] - faces[n]);
        let residual = (after - before) - boundary;
        if !residual.is_finite() {
            return Err(invalid("transient liquid volume balance is unrepresentable"));
        }
        if cancelled() {
            return Err(invalid("transient liquid saturation step cancelled"));
        }
        Ok(LiquidSaturationCandidate {
            state: next,
            step: LiquidSaturationStep {
                duration_s: input.duration_s,
                face_flux_m_s: faces,
                liquid_pressure_pa_before: before_pressures,
                liquid_pressure_pa_after: after_pressures,
                liquid_volume_per_area_before_m: before,
                liquid_volume_per_area_after_m: after,
                boundary_volume_per_area_m: boundary,
                balance_residual_m: residual,
            },
        })
    }
}

fn admit_state(
    state: &WetPoreLiquidState,
    temperature_k: f64,
    liquid_pressure_pa: f64,
    speed_m_s: f64,
) -> Result<(), MaterialError> {
    let support = state.support();
    if temperature_k < support.temperature_k.0
        || temperature_k > support.temperature_k.1
        || liquid_pressure_pa < support.liquid_pressure_pa.0
        || liquid_pressure_pa > support.liquid_pressure_pa.1
        || !speed_m_s.is_finite()
        || speed_m_s > support.maximum_superficial_speed_m_s
    {
        return Err(invalid(
            "transient liquid state lies outside calibration temperature/liquid-pressure/speed support",
        ));
    }
    Ok(())
}

fn pressures(
    states: &[WetPoreLiquidState],
    gas_pressure_pa: f64,
) -> Result<Vec<Option<f64>>, MaterialError> {
    states
        .iter()
        .map(|state| match state.capillary_pressure() {
            LiquidCapillaryPressure::Finite(pc) => {
                let p = gas_pressure_pa - pc;
                if !positive(p) {
                    Err(invalid("capillary pressure exceeds prescribed gas pressure"))
                } else {
                    Ok(Some(p))
                }
            }
            LiquidCapillaryPressure::DryLimit => Ok(None),
        })
        .collect()
}

fn harmonic(a: f64, b: f64) -> Result<f64, MaterialError> {
    if a == 0.0 || b == 0.0 {
        return Ok(0.0);
    }
    let h = 2.0 / (1.0 / a + 1.0 / b);
    if positive(h) { Ok(h) } else { Err(invalid("harmonic liquid mobility is unrepresentable")) }
}

fn face_flux(
    mobility_m2_pa_s: f64,
    left_pressure_pa: f64,
    right_pressure_pa: f64,
    distance_m: f64,
    body_force_pressure_gradient_pa_m: f64,
) -> Result<f64, MaterialError> {
    let gradient = (right_pressure_pa - left_pressure_pa) / distance_m;
    let q = -mobility_m2_pa_s * (gradient - body_force_pressure_gradient_pa_m);
    if q.is_finite() { Ok(q) } else { Err(invalid("transient liquid Darcy face flux is unrepresentable")) }
}

fn compensated_sum(values: impl Iterator<Item = f64>) -> Result<f64, MaterialError> {
    let mut sum = 0.0;
    let mut c = 0.0;
    for value in values {
        let y = value - c;
        let t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    if sum.is_finite() { Ok(sum) } else { Err(invalid("liquid storage sum is unrepresentable")) }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn law() -> WetPoreLiquidLaw {
        WetPoreLiquidLaw::try_new(
            0.1,
            0.1,
            0.5,
            0.5,
            10_000.0,
            super::super::WetPoreLiquidSupport {
                liquid_saturation: (0.1, 0.9),
                temperature_k: (280.0, 330.0),
                liquid_pressure_pa: (50_000.0, 200_000.0),
                maximum_superficial_speed_m_s: 0.1,
            },
            "synthetic transient retention fixture",
        )
        .unwrap()
    }
    fn column(s: [f64; 2]) -> LiquidSaturationColumn {
        LiquidSaturationColumn::new(
            vec![
                LiquidSaturationCell {
                    thickness_m: 0.01,
                    porosity: 0.4,
                    intrinsic_permeability_m2: 2e-10,
                    law: law(),
                },
                LiquidSaturationCell {
                    thickness_m: 0.01,
                    porosity: 0.4,
                    intrinsic_permeability_m2: 2e-10,
                    law: law(),
                },
            ],
            s.to_vec(),
        )
        .unwrap()
    }
    fn input(boundary: [f64; 2], dt: f64) -> LiquidSaturationStepInput {
        LiquidSaturationStepInput {
            duration_s: dt,
            temperature_k: 300.0,
            gas_pressure_pa: 150_000.0,
            dynamic_viscosity_pa_s: 1e-3,
            body_force_pressure_gradient_pa_m: 0.0,
            boundary_liquid_pressure_pa: boundary,
        }
    }

    #[test]
    fn uniform_equilibrium_is_exactly_stationary() {
        let state = column([0.5, 0.5]);
        let p = pressures(&state.liquid_states().unwrap(), 150_000.0).unwrap()[0].unwrap();
        let candidate = state.step_candidate(input([p, p], 1.0), || false).unwrap();
        assert_eq!(candidate.state().saturations(), state.saturations());
        assert_eq!(candidate.step().face_flux_m_s, vec![0.0; 3]);
        assert_eq!(candidate.step().balance_residual_m, 0.0);
    }

    #[test]
    fn boundary_infiltration_conserves_liquid_volume() {
        let state = column([0.5, 0.5]);
        let p = pressures(&state.liquid_states().unwrap(), 150_000.0).unwrap()[0].unwrap();
        let candidate = state
            .step_candidate(input([p + 500.0, p], 0.1), || false)
            .unwrap();
        assert!(candidate.state().saturations()[0] > 0.5);
        assert!(candidate.step().face_flux_m_s[0] > 0.0);
        assert!(
            candidate.step().balance_residual_m.abs()
                <= 32.0 * f64::EPSILON * candidate.step().liquid_volume_per_area_after_m
        );
    }

    #[test]
    fn dry_seal_and_oversized_explicit_step_refuse_without_mutation() {
        let dry = column([0.1, 0.5]);
        let saved = dry.clone();
        let p = pressures(&dry.liquid_states().unwrap(), 150_000.0).unwrap()[1].unwrap();
        let candidate = dry.step_candidate(input([150_000.0, p], 1.0), || false).unwrap();
        assert_eq!(candidate.step().face_flux_m_s[0], 0.0);
        assert_eq!(candidate.step().face_flux_m_s[1], 0.0);
        assert_eq!(dry, saved);

        let wet = column([0.5, 0.5]);
        let p = pressures(&wet.liquid_states().unwrap(), 150_000.0).unwrap()[0].unwrap();
        assert!(wet.step_candidate(input([p + 100_000.0, p], 1e9), || false).is_err());
        assert_eq!(wet.saturations(), &[0.5, 0.5]);
    }
}
