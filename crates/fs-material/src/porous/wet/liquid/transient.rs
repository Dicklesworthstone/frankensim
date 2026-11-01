//! Conservative one-dimensional evolution of prescribed wetting-liquid storage.
//!
//! This owner closes the gap between the instantaneous van Genuchten--Mualem
//! state and callers that need a changing saturation field. Each finite volume
//! owns porosity, thickness and intrinsic normal permeability. Gas pressure and
//! temperature are prescribed over one step; liquid pressure follows
//! p_l = p_g - Pc(S_l). Darcy face flux uses harmonic hydraulic mobility
//! k*kr/mu and an explicit signed body-force pressure gradient.
//!
//! Explicit and fully implicit backward-Euler candidates share the same source
//! samples, half-cell series resistances and conservative face-flux owner.
//! Boundaries admit sealed faces, signed superficial flux, or external pressure
//! through an explicit surface hydraulic resistance; original material limits
//! apply to the actual material-side face pressure and speed.
//! Neither clips saturation or opens a dry seal with an epsilon. Temperature
//! and gas pressure remain prescribed; no coupled gas or caloric solve is added.

use super::{
    LiquidCapillaryPressure, WetPoreLiquidLaw, WetPoreLiquidState, fraction, invalid, positive,
};
use crate::MaterialError;

mod boundary;
pub use boundary::{LiquidSaturationBoundary, LiquidSaturationBoundaryStepInput};
mod implicit;
pub use implicit::{
    LiquidSaturationImplicitControls, LiquidSaturationImplicitError, LiquidSaturationImplicitReport,
};

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
            || !positive(self.capacity_m())
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

/// Prescribed thermohydraulic coordinates with fixed liquid pressures at both
/// ends. This is the zero-surface-resistance specialization of
/// LiquidSaturationBoundaryStepInput.
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
    pub const fn state(&self) -> &LiquidSaturationColumn {
        &self.state
    }
    pub const fn step(&self) -> &LiquidSaturationStep {
        &self.step
    }
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
        Ok(Self {
            cells,
            saturations,
            states,
            elapsed_s: 0.0,
        })
    }

    pub fn cells(&self) -> &[LiquidSaturationCell] {
        &self.cells
    }
    pub fn saturations(&self) -> &[f64] {
        &self.saturations
    }
    /// Cached constitutive states for allocation-free hydraulic rebinding.
    pub fn states(&self) -> &[WetPoreLiquidState] {
        &self.states
    }
    pub const fn elapsed_s(&self) -> f64 {
        self.elapsed_s
    }

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
        cancelled: impl FnMut() -> bool,
    ) -> Result<LiquidSaturationCandidate, MaterialError> {
        self.explicit_candidate(input.into(), None, cancelled)
    }

    /// Advance one explicit interval with typed boundary conditions and one
    /// fixed original-to-final liquid-volume budget per transverse area [m].
    /// Internal face transfers cancel; positive left and negative right flux
    /// add liquid. The original pressure-only entrypoint retains its existing
    /// report-only volume accounting through this same step owner.
    pub fn step_candidate_with_boundaries(
        &self,
        input: LiquidSaturationBoundaryStepInput,
        volume_balance_tolerance_m: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<LiquidSaturationCandidate, MaterialError> {
        if !positive(volume_balance_tolerance_m) {
            return Err(invalid(
                "explicit liquid step needs a positive finite volume budget",
            ));
        }
        self.explicit_candidate(input, Some(volume_balance_tolerance_m), cancelled)
    }

    fn explicit_candidate(
        &self,
        input: LiquidSaturationBoundaryStepInput,
        volume_balance_tolerance_m: Option<f64>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<LiquidSaturationCandidate, MaterialError> {
        input.validate()?;
        if cancelled() {
            return Err(invalid("transient liquid saturation step cancelled"));
        }
        let elapsed_s = self.elapsed_s + input.duration_s;
        if !elapsed_s.is_finite() || elapsed_s <= self.elapsed_s {
            return Err(invalid("transient liquid time must advance representably"));
        }
        let n = self.cells.len();
        let flow = hydraulics(&self.cells, &self.states, input, || {
            if cancelled() {
                Err(invalid("transient liquid saturation step cancelled"))
            } else {
                Ok(())
            }
        })?;
        let before_pressures = flow.pressures;
        let faces = flow.faces;

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
            elapsed_s,
        };
        let after_pressures = pressures(&next.states, input.gas_pressure_pa)?;
        for i in 0..n {
            let state = &next.states[i];
            let speed = faces[i].abs().max(faces[i + 1].abs());
            if state.liquid_relative_permeability() > 0.0 {
                let p = after_pressures[i]
                    .ok_or_else(|| invalid("accepted mobile liquid pressure is missing"))?;
                admit_state(state, input.temperature_k, p, speed)?;
            } else if speed != 0.0 {
                return Err(invalid("nonzero liquid flux crosses an accepted dry seal"));
            }
        }
        input.admit_faces(&next.cells, &next.states, &after_pressures, &faces)?;
        let before = self.liquid_volume_per_area_m()?;
        let after = next.liquid_volume_per_area_m()?;
        let boundary = input.duration_s * (faces[0] - faces[n]);
        let change = compensated_sum(
            self.cells
                .iter()
                .enumerate()
                .map(|(i, cell)| cell.capacity_m() * (next.saturations[i] - self.saturations[i])),
        )?;
        let residual = change - boundary;
        if !residual.is_finite() {
            return Err(invalid(
                "transient liquid volume balance is unrepresentable",
            ));
        }
        if volume_balance_tolerance_m.is_some_and(|budget| residual.abs() > budget) {
            return Err(invalid(
                "explicit liquid column exceeds its fixed volume-balance budget",
            ));
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
                    Err(invalid(
                        "capillary pressure exceeds prescribed gas pressure",
                    ))
                } else {
                    Ok(Some(p))
                }
            }
            LiquidCapillaryPressure::DryLimit => Ok(None),
        })
        .collect()
}

/// Half-cell resistances in series, including unequal widths. The returned
/// resistance fractions also differentiate this exact mobility in the implicit
/// owner. A zero intrinsic/relative mobility remains an exact seal.
fn interface_mobility(
    a: f64,
    b: f64,
    left_half: f64,
    right_half: f64,
) -> Result<(f64, [f64; 2]), MaterialError> {
    if a == 0.0 || b == 0.0 {
        return Ok((0.0, [0.0; 2]));
    }
    let distance = left_half + right_half;
    if !positive(left_half) || !positive(right_half) || !positive(distance) {
        return Err(invalid("liquid face distance is unrepresentable"));
    }
    let left = left_half / distance;
    let right = right_half / distance;
    let (mobility, shares) = if a <= b {
        let other = right * (a / b);
        let denominator = left + other;
        (a / denominator, [left / denominator, other / denominator])
    } else {
        let other = left * (b / a);
        let denominator = other + right;
        (b / denominator, [other / denominator, right / denominator])
    };
    if !positive(mobility)
        || shares
            .iter()
            .any(|s| !s.is_finite() || *s < 0.0 || *s > 1.0)
    {
        return Err(invalid("series liquid mobility is unrepresentable"));
    }
    Ok((mobility, shares))
}

struct Hydraulics {
    pressures: Vec<Option<f64>>,
    mobility: Vec<f64>,
    faces: Vec<f64>,
    boundaries: [boundary::BoundaryFlow; 2],
}

fn hydraulic_mobility(
    cell: &LiquidSaturationCell,
    state: &WetPoreLiquidState,
    dynamic_viscosity_pa_s: f64,
) -> Result<f64, MaterialError> {
    let value = cell.intrinsic_permeability_m2 * state.liquid_relative_permeability()
        / dynamic_viscosity_pa_s;
    if !value.is_finite()
        || value < 0.0
        || (value == 0.0
            && cell.intrinsic_permeability_m2 > 0.0
            && state.liquid_relative_permeability() > 0.0)
    {
        Err(invalid(
            "transient liquid hydraulic mobility is unrepresentable",
        ))
    } else {
        Ok(value)
    }
}

/// One original source/flux owner for both coefficient time levels. The generic
/// error channel retains the caller's typed cancellation without inspecting text.
fn hydraulics<E: From<MaterialError>>(
    cells: &[LiquidSaturationCell],
    states: &[WetPoreLiquidState],
    input: LiquidSaturationBoundaryStepInput,
    mut checkpoint: impl FnMut() -> Result<(), E>,
) -> Result<Hydraulics, E> {
    checkpoint()?;
    let n = cells.len();
    let pressures = pressures(states, input.gas_pressure_pa)?;
    let mut mobility = Vec::with_capacity(n);
    for (i, (cell, state)) in cells.iter().zip(states).enumerate() {
        if i % 128 == 0 {
            checkpoint()?;
        }
        let support = state.support();
        if input.temperature_k < support.temperature_k.0
            || input.temperature_k > support.temperature_k.1
        {
            return Err(
                invalid("transient liquid temperature lies outside calibration support").into(),
            );
        }
        let value = hydraulic_mobility(cell, state, input.dynamic_viscosity_pa_s)?;
        if let Some(p) = pressures[i] {
            admit_state(state, input.temperature_k, p, 0.0)?;
        } else if value != 0.0 {
            return Err(invalid("mobile liquid requires a finite liquid pressure").into());
        }
        mobility.push(value);
    }
    let mut faces = vec![0.0; n + 1];
    let left_boundary = input.boundaries[0].flow(
        mobility[0],
        pressures[0],
        0.5 * cells[0].thickness_m,
        input.body_force_pressure_gradient_pa_m,
        true,
    )?;
    faces[0] = left_boundary.flux;
    for i in 0..n - 1 {
        if i % 128 == 0 {
            checkpoint()?;
        }
        let left = 0.5 * cells[i].thickness_m;
        let right = 0.5 * cells[i + 1].thickness_m;
        let (m, _) = interface_mobility(mobility[i], mobility[i + 1], left, right)?;
        if m > 0.0 {
            faces[i + 1] = face_flux(
                m,
                pressures[i].expect("admitted mobile pressure"),
                pressures[i + 1].expect("admitted mobile pressure"),
                left + right,
                input.body_force_pressure_gradient_pa_m,
            )?;
        }
    }
    let right_boundary = input.boundaries[1].flow(
        mobility[n - 1],
        pressures[n - 1],
        0.5 * cells[n - 1].thickness_m,
        input.body_force_pressure_gradient_pa_m,
        false,
    )?;
    faces[n] = right_boundary.flux;
    for (i, state) in states.iter().enumerate() {
        if i % 128 == 0 {
            checkpoint()?;
        }
        let speed = faces[i].abs().max(faces[i + 1].abs());
        if let Some(p) = pressures[i] {
            admit_state(state, input.temperature_k, p, speed)?;
        } else if speed != 0.0 {
            return Err(invalid("nonzero liquid flux crosses a dry seal").into());
        }
    }
    input.admit_faces(cells, states, &pressures, &faces)?;
    checkpoint()?;
    Ok(Hydraulics {
        pressures,
        mobility,
        faces,
        boundaries: [left_boundary, right_boundary],
    })
}

fn face_flux(
    mobility_m2_pa_s: f64,
    left_pressure_pa: f64,
    right_pressure_pa: f64,
    distance_m: f64,
    body_force_pressure_gradient_pa_m: f64,
) -> Result<f64, MaterialError> {
    if !positive(distance_m) {
        return Err(invalid("liquid face distance is unrepresentable"));
    }
    let gradient = (right_pressure_pa - left_pressure_pa) / distance_m;
    let q = -mobility_m2_pa_s * (gradient - body_force_pressure_gradient_pa_m);
    if q.is_finite() {
        Ok(q)
    } else {
        Err(invalid(
            "transient liquid Darcy face flux is unrepresentable",
        ))
    }
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
    if sum.is_finite() {
        Ok(sum)
    } else {
        Err(invalid("liquid storage sum is unrepresentable"))
    }
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
        let candidate = dry
            .step_candidate(input([150_000.0, p], 1.0), || false)
            .unwrap();
        assert_eq!(candidate.step().face_flux_m_s[0], 0.0);
        assert_eq!(candidate.step().face_flux_m_s[1], 0.0);
        assert_eq!(dry, saved);

        let wet = column([0.5, 0.5]);
        let p = pressures(&wet.liquid_states().unwrap(), 150_000.0).unwrap()[0].unwrap();
        assert!(
            wet.step_candidate(input([p + 100_000.0, p], 1e9), || false)
                .is_err()
        );
        assert_eq!(wet.saturations(), &[0.5, 0.5]);
    }
}
