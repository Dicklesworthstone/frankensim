//! Continuous thermal impedance on the original enthalpy graph.
//!
//! Peak phasors use exp(-i omega t). With extensive stored energy E=m*h,
//! K=-d(net inward watts)/dE, solve (K-i*omega*I) delta E=delta P and recover
//! delta T=(dT/dh)/m * delta E. Thus T/P is in K/W; neither a fitted heat
//! capacity nor a time-discretized thermal RC surrogate replaces the source.
//! A temperature/heat-rate product has units K W, not physical power.
//!
//! Preparation requires steady cell heat balance. External clocks, baths and
//! prescribed geometry are frozen; native density-geometry tangents remain.
//! Callers admit differentiable branches and small perturbations within them.
//! Latent plateaus retain energy response with zero dT/dh. One-sided corners
//! do not establish two-sided harmonic response. No history advances, and no
//! stability, experimental or continuum-error certificate is implied.

use super::*;
use fs_la::eigen_complex::{LuComplex, LuComplexError, lu_complex_with_cancel, with_lu_workers};
use fs_math::c64::C64;

#[path = "harmonic/conservation.rs"]
mod conservation;
/// Closed component response with conserved original stored energy.
#[path = "harmonic/isolated.rs"]
pub mod isolated;

/// Added dense resources, excluding native graph maps and allocator overhead.
/// Native preparation keeps its graph budgets; at most 1024 cells are admitted.
#[derive(Clone, Copy, Debug)]
pub struct ThermalImpedanceBudget {
    pub maximum_cells: usize,
    /// Conservative added numerical payload: 64*n*n + 256*n bytes.
    pub maximum_dense_bytes: usize,
    /// Bound n*n*(n+16) for one preparation, factor or drive.
    pub maximum_dense_work: u64,
    /// Maximum absolute original net inward heat per cell at the operating point.
    pub equilibrium_tolerance_w: f64,
    /// Componentwise backward residual in the original continuous equations.
    pub relative_residual_tolerance: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum ThermalImpedanceError {
    Implicit(ImplicitConductionError),
    Linear(LuComplexError),
    NotEquilibrium { cell: CellId, net_power_w: f64 },
    Residual { relative: f64, tolerance: f64 },
}
impl From<ImplicitConductionError> for ThermalImpedanceError {
    fn from(error: ImplicitConductionError) -> Self {
        Self::Implicit(error)
    }
}
impl From<DepositionError> for ThermalImpedanceError {
    fn from(error: DepositionError) -> Self {
        Self::Implicit(error.into())
    }
}
impl From<LuComplexError> for ThermalImpedanceError {
    fn from(error: LuComplexError) -> Self {
        Self::Linear(error)
    }
}
impl core::fmt::Display for ThermalImpedanceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "thermal impedance: {self:?}")
    }
}
impl std::error::Error for ThermalImpedanceError {}

type Result<T> = core::result::Result<T, ThermalImpedanceError>;
fn invalid(message: &'static str) -> ThermalImpedanceError {
    ImplicitConductionError::Invalid(message).into()
}
fn complex(value: C64) -> Result<C64> {
    finite(value.re)?;
    finite(value.im)?;
    Ok(value)
}

/// Immutable original-equilibrium tangent, reusable at many frequencies.
pub struct ThermalImpedanceModel {
    ids: Vec<CellId>,
    time_s: f64,
    temperatures_k: Vec<f64>,
    temperature_per_energy_k_j: Vec<f64>,
    outward_tangent_per_s: Vec<f64>,
    budget: ThermalImpedanceBudget,
}

/// One row-equilibrated factor, reused for all driving and mutual impedances.
pub struct ThermalImpedanceFactor<'a> {
    model: &'a ThermalImpedanceModel,
    omega_rad_s: f64,
    row_scales: Vec<f64>,
    lu: LuComplex,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ThermalImpedanceResponse {
    /// Extensive original stored-energy perturbations, in cell_ids() order.
    pub energy_j: Vec<C64>,
    /// Temperature perturbations, including exactly zero on a latent plateau.
    pub temperature_k: Vec<C64>,
    pub maximum_residual_w: f64,
    pub relative_residual: f64,
}

impl ImplicitConductionChain<'_> {
    /// Linearize native conduction, contacts, boundaries, source feedback and
    /// prescribed whole-enthalpy upwind transport at this frozen operating clock.
    /// The directed advective tangent remains nonsymmetric. Inlet data/rates are
    /// held constant; heat-drive phasors do not represent mass-flow modulation.
    pub fn prepare_thermal_impedance(
        &self,
        thermal: &DepositionSimulation<'_>,
        budget: ThermalImpedanceBudget,
        mut checkpoint: impl FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<ThermalImpedanceModel> {
        checkpoint()?;
        let n = self.ids.len();
        let square = n as u128 * n as u128;
        if n == 0
            || n > 1024
            || n > budget.maximum_cells
            || 64 * square + 256 * n as u128 > budget.maximum_dense_bytes as u128
            || square * (n as u128 + 16) > budget.maximum_dense_work as u128
            || !nonnegative(budget.equilibrium_tolerance_w)
            || !positive(budget.relative_residual_tolerance)
            || budget.relative_residual_tolerance >= 1.0
        {
            return Err(invalid("invalid thermal impedance resources or tolerances"));
        }
        self.admit(thermal, &mut checkpoint)?;
        let endpoint = self.stationary_endpoint(thermal, &mut checkpoint)?;
        let rates = self.stationary_heat_rates(thermal, &endpoint, &mut checkpoint)?;
        let mut slopes = Vec::with_capacity(n);
        for (i, &id) in self.ids.iter().enumerate() {
            checkpoint()?;
            let net_power_w = rates[i];
            if net_power_w.abs() > budget.equilibrium_tolerance_w {
                return Err(ThermalImpedanceError::NotEquilibrium {
                    cell: id,
                    net_power_w,
                });
            }
            let cell = &thermal.cells[&id];
            self.laws[&cell.material_id]
                .admit_harmonic_operating_point(cell.specific_h_j_kg)?;
            let [lo, hi] = self.laws[&cell.material_id].enthalpy_bounds_j_kg();
            if cell.specific_h_j_kg <= lo || cell.specific_h_j_kg >= hi {
                return Err(invalid(
                    "harmonic perturbations need an interior caloric operating point",
                ));
            }
            let slope = endpoint.states[i].temperature_slope_kg_k_j;
            let per_energy = finite(slope / cell.mass_kg)?;
            if slope != 0.0 && per_energy == 0.0 {
                return Err(invalid("caloric energy tangent is not representable"));
            }
            slopes.push(per_energy);
        }
        let mut workspace = ImplicitConductionWorkspace::new(n)?;
        self.assemble_tangent(
            thermal,
            &endpoint,
            1.0,
            0.0,
            &mut workspace,
            &mut checkpoint,
        )?;
        let mut tangent = vec![0.0; n * n];
        for (i, &diagonal) in workspace.diagonal.iter().enumerate() {
            checkpoint()?;
            tangent[i * n + i] = diagonal;
        }
        if let Some(entries) = &workspace.graph_entries {
            for &(row, column, value) in entries {
                checkpoint()?;
                let entry = &mut tangent[row * n + column];
                *entry = finite(*entry + value)?;
            }
        } else {
            for i in 0..n - 1 {
                checkpoint()?;
                tangent[(i + 1) * n + i] = workspace.lower[i];
                tangent[i * n + i + 1] = workspace.upper[i];
            }
        }
        checkpoint()?;
        // The operating temperatures are the ones the stationary rates were taken at.
        let temperatures_k = self
            .ids
            .iter()
            .map(|&id| thermal.temperature_k(id))
            .collect::<core::result::Result<Vec<_>, _>>()?;
        Ok(ThermalImpedanceModel {
            ids: self.ids.clone(),
            time_s: thermal.time_s,
            temperatures_k,
            temperature_per_energy_k_j: slopes,
            outward_tangent_per_s: tangent,
            budget,
        })
    }
}

impl ThermalImpedanceModel {
    pub fn cell_ids(&self) -> &[CellId] {
        &self.ids
    }
    pub fn operating_time_s(&self) -> f64 {
        self.time_s
    }
    pub fn operating_temperatures_k(&self) -> &[f64] {
        &self.temperatures_k
    }

    /// Factor once; DC nullspaces refuse without artificial grounding/leakage.
    pub fn factor(
        &self,
        omega_rad_s: f64,
        mut checkpoint: impl FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<ThermalImpedanceFactor<'_>> {
        checkpoint()?;
        if !nonnegative(omega_rad_s) {
            return Err(invalid(
                "thermal angular frequency must be finite and nonnegative",
            ));
        }
        let n = self.ids.len();
        if omega_rad_s == 0.0 {
            self.admit_dc(&mut checkpoint)?;
        }
        let mut row_scales = Vec::with_capacity(n);
        let mut matrix = vec![C64::ZERO; n * n];
        for i in 0..n {
            checkpoint()?;
            let row = &self.outward_tangent_per_s[i * n..(i + 1) * n];
            let scale = row.iter().fold(omega_rad_s, |s, value| s.max(value.abs()));
            if scale == 0.0 {
                return Err(LuComplexError::Singular { index: i }.into());
            }
            row_scales.push(scale);
            for j in 0..n {
                matrix[i * n + j] = C64::new(
                    row[j] / scale,
                    if i == j { -omega_rad_s / scale } else { 0.0 },
                );
            }
        }
        self.factor_scaled_matrix(omega_rad_s, row_scales, matrix, &mut checkpoint)
    }

    // Both ordinary and closed-energy response use this original native LU.
    fn factor_scaled_matrix(
        &self,
        omega_rad_s: f64,
        row_scales: Vec<f64>,
        matrix: Vec<C64>,
        checkpoint: &mut impl FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<ThermalImpedanceFactor<'_>> {
        let n = self.ids.len();
        let mut refusal = None;
        // Keep the declared dense budget independent of an enclosing thread
        // budget; this bounded preparation uses the native serial factor.
        let factor = with_lu_workers(1, || {
            lu_complex_with_cancel(&matrix, n, || match checkpoint() {
                Ok(()) => false,
                Err(error) => {
                    refusal = Some(error);
                    true
                }
            })
        });
        if let Some(error) = refusal {
            return Err(error.into());
        }
        let lu = factor?;
        checkpoint()?;
        Ok(ThermalImpedanceFactor {
            model: self,
            omega_rad_s,
            row_scales,
            lu,
        })
    }

    // Detect a component's left constant energy null vector before LU; a
    // rounded last pivot must not manufacture finite DC impedance.
    fn admit_dc(
        &self,
        checkpoint: &mut impl FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<()> {
        let n = self.ids.len();
        let mut seen = vec![false; n];
        for first in 0..n {
            if seen[first] {
                continue;
            }
            seen[first] = true;
            let mut pending = vec![first];
            let mut anchored = false;
            while let Some(column) = pending.pop() {
                checkpoint()?;
                let (mut sum, mut magnitude) = (Sum::default(), Sum::default());
                for row in 0..n {
                    let value = self.outward_tangent_per_s[row * n + column];
                    sum.add(value)?;
                    magnitude.add(value.abs())?;
                    if !seen[row]
                        && (value != 0.0 || self.outward_tangent_per_s[column * n + row] != 0.0)
                    {
                        seen[row] = true;
                        pending.push(row);
                    }
                }
                anchored |= sum.total()?.abs() > 64.0 * f64::EPSILON * magnitude.total()?;
            }
            if !anchored {
                return Err(invalid(
                    "DC thermal impedance has an unanchored energy nullspace",
                ));
            }
        }
        Ok(())
    }
}

impl ThermalImpedanceFactor<'_> {
    pub fn omega_rad_s(&self) -> f64 {
        self.omega_rad_s
    }

    /// Additional peak watts in cell order, with original-row residual checks.
    /// Reuses the LU; nominal sources are not added again.
    pub fn response(
        &self,
        heat_w: &[C64],
        checkpoint: impl FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<ThermalImpedanceResponse> {
        self.response_with_rhs(heat_w, heat_w, checkpoint)
    }

    // Constrained component rows change only the solved RHS. Admission below
    // always rechecks every original cell equation against the actual heat.
    fn response_with_rhs(
        &self,
        heat_w: &[C64],
        rhs: &[C64],
        mut checkpoint: impl FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<ThermalImpedanceResponse> {
        checkpoint()?;
        let n = self.model.ids.len();
        if heat_w.len() != n || rhs.len() != n {
            return Err(invalid("one heat phasor per original cell required"));
        }
        let mut energy = Vec::with_capacity(n);
        for ((&value, &heat), &scale) in rhs.iter().zip(heat_w).zip(&self.row_scales) {
            checkpoint()?;
            complex(value)?;
            complex(heat)?;
            energy.push(complex(C64::new(value.re / scale, value.im / scale))?);
        }
        let mut refusal = None;
        let solved = self
            .lu
            .solve_with_cancel(&mut energy, || match checkpoint() {
                Ok(()) => false,
                Err(error) => {
                    refusal = Some(error);
                    true
                }
            });
        if let Some(error) = refusal {
            return Err(error.into());
        }
        solved?;
        let mut temperatures = Vec::with_capacity(n);
        let (mut maximum, mut relative) = (0.0_f64, 0.0_f64);
        for i in 0..n {
            checkpoint()?;
            complex(energy[i])?;
            temperatures.push(complex(
                energy[i].scale(self.model.temperature_per_energy_k_j[i]),
            )?);
            let storage = complex(C64::new(
                self.omega_rad_s * energy[i].im,
                -self.omega_rad_s * energy[i].re,
            ))?;
            let (mut real, mut imag, mut scale) = (Sum::default(), Sum::default(), Sum::default());
            real.add(storage.re)?;
            imag.add(storage.im)?;
            scale.add(storage.abs())?;
            for (j, &e) in energy.iter().enumerate() {
                let flux = complex(e.scale(self.model.outward_tangent_per_s[i * n + j]))?;
                real.add(flux.re)?;
                imag.add(flux.im)?;
                scale.add(flux.abs())?;
            }
            real.add(-heat_w[i].re)?;
            imag.add(-heat_w[i].im)?;
            scale.add(heat_w[i].abs())?;
            let residual = finite(C64::new(real.total()?, imag.total()?).abs())?;
            let scale = scale.total()?;
            maximum = maximum.max(residual);
            relative = relative.max(if scale == 0.0 { 0.0 } else { residual / scale });
        }
        if relative > self.model.budget.relative_residual_tolerance {
            return Err(ThermalImpedanceError::Residual {
                relative,
                tolerance: self.model.budget.relative_residual_tolerance,
            });
        }
        checkpoint()?;
        Ok(ThermalImpedanceResponse {
            energy_j: energy,
            temperature_k: temperatures,
            maximum_residual_w: maximum,
            relative_residual: relative,
        })
    }

    /// One physical driving/mutual impedance column [K/W], from unit heat at
    /// the named original cell. Repeated columns share this frequency's LU.
    pub fn impedance_column(
        &self,
        source: CellId,
        mut checkpoint: impl FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<Vec<C64>> {
        checkpoint()?;
        // Both native owners take sorted IDs from the original BTreeMap.
        let column = self
            .model
            .ids
            .binary_search(&source)
            .map_err(|_| invalid("unknown thermal impedance source cell"))?;
        let mut heat = vec![C64::ZERO; self.model.ids.len()];
        heat[column] = C64::ONE;
        Ok(self.response(&heat, checkpoint)?.temperature_k)
    }
}
