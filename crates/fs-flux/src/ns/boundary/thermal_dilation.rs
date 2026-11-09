//! Caloric thermal expansion on the native BDM1/P0 Stokes constraint.
//!
//! At fixed thermodynamic pressure, mass continuity gives
//! `div(u) = (1/v) Dv/Dt = rho(h) (dv/dh) Dh/Dt`, where `v=1/rho`.
//! The original material chart supplies both factors, including its actual
//! sensible-density interpolation and variable heat capacity. The input rate
//! is a MATERIAL derivative, not an Eulerian local heating increment.

/// Original measured Cp/free-expansion sources on the same native constraint.
#[path = "thermal_dilation/isobaric.rs"]
pub mod isobaric;

use super::{StokesAssembly, StokesBoundary, StokesBoundaryError, StokesSolution};
use crate::ns::{FluxParams, FluxSystem};
use fs_material::phase::{
    EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState, PhaseStateError, SolidLiquidPhase,
};
use fs_math::dd::Dd;

/// One cell's original isobaric liquid chart and accepted caloric rate.
#[derive(Clone, Copy)]
pub struct CaloricDilationCell<'a> {
    pub caloric: &'a EquilibriumEnthalpyPhaseCurve,
    /// Current specific enthalpy in this chart's original datum [J/kg].
    pub specific_enthalpy_j_kg: f64,
    /// Material derivative Dh/Dt [W/kg]. It includes the thermal owner's
    /// heat/work convention; pressure work must not be added a second time.
    pub material_specific_enthalpy_rate_w_kg: f64,
}

/// Source identities and the exact chart tangent used by the constraint.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CaloricDilationResponse {
    pub state: EquilibriumPhaseState,
    pub material_specific_enthalpy_rate_w_kg: f64,
    pub specific_volume_derivative_m3_j: f64,
    /// Prescribed cellwise divergence [1/s]; negative means contraction.
    pub dilation_per_s: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum CaloricDilationError {
    Flow(StokesBoundaryError),
    Material {
        cell: usize,
        source: PhaseStateError,
    },
    Invalid(&'static str),
    /// A solid or coexistence state needs a separately admitted mobility law.
    NonLiquid {
        cell: usize,
    },
    /// Actual returned edge moments failed the independent physical budget.
    DilationDefect {
        cell: usize,
        expected_per_s: f64,
        actual_per_s: f64,
        maximum_defect_per_s: f64,
    },
}

impl core::fmt::Display for CaloricDilationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "caloric Stokes dilation: {self:?}")
    }
}
impl std::error::Error for CaloricDilationError {}
impl From<StokesBoundaryError> for CaloricDilationError {
    fn from(error: StokesBoundaryError) -> Self {
        Self::Flow(error)
    }
}

/// The original native operator with its thermodynamic continuity loads.
pub struct CaloricDilationAssembly<'a, 'm> {
    native: StokesAssembly<'a, 'm>,
    source: Vec<CaloricDilationResponse>,
}

/// Native pressure/velocity/reactions and independently measured dilation.
pub struct CaloricDilationSolution {
    pub flow: StokesSolution,
    pub source: Vec<CaloricDilationResponse>,
    pub maximum_dilation_defect_per_s: f64,
}

impl<'a, 'm> CaloricDilationAssembly<'a, 'm> {
    /// Original matrix and RHS remain available for checked native adjoints.
    /// They are not derivatives of the upstream thermodynamic source model.
    #[must_use]
    pub const fn native(&self) -> &StokesAssembly<'a, 'm> {
        &self.native
    }

    #[must_use]
    pub fn source(&self) -> &[CaloricDilationResponse] {
        &self.source
    }

    /// Solve the original mixed system, then read back EVERY cell's flux,
    /// including the continuity row replaced by a sealed component's gauge.
    /// The allowance is absolute in 1/s; no velocity, mass or density repair
    /// is performed if the actual source constraint fails it.
    pub fn solve(
        &self,
        maximum_dilation_defect_per_s: f64,
    ) -> Result<CaloricDilationSolution, CaloricDilationError> {
        let (flow, maximum_defect) = solve_dilation(
            &self.native,
            self.source.iter().map(|source| source.dilation_per_s),
            maximum_dilation_defect_per_s,
        )?;
        Ok(CaloricDilationSolution {
            flow,
            source: self.source.clone(),
            maximum_dilation_defect_per_s: maximum_defect,
        })
    }
}

fn solve_dilation(
    native: &StokesAssembly<'_, '_>,
    dilation: impl Iterator<Item = f64>,
    maximum_dilation_defect_per_s: f64,
) -> Result<(StokesSolution, f64), CaloricDilationError> {
    if !maximum_dilation_defect_per_s.is_finite() || maximum_dilation_defect_per_s < 0.0 {
        return Err(CaloricDilationError::Invalid(
            "finite nonnegative physical dilation allowance required",
        ));
    }
    let flow = native.solve()?;
    let mesh = native.system.mesh;
    let mut maximum_defect = 0.0_f64;
    for (cell, expected) in dilation.enumerate() {
        let mut outflow = Dd::from_f64(0.0);
        for &(edge, sign) in &mesh.tri_edges[cell] {
            outflow = outflow + Dd::from_f64(sign * mesh.edges[edge].len * flow.field.x[2 * edge]);
        }
        let actual = (outflow.hi + outflow.lo) / mesh.areas[cell];
        let defect = (actual - expected).abs();
        if !defect.is_finite() || defect > maximum_dilation_defect_per_s {
            return Err(CaloricDilationError::DilationDefect {
                cell,
                expected_per_s: expected,
                actual_per_s: actual,
                maximum_defect_per_s: maximum_dilation_defect_per_s,
            });
        }
        maximum_defect = maximum_defect.max(defect);
    }
    Ok((flow, maximum_defect))
}

impl<'m> FluxSystem<'m> {
    /// Quasistatic liquid flow driven by its original caloric expansion law.
    ///
    /// Each triangle receives `rho * d(1/rho)/dh * Dh/Dt` as its P0
    /// divergence. Positive and negative expansion and nonlinear Cp/rho laws
    /// keep the SAME source chart; its right-sided interior-knot and left-sided
    /// upper-endpoint tangent convention is unchanged. Outward rates at the
    /// caloric domain boundary and cooling at an interior knot refuse: the
    /// right-sided tangent there cannot supply a left-sided cooling derivative.
    /// Only liquid states are admitted. No mass source, independent beta fit or density change is
    /// invented by the flow operator.
    ///
    /// The existing symmetric `-div(2 nu D(u)) + grad(p) = f` operator,
    /// one-sided viscosity, boundary tractions, pressure gauges and native
    /// solver are reused. For this nonzero-divergence lane, `p` is the pressure
    /// multiplying `-I` in that declared stress; any separately modelled bulk
    /// stress must be included consistently by the composing physical owner.
    /// For a 3-D Newtonian fluid with planar velocity, physical kinematic
    /// pressure is `p + (zeta/rho_reference - 2*nu/3)*dilation`; the native
    /// total-traction boundary already acts on the complete declared stress.
    /// Viscosity, force and traction use the existing fixed-reference-density
    /// normalization. `params.nu` is unused, as in the native variable-viscosity
    /// entry. Temperature-dependent viscosity remains caller/source-owned.
    ///
    /// Every component with prescribed boundary velocities must satisfy
    /// `integral_boundary u.n = integral_volume dilation`. A rigid sealed
    /// heated cavity with nonzero net expansion refuses: mechanical pressure
    /// is not an EOS that can compress it. A traction boundary instead lets the
    /// actual native flow carry the demanded volume out, or draw volume in.
    ///
    /// This is an instantaneous fixed-geometry Stokes response, not a time
    /// update of material mass/enthalpy, a free-surface remap, or a closed
    /// thermal-mechanical energy solve. Local sources are P0 approximations;
    /// finite-time motion, material advection and source/mesh refinement remain
    /// with the composing owner. Work associated with pressure and dilation
    /// is not silently converted into heat.
    #[allow(clippy::too_many_arguments)]
    pub fn assemble_stokes_caloric_dilation<'a>(
        &'a self,
        params: FluxParams,
        viscosity: &dyn Fn(usize, [f64; 2]) -> f64,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[StokesBoundary<'_>],
        cells: &[CaloricDilationCell<'_>],
    ) -> Result<CaloricDilationAssembly<'a, 'm>, CaloricDilationError> {
        if cells.len() != self.mesh.tris.len() {
            return Err(CaloricDilationError::Invalid(
                "one original caloric state and material enthalpy rate per triangle required",
            ));
        }
        let mut source = Vec::with_capacity(cells.len());
        let mut cell_outflow = Vec::with_capacity(cells.len());
        for (cell, input) in cells.iter().enumerate() {
            let rate = input.material_specific_enthalpy_rate_w_kg;
            if !rate.is_finite() {
                return Err(CaloricDilationError::Invalid(
                    "nonfinite material enthalpy rate",
                ));
            }
            let material_error = |source| CaloricDilationError::Material { cell, source };
            let state = input
                .caloric
                .state_at_specific_enthalpy(input.specific_enthalpy_j_kg)
                .map_err(material_error)?;
            if state.phase() != SolidLiquidPhase::Liquid {
                return Err(CaloricDilationError::NonLiquid { cell });
            }
            let knots = input.caloric.knots();
            let knot =
                knots.partition_point(|k| k.specific_enthalpy_j_kg < input.specific_enthalpy_j_kg);
            if knot < knots.len()
                && knots[knot].specific_enthalpy_j_kg == input.specific_enthalpy_j_kg
                && ((rate < 0.0 && knot + 1 < knots.len())
                    || (rate > 0.0 && knot + 1 == knots.len()))
            {
                return Err(CaloricDilationError::Invalid(
                    "caloric rate requires an unsupported outward or left-sided knot tangent",
                ));
            }
            let slope = input
                .caloric
                .specific_volume_derivative_at_specific_enthalpy(input.specific_enthalpy_j_kg)
                .map_err(material_error)?;
            let coefficient = state.bulk_density_kg_m3() * slope;
            let dilation_per_s = coefficient * rate;
            let outflow = self.mesh.areas[cell] * dilation_per_s;
            if !coefficient.is_finite()
                || !dilation_per_s.is_finite()
                || !outflow.is_finite()
                || (slope != 0.0 && coefficient == 0.0)
                || (coefficient != 0.0 && rate != 0.0 && dilation_per_s == 0.0)
                || (dilation_per_s != 0.0 && outflow == 0.0)
            {
                return Err(CaloricDilationError::Invalid(
                    "unrepresentable thermal dilation",
                ));
            }
            source.push(CaloricDilationResponse {
                state,
                material_specific_enthalpy_rate_w_kg: rate,
                specific_volume_derivative_m3_j: slope,
                dilation_per_s,
            });
            cell_outflow.push(outflow);
        }
        let native = self.assemble_stokes_dilation_native(
            params,
            viscosity,
            force,
            boundaries,
            &cell_outflow,
        )?;
        Ok(CaloricDilationAssembly { native, source })
    }

    fn assemble_stokes_dilation_native<'a>(
        &'a self,
        params: FluxParams,
        viscosity: &dyn Fn(usize, [f64; 2]) -> f64,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[StokesBoundary<'_>],
        cell_outflow: &[f64],
    ) -> Result<StokesAssembly<'a, 'm>, CaloricDilationError> {
        let admitted =
            self.admit_stokes_problem_with_dilation(params, boundaries, &[], Some(cell_outflow))?;
        let mut native =
            self.assemble_stokes_admitted_viscosity(params, viscosity, force, &admitted)?;
        for (cell, &outflow) in cell_outflow.iter().enumerate() {
            let component = self.pressure_component[cell];
            if outflow != 0.0
                && (admitted.has_traction[component] || self.pressure_roots[component] != cell)
            {
                // Native continuity row is -integral div(u), not +div(u).
                native.rhs[self.n_u + cell] = -outflow;
            }
        }
        Ok(native)
    }
}
