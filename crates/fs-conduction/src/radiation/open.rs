//! Diffuse-gray exchange with an explicitly black external environment.
//!
//! A surface row may sum to less than one: its deficit is the view of the
//! black bath, not heat to redistribute among the retained surfaces. Complete
//! the reciprocal system with a black bath of effective exchange area
//! `A_b = sum_i A_i F_ib` and `F_bi = A_i F_ib/A_b`. Its radiosity is fixed by
//! its temperature. This algebraic completion does not invent a mesh surface,
//! close an aperture, or infer the temperature of an unassigned opaque wall.
//!
//! The original shifted radiosity kernel and prepared factorization remain
//! the numerical owners. A patch is isothermal at its supplied temperature;
//! a solid-field binding uses its area mean. The bath is diffuse, isotropic
//! and black, with no spectral, participating-medium or directional-source
//! claim. See the net-radiation/view-factor equations in the MOOSE Heat
//! Transfer module: <https://mooseframework.inl.gov/modules/heat_transfer/>.

use std::collections::BTreeSet;

use fs_exec::Cx;
use fs_matdb::UncertaintyModel;

use super::{
    CoupledRadiationConfig, CoupledRadiationReport, CoupledRadiationSolution,
    MAX_ENCLOSURE_SURFACES, RadiationSurface, RadiosityReport, ViewFactorEvidence,
    ViewFactorMatrix, ViewFactorTolerance, enclosure_exchange as exchange,
    prepared::DenseRadiosityFactorization, radiation_error, require_temperature,
    solve_radiating_surfaces, solve_radiosity_system,
};
use crate::{
    ConductionError, ConductionProblem, ConductionSolution, RobinFlux, SolveConfig,
    ThermalBoundary, ThermalInterfaces,
};

/// Fixed reciprocal surface geometry with one explicitly black open bath.
/// Matrix preparation is bounded/cancellable and performed once. Temperature
/// evaluations thereafter reuse the factorization with O(n²) work and storage.
#[derive(Debug)]
pub struct OpenDiffuseEnclosure {
    surfaces: Vec<RadiationSurface>,
    complete: ViewFactorMatrix,
    factors: DenseRadiosityFactorization,
}

/// Original radiosity equations for physical surfaces FOLLOWED by the bath.
/// The last row is an algebraic black reservoir, never a solid thermal degree
/// of freedom. Its heat is independently evaluated, not minus a forced sum.
#[derive(Debug, Clone, PartialEq)]
pub struct OpenRadiosityReport {
    pub exchange: RadiosityReport,
}

impl OpenRadiosityReport {
    /// Signed outward powers of the physical surfaces, in declaration order.
    pub fn surface_heat_w(&self) -> &[f64] {
        let n = self.exchange.net_outward_heat_w.len() - 1;
        &self.exchange.net_outward_heat_w[..n]
    }

    /// Heat leaving the ambient reservoir, W. Negative means the ambient
    /// absorbs energy from the solids. It is separate from convective heat.
    pub fn ambient_heat_w(&self) -> f64 {
        self.exchange.net_outward_heat_w[self.exchange.net_outward_heat_w.len() - 1]
    }
}

/// A frozen radiation iterate composed with the ORIGINAL natural boundary.
/// This does not advance the solid or assert nonlinear convergence. A coupled
/// owner must repeat until its actual temperature/residual/heat criteria pass.
#[derive(Debug, Clone)]
pub struct OpenRadiationBoundary {
    pub boundary: ThermalBoundary,
    pub radiation: OpenRadiosityReport,
}

/// Nonlinear solid equilibrium with the original external heat channels.
#[derive(Debug, Clone, PartialEq)]
pub struct OpenEnclosureSolution {
    pub conduction: ConductionSolution,
    /// Applied and fresh final-state radiosity. Every exchange vector contains
    /// the physical surfaces followed by the explicit black ambient reservoir.
    pub radiation: CoupledRadiationReport,
    pub combined_boundary: ThermalBoundary,
    pub convective_robin_fluxes: Vec<RobinFlux>,
    /// Original convection only, W outward; this is the heat delivered to air.
    pub convective_out_w: f64,
    /// Original authored Neumann heat, excluding radiation, W outward.
    pub prescribed_neumann_out_w: f64,
    /// Sum of fresh physical-surface radiative heat, W outward.
    pub radiative_out_w: f64,
    /// Independently evaluated heat absorbed by the ambient bath, W.
    pub ambient_absorbed_w: f64,
}

impl OpenDiffuseEnclosure {
    /// Admit reciprocal physical-surface factors whose row deficits see one
    /// black environment, and prepare the shared radiosity matrix once.
    ///
    /// At most 255 physical surfaces leave room for the bath in the existing
    /// 256-row resource cap. Negative/overfull rows, absent escape, overlapping
    /// traces, mixed meshes and nonreciprocal exchange refuse. Tolerances never
    /// authorize clipping a negative bath view or redistributing an escape.
    pub fn new(
        cx: &Cx<'_>,
        surfaces: Vec<RadiationSurface>,
        surface_factors: Vec<Vec<f64>>,
        evidence: ViewFactorEvidence,
        tolerance: ViewFactorTolerance,
    ) -> Result<Self, ConductionError> {
        poll(cx, 0)?;
        let n = surfaces.len();
        if n == 0
            || n >= MAX_ENCLOSURE_SURFACES
            || surface_factors.len() != n
            || surface_factors.iter().any(|row| row.len() != n)
        {
            return Err(invalid(
                "open exchange requires 1..=255 surfaces and a square surface matrix",
            ));
        }
        let mut names = BTreeSet::new();
        let mut slots = BTreeSet::new();
        let mesh_identity = surfaces[0].mesh_identity;
        let mut areas = Vec::with_capacity(n + 1);
        let mut completed = Vec::with_capacity(n + 1);
        let mut ambient_exchange = Vec::with_capacity(n);
        let mut ambient_area = 0.0;
        for (i, (surface, row)) in surfaces.iter().zip(surface_factors).enumerate() {
            poll(cx, i)?;
            if surface.mesh_identity != mesh_identity
                || !names.insert(surface.name())
                || surface.face_slots().iter().any(|&slot| !slots.insert(slot))
            {
                return Err(invalid(
                    "open surfaces need distinct names and disjoint traces on one exact mesh",
                ));
            }
            if row
                .iter()
                .any(|&f| !f.is_finite() || !(0.0..=1.0).contains(&f))
            {
                return Err(invalid("surface view factors must be finite and in [0,1]"));
            }
            let sum: f64 = row.iter().sum();
            if sum > 1.0 {
                return Err(invalid(
                    "surface view-factor row exceeds one; no negative ambient view is admitted",
                ));
            }
            let escape = 1.0 - sum;
            let exchange_area = finite(surface.area_m2() * escape)?;
            if escape > 0.0 && exchange_area == 0.0 {
                return Err(invalid("positive ambient exchange area underflows"));
            }
            ambient_area = finite(ambient_area + exchange_area)?;
            ambient_exchange.push(exchange_area);
            areas.push(surface.area_m2());
            let mut full = row;
            full.push(escape);
            completed.push(full);
        }
        if ambient_area <= 0.0 {
            return Err(invalid(
                "open exchange has no ambient view; use GrayDiffuseEnclosure for a closed system",
            ));
        }
        let mut ambient_row: Vec<_> = ambient_exchange
            .iter()
            .map(|area| area / ambient_area)
            .collect();
        for (&area, &factor) in ambient_exchange.iter().zip(&ambient_row) {
            if area > 0.0 && factor == 0.0 {
                return Err(invalid("positive reciprocal ambient view underflows"));
            }
        }
        ambient_row.push(0.0);
        completed.push(ambient_row);
        areas.push(ambient_area);
        let complete = ViewFactorMatrix::admit(areas, completed, evidence, tolerance)?;
        let mut matrix = vec![vec![0.0; n + 1]; n + 1];
        for (i, row) in matrix.iter_mut().enumerate() {
            poll(cx, i)?;
            let epsilon = if i == n {
                1.0
            } else {
                surfaces[i].emissivity().value()
            };
            for (j, entry) in row.iter_mut().enumerate() {
                *entry =
                    if i == j { 1.0 } else { 0.0 } - (1.0 - epsilon) * complete.factors()[i][j];
            }
        }
        let factors = DenseRadiosityFactorization::new(cx, &matrix)?;
        poll(cx, n)?;
        Ok(Self {
            surfaces,
            complete,
            factors,
        })
    }

    pub fn surfaces(&self) -> &[RadiationSurface] {
        &self.surfaces
    }

    /// Completed reciprocal system; the final row is the explicit black bath.
    /// Its area is an effective exchange area, not a meshed reservoir area.
    pub fn view_factors(&self) -> &ViewFactorMatrix {
        &self.complete
    }

    /// Evaluate changed surface AND ambient temperatures with the same matrix.
    /// Only the solid temperatures query emissivity validity; the bath has its
    /// own absolute temperature and is not constrained by a surface material.
    pub fn solve(
        &self,
        cx: &Cx<'_>,
        surface_temperatures_k: &[f64],
        ambient_temperature_k: f64,
    ) -> Result<OpenRadiosityReport, ConductionError> {
        poll(cx, 0)?;
        let n = self.surfaces.len();
        if surface_temperatures_k.len() != n {
            return Err(invalid(
                "supply one temperature per physical surface, followed separately by the bath temperature",
            ));
        }
        require_temperature(
            "open-ambient",
            "black ambient temperature",
            ambient_temperature_k,
        )?;
        for (i, (surface, &temperature)) in
            self.surfaces.iter().zip(surface_temperatures_k).enumerate()
        {
            poll(cx, i)?;
            surface.emissivity.validate_temperature(
                surface.name(),
                "open surface temperature",
                temperature,
            )?;
        }
        let mut temperatures = surface_temperatures_k.to_vec();
        temperatures.push(ambient_temperature_k);
        let exchange = solve_radiosity_system(
            cx,
            &temperatures,
            |i| {
                if i == n {
                    1.0
                } else {
                    self.surfaces[i].emissivity.value
                }
            },
            |i| self.complete.areas_m2()[i],
            &self.complete,
            self.surfaces
                .iter()
                .all(|s| !matches!(s.emissivity.uncertainty, UncertaintyModel::Unstated)),
            Some(&self.factors),
        )?;
        for &value in exchange
            .radiosity_w_m2
            .iter()
            .chain(&exchange.irradiation_w_m2)
        {
            if !value.is_finite() || value < 0.0 {
                return Err(invalid("nonphysical open-system radiosity or irradiation"));
            }
        }
        for &value in exchange
            .net_outward_flux_w_m2
            .iter()
            .chain(&exchange.net_outward_heat_w)
        {
            finite(value)?;
        }
        finite(exchange.enclosure_energy_closure_w)?;
        finite(exchange.enclosure_energy_scale_w)?;
        finite(exchange.linear_residual_max_w_m2)?;
        poll(cx, n)?;
        Ok(OpenRadiosityReport { exchange })
    }

    /// Compose a frozen radiative iterate with the original uniform natural
    /// traces for an existing conduction solve. Neumann heat is ADDED; Robin
    /// heat is represented by `T_ref-q_rad/h`, an algebraic reference rather
    /// than a fluid temperature. Callers retain the original convection law
    /// for fluid heat accounting and audit the returned radiative heat apart.
    pub fn boundary_with_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        temperature: &[f64],
        ambient_temperature_k: f64,
    ) -> Result<OpenRadiationBoundary, ConductionError> {
        let patches = exchange::bind_surface_traces(cx, problem, &self.surfaces)?;
        let mut means = Vec::with_capacity(self.surfaces.len());
        for (i, surface) in self.surfaces.iter().enumerate() {
            poll(cx, i)?;
            means.push(surface.mean_temperature(problem.mesh, temperature)?);
        }
        let radiation = self.solve(cx, &means, ambient_temperature_k)?;
        let boundary = exchange::shifted_boundary(problem, &patches, &radiation.exchange)?;
        poll(cx, self.surfaces.len())?;
        Ok(OpenRadiationBoundary {
            boundary,
            radiation,
        })
    }
}

/// Solve coupled steady conduction and open-cavity radiation on 1..=64 named
/// natural traces. Reflections, fixed contact and heterogeneous k(T) participate
/// in every accepted solid solve. The SAME existing enclosure fixed point and
/// fresh nonlinear residual gate govern closed and open radiation.
///
/// Original convection, prescribed flux and heat absorbed by the bath are
/// separate channels. An initial conduction solution must be thermally anchored
/// before radiation is applied, as in the existing closed-enclosure driver.
/// This is fixed-emissivity patch-mean radiation, not a pointwise radiation
/// field, a radiation-only nullspace remedy, or a coupled adjoint.
#[allow(clippy::too_many_arguments)]
pub fn solve_with_open_gray_diffuse_enclosure(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    interfaces: Option<&ThermalInterfaces>,
    enclosure: &OpenDiffuseEnclosure,
    ambient_temperature_k: f64,
    conduction_config: SolveConfig,
    coupling_config: CoupledRadiationConfig,
) -> Result<OpenEnclosureSolution, ConductionError> {
    poll(cx, 0)?;
    require_temperature(
        "open-ambient",
        "black ambient temperature",
        ambient_temperature_k,
    )?;
    let patches = exchange::bind_surface_traces(cx, problem, &enclosure.surfaces)?;
    let prescribed_neumann_out_w = exchange::prescribed_neumann_heat(cx, problem)?;
    let apply_boundary =
        |report: &RadiosityReport| exchange::shifted_boundary(problem, &patches, report);
    let evaluate = |temperatures: &[f64]| {
        enclosure
            .solve(cx, temperatures, ambient_temperature_k)
            .map(|report| report.exchange)
    };
    let (
        CoupledRadiationSolution {
            conduction,
            radiation,
        },
        combined_boundary,
    ) = solve_radiating_surfaces(
        cx,
        problem,
        interfaces,
        &enclosure.surfaces,
        conduction_config,
        coupling_config,
        &apply_boundary,
        &evaluate,
        true,
    )?;
    let (convective_robin_fluxes, convective_out_w) =
        exchange::convection(cx, problem, &conduction.report.robin_fluxes)?;
    let heat = &radiation.equilibrium.nonlinear_radiosity.net_outward_heat_w;
    let radiative_out_w = finite(heat[..enclosure.surfaces.len()].iter().sum())?;
    let ambient_absorbed_w = -heat[enclosure.surfaces.len()];
    poll(cx, enclosure.surfaces.len())?;
    Ok(OpenEnclosureSolution {
        conduction,
        radiation,
        combined_boundary,
        convective_robin_fluxes,
        convective_out_w,
        prescribed_neumann_out_w,
        radiative_out_w,
        ambient_absorbed_w,
    })
}

fn finite(value: f64) -> Result<f64, ConductionError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("nonfinite open-radiation arithmetic"))
    }
}

fn poll(cx: &Cx<'_>, at: usize) -> Result<(), ConductionError> {
    cx.checkpoint().map_err(|_| ConductionError::Cancelled {
        stage: "open-radiation",
        at,
    })
}

fn invalid(what: &str) -> ConductionError {
    radiation_error(
        "open-enclosure",
        what,
        "retain reciprocal surface factors, explicit black ambient exposure, physical material support and finite heat accounts",
    )
}
