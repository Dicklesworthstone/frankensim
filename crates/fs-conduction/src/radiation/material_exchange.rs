//! A closed material enclosure as one coupled heat source of the implicit
//! thermal graph, with its exact enthalpy tangent and a declared passive
//! exchange.

use fs_material::state_point::MaterialPropertySelection;
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_thermal::enthalpy::conduction::implicit::coupled_source::{
    ImplicitCoupledHeatSourceLaw, ImplicitCoupledHeatSourceState,
};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation, ThermalLaw};

use super::pointwise::material::source::{DerivativeMode, MaterialEmissivitySource};
use super::prepared::PreparedRadiosity;
use super::{
    ConductionError, GrayDiffuseEnclosure, RadiosityPoll, RadiosityReport,
    STEFAN_BOLTZMANN_W_M2_K4, SurfaceEmissivity,
};
use crate::transient::deposition::radiation::{MaterialRadiation, RadiationCellError};

const SWITCHED: &str = "material radiation emissivity selection would leave its attached claim";
const EMISSIVITY: &str = "material radiation emissivity query refused at the trial state";
const RADIOSITY: &str = "material radiation radiosity solve refused";
const EQUATIONS: &str =
    "original radiosity equations or closed-enclosure power balance did not close";
const NONFINITE: &str = "accumulated radiation cell power is nonfinite";
const PARTITION: &str = "surface-to-cell power partition did not close";

/// The enclosure of a [`MaterialRadiation`] binding as one
/// [`ImplicitCoupledHeatSourceLaw`] over its radiating cells, in `CellId`
/// order ([`Self::cells`]).
///
/// Each query resolves every surface's emissivity at its cell's specific
/// enthalpy, solves the original radiosity equations once and applies
/// [`MaterialRadiation::resolve`]'s gates within the binding's tolerances:
/// the radiosity residual, the enclosure's power closure and the
/// surface-to-cell partition. The emissivity comes from the claim that the
/// binding's selection picks at attachment. A query is refused where
/// `resolve`'s per-query selection would pick another claim, or refuse.
/// Otherwise the cell powers are `resolve`'s at the same state.
///
/// The tangent is exact, with no differences. For each cell `c`, let
/// `r_s = (de_s/dh_c) (B_s - G_s) + e_s dB_s/dh_c` on `c`'s surfaces and zero
/// elsewhere. Here `B = sigma T^4` is the emission, `G = F J` the
/// irradiation, and `de_s/dh_c` is the emissivity source's own segment slope
/// (one-sided at a knot). Then `dJ = M^-1 r` on the same factors, and the
/// derivative of the outward heat `A_s e_s (B_s - G_s)` is
/// `A_s (r_s - e_s (F dJ)_s)`. The tangent steers Newton only; acceptance
/// stays on the residuals.
///
/// The enclosure is closed, so the powers are declared a passive exchange
/// within the binding's power tolerance
/// ([`ImplicitCoupledHeatSourceLaw::passive_exchange_tolerance_w`]).
///
/// Refusals: a material query refused at a trial state is
/// [`DepositionError::MaterialDomain`], so the line search backs off. A
/// switching selection, a refused radiosity solve and `resolve`'s gates are
/// [`DepositionError::InvalidInput`], each with its own message. A
/// checkpoint's own error is returned unchanged.
pub struct RadiationExchange<'a> {
    enclosure: &'a GrayDiffuseEnclosure,
    surfaces: Vec<ExchangeSurface<'a>>,
    cells: Vec<CellId>,
    calorics: Vec<&'a dyn ThermalLaw>,
    dependencies: Vec<(usize, usize)>,
    power_tolerance_w: f64,
    radiosity_tolerance_w_m2: f64,
}

struct ExchangeSurface<'a> {
    source: MaterialEmissivitySource<'a>,
    material: ThermalMaterialPoint<'a>,
    /// The binding's selection, or `None` when it pins its claim.
    selection: Option<MaterialPropertySelection>,
    cell: usize,
}

impl<'a> MaterialRadiation<'a> {
    /// This enclosure as an implicit coupled heat source ([`RadiationExchange`]).
    ///
    /// The admission is `resolve`'s, made once here against `thermal`: complete
    /// bindings, active cells with their original caloric objects, and traces
    /// within each cell's exposed area. A later event that changes an exposed
    /// area needs a fresh attachment. Each surface's emissivity source is
    /// pinned to the claim its binding's selection picks at the binding's
    /// material point.
    pub fn exchange(
        &self,
        thermal: &DepositionSimulation<'_>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<RadiationExchange<'a>, RadiationCellError> {
        checkpoint()?;
        let ordered = self.admitted(thermal, &mut checkpoint)?;
        let mut cells: Vec<CellId> = ordered.iter().map(|binding| binding.cell).collect();
        cells.sort_unstable();
        cells.dedup();
        let mut calorics: Vec<Option<&'a dyn ThermalLaw>> = vec![None; cells.len()];
        let mut surfaces = Vec::with_capacity(ordered.len());
        for (index, binding) in ordered.into_iter().enumerate() {
            checkpoint()?;
            let cell = cells.partition_point(|&id| id < binding.cell);
            calorics[cell] = Some(binding.curve as &dyn ThermalLaw);
            let selection = match &binding.selection {
                MaterialPropertySelection::PinnedByProperty(_) => None,
                other => Some(other.clone()),
            };
            surfaces.push(ExchangeSurface {
                source: MaterialEmissivitySource::prepare(
                    self.enclosure.surfaces()[index].name(),
                    binding.curve,
                    binding.material,
                    binding.selection.clone(),
                )?,
                material: binding.material,
                selection,
                cell,
            });
        }
        let n = cells.len();
        Ok(RadiationExchange {
            enclosure: self.enclosure,
            surfaces,
            cells,
            calorics: calorics.into_iter().flatten().collect(),
            dependencies: (0..n)
                .flat_map(|row| (0..n).map(move |col| (row, col)))
                .collect(),
            power_tolerance_w: self.power_tolerance_w,
            radiosity_tolerance_w_m2: self.radiosity_tolerance_w_m2,
        })
    }
}

impl RadiationExchange<'_> {
    /// The radiating cells in `CellId` order: the attachment's cell order.
    #[must_use]
    pub fn cells(&self) -> &[CellId] {
        &self.cells
    }

    fn query(
        &self,
        h: &[f64],
        poll: &mut RadiosityPoll<'_>,
    ) -> Result<ImplicitCoupledHeatSourceState, DepositionError> {
        if h.len() != self.cells.len() {
            return Err(DepositionError::InvalidInput(
                "one specific enthalpy per radiating cell",
            ));
        }
        let n = self.surfaces.len();
        let mut emissivities = Vec::with_capacity(n);
        let mut temperatures = Vec::with_capacity(n);
        let mut slopes = Vec::with_capacity(n);
        for (index, surface) in self.surfaces.iter().enumerate() {
            poll("radiation-exchange-surface", index).map_err(|e| refusal(&e, EMISSIVITY))?;
            let node = surface
                .source
                .resolve(index, h[surface.cell], None, DerivativeMode::Primal)
                .map_err(|e| refusal(&e, EMISSIVITY))?;
            if let Some(selection) = &surface.selection {
                let (selected, _) = SurfaceEmissivity::from_thermal_state(
                    self.enclosure.surfaces()[index].name(),
                    surface.material,
                    node.state.phase,
                    selection.clone(),
                )
                .map_err(|_| DepositionError::InvalidInput(SWITCHED))?;
                if selected.receipt().selected != surface.source.claim() {
                    return Err(DepositionError::InvalidInput(SWITCHED));
                }
            }
            temperatures.push(node.state.phase.temperature_k());
            slopes.push((node.depsilon_dh, node.dt_dh));
            emissivities.push(node.state.emissivity);
        }
        let enclosure = self
            .enclosure
            .with_emissivities(emissivities)
            .map_err(|e| refusal(&e, RADIOSITY))?;
        let prepared =
            PreparedRadiosity::new_polled(poll, &enclosure).map_err(|e| refusal(&e, RADIOSITY))?;
        let report = prepared
            .solve_polled(poll, &temperatures)
            .map_err(|e| refusal(&e, RADIOSITY))?;
        let powers_w = self.powers(&report)?;
        let cells = self.cells.len();
        let factors = &enclosure.view_factors.factors;
        let mut derivatives = vec![0.0; cells * cells];
        for column in 0..cells {
            let mut rhs = vec![0.0; n];
            for (s, surface) in self.surfaces.iter().enumerate() {
                if surface.cell == column {
                    let emissivity = enclosure.surfaces()[s].emissivity().value();
                    let (depsilon_dh, dt_dh) = slopes[s];
                    let t = temperatures[s];
                    let emission_slope = 4.0 * STEFAN_BOLTZMANN_W_M2_K4 * t * t * t * dt_dh;
                    rhs[s] = finite(
                        depsilon_dh * report.net_outward_flux_w_m2[s] / emissivity
                            + emissivity * emission_slope,
                    )?;
                }
            }
            if rhs.iter().all(|&value| value == 0.0) {
                continue;
            }
            let dj = prepared
                .solve_rhs_polled(poll, &rhs)
                .map_err(|e| refusal(&e, RADIOSITY))?;
            for (s, surface) in self.surfaces.iter().enumerate() {
                let irradiation: f64 = factors[s].iter().zip(&dj).map(|(f, d)| f * d).sum();
                let emissivity = enclosure.surfaces()[s].emissivity().value();
                let outward =
                    enclosure.surfaces()[s].area_m2() * (rhs[s] - emissivity * irradiation);
                let entry = &mut derivatives[surface.cell * cells + column];
                *entry = finite(*entry - outward)?;
            }
        }
        Ok(ImplicitCoupledHeatSourceState {
            powers_w,
            power_enthalpy_derivatives_w_kg_j: derivatives,
        })
    }

    /// `resolve`'s gates and cell sums, in its order.
    fn powers(&self, report: &RadiosityReport) -> Result<Vec<f64>, DepositionError> {
        if report
            .net_outward_heat_w
            .iter()
            .chain(&report.radiosity_w_m2)
            .chain(&report.irradiation_w_m2)
            .any(|v| !v.is_finite())
            || !report.enclosure_energy_closure_w.is_finite()
            || !nonnegative(report.linear_residual_max_w_m2)
            || report.linear_residual_max_w_m2 > self.radiosity_tolerance_w_m2
            || report.enclosure_energy_closure_w.abs() > self.power_tolerance_w
        {
            return Err(DepositionError::InvalidInput(EQUATIONS));
        }
        let mut powers = vec![0.0; self.cells.len()];
        for (surface, &outward) in self.surfaces.iter().zip(&report.net_outward_heat_w) {
            powers[surface.cell] -= outward;
            if !powers[surface.cell].is_finite() {
                return Err(DepositionError::InvalidInput(NONFINITE));
            }
        }
        let total: f64 = powers.iter().copied().sum();
        if !total.is_finite()
            || total.abs() > self.power_tolerance_w
            || (total + report.enclosure_energy_closure_w).abs() > self.power_tolerance_w
        {
            return Err(DepositionError::InvalidInput(PARTITION));
        }
        Ok(powers)
    }
}

impl ImplicitCoupledHeatSourceLaw for RadiationExchange<'_> {
    fn caloric_laws(&self) -> &[&dyn ThermalLaw] {
        &self.calorics
    }

    fn power_dependencies(&self) -> &[(usize, usize)] {
        &self.dependencies
    }

    fn state_at_specific_enthalpies_and_time(
        &self,
        specific_h_j_kg: &[f64],
        _time_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitCoupledHeatSourceState, DepositionError> {
        let mut stopped = None;
        let result = self.query(specific_h_j_kg, &mut |stage, at| {
            checkpoint().map_err(|error| {
                stopped = Some(error);
                ConductionError::Cancelled { stage, at }
            })
        });
        result.map_err(|error| stopped.take().unwrap_or(error))
    }

    fn passive_exchange_tolerance_w(&self) -> Option<f64> {
        Some(self.power_tolerance_w)
    }
}

fn refusal(error: &ConductionError, message: &'static str) -> DepositionError {
    match error {
        ConductionError::Cancelled { .. } => DepositionError::Cancelled,
        ConductionError::MaterialQuery { .. } => DepositionError::MaterialDomain,
        _ => DepositionError::InvalidInput(message),
    }
}

fn finite(value: f64) -> Result<f64, DepositionError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DepositionError::Arithmetic)
    }
}

fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
