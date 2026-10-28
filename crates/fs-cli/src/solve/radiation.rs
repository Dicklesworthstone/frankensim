//! Canonical-project lowering for the shared area-mean surface-radiation model.
//!
//! The conduction report accounts for convection plus radiation. Only the
//! original convection rows may feed the stream-wise air energy equation.

use fs_conduction::{
    AmbientRadiationConfig, AmbientRadiationPatch, AmbientRadiationReport, ConductionError,
    ConductionProblem, ConductionSolution, RobinFlux, SolveConfig, SurfaceEmissivity,
    ThermalInterfaces,
};
use fs_exec::Cx;
use fs_conduction::radiation::GrayDiffuseEnclosure;
use fs_conduction::transient::enthalpy::enclosure::EnthalpyEnclosureReport;
use fs_matdb::SelectionPolicy;
use fs_project::ProjectSpec;
use fs_project::spec::{ConductionSetup, ThermalBoundaryCondition};

use super::{SolveRefusal, conduction_error, parse_claim_id};
use crate::cards::CardPackSet;
use crate::import::json_string;

pub(super) mod adaptive;
mod enclosure;
mod pointwise;
use pointwise::PointwiseEndpoint;
pub(super) use pointwise::receipt as pointwise_receipt;

pub(super) struct LoweredRadiation {
    pub(super) patches: Vec<AmbientRadiationPatch>,
    pub(super) config: AmbientRadiationConfig,
    surfaces: Vec<SurfaceSource>,
    enclosure: Option<enclosure::LoweredEnclosure>,
    pointwise: bool,
}

struct SurfaceSource {
    name: String,
    target: String,
    emissivity: SurfaceEmissivity,
}

impl LoweredRadiation {
    pub(super) fn is_pointwise(&self) -> bool {
        self.pointwise
    }

    /// Pointwise emission is inside one Newton solve, without a radiation fixed point.
    pub(super) fn max_solid_solves(&self) -> usize {
        if self.pointwise { 1 } else { self.config.max_iterations }
    }

    pub(super) fn is_enclosure(&self) -> bool {
        self.enclosure.is_some()
    }

    pub(super) fn surface_count(&self) -> usize {
        self.surfaces.len()
    }

    /// Reuse the steady producer's exact face partition and area admission.
    pub(super) fn bind_enclosure(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
    ) -> Result<Option<GrayDiffuseEnclosure>, ConductionError> {
        self.enclosure.as_ref()
            .map(|_| enclosure::bind(cx, problem, self))
            .transpose()
    }

    /// Immutable user-facing identities in the same order as the bound patches.
    pub(super) fn surface_labels(&self) -> impl ExactSizeIterator<Item = (&str, &str)> {
        self.surfaces
            .iter()
            .map(|source| (source.name.as_str(), source.target.as_str()))
    }
}

pub(super) fn pointwise_requested(spec: &ProjectSpec) -> bool {
    use fs_project::RadiationSpatialModel::Pointwise;
    spec.cooling.as_ref().and_then(|cooling| cooling.conduction.as_ref())
        .is_some_and(|setup| {
            setup.radiation.as_ref().is_some_and(|row| row.spatial_model == Some(Pointwise))
                || setup.enclosure_radiation.as_ref()
                    .is_some_and(|row| row.spatial_model == Some(Pointwise))
        })
}

pub(super) fn lower(
    spec: &ProjectSpec,
    setup: &ConductionSetup,
    cards: &CardPackSet,
) -> Result<Option<LoweredRadiation>, SolveRefusal> {
    if pointwise_requested(spec) {
        super::transient::admit(spec)?;
    }
    if setup.enclosure_radiation.is_some() {
        return enclosure::lower(spec, setup, cards).map(Some);
    }
    let Some(declaration) = setup.radiation.as_ref() else {
        return Ok(None);
    };
    let mut patches = Vec::with_capacity(declaration.surfaces.len());
    let mut surfaces = Vec::with_capacity(declaration.surfaces.len());
    for surface in &declaration.surfaces {
        if !setup.boundaries.iter().any(|boundary| {
            boundary.target == surface.target && boundary.condition.supports_ambient_radiation()
        }) {
            return Err(conduction_error(
                "cli-solve-conduction-radiation-boundary",
                format!(
                    "radiating surface `{}` target `{}` has no convection or explicit zero-flux boundary",
                    surface.name, surface.target
                ),
                "declare convection or exactly zero outward heat flux on this exterior target; no fictitious convection is needed for radiation-only exchange",
            ));
        }
        let emissivity = resolve_emissivity(
            &surface.name,
            &surface.card,
            surface.claim.as_deref(),
            surface.query_temperature.value,
            cards,
        )?;
        patches.push(
            AmbientRadiationPatch::new(
                &surface.target,
                emissivity.clone(),
                surface.reservoir_temperature.value,
            )
            .map_err(|error| {
                conduction_error(
                    "cli-solve-conduction-radiation-surface",
                    format!(
                        "surface `{}` radiation declaration refused: {error}",
                        surface.name
                    ),
                    "check the exterior target and positive absolute reservoir temperature",
                )
            })?,
        );
        surfaces.push(SurfaceSource {
            name: surface.name.clone(),
            target: surface.target.clone(),
            emissivity,
        });
    }
    Ok(Some(LoweredRadiation {
        patches,
        config: AmbientRadiationConfig {
            max_iterations: declaration.max_iterations as usize,
            temperature_tolerance_k: declaration.temperature_tolerance.value,
            balance_tolerance_w: declaration.heat_tolerance.value,
            balance_relative_tolerance: 0.0,
            relaxation: declaration.relaxation,
        },
        surfaces,
        enclosure: None,
        pointwise: declaration.spatial_model == Some(fs_project::RadiationSpatialModel::Pointwise),
    }))
}

/// Both native radiation models resolve the same immutable emissivity claim.
fn resolve_emissivity(
    name: &str,
    identity: &str,
    claim: Option<&str>,
    temperature_k: f64,
    cards: &CardPackSet,
) -> Result<SurfaceEmissivity, SolveRefusal> {
    let library = cards.library();
    let card = library.material(identity).ok_or_else(|| {
        conduction_error(
            "cli-solve-conduction-radiation-card",
            format!("surface {name} references absent emissivity card {identity}"),
            "supply the exact immutable material pack containing the surface emissivity claim",
        )
    })?;
    match claim {
        Some(pin) => SurfaceEmissivity::from_card_pinned(
            name, card, temperature_k,
            parse_claim_id(pin, &format!("radiating surface {name}"))?,
        ),
        None => SurfaceEmissivity::from_card(
            name, card, temperature_k, SelectionPolicy::SingleClaimOnly,
        ),
    }.map_err(|error| conduction_error(
        "cli-solve-conduction-radiation-card",
        format!("surface {name} emissivity query refused: {error}"),
        "provide an in-domain dimensionless hemispherical-total-emissivity claim in (0, 1], with an exact claim pin when selection is ambiguous",
    ))
}

/// Last accepted enclosure endpoint for either declared storage coordinate. The model is bound once
/// per spatial rung; only endpoint evidence is retained between time steps.
pub(super) struct EnclosureEndpoint {
    pub(super) combined_boundary: fs_conduction::ThermalBoundary,
    pub(super) convective_robin_fluxes: Vec<RobinFlux>,
    pub(super) convective_out_w: f64,
    pub(super) report: EnthalpyEnclosureReport,
    pub(super) time_integration: &'static str,
    pub(super) dt_s: f64,
    pub(super) physical_residual_norm_j: f64,
    pub(super) physical_residual_tolerance_j: f64,
    pub(super) physical_energy_residual_j: f64,
    pub(super) decomposition_residual_w: f64,
}

pub(super) struct SolidSolution {
    pub(super) conduction: ConductionSolution,
    convective_fluxes: Option<Vec<RobinFlux>>,
    convective_out: Option<f64>,
    radiation: Option<AmbientRadiationReport>,
    enclosure: Option<enclosure::EnclosureResult>,
    pointwise: Option<PointwiseEndpoint>,
    /// The exact combined convective + radiative Robin partition behind a
    /// radiating field; the roundoff bound must see this operator.
    pub(super) combined_boundary: Option<fs_conduction::ThermalBoundary>,
}

impl SolidSolution {
    /// Wrap a conduction producer that carries no radiative split or tangent.
    pub(super) fn from_conduction(conduction: ConductionSolution) -> Self {
        Self {
            conduction,
            convective_fluxes: None,
            convective_out: None,
            radiation: None,
            enclosure: None,
            pointwise: None,
            combined_boundary: None,
        }
    }

    /// Retain the exact accepted radiative endpoint and its original convection
    /// split. A transient producer supplies the same physical boundary evidence
    /// as the steady producer; storage stays in its separate time receipt.
    pub(super) fn from_radiative_endpoint(
        conduction: ConductionSolution,
        combined_boundary: fs_conduction::ThermalBoundary,
        convective_robin_fluxes: Vec<RobinFlux>,
        convective_out_w: f64,
        radiation: AmbientRadiationReport,
    ) -> Self {
        Self {
            conduction,
            convective_fluxes: Some(convective_robin_fluxes),
            convective_out: Some(convective_out_w),
            radiation: Some(radiation),
            enclosure: None,
            pointwise: None,
            combined_boundary: Some(combined_boundary),
        }
    }

    pub(super) fn from_enclosure_endpoint(
        conduction: ConductionSolution,
        model: GrayDiffuseEnclosure,
        endpoint: EnclosureEndpoint,
    ) -> Self {
        let enclosure = enclosure::EnclosureResult::from_transient(model, &endpoint);
        Self {
            conduction,
            convective_fluxes: Some(endpoint.convective_robin_fluxes),
            convective_out: Some(endpoint.convective_out_w),
            radiation: None,
            enclosure: Some(enclosure),
            pointwise: None,
            combined_boundary: Some(endpoint.combined_boundary),
        }
    }

    pub(super) fn from_pointwise_endpoint(
        conduction: ConductionSolution,
        ambient: Option<fs_conduction::PointwiseRadiationReport>,
        enclosure: Option<fs_conduction::PointwiseEnclosureReport>,
    ) -> Result<Self, SolveRefusal> {
        let endpoint = PointwiseEndpoint::new(ambient, enclosure)?;
        // Original convection is unmodified by the pointwise owner. The common
        // energy row includes radiation, so retain its independently summed
        // convection account rather than inferring it from that combined row.
        let convective_out_w = conduction.report.robin_fluxes.iter()
            .try_fold(0.0, |sum, row| {
                let value = sum + row.heat_rate_w;
                if value.is_finite() { Ok(value) } else {
                    Err(conduction_error("cli-solve-conduction-radiation-nonfinite",
                        "pointwise convection sum is nonfinite", "inspect the physical heat inputs"))
                }
            })?;
        Ok(Self {
            conduction,
            convective_fluxes: None,
            convective_out: Some(convective_out_w),
            radiation: None,
            enclosure: None,
            pointwise: Some(endpoint),
            combined_boundary: None,
        })
    }

    pub(super) fn enclosure_model(&self) -> Option<&GrayDiffuseEnclosure> {
        self.enclosure.as_ref().map(enclosure::EnclosureResult::model)
    }

    /// Retain the producer's immutable laws, never infer emissivity from a
    /// rounded receipt or a difference of combined coefficients.
    pub(super) fn radiation_patches(&self) -> Vec<AmbientRadiationPatch> {
        self.radiation.as_ref().map_or_else(Vec::new, |report| {
            report.patches.iter().map(|row| row.patch.clone()).collect()
        })
    }

    pub(super) fn convective_robin_fluxes(&self) -> &[RobinFlux] {
        self.convective_fluxes
            .as_deref()
            .unwrap_or(&self.conduction.report.robin_fluxes)
    }

    pub(super) fn convective_out_w(&self) -> f64 {
        self.convective_out
            .unwrap_or(self.conduction.report.energy.robin_out_w)
    }

    pub(super) fn radiation_receipt(
        &self,
        lowering: Option<&LoweredRadiation>,
    ) -> Result<Option<String>, SolveRefusal> {
        if let Some(endpoint) = &self.pointwise {
            return endpoint.receipt(lowering.ok_or_else(|| conduction_error(
                "cli-solve-conduction-pointwise-radiation",
                "pointwise endpoint lost its declaration", "retain the original radiation model",
            ))?, self.convective_out_w()).map(Some);
        }
        if let Some(report) = &self.enclosure {
            return report
                .receipt(
                    lowering.expect("enclosure solution retains its declaration"),
                    self.convective_out_w(),
                )
                .map(Some);
        }
        let Some(report) = &self.radiation else {
            return Ok(None);
        };
        let lowering = lowering.expect("a radiation solution has a lowering");
        let mut rows = Vec::with_capacity(report.patches.len());
        for row in &report.patches {
            let source = lowering
                .surfaces
                .iter()
                .find(|source| source.target == row.patch.region())
                .expect("the backend retains the declared patches");
            rows.push(format!(
                "{{\"name\":{},\"target\":{},\"card\":{},\"material_state\":{},\"emissivity_receipt\":{},\"emissivity\":{},\"query_temperature_k\":{},\"ambient_temperature_k\":{},\"area_m2\":{},\"mean_temperature_k\":{},\"driving_temperature_k\":{},\"coefficient_w_m2_k\":{},\"applied_heat_w\":{},\"nonlinear_heat_w\":{},\"heat_mismatch_w\":{},\"heat_tolerance_w\":{}}}",
                json_string(&source.name), json_string(&source.target),
                json_string(&source.emissivity.card_identity().to_hex()),
                json_string(source.emissivity.material_state()),
                json_string(&source.emissivity.receipt().content_hash().to_hex()),
                num(source.emissivity.value())?, num(source.emissivity.temperature_k())?,
                num(row.patch.ambient_temperature_k())?, num(row.area_m2)?,
                num(row.mean_surface_temperature_k)?, num(row.driving_temperature_k)?,
                num(row.applied_coefficient_w_m2_k)?, num(row.applied_heat_w)?,
                num(row.nonlinear_heat_w)?, num(row.heat_mismatch_w)?, num(row.heat_tolerance_w)?,
            ));
        }
        let convective = self.convective_out_w();
        let radiative = report.applied_radiation_out_w;
        let denominator = convective.abs() + radiative.abs();
        let share = if denominator > 0.0 {
            num(radiative.abs() / denominator)?
        } else {
            "null".to_string()
        };
        Ok(Some(format!(
            "{{\"model\":\"area-mean-gray-surface-to-black-reservoir\",\"surfaces\":[{}],\"solid_solves\":{},\"solid_iterations\":{},\"krylov_iterations\":{},\"max_temperature_change_k\":{},\"max_heat_mismatch_w\":{},\"convective_out_w\":{},\"radiative_out_w\":{},\"nonlinear_radiative_out_w\":{},\"decomposition_residual_w\":{},\"radiative_magnitude_share\":{},\"controls\":{{\"max_iterations\":{},\"temperature_tolerance_k\":{},\"heat_tolerance_w\":{},\"relative_heat_tolerance\":{},\"relaxation\":{}}},\"authority\":\"Estimated\",\"no_claim\":{}}}",
            rows.join(","),
            report.iterations,
            report.solid_iterations,
            report.krylov_iterations,
            num(report.max_temperature_change_k)?,
            num(report.max_heat_mismatch_w)?,
            num(convective)?,
            num(radiative)?,
            num(report.nonlinear_radiation_out_w)?,
            num(report.decomposition_residual_w)?,
            share,
            report.config.max_iterations,
            num(report.config.temperature_tolerance_k)?,
            num(report.config.balance_tolerance_w)?,
            num(report.config.balance_relative_tolerance)?,
            num(report.config.relaxation)?,
            json_string(RADIATION_NO_CLAIM),
        )))
    }
}

fn num(value: f64) -> Result<String, SolveRefusal> {
    super::canonical_f64(value).ok_or_else(|| {
        conduction_error(
            "cli-solve-conduction-radiation-nonfinite",
            "the radiation report contains a non-finite number",
            "report the solver defect; non-finite radiation evidence is never published",
        )
    })
}

const RADIATION_NO_CLAIM: &str = "card-backed emissivity is fixed at the declared query temperature; each patch radiates epsilon sigma A (area-mean temperature^4 - reservoir temperature^4) to a black isothermal reservoir with view factor one; the pointwise Robin trace uses a converged secant coefficient; no enclosure, shadowing, participating medium, temperature-dependent emissivity, pointwise fourth-power integral, or experimental model validation is claimed; the radiative magnitude share is |radiation|/(|radiation|+|convection|), not a source-power fraction";

/// A paired physical-model comparison, never an uncertainty bound.
#[derive(Debug)]
pub(super) enum RadiationSensitivity {
    Measured { on_k: f64, off_k: f64 },
    Unavailable { reason: String },
}

impl RadiationSensitivity {
    pub(super) fn json(&self) -> Result<String, SolveRefusal> {
        Ok(match self {
            Self::Measured { on_k, off_k } => format!(
                "{{\"state\":\"measured\",\"method\":\"paired-base-fidelity-radiation-on-off\",\"radiation_on_k\":{},\"radiation_off_k\":{},\"delta_on_minus_off_k\":{},\"absolute_difference_k\":{},\"authority\":\"Estimated\",\"no_claim\":{}}}",
                num(*on_k)?,
                num(*off_k)?,
                num(on_k - off_k)?,
                num((on_k - off_k).abs())?,
                json_string(
                    "a sensitivity between two declared physical models on the same base mesh; not a bound on model-form error or a validation of either model"
                ),
            ),
            Self::Unavailable { reason } => format!(
                "{{\"state\":\"no-data\",\"reason\":{}}}",
                json_string(reason),
            ),
        })
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn measure_sensitivity(
    ledger: &fs_ledger::Ledger,
    spec: &ProjectSpec,
    cards: &CardPackSet,
    context: &super::StageContext,
    run: super::SolveRunId,
    work: super::EvidenceWork<'_>,
    resume: bool,
    available_wall_s: f64,
    propagated_nominal: Option<f64>,
) -> Result<Option<RadiationSensitivity>, SolveRefusal> {
    let Some(setup) = spec
        .cooling
        .as_ref()
        .and_then(|cooling| cooling.conduction.as_ref())
    else {
        return Ok(None);
    };
    if setup.radiation.is_none() && setup.enclosure_radiation.is_none() {
        return Ok(None);
    }
    if setup.transient.is_none()
        && setup.boundaries.iter().all(|boundary| {
            matches!(
                boundary.condition,
                ThermalBoundaryCondition::HeatFlux { .. }
            )
        })
    {
        // Removing the only temperature anchor has no uniquely defined
        // steady comparison. Do not spend the nominal solve's budget trying
        // a singular off-model or invent an ambient convection reference.
        return Ok(Some(RadiationSensitivity::Unavailable {
            reason: "radiation-off steady problem has no temperature anchor; no unique thermal equilibrium can be compared".into(),
        }));
    }
    let Some(region) = super::temperature_maximum_region(spec) else {
        return Ok(None);
    };
    let mut base = super::base_fidelity(spec);
    let solve = |project: &ProjectSpec| {
        let product = super::conduction_solve_receipt(
            ledger,
            project,
            cards,
            context,
            run,
            work,
            resume,
            available_wall_s,
            None,
            1.0,
            0.0,
            0.0,
        )?;
        super::region_maximum(&product.qoi_inputs, region, work)
    };
    let on = match propagated_nominal {
        Some(value) => Some(value),
        None => solve(&base)?,
    };
    base.cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap()
        .radiation = None;
    base.cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap()
        .enclosure_radiation = None;
    let comparison = match solve(&base) {
        Ok(Some(off_k)) => match on {
            Some(on_k) => RadiationSensitivity::Measured { on_k, off_k },
            None => RadiationSensitivity::Unavailable {
                reason: "the radiation-on solve has no requested region maximum".to_string(),
            },
        },
        Ok(None) => RadiationSensitivity::Unavailable {
            reason: "the radiation-off solve has no requested region maximum".to_string(),
        },
        Err(error)
            if matches!(
                error.code,
                "cli-solve-cancelled" | "cli-solve-work-envelope"
            ) =>
        {
            return Err(error);
        }
        Err(error) => RadiationSensitivity::Unavailable {
            reason: format!(
                "radiation-off solve refused: {} ({})",
                error.what, error.code
            ),
        },
    };
    Ok(Some(comparison))
}

pub(super) fn solve(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    interfaces: Option<&ThermalInterfaces>,
    config: SolveConfig,
    radiation: Option<&LoweredRadiation>,
) -> Result<SolidSolution, ConductionError> {
    if let Some(radiation) = radiation {
        if radiation.enclosure.is_some() {
            return enclosure::solve(cx, problem, interfaces, config, radiation);
        }
        let solved = fs_conduction::solve_with_ambient_radiation(
            cx,
            problem,
            interfaces,
            &radiation.patches,
            config,
            radiation.config,
        )?;
        Ok(SolidSolution::from_radiative_endpoint(
            solved.conduction,
            solved.combined_boundary,
            solved.convective_robin_fluxes,
            solved.convective_out_w,
            solved.radiation,
        ))
    } else {
        let conduction = match interfaces {
            Some(interfaces) => {
                fs_conduction::solve_with_interfaces(cx, problem, interfaces, config)
            }
            None => fs_conduction::solve(cx, problem, config),
        }?;
        Ok(SolidSolution {
            conduction,
            convective_fluxes: None,
            convective_out: None,
            radiation: None,
            enclosure: None,
            pointwise: None,
            combined_boundary: None,
        })
    }
}
