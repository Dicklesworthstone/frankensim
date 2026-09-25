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
use fs_matdb::SelectionPolicy;
use fs_project::ProjectSpec;
use fs_project::spec::{ConductionSetup, ThermalBoundaryCondition};

use super::{SOLVER_FIDELITY_ADAPTIVE, SolveRefusal, conduction_error, parse_claim_id};
use crate::cards::CardPackSet;
use crate::import::json_string;

pub(super) struct LoweredRadiation {
    patches: Vec<AmbientRadiationPatch>,
    config: AmbientRadiationConfig,
    surfaces: Vec<SurfaceSource>,
}

struct SurfaceSource {
    name: String,
    target: String,
    emissivity: SurfaceEmissivity,
}

pub(super) fn lower(
    spec: &ProjectSpec,
    setup: &ConductionSetup,
    cards: &CardPackSet,
) -> Result<Option<LoweredRadiation>, SolveRefusal> {
    let Some(declaration) = setup.radiation.as_ref() else {
        return Ok(None);
    };
    if spec
        .solver
        .as_ref()
        .is_some_and(|solver| solver.fidelity == SOLVER_FIDELITY_ADAPTIVE)
    {
        return Err(conduction_error(
            "cli-solve-conduction-radiation-adaptive",
            "adaptive goals do not yet include the total tangent of the area-mean radiation law",
            "use base or ladder fidelity for declared radiation; an adjoint of frozen Robin coefficients cannot measure this model's goal error",
        ));
    }
    let library = cards.library();
    let mut patches = Vec::with_capacity(declaration.surfaces.len());
    let mut surfaces = Vec::with_capacity(declaration.surfaces.len());
    for surface in &declaration.surfaces {
        if !setup.boundaries.iter().any(|boundary| {
            boundary.target == surface.target
                && matches!(
                    boundary.condition,
                    ThermalBoundaryCondition::Convection { .. }
                        | ThermalBoundaryCondition::AirflowConvection { .. }
                )
        }) {
            return Err(conduction_error(
                "cli-solve-conduction-radiation-boundary",
                format!(
                    "radiating surface `{}` target `{}` has no conventional convection boundary",
                    surface.name, surface.target
                ),
                "declare convection or airflow-convection on this exterior target; this radiation model augments that boundary",
            ));
        }
        let card = library.material(&surface.card).ok_or_else(|| {
            conduction_error(
                "cli-solve-conduction-radiation-card",
                format!(
                    "surface `{}` references absent emissivity card {}",
                    surface.name, surface.card
                ),
                "supply the exact immutable material pack containing the surface emissivity claim",
            )
        })?;
        let emissivity = match surface.claim.as_deref() {
            Some(pin) => SurfaceEmissivity::from_card_pinned(
                &surface.name, card, surface.query_temperature.value,
                parse_claim_id(pin, &format!("radiating surface `{}`", surface.name))?,
            ),
            None => SurfaceEmissivity::from_card(
                &surface.name, card, surface.query_temperature.value, SelectionPolicy::SingleClaimOnly,
            ),
        }.map_err(|error| conduction_error(
            "cli-solve-conduction-radiation-card",
            format!("surface `{}` emissivity query refused: {error}", surface.name),
            "provide an in-domain dimensionless hemispherical-total-emissivity claim in (0, 1], with an exact claim pin when selection is ambiguous",
        ))?;
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
    }))
}

pub(super) struct SolidSolution {
    pub(super) conduction: ConductionSolution,
    convective_fluxes: Option<Vec<RobinFlux>>,
    convective_out: Option<f64>,
    radiation: Option<AmbientRadiationReport>,
}

impl SolidSolution {
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

const RADIATION_NO_CLAIM: &str = "card-backed emissivity is fixed at the declared query temperature; each patch radiates epsilon sigma A (area-mean temperature^4 - reservoir temperature^4) to a black isothermal reservoir with view factor one; the pointwise Robin trace uses a converged secant coefficient; no enclosure, shadowing, participating medium, temperature-dependent emissivity, pointwise fourth-power integral, radiative adjoint, or experimental model validation is claimed; the radiative magnitude share is |radiation|/(|radiation|+|convection|), not a source-power fraction";

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
    if !spec
        .cooling
        .as_ref()
        .and_then(|cooling| cooling.conduction.as_ref())
        .is_some_and(|setup| setup.radiation.is_some())
    {
        return Ok(None);
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
        let solved = fs_conduction::solve_with_ambient_radiation(
            cx,
            problem,
            interfaces,
            &radiation.patches,
            config,
            radiation.config,
        )?;
        Ok(SolidSolution {
            conduction: solved.conduction,
            convective_fluxes: Some(solved.convective_robin_fluxes),
            convective_out: Some(solved.convective_out_w),
            radiation: Some(solved.radiation),
        })
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
        })
    }
}
