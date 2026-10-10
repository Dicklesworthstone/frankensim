//! Opt-in continuum enclosure of a native selected-region volume mean.
//!
//! This never reinterprets the decision maximum as a bounded functional. The
//! original final field is passed to the existing equilibrated primal/dual
//! producer. Only a genuinely global scalar source law can supply unbounded
//! coefficient validity; neither flat sampled data nor a small nodal range do.

use super::{
    EvidenceWork, ProjectSpec, RungSolved, SolveRefusal, conduction_error,
    temperature_maximum_region,
};
use crate::cards::CardPackSet;
use fs_exec::Cx;
use std::collections::BTreeMap;

pub(super) const OUTPUT: &str = "temperature-volume-mean-bound";
#[cfg(feature = "thermal-verification")]
const SCOPE: &str = "Verified numerical interval for the region volume-mean temperature of the nominal linear conduction PDE on the exact published polyhedral mesh domain, with original fixed Dirichlet/Neumann/Robin data and globally declared scalar conductivity. Independent equilibrated primal and dual fluxes include discretization, algebraic residual and outward numerical evaluation. Original material cards and exact selected property receipts are retained; their declarations are not experimental validation. This is not a maximum-temperature bound, material or boundary uncertainty bound, CAD/as-built geometry certificate, nonlinear/coupled/contact bound, or a complete engineering decision. The existing temperature-max requirement and its uncertainty budget are unchanged.";

fn bad(message: impl Into<String>) -> SolveRefusal {
    conduction_error(
        "cli-solve-mean-bound",
        message,
        "request a steady no-contact linear model with fixed thermal boundaries and globally declared scalar conductivity; retain bounded material spans when no continuum temperature-range proof exists",
    )
}

fn unavailable() -> SolveRefusal {
    conduction_error(
        "cli-solve-mean-bound-feature",
        "temperature-volume-mean-bound requires the opt-in thermal-verification feature",
        "build frankensim with --features fs-cli/thermal-verification, then rerun the unchanged project",
    )
}

pub(super) fn requested(spec: &ProjectSpec) -> Result<bool, SolveRefusal> {
    let rows: Vec<_> = spec
        .outputs
        .as_deref()
        .unwrap_or(&[])
        .iter()
        .filter(|row| row.name == OUTPUT)
        .collect();
    if rows.is_empty() {
        return Ok(false);
    }
    if rows.len() != 1 || rows[0].kind != "report" || rows[0].region.is_some() {
        return Err(bad(
            "declare exactly one temperature-volume-mean-bound report without :region; it uses the existing temperature-max requirement's volume region",
        ));
    }
    if temperature_maximum_region(spec).is_none() {
        return Err(bad(
            "the region-mean bound requires the existing temperature-max requirement and its declared volume region",
        ));
    }
    if !cfg!(feature = "thermal-verification") {
        return Err(unavailable());
    }
    let setup = spec
        .cooling
        .as_ref()
        .and_then(|c| c.conduction.as_ref())
        .ok_or_else(|| bad("the volume-mean bound requires a conduction model"))?;
    if setup.transient.is_some()
        || setup.radiation.is_some()
        || setup.boundaries.iter().any(|b| {
            matches!(
                b.condition,
                fs_project::ThermalBoundaryCondition::NaturalConvection { .. }
                    | fs_project::ThermalBoundaryCondition::AirflowConvection { .. }
            )
        })
    {
        return Err(bad(
            "the volume-mean bound requires steady fixed thermal boundary laws; transient, radiation and natural/forced-convection feedback need their own continuum bounds",
        ));
    }
    Ok(true)
}

pub(super) fn extract(
    cx: &Cx<'_>,
    spec: &ProjectSpec,
    cards: &CardPackSet,
    solved: &RungSolved,
    ids: &BTreeMap<String, u32>,
    work: EvidenceWork<'_>,
) -> Result<String, SolveRefusal> {
    #[cfg(feature = "thermal-verification")]
    {
        enabled::extract(cx, spec, cards, solved, ids, work)
    }
    #[cfg(not(feature = "thermal-verification"))]
    {
        let _ = (cx, spec, cards, solved, ids, work);
        Err(unavailable())
    }
}

#[cfg(feature = "thermal-verification")]
mod enabled {
    use super::super::{canonical_f64, json_string};
    use super::*;
    use fs_conduction::verification::region::{
        GoalResidualLimits, RegionMeanConfig, bound_temperature_region_mean,
    };
    use fs_conduction::verification::{ConductionBoundError, FluxBudget, TetError};
    use fs_conduction::{
        ConductionProblem, ConductivityModel, ConductivityTable, ElementMaterials, MaterialId,
        MaterialTable,
    };

    fn poll(cx: &Cx<'_>, work: EvidenceWork<'_>) -> Result<(), SolveRefusal> {
        if work.is_requested() || cx.checkpoint().is_err() {
            return Err(conduction_error(
                "cli-solve-cancelled",
                "the requested continuum volume-mean bound was cancelled",
                "resume the completed pipeline prefix; no partial bound was published",
            ));
        }
        Ok(())
    }
    fn number(value: f64) -> Result<String, SolveRefusal> {
        canonical_f64(value).ok_or_else(|| bad("the mean bound contains nonfinite arithmetic"))
    }
    fn interval(lo: f64, hi: f64) -> Result<String, SolveRefusal> {
        if lo > hi {
            return Err(bad("the mean bound interval has reversed endpoints"));
        }
        Ok(format!(
            "{{\"lower\":{},\"upper\":{}}}",
            number(lo)?,
            number(hi)?
        ))
    }
    fn lower(error: fs_conduction::ConductionError) -> SolveRefusal {
        bad(format!("original material/thermal model refused: {error}"))
    }

    pub(super) fn extract(
        cx: &Cx<'_>,
        spec: &ProjectSpec,
        cards: &CardPackSet,
        solved: &RungSolved,
        ids: &BTreeMap<String, u32>,
        work: EvidenceWork<'_>,
    ) -> Result<String, SolveRefusal> {
        poll(cx, work)?;
        if !requested(spec)? {
            return Err(bad("the volume-mean bound was not requested"));
        }
        let data = solved
            .adjoint_data
            .as_ref()
            .ok_or_else(|| bad("final native operator is absent"))?;
        if solved.interface_pair_count != 0
            || data
                .interfaces
                .as_ref()
                .is_some_and(|i| i.surface_count() != 0)
            || !data.air_paths.is_empty()
            || data.natural_goal.is_some()
            || solved.natural_fragment.is_some()
            || solved.conjugate_fragment.is_some()
            || data.radiating_boundary.is_some()
            || !data.radiation_patches.is_empty()
        {
            return Err(bad(
                "this mean-bound product rung requires no contacts and fixed thermal boundaries; a frozen coupled/nonlinear operator cannot certify the original PDE",
            ));
        }
        let n = solved.mesh.vertex_count();
        let ne = solved.mesh.element_count();
        if solved.labels.len() != ne || solved.solution.temperature.len() != n {
            return Err(bad("published mesh, field and region labels disagree"));
        }
        // Conservative additional-work admission, separate from the primal's
        // mesh allowance. The lower producer also enforces its cell/iteration
        // bounds; this does not assert exact allocator peak accounting.
        let memory = spec.budgets.as_ref().map_or(0, |b| b.memory_bytes);
        let needed = ne
            .checked_mul(32_768)
            .and_then(|v| n.checked_mul(8_192).and_then(|w| v.checked_add(w)))
            .and_then(|v| u64::try_from(v).ok());
        if needed.is_none_or(|bytes| bytes > memory / 2) || ne > 65_536 {
            return Err(bad(
                "the declared memory budget does not admit the additional primal/dual flux reconstruction",
            ));
        }
        let region =
            temperature_maximum_region(spec).ok_or_else(|| bad("missing volume region"))?;
        let region_id = *ids
            .get(region)
            .ok_or_else(|| bad("the requested volume region has no mesh label"))?;
        let cells: Vec<_> = solved
            .labels
            .iter()
            .enumerate()
            .filter_map(|(e, &label)| (label == region_id).then_some(e))
            .collect();
        if cells.is_empty() {
            return Err(bad("the requested volume region contains no tetrahedra"));
        }

        let library = cards.library();
        let mut models = Vec::with_capacity(ids.len());
        let mut material_rows = Vec::with_capacity(ids.len());
        for (name, &id) in ids {
            poll(cx, work)?;
            let binding = spec
                .materials
                .as_deref()
                .unwrap_or(&[])
                .iter()
                .find(|binding| binding.region == *name)
                .ok_or_else(|| bad("an original region material binding is absent"))?;
            let card = library
                .material(&binding.card)
                .ok_or_else(|| bad("the exact original material card is absent"))?;
            let native = data.materials.table().get(MaterialId(id)).map_err(lower)?;
            if native.is_temperature_dependent() {
                return Err(bad(
                    "temperature-dependent conductivity requires its own nonlinear continuum bound",
                ));
            }
            let native_receipts = native.receipts();
            let selected = native_receipts
                .first()
                .ok_or_else(|| bad("native conductivity has no source receipt"))?
                .selected;
            if native_receipts.iter().any(|r| r.selected != selected) {
                return Err(bad(
                    "native conductivity changed source claims across its sampled span",
                ));
            }
            let constant = ConductivityTable::from_unbounded_constant_claim(
                card.claims(),
                fs_project::THERMAL_CONDUCTIVITY_PROPERTY,
                binding.temp_lo.value,
                selected,
            )
            .map_err(lower)?;
            let receipt = &constant.receipts()[0];
            let receipt_bytes = receipt
                .to_bytes()
                .map_err(|e| bad(format!("constant source receipt refused: {e}")))?;
            let receipt_hex = receipt_bytes
                .iter()
                .map(|b| format!("{b:02x}"))
                .collect::<String>();
            let conductivity = constant.eval(binding.temp_lo.value).map_err(lower)?;
            material_rows.push(format!(
                concat!(
                    "{{\"region\":{},\"card\":{},\"selected_claim\":{},",
                    "\"receipt\":{},\"receipt_bytes_hex\":{},\"conductivity_w_m_k\":{}}}"
                ),
                json_string(name),
                json_string(&binding.card),
                json_string(&selected.0.to_hex()),
                json_string(&receipt.content_hash().to_hex()),
                json_string(&receipt_hex),
                number(conductivity)?
            ));
            let model = ConductivityModel::isotropic(constant);
            if model.tensor_at(binding.temp_lo.value).map_err(lower)?
                != native.tensor_at(binding.temp_lo.value).map_err(lower)?
            {
                return Err(bad(
                    "the globally constant source law does not equal the original native conductivity",
                ));
            }
            models.push((MaterialId(id), model));
        }
        let materials = ElementMaterials::new(
            MaterialTable::new(models).map_err(lower)?,
            data.materials.of_element().to_vec(),
        )
        .map_err(lower)?;
        for (e, tet) in solved.mesh.complex().tets.iter().enumerate() {
            poll(cx, work)?;
            let temperature = tet
                .iter()
                .map(|&v| solved.solution.temperature[v as usize] / 4.0)
                .sum();
            if materials
                .model_for(e)
                .map_err(lower)?
                .tensor_at(temperature)
                .map_err(lower)?
                != data
                    .materials
                    .model_for(e)
                    .map_err(lower)?
                    .tensor_at(temperature)
                    .map_err(lower)?
            {
                return Err(bad(
                    "the verifier and published native element conductivity differ",
                ));
            }
        }
        let fallback = materials.model_for(0).map_err(lower)?;
        let problem = ConductionProblem {
            mesh: &solved.mesh,
            boundary: &data.boundary,
            material: fallback,
            element_materials: Some(&materials),
            source: &data.source,
        };
        let result = bound_temperature_region_mean(
            cx,
            problem,
            &[],
            &solved.solution.temperature,
            &cells,
            RegionMeanConfig {
                dual: data.linear,
                residual_limits: GoalResidualLimits {
                    max_rows: n,
                    max_nonzeros: usize::try_from(memory / 64).unwrap_or(usize::MAX),
                },
                flux: FluxBudget {
                    max_cells: ne,
                    max_iterations: 128,
                },
            },
        )
        .map_err(|error| match error {
            ConductionBoundError::Verification(TetError::Cancelled)
            | ConductionBoundError::Conduction(fs_conduction::ConductionError::Cancelled {
                ..
            }) => conduction_error(
                "cli-solve-cancelled",
                "continuum mean verification cancelled",
                "resume the accepted pipeline prefix; no partial bound was published",
            ),
            error => bad(format!("continuum region-mean producer refused: {error}")),
        })?;
        poll(cx, work)?;
        let bound = &result.bound;
        let difference = bound.enclosure.sub(bound.candidate_mean);
        let error_upper = difference.lo.abs().max(difference.hi.abs());
        let candidate = bound.candidate_mean.lo * 0.5 + bound.candidate_mean.hi * 0.5;
        Ok(format!(
            concat!(
                "{{\"schema\":\"frankensim.cli.volume-mean-bound.v1\",",
                "\"name\":{},\"functional\":\"region-volume-mean-temperature\",\"region\":{},\"unit\":\"K\",",
                "\"value_k\":{},\"candidate_mean_k\":{},\"enclosure_k\":{},\"error_upper_k\":{},",
                "\"authority\":\"Verified\",\"producer\":\"fs-conduction::verification::region\",",
                "\"selected_cells\":{},\"domain_cells\":{},\"region_volume_m3\":{},",
                "\"primal_energy_error_upper\":{},\"dual_energy_error_upper\":{},",
                "\"residual_correction_integral\":{},\"remainder_integral_upper\":{},",
                "\"dual_iterations\":{},\"materials\":[{}],\"scope\":{}}}"
            ),
            json_string(OUTPUT),
            json_string(region),
            number(candidate)?,
            interval(bound.candidate_mean.lo, bound.candidate_mean.hi)?,
            interval(bound.enclosure.lo, bound.enclosure.hi)?,
            number(error_upper)?,
            cells.len(),
            ne,
            interval(bound.region_volume.lo, bound.region_volume.hi)?,
            number(bound.integral.primal.energy_error_upper)?,
            number(bound.integral.dual.energy_error_upper)?,
            interval(
                bound.integral.residual_correction.lo,
                bound.integral.residual_correction.hi
            )?,
            number(bound.integral.remainder_upper)?,
            result.dual_analysis.dual_iterations,
            material_rows.join(","),
            json_string(SCOPE)
        ))
    }
}
