//! Opt-in continuum enclosure of a native selected-region volume mean.
//!
//! This never reinterprets the decision maximum as a bounded functional. The
//! original final field is passed to the existing equilibrated primal/dual
//! producer. Only a genuinely global scalar source law can supply unbounded
//! coefficient validity; neither flat sampled data nor a small nodal range do.
//! Finite scalar-constant cards require a separate whole-continuum slab proof.

use super::{
    EvidenceWork, ProjectSpec, RungSolved, SolveRefusal, conduction_error,
    temperature_maximum_region,
};
use crate::cards::CardPackSet;
use fs_exec::Cx;
use std::collections::BTreeMap;

pub(super) const OUTPUT: &str = "temperature-volume-mean-bound";
#[cfg(feature = "thermal-verification")]
mod material_range;
#[cfg(feature = "thermal-verification")]
const SCOPE: &str = "Verified numerical interval for the region volume-mean temperature of the nominal linear conduction PDE on the exact published polyhedral mesh domain, with original fixed Dirichlet/Neumann/Robin data and constant scalar conductivity. Material validity is either globally declared or independently established by the retained whole-continuum slab comparison; layered contacts include their original finite-resistance jumps in that proof. Original matching finite-resistance contacts retain separate temperature traces and constant scalar resistance whose source domain is global or established by the same whole-continuum comparison; their resistance-weighted jumps enter the continuum majorant. Independent equilibrated primal and dual fluxes include discretization, algebraic residual and outward numerical evaluation. Original material and interface cards and exact selected property receipts are retained; their declarations are not experimental validation. The primary enclosure bounds the mean, not the error of the selected maximum. This is not a material/contact/boundary uncertainty bound, CAD/as-built geometry certificate, nonlinear/nonmatching-contact/convection-feedback bound, or a complete engineering decision. The existing temperature-max requirement and its uncertainty budget are unchanged.";

fn bad(message: impl Into<String>) -> SolveRefusal {
    conduction_error(
        "cli-solve-mean-bound",
        message,
        "request steady fixed thermal boundaries and constant scalar conductivity; finite spans require an admitted slab or layered-contact continuum range with all original scalar source domains retained",
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
        || setup.enclosure_radiation.is_some()
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
        GoalResidualLimits, RegionMeanConfig,
    };
    use fs_conduction::verification::{ConductionBoundError, FluxBudget, TetError};
    use fs_conduction::{
        ConductionProblem, ElementMaterials, MaterialId,
        MaterialTable,
    };

    fn poll(cx: &Cx<'_>, work: EvidenceWork<'_>) -> Result<(), SolveRefusal> {
        if work.is_requested() || cx.checkpoint().is_err() {
            return Err(conduction_error(
                "cli-solve-cancelled",
                "the requested continuum volume-mean bound was cancelled",
                "resume the accepted pipeline prefix; no partial bound was published",
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
    fn span_json(span: fs_conduction::TemperatureSpan) -> Result<String, SolveRefusal> {
        match span {
            fs_conduction::TemperatureSpan::Unbounded => Ok("null".into()),
            fs_conduction::TemperatureSpan::Sampled { low, high } => interval(low, high),
        }
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
        if !data.air_paths.is_empty()
            || data.natural_goal.is_some()
            || solved.natural_fragment.is_some()
            || solved.conjugate_fragment.is_some()
            || data.radiating_boundary.is_some()
            || !data.radiation_patches.is_empty()
        {
            return Err(bad(
                "this mean-bound product rung requires fixed thermal boundaries; a frozen coupled/nonlinear operator cannot certify the original PDE",
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
        let pair_count = data
            .interface_surfaces
            .iter()
            .try_fold(0usize, |count, surface| {
                count.checked_add(surface.face_pairs().len())
            });
        if pair_count != Some(solved.interface_pair_count)
            || data.interfaces.as_ref().map_or(0, |i| i.surface_count())
                != data.interface_surfaces.len()
        {
            return Err(bad(
                "the original contact declarations and published operator disagree",
            ));
        }
        let mut contact_rows = Vec::with_capacity(data.interface_surfaces.len());
        let mut contact_span = fs_conduction::TemperatureSpan::Unbounded;
        for surface in &data.interface_surfaces {
            poll(cx, work)?;
            let resistance = surface.resistance().ok_or_else(|| {
                bad("mapped contact resistance requires its own source-validity admission")
            })?;
            let binding = spec
                .interface_cards
                .as_deref()
                .unwrap_or(&[])
                .iter()
                .find(|binding| binding.interface == surface.name())
                .ok_or_else(|| bad("an original contact card binding is absent"))?;
            let card = library
                .interface(&binding.card)
                .ok_or_else(|| bad("the exact original contact card is absent"))?;
            if card.content_hash() != resistance.card_identity()
                || resistance.value_origin() != &fs_conduction::ResistanceValueOrigin::CardClaim
                || data
                    .interfaces
                    .as_ref()
                    .is_none_or(|i| i.surface_is_nonmatching(surface.name()))
            {
                return Err(bad(
                    "contact verification requires the original matching card-backed resistance",
                ));
            }
            let receipt = resistance.receipt();
            let claim = card
                .claims()
                .claim(receipt.selected)
                .ok_or_else(|| bad("the selected original contact claim is absent"))?;
            // Retain the selected source's actual domain. A finite domain
            // requires a whole-continuum proof below, never nodal extrema.
            let domain = contact_domain(claim, resistance.value_m2_k_per_w())?;
            contact_span = contact_span.intersect(domain);
            card.claims().verify_receipt(receipt).map_err(|error| {
                bad(format!("original contact source receipt refused: {error}"))
            })?;
            let receipt_bytes = receipt
                .to_bytes()
                .map_err(|error| bad(format!("contact source receipt refused: {error}")))?;
            let receipt_hex = receipt_bytes
                .iter()
                .map(|b| format!("{b:02x}"))
                .collect::<String>();
            let pairs = surface
                .face_pairs()
                .iter()
                .map(|pair| format!("[{},{}]", pair.side_a, pair.side_b))
                .collect::<Vec<_>>()
                .join(",");
            contact_rows.push(format!(
                concat!(
                    "{{\"interface\":{},\"card\":{},\"selected_claim\":{},",
                    "\"receipt\":{},\"receipt_bytes_hex\":{},\"resistance_m2_k_per_w\":{},",
                    "\"face_pairs\":[{}],\"temperature_domain_k\":{}}}"
                ),
                json_string(surface.name()),
                json_string(&binding.card),
                json_string(&receipt.selected.0.to_hex()),
                json_string(&receipt.content_hash().to_hex()),
                json_string(&receipt_hex),
                number(resistance.value_m2_k_per_w())?,
                pairs, span_json(domain)?,
            ));
        }
        let mut models = Vec::with_capacity(ids.len());
        let mut finite_material = false;
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
            let (model, receipt) = material_range::constant_model(
                card.claims(), native, binding.temp_lo.value, selected,
            )?;
            let receipt_bytes = receipt
                .to_bytes()
                .map_err(|e| bad(format!("constant source receipt refused: {e}")))?;
            let receipt_hex = receipt_bytes
                .iter()
                .map(|b| format!("{b:02x}"))
                .collect::<String>();
            let conductivity = model.tensor_at(binding.temp_lo.value).map_err(lower)?[0][0];
            finite_material |= matches!(model.temperature_span(), fs_conduction::TemperatureSpan::Sampled { .. });
            let temperature_domain = match model.temperature_span() {
                fs_conduction::TemperatureSpan::Unbounded => "null".to_string(),
                fs_conduction::TemperatureSpan::Sampled { low, high } => interval(low, high)?,
            };
            material_rows.push(format!(
                concat!(
                    "{{\"region\":{},\"card\":{},\"selected_claim\":{},",
                    "\"receipt\":{},\"receipt_bytes_hex\":{},\"conductivity_w_m_k\":{},\"temperature_domain_k\":{}}}"
                ),
                json_string(name),
                json_string(&binding.card),
                json_string(&selected.0.to_hex()),
                json_string(&receipt.content_hash().to_hex()),
                json_string(&receipt_hex),
                number(conductivity)?, temperature_domain
            ));
            if model.tensor_at(binding.temp_lo.value).map_err(lower)?
                != native.tensor_at(binding.temp_lo.value).map_err(lower)?
            {
                return Err(bad(
                    "the retained constant source law does not equal the original native conductivity",
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
        let config = RegionMeanConfig {
            dual: data.linear,
            residual_limits: GoalResidualLimits {
                max_rows: n,
                max_nonzeros: usize::try_from(memory / 64).unwrap_or(usize::MAX),
            },
            flux: FluxBudget { max_cells: ne, max_iterations: 128 },
        };
        let needs_contact_proof = !data.interface_surfaces.is_empty()
            && (finite_material || matches!(contact_span, fs_conduction::TemperatureSpan::Sampled { .. }));
        let (result, material_domain, contact_proof) = if needs_contact_proof {
            fs_conduction::verification::region::material_range::
                bound_temperature_region_mean_with_contact_slab_validity(
                    cx, problem, &data.interface_surfaces, &solved.solution.temperature, &cells,
                    &[(0,false),(0,true),(1,false),(1,true),(2,false),(2,true)], contact_span, config,
                ).map(|checked| {
                    let validity = (checked.continuum_temperature_range_k, checked.material_temperature_span,
                        checked.contact_temperature_span, checked.axis, checked.reverse, checked.comparison_attempts);
                    (checked.region, None, Some(validity))
                })
        } else {
            material_range::bound(cx, problem, &data.interface_surfaces,
                &solved.solution.temperature, &cells, config)
                .map(|(result, proof)| (result, proof, None))
        }
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
        let material_validity = if let Some((range, materials, contacts, axis, reverse, attempts)) = contact_proof {
            format!(
                "{{\"method\":\"continuum-layered-contact-comparison\",\"continuum_temperature_range_k\":{},\"retained_material_domain_k\":{},\"retained_contact_domain_k\":{},\"axis\":{},\"reverse\":{},\"comparison_attempts\":{}}}",
                interval(range[0], range[1])?, span_json(materials)?, span_json(contacts)?,
                axis, reverse, attempts,
            )
        } else { match material_domain {
            None => "{\"method\":\"global-constant-source\"}".to_string(),
            Some(proof) => format!(
                "{{\"method\":\"continuum-slab-comparison\",\"continuum_temperature_range_k\":{},\"retained_material_domain_k\":{},\"axis\":{},\"reverse\":{},\"two_sided\":{},\"comparison_attempts\":{}}}",
                interval(proof.range[0], proof.range[1])?, interval(proof.span[0], proof.span[1])?,
                proof.axis, proof.reverse, proof.two_sided, proof.attempts,
            ),
        }};
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
                "\"dual_iterations\":{},\"materials\":[{}],\"contacts\":[{}],\"material_validity\":{},\"scope\":{}}}"
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
            contact_rows.join(","),
            material_validity,
            json_string(SCOPE)
        ))
    }
}

#[cfg(feature = "thermal-verification")]
fn contact_domain(
    claim: &fs_matdb::PropertyClaim, actual_resistance: f64,
) -> Result<fs_conduction::TemperatureSpan, SolveRefusal> {
    use fs_conduction::{TemperatureSpan, AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS};
    let fs_matdb::PropertyValue::Scalar { value, dims } = &claim.value else {
        return Err(bad("contact bounds require the original scalar resistance claim"));
    };
    if *dims != AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS || *value != actual_resistance
        || !value.is_finite() || *value <= 0.0
        || claim.interpolation != fs_matdb::InterpolationPolicy::ConstantWithinValidity
        || claim.validity.has_typed_axes()
    {
        return Err(bad("contact bound does not match an original constant scalar resistance law"));
    }
    if claim.validity.bounds().is_empty() { return Ok(TemperatureSpan::Unbounded); }
    let Some((low, high)) = claim.validity.bound("T") else {
        return Err(bad("finite contact validity requires a temperature-only source domain"));
    };
    if claim.validity.bounds().len() != 1 || !low.is_finite() || !high.is_finite() || low > high {
        return Err(bad("finite contact validity requires one finite nonempty temperature domain"));
    }
    Ok(TemperatureSpan::Sampled { low, high })
}

#[cfg(all(test, feature = "thermal-verification"))]
mod contact_domain_tests {
    use super::*;
    use fs_evidence::ValidityDomain;
    use fs_matdb::{InterpolationPolicy, PropertyClaim, PropertyKey, PropertyValue, Provenance, UncertaintyModel};

    fn claim(validity: ValidityDomain) -> PropertyClaim {
        let dims = fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS;
        PropertyClaim {
            key: PropertyKey::new(fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY, dims),
            value: PropertyValue::Scalar { value: 0.2, dims }, validity,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            uncertainty: UncertaintyModel::Unstated, observations: Vec::new(),
            provenance: Provenance { source: "synthetic contact domain fixture".into(),
                license: "internal-test-use".into(), artifact: None },
        }
    }
    #[test]
    fn exact_contact_source_keeps_its_finite_domain_or_global_declaration() {
        let source = claim(ValidityDomain::unconstrained().with("T", 250.0, 350.0));
        assert_eq!(contact_domain(&source, 0.2).unwrap(),
            fs_conduction::TemperatureSpan::Sampled { low: 250.0, high: 350.0 });
        assert_eq!(source.validity.bound("T"), Some((250.0,350.0)));
        assert_eq!(contact_domain(&claim(ValidityDomain::unconstrained()), 0.2).unwrap(),
            fs_conduction::TemperatureSpan::Unbounded);
        assert!(contact_domain(&source, 0.3).is_err(), "a changed numerical resistance is not the selected claim");
    }
    #[test]
    fn other_state_axes_flat_curves_and_wrong_dimensions_do_not_become_contact_proofs() {
        for domain in [ValidityDomain::unconstrained().with("pressure", 0.0, 1.0),
            ValidityDomain::unconstrained().with("T",250.0,350.0).with("pressure",0.0,1.0)] {
            assert!(contact_domain(&claim(domain),0.2).is_err());
        }
        let mut source = claim(ValidityDomain::unconstrained().with("T",250.0,350.0));
        source.value = PropertyValue::Curve { abscissa:"T".into(),
            abscissa_dims:fs_qty::Dims([0,0,0,1,0,0]), knots:vec![(250.0,0.2),(350.0,0.2)],
            dims:fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS };
        assert!(contact_domain(&source,0.2).is_err());
        source.value = PropertyValue::Scalar { value:0.2, dims:fs_qty::Dims::NONE };
        assert!(contact_domain(&source,0.2).is_err());
    }
}
