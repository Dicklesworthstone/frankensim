//! Keep source-card validity and independently proved continuum validity distinct.
use super::{SolveRefusal, bad};
use fs_conduction::{ConductionProblem, ConductivityModel, ConductivityTable,
    InterfaceSurface, TemperatureSpan};
use fs_conduction::verification::{ConductionBoundError, TetError};
use fs_conduction::verification::region::{RegionMeanConfig, RegionMeanFieldBound,
    bound_temperature_region_mean};
use fs_conduction::verification::region::material_range::
    bound_temperature_region_mean_with_slab_directions;
use fs_exec::Cx;
use fs_matdb::{ClaimId, ClaimSet, InterpolationPolicy, PropertyUsageReceipt, PropertyValue};

pub(super) fn constant_model(
    claims: &ClaimSet, native: &ConductivityModel, query: f64, selected: ClaimId,
) -> Result<(ConductivityModel, PropertyUsageReceipt), SolveRefusal> {
    let claim = claims.claim(selected).ok_or_else(|| bad("selected original material claim is absent"))?;
    if claim.validity.bounds().is_empty() && !claim.validity.has_typed_axes() {
        // Preserve the original global-source path, including its exact receipt.
        let constant = ConductivityTable::from_unbounded_constant_claim(claims,
            fs_project::THERMAL_CONDUCTIVITY_PROPERTY, query, selected)
            .map_err(|e| bad(format!("global constant source refused: {e}")))?;
        let receipt = constant.receipts()[0].clone();
        return Ok((ConductivityModel::isotropic(constant), receipt));
    }
    let PropertyValue::Scalar { value, dims } = &claim.value else {
        return Err(bad("finite spans require an original scalar constant claim, not flat sampled data"));
    };
    if *dims != fs_conduction::CONDUCTIVITY_DIMS
        || claim.interpolation != InterpolationPolicy::ConstantWithinValidity
        || claim.validity.has_typed_axes() || claim.validity.bounds().len() != 1
        || !(query.is_finite() && query > 0.0) || native.is_temperature_dependent()
    {
        return Err(bad("finite spans require constant scalar conductivity with only a temperature validity axis"));
    }
    let (source_low, source_high) = claim.validity.bound("T")
        .ok_or_else(|| bad("finite source validity is not a temperature domain"))?;
    let TemperatureSpan::Sampled { low, high } = native.temperature_span() else {
        return Err(bad("finite source validity was not retained in the original native material"));
    };
    if ![low, high, source_low, source_high].iter().all(|x| x.is_finite())
        || low > high || low < source_low || high > source_high
    {
        return Err(bad("native temperature span exceeds its original material claim validity"));
    }
    let tensor = native.tensor_at(query).map_err(|e| bad(format!("native material refused: {e}")))?;
    if tensor != [[*value,0.0,0.0],[0.0,*value,0.0],[0.0,0.0,*value]] {
        return Err(bad("finite constant source does not equal the original native conductivity"));
    }
    let receipts = native.receipts();
    for receipt in &receipts {
        if receipt.selected != selected { return Err(bad("native source changed across its retained span")); }
        claims.verify_receipt(receipt).map_err(|e| bad(format!("original material receipt refused: {e}")))?;
    }
    let receipt = receipts.first().ok_or_else(|| bad("finite constant material has no original receipt"))?;
    Ok((native.clone(), (**receipt).clone()))
}

pub(super) struct MaterialDomainProof {
    pub(super) range: [f64; 2],
    pub(super) span: [f64; 2],
    pub(super) axis: usize,
    pub(super) reverse: bool,
    pub(super) attempts: usize,
}

pub(super) fn bound(
    cx: &Cx<'_>, problem: ConductionProblem<'_>, contacts: &[InterfaceSurface],
    temperature: &[f64], cells: &[usize], config: RegionMeanConfig,
) -> Result<(RegionMeanFieldBound, Option<MaterialDomainProof>), ConductionBoundError> {
    let mut finite = false;
    for e in 0..problem.mesh.element_count() {
        cx.checkpoint().map_err(|_| ConductionBoundError::Verification(TetError::Cancelled))?;
        let model = match problem.element_materials {
            Some(materials) => materials.model_for(e)?, None => problem.material,
        };
        finite |= matches!(model.temperature_span(), TemperatureSpan::Sampled { .. });
    }
    if !finite {
        return bound_temperature_region_mean(cx, problem, contacts, temperature, cells, config)
            .map(|bound| (bound, None));
    }
    if !contacts.is_empty() {
        return Err(ConductionBoundError::Verification(TetError::Unsupported(
            "finite spans with contact need a contact-aware continuum temperature proof")));
    }
    // At most six deterministic geometry/barrier checks, then ONE dual/majorant.
    // Numerical, budget and cancellation errors never retry as another axis.
    let checked = bound_temperature_region_mean_with_slab_directions(cx, problem,
        temperature, cells, &[(0,false),(0,true),(1,false),(1,true),(2,false),(2,true)], config)?;
    let TemperatureSpan::Sampled { low, high } = checked.material_temperature_span else {
        return Err(ConductionBoundError::Verification(TetError::Invalid("missing original finite material span")));
    };
    Ok((checked.region, Some(MaterialDomainProof {
        range: checked.continuum_temperature_range_k, span: [low, high],
        axis: checked.axis, reverse: checked.reverse, attempts: checked.comparison_attempts,
    })))
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_evidence::ValidityDomain;
    use fs_matdb::{PropertyClaim, PropertyKey, Provenance, UncertaintyModel};

    fn source(value: PropertyValue, validity: ValidityDomain, interpolation: InterpolationPolicy)
        -> (ClaimSet, ClaimId)
    {
        let mut claims=ClaimSet::new();
        let id=claims.insert_claim(PropertyClaim {
            key:PropertyKey::new(fs_project::THERMAL_CONDUCTIVITY_PROPERTY,fs_conduction::CONDUCTIVITY_DIMS),
            value,validity,interpolation,uncertainty:UncertaintyModel::Unstated,
            observations:Vec::new(),provenance:Provenance {
                source:"synthetic finite-domain claim".into(),license:"internal-test-use".into(),artifact:None,
            },
        }).unwrap();
        (claims,id)
    }
    fn scalar()->PropertyValue {
        PropertyValue::Scalar {value:10.0,dims:fs_conduction::CONDUCTIVITY_DIMS}
    }
    #[test]
    fn original_finite_receipt_and_narrower_declared_span_are_preserved() {
        let (claims,id)=source(scalar(),ValidityDomain::unconstrained().with("T",250.0,400.0),
            InterpolationPolicy::ConstantWithinValidity);
        let native=ConductivityModel::isotropic(ConductivityTable::from_claims_pinned(
            &claims,fs_project::THERMAL_CONDUCTIVITY_PROPERTY,&[290.0,310.0],id).unwrap());
        let original=native.clone();
        let (retained,receipt)=constant_model(&claims,&native,290.0,id).unwrap();
        assert_eq!(retained,native); assert_eq!(native,original);
        assert_eq!(retained.temperature_span(),TemperatureSpan::Sampled {low:290.0,high:310.0});
        assert_eq!(&receipt,native.receipts()[0]);
        claims.verify_receipt(&receipt).unwrap();
        assert!(retained.tensor_at(320.0).is_err(),"the card's wider domain does not widen the authored solve span");
    }
    #[test]
    fn flat_curves_other_axes_and_an_unbounded_native_copy_do_not_establish_validity() {
        let native=ConductivityModel::isotropic_declared(10.0).unwrap();
        for (value,domain,interpolation) in [
            (scalar(),ValidityDomain::unconstrained().with("T",250.0,400.0),InterpolationPolicy::ConstantWithinValidity),
            (scalar(),ValidityDomain::unconstrained().with("pressure",0.0,1.0),InterpolationPolicy::ConstantWithinValidity),
            (PropertyValue::Curve {abscissa:"T".into(),abscissa_dims:fs_qty::Dims([0,0,0,1,0,0]),
                knots:vec![(250.0,10.0),(400.0,10.0)],dims:fs_conduction::CONDUCTIVITY_DIMS},
                ValidityDomain::unconstrained().with("T",250.0,400.0),InterpolationPolicy::LinearInside),
        ] {
            let (claims,id)=source(value,domain,interpolation);
            assert!(constant_model(&claims,&native,300.0,id).is_err());
        }
    }
}
