//! Exact original-source admission and explicit material-point mapping.

use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_matdb::PropertyKey;
use fs_material::hardness::plastic_history::ACCUMULATED_EQUIVALENT_PLASTIC_STRAIN_AXIS as ALPHA;
use fs_material::state_point::{
    POISSON_RATIO_PROPERTY, ResolvedInterfaceStatePoint, YOUNG_MODULUS_PROPERTY,
};
use fs_qty::{Dims, QuantitySpec};

use super::super::{validate_envelope, validate_fixed_normal};
use super::{
    BoundNormalInterfaceState, ContactHardnessField, FiniteGapPlasticModel,
    HistoryHardnessContactError, PlasticContactEnvelope, checkpoint, invalid,
};

pub(super) fn admit_sources(
    normal: &BoundNormalInterfaceState,
    coefficients: &ResolvedInterfaceStatePoint,
    keys: &[PropertyKey; 2],
    hardness: &[ContactHardnessField; 2],
    envelope: PlasticContactEnvelope,
    maximum_recession_m: [f64; 2],
) -> Result<([f64; 2], [usize; 2]), HistoryHardnessContactError> {
    validate_fixed_normal(normal)?;
    validate_envelope(envelope)?;
    let interface = normal.interface_state();
    if coefficients.card_identity() != interface.card_identity()
        || coefficients.surface_materials() != interface.surface_materials()
        || coefficients.query_point() != interface.query_point()
        || coefficients.axis_quantities() != interface.axis_quantities()
    {
        return Err(invalid(
            "wear coefficients must retain the original ordered interface and fixed query",
        ));
    }
    let surfaces = [
        normal.elastic().surface_a_state(),
        normal.elastic().surface_b_state(),
    ];
    let mut k = [0.0; 2];
    let mut selected = [0; 2];
    for side in 0..2 {
        let surface = surfaces[side];
        let bulk = surface.resolved();
        let source = hardness[side].hardness_source();
        let h = source.resolved();
        if hardness[side].is_empty()
            || h.card_identity() != bulk.card_identity()
            || h.material() != bulk.material()
            || bulk
                .properties()
                .iter()
                .any(|p| !p.is_static_at_fixed_external_conditions())
        {
            return Err(invalid(
                "mapped hardness requires nonempty original ordered material sources and fixed elasticity",
            ));
        }
        match &hardness[side] {
            ContactHardnessField::Fixed(_) => {
                if h.query_point() != bulk.query_point()
                    || h.axis_quantities() != bulk.axis_quantities()
                    || !source.property().is_static_at_fixed_external_conditions()
                {
                    return Err(invalid(
                        "fixed hardness must retain the exact original surface query and test",
                    ));
                }
            }
            ContactHardnessField::Accepted(field) => {
                let mechanical = field.mechanical_source();
                if mechanical.card_identity() != bulk.card_identity()
                    || mechanical.material() != bulk.material()
                    || mechanical.query_point() != bulk.query_point()
                    || mechanical.axis_quantities() != bulk.axis_quantities()
                    || !h
                        .query_point()
                        .iter()
                        .filter(|(name, _)| name.as_str() != ALPHA)
                        .eq(bulk.query_point().iter())
                    || !h
                        .axis_quantities()
                        .iter()
                        .filter(|(name, _)| name.as_str() != ALPHA)
                        .eq(bulk.axis_quantities().iter())
                    || h.axis_quantities().get(ALPHA)
                        != Some(&QuantitySpec::dimensional(Dims::NONE))
                    || !source
                        .property()
                        .has_only_owned_and_fixed_external_conditions(
                            ALPHA,
                            QuantitySpec::dimensional(Dims::NONE),
                        )
                {
                    return Err(invalid(
                        "accepted history may add only its original dimensionless alpha to the fixed surface environment",
                    ));
                }
                // The native J2/cyclic source owns canonical E and nu. Equal
                // numerical moduli from unrelated source roles are insufficient.
                // A different elastic parameterization needs an explicit adapter.
                if surface.elastic_state().derivation().is_some() {
                    return Err(invalid(
                        "native-history contact requires original canonical Young/Poisson elastic roles",
                    ));
                }
                for (name, value) in [
                    (YOUNG_MODULUS_PROPERTY, surface.young_modulus_pa()),
                    (POISSON_RATIO_PROPERTY, surface.poisson_ratio()),
                ] {
                    let selected = mechanical
                        .property(name)
                        .ok_or_else(|| invalid("missing original native elastic source role"))?;
                    if Some(selected) != bulk.property(name)
                        || !selected.is_static_at_fixed_external_conditions()
                        || selected.value_si().to_bits() != value.to_bits()
                    {
                        return Err(invalid(
                            "native-history and half-space elasticity must retain identical original selected claims",
                        ));
                    }
                }
            }
        }
        let coefficient = coefficients
            .property_by_key(&keys[side])
            .ok_or_else(|| invalid("missing exact original wear coefficient"))?;
        if keys[side].quantity() != QuantitySpec::dimensional(Dims::NONE)
            || !coefficient.is_static_at_fixed_external_conditions()
            || !coefficient.value_si().is_finite()
            || coefficient.value_si() < 0.0
            || !maximum_recession_m[side].is_finite()
            || maximum_recession_m[side] < 0.0
            || maximum_recession_m[side] >= envelope.half_space_depth_m[side]
        {
            return Err(invalid(
                "fixed nonnegative k and finite recession support within each homogeneous material column required",
            ));
        }
        k[side] = coefficient.value_si();
        selected[side] = coefficients
            .properties()
            .iter()
            .position(|p| p == coefficient)
            .ok_or_else(|| invalid("wear coefficient must be retained by its original bundle"))?;
    }
    Ok((k, selected))
}

pub(super) fn mapped_caps(
    fields: &[ContactHardnessField; 2],
    mapping: &[[usize; 2]],
    cx: &Cx<'_>,
) -> Result<Vec<f64>, HistoryHardnessContactError> {
    let mut caps = Vec::with_capacity(mapping.len());
    for (cell, points) in mapping.iter().enumerate() {
        if cell % 256 == 0 {
            checkpoint(cx)?;
        }
        let a = fields[0]
            .hardness_pa(points[0])
            .ok_or_else(|| invalid("body A material-point map index outside original field"))?;
        let b = fields[1]
            .hardness_pa(points[1])
            .ok_or_else(|| invalid("body B material-point map index outside original field"))?;
        caps.push(a.min(b));
    }
    checkpoint(cx)?;
    Ok(caps)
}

pub(super) fn admit_successor(
    old: &[ContactHardnessField; 2],
    next: &[ContactHardnessField; 2],
    cx: &Cx<'_>,
) -> Result<(), HistoryHardnessContactError> {
    for side in 0..2 {
        if old[side].source_identity() != next[side].source_identity()
            || old[side].len() != next[side].len()
        {
            return Err(invalid(
                "history refresh must retain the original source, test, native field count and order",
            ));
        }
        match (&old[side], &next[side]) {
            (ContactHardnessField::Fixed(_), ContactHardnessField::Fixed(_)) => {}
            (ContactHardnessField::Accepted(before), ContactHardnessField::Accepted(after)) => {
                for (index, (a, b)) in before.samples().iter().zip(after.samples()).enumerate() {
                    if index % 256 == 0 {
                        checkpoint(cx)?;
                    }
                    if b.alpha() < a.alpha() {
                        return Err(invalid(
                            "accepted accumulated plastic strain cannot decrease on field refresh",
                        ));
                    }
                }
            }
            _ => {
                return Err(invalid(
                    "field refresh cannot replace a fixed source by a native-history source",
                ));
            }
        }
    }
    checkpoint(cx)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn identity(
    normal: &BoundNormalInterfaceState,
    coefficients: &ResolvedInterfaceStatePoint,
    selected: &[usize; 2],
    fields: &[ContactHardnessField; 2],
    mapping: &[[usize; 2]],
    model: &FiniteGapPlasticModel,
    envelope: PlasticContactEnvelope,
    maximum_recession_m: [f64; 2],
    cx: &Cx<'_>,
) -> Result<ContentHash, HistoryHardnessContactError> {
    let mut hash = DomainHasher::new("org.frankensim.fs-contact.source-history-hardness-wear.v1");
    hash.update(normal.identity().as_bytes());
    hash.update(coefficients.identity().as_bytes());
    hash.update(model.identity().as_bytes());
    for side in 0..2 {
        hash.update(&(selected[side] as u64).to_le_bytes());
        hash.update(fields[side].source_identity().as_bytes());
        hash.update(&(fields[side].len() as u64).to_le_bytes());
        match &fields[side] {
            ContactHardnessField::Fixed(source) => {
                hash.update(&[0]);
                hash.update(source.identity().as_bytes());
            }
            ContactHardnessField::Accepted(field) => {
                hash.update(&[1]);
                hash.update(field.identity().as_bytes());
            }
        }
    }
    hash.update(&(mapping.len() as u64).to_le_bytes());
    for (cell, points) in mapping.iter().enumerate() {
        if cell % 256 == 0 {
            checkpoint(cx)?;
        }
        for point in points {
            hash.update(&(*point as u64).to_le_bytes());
        }
    }
    for value in [
        envelope.maximum_approach_m,
        envelope.maximum_plastic_depth_m,
        envelope.half_space_depth_m[0],
        envelope.half_space_depth_m[1],
        envelope.maximum_patch_to_depth[0],
        envelope.maximum_patch_to_depth[1],
        maximum_recession_m[0],
        maximum_recession_m[1],
    ] {
        hash.update(&value.to_bits().to_le_bytes());
    }
    checkpoint(cx)?;
    Ok(hash.finalize())
}
