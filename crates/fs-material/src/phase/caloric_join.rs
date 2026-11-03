//! Join original sensible-heat spans without erasing their constitutive source.
//!
//! The common temperature is an exact caloric node, not an interpolated
//! temperature or a new material state. Only the enthalpy datum changes.
use super::*;

impl ResolvedHeatCapacityCurve {
    /// Join adjacent, continuously sourced spans of the SAME single phase.
    ///
    /// `self` ends at the reference temperature; `upper` starts there. The
    /// returned chart assigns exactly `reference_specific_enthalpy_j_kg` to
    /// that node. Both input charts keep their original values and identities.
    /// Their old enthalpy datums may differ.
    ///
    /// The seam must have identical temperature, density, phase and original
    /// Cp/density usage receipts. This rejects different source claims, keys,
    /// query contexts and material revisions even when numerical values agree.
    /// This is NOT a phase transition, a material blend or a source-switching
    /// policy. It does not bridge a temperature gap or extrapolate a source.
    ///
    /// Every native enthalpy knot and the affine sensible-density law survive.
    /// Original receipts remain in increasing temperature order, Cp then
    /// density, with the common pair retained once. The returned plain phase
    /// curve therefore remains usable by source-checking thermal consumers.
    ///
    /// Joining preserves each span's interpolation approximation; it does not
    /// improve its error allowance or identify experimental material behavior.
    /// Rebasing adds floating-point roundoff, not latent heat or physical work.
    /// Unrepresentable arithmetic, collapsed enthalpy intervals and a combined
    /// chart exceeding the existing source-knot budget refuse atomically.
    pub fn try_join_at_reference(
        &self,
        upper: &Self,
        reference_specific_enthalpy_j_kg: f64,
    ) -> Result<Self, PhaseStateError> {
        let invalid = |what| PhaseStateError::InvalidCurve { what };
        let lower = &self.curve;
        let upper = &upper.curve;
        if !reference_specific_enthalpy_j_kg.is_finite()
            || lower.material_card_identity != upper.material_card_identity
            || lower.density_interpolation != PhaseDensityInterpolation::LinearSensibleDensity
            || upper.density_interpolation != PhaseDensityInterpolation::LinearSensibleDensity
        {
            return Err(invalid("caloric join needs a finite datum and one original sensible material"));
        }
        let low_seam = lower.knots.last().ok_or_else(|| invalid("empty lower caloric span"))?;
        let high_seam = upper.knots.first().ok_or_else(|| invalid("empty upper caloric span"))?;
        if low_seam.temperature_k.to_bits() != high_seam.temperature_k.to_bits()
            || low_seam.bulk_density_kg_m3.to_bits() != high_seam.bulk_density_kg_m3.to_bits()
            || low_seam.liquid_mass_fraction.to_bits() != high_seam.liquid_mass_fraction.to_bits()
        {
            return Err(invalid("caloric spans must share their exact temperature, density and phase"));
        }
        let phase = match low_seam.liquid_mass_fraction {
            0.0 => SolidLiquidPhase::Solid,
            1.0 => SolidLiquidPhase::Liquid,
            _ => return Err(invalid("caloric joining cannot introduce a coexistence interval")),
        };
        let a = &lower.source_receipts;
        let b = &upper.source_receipts;
        if a.len() < 4 || b.len() < 4 || a.len() % 2 != 0 || b.len() % 2 != 0
            || a[a.len() - 2..] != b[..2]
        {
            return Err(invalid("caloric spans need identical original Cp/density source usages at the seam"));
        }
        let count = lower.knots.len().checked_add(upper.knots.len() - 1)
            .filter(|&count| count <= MAX_HEAT_CAPACITY_ENTHALPY_KNOTS)
            .ok_or_else(|| invalid("joined caloric chart exceeds the source-knot budget"))?;
        let receipt_count = a.len().checked_add(b.len() - 2)
            .ok_or_else(|| invalid("joined caloric source count overflows"))?;
        let mut knots = Vec::with_capacity(count);
        let rebase = |knot: &EnthalpyPhaseKnot, datum: f64| {
            let delta = knot.specific_enthalpy_j_kg - datum;
            let enthalpy = reference_specific_enthalpy_j_kg + delta;
            if !delta.is_finite() || !enthalpy.is_finite() {
                return Err(invalid("joined caloric enthalpy is not representable"));
            }
            Ok(EnthalpyPhaseKnot { specific_enthalpy_j_kg: enthalpy, ..*knot })
        };
        for knot in &lower.knots[..lower.knots.len() - 1] {
            knots.push(rebase(knot, low_seam.specific_enthalpy_j_kg)?);
        }
        // Assign the requested datum directly, including its signed-zero bits.
        knots.push(EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: reference_specific_enthalpy_j_kg,
            ..*low_seam
        });
        for knot in &upper.knots[1..] {
            knots.push(rebase(knot, high_seam.specific_enthalpy_j_kg)?);
        }
        // Native admission detects collapsed/reversed intervals after rebasing.
        let mut curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(
            lower.material_card_identity, phase, knots,
        )?;
        curve.retain_sensible_density_law();
        let mut receipts = Vec::with_capacity(receipt_count);
        receipts.extend_from_slice(a);
        receipts.extend_from_slice(&b[2..]);
        // Same source-bound identity construction as try_new, over the actual
        // output knots/density law and the actual retained query usages.
        let mut identity = DomainHasher::new(
            "org.frankensim.fs-material.sourced-heat-capacity-curve.v1",
        );
        identity.update(curve.identity().as_bytes());
        identity.update(&(receipts.len() as u64).to_le_bytes());
        for receipt in &receipts {
            identity.update(receipt.content_hash().as_bytes());
        }
        curve.identity = identity.finalize();
        curve.source_receipts = receipts;
        Ok(Self { curve })
    }
}

#[cfg(test)]
mod tests;
