//! Exact held-reference pullback of the original objective strain derivative.
//!
//! The native correction stays in elastic C, preserving anisotropic material
//! axes. Both endpoints share one G: W0 = Jg We, P0 = Jg Pe G^-T.
//! Consequently P0:Fdelta = Jg Pe:Fedelta and P0 Fmid^T = Jg Pe Femid^T.
//! No new discrete-gradient kernel or midpoint inverse is introduced.

use super::{
    DiscreteGradientPoint, DiscreteGradientSettings, EnergyPart, PairResponse, Potential, error,
};
use crate::MaterialError;
use crate::hyper::Hyperelastic;
use crate::hyper::reference::ReferenceDeformation;
use fs_math::dd::Dd;

#[derive(Clone, Copy, Debug)]
pub(super) struct ReferenceBinding<'a> {
    pub(super) reference: &'a ReferenceDeformation,
    pub(super) accepted_total_f: [f64; 9],
}

fn elastic_allowance(original: f64, reference_volume: f64) -> Result<f64, MaterialError> {
    let elastic = original / reference_volume;
    if !original.is_finite()
        || original < 0.0
        || !elastic.is_finite()
        || (original > 0.0 && elastic == 0.0)
    {
        return Err(error(
            "original work allowance cannot be represented per relaxed volume",
        ));
    }
    Ok(elastic)
}

impl<'a> DiscreteGradientPoint<'a> {
    /// Original stored-energy change between two admitted material points.
    ///
    /// Both points must select the same energy part. A change of reference G
    /// is explicit and includes its relaxed-volume factor. Equal native laws
    /// retain their original finite-pair potential increment independently of
    /// absolute storage. Different laws retain their own admitted endpoint
    /// potentials; that fallback cannot recover a parameter change lost in
    /// either scalar endpoint energy.
    ///
    /// At fixed total F, the result is signed external work from the changed
    /// reference or material source. It is not a strain work, discrete stress,
    /// dissipation law or authorization to change an accepted history.
    pub fn energy_increment_to(&self, after: &Self) -> Result<f64, MaterialError> {
        if self.part != after.part {
            return Err(error(
                "stored-energy increment requires the same energy part",
            ));
        }
        let elastic_increment = if core::ptr::eq(self.law, after.law) || self.law == after.law {
            Potential::Material(self.law, self.part).increment(&self.accepted_f, &after.accepted_f)
        } else {
            (Dd::from_f64(after.energy_before_pa) - Dd::from_f64(self.energy_before_pa)).to_f64()
        };
        let identity = ReferenceDeformation::identity();
        let before_reference = self.reference.map_or(&identity, |r| r.reference);
        let after_reference = after.reference.map_or(&identity, |r| r.reference);
        let increment = before_reference.original_density_increment_to(
            after_reference,
            self.energy_before_pa,
            after.energy_before_pa,
            elastic_increment,
        );
        if !increment.is_finite() {
            return Err(error("stored-energy increment is unrepresentable"));
        }
        Ok(increment)
    }

    /// Bind an admitted elastic law to one held reference deformation G.
    ///
    /// Accepted/trial deformations, energies, stresses, tangents and volume
    /// work use the ORIGINAL material reference. Constitutive stretch limits
    /// and the small-increment criterion apply to Fe = F G^-1 and elastic C.
    /// Identity is the original constructor, including floating-point order.
    ///
    /// G must remain fixed over this pair of endpoints. A changing prescribed
    /// reference requires its own signed external-work account. This elastic
    /// owner does not integrate viscous history or its incremental potential.
    pub fn try_new_with_reference(
        law: &'a Hyperelastic,
        part: EnergyPart,
        reference: &'a ReferenceDeformation,
        accepted_total_f: [f64; 9],
        settings: DiscreteGradientSettings,
    ) -> Result<Self, MaterialError> {
        if *reference == ReferenceDeformation::identity() {
            return Self::try_new(law, part, accepted_total_f, settings);
        }
        elastic_allowance(
            settings.absolute_work_tolerance_pa,
            reference.volume_ratio(),
        )?;
        let mut point = Self::try_new(
            law,
            part,
            reference.elastic_deformation(&accepted_total_f),
            settings,
        )?;
        point.reference = Some(ReferenceBinding {
            reference,
            accepted_total_f,
        });
        let volume = point.original_density(point.volume_ratio_before);
        if !point.energy_before_pa().is_finite() || !volume.is_finite() || volume <= 0.0 {
            return Err(error(
                "accepted original-reference storage or volume is unrepresentable",
            ));
        }
        Ok(point)
    }

    pub(super) fn original_density(&self, elastic: f64) -> f64 {
        self.reference
            .map_or(elastic, |r| r.reference.volume_ratio() * elastic)
    }

    pub(super) fn evaluate(
        &self,
        next_f: &[f64; 9],
        potential: Potential<'_>,
        before: f64,
        absolute_tolerance: f64,
        with_tangent: bool,
    ) -> Result<PairResponse, MaterialError> {
        let Some(binding) = self.reference else {
            return self.evaluate_elastic(
                next_f,
                potential,
                before,
                absolute_tolerance,
                with_tangent,
            );
        };
        let reference = binding.reference;
        let jg = reference.volume_ratio();
        let elastic = self.evaluate_elastic(
            &reference.elastic_deformation(next_f),
            potential,
            before,
            elastic_allowance(absolute_tolerance, jg)?,
            with_tangent,
        )?;
        let before = jg * before;
        let after = jg * elastic.value_after;
        let increment = jg * elastic.increment;
        let gradient = reference.pullback_piola(&elastic.gradient);
        let tangent = if with_tangent {
            reference.pullback_tangent(&elastic.tangent)
        } else {
            [[0.0; 9]; 9]
        };
        // Gate the actual original-volume contraction, including transformation
        // roundoff. An elastic-volume work check alone cannot attest this value.
        let work = (0..9).fold(0.0, |sum, i| {
            gradient[i].mul_add(next_f[i] - binding.accepted_total_f[i], sum)
        });
        let defect = work - increment;
        let allowance = absolute_tolerance
            + self.settings.relative_work_tolerance * before.abs().max(after.abs()).max(work.abs());
        if ![before, after, work, defect, allowance]
            .iter()
            .all(|v| v.is_finite())
            || gradient
                .iter()
                .chain(tangent.iter().flatten())
                .any(|v| !v.is_finite())
            || (matches!(potential, Potential::Volume) && after <= 0.0)
        {
            return Err(error(
                "original-reference stress, tangent, volume or work is unrepresentable",
            ));
        }
        if defect.abs() > allowance {
            return Err(error(&format!(
                "original-reference work defect {defect:e} exceeds {allowance:e}"
            )));
        }
        Ok(PairResponse {
            value_after: after,
            increment,
            gradient,
            tangent,
            work,
            defect,
            limit: elastic.limit,
        })
    }
}
