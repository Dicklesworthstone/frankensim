//! Preserve original bending and torsional Prony laws in the same bank owner.
//!
//! Source histories precede regional body histories and explicit layered
//! additions. All use the original exact mechanical and material-memory owner.
use super::*;
use fs_scenario::acoustic::PronyBending;
use fs_material::section::torsion::ConcentricPronyTorsionalSection;

impl MaterialStringSource {
    pub(in crate::bowed_string::runtime::reciprocal::bank) fn torsional_relaxation(&self)
        -> Option<&ConcentricPronyTorsionalSection>
    {
        self.torsion_source().and_then(|twist| twist.relaxation())
    }

    pub(in crate::bowed_string::runtime::reciprocal::bank) fn bending_relaxation(&self)
        -> Option<PronyBending>
    {
        match self {
            Self::WithPolarizations(p) => p.source().bending_relaxation(),
            Self::WithTorsion(t) => t.bending().bending_relaxation(),
            Self::Homogeneous(s) => s.string().relaxation_bending,
            Self::Wound(s) => s.string().relaxation_bending,
            Self::Concentric { .. } | Self::ConcentricModalLoss { .. } | Self::Measured(_) => None,
            Self::ConcentricProny { section, .. } =>
                Some(crate::material_assembly::concentric::concentric_prony_bending(section)),
        }
    }

    // The homogeneous adapter admits eight arms; a source-bound section
    // already admits at most 128 across its physical layers. The complete
    // bank's 256-state and setup budgets still bound the projected histories.
    fn bending_relaxation_arm_limit(&self) -> usize {
        match self {
            Self::WithPolarizations(p) => p.source().bending_relaxation_arm_limit(),
            Self::WithTorsion(t) => t.bending().bending_relaxation_arm_limit(),
            Self::ConcentricProny { .. } =>
                fs_material::section::relaxation::MAX_SECTION_RELAXATION_BRANCHES,
            _ => 8,
        }
    }

    // Other cards/impact consumers retain `lower` and its memory refusal.
    pub(in crate::bowed_string::runtime::reciprocal::bank) fn lower_for_bank(
        &self, cx: &Cx<'_>, gas: &GasState, rate: u32,
    ) -> Result<BowedStringCard, ReciprocalStringError> {
        match self {
            Self::WithPolarizations(p) => p.source().lower_for_bank(cx, gas, rate),
            Self::WithTorsion(t) => t.bending().lower_for_bank(cx, gas, rate),
            Self::Homogeneous(s) if s.string().relaxation_bending.is_some() =>
                BowedStringCard::material_relaxation_backbone_with_gas(cx, s, gas, rate)
                    .map_err(ReciprocalStringError::Material),
            Self::Wound(s) if s.string().relaxation_bending.is_some() =>
                BowedStringCard::wound_relaxation_backbone_with_gas(cx, s, gas, rate)
                    .map_err(ReciprocalStringError::Material),
            Self::ConcentricProny { template, section, prestress } =>
                BowedStringCard::concentric_relaxation_backbone_with_gas(
                    cx, template.clone(), section, *prestress, gas, rate)
                    .map_err(ReciprocalStringError::Material),
            _ => self.lower(cx, gas, rate),
        }
    }
}

impl MaterialStringBankSource {
    pub(in crate::bowed_string::runtime::reciprocal::bank) fn validate_source_memory_frequencies(
        &self, system: &FrictionModalSystem, gate: &CancelGate,
    ) -> Result<(), ReciprocalStringError> {
        let high = system.linear_frequency_upper_bound_rad_s();
        // Bending and torsion histories are original first-order states in the
        // exact coupled owner, including when nonlinear friction is attached.
        // Their source support is independent of the midpoint dt/tau cap.
        for member in &self.strings {
            poll(gate)?;
            if let Some(law) = member.source.bending_relaxation() {
                let (lo, hi) = law.omega_band_rad_s;
                if !lo.is_finite() || !hi.is_finite() || lo != 0.0 || hi < high
                    || law.branches.iter().any(|b| !b.relaxing_stiffness_n_m2.is_finite()
                        || b.relaxing_stiffness_n_m2 < 0.0
                        || !b.relaxation_time_s.is_finite() || b.relaxation_time_s <= 0.0) {
                    return Err(invalid("string material memory requires positive finite times and source support from DC through the complete instantaneous coupled band"));
                }
            }
            if let Some(law) = member.source.torsional_relaxation() {
                let (lo, hi) = law.omega_band_rad_s();
                if lo != 0.0 || hi < high {
                    return Err(invalid("torsional material memory requires source support from DC through the complete instantaneous coupled band"));
                }
            }
        }
        Ok(())
    }

    pub(in crate::bowed_string::runtime::reciprocal::bank) fn source_memory_count(
        &self, cards: &[BowedStringCard], gate: &CancelGate,
    ) -> Result<usize, ReciprocalStringError> {
        if cards.len() != self.strings.len() {
            return Err(invalid("source memory requires its original complete string layout"));
        }
        let mut count = 0usize;
        for (member, card) in self.strings.iter().zip(cards) {
            poll(gate)?;
            if let Some(law) = member.source.bending_relaxation() {
                if law.branches.len() > member.source.bending_relaxation_arm_limit() {
                    return Err(invalid("string memory exceeds the bounded source-arm count"));
                }
                let positive = law.branches.iter().filter(|b| b.relaxing_stiffness_n_m2 > 0.0).count();
                let planes = if matches!(&member.source, MaterialStringSource::WithPolarizations(_)) { 2 } else { 1 };
                count = positive.checked_mul(planes).and_then(|n| n.checked_mul(card.mode_count))
                    .and_then(|n| n.checked_add(count))
                    .ok_or_else(|| invalid("complete string material history extent overflow"))?;
            }
            if let Some(twist) = member.source.torsion_source() {
                if let Some(law) = twist.relaxation() {
                    let positive = law.branches().iter()
                        .filter(|arm| arm.relaxing_stiffness_n_m2 > 0.0).count();
                    count = positive.checked_mul(twist.mode_count())
                        .and_then(|n| n.checked_add(count))
                        .ok_or_else(|| invalid("complete torsional material history extent overflow"))?;
                }
            }
        }
        if count > 256 {
            return Err(invalid("source string memory exceeds the original state ceiling"));
        }
        let substrate = self.material_shell_substrate_relaxation()
            .map_or(0, |image| image.relaxation_times_s().len());
        count = count.checked_add(substrate)
            .ok_or_else(|| invalid("complete source substrate history extent overflow"))?;
        Ok(count)
    }

    fn source_memory_branches(
        &self, cards: &[BowedStringCard], system: &FrictionModalSystem, gate: &CancelGate,
    ) -> Result<Vec<RelaxationBranch>, ReciprocalStringError> {
        if cards.len() > system.components().len() {
            return Err(invalid("source memory lost its original string components"));
        }
        let dimension = system.mode_count().checked_mul(2)
            .ok_or_else(|| invalid("string memory state extent overflow"))?;
        let count = self.source_memory_count(cards, gate)?;
        if count.checked_mul(dimension).is_none_or(|n| n > self.options.network.max_setup_terms) {
            return Err(invalid("source string memory exceeds the original setup allowance"));
        }
        let mut branches = Vec::with_capacity(count);
        let mut start = 0;
        for (index, (member, card)) in self.strings.iter().zip(cards).enumerate() {
            poll(gate)?;
            let component = &system.components()[index];
            if component.modes().len() != member.source.modes() {
                return Err(invalid("string memory lost its original mechanical coordinate ordering"));
            }
            if let Some(law) = member.source.bending_relaxation() {
                let planes = if matches!(&member.source, MaterialStringSource::WithPolarizations(_)) { 2 } else { 1 };
                for arm in law.branches.iter().filter(|b| b.relaxing_stiffness_n_m2 > 0.0) {
                    for plane in 0..planes {
                        for mode in 0..card.mode_count {
                            poll(gate)?;
                            let wave_number = (mode + 1) as f64 * core::f64::consts::PI / card.length_m;
                            let stiffness = arm.relaxing_stiffness_n_m2 * wave_number.powi(4)
                                / card.linear_density_kg_m;
                            if !stiffness.is_finite() || stiffness <= 0.0 {
                                return Err(invalid("source string bending increment is unrepresentable in the actual mass basis"));
                            }
                            let mut projection = vec![0.0; dimension];
                            projection[2*(start + plane*card.mode_count + mode)] = 1.0;
                            branches.push(RelaxationBranch { projection, stiffness,
                                relaxation_time_s: arm.relaxation_time_s });
                        }
                    }
                }
            }
            start += component.modes().len();
        }
        // Keep every original bending history first, then append source twist
        // histories in string/layer/arm/mode order. Body and explicit bending
        // additions remain after this source block in the same accepted state.
        start = 0;
        for (index, (member, card)) in self.strings.iter().zip(cards).enumerate() {
            poll(gate)?;
            let component = &system.components()[index];
            if let Some(twist) = member.source.torsion_source() {
                if let Some(law) = twist.relaxation() {
                    let offset = component.modes().len().checked_sub(twist.mode_count())
                        .ok_or_else(|| invalid("torsional memory lost its original twist coordinates"))?;
                    let inertia = law.torsional().rotary_inertia_kg_m();
                    for arm in law.branches().iter().filter(|arm| arm.relaxing_stiffness_n_m2 > 0.0) {
                        for mode in 0..twist.mode_count() {
                            poll(gate)?;
                            let wave_number = (mode + 1) as f64 * core::f64::consts::PI / card.length_m;
                            // Integral of twist-gradient storage in the actual
                            // rotary mass basis: delta_GJ * k^2 / sum(rho J).
                            let stiffness = arm.relaxing_stiffness_n_m2
                                * wave_number.powi(2) / inertia;
                            if !stiffness.is_finite() || stiffness <= 0.0 {
                                return Err(invalid("torsional memory stiffness is unrepresentable in its rotary mass basis"));
                            }
                            let mut projection = vec![0.0; dimension];
                            projection[2 * (start + offset + mode)] = 1.0;
                            branches.push(RelaxationBranch { projection, stiffness,
                                relaxation_time_s: arm.relaxation_time_s });
                        }
                    }
                }
            }
            start += component.modes().len();
        }
        // Explicit substrate memory is a source history, before optional
        // coating/body additions. The same ordering survives cold rebuilds.
        branches.extend(self.shell_substrate_memory_branches(system, gate)?);
        Ok(branches)
    }

    // Rebuild one combined owner after a cold topology or gas-inertia change.
    // Accepted energy-normalized strains z do not rescale with q's mass basis.
    // Explicit layered selections retain their existing stale-basis refusal.
    pub(in crate::bowed_string::runtime::reciprocal::bank) fn restore_material_memory(
        &self, mut system: FrictionModalSystem, cards: &[BowedStringCard],
        memory: Option<Vec<f64>>, gate: &CancelGate,
    ) -> Result<FrictionModalSystem, ReciprocalStringError> {
        poll(gate)?;
        if self.material_chart_relaxation.is_none() && self.string_bending_relaxation.is_none()
            && self.longitudinal_waves.as_ref().and_then(|axial| axial.material_relaxation()).is_none()
            && self.material_chart_tensor_relaxation.is_none()
            && self.material_shell_coating_relaxation.is_none()
            && self.material_shell_substrate_relaxation().is_none()
            && self.spatial_cavity.as_ref().and_then(|cavity|cavity.thermal_loss()).is_none()
            && self.strings.iter().all(|m| m.source.bending_relaxation().is_none()
                && m.source.torsional_relaxation().is_none()) {
            if memory.as_ref().is_some_and(|z| !z.is_empty()) {
                return Err(invalid("material history has no original source binding"));
            }
            return self.restore_body_radiation(system, gate);
        }
        let mut branches = self.source_memory_branches(cards, &system, gate)?;
        branches.extend(self.body_memory_branches(&system, gate)?);
        if let Some(image) = &self.string_bending_relaxation {
            branches.extend(image.branches(&system, gate)?);
        }
        if let Some(cavity) = &self.spatial_cavity {
            branches.extend(cavity.thermal_memory_branches(self, &system, gate)?);
        }
        if let Some(image) = self.longitudinal_waves.as_ref()
            .and_then(|axial| axial.material_relaxation()) {
            branches.extend(image.branches(&system, gate)?);
        }
        system = system.with_material_relaxation(branches, memory, gate)
            .map_err(ReciprocalStringError::Friction)?;
        system=self.restore_body_tensor_constraints(system,cards,gate)?;
        self.validate_body_memory_frequencies(&system, gate)?;
        self.validate_material_bridge_frequencies(&system, gate)?;
        poll(gate)?;
        self.restore_body_radiation(system, gate)
    }
}

