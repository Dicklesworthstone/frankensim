//! Exact restart data, not an archived solver cache or an authenticated history.
//! Recovery re-solves the original baseline and accepted field and compares all
//! stored numerical evidence bitwise before publishing an optimizer. Recovery
//! work has its own caller-owned allowance; study work is never refunded.

use super::*;

const MAGIC: &[u8] = b"fs-topols/projected-multiload/checkpoint/1\n";
const MAX_BYTES: usize = 4 * 1024 * 1024;

fn word(out: &mut Vec<u8>, value: u64) { out.extend_from_slice(&value.to_le_bytes()); }
fn number(out: &mut Vec<u8>, value: f64) { word(out, value.to_bits()); }
fn count(out: &mut Vec<u8>, value: usize) { word(out, value as u64); }
fn edge_tag(edge: DesignBoxEdge) -> u64 {
    match edge {
        DesignBoxEdge::Left => 0, DesignBoxEdge::Right => 1,
        DesignBoxEdge::Bottom => 2, DesignBoxEdge::Top => 3,
    }
}

fn state_bytes(out: &mut Vec<u8>, state: &MultiLoadProjectedState) {
    number(out, state.objective);
    number(out, state.volume);
    word(out, state.snapshot);
    word(out, state.active_case.map_or(u64::MAX, |index| index as u64));
    for &value in &state.case_compliances { number(out, value); }
}

fn stress_bytes(out: &mut Vec<u8>, stress: Option<&RobustSampledStressEvaluation>) {
    let Some(stress) = stress else { word(out, 0); return; };
    word(out, 1);
    for &value in &stress.case_compliances { number(out, value); }
    for &value in &stress.case_sampled_max_von_mises { number(out, value); }
    for &point in &stress.case_max_locations {
        number(out, point[0]); number(out, point[1]);
    }
    for &value in &stress.case_sample_counts { count(out, value); }
    number(out, stress.weighted_sum_compliance);
    number(out, stress.worst_weighted_compliance);
    number(out, stress.objective);
    number(out, stress.worst_sampled_von_mises);
    count(out, stress.worst_stress_case);
    number(out, stress.volume);
    word(out, stress.snapshot);
}

struct Reader<'a> { bytes: &'a [u8], at: usize }
impl Reader<'_> {
    fn word(&mut self) -> Result<u64, CutFemError> {
        let raw = self.bytes.get(self.at..self.at + 8)
            .ok_or_else(|| invalid("truncated projected checkpoint"))?;
        let mut array = [0; 8];
        array.copy_from_slice(raw);
        self.at += 8;
        Ok(u64::from_le_bytes(array))
    }
    fn number(&mut self) -> Result<f64, CutFemError> { Ok(f64::from_bits(self.word()?)) }
    fn count(&mut self, cap: usize) -> Result<usize, CutFemError> {
        let value = usize::try_from(self.word()?)
            .map_err(|_| invalid("checkpoint count exceeds platform address space"))?;
        if value > cap { return Err(invalid("projected checkpoint count exceeds its bound")); }
        Ok(value)
    }
    fn field(&mut self, n: usize) -> Result<GridSdf, CutFemError> {
        let nodes = (n + 1) * (n + 1);
        if self.bytes.len() - self.at < nodes * 8 {
            return Err(invalid("truncated projected checkpoint field"));
        }
        let mut phi = GridSdf::from_fn(n, &|_, _| 0.0);
        for value in phi.nodes_mut() {
            *value = self.number()?;
            if !value.is_finite() { return Err(invalid("non-finite checkpoint node")); }
        }
        Ok(phi)
    }
}

fn same_field(left: &GridSdf, right: &GridSdf) -> bool {
    left.n() == right.n() && left.nodes().iter().zip(right.nodes())
        .all(|(a, b)| a.to_bits() == b.to_bits())
}

impl MultiLoadProjectedOptimizer {
    /// Exact original feasible baseline geometry, preserved across recovery.
    #[must_use]
    pub fn baseline_geometry(&self) -> &GridSdf { &self.baseline_geometry }

    /// Immutable declaration, including the ORIGINAL total accepted-update goal.
    #[must_use]
    pub fn settings(&self) -> OptimizeSettings { self.kernel.settings }

    /// Immutable authored load scenarios, in their original order.
    #[must_use]
    pub fn load_cases(&self) -> &[RobustLoadCase] { &self.kernel.load_cases }

    /// Immutable objective aggregation policy.
    #[must_use]
    pub fn aggregate(&self) -> RobustAggregate { self.kernel.aggregate }

    /// Immutable candidate and total study-work policy.
    #[must_use]
    pub fn controls(&self) -> MultiLoadProjectedSettings { self.controls }

    /// Prescribed row-major nodal values in canonical index order.
    #[must_use]
    pub fn fixed_nodes(&self) -> &[(usize, f64)] { &self.fixed }

    /// Immutable numerical area policy.
    #[must_use]
    pub fn projection_settings(&self) -> VolumeProjectionSettings { self.projection }

    /// Deterministic restart bytes for the last accepted state. Includes exact
    /// controls, loads, fixed nodes, area/stress policy, baseline and current
    /// fields/evidence, global ordinal, AL search state and spent study work.
    /// Does not include transient candidate buffers or claim execution history
    /// authenticity. Transport must provide byte integrity (the CLI hashes it).
    /// Future evolution changes require an explicit checkpoint-version change.
    #[must_use]
    pub(super) fn checkpoint_v1_bytes(&self) -> Vec<u8> {
        let mut out = MAGIC.to_vec();
        let s = self.kernel.settings;
        word(&mut out, u64::from(s.level)); number(&mut out, s.volfrac);
        count(&mut out, s.iterations);
        for v in [s.band_cells, s.move_cells, s.ell0, s.mu_al, s.sobolev_alpha] { number(&mut out, v); }
        count(&mut out, s.nucleation_period);
        for v in [s.hole_radius_cells, s.youngs, s.poisson] { number(&mut out, v); }
        count(&mut out, self.controls.max_candidates);
        number(&mut out, self.controls.contraction);
        number(&mut out, self.controls.min_relative_improvement);
        count(&mut out, self.controls.max_solves);
        for v in [self.projection.target, self.projection.tolerance, self.projection.max_shift] {
            number(&mut out, v);
        }
        count(&mut out, self.projection.max_evaluations);
        word(&mut out, match self.kernel.aggregate {
            RobustAggregate::WeightedSum => 0, RobustAggregate::WorstWeightedCase => 1,
        });
        count(&mut out, self.kernel.load_cases.len());
        for case in &self.kernel.load_cases {
            word(&mut out, edge_tag(case.edge()));
            for v in case.interval().into_iter().chain(case.traction()).chain([case.weight()]) {
                number(&mut out, v);
            }
        }
        count(&mut out, self.fixed.len());
        for &(index, value) in &self.fixed { count(&mut out, index); number(&mut out, value); }
        match self.stress_limit {
            Some(limit) => {
                word(&mut out, 1); number(&mut out, limit.max_von_mises);
                number(&mut out, limit.absolute_tolerance);
            }
            None => word(&mut out, 0),
        }
        count(&mut out, self.next_iteration); number(&mut out, self.ell);
        count(&mut out, self.solves_started);
        for field in [&self.baseline_geometry, &self.current.phi] {
            for &value in field.nodes() { number(&mut out, value); }
        }
        state_bytes(&mut out, &self.baseline);
        state_bytes(&mut out, &self.current());
        stress_bytes(&mut out, self.baseline_stress.as_ref());
        stress_bytes(&mut out, self.current_stress.as_ref());
        out
    }

    /// Rehydrate without projecting a different geometry or resetting the study.
    ///
    /// Exactly two complete load families are reserved from the independent
    /// `remaining_recovery_solves` allowance before any PDE work. Every actual
    /// solve start decrements that allowance, including a solve that refuses;
    /// it is not refunded on a later evidence mismatch. Study work and limits
    /// remain exactly as recorded. This is synchronous recovery, not a bounded
    /// wall-time promise. Requires the same numerical implementation/platform.
    ///
    /// # Errors
    /// Refuses unsupported/malformed/bounded-count data before PDE work; refuses
    /// nonfeasible fields rather than silently repairing them; refuses any
    /// bitwise evidence disagreement after canonical solves/stress sampling.
    pub fn restore_checkpoint(
        bytes: &[u8], remaining_recovery_solves: &mut usize,
    ) -> Result<Self, CutFemError> {
        if bytes.len() > MAX_BYTES {
            return Err(invalid("unsupported or oversized projected checkpoint"));
        }
        let original_bytes = bytes;
        let (bytes, restoration_reduction, restoration_updates) = restoration::checkpoint_payload(bytes)?;
        if !bytes.starts_with(MAGIC) {
            return Err(invalid("unsupported or oversized projected checkpoint"));
        }
        let mut r = Reader { bytes, at: MAGIC.len() };
        let level = u32::try_from(r.count(8)?)
            .map_err(|_| invalid("invalid checkpoint level"))?;
        if level == 0 { return Err(invalid("checkpoint level must be positive")); }
        let settings = OptimizeSettings {
            level, volfrac: r.number()?, iterations: r.count(10_000)?,
            band_cells: r.number()?, move_cells: r.number()?, ell0: r.number()?,
            mu_al: r.number()?, sobolev_alpha: r.number()?,
            nucleation_period: r.count(usize::MAX)?, hole_radius_cells: r.number()?,
            youngs: r.number()?, poisson: r.number()?,
        };
        let controls = MultiLoadProjectedSettings {
            max_candidates: r.count(64)?, contraction: r.number()?,
            min_relative_improvement: r.number()?, max_solves: r.count(usize::MAX)?,
        };
        let projection = VolumeProjectionSettings {
            target: r.number()?, tolerance: r.number()?, max_shift: r.number()?,
            max_evaluations: r.count(128)?,
        };
        let aggregate = match r.word()? {
            0 => RobustAggregate::WeightedSum, 1 => RobustAggregate::WorstWeightedCase,
            _ => return Err(invalid("unknown checkpoint aggregate")),
        };
        let case_count = r.count(64)?;
        let mut cases = Vec::with_capacity(case_count);
        for _ in 0..case_count {
            let edge = match r.word()? {
                0 => DesignBoxEdge::Left, 1 => DesignBoxEdge::Right,
                2 => DesignBoxEdge::Bottom, 3 => DesignBoxEdge::Top,
                _ => return Err(invalid("unknown checkpoint load edge")),
            };
            cases.push(RobustLoadCase::new(edge, r.number()?, r.number()?,
                [r.number()?, r.number()?], r.number()?)?);
        }
        let n = 1usize << level;
        let fixed_count = r.count((n + 1) * (n + 1))?;
        let mut fixed = Vec::with_capacity(fixed_count);
        for _ in 0..fixed_count { fixed.push((r.count((n + 1) * (n + 1) - 1)?, r.number()?)); }
        let stress_limit = match r.word()? {
            0 => None, 1 => Some(SampledStressLimit::new(r.number()?, r.number()?)?),
            _ => return Err(invalid("unknown checkpoint stress policy")),
        };
        let ordinal = r.count(settings.iterations)?;
        let ell = r.number()?;
        let spent = r.count(controls.max_solves)?;
        let baseline_geometry = r.field(n)?;
        let geometry = r.field(n)?;
        // Evidence has fixed size derived from the admitted case/stress policy;
        // reject truncation/trailing data before spending a single solve.
        let state_len = 32 + 8 * case_count;
        let stress_len = if stress_limit.is_some() { 64 + 40 * case_count } else { 8 };
        if r.at + 2 * (state_len + stress_len) != bytes.len() {
            return Err(invalid("checkpoint evidence length mismatch"));
        }
        let minimum_spent = (ordinal + 1) * case_count;
        if case_count == 0 || settings.iterations == 0 || controls.max_candidates == 0
            || restoration_updates > ordinal
            || (restoration_reduction.is_some() && stress_limit.is_none())
            || projection.target.to_bits() != settings.volfrac.to_bits()
            || !(settings.move_cells.is_finite() && settings.move_cells > 0.0)
            || !(controls.contraction.is_finite() && controls.contraction > 0.0 && controls.contraction < 1.0)
            || !(controls.min_relative_improvement.is_finite() && (0.0..1.0).contains(&controls.min_relative_improvement))
            || spent < minimum_spent || !(ell.is_finite() && ell >= 0.0)
            || (ordinal == 0 && (ell.to_bits() != settings.ell0.to_bits()
                || !same_field(&baseline_geometry, &geometry)))
        {
            return Err(invalid("inconsistent projected checkpoint controls or continuation state"));
        }
        validate(&baseline_geometry, &cases, settings)?;
        validate(&geometry, &cases, settings)?;
        material(settings)?;
        for field in [&baseline_geometry, &geometry] {
            let mut checked = field.clone();
            project_material_volume(&mut checked, level, &fixed, projection)?;
            if !same_field(field, &checked) {
                return Err(invalid("checkpoint field violates fixed nodes or numerical area"));
            }
        }
        if *remaining_recovery_solves < 2 * case_count {
            return Err(invalid("checkpoint recovery needs two complete independent load families"));
        }
        let kernel = Kernel::new(&geometry, &cases, settings, aggregate)?;
        let mut evaluate = |field| -> Result<MultiState, CutFemError> {
            match kernel.evaluate_controlled(field, |_, complete| {
                if !complete { *remaining_recovery_solves -= 1; }
                ControlFlow::<Infallible>::Continue(())
            })? {
                ControlFlow::Continue(state) => Ok(state),
                ControlFlow::Break(never) => match never {},
            }
        };
        let baseline = evaluate(baseline_geometry.clone())?;
        let current = evaluate(geometry)?;
        if (baseline.volume - projection.target).abs() > projection.tolerance
            || (current.volume - projection.target).abs() > projection.tolerance
            || (restoration_reduction.is_none() && ordinal > 0 && !(current.objective < baseline.objective))
        {
            return Err(invalid("checkpoint independent feasibility/descent check failed"));
        }
        let mut restored = Self {
            kernel, current, baseline: MultiLoadProjectedState::of(&baseline), baseline_geometry,
            fixed, projection, controls, next_iteration: ordinal, ell, solves_started: spent,
            stress_limit, baseline_stress: None, current_stress: None,
            restoration_reduction, restoration_updates,
        };
        if stress_limit.is_some() {
            let sample = |state: &MultiState| -> Result<RobustSampledStressEvaluation, CutFemError> {
                match restored.sample_stress_controlled(state, |_, _| ControlFlow::<Infallible>::Continue(()))? {
                    ControlFlow::Continue(stress) => Ok(stress),
                    ControlFlow::Break(never) => match never {},
                }
            };
            let baseline_stress = sample(&baseline)?;
            let current_stress = sample(&restored.current)?;
            if restoration_reduction.is_none() {
                stress::require_feasible(stress_limit, Some(&baseline_stress))?;
                stress::require_feasible(stress_limit, Some(&current_stress))?;
            }
            restored.baseline_stress = Some(baseline_stress);
            restored.current_stress = Some(current_stress);
        }
        restored.check_restoration_history()?;
        if restored.checkpoint_bytes() != original_bytes {
            return Err(invalid("checkpoint evidence differs from independent replay; no optimizer restored"));
        }
        Ok(restored)
    }
}

#[cfg(test)]
mod tests;
