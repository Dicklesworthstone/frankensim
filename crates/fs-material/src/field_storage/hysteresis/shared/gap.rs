//! Uniform series-gap storage and exact source active-set path integration.
//!
//! In reference-linkage coordinates, s=M(Phi)+(k+R(g))*Phi and
//! R(g)=g/(mu*A*Nref^2). Along an affine linkage/gap chord,
//! Phi'=(delta(s)-delta(R)*Phi)/(M'+k+R). The positive denominator
//! preserves one flux direction. Each original play span has affine M(Phi),
//! hence a linear-over-linear Phi(t). Its first two analytic moments supply
//! winding work and the always-attractive force without a midpoint-force fit.

use super::super::branch_response;
use super::{
    FieldStorageError, HystereticCoreStep, SharedHystereticCore, SharedHystereticState,
    SharedHystereticStep, finite, invalid, positive,
};
use crate::field_storage::{MAGNETIC_PERMEABILITY_PROPERTY, ResolvedInductor};
use fs_matdb::EvaluationDecision;

/// Original constant gap-medium source, preparation geometry and open limits.
/// The accepted gap is in `SharedHystereticState`; this source retains its
/// original external conditions. Area, turns and permeability never evolve here.
#[derive(Clone, Debug, PartialEq)]
pub struct SharedHystereticGap {
    medium: ResolvedInductor,
    limits_m: [f64; 2],
    reluctance_per_m: f64,
}
impl SharedHystereticGap {
    pub(super) fn bind(
        medium: ResolvedInductor,
        reference_turns: f64,
        limits_m: [f64; 2],
    ) -> Result<Self, FieldStorageError> {
        if f64::from(medium.turns()) != reference_turns.abs()
            || limits_m.iter().any(|&x| !positive(x))
            || !(limits_m[0] < medium.path_length_m() && medium.path_length_m() < limits_m[1])
        {
            return Err(invalid(
                "series gap needs matching reference turns and positive open geometric limits",
            ));
        }
        let property = medium
            .material()
            .property(MAGNETIC_PERMEABILITY_PROPERTY)
            .ok_or_else(|| invalid("series gap lost its original permeability source"))?;
        if property.answer().receipt.decision != EvaluationDecision::ConstantWithinValidity
            || !property.is_static_at_fixed_external_conditions()
        {
            return Err(invalid(
                "series-gap motion needs a source-constant static permeability",
            ));
        }
        let turns = f64::from(medium.turns());
        let scale = finite(property.value_si() * medium.area_m2() * turns * turns)?;
        let reluctance_per_m = scale.recip();
        if !positive(scale)
            || !positive(reluctance_per_m)
            || limits_m.iter().any(|&g| !positive(g * reluctance_per_m))
        {
            return Err(invalid("series-gap reluctance is unrepresentable"));
        }
        Ok(Self {
            medium,
            limits_m,
            reluctance_per_m,
        })
    }
    pub const fn medium(&self) -> &ResolvedInductor {
        &self.medium
    }
    pub const fn limits_m(&self) -> [f64; 2] {
        self.limits_m
    }
    /// dR/dg = 1/(mu*A*Nref^2), in reciprocal henry-metres.
    pub const fn reluctance_per_m(&self) -> f64 {
        self.reluctance_per_m
    }
    pub(super) fn reluctance(&self, gap_m: f64) -> Result<f64, FieldStorageError> {
        if !gap_m.is_finite() || gap_m <= self.limits_m[0] || gap_m >= self.limits_m[1] {
            return Err(invalid(
                "shared magnetic gap is outside its open geometric interval",
            ));
        }
        finite(self.reluctance_per_m * gap_m)
    }
}

/// Original series-gap energy and force, positive displacement toward opening.
#[derive(Clone, Debug, PartialEq)]
pub struct SharedHystereticGapStep {
    pub gap_before_m: f64,
    pub gap_after_m: f64,
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    pub stored_change_j: f64,
    pub force_n: f64,
    /// Exact affine-chord mean force; never an endpoint-product force fit.
    pub mean_force_n: f64,
    /// Work delivered TO native mechanics: mean_force * delta(gap) [J].
    pub mechanical_work_j: f64,
    /// Original material active-set spans integrated, bounded by 4*branches+1.
    pub path_spans: usize,
}

pub(super) fn candidate(
    model: &SharedHystereticCore,
    previous: &SharedHystereticState,
    linkages: &[f64],
    gap_m: f64,
) -> Result<SharedHystereticStep, FieldStorageError> {
    if previous.geometry != model.geometry
        || linkages.len() != model.windings().len()
        || linkages.iter().any(|x| !x.is_finite())
    {
        return Err(invalid(
            "shared gap step requires original geometry and complete finite linkages",
        ));
    }
    let gap = model
        .gap()
        .ok_or_else(|| invalid("shared magnetic core has no physical gap"))?;
    let old_gap = previous
        .gap_m
        .ok_or_else(|| invalid("shared gap history is incomplete"))?;
    let (r0, r1) = (gap.reluctance(old_gap)?, gap.reluctance(gap_m)?);
    let dg = finite(gap_m - old_gap)?;
    let dr = finite(gap.reluctance_per_m * dg)?;
    let old_core = model.core.hold_candidate(previous.core_state())?;
    let mut s0 = 0.0;
    let mut s1 = 0.0;
    let mut ds = 0.0;
    for ((&w, &a), &b) in model
        .geometry
        .drive_weights
        .iter()
        .zip(previous.flux_linkages_wb())
        .zip(linkages)
    {
        s0 = finite(w.mul_add(a, s0))?;
        s1 = finite(w.mul_add(b, s1))?;
        ds = finite(w.mul_add(finite(b - a)?, ds))?;
    }
    let driving_change = finite((-dr).mul_add(old_core.flux_linkage_wb, ds))?;
    let core = if driving_change == 0.0 {
        old_core.clone()
    } else {
        let loading = finite(model.geometry.loading_per_h + r1)?;
        let target = finite((s1 / loading) / model.core.linkage_per_t)?;
        let coefficient = (model.core.current_per_field / model.core.linkage_per_t) / loading;
        if !positive(coefficient) {
            return Err(invalid("series-gap load line is unrepresentable"));
        }
        model.core.finish(model.core.law.candidate_linear_flux(
            previous.core_state(),
            target,
            coefficient,
        )?)?
    };
    let state = model.state_at(&core, linkages.to_vec(), Some(gap_m))?;
    let (p0, p1) = (old_core.flux_linkage_wb, core.flux_linkage_wb);
    let dp = finite(p1 - p0)?;
    if dp != 0.0 && driving_change.signum() != dp.signum() {
        return Err(invalid(
            "series-gap inverse changed the common monotone flux direction",
        ));
    }
    let gap_before = finite((0.5 * r0 * p0) * p0)?;
    let gap_after = finite((0.5 * r1 * p1) * p1)?;
    // Exact endpoint identity, with delta(R) taken from conserved displacement.
    let gap_change =
        finite(f64::midpoint(r0 * p0, r1 * p1).mul_add(dp, finite((0.5 * p0 * p1) * dr)?))?;
    let mut means: Vec<_> = previous
        .currents_a()
        .iter()
        .zip(state.currents_a())
        .map(|(&a, &b)| f64::midpoint(a, b))
        .collect();
    let mut leakage_change = 0.0;
    let mut mean_mmf = 0.0;
    for (((w, &i), &a), &b) in model
        .windings()
        .iter()
        .zip(&means)
        .zip(previous.flux_linkages_wb())
        .zip(linkages)
    {
        leakage_change =
            finite(i.mul_add(finite((-w.turns_ratio).mul_add(dp, b - a))?, leakage_change))?;
        mean_mmf = finite(w.turns_ratio.mul_add(i, mean_mmf))?;
    }
    let moments = path_moments(model, previous, &old_core, &core, s0, ds, r0, dr)?;
    let correction = if dg == 0.0 && ds != 0.0 {
        // Preserve the stable fixed-geometry integration-by-parts path.
        finite(mean_mmf.mul_add(dp, -finite(core.magnetic_work_j + gap_change)?) / ds)?
    } else {
        finite(moments.mean - f64::midpoint(p0, p1))?
    };
    for (&w, i) in model.geometry.drive_weights.iter().zip(&mut means) {
        *i = finite((-w).mul_add(correction, *i))?;
    }
    let force = finite((-0.5 * gap.reluctance_per_m * p1) * p1)?;
    let mean_force = finite(-0.5 * gap.reluctance_per_m * moments.square)?;
    let mechanical_work = finite(mean_force * dg)?;
    let core_change = finite(model.core.volume_m3 * core.material.energy_change_j_m3)?;
    let stored_change = finite(leakage_change + core_change + gap_change)?;
    let work = finite(stored_change + core.core_heat_j + mechanical_work)?;
    let mut paired = 0.0;
    // Endpoint storage also bounds representational error in the loaded source
    // inverse. Native circuit/mechanical owners impose their own absolute gates.
    let mut scale = finite(
        previous.stored_energy_j()
            + state.stored_energy_j()
            + leakage_change.abs()
            + core_change.abs()
            + gap_change.abs()
            + core.core_heat_j.abs()
            + mechanical_work.abs()
            + work.abs(),
    )?;
    for ((&i, &a), &b) in means.iter().zip(previous.flux_linkages_wb()).zip(linkages) {
        let term = finite(i * (b - a))?;
        paired = finite(paired + term)?;
        scale = finite(scale + term.abs())?;
    }
    if (paired - work).abs() > 4096.0 * f64::EPSILON * scale || mean_force > 0.0 {
        return Err(invalid(
            "shared gap winding, stored, pinning and mechanical work do not close",
        ));
    }
    Ok(SharedHystereticStep {
        stored_energy_j: state.stored_energy_j(),
        state,
        core,
        mean_currents_a: means,
        stored_change_j: stored_change,
        magnetic_work_j: work,
        gap: Some(SharedHystereticGapStep {
            gap_before_m: old_gap,
            gap_after_m: gap_m,
            stored_before_j: gap_before,
            stored_after_j: gap_after,
            stored_change_j: gap_change,
            force_n: force,
            mean_force_n: mean_force,
            mechanical_work_j: mechanical_work,
            path_spans: moments.spans,
        }),
    })
}

struct Moments {
    mean: f64,
    square: f64,
    spans: usize,
}

// Positive slope sums are recomputed through a balanced tree after each event.
// Unlike subtracting susceptibilities from a running sum, this retains the
// small positive background when large branches enter/leave saturation.
struct Slopes {
    leaf: usize,
    values: Vec<f64>,
}
impl Slopes {
    fn new(values: impl Iterator<Item = f64>, count: usize) -> Self {
        let leaf = count.max(1).next_power_of_two();
        let mut data = vec![0.0; 2 * leaf];
        for (i, value) in values.enumerate() {
            data[leaf + i] = value;
        }
        for i in (1..leaf).rev() {
            data[i] = data[2 * i] + data[2 * i + 1];
        }
        Self { leaf, values: data }
    }
    fn set(&mut self, index: usize, value: f64) {
        let mut at = self.leaf + index;
        self.values[at] = value;
        while at > 1 {
            at /= 2;
            self.values[at] = self.values[2 * at] + self.values[2 * at + 1];
        }
    }
}

fn add(sum: &mut f64, correction: &mut f64, value: f64) -> Result<(), FieldStorageError> {
    let y = finite(value - *correction)?;
    let next = finite(*sum + y)?;
    *correction = finite((next - *sum) - y)?;
    *sum = next;
    Ok(())
}

// Mean shape f=(1+z)u/(1+zu) and its square, with positive denominator ratio.
// Series around z=0 avoid cancellation in logarithmic divided differences.
fn shape_moments(z: f64) -> Result<(f64, f64), FieldStorageError> {
    if !z.is_finite() || z < 0.0 {
        return Err(invalid("gap path denominator ratio is unrepresentable"));
    }
    let (mut first, mut second) = (0.5, 1.0 / 3.0);
    if z <= 0.25 {
        let mut power = z;
        for j in 1..=32 {
            let n = f64::from(j);
            first += power / ((n + 1.0) * (n + 2.0));
            second += 4.0 * power / ((n + 1.0) * (n + 2.0) * (n + 3.0));
            power *= -z;
        }
    } else {
        let log_over = z.ln_1p() / z;
        let ratio = 1.0 + z.recip();
        first = ratio * (1.0 - log_over);
        second = ratio * ratio * (1.0 - 2.0 * log_over + 1.0 / (1.0 + z));
    }
    if !first.is_finite()
        || !second.is_finite()
        || first < 0.5
        || first > 1.0
        || second < 1.0 / 3.0
        || second > 1.0
    {
        return Err(invalid("gap path moments are unrepresentable"));
    }
    Ok((first, second))
}

fn span_moments(
    mut a: f64,
    mut b: f64,
    mut da: f64,
    mut db: f64,
) -> Result<(f64, f64), FieldStorageError> {
    if !positive(da) || !positive(db) {
        return Err(invalid("gap path lost its positive load line"));
    }
    if db < da {
        std::mem::swap(&mut a, &mut b);
        std::mem::swap(&mut da, &mut db);
    }
    let z = finite((db - da) / da)?;
    let (m1, m2) = shape_moments(z)?;
    let delta = finite(b - a)?;
    let mean = finite(delta.mul_add(m1, a))?;
    let variance = m2 - m1 * m1;
    if variance < -64.0 * f64::EPSILON {
        return Err(invalid("gap path has negative flux variance"));
    }
    let square = finite((delta * delta).mul_add(variance.max(0.0), mean * mean))?;
    Ok((mean, square))
}

#[allow(clippy::too_many_arguments)]
fn path_moments(
    model: &SharedHystereticCore,
    previous: &SharedHystereticState,
    old: &HystereticCoreStep,
    next: &HystereticCoreStep,
    s0: f64,
    ds: f64,
    r0: f64,
    dr: f64,
) -> Result<Moments, FieldStorageError> {
    let (p0, p1) = (old.flux_linkage_wb, next.flux_linkage_wb);
    if p0 == p1 {
        return Ok(Moments {
            mean: p0,
            square: finite(p0 * p0)?,
            spans: 0,
        });
    }
    let law = &model.core.law;
    let parameters = &law.parameters;
    let old_z = previous.core_state().polarizations_t();
    let (h0, h1) = (
        old.material.state.field_a_m(),
        next.material.state.field_a_m(),
    );
    let (low, high) = (h0.min(h1), h0.max(h1));
    let mut events = Vec::with_capacity(4 * parameters.branches.len());
    for (index, (branch, &z)) in parameters.branches.iter().zip(old_z).enumerate() {
        let chi = branch.susceptibility_h_m;
        let hc = branch.pinning_field_a_m;
        let saturation = branch.saturation_polarization_t / chi;
        for field in [
            z / chi - hc,
            z / chi + hc,
            -saturation - hc,
            saturation + hc,
        ] {
            if field > low && field < high {
                events.push((field, index));
            }
        }
    }
    events.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
    events.dedup();
    let first = events.first().map_or(high, |event| event.0);
    let sample = f64::midpoint(low, first);
    let mut slopes = Slopes::new(
        parameters
            .branches
            .iter()
            .zip(old_z)
            .map(|(branch, &z)| branch_response(branch, z, sample).1),
        parameters.branches.len(),
    );
    let mut field = low;
    let mut density = if low == h0 {
        old.material.state.flux_density_t()
    } else {
        next.material.state.flux_density_t()
    };
    let mut density_correction = 0.0;
    let mut phi = if low == h0 { p0 } else { p1 };
    let mut time = if low == h0 { 0.0 } else { 1.0 };
    let mut cursor = 0;
    let mut result = Moments {
        mean: 0.0,
        square: 0.0,
        spans: 0,
    };
    let (mut mean_correction, mut square_correction, mut covered) = (0.0, 0.0, 0.0);
    while field < high {
        let end = events.get(cursor).map_or(high, |event| event.0);
        let slope = finite(parameters.background_permeability_h_m + slopes.values[1])?;
        let q = finite((model.core.current_per_field / model.core.linkage_per_t) / slope)?;
        if !positive(q) {
            return Err(invalid(
                "gap source active-set reluctance is unrepresentable",
            ));
        }
        add(
            &mut density,
            &mut density_correction,
            finite(slope * (end - field))?,
        )?;
        let end_phi = if end == h0 {
            p0
        } else if end == h1 {
            p1
        } else {
            finite(model.core.linkage_per_t * density)?
        };
        let end_time = if end == h0 {
            0.0
        } else if end == h1 {
            1.0
        } else {
            let numerator = finite(
                (model.geometry.loading_per_h + r0)
                    .mul_add(end_phi, model.core.current_per_field * end - s0),
            )?;
            let denominator = finite((-dr).mul_add(end_phi, ds))?;
            let value = finite(numerator / denominator)?;
            if !(0.0..=1.0).contains(&value) {
                return Err(invalid(
                    "source threshold cannot be resolved on the shared gap chord",
                ));
            }
            value
        };
        let width = finite((end_time - time).abs())?;
        if width > 0.0 {
            let base = finite(model.geometry.loading_per_h + r0 + q)?;
            let (mean, square) = span_moments(
                phi,
                end_phi,
                finite(dr.mul_add(time, base))?,
                finite(dr.mul_add(end_time, base))?,
            )?;
            add(
                &mut result.mean,
                &mut mean_correction,
                finite(width * mean)?,
            )?;
            add(
                &mut result.square,
                &mut square_correction,
                finite(width * square)?,
            )?;
            covered = finite(covered + width)?;
            result.spans += 1;
        }
        field = end;
        phi = end_phi;
        time = end_time;
        let start = cursor;
        while cursor < events.len() && events[cursor].0 == end {
            cursor += 1;
        }
        let following = events.get(cursor).map_or(high, |event| event.0);
        let sample = f64::midpoint(end, following);
        for &(_, index) in &events[start..cursor] {
            slopes.set(
                index,
                branch_response(&parameters.branches[index], old_z[index], sample).1,
            );
        }
    }
    if (covered - 1.0).abs() > 256.0 * f64::EPSILON * (result.spans as f64).max(1.0) {
        return Err(invalid(
            "source active-set spans do not cover the shared gap chord",
        ));
    }
    Ok(result)
}
