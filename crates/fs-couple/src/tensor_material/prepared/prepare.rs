//! Bounded source assembly and immutable native flow/positive-loss caches.
use super::*;
use fs_phs::linear_exact::{LinearSplitBudget, MAX_LINEAR_EXPONENTIAL_BLOCK_DIM, split_quadratic};
use fs_phs::{PortHamiltonian, QuadraticStorage};

impl PreparedTensorBody {
    /// Retain original consistent inertia and every projected source point.
    /// `heat_groups[p]` names a channel for point p; labels must cover 0..groups
    /// without gaps. Grouping changes neither mechanics nor material history.
    ///
    /// Equal source times share y_tau exactly, including distinct weak and
    /// rank-deficient branches. Every six-row source factor remains available
    /// for original storage and positive heat. Dense limits refuse the entire
    /// request; they never reduce quadrature, material memory or coordinates.
    #[allow(clippy::too_many_arguments)]
    pub fn from_mass_matrix(
        coordinates: usize,
        mass_kg: Vec<f64>,
        points: Vec<ProjectedTensorProny>,
        heat_groups: Vec<usize>,
        dt_s: f64,
        controls: PreparedTensorControls,
        gate: &CancelGate,
    ) -> Result<Self> {
        Self::prepare_sources(
            coordinates,
            mass_kg,
            points,
            heat_groups,
            dt_s,
            None,
            controls,
            gate,
        )
    }

    /// Prepare continuous F(s)=F_start+(F_end-F_start)*s, s=t/h. The positive
    /// `force_scale_n` conditions additional force records only; it changes no
    /// physical source amplitude, inertia or material storage. Native capacity
    /// includes these records and a separately bounded one-way moment observer.
    #[allow(clippy::too_many_arguments)]
    pub fn from_mass_matrix_affine(
        coordinates: usize,
        mass_kg: Vec<f64>,
        points: Vec<ProjectedTensorProny>,
        heat_groups: Vec<usize>,
        dt_s: f64,
        force_scale_n: f64,
        controls: PreparedTensorControls,
        gate: &CancelGate,
    ) -> Result<Self> {
        Self::prepare_sources(
            coordinates,
            mass_kg,
            points,
            heat_groups,
            dt_s,
            Some(force_scale_n),
            controls,
            gate,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn prepare_sources(
        coordinates: usize,
        mass_kg: Vec<f64>,
        points: Vec<ProjectedTensorProny>,
        heat_groups: Vec<usize>,
        dt_s: f64,
        force_scale_n: Option<f64>,
        controls: PreparedTensorControls,
        gate: &CancelGate,
    ) -> Result<Self> {
        poll(gate)?;
        let n = coordinates;
        let source_blocks = usize::from(force_scale_n.is_some());
        if n == 0
            || n > 128
            || !dt_s.is_finite()
            || dt_s <= 0.0
            || controls.maximum_steps == 0
            || controls.maximum_steps > (1_u64 << 53)
            || !controls.maximum_relative_operator_error.is_finite()
            || !(0.0..1.0).contains(&controls.maximum_relative_operator_error)
            || !controls.maximum_work_error_j.is_finite()
            || controls.maximum_work_error_j < 0.0
            || points.is_empty()
            || heat_groups.len() != points.len()
            || mass_kg.len() != n * n
            || !finite(&mass_kg)
            || force_scale_n.is_some_and(|scale| {
                !scale.is_finite()
                    || scale <= 0.0
                    || !scale.recip().is_finite()
                    || !dt_s.recip().is_finite()
            })
        {
            return Err(invalid(
                "finite complete source, clock and controls required",
            ));
        }
        bound(n as u128, controls.maximum_coordinates as u128)?;
        bound(
            points.len() as u128,
            controls.maximum_material_points as u128,
        )?;
        let groups = heat_groups
            .iter()
            .copied()
            .max()
            .unwrap()
            .checked_add(1)
            .ok_or_else(|| invalid("heat channel extent overflow"))?;
        bound(groups as u128, controls.maximum_heat_groups as u128)?;
        if groups > points.len() {
            return Err(invalid("heat groups must be contiguous and nonempty"));
        }
        // Even the preliminary grouping scratch is admitted BEFORE allocation.
        // A fixed-capacity sorted vector needs no unbounded tree nodes and keeps
        // the exact positive source-time ordering without changing arithmetic.
        bound(
            ((2 + source_blocks) * n) as u128,
            controls.maximum_state_coordinates as u128,
        )?;
        bound(
            (2 * (2 + source_blocks) * n) as u128,
            MAX_LINEAR_EXPONENTIAL_BLOCK_DIM as u128,
        )?;
        let source_branches = points
            .iter()
            .try_fold(0usize, |sum, point| sum.checked_add(point.branches().len()))
            .ok_or_else(|| invalid("source branch count overflow"))?;
        let time_capacity = (controls
            .maximum_state_coordinates
            .min(MAX_LINEAR_EXPONENTIAL_BLOCK_DIM / 2)
            / n)
            .saturating_sub(2 + source_blocks)
            .min(source_branches);
        bound(
            128 + 9 * groups as u128 + 8 * time_capacity as u128,
            controls.maximum_preparation_bytes as u128,
        )?;
        bound(
            8 * (points.len() as u128 + source_branches as u128 + groups as u128 + 1),
            controls.maximum_preparation_multiply_terms as u128,
        )?;
        let mut seen = allocate(groups, false)?;
        let mut row_counts = allocate(groups, 0usize)?;
        let mut times = Vec::<f64>::new();
        times
            .try_reserve_exact(time_capacity)
            .map_err(|_| invalid("source-time allocation refused"))?;
        let (mut branches, mut kelvin_points, mut compatible_points) = (0usize, 0usize, 0usize);
        for (point, &group) in points.iter().zip(&heat_groups) {
            poll(gate)?;
            if point.coordinates() != n || point.free_strain_mandel().iter().any(|x| *x != 0.0) {
                return Err(invalid(
                    "complete matching source coordinates and zero eigenstrain required",
                ));
            }
            compatible_points += usize::from(point.compatible_strain_columns().is_some());
            seen[group] = true;
            if let Some(rows) = point.viscosity_gram_rows() {
                kelvin_points += 1;
                row_counts[group] = row_counts[group]
                    .checked_add(rows.len())
                    .ok_or_else(|| invalid("source row count overflow"))?;
            }
            for branch in point.branches() {
                poll(gate)?;
                let tau = branch.relaxation_time_s();
                if !tau.is_finite() || tau <= 0.0 || !tau.recip().is_finite() || tau.recip() <= 0.0
                {
                    return Err(invalid(
                        "positive representable original relaxation rate required",
                    ));
                }
                if let Err(at) = times.binary_search_by_key(&tau.to_bits(), |t| t.to_bits()) {
                    let states = (3usize + source_blocks + times.len())
                        .checked_mul(n)
                        .ok_or_else(|| invalid("complete history extent overflow"))?;
                    bound(states as u128, controls.maximum_state_coordinates as u128)?;
                    bound(2 * states as u128, MAX_LINEAR_EXPONENTIAL_BLOCK_DIM as u128)?;
                    // The capacity was proved from these same state limits.
                    if times.len() == time_capacity {
                        return Err(invalid("complete source-time capacity exhausted"));
                    }
                    times.insert(at, tau);
                }
                branches = branches
                    .checked_add(1)
                    .ok_or_else(|| invalid("branch count overflow"))?;
                row_counts[group] = row_counts[group]
                    .checked_add(branch.rows().len())
                    .ok_or_else(|| invalid("source row count overflow"))?;
            }
        }
        if seen.iter().any(|x| !x) {
            return Err(invalid("heat groups must be contiguous and nonempty"));
        }
        let physical_states = (2 + times.len()) * n;
        let states = physical_states + source_blocks * n;
        let dimension = states + n + 1;
        bound(states as u128, controls.maximum_state_coordinates as u128)?;
        bound(2 * states as u128, MAX_LINEAR_EXPONENTIAL_BLOCK_DIM as u128)?;
        // The unchanged base split and the extended/native loss each have their
        // own complete dimension proposition; no cap is silently widened.
        bound(
            (2 * (3 * n + 1)) as u128,
            MAX_LINEAR_EXPONENTIAL_BLOCK_DIM as u128,
        )?;
        bound(dimension as u128, MAX_LINEAR_EXPONENTIAL_BLOCK_DIM as u128)?;
        if force_scale_n.is_some() {
            bound(
                (states + 2 * n) as u128,
                MAX_LINEAR_EXPONENTIAL_BLOCK_DIM as u128,
            )?;
        }
        let rows = row_counts
            .iter()
            .try_fold(0usize, |a, b| a.checked_add(*b))
            .ok_or_else(|| invalid("total source row count overflow"))?;
        admit_resources(
            n,
            states,
            points.len(),
            branches,
            kelvin_points,
            compatible_points,
            groups,
            rows,
            force_scale_n.is_some(),
            controls,
        )?;
        let budget = DenseExponentialBudget {
            max_dimension: MAX_LINEAR_EXPONENTIAL_BLOCK_DIM,
            max_dense_bytes: controls.maximum_preparation_bytes,
            max_multiply_terms: controls.maximum_preparation_multiply_terms,
        };
        let inertia = Inertia::consistent(n, mass_kg, gate)?;
        let mass = match &inertia {
            Inertia::Consistent { mass, .. } => mass,
            Inertia::Diagonal(_) => return Err(invalid("consistent inertia owner required")),
        };
        for i in 0..n {
            poll(gate)?;
            for j in 0..n {
                let product = (0..n).fold(0.0, |sum, k| {
                    mass[i * n + k].mul_add(inertia.inverse_entry(k, j), sum)
                });
                let defect = (product - if i == j { 1.0 } else { 0.0 }).abs();
                if !defect.is_finite() || defect > controls.maximum_relative_operator_error {
                    return Err(PreparedTensorError::Operator {
                        defect,
                        allowance: controls.maximum_relative_operator_error,
                    });
                }
            }
        }
        let mut equilibrium = zeros(n * n)?;
        let mut viscosity = zeros(n * n)?;
        let mut equilibrium_roundoff = zeros(n * n)?;
        let mut viscosity_roundoff = zeros(n * n)?;
        let mut branch_stiffness = zeros(times.len() * n * n)?;
        let mut branch_roundoff = zeros(times.len() * n * n)?;
        let mut point_times = Vec::with_capacity(points.len());
        for point in &points {
            poll(gate)?;
            for i in 0..n * n {
                add(
                    &mut equilibrium[i],
                    &mut equilibrium_roundoff[i],
                    point.equilibrium_stiffness()[i],
                );
                add(
                    &mut viscosity[i],
                    &mut viscosity_roundoff[i],
                    point.viscosity_matrix()[i],
                );
            }
            let mut indices = Vec::with_capacity(point.branches().len());
            for branch in point.branches() {
                poll(gate)?;
                let index = times
                    .binary_search_by_key(&branch.relaxation_time_s().to_bits(), |t| t.to_bits())
                    .map_err(|_| invalid("original source time missing"))?;
                indices.push(index);
                for row in branch.rows() {
                    for i in 0..n {
                        for j in 0..n {
                            let at = index * n * n + i * n + j;
                            add(
                                &mut branch_stiffness[at],
                                &mut branch_roundoff[at],
                                row[i] * row[j],
                            );
                        }
                    }
                }
            }
            point_times.push(indices);
        }
        if !finite(&equilibrium) || !finite(&viscosity) || !finite(&branch_stiffness) {
            return Err(invalid("original source assembly overflow"));
        }
        let base = 2 * n;
        let mut q = zeros(base * base)?;
        let mut j = zeros(base * base)?;
        let mut r = zeros(base * base)?;
        let mut g = zeros(base * n)?;
        for i in 0..n {
            poll(gate)?;
            j[i * base + n + i] = 1.0;
            j[(n + i) * base + i] = -1.0;
            g[(n + i) * n + i] = 1.0;
            for k in 0..n {
                q[i * base + k] = equilibrium[i * n + k];
                q[(n + i) * base + n + k] = inertia.inverse_entry(i, k);
                r[(n + i) * base + n + k] = viscosity[i * n + k];
            }
        }
        let storage = QuadraticStorage::new(q, base)?;
        let system = PortHamiltonian::new(base, n, j, r, g, Box::new(storage))?;
        poll(gate)?;
        let split = split_quadratic(
            &system,
            LinearSplitBudget {
                max_modes: base,
                max_omega_dt: core::f64::consts::PI,
            },
        )?;
        poll(gate)?;
        let mut generator = zeros(states * states)?;
        let mut input = zeros(states * n)?;
        for i in 0..n {
            poll(gate)?;
            if force_scale_n.is_some() {
                // The held input is the interval force INCREMENT [N]. The
                // force record starts at F0 and evolves with derivative dF/h.
                input[(physical_states + i) * n + i] = dt_s.recip();
                generator[(n + i) * states + physical_states + i] = 1.0;
            } else {
                input[(n + i) * n + i] = 1.0;
            }
            for k in 0..n {
                generator[i * states + n + k] = inertia.inverse_entry(i, k);
                generator[(n + i) * states + k] = -equilibrium[i * n + k];
                generator[(n + i) * states + n + k] = (0..n).fold(0.0, |sum, l| {
                    (-viscosity[i * n + l]).mul_add(inertia.inverse_entry(l, k), sum)
                });
                for (t, &tau) in times.iter().enumerate() {
                    let y = (2 + t) * n;
                    generator[(n + i) * states + y + k] = -branch_stiffness[t * n * n + i * n + k];
                    generator[(y + i) * states + n + k] = inertia.inverse_entry(i, k);
                    if i == k {
                        generator[(y + i) * states + y + k] = -tau.recip();
                    }
                }
            }
        }
        if !finite(&generator) {
            return Err(invalid("hereditary generator overflow"));
        }
        let mut heat_scales = zeros(states)?;
        for i in 0..n {
            heat_scales[i] = scale(equilibrium[i * n + i])?;
            heat_scales[n + i] = scale(inertia.inverse_entry(i, i))?;
            for t in 0..times.len() {
                heat_scales[(2 + t) * n + i] = scale(branch_stiffness[t * n * n + i * n + i])?;
            }
            if let Some(scale) = force_scale_n {
                heat_scales[physical_states + i] = scale.recip();
            }
        }
        let flow = if force_scale_n.is_some() {
            let mut work = zeros(states * states)?;
            for i in 0..n {
                poll(gate)?;
                for j in 0..n {
                    let value = 0.5 * inertia.inverse_entry(i, j);
                    work[(physical_states + i) * states + n + j] = value;
                    work[(n + j) * states + physical_states + i] = value;
                }
            }
            // This observes the original F(t)^T M^-1 p(t), independently of
            // energy or the separately prepared linear displacement moment.
            PreparedQuadraticFlow::new_extended_input_controlled(
                &system,
                &split,
                &generator,
                &heat_scales[base..],
                &input,
                dt_s,
                &[work],
                &[],
                budget,
                || gate.is_requested(),
            )?
        } else {
            PreparedQuadraticFlow::new_extended_controlled(
                &system,
                &split,
                &generator,
                &heat_scales[base..],
                dt_s,
                &[],
                &[],
                budget,
                || gate.is_requested(),
            )?
        };
        // Independent positive-loss caches use a diagonal similarity only.
        // Source rows are changed contragrediently; physical squares stay exact.
        for i in 0..states {
            poll(gate)?;
            for j in 0..states {
                generator[i * states + j] =
                    (generator[i * states + j] / heat_scales[j]) * heat_scales[i];
            }
            for p in 0..n {
                input[i * n + p] *= heat_scales[i];
            }
        }
        if !finite(&generator) || !finite(&input) {
            return Err(invalid("conditioned original heat trajectory overflow"));
        }
        let affine = force_scale_n
            .map(|force_scale_n| {
                super::affine::prepare_observer(
                    n,
                    &generator,
                    &input,
                    &heat_scales,
                    dt_s,
                    force_scale_n,
                    budget,
                    gate,
                )
            })
            .transpose()?;
        let mut sources = Vec::with_capacity(groups);
        for &count in &row_counts {
            sources.push(zeros(count * states)?);
        }
        let mut offsets = vec![0usize; groups];
        for ((point, times), &group) in points.iter().zip(&point_times).zip(&heat_groups) {
            poll(gate)?;
            if let Some(rows) = point.viscosity_gram_rows() {
                for row in rows {
                    let at = offsets[group] * states;
                    for j in 0..n {
                        sources[group][at + n + j] = (0..n).fold(0.0, |sum, k| {
                            row[k].mul_add(inertia.inverse_entry(k, j), sum)
                        }) / heat_scales[n + j];
                    }
                    offsets[group] += 1;
                }
            }
            for (branch, &time) in point.branches().iter().zip(times) {
                let root_rate = branch.relaxation_time_s().recip().sqrt();
                for row in branch.rows() {
                    let at = offsets[group] * states + (2 + time) * n;
                    for j in 0..n {
                        let value = (row[j] / heat_scales[(2 + time) * n + j]) * root_rate;
                        if row[j] != 0.0 && value == 0.0 {
                            return Err(invalid(
                                "source heat scaling lost a nonzero material factor",
                            ));
                        }
                        sources[group][at + j] = value;
                    }
                    offsets[group] += 1;
                }
            }
        }
        let mut heat = Vec::with_capacity(groups);
        for (rows, count) in sources.into_iter().zip(row_counts) {
            poll(gate)?;
            let factor = compress(
                rows,
                count,
                states,
                controls.maximum_relative_operator_error,
                gate,
            )?;
            heat.push(PreparedSourceLoss::new_controlled(
                states,
                n,
                &generator,
                &input,
                &factor,
                dt_s,
                budget,
                || gate.is_requested(),
            )?);
        }
        poll(gate)?;
        Ok(Self {
            inertia,
            points,
            heat_groups,
            times,
            point_times,
            equilibrium,
            flow,
            heat,
            heat_scales,
            affine,
            states: physical_states,
            dt_s,
            controls,
        })
    }
}

fn allocate<T: Copy>(count: usize, value: T) -> Result<Vec<T>> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(count)
        .map_err(|_| invalid("grouping allocation refused"))?;
    result.resize(count, value);
    Ok(result)
}

fn scale(diagonal: f64) -> Result<f64> {
    if !diagonal.is_finite() || diagonal < 0.0 {
        return Err(invalid("invalid original storage diagonal"));
    }
    let value = if diagonal == 0.0 {
        1.0
    } else {
        diagonal.sqrt()
    };
    if !value.is_finite() || value <= 0.0 || !value.recip().is_finite() {
        return Err(invalid("source coordinate scale is unrepresentable"));
    }
    Ok(value)
}

/// Conservative whole-owner simultaneous storage and product envelopes. Native
/// flow uses 624*d^3, each loss uses 2048*(2*s)^3; every channel is counted.
/// Additional terms include source storage, native base eigensolves, full source
/// assembly, tall QR and original-Gram admission. Warm work contains no factor,
/// eigensolve, time quadrature, exponential or iteration.
#[allow(clippy::too_many_arguments)]
fn admit_resources(
    n: usize,
    s: usize,
    p: usize,
    b: usize,
    kv: usize,
    compatible: usize,
    g: usize,
    rows: usize,
    affine: bool,
    controls: PreparedTensorControls,
) -> Result<()> {
    let (n, s, p, b, kv, compatible, g, rows) = (
        n as u128,
        s as u128,
        p as u128,
        b as u128,
        kv as u128,
        compatible as u128,
        g as u128,
        rows as u128,
    );
    let d = s + n + 1;
    let source = 8 * (2 * p * n * n + 6 * (p + b + kv + compatible) * n + 256 * p + 72 * b)
        + 1024 * p
        + 512 * b;
    let mut bytes = source
        + 8 * (36 * d * d
            + g * (132 * s * s + s * n)
            + 32 * s * s
            + 8 * s * n
            + 4 * rows * s
            + 256 * (s + p + b));
    let mut work = 624 * d * d * d
        + 2048 * g * (2 * s) * (2 * s) * (2 * s)
        + 128 * (2 * n) * (2 * n) * (2 * n) * (2 * n)
        + 16 * (rows + g * s) * s * s
        + 64 * (p + b + 1) * n * n;
    // Assumed-strain points also contract and admit their independent original
    // compatible map. It adds storage and endpoint work, never material energy.
    let mut warm =
        64 * s * s + g * (16 * s * s + s * n) + 64 * (p + b + compatible + 1) * n + 16 * n * n;
    if affine {
        let observer = s + 2 * n;
        // Original force-work Gramian plus one native one-way integral
        // exponential and its retained linear maps. No n-fold Gramian cache.
        bytes += 8 * (8 * observer * observer + 2 * d * d + 4 * n * (s + n) + 256 * (s + n));
        work += 80 * observer * observer * observer + 32 * observer * observer;
        warm += 8 * d * d + 16 * n * (s + n);
    }
    bound(bytes, controls.maximum_preparation_bytes as u128)?;
    bound(work, controls.maximum_preparation_multiply_terms as u128)?;
    bound(warm, controls.maximum_step_multiply_terms as u128)
}

/// Native scaled Householder QR preserves the complete positive row space.
/// No relative rank decision is made. The returned factor is admitted against
/// the ORIGINAL normalized source Gram, including weak off-diagonal directions.
fn compress(
    mut rows: Vec<f64>,
    height: usize,
    n: usize,
    tolerance: f64,
    gate: &CancelGate,
) -> Result<Vec<f64>> {
    if !finite(&rows) {
        return Err(invalid("original heat source overflow"));
    }
    if height <= n {
        return Ok(rows);
    }
    let mut scales = zeros(n)?;
    for row in rows.chunks_exact(n) {
        poll(gate)?;
        for j in 0..n {
            scales[j] = scales[j].max(row[j].abs());
        }
    }
    for row in rows.chunks_exact_mut(n) {
        poll(gate)?;
        for j in 0..n {
            if scales[j] != 0.0 {
                let before = row[j];
                row[j] /= scales[j];
                if before != 0.0 && row[j] == 0.0 {
                    return Err(invalid("heat factor equilibration lost a nonzero entry"));
                }
            }
        }
    }
    poll(gate)?;
    let qr = fs_la::factor::qr(&rows, height, n);
    poll(gate)?;
    for i in 0..n {
        poll(gate)?;
        let left = (0..height)
            .map(|k| rows[k * n + i] * rows[k * n + i])
            .sum::<f64>()
            .sqrt();
        for j in 0..=i {
            let original = (0..height)
                .map(|k| rows[k * n + i] * rows[k * n + j])
                .sum::<f64>();
            let actual = (0..n).map(|k| qr.r(k, i) * qr.r(k, j)).sum::<f64>();
            let right = (0..height)
                .map(|k| rows[k * n + j] * rows[k * n + j])
                .sum::<f64>()
                .sqrt();
            let allowance = tolerance * left * right;
            let defect = (actual - original).abs();
            if !defect.is_finite() || defect > allowance {
                return Err(PreparedTensorError::Operator { defect, allowance });
            }
        }
    }
    let mut factor = zeros(n * n)?;
    for i in 0..n {
        poll(gate)?;
        for j in i..n {
            let value = qr.r(i, j) * scales[j];
            if !value.is_finite() || qr.r(i, j) != 0.0 && scales[j] != 0.0 && value == 0.0 {
                return Err(invalid("source heat factor scale is unrepresentable"));
            }
            factor[i * n + j] = value;
        }
    }
    Ok(factor)
}
