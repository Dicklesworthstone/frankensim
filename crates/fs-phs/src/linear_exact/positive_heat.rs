//! Preserve the positive factor of the physical resistance integral.
//!
//! Expanding W into z^T W z can erase the heat of a tiny vibration about a
//! force-supported equilibrium. We instead retain W = L^T L and evaluate
//! ||L z||^2, using initial loss effort and the actual force imbalance as
//! coordinates. Positive Gauss weights on a norm-bounded short interval and
//! QR doubling preserve this factor without clipping heat or resistance modes.
use super::{LinearSplit, checked_eigh, exponential, finite, product, refuse};
use crate::interconnection::accumulation::Sum;
use crate::{PhsError, PortHamiltonian, PreparedStepError};
use fs_la::factor::qr;
use fs_math::{det, eft::two_sum};

pub(super) struct PositiveHeat {
    n: usize,
    m: usize,
    d: usize,
    forward: Vec<f64>,
    resistance_root: Vec<f64>,
    rate: f64,
    root: Vec<f64>,
    /// Doublings of the Gramian chain that integrated the dissipated heat (`scaling_doublings`);
    /// the step's energy-balance gate reads it.
    doublings: u32,
}

impl PositiveHeat {
    /// Same physical generator and energy transform as PreparedQuadraticFlow.
    /// The caller has already admitted the complete state/port capacity.
    pub(super) fn new(
        sys: &PortHamiltonian,
        split: &LinearSplit,
        dt: f64,
        resistance: &[f64],
    ) -> Result<Self, PhsError> {
        if resistance.len() != split.n * split.n || !finite(resistance) {
            return Err(refuse("positive heat requires complete finite resistance"));
        }
        // Passive resistance graphs have an explicit positive decomposition.
        // Their exact effort nullspaces must not depend on a rounded eigenvalue.
        let root = match &sys.resistance_witness {
            Some(witness) if witness.matrix == resistance => witness.root.clone(),
            _ => factor_resistance(resistance, split.n)?,
        };
        Self::from_factor(sys, split, dt, resistance, root)
    }

    // A physical port already supplies its rank-one factor. Keep that factor
    // instead of reconstructing its exact nullspace through eigendecomposition.
    // The same entrywise reconstruction admission still checks the supplied R.
    pub(super) fn from_factor(
        sys: &PortHamiltonian,
        split: &LinearSplit,
        dt: f64,
        resistance: &[f64],
        resistance_root: Vec<f64>,
    ) -> Result<Self, PhsError> {
        let (n, m) = (split.n, split.m);
        let d = 2 * n;
        if resistance.len() != n * n
            || resistance_root.len() != n * n
            || !finite(resistance)
            || !finite(&resistance_root)
            || !dt.is_finite()
            || dt <= 0.0
        {
            return Err(refuse(
                "positive heat requires finite resistance and a positive clock",
            ));
        }
        check_resistance_root(resistance, &resistance_root, n)?;
        let mut forward = vec![0.0; n * n];
        let mut inverse = vec![0.0; n * n];
        for i in 0..n {
            let scale = split.coordinate_scale(i);
            forward[i * n + i] = scale;
            inverse[i * n + i] = 1.0 / scale;
        }
        let scale = resistance.iter().fold(0.0_f64, |s, r| s.max(r.abs()));
        if scale == 0.0 {
            return Ok(Self {
                n,
                m,
                d,
                forward,
                resistance_root: Vec::new(),
                rate: 1.0,
                root: Vec::new(),
                doublings: 0,
            });
        }
        let mut original = vec![0.0; n * n];
        for i in 0..n {
            for k in 0..n {
                for j in 0..n {
                    original[i * n + j] +=
                        (sys.j[i * n + k] - sys.r[i * n + k]) * split.quadratic.matrix[k * n + j];
                }
            }
        }
        let balanced = product(&product(&forward, &original, n), &inverse, n);
        let effort = product(&split.quadratic.matrix, &inverse, n);
        let loss_map = product(&resistance_root, &effort, n);
        let (rate, root, doublings) = source_integral_counted(&balanced, &loss_map, n, dt)?;
        Ok(Self {
            n,
            m,
            d,
            forward,
            resistance_root,
            rate,
            root,
            doublings,
        })
    }

    /// Doublings of the plain f64 chain behind the dissipated heat of one step.
    pub(super) fn doublings(&self) -> u32 {
        self.doublings
    }

    pub(super) fn evaluate(
        &self,
        sys: &PortHamiltonian,
        state: &[f64],
        input: &[f64],
    ) -> Result<f64, PhsError> {
        self.evaluate_into(
            sys,
            state,
            input,
            &mut vec![0.0; self.n],
            &mut vec![0.0; self.n],
            &mut vec![0.0; self.d],
        )
    }

    /// No allocation, exponential, equilibrium solve or expanded heat form.
    pub(super) fn evaluate_into(
        &self,
        sys: &PortHamiltonian,
        state: &[f64],
        input: &[f64],
        effort: &mut [f64],
        derivative: &mut [f64],
        coordinates: &mut [f64],
    ) -> Result<f64, PhsError> {
        self.evaluate_shifted_into(sys, state, input, 0.0, effort, derivative, coordinates)
    }

    // Differentiate the SAME retained positive factor, at fixed initial state.
    // Only the original force-imbalance coordinates depend on the held input:
    // dz/du_p = [0; forward * G[:,p] / rate]. Do not expand L^T L or obtain this
    // derivative from two heats; both erase weak losses near static support.
    // `gradient` is private candidate scratch and may change on refusal.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn input_gradient_into_controlled(
        &self,
        sys: &PortHamiltonian,
        state: &[f64],
        input: &[f64],
        effort: &mut [f64],
        derivative: &mut [f64],
        coordinates: &mut [f64],
        gradient: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64, PreparedStepError> {
        let mut poll = || {
            if cancelled() {
                Err(PreparedStepError::Cancelled)
            } else {
                Ok(())
            }
        };
        poll()?;
        if gradient.len() != self.m {
            return Err(refuse("positive heat derivative requires every input column").into());
        }
        let heat = self.evaluate_into(sys, state, input, effort, derivative, coordinates)?;
        poll()?;
        if self.root.is_empty() {
            gradient.fill(0.0);
            return Ok(heat);
        }
        let n = self.n;
        for (port, value) in gradient.iter_mut().enumerate() {
            poll()?;
            for (i, delta) in derivative.iter_mut().enumerate() {
                *delta = (0..n)
                    .map(|j| self.forward[i * n + j] * sys.g[j * self.m + port])
                    .sum::<f64>()
                    / self.rate;
            }
            if !finite(derivative) {
                return Err(refuse("positive heat input direction is unrepresentable").into());
            }
            let mut contraction = Sum::default();
            for row in self.root.chunks_exact(self.d) {
                poll()?;
                let loss = row
                    .iter()
                    .zip(coordinates.iter())
                    .map(|(a, z)| a * z)
                    .sum::<f64>();
                let delta = row[n..]
                    .iter()
                    .zip(derivative.iter())
                    .map(|(a, z)| a * z)
                    .sum::<f64>();
                contraction.product(loss, delta).ok_or_else(|| {
                    refuse("positive heat derivative contraction is unrepresentable")
                })?;
            }
            *value = 2.0
                * contraction
                    .value()
                    .ok_or_else(|| refuse("positive heat input derivative is unrepresentable"))?;
            if !value.is_finite() {
                return Err(refuse("positive heat input derivative is unrepresentable").into());
            }
        }
        poll()?;
        Ok(heat)
    }

    // Integrate the same positive square after subtracting one HELD scalar
    // from its first factor row. Its derivative is unchanged. This preserves
    // nearly matched outgoing wave energy without subtracting large energies.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn evaluate_shifted_into(
        &self,
        sys: &PortHamiltonian,
        state: &[f64],
        input: &[f64],
        first_shift: f64,
        effort: &mut [f64],
        derivative: &mut [f64],
        coordinates: &mut [f64],
    ) -> Result<f64, PhsError> {
        let n = self.n;
        if state.len() != n
            || input.len() != self.m
            || effort.len() != n
            || derivative.len() != n
            || coordinates.len() != self.d
            || !finite(state)
            || !finite(input)
            || !first_shift.is_finite()
        {
            return Err(refuse(
                "positive heat requires complete finite state/input/scratch",
            ));
        }
        if self.root.is_empty() {
            if first_shift != 0.0 {
                return Err(refuse(
                    "constant shifted loss requires explicit interval power",
                ));
            }
            return Ok(0.0);
        }
        sys.storage.gradient(state, effort);
        self.evaluate_effort_shifted_into(sys, effort, input, first_shift, derivative, coordinates)
    }

    /// Integrate the original positive resistance square for
    /// `x' = (J-R)(Q x + held_effort) + G u`, with Q from preparation.
    /// The caller supplies `Q*x_initial + held_effort`, not grad(H) at the
    /// initial endpoint. For a held discrete gradient of the nonlinear
    /// remainder, the initial derivative still evolves as f'=((J-R)Q)f;
    /// the prepared loss factor and its complete physical nullspace therefore
    /// remain valid. No heat is obtained by subtracting storage or port work.
    ///
    /// This path does not evaluate storage, allocate scratch, or alter any
    /// factor. Inputs and both scratch extents must match the prepared owner.
    pub(super) fn evaluate_with_effort_into(
        &self,
        sys: &PortHamiltonian,
        initial_effort: &[f64],
        input: &[f64],
        derivative: &mut [f64],
        coordinates: &mut [f64],
    ) -> Result<f64, PhsError> {
        self.evaluate_effort_shifted_into(sys, initial_effort, input, 0.0, derivative, coordinates)
    }

    fn evaluate_effort_shifted_into(
        &self,
        sys: &PortHamiltonian,
        effort: &[f64],
        input: &[f64],
        first_shift: f64,
        derivative: &mut [f64],
        coordinates: &mut [f64],
    ) -> Result<f64, PhsError> {
        let n = self.n;
        if sys.n != n
            || sys.m != self.m
            || effort.len() != n
            || input.len() != self.m
            || derivative.len() != n
            || coordinates.len() != self.d
            || !finite(effort)
            || !finite(input)
            || !first_shift.is_finite()
        {
            return Err(refuse(
                "positive heat requires complete finite initial effort/input/scratch",
            ));
        }
        if self.root.is_empty() {
            if first_shift != 0.0 {
                return Err(refuse(
                    "constant shifted loss requires explicit interval power",
                ));
            }
            return Ok(0.0);
        }
        let root_rate = det::sqrt(self.rate);
        sys.vector_field_from_effort_into(effort, input, derivative)?;
        for i in 0..n {
            let loss = (0..n)
                .map(|j| self.resistance_root[i * n + j] * effort[j])
                .sum::<f64>();
            coordinates[i] = (loss - if i == 0 { first_shift } else { 0.0 }) / root_rate;
        }
        for i in 0..n {
            coordinates[n + i] = (0..n)
                .map(|j| self.forward[i * n + j] * derivative[j])
                .sum::<f64>()
                / self.rate;
        }
        if !finite(coordinates) {
            return Err(refuse("positive heat initial effort/rate is nonfinite"));
        }
        let mut value = 0.0;
        for row in self.root.chunks_exact(self.d) {
            let loss: f64 = row.iter().zip(coordinates.iter()).map(|(a, b)| a * b).sum();
            value += loss * loss;
        }
        if !value.is_finite() {
            return Err(refuse("physical resistance heat is unrepresentable"));
        }
        Ok(value)
    }
}

/// Common positive integral for an unchanged linear generator and original
/// squared output. For held input, f=x' obeys f'=A f even when A is singular.
/// With a=S*x/sqrt(rate), b=f/rate, a'=sqrt(rate)*S*b and b'=A*b.
/// The initial source/internal-force cancellation occurs before squaring.
pub(super) fn source_integral(
    flow: &[f64], loss_map: &[f64], n: usize, dt: f64,
) -> Result<(f64, Vec<f64>), PhsError> {
    let (rate, root, _) = source_integral_counted(flow, loss_map, n, dt)?;
    Ok((rate, root))
}

/// `source_integral` with the doublings of its chain, for the step whose energy balance reads them
/// (`PositiveHeat`); the source-loss cache has no balance to gate and keeps the two-value form.
pub(super) fn source_integral_counted(
    flow: &[f64], loss_map: &[f64], n: usize, dt: f64,
) -> Result<(f64, Vec<f64>, u32), PhsError> {
    let d = 2 * n;
    let mut rate = 1.0_f64;
    for i in 0..n {
        rate = rate
            .max((0..n).map(|j| flow[i * n + j].abs()).sum())
            .max((0..n).map(|j| flow[j * n + i].abs()).sum());
    }
    if !rate.is_finite() {
        return Err(refuse("positive heat coordinate rate is nonfinite"));
    }
    let root_rate = det::sqrt(rate);
    let mut generator = vec![0.0; d * d];
    let mut factor = vec![0.0; d * d];
    for i in 0..n {
        factor[i * d + i] = root_rate;
        for j in 0..n {
            generator[i * d + n + j] = root_rate * loss_map[i * n + j];
            generator[(n + i) * d + n + j] = flow[i * n + j];
        }
    }
    Ok((
        rate,
        integral_root(&generator, &factor, d, dt)?,
        scaling_doublings(&generator, d, dt)?,
    ))
}

// Exact zero entries separate independent resistance components. A general
// full-C material block must not force a Foster graph's nullspace through an
// eigensolve, or set the absolute scale of a physically weak thermal block.
pub(crate) fn factor_resistance(resistance: &[f64], n: usize) -> Result<Vec<f64>, PhsError> {
    let mut root = vec![0.0; n * n];
    let mut seen = vec![false; n];
    for start in 0..n {
        if seen[start] {
            continue;
        }
        let mut indices = vec![start];
        seen[start] = true;
        let mut next = 0;
        while next < indices.len() {
            let i = indices[next];
            next += 1;
            for j in 0..n {
                if !seen[j] && (resistance[i * n + j] != 0.0 || resistance[j * n + i] != 0.0) {
                    indices.push(j);
                    seen[j] = true;
                }
            }
        }
        indices.sort_unstable();
        let count = indices.len();
        let matrix: Vec<_> = indices
            .iter()
            .flat_map(|&i| indices.iter().map(move |&j| resistance[i * n + j]))
            .collect();
        let factor = if let Some(factor) = graph_root(&matrix, count)? {
            factor
        } else {
            fs_la::eigen::admit_jacobi_eigh(count)
                .map_err(|_| refuse("physical resistance factor exceeds Jacobi work capacity"))?;
            let scale = matrix.iter().fold(0.0_f64, |s, r| s.max(r.abs()));
            let normalized: Vec<_> = matrix.iter().map(|r| r / scale).collect();
            let (values, vectors) = checked_eigh(&normalized, count)?;
            if values.iter().any(|&value| value < 0.0) {
                return Err(refuse(
                    "physical heat resistance has an unresolved negative eigenvalue",
                ));
            }
            let mut factor = vec![0.0; count * count];
            for i in 0..count {
                let weight = det::sqrt(values[i]) * det::sqrt(scale);
                for j in 0..count {
                    factor[i * count + j] = weight * vectors[j * count + i];
                }
            }
            factor
        };
        check_resistance_root(&matrix, &factor, count)?;
        for (i, &row) in indices.iter().enumerate() {
            for (j, &column) in indices.iter().enumerate() {
                root[row * n + column] = factor[i * count + j];
            }
        }
    }
    Ok(root)
}

// Exact signed diagonal dominance supplies physical edge/ground squares:
// |R_ij| (e_i + sign(R_ij) e_j)^2 and (R_ii-sum_j |R_ij|) e_i^2.
// Every positive coefficient survives at its own absolute scale. No tolerance
// selects rank, rescues a negative remainder or changes the original R.
fn graph_root(resistance: &[f64], n: usize) -> Result<Option<Vec<f64>>, PhsError> {
    let mut grounds = vec![0.0; n];
    let mut count = 0;
    for i in 0..n {
        for j in 0..n {
            if resistance[i * n + j] != resistance[j * n + i] {
                return Ok(None);
            }
            if i != j {
                if i < j && resistance[i * n + j] != 0.0 {
                    count += 1;
                }
            }
        }
        let Some(ground) = nonnegative_ground(resistance, n, i) else {
            return Ok(None);
        };
        grounds[i] = ground;
        if grounds[i] > 0.0 {
            count += 1;
        }
    }
    // Forests and diagonal damping already fit in n rows. Denser passive
    // graphs use the existing QR owner to retain every edge in n coordinates.
    let mut rows = vec![0.0; count.max(n) * n];
    let mut row = 0;
    for i in 0..n {
        if grounds[i] > 0.0 {
            rows[row * n + i] = det::sqrt(grounds[i]);
            row += 1;
        }
        for j in i + 1..n {
            let value = resistance[i * n + j];
            if value != 0.0 {
                let weight = det::sqrt(value.abs());
                rows[row * n + i] = weight;
                rows[row * n + j] = if value < 0.0 { -weight } else { weight };
                row += 1;
            }
        }
    }
    let root = if count > n {
        compress(&rows, count, n)?
    } else {
        rows
    };
    Ok(Some(root))
}

// Grow an exact floating-point expansion of R_ii-sum_j |R_ij|. A rounded
// off-diagonal sum can otherwise conceal a negative ground below one ulp of
// the diagonal. Nonoverlapping expansion components are ordered by magnitude,
// so the last nonzero component supplies the exact represented residual sign.
fn nonnegative_ground(resistance: &[f64], n: usize, i: usize) -> Option<f64> {
    let mut expansion = vec![resistance[i * n + i]];
    for j in 0..n {
        if i == j || resistance[i * n + j] == 0.0 {
            continue;
        }
        let mut value = -resistance[i * n + j].abs();
        let mut next = Vec::with_capacity(expansion.len() + 1);
        for &part in &expansion {
            let (sum, error) = two_sum(value, part);
            if !sum.is_finite() || !error.is_finite() {
                return None;
            }
            if error != 0.0 {
                next.push(error);
            }
            value = sum;
        }
        if value != 0.0 {
            next.push(value);
        }
        expansion = next;
    }
    if expansion.last().is_some_and(|&value| value < 0.0) {
        return None;
    }
    let ground: f64 = expansion.iter().sum();
    (ground.is_finite() && ground >= 0.0).then_some(ground)
}

/// Constant of the normwise allowance for a factor of a matrix outside passive-graph form, in
/// units of n eps max|R_kl| (see `check_resistance_root`): the measured deviation is 0.017 of it.
const GENERAL_ROOT_BACKWARD_ERROR: f64 = 1.0;

/// Whether `graph_root` would decompose this matrix: symmetric, with every ground nonnegative.
fn passive_graph_form(resistance: &[f64], n: usize) -> bool {
    for i in 0..n {
        for j in 0..n {
            if resistance[i * n + j] != resistance[j * n + i] {
                return false;
            }
        }
        if nonnegative_ground(resistance, n, i).is_none() {
            return false;
        }
    }
    true
}

// Both graph and eigen factors must resolve each represented nonzero entry;
// a large unrelated diagonal or cancellation inside QR cannot hide its loss.
// At structural zeros only, use the local product sum as the roundoff scale.
pub(crate) fn check_resistance_root(resistance: &[f64], root: &[f64], n: usize) -> Result<(), PhsError> {
    // A matrix in passive-graph form (symmetric, every ground R_ii - sum_j |R_ij| nonnegative) has
    // an exact decomposition into edge and ground squares (`graph_root`) and keeps the
    // componentwise bound below. A symmetric matrix with a nonnegative diagonal (necessary for a
    // resistance) that is not in that form is factored by a symmetric eigensolve, which is
    // backward stable in the NORM: the computed eigenpairs are exact for R + E with
    // ||E||_2 <= p(n) u ||R||_2 (Golub and Van Loan, Matrix Computations, 4th ed., sections 8.3 and
    // 8.5 for symmetric QR and Jacobi, as stated to this change: the pages were not read here). A
    // componentwise bound asks of such a factor an entry it cannot give: at a structural zero of
    // the n = 17 laminar resistance (R_11 = 4.5e5) the Gram entry was 2.887e-11 off against an
    // allowance of 7.7e-12. For that class the allowance is raised to
    // GENERAL_ROOT_BACKWARD_ERROR * n * eps * max|R_kl|, 59 times that measured deviation (1.7e-9).
    // It never makes the bound stricter, and an entry the class cannot resolve below that level is
    // no longer certified there (a weak coupling under n eps max|R| in an eigen-factored component).
    // A negative diagonal or an asymmetry stays under the componentwise bound that refuses it.
    let symmetric = (0..n).all(|i| (0..n).all(|j| resistance[i * n + j] == resistance[j * n + i]));
    let general = symmetric
        && (0..n).all(|i| resistance[i * n + i] >= 0.0)
        && !passive_graph_form(resistance, n);
    let largest = resistance.iter().fold(0.0_f64, |s, r| s.max(r.abs()));
    for i in 0..n {
        for j in i..n {
            let mut gram = 0.0;
            let mut products = 0.0;
            for k in 0..n {
                let product = root[k * n + i] * root[k * n + j];
                gram += product;
                products += product.abs();
            }
            let magnitude = if resistance[i * n + j] != 0.0 {
                resistance[i * n + j].abs()
            } else {
                products
            };
            let mut allowance = 128.0 * n as f64 * f64::EPSILON * magnitude;
            if general {
                allowance =
                    allowance.max(GENERAL_ROOT_BACKWARD_ERROR * n as f64 * f64::EPSILON * largest);
            }
            if !gram.is_finite()
                || !magnitude.is_finite()
                || (gram - resistance[i * n + j]).abs() > allowance
            {
                return Err(refuse(
                    "physical resistance factor exceeds componentwise roundoff",
                ));
            }
        }
    }
    Ok(())
}

// QR's factors retain exact zero rows/columns without asking for numerical
// rank. A common scale protects its squared column norms from tiny absolute
// loss magnitudes; it does not discard small directions relative to others.
pub(crate) fn compress(rows: &[f64], height: usize, n: usize) -> Result<Vec<f64>, PhsError> {
    Ok(compress_paired(rows, None, height, n)?.0)
}

// Identical QR rotations may also act on a second observation. Keeping the
// two factors preserves a signed cross integral without expanding products
// that would lose small work under destructive source interference.
fn compress_paired(
    rows: &[f64], paired: Option<&[f64]>, height: usize, n: usize,
) -> Result<(Vec<f64>, Option<Vec<f64>>), PhsError> {
    let scale = rows.iter().fold(0.0_f64, |s, x| s.max(x.abs()));
    if !scale.is_finite() || paired.is_some_and(|p| p.len() != height * n || !finite(p)) {
        return Err(refuse("positive heat factor is nonfinite"));
    }
    if scale == 0.0 {
        return Ok((vec![0.0; n * n], paired.map(|_| vec![0.0; n * n])));
    }
    let normalized: Vec<_> = rows.iter().map(|x| x / scale).collect();
    let decomposition = qr(&normalized, height, n);
    let mut root = vec![0.0; n * n];
    for i in 0..n {
        for j in i..n {
            root[i * n + j] = scale * decomposition.r(i, j);
        }
    }
    if !finite(&root) {
        return Err(refuse("positive heat QR factor is nonfinite"));
    }
    let other = if let Some(paired) = paired {
        let mut result = vec![0.0; n * n];
        let mut column = vec![0.0; height];
        for j in 0..n {
            for i in 0..height {
                column[i] = paired[i * n + j];
            }
            decomposition.apply_qt(&mut column);
            for i in 0..n {
                result[i * n + j] = column[i];
            }
        }
        if !finite(&result) {
            return Err(refuse("paired power factor is nonfinite"));
        }
        Some(result)
    } else {
        None
    };
    Ok((root, other))
}

pub(super) fn integral_root(
    generator: &[f64],
    factor: &[f64],
    n: usize,
    dt: f64,
) -> Result<Vec<f64>, PhsError> {
    Ok(integral_factors(generator, factor, None, n, dt)?.0)
}

pub(super) fn integral_cross_root(
    generator: &[f64], left: &[f64], right: &[f64], n: usize, dt: f64,
) -> Result<(Vec<f64>, Vec<f64>), PhsError> {
    let (left, right) = integral_factors(generator, left, Some(right), n, dt)?;
    Ok((left, right.expect("paired integral retains its second factor")))
}

/// Doublings of the positive-heat Gramian chain for `generator` over `dt`: the scaling loop of
/// `integral_factors`, kept in one place so that the energy-balance gate reads the count the
/// chain used. The chain is plain f64 and its rounding grows as 2^s eps (see `HEAT_CHAIN_ROUNDING`
/// in linear_exact.rs).
pub(super) fn scaling_doublings(generator: &[f64], n: usize, dt: f64) -> Result<u32, PhsError> {
    let mut norm = 0.0_f64;
    for i in 0..n {
        let row: f64 = (0..n).map(|j| (dt * generator[i * n + j]).abs()).sum();
        let column: f64 = (0..n).map(|j| (dt * generator[j * n + i]).abs()).sum();
        norm = norm.max(row).max(column);
    }
    if !norm.is_finite() {
        return Err(refuse("positive heat operator is nonfinite"));
    }
    let mut scale = 1.0;
    let mut doublings = 0_u32;
    while norm * scale > 0.25 {
        if doublings == 48 {
            return Err(refuse("positive heat scaling exceeds 48 doublings"));
        }
        scale *= 0.5;
        doublings += 1;
    }
    Ok(doublings)
}

fn integral_factors(
    generator: &[f64], factor: &[f64], paired: Option<&[f64]>, n: usize, dt: f64,
) -> Result<(Vec<f64>, Option<Vec<f64>>), PhsError> {
    if !finite(factor) {
        return Err(refuse("positive heat operator is nonfinite"));
    }
    let doublings = scaling_doublings(generator, n, dt)?;
    let mut scale = 1.0;
    for _ in 0..doublings {
        scale *= 0.5;
    }
    let h = dt * scale;
    let a: Vec<_> = generator.iter().map(|x| h * x).collect();
    let mut transition = exponential(&a, n)?;
    // Eight-point Gauss on ||A h|| <= 1/4: the first omitted degree is 16.
    // The exponential product's truncation is below f64 roundoff at this norm;
    // quadrature weights remain positive, unlike an expanded Lyapunov series.
    const NODES: [f64; 8] = [
        -0.9602898564975363,
        -0.7966664774136267,
        -0.525532409916329,
        -0.1834346424956498,
        0.1834346424956498,
        0.525532409916329,
        0.7966664774136267,
        0.9602898564975363,
    ];
    const WEIGHTS: [f64; 8] = [
        0.1012285362903763,
        0.2223810344533745,
        0.3137066458778873,
        0.362683783378362,
        0.362683783378362,
        0.3137066458778873,
        0.2223810344533745,
        0.1012285362903763,
    ];
    let mut samples = vec![0.0; 8 * n * n];
    let mut paired_samples = paired.map(|_| vec![0.0; 8 * n * n]);
    for k in 0..8 {
        let node = 0.5 * (1.0 + NODES[k]);
        let at: Vec<_> = a.iter().map(|x| node * x).collect();
        let map = exponential(&at, n)?;
        let output = product(factor, &map, n);
        let weight = det::sqrt(h) * det::sqrt(0.5 * WEIGHTS[k]);
        for (target, value) in samples[k * n * n..(k + 1) * n * n].iter_mut().zip(output) {
            *target = weight * value;
        }
        if let (Some(factor), Some(samples)) = (paired, paired_samples.as_mut()) {
            let output = product(factor, &map, n);
            for (target, value) in samples[k * n * n..(k + 1) * n * n].iter_mut().zip(output) {
                *target = weight * value;
            }
        }
    }
    let (mut root, mut other) = compress_paired(&samples, paired_samples.as_deref(), 8 * n, n)?;
    for _ in 0..doublings {
        // W(2h)=W(h)+E(h)^T W(h) E(h): QR of [L; L E], never a
        // Gram matrix or exp(-A h). Stiff passive relaxation stays bounded.
        let advanced = product(&root, &transition, n);
        let mut stacked = Vec::with_capacity(2 * n * n);
        stacked.extend_from_slice(&root);
        stacked.extend_from_slice(&advanced);
        let paired_stacked = other.as_ref().map(|other| {
            let advanced = product(other, &transition, n);
            let mut stacked = Vec::with_capacity(2 * n * n);
            stacked.extend_from_slice(other);
            stacked.extend_from_slice(&advanced);
            stacked
        });
        (root, other) = compress_paired(&stacked, paired_stacked.as_deref(), 2 * n, n)?;
        transition = product(&transition, &transition, n);
        if !finite(&transition) {
            return Err(refuse("positive heat transition overflows"));
        }
    }
    Ok((root, other))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{QuadraticStorage, linear_exact::{LinearSplitBudget, split_quadratic}};

    #[test]
    fn a_matrix_outside_graph_form_is_held_to_the_normwise_bound_of_its_eigen_factor() {
        // Ground 1 - 2 < 0: not a passive graph, so the factor comes from the eigensolve.
        let (r, n) = (vec![1.0, 2.0, 2.0, 5.0], 2);
        assert!(!passive_graph_form(&r, n));
        let root = factor_resistance(&r, n).unwrap();
        check_resistance_root(&r, &root, n).unwrap();
        // One factor entry moved by 1e-6 relative is still refused.
        for k in 0..n * n {
            if root[k].abs() > 1e-3 {
                let mut bent = root.clone();
                bent[k] *= 1.0 + 1e-6;
                assert!(check_resistance_root(&r, &bent, n).is_err(), "entry {k}");
            }
        }
    }

    #[test]
    fn a_negative_diagonal_stays_under_the_componentwise_bound() {
        // A large unrelated entry cannot conceal unsupported negative tiny heat: 1e-200 is far below
        // n eps max|R|, so the normwise allowance must not reach this matrix.
        let (r, n) = (vec![1.0, 0.0, 0.0, -1e-200], 2);
        assert!(check_resistance_root(&r, &[1.0, 0.0, 0.0, 0.0], n).is_err());
    }

    #[test]
    fn a_passive_graph_keeps_the_componentwise_bound() {
        let (r, n) = (vec![3.0, -1.0, -1.0, 3.0], 2);
        assert!(passive_graph_form(&r, n));
        let root = factor_resistance(&r, n).unwrap();
        check_resistance_root(&r, &root, n).unwrap();
        let mut bent = root.clone();
        bent[0] *= 1.0 + 1e-6;
        assert!(check_resistance_root(&r, &bent, n).is_err());
    }

    #[test]
    fn held_effort_offset_integrates_the_original_resistance_square() {
        // x'=-r*(q*x+fbar)+u, e=q*x+fbar. The physical effort follows
        // e(t)=u/r+(e0-u/r)*exp(-r*q*t), regardless of the storage gradient
        // at x0. Integrate this analytical square, not an energy remainder.
        let (q, r, dt, input, effort0) = (2.0_f64, 0.3_f64, 0.7_f64, 0.5_f64, -2.0_f64);
        let sys = PortHamiltonian::new(1, 1, vec![0.0], vec![r], vec![1.0],
            Box::new(QuadraticStorage::new(vec![q], 1).unwrap())).unwrap();
        let split = split_quadratic(&sys, LinearSplitBudget {
            max_modes: 1, max_omega_dt: core::f64::consts::PI,
        }).unwrap();
        let heat = PositiveHeat::new(&sys, &split, dt, &[r]).unwrap();
        let (mut derivative, mut coordinates) = ([0.0], [0.0; 2]);
        let actual = heat.evaluate_with_effort_into(&sys, &[effort0], &[input],
            &mut derivative, &mut coordinates).unwrap();
        let equilibrium = input/r;
        let deviation = effort0-equilibrium;
        let rate = r*q;
        let first = -(-rate*dt).exp_m1()/rate;
        let second = -(-2.0*rate*dt).exp_m1()/(2.0*rate);
        let expected = r*(equilibrium*equilibrium*dt
            + 2.0*equilibrium*deviation*first + deviation*deviation*second);
        assert!(actual >= 0.0);
        assert!((actual-expected).abs() < 2e-12*expected);
        assert_eq!(actual.to_bits(), heat.evaluate(&sys, &[effort0/q], &[input]).unwrap().to_bits());
        assert!((actual-heat.evaluate(&sys, &[1.0], &[input]).unwrap()).abs() > 0.1*actual,
            "endpoint storage gradient must not replace the held nonlinear effort");
        assert_eq!(heat.evaluate_with_effort_into(&sys, &[0.0], &[0.0],
            &mut derivative, &mut coordinates).unwrap(), 0.0);
        for bad in [f64::NAN, f64::INFINITY] {
            assert!(heat.evaluate_with_effort_into(&sys, &[bad], &[input],
                &mut derivative, &mut coordinates).is_err());
        }
        assert!(heat.evaluate_with_effort_into(&sys, &[], &[input],
            &mut derivative, &mut coordinates).is_err());
    }
}
