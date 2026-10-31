//! Closed, finite-bath ideal Gibbs--Donnan equilibrium for mobile ions.
//!
//! Every declared ion permeates and reaches electrochemical equilibrium. The
//! two volumes and EACH species' total inventory remain finite. Fixed charge
//! is an amount of elementary-charge equivalents [mol], not an ion reservoir
//! or an osmotic particle count. Independent trapped osmotic particles are
//! supplied separately. No counterions, salt dissociation, or membrane charge
//! are invented; both compartments must admit electroneutrality.
//!
//! With `psi=F*(phi_left-phi_right)/(RT)`, equilibrium gives
//! `n_left_i=N_i*logistic(ln(V_left/V_right)-z_i*psi)`. The left charge residual
//! decreases monotonically with psi. Equal-valence ions share one evaluation
//! during the bounded scalar solve, while their individual amounts remain
//! conserved. Solvent force uses the FULL ionic plus trapped osmolarity.
//!
//! This owner assumes ideal dilute ions with negligible solute partial volume,
//! no association, binding, imposed voltage/current, or double-layer storage.
//! It is a fast-ion equilibrium limit, not transient Nernst--Planck transport.
//! Nonideal solvent pressure must not be combined with this ideal ion closure.
//!
//! Primary finite-bath balance and electrochemical-equilibrium formulation:
//! Waniewski et al. (2021), equations 7--9 and 34--37,
//! <https://doi.org/10.1038/s41598-021-00899-y>. The more general partial-volume
//! and activity corrections intentionally excluded here are described by
//! Waniewski (2023), <https://doi.org/10.1038/s41598-023-47592-w>.

use fs_math::det;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

/// Faraday constant [C/mol], from the exact SI elementary charge and Avogadro constant.
pub const FARADAY_C_PER_MOL: f64 = 1.602_176_634e-19 * 6.022_140_76e23;
/// Bounded ionic layout and stack workspace.
pub const MAX_DONNAN_SPECIES: usize = 64;

/// Physical admission, representability, work, or cancellation failure.
#[derive(Clone, Debug, PartialEq)]
pub enum DonnanError {
    /// Invalid physical input or control.
    Invalid(&'static str),
    /// Caller-owned inventory/output shape differs from the prepared layout.
    Shape { expected: usize, actual: usize },
    /// Finite inventories and fixed charges are not globally electroneutral.
    NotElectroneutral { residual_mol: f64 },
    /// Available finite ions cannot neutralize each compartment's fixed charge.
    InfeasibleCharge,
    /// A required finite physical result is not representable.
    Arithmetic(&'static str),
    /// The bounded scalar solve did not meet the charge tolerance.
    NonConvergence,
    /// Caller cancelled before atomic publication.
    Cancelled,
}
impl core::fmt::Display for DonnanError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Donnan equilibrium: {self:?}")
    }
}
impl std::error::Error for DonnanError {}
type Result<T> = core::result::Result<T, DonnanError>;

/// Explicit charge accuracy and maximum finite-potential residual evaluations.
#[derive(Clone, Copy, Debug)]
pub struct DonnanControls {
    /// Absolute residual allowance [mol of elementary-charge equivalents].
    pub absolute_charge_tolerance_mol: f64,
    /// Relative allowance against total absolute mobile and fixed charge.
    pub relative_charge_tolerance: f64,
    /// Includes initial evaluation, bracket expansion, and safeguarded Newton steps.
    pub maximum_iterations: usize,
}
impl Default for DonnanControls {
    fn default() -> Self {
        Self {
            absolute_charge_tolerance_mol: 0.0,
            relative_charge_tolerance: 4e-13,
            maximum_iterations: 128,
        }
    }
}

/// Finite inventory and volume problem on the prepared ion ordering.
#[derive(Clone, Copy, Debug)]
pub struct DonnanProblem<'a> {
    /// Common absolute temperature [K].
    pub temperature_k: f64,
    /// Positive left and right solvent volumes [m³].
    pub volumes_m3: [f64; 2],
    /// Total moles of each mobile ion in BOTH compartments.
    pub mobile_totals_mol: &'a [f64],
    /// Signed, immobile charge amounts [mol of elementary-charge equivalents].
    pub fixed_charge_mol: [f64; 2],
    /// Independent retained osmotic particles [mol], unrelated to charge valence.
    pub trapped_particles_mol: [f64; 2],
}

/// Potential diagnostic; counterion-only limits need no fictitious salt floor.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum DonnanPotential {
    /// A finite left-minus-right potential difference.
    Finite {
        dimensionless_left_minus_right: f64,
        voltage_left_minus_right_v: f64,
    },
    /// Exact limiting partition requires psi tending to negative infinity.
    NegativeInfinite,
    /// Exact limiting partition requires psi tending to positive infinity.
    PositiveInfinite,
    /// No mobile ions and no fixed charges constrain the potential.
    Indeterminate,
}

/// Accepted finite partition and coherent-SI solvent diagnostics.
#[derive(Clone, Copy, Debug)]
pub struct DonnanReport {
    pub potential: DonnanPotential,
    /// Mobile ionic particle amounts only; each ion counts once.
    pub mobile_osmoles_mol: [f64; 2],
    /// `RT*(mobile_osmoles+trapped_particles)/volume` [Pa].
    pub osmotic_pressure_pa: [f64; 2],
    /// Actual mobile plus fixed charge residual in each compartment [mol].
    pub charge_residual_mol: [f64; 2],
    /// Scalar charge residual evaluations performed.
    pub iterations: usize,
    /// `d(Pi_left-Pi_right)/dx` [Pa/m³], x water transferred left to right.
    /// Nonnegative at fixed total volume/inventories/charges. `None` means
    /// this derivative alone is unrepresentable; the finite partition is valid.
    pub osmotic_pressure_drop_slope_pa_per_m3: Option<f64>,
}
impl DonnanReport {
    /// Full solvent driving pressure, positive for left-to-right water flow.
    /// Both mechanical pressures must use the same gauge. Fixed charge is not
    /// counted a second time as osmotic particles or an ionic reflection factor.
    pub fn solvent_driving_pressure_pa(&self, pressure_pa: [f64; 2]) -> Result<f64> {
        if pressure_pa.iter().any(|p| !p.is_finite()) {
            return Err(DonnanError::Invalid("mechanical pressures must be finite"));
        }
        finite(
            (pressure_pa[0] - pressure_pa[1])
                - (self.osmotic_pressure_pa[0] - self.osmotic_pressure_pa[1]),
            "solvent driving pressure is not representable",
        )
    }
}

/// Prepared ion valences; allocations occur only at preparation.
#[derive(Clone, Debug)]
pub struct IdealDonnan {
    valences: Box<[i32]>,
    group_for_species: Box<[usize]>,
    group_valences: Box<[i32]>,
}
impl IdealDonnan {
    /// Prepare 0--64 explicitly charged species. Zero valence is not an ion;
    /// separately retained neutral osmoles belong in `trapped_particles_mol`.
    pub fn new(valences: Vec<i32>) -> Result<Self> {
        if valences.len() > MAX_DONNAN_SPECIES || valences.contains(&0) {
            return Err(DonnanError::Invalid(
                "Donnan layout needs at most 64 nonzero ion valences",
            ));
        }
        let mut groups = Vec::new();
        let mut indices = Vec::with_capacity(valences.len());
        for &z in &valences {
            let group = groups
                .iter()
                .position(|&existing| existing == z)
                .unwrap_or_else(|| {
                    groups.push(z);
                    groups.len() - 1
                });
            indices.push(group);
        }
        Ok(Self {
            valences: valences.into_boxed_slice(),
            group_for_species: indices.into_boxed_slice(),
            group_valences: groups.into_boxed_slice(),
        })
    }
    pub fn valences(&self) -> &[i32] {
        &self.valences
    }
    /// Immutable preparation allocation, excluding `size_of::<Self>()`.
    pub fn heap_bytes(&self) -> usize {
        core::mem::size_of_val(self.valences.as_ref())
            + core::mem::size_of_val(self.group_for_species.as_ref())
            + core::mem::size_of_val(self.group_valences.as_ref())
    }

    /// Equilibrate finite mobile ions while conserving every species. Failure,
    /// exhausted work, or cancellation leaves BOTH output slices unchanged.
    /// Every solve uses bounded stack buffers and no allocation.
    pub fn equilibrate_into(
        &self,
        problem: DonnanProblem<'_>,
        controls: DonnanControls,
        left_amounts_mol: &mut [f64],
        right_amounts_mol: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DonnanReport> {
        let count = self.valences.len();
        for actual in [
            problem.mobile_totals_mol.len(),
            left_amounts_mol.len(),
            right_amounts_mol.len(),
        ] {
            if actual != count {
                return Err(DonnanError::Shape {
                    expected: count,
                    actual,
                });
            }
        }
        let rt = admit(problem, controls)?;
        poll(&mut cancelled)?;
        let mut grouped = [Accumulator::default(); MAX_DONNAN_SPECIES];
        let mut global = Accumulator::default();
        let mut scale = Accumulator::default();
        let mut minimum = Accumulator::default();
        let mut maximum = Accumulator::default();
        for q in problem.fixed_charge_mol {
            global.add(q);
            scale.add(q.abs());
        }
        minimum.add(problem.fixed_charge_mol[0]);
        maximum.add(problem.fixed_charge_mol[0]);
        for (i, (&z, &n)) in self
            .valences
            .iter()
            .zip(problem.mobile_totals_mol)
            .enumerate()
        {
            let charge = finite(f64::from(z) * n, "mobile charge is not representable")?;
            global.add(charge);
            scale.add(charge.abs());
            if z < 0 {
                minimum.add(charge);
            } else {
                maximum.add(charge);
            }
            grouped[self.group_for_species[i]].add(n);
        }
        let tolerance = finite(
            controls.absolute_charge_tolerance_mol
                + controls.relative_charge_tolerance * scale.finish()?,
            "charge tolerance is not representable",
        )?;
        let global_charge = global.finish()?;
        if global_charge.abs() > tolerance {
            return Err(DonnanError::NotElectroneutral {
                residual_mol: global_charge,
            });
        }
        let (minimum, maximum) = (minimum.finish()?, maximum.finish()?);
        if minimum > 0.0 || maximum < 0.0 {
            return Err(DonnanError::InfeasibleCharge);
        }
        let mut totals = [0.0; MAX_DONNAN_SPECIES];
        for (out, group) in totals
            .iter_mut()
            .zip(grouped)
            .take(self.group_valences.len())
        {
            *out = group.finish()?;
        }
        let log_volume_ratio = det::ln(problem.volumes_m3[0]) - det::ln(problem.volumes_m3[1]);
        let mut iterations = 0;
        let potential = if minimum == 0.0 && maximum == 0.0 {
            DonnanPotential::Indeterminate
        } else if minimum == 0.0 {
            DonnanPotential::PositiveInfinite
        } else if maximum == 0.0 {
            DonnanPotential::NegativeInfinite
        } else {
            let largest_valence = self
                .group_valences
                .iter()
                .zip(&totals)
                .filter(|&(_, &n)| n > 0.0)
                .map(|(&z, _)| f64::from(z).abs())
                .fold(1.0, f64::max);
            let mut evaluate = |psi| {
                poll(&mut cancelled)?;
                if iterations >= controls.maximum_iterations {
                    return Err(DonnanError::NonConvergence);
                }
                iterations += 1;
                charge_at(
                    &self.group_valences,
                    &totals,
                    log_volume_ratio,
                    problem.fixed_charge_mol[0],
                    psi,
                )
            };
            let psi = solve_potential(0.25 * tolerance, 1.0 / largest_valence, &mut evaluate)?;
            DonnanPotential::Finite {
                dimensionless_left_minus_right: psi,
                voltage_left_minus_right_v: finite(
                    psi * (rt / FARADAY_C_PER_MOL),
                    "Donnan voltage is not representable",
                )?,
            }
        };
        let mut left = [0.0; MAX_DONNAN_SPECIES];
        let mut right = [0.0; MAX_DONNAN_SPECIES];
        let mut amount = [Accumulator::default(); 2];
        let mut charge = [Accumulator::default(); 2];
        charge[0].add(problem.fixed_charge_mol[0]);
        charge[1].add(problem.fixed_charge_mol[1]);
        for (i, (&z, &n)) in self
            .valences
            .iter()
            .zip(problem.mobile_totals_mol)
            .enumerate()
        {
            (left[i], right[i]) = partition(n, z, log_volume_ratio, potential);
            for (side, value) in [left[i], right[i]].into_iter().enumerate() {
                amount[side].add(value);
                charge[side].add(f64::from(z) * value);
            }
        }
        let mobile_osmoles_mol = [amount[0].finish()?, amount[1].finish()?];
        let charge_residual_mol = [charge[0].finish()?, charge[1].finish()?];
        if charge_residual_mol.iter().any(|q| q.abs() > tolerance) {
            return Err(DonnanError::NonConvergence);
        }
        let mut osmotic_pressure_pa = [0.0; 2];
        for side in 0..2 {
            let particles = finite(
                mobile_osmoles_mol[side] + problem.trapped_particles_mol[side],
                "total osmotic particles are not representable",
            )?;
            let concentration = finite(
                particles / problem.volumes_m3[side],
                "Donnan osmolarity is not representable",
            )?;
            osmotic_pressure_pa[side] = finite(
                rt * concentration,
                "Donnan osmotic pressure is not representable",
            )?;
        }
        let report = DonnanReport {
            potential,
            mobile_osmoles_mol,
            osmotic_pressure_pa,
            charge_residual_mol,
            iterations,
            osmotic_pressure_drop_slope_pa_per_m3: pressure_slope(
                &self.group_valences,
                &totals,
                problem,
                potential,
                log_volume_ratio,
                rt,
            ),
        };
        poll(&mut cancelled)?;
        left_amounts_mol.copy_from_slice(&left[..count]);
        right_amounts_mol.copy_from_slice(&right[..count]);
        Ok(report)
    }
}

fn admit(problem: DonnanProblem<'_>, controls: DonnanControls) -> Result<f64> {
    if !positive(problem.temperature_k)
        || problem.volumes_m3.iter().any(|&v| !positive(v))
        || problem
            .mobile_totals_mol
            .iter()
            .chain(&problem.trapped_particles_mol)
            .any(|&n| !nonnegative(n))
        || problem.fixed_charge_mol.iter().any(|q| !q.is_finite())
    {
        return Err(DonnanError::Invalid(
            "Donnan inputs require finite temperature, positive volumes, nonnegative moles and finite charges",
        ));
    }
    if !nonnegative(controls.absolute_charge_tolerance_mol)
        || !positive(controls.relative_charge_tolerance)
        || controls.relative_charge_tolerance >= 1.0
        || controls.maximum_iterations == 0
        || controls.maximum_iterations > 4096
    {
        return Err(DonnanError::Invalid(
            "invalid Donnan charge tolerance or iteration budget",
        ));
    }
    finite(R * problem.temperature_k, "Donnan RT is not representable")
}

#[derive(Clone, Copy, Default)]
struct Accumulator {
    value: f64,
    correction: f64,
}
impl Accumulator {
    fn add(&mut self, value: f64) {
        let next = self.value + value;
        self.correction += if self.value.abs() >= value.abs() {
            (self.value - next) + value
        } else {
            (value - next) + self.value
        };
        self.value = next;
    }
    fn finish(self) -> Result<f64> {
        finite(
            self.value + self.correction,
            "Donnan sum is not representable",
        )
    }
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn nonnegative(x: f64) -> bool {
    x.is_finite() && x >= 0.0
}
fn finite(x: f64, reason: &'static str) -> Result<f64> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(DonnanError::Arithmetic(reason))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(DonnanError::Cancelled)
    } else {
        Ok(())
    }
}

// Preserve the small minority amount, then debit that SAME amount from the
// majority. A logarithmic fallback recovers a representable minority even if
// its fraction alone would underflow. No concentration or salt floor is used.
fn split_amount(total: f64, log_odds: f64) -> (f64, f64) {
    if total == 0.0 {
        return (0.0, 0.0);
    }
    let attenuation = det::exp(-log_odds.abs());
    let mut minority = total * (attenuation / (1.0 + attenuation));
    if minority == 0.0 {
        minority = det::exp(det::ln(total) - log_odds.abs() - det::ln(1.0 + attenuation));
    }
    if log_odds >= 0.0 {
        (total - minority, minority)
    } else {
        (minority, total - minority)
    }
}
fn partition(total: f64, z: i32, log_volume_ratio: f64, potential: DonnanPotential) -> (f64, f64) {
    match potential {
        DonnanPotential::Finite {
            dimensionless_left_minus_right: psi,
            ..
        } => split_amount(total, log_volume_ratio - f64::from(z) * psi),
        DonnanPotential::PositiveInfinite => {
            if z < 0 {
                (total, 0.0)
            } else {
                (0.0, total)
            }
        }
        DonnanPotential::NegativeInfinite => {
            if z > 0 {
                (total, 0.0)
            } else {
                (0.0, total)
            }
        }
        DonnanPotential::Indeterminate => (0.0, total),
    }
}

#[derive(Clone, Copy)]
struct ChargePoint {
    psi: f64,
    residual: f64,
    negative_derivative: f64,
}
fn charge_at(
    valences: &[i32],
    totals: &[f64],
    log_volume_ratio: f64,
    fixed_left: f64,
    psi: f64,
) -> Result<ChargePoint> {
    let mut charge = Accumulator::default();
    charge.add(fixed_left);
    let mut slope = 0.0;
    for (&z, &n) in valences.iter().zip(totals) {
        if n == 0.0 {
            continue;
        }
        let (left, right) = split_amount(n, log_volume_ratio - f64::from(z) * psi);
        charge.add(f64::from(z) * left);
        let weight = left.min(right) * (left.max(right) / n);
        slope += f64::from(z) * f64::from(z) * weight;
    }
    Ok(ChargePoint {
        psi,
        residual: charge.finish()?,
        negative_derivative: slope,
    })
}
fn solve_potential(
    tolerance: f64,
    mut distance: f64,
    evaluate: &mut impl FnMut(f64) -> Result<ChargePoint>,
) -> Result<f64> {
    let zero = evaluate(0.0)?;
    if zero.residual.abs() <= tolerance {
        return Ok(0.0);
    }
    let direction = if zero.residual > 0.0 { 1.0 } else { -1.0 };
    let mut previous = zero;
    let (mut lo, mut hi) = loop {
        let point = evaluate(direction * distance)?;
        if point.residual.abs() <= tolerance {
            return Ok(point.psi);
        }
        if direction * point.residual < 0.0 {
            break if direction > 0.0 {
                (previous, point)
            } else {
                (point, previous)
            };
        }
        previous = point;
        distance *= 2.0;
        if !distance.is_finite() {
            return Err(DonnanError::NonConvergence);
        }
    };
    let mut bisect_next = false;
    loop {
        let width = hi.psi - lo.psi;
        let best = if lo.residual.abs() < hi.residual.abs() {
            lo
        } else {
            hi
        };
        let proposed = best.psi + best.residual / best.negative_derivative;
        let interpolate =
            !bisect_next && proposed.is_finite() && proposed > lo.psi && proposed < hi.psi;
        let psi = if interpolate {
            proposed
        } else {
            lo.psi + 0.5 * width
        };
        if psi == lo.psi || psi == hi.psi {
            return Err(DonnanError::NonConvergence);
        }
        let point = evaluate(psi)?;
        if point.residual.abs() <= tolerance {
            return Ok(psi);
        }
        if point.residual > 0.0 {
            lo = point;
        } else {
            hi = point;
        }
        bisect_next = interpolate && hi.psi - lo.psi > 0.5 * width;
    }
}

fn pressure_slope(
    valences: &[i32],
    totals: &[f64],
    problem: DonnanProblem<'_>,
    potential: DonnanPotential,
    log_volume_ratio: f64,
    rt: f64,
) -> Option<f64> {
    let [vl, vr] = problem.volumes_m3;
    let mut curvature = 0.0;
    let mut first = Accumulator::default();
    let mut second = 0.0;
    // A common valence scale cancels from S1²/S2. It prevents squared
    // multivalent charges from overflowing an otherwise finite derivative.
    let charge_scale = valences
        .iter()
        .map(|&z| f64::from(z).abs())
        .fold(1.0, f64::max);
    for (&z, &n) in valences.iter().zip(totals) {
        if n == 0.0 {
            continue;
        }
        let (left, right) = partition(n, z, log_volume_ratio, potential);
        let difference = left / vl - right / vr;
        curvature += square_over(difference.abs(), n);
        let weight = left.min(right) * (left.max(right) / n);
        let charge = f64::from(z) / charge_scale;
        first.add(charge * weight);
        second += charge * charge * weight;
    }
    let first = first.finish().ok()?;
    if second > 0.0 && first != 0.0 {
        let inverse_volume = 1.0 / vl + 1.0 / vr;
        let term = square_over(first.abs(), second) * inverse_volume * inverse_volume;
        curvature += if term.is_normal() {
            term
        } else {
            let log_left = -det::ln(vl);
            let log_right = -det::ln(vr);
            let maximum = log_left.max(log_right);
            let log_inverse =
                maximum + det::ln(det::exp(log_left - maximum) + det::exp(log_right - maximum));
            det::exp(2.0 * (det::ln(first.abs()) + log_inverse) - det::ln(second))
        };
    }
    for (trapped, volume) in problem
        .trapped_particles_mol
        .into_iter()
        .zip(problem.volumes_m3)
    {
        if trapped > 0.0 {
            curvature += (trapped / volume) / volume;
        }
    }
    let slope = rt * curvature;
    if nonnegative(slope) {
        Some(slope)
    } else {
        None
    }
}

fn square_over(value: f64, denominator: f64) -> f64 {
    if value == 0.0 {
        return 0.0;
    }
    let ratio = value / denominator;
    let result = value * ratio;
    if ratio.is_normal() && result.is_normal() {
        result
    } else {
        det::exp(2.0 * det::ln(value) - det::ln(denominator))
    }
}

#[cfg(test)]
mod tests;
