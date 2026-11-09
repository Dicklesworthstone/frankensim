//! Reversible shared-core storage on the ORIGINAL winding flux linkages.
//!
//! For signed turns N, positive physical leakage L_l and source MMF F(Phi),
//! lambda_i = L_li i_i + N_i Phi and F(Phi) = sum N_i i_i. Eliminating currents
//! gives F(Phi) + A Phi = s, A=sum N_i^2/L_li, s=sum N_i lambda_i/L_li.
//! The original piecewise-affine source therefore needs only one prepared
//! scalar inverse, not a second nonlinear solver. H=integral F dPhi plus
//! sum L_li i_i^2/2. Its gradient is current, including mutual saturation.
//!
//! Equivalently H=lambda^T D lambda/2 - V(s), V'=Phi. The SAME scalar integral
//! supplies exact chord work and its analytic Newton action. Energy queries use
//! the positive core/leakage decomposition, avoiding cancellation of the two
//! large terms in that equivalent expression. No extrapolation, invented
//! leakage, hysteresis, eddy loss, motion or thermal law is inferred.
//!
//! Energy/flux conventions: MIT 6.685 Class Notes 1,
//! <https://web.mit.edu/6.685/www/chapter1.pdf>.

use std::sync::Arc;

use fs_math::dd::Dd;

use super::super::scalar::ScalarCurve;
use super::{CircuitError, curve_error};

/// Bounds preparation and each allocation-free native constitutive callback.
pub const MAX_SHARED_CORE_WINDINGS: usize = 64;

fn invalid(what: &'static str) -> CircuitError {
    CircuitError::Invalid { what }
}

fn finite(value: f64) -> Result<f64, CircuitError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(
            "shared-core magnetic arithmetic is not representable",
        ))
    }
}

fn dd(value: f64) -> Dd {
    Dd::from_f64(value)
}

/// One immutable, bounded common-flux path and its physical winding leakage.
///
/// Original source knots are (common flux Phi [Wb], magnetomotive force [A
/// turns]); their integral is core-plus-declared-gap energy [J]. Both polarities
/// and an explicit (0,0) knot are required. Turns include winding orientation.
/// Leakage is strictly positive and supplied by the caller: the singular
/// perfect-coupling limit requires constrained coordinates and is not repaired.
/// Immutable tables are shared across clones and repeated circuit assemblies.
#[derive(Clone, Debug, PartialEq)]
pub struct SharedCoreInductance {
    turns: Arc<[f64]>,
    leakage: Arc<[f64]>,
    weights: Arc<[Dd]>,
    reluctance: Dd,
    core: ScalarCurve,
    slopes: Arc<[Dd]>,
    drive_knots: Arc<[Dd]>,
}

#[derive(Clone, Copy)]
struct CorePoint {
    drive: Dd,
    phi: Dd,
    segment: usize,
}

impl SharedCoreInductance {
    pub fn try_new(
        turns: Vec<f64>,
        leakage_inductance_h: Vec<f64>,
        core_flux_mmf_knots: Vec<(f64, f64)>,
    ) -> Result<Self, CircuitError> {
        let n = turns.len();
        if !(1..=MAX_SHARED_CORE_WINDINGS).contains(&n)
            || leakage_inductance_h.len() != n
            || turns.iter().any(|v| !v.is_finite() || *v == 0.0)
            || leakage_inductance_h
                .iter()
                .any(|v| !v.is_finite() || *v <= 0.0 || !v.recip().is_finite() || v.recip() <= 0.0)
        {
            return Err(invalid(
                "shared core requires bounded nonzero signed turns and positive physical leakage",
            ));
        }
        let core = ScalarCurve::try_new(core_flux_mmf_knots).map_err(curve_error)?;
        let mut weights = Vec::with_capacity(n);
        let mut reluctance = Dd::ZERO;
        for (&turn, &leakage) in turns.iter().zip(&leakage_inductance_h) {
            let weight = dd(turn) / dd(leakage);
            let contribution = dd(turn) * weight;
            if finite(weight.to_f64())? == 0.0 || finite(contribution.to_f64())? <= 0.0 {
                return Err(invalid("shared-core winding reduction is unrepresentable"));
            }
            weights.push(weight);
            reluctance = reluctance + contribution;
        }
        finite(reluctance.to_f64())?;
        let pivot_floor = 16.0 * f64::EPSILON * n as f64;
        let mut slopes = Vec::with_capacity(core.knots().len() - 1);
        for pair in core.knots().windows(2) {
            let slope = (dd(pair[1].1) - dd(pair[0].1)) / (dd(pair[1].0) - dd(pair[0].0));
            let total = slope + reluctance;
            finite(total.to_f64())?;
            // If core stiffness disappears beside leakage elimination, a
            // represented inverse would lose its positive common-flux mode.
            if (slope / total).to_f64() <= pivot_floor {
                return Err(invalid(
                    "shared-core leakage/common-flux separation is numerically unresolved",
                ));
            }
            slopes.push(slope);
        }
        let drive_knots = core
            .knots()
            .iter()
            .map(|&(phi, mmf)| {
                let drive = reluctance * dd(phi) + dd(mmf);
                finite(drive.to_f64())?;
                Ok(drive)
            })
            .collect::<Result<Vec<_>, CircuitError>>()?;
        Ok(Self {
            turns: turns.into(),
            leakage: leakage_inductance_h.into(),
            weights: weights.into(),
            reluctance,
            core,
            slopes: slopes.into(),
            drive_knots: drive_knots.into(),
        })
    }

    pub fn turns(&self) -> &[f64] {
        &self.turns
    }
    pub fn leakage_inductance_h(&self) -> &[f64] {
        &self.leakage
    }
    pub fn core_flux_mmf_knots(&self) -> &[(f64, f64)] {
        self.core.knots()
    }
    /// Closed physical common-flux support [Wb]. Individual leakage linkages
    /// are not bounded by this interval: opposing currents can cancel core MMF.
    pub fn common_flux_range_wb(&self) -> (f64, f64) {
        self.core.range()
    }

    fn admit_vector(&self, values: &[f64]) -> Result<(), CircuitError> {
        if values.len() != self.turns.len() || values.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "complete finite shared-core winding vector required",
            ));
        }
        Ok(())
    }

    fn drive(&self, flux: impl Fn(usize) -> f64) -> Result<Dd, CircuitError> {
        let mut sum = Dd::ZERO;
        for (i, &weight) in self.weights.iter().enumerate() {
            sum = sum + weight * dd(finite(flux(i))?);
        }
        finite(sum.to_f64())?;
        Ok(sum)
    }

    fn slope(&self, segment: usize) -> Dd {
        self.slopes[segment]
    }

    fn knot_drive(&self, knot: usize) -> Dd {
        self.drive_knots[knot]
    }

    // Keep Phi's low word through leakage subtraction. Rounding it first can
    // erase magnetizing current long before a coupled matrix is singular.
    fn point(&self, flux: impl Fn(usize) -> f64) -> Result<CorePoint, CircuitError> {
        let drive = self.drive(flux)?;
        let last = self.drive_knots.len() - 1;
        if drive.lt(self.drive_knots[0]) || self.drive_knots[last].lt(drive) {
            return Err(invalid(
                "common flux is outside the shared-core source support",
            ));
        }
        let knots = self.core.knots();
        let segment = self
            .drive_knots
            .partition_point(|&knot| !drive.lt(knot))
            .saturating_sub(1)
            .min(knots.len() - 2);
        let anchor = if drive.to_f64() < 0.0 {
            segment + 1
        } else {
            segment
        };
        let phi = if drive == self.knot_drive(segment) {
            dd(knots[segment].0)
        } else if drive == self.knot_drive(segment + 1) {
            dd(knots[segment + 1].0)
        } else {
            dd(knots[anchor].0)
                + (drive - self.knot_drive(anchor)) / (self.reluctance + self.slope(segment))
        };
        finite(phi.to_f64())?;
        let (lo, hi) = self.core.range();
        if phi.lt(dd(lo)) || dd(hi).lt(phi) {
            return Err(invalid(
                "common flux is outside the shared-core source support",
            ));
        }
        Ok(CorePoint {
            drive,
            phi,
            segment,
        })
    }

    fn current(&self, i: usize, linkage: f64, phi: Dd) -> Result<f64, CircuitError> {
        finite(((dd(linkage) - dd(self.turns[i]) * phi) / dd(self.leakage[i])).to_f64())
    }

    fn core_tail(&self, point: CorePoint) -> Result<Dd, CircuitError> {
        let rounded = point.phi.to_f64();
        let tail = point.phi - dd(rounded);
        Ok(dd(self.core.effort(rounded).map_err(curve_error)?) * tail
            + dd(0.5) * self.slope(point.segment) * tail * tail)
    }

    fn energy(&self, flux: impl Fn(usize) -> f64) -> Result<f64, CircuitError> {
        let point = self.point(&flux)?;
        let mut energy = dd(self.core.energy(point.phi.to_f64()).map_err(curve_error)?)
            + self.core_tail(point)?;
        for i in 0..self.turns.len() {
            let leakage_flux = dd(flux(i)) - dd(self.turns[i]) * point.phi;
            energy = energy + dd(0.5) * leakage_flux * (leakage_flux / dd(self.leakage[i]));
        }
        let energy = finite(energy.to_f64())?;
        if energy < 0.0 || (energy == 0.0 && (0..self.turns.len()).any(|i| flux(i) != 0.0)) {
            return Err(invalid(
                "nonzero shared-core stored energy is unrepresentable",
            ));
        }
        Ok(energy)
    }

    pub fn common_flux_wb(&self, flux_linkages_wb: &[f64]) -> Result<f64, CircuitError> {
        self.admit_vector(flux_linkages_wb)?;
        Ok(self.point(|i| flux_linkages_wb[i])?.phi.to_f64())
    }

    /// Actual endpoint currents, in the original declared winding order.
    pub fn currents_a(&self, flux_linkages_wb: &[f64]) -> Result<Vec<f64>, CircuitError> {
        self.admit_vector(flux_linkages_wb)?;
        let phi = self.point(|i| flux_linkages_wb[i])?.phi;
        (0..self.turns.len())
            .map(|i| self.current(i, flux_linkages_wb[i], phi))
            .collect()
    }

    /// Initialize all original linkages from currents using the SAME source.
    /// Cancellation of signed ampere-turns occurs before core-domain admission.
    /// If rounding the resulting linkages leaves the closed source domain or
    /// makes their energy unrepresentable, initialization refuses.
    pub fn flux_linkages_wb(&self, currents_a: &[f64]) -> Result<Vec<f64>, CircuitError> {
        self.admit_vector(currents_a)?;
        let mut mmf = Dd::ZERO;
        for (&turns, &current) in self.turns.iter().zip(currents_a) {
            mmf = mmf + dd(turns) * dd(current);
        }
        finite(mmf.to_f64())?;
        let knots = self.core.knots();
        if mmf.lt(dd(knots[0].1)) || dd(knots[knots.len() - 1].1).lt(mmf) {
            return Err(invalid(
                "combined winding MMF is outside the shared-core source",
            ));
        }
        let segment = knots
            .partition_point(|&(_, f)| !mmf.lt(dd(f)))
            .saturating_sub(1)
            .min(knots.len() - 2);
        let anchor = if mmf.to_f64() < 0.0 {
            segment + 1
        } else {
            segment
        };
        let phi = if mmf == dd(knots[segment].1) {
            dd(knots[segment].0)
        } else if mmf == dd(knots[segment + 1].1) {
            dd(knots[segment + 1].0)
        } else {
            dd(knots[anchor].0) + (mmf - dd(knots[anchor].1)) / self.slope(segment)
        };
        let linkages = self
            .turns
            .iter()
            .zip(self.leakage.iter())
            .zip(currents_a)
            .map(|((&n, &l), &i)| finite((dd(l) * dd(i) + dd(n) * phi).to_f64()))
            .collect::<Result<Vec<_>, _>>()?;
        self.energy(|i| linkages[i])?;
        Ok(linkages)
    }

    /// Positive core energy plus positive leakage energy [J]. Neither term is
    /// an extra store to add again to the circuit's Hamiltonian.
    pub fn stored_energy_j(&self, flux_linkages_wb: &[f64]) -> Result<f64, CircuitError> {
        self.admit_vector(flux_linkages_wb)?;
        self.energy(|i| flux_linkages_wb[i])
    }

    /// Reciprocal d(lambda)/d(i) [H], row-major in winding order. This is
    /// diag(L_leak)+N N^T/F'(Phi), not a matrix of secant inductances.
    /// Exact source knots select the right segment, except the final endpoint;
    /// a low-word offset beside a knot keeps its actual source segment. Native
    /// continuous small-signal actions separately refuse nonsmooth knots/edges.
    pub fn differential_inductance_h(
        &self,
        flux_linkages_wb: &[f64],
    ) -> Result<Vec<f64>, CircuitError> {
        self.admit_vector(flux_linkages_wb)?;
        let point = self.point(|i| flux_linkages_wb[i])?;
        let permeance = Dd::ONE / self.slope(point.segment);
        let n = self.turns.len();
        let mut matrix = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..=i {
                let mut entry = dd(self.turns[i]) * permeance * dd(self.turns[j]);
                if i == j {
                    entry = entry + dd(self.leakage[i]);
                }
                let entry = finite(entry.to_f64())?;
                matrix[i * n + j] = entry;
                matrix[j * n + i] = entry;
            }
        }
        Ok(matrix)
    }

    fn mean_flux(
        &self,
        a: impl Fn(usize) -> f64,
        b: impl Fn(usize) -> f64,
    ) -> Result<Dd, CircuitError> {
        self.mean_points(self.point(a)?, self.point(b)?)
    }

    // Integrate the original affine inverse spans relative to the lower flux,
    // without subtracting two large inverse-potential values. Low words survive
    // the final leakage subtraction even in a tightly coupled common mode.
    fn mean_points(&self, a: CorePoint, b: CorePoint) -> Result<Dd, CircuitError> {
        if a.segment == b.segment || a.drive == b.drive {
            return Ok(dd(0.5) * a.phi + dd(0.5) * b.phi);
        }
        let (lo, hi) = if a.drive.lt(b.drive) { (a, b) } else { (b, a) };
        let mut position = lo.drive;
        let mut phi = lo.phi;
        let mut relative_integral = Dd::ZERO;
        for segment in lo.segment..=hi.segment {
            let (end, end_phi) = if segment == hi.segment {
                (hi.drive, hi.phi)
            } else {
                (
                    self.knot_drive(segment + 1),
                    dd(self.core.knots()[segment + 1].0),
                )
            };
            relative_integral = relative_integral
                + dd(0.5) * ((phi - lo.phi) + (end_phi - lo.phi)) * (end - position);
            position = end;
            phi = end_phi;
        }
        let mean = lo.phi + relative_integral / (hi.drive - lo.drive);
        finite(mean.to_f64())?;
        Ok(mean)
    }

    fn mean_current(&self, i: usize, a: f64, b: f64, mean_phi: Dd) -> Result<f64, CircuitError> {
        finite(
            ((dd(0.5) * dd(a) + dd(0.5) * dd(b) - dd(self.turns[i]) * mean_phi)
                / dd(self.leakage[i]))
            .to_f64(),
        )
    }

    /// Exact work-conjugate chord currents. These are not generally a temporal
    /// average or the arithmetic mean of the endpoint currents after saturation.
    pub fn discrete_currents_a(
        &self,
        before: &[f64],
        after: &[f64],
    ) -> Result<Vec<f64>, CircuitError> {
        self.admit_vector(before)?;
        self.admit_vector(after)?;
        let phi = self.mean_flux(|i| before[i], |i| after[i])?;
        (0..self.turns.len())
            .map(|i| self.mean_current(i, before[i], after[i], phi))
            .collect()
    }

    // Difference of the ORIGINAL core and leakage energies, independently of
    // discrete currents. Compensated products retain small work beside bias.
    fn energy_change(
        &self,
        a: impl Fn(usize) -> f64,
        b: impl Fn(usize) -> f64,
    ) -> Result<Dd, CircuitError> {
        let pa = self.point(&a)?;
        let pb = self.point(&b)?;
        let mut result = dd(self
            .core
            .energy_change(pa.phi.to_f64(), pb.phi.to_f64())
            .map_err(curve_error)?)
            + self.core_tail(pb)?
            - self.core_tail(pa)?;
        for i in 0..self.turns.len() {
            let da = dd(a(i)) - dd(self.turns[i]) * pa.phi;
            let db = dd(b(i)) - dd(self.turns[i]) * pb.phi;
            result = result + (dd(0.5) * da + dd(0.5) * db) * ((db - da) / dd(self.leakage[i]));
        }
        finite(result.to_f64())?;
        Ok(result)
    }

    pub fn energy_change_j(&self, before: &[f64], after: &[f64]) -> Result<f64, CircuitError> {
        self.admit_vector(before)?;
        self.admit_vector(after)?;
        Ok(self.energy_change(|i| before[i], |i| after[i])?.to_f64())
    }
}

/// Complete common-core storage for distinct ORIGINAL graph inductor branches.
/// Branch order selects the signed turn/leakage order. Native assembly admits
/// ownership and graph types; the law itself contains no graph or time state.
#[derive(Clone, Debug, PartialEq)]
pub struct MagneticStorageBlock {
    pub branches: Vec<usize>,
    pub law: SharedCoreInductance,
}

/// Prepared state-index binding; every callback borrows original state directly.
#[derive(Clone)]
pub(in crate::circuit) struct BoundMagneticStorage {
    pub coordinates: Vec<usize>,
    pub law: SharedCoreInductance,
}

impl BoundMagneticStorage {
    fn admit(&self, state: &[f64]) -> bool {
        self.coordinates
            .iter()
            .all(|&i| state.get(i).is_some_and(|v| v.is_finite()))
    }

    pub(in crate::circuit) fn energy(&self, state: &[f64]) -> f64 {
        if !self.admit(state) {
            return f64::NAN;
        }
        self.law
            .energy(|i| state[self.coordinates[i]])
            .unwrap_or(f64::NAN)
    }

    pub(in crate::circuit) fn gradient(&self, state: &[f64], out: &mut [f64]) -> bool {
        if !self.admit(state) || !self.coordinates.iter().all(|&i| i < out.len()) {
            return false;
        }
        let Ok(point) = self.law.point(|i| state[self.coordinates[i]]) else {
            return false;
        };
        for (i, &coordinate) in self.coordinates.iter().enumerate() {
            let Ok(current) = self.law.current(i, state[coordinate], point.phi) else {
                return false;
            };
            out[coordinate] = current;
        }
        true
    }

    pub(in crate::circuit) fn mean(&self, a: &[f64], b: &[f64], out: &mut [f64]) -> bool {
        if !self.admit(a) || !self.admit(b) || !self.coordinates.iter().all(|&i| i < out.len()) {
            return false;
        }
        let Ok(phi) = self
            .law
            .mean_flux(|i| a[self.coordinates[i]], |i| b[self.coordinates[i]])
        else {
            return false;
        };
        for (i, &coordinate) in self.coordinates.iter().enumerate() {
            let Ok(current) = self.law.mean_current(i, a[coordinate], b[coordinate], phi) else {
                return false;
            };
            out[coordinate] = current;
        }
        true
    }

    pub(in crate::circuit) fn remainder(&self, a: &[f64], b: &[f64], mid: &[f64]) -> f64 {
        if !self.admit(a) || !self.admit(b) || !self.admit(mid) {
            return f64::NAN;
        }
        let Ok(mut change) = self
            .law
            .energy_change(|i| a[self.coordinates[i]], |i| b[self.coordinates[i]])
        else {
            return f64::NAN;
        };
        let Ok(point) = self.law.point(|i| mid[self.coordinates[i]]) else {
            return f64::NAN;
        };
        for (i, &coordinate) in self.coordinates.iter().enumerate() {
            let Ok(current) = self.law.current(i, mid[coordinate], point.phi) else {
                return f64::NAN;
            };
            change = change - dd(current) * (dd(b[coordinate]) - dd(a[coordinate]));
        }
        change.to_f64()
    }

    fn action(&self, direction: &[f64], out: &mut [f64], slope: Dd, diagonal_scale: f64) -> bool {
        if !self.admit(direction) || !self.coordinates.iter().all(|&i| i < out.len()) {
            return false;
        }
        let Ok(drive) = self.law.drive(|i| direction[self.coordinates[i]]) else {
            return false;
        };
        for (i, &coordinate) in self.coordinates.iter().enumerate() {
            let value = ((dd(diagonal_scale) * dd(direction[coordinate])
                - dd(self.law.turns[i]) * slope * drive)
                / dd(self.law.leakage[i]))
            .to_f64();
            if !value.is_finite() {
                return false;
            }
            out[coordinate] = value;
        }
        true
    }

    pub(in crate::circuit) fn continuous_action(
        &self,
        state: &[f64],
        direction: &[f64],
        out: &mut [f64],
    ) -> bool {
        if !self.admit(state) {
            return false;
        }
        let Ok(point) = self.law.point(|i| state[self.coordinates[i]]) else {
            return false;
        };
        let k = point.segment;
        let last = self.law.drive_knots.len() - 1;
        if point.drive == self.law.drive_knots[0]
            || point.drive == self.law.drive_knots[last]
            || (k > 0
                && point.drive == self.law.drive_knots[k]
                && self.law.slopes[k - 1] != self.law.slopes[k])
        {
            return false;
        }
        let slope = Dd::ONE / (self.law.reluctance + self.law.slope(point.segment));
        self.action(direction, out, slope, 1.0)
    }

    pub(in crate::circuit) fn end_action(
        &self,
        a: &[f64],
        b: &[f64],
        direction: &[f64],
        out: &mut [f64],
    ) -> bool {
        if !self.admit(a) || !self.admit(b) {
            return false;
        }
        let (Ok(pa), Ok(pb)) = (
            self.law.point(|i| a[self.coordinates[i]]),
            self.law.point(|i| b[self.coordinates[i]]),
        ) else {
            return false;
        };
        let slope = if pa.segment == pb.segment || pa.drive == pb.drive {
            dd(0.5) / (self.law.reluctance + self.law.slope(pb.segment))
        } else {
            let Ok(mean) = self.law.mean_points(pa, pb) else {
                return false;
            };
            (pb.phi - mean) / (pb.drive - pa.drive)
        };
        self.action(direction, out, slope, 0.5)
    }
}
