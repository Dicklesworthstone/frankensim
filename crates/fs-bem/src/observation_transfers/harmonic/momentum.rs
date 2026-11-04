//! Cycle-mean acoustic momentum through an authored oriented surface.
//!
//! The original pressure and all three particle velocities give the inviscid
//! radiation stress, without a plane-wave or far-field approximation:
//! `S_ij = (|p|²/(4 rho c²) - rho |u|²/4) delta_ij
//!           + rho Re(u_i conj(u_j))/2` for peak exp(-i omega t) phasors.
//! Force-rate and moment-rate fluxes are `sum S dA` and
//! `sum (x - origin) cross (S dA)`. Source columns are combined before products;
//! independent coherence factors are combined afterward.
//!
//! This is the second-order, complete-cycle momentum flux of linear acoustics
//! in a homogeneous, inviscid, nondispersive fluid at rest. With a closed,
//! outward surface in that fluid enclosing all relevant sources/scatterers,
//! its negative is the cycle-mean reaction on the enclosed system. An open
//! patch, incomplete incident field, or selected receivers do not establish
//! such a force balance. No closure, enclosure, spatial convergence, streaming,
//! viscous wall stress, finite-time field storage or moving-medium model is
//! inferred. Torque includes the full near-field velocity, not just radial
//! far-field power. See the stress-integral derivation in
//! <https://doi.org/10.1121/10.0002491>.

use super::*;
use fs_math::dd::Dd;

mod welch;
pub use welch::{HarmonicMomentumWelchLine, PreparedWelchSurfaceMomentum};

/// Signed outgoing momentum rates through the declared oriented quadrature.
/// Signs reverse when every area vector reverses. These are whole-cycle means.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct HarmonicMomentumFlux {
    /// Linear momentum per second [N], in the original receiver frame.
    pub linear_momentum_rate_n: [f64; 3],
    /// Angular momentum per second [N m], about the explicitly supplied origin.
    pub angular_momentum_rate_nm: [f64; 3],
}

/// Borrowed, allocation-free momentum observer of the original BEM fields.
///
/// The geometry and fields stay borrowed and immutable. Repeated evaluations
/// perform no solve, kernel quadrature, covariance assembly or heap allocation.
/// Repeated/zero/signed area entries retain their authored meaning; surface
/// normals and closed geometry are not inferred from area-vector cancellation.
pub struct PreparedHarmonicSurfaceMomentum<'a> {
    field: &'a PreparedHarmonicObservation,
    quadrature: &'a [HarmonicSurfaceQuadrature],
    moment_origin_m: [f64; 3],
}

impl PreparedHarmonicObservation {
    /// Prepare cycle-mean momentum transport on existing physical receivers.
    /// All entries are checked, including silent and zero-area entries.
    pub fn momentum_surface<'a>(
        &'a self,
        quadrature: &'a [HarmonicSurfaceQuadrature],
        moment_origin_m: [f64; 3],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedHarmonicSurfaceMomentum<'a>, TransferError> {
        poll(&mut cancelled)?;
        if !self.radiation_momentum {
            return Err(TransferError::Invalid(
                "acoustic radiation momentum requires a nondispersive exterior",
            ));
        }
        cap(quadrature.len(), MAX_SURFACE_QUADRATURE_POINTS)?;
        if moment_origin_m.iter().any(|x| !x.is_finite()) {
            return Err(TransferError::Invalid(
                "acoustic momentum origin is invalid",
            ));
        }
        for sample in quadrature {
            poll(&mut cancelled)?;
            if sample.receiver >= self.receivers.len()
                || sample.oriented_area_m2.iter().any(|a| !a.is_finite())
                || self.receivers[sample.receiver]
                    .position_m
                    .iter()
                    .any(|x| !x.is_finite())
            {
                return Err(TransferError::Invalid(
                    "acoustic momentum quadrature is invalid",
                ));
            }
        }
        poll(&mut cancelled)?;
        Ok(PreparedHarmonicSurfaceMomentum {
            field: self,
            quadrature,
            moment_origin_m,
        })
    }
}

impl PreparedHarmonicSurfaceMomentum<'_> {
    pub fn angular_frequency_rad_s(&self) -> f64 {
        self.field.omega
    }
    pub fn sources(&self) -> &[SourceId] {
        &self.field.sources
    }
    pub fn physical_receivers(&self) -> &[ReceiverSpec] {
        &self.field.receivers
    }
    pub fn quadrature(&self) -> &[HarmonicSurfaceQuadrature] {
        self.quadrature
    }
    pub fn moment_origin_m(&self) -> [f64; 3] {
        self.moment_origin_m
    }
    pub fn producer_hash(&self) -> ContentHash {
        self.field.producer_hash
    }
    pub fn producer_uncertainty(&self) -> ProducerUncertainty {
        self.field.producer_uncertainty()
    }

    /// One coherent drive, retaining incident/scattered and cross-source terms.
    pub fn evaluate_coherent(
        &self,
        coefficients: &[C64],
        budget: HarmonicWorkBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<HarmonicMomentumFlux, TransferError> {
        self.evaluate_coherence_modes(coefficients, 1, budget, cancelled)
    }

    /// Source peak covariance is B B*, with B stored mode-major/source-minor.
    /// Zero modes are explicit silence. Work is Q M (4 S + 7) complex
    /// transfer/quadratic products, with fixed-size stress/torque contractions.
    /// Only quadrature-selected receivers are evaluated. No candidate escapes
    /// on cancellation, budget refusal or an unrepresentable physical product.
    pub fn evaluate_coherence_modes(
        &self,
        factors: &[C64],
        modes: usize,
        budget: HarmonicWorkBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HarmonicMomentumFlux, TransferError> {
        self.field.factors(factors, modes)?;
        self.work(modes, 4 * self.field.sources.len() + 7, budget)?;
        poll(&mut cancelled)?;
        let mut sum = MomentumSum::default();
        for sample in self.quadrature {
            poll(&mut cancelled)?;
            for factor in factors.chunks_exact(self.field.sources.len()) {
                poll(&mut cancelled)?;
                self.add(
                    &mut sum,
                    sample,
                    self.field.field(factor, sample.receiver)?,
                    1.0,
                )?;
            }
        }
        poll(&mut cancelled)?;
        sum.finish()
    }

    /// Independent source peak second moments in each coordinate's squared SI
    /// units. This uses Q S * 7 quadratic products, without a diagonal factor
    /// matrix or square roots of source powers. Silent sources cost their
    /// admitted work allowance but skip physical multiplication.
    pub fn evaluate_independent_sources(
        &self,
        peak_second_moments: &[f64],
        budget: HarmonicWorkBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HarmonicMomentumFlux, TransferError> {
        self.field.powers(peak_second_moments)?;
        self.work(self.field.sources.len(), 7, budget)?;
        poll(&mut cancelled)?;
        let mut sum = MomentumSum::default();
        for sample in self.quadrature {
            poll(&mut cancelled)?;
            for (source, &power) in peak_second_moments.iter().enumerate() {
                poll(&mut cancelled)?;
                if power != 0.0 {
                    let field =
                        self.field.transfers[source * self.field.receivers.len() + sample.receiver];
                    self.add(&mut sum, sample, field, power)?;
                }
            }
        }
        poll(&mut cancelled)?;
        sum.finish()
    }

    fn work(
        &self,
        modes: usize,
        products: usize,
        budget: HarmonicWorkBudget,
    ) -> Result<(), TransferError> {
        cap(budget.maximum_products, MAX_PRODUCTS)?;
        if budget.maximum_products == 0 {
            return Err(TransferError::Invalid("empty harmonic work allowance"));
        }
        let work = self
            .quadrature
            .len()
            .checked_mul(modes)
            .and_then(|n| n.checked_mul(products))
            .ok_or(TransferError::Invalid("acoustic momentum work overflow"))?;
        cap(work, budget.maximum_products)
    }

    fn add(
        &self,
        sum: &mut MomentumSum,
        sample: &HarmonicSurfaceQuadrature,
        field: [C64; 4],
        power: f64,
    ) -> Result<(), TransferError> {
        let stress = stress(field, self.field.medium, power)?;
        let area = sample.oriented_area_m2.map(Dd::from_f64);
        let position = self.field.receivers[sample.receiver].position_m;
        let arm = core::array::from_fn::<_, 3, _>(|i| {
            Dd::from_f64(position[i]) - Dd::from_f64(self.moment_origin_m[i])
        });
        let traction = core::array::from_fn::<_, 3, _>(|i| {
            (0..3).fold(Dd::ZERO, |v, j| v + stress[i][j] * area[j])
        });
        for i in 0..3 {
            sum.linear[i] = sum.linear[i] + traction[i];
            let j = (i + 1) % 3;
            let k = (i + 2) % 3;
            sum.angular[i] = sum.angular[i] + arm[j] * traction[k] - arm[k] * traction[j];
        }
        if sum
            .linear
            .iter()
            .chain(&sum.angular)
            .any(|&x| !dd_finite(x))
        {
            return Err(TransferError::Invalid(
                "acoustic momentum surface product is unrepresentable",
            ));
        }
        Ok(())
    }
}

#[derive(Default)]
struct MomentumSum {
    linear: [Dd; 3],
    angular: [Dd; 3],
}
impl MomentumSum {
    fn finish(self) -> Result<HarmonicMomentumFlux, TransferError> {
        if self
            .linear
            .iter()
            .chain(&self.angular)
            .any(|&x| !dd_finite(x))
        {
            return Err(TransferError::Invalid(
                "acoustic momentum sum is unrepresentable",
            ));
        }
        Ok(HarmonicMomentumFlux {
            linear_momentum_rate_n: self.linear.map(Dd::to_f64),
            angular_momentum_rate_nm: self.angular.map(Dd::to_f64),
        })
    }
}

// Retain cancellation in the isotropic term and across oriented patches in DD.
// Form pressure/c before squaring: this avoids a needless c² intermediate.
fn stress(field: [C64; 4], medium: Medium, power: f64) -> Result<[[Dd; 3]; 3], TransferError> {
    let density = Dd::from_f64(medium.density);
    let speed = Dd::from_f64(medium.sound_speed);
    let weight = Dd::from_f64(power);
    let half = Dd::from_f64(0.5);
    let quarter = Dd::from_f64(0.25);
    let pressure = [
        Dd::from_f64(field[0].re) / speed,
        Dd::from_f64(field[0].im) / speed,
    ];
    let pressure_energy = (pressure[0] * pressure[0] + pressure[1] * pressure[1]) / density;
    let velocity = [field[1], field[2], field[3]].map(|z| [Dd::from_f64(z.re), Dd::from_f64(z.im)]);
    let mut products = [[Dd::ZERO; 3]; 3];
    let mut kinetic = Dd::ZERO;
    for i in 0..3 {
        for j in i..3 {
            let value =
                density * (velocity[i][0] * velocity[j][0] + velocity[i][1] * velocity[j][1]);
            products[i][j] = value;
            products[j][i] = value;
        }
        kinetic = kinetic + products[i][i];
    }
    let isotropic = (pressure_energy - kinetic) * quarter;
    for (i, row) in products.iter_mut().enumerate() {
        for (j, value) in row.iter_mut().enumerate() {
            *value = (*value * half + if i == j { isotropic } else { Dd::ZERO }) * weight;
            if !dd_finite(*value) {
                return Err(TransferError::Invalid(
                    "acoustic radiation stress is unrepresentable",
                ));
            }
        }
    }
    Ok(products)
}

fn dd_finite(value: Dd) -> bool {
    value.hi.is_finite() && value.lo.is_finite()
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), TransferError> {
    if cancelled() {
        Err(TransferError::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests;
