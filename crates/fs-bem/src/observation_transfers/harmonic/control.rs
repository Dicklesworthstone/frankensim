//! Complex source control on the original retained physical BEM field map.
//!
//! This differentiates source coordinates with the solved operator held fixed.
//! It neither differentiates geometry/materials nor invokes reciprocity to
//! replace the Hermitian adjoint. Complex boundary/propagation phases remain
//! those of the original transfer bank, including for nonreciprocal maps.

use super::*;
use fs_math::dd::Dd;

type Parts = [[Dd; 2]; 4];

/// A real objective's complex field covector at one original receiver.
/// Convention: dJ = Re(conj(pressure)*dp + sum conj(velocity_i)*du_i).
/// These are dual coordinates, with units objective/Pa and objective/(m/s),
/// not a physical pressure/velocity phasor. No peak/RMS factor is inferred.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct HarmonicFieldCovector {
    pub pressure: C64,
    pub particle_velocity: [C64; 3],
}

/// Admission for a fixed-frequency acoustic-energy control objective.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicControlLimits {
    pub maximum_quadrature_points: usize,
    /// Conservative complex-product equivalents: R*(8*S+16), with four
    /// components per field and original source/receiver caps. Zero/repeated
    /// volume entries count in geometry admission; zero receiver weights skip
    /// arithmetic but do not reduce this admitted whole-query allowance.
    pub maximum_products: usize,
    /// Logical retained receiver weights/covectors and source-gradient output.
    /// Original fields, targets and quadrature remain borrowed. Allocator
    /// overhead and bounded stack temporaries are excluded.
    pub maximum_workspace_bytes: usize,
}

/// Objective and Cartesian complex source gradient from the same evaluation.
/// For source a=x+i*y, gradient.re=dJ/dx and gradient.im=dJ/dy. Equivalently,
/// dJ=Re(sum conj(gradient_s)*da_s). Entries retain original SourceId order
/// and each source coordinate's physical SI unit; no factor of two is hidden.
pub struct HarmonicEnergyGradient<'a> {
    pub mean_residual_energy_j: f64,
    pub source_gradient: &'a [C64],
}

/// Reusable source objective with actual pressure AND particle velocity.
///
/// For the homogeneous nondispersive resting exterior, the objective is
/// J = sum_r V_r (|p_r-p_target,r|²/(rho*c²)
///                  + rho*|u_r-u_target,r|²) / 4.
/// The factor 1/4 is the cycle-mean energy of peak phasors. None targets a
/// quiet field; supplied targets use the same frequency, receiver frame and
/// common phase reference. An arbitrary target defines an energy-metric
/// mismatch; it does not establish a physically realizable target field.
///
/// Both forward and adjoint applications cost O(sources*receivers). Repeated
/// volume weights coalesce once. Queries and Hessian-vector products allocate
/// nothing, form no dense Hessian and repeat no boundary solve. DD products
/// retain source and target cancellation until the complete residual and
/// gradient have been formed. No relative error certificate near every null,
/// spatial cubature accuracy or geometry/material derivative is inferred.
pub struct PreparedHarmonicEnergyControl<'a> {
    field: &'a PreparedHarmonicObservation,
    quadrature: &'a [HarmonicVolumeQuadrature],
    targets: Option<&'a [HarmonicField]>,
    weights: Vec<Dd>,
    covectors: Vec<Parts>,
    gradient: Vec<C64>,
    products: usize,
    workspace_bytes: usize,
}

impl PreparedHarmonicObservation {
    /// Apply the exact Hermitian adjoint of this retained source-to-field map.
    /// This is a source-coordinate pullback, not a new boundary solve or a
    /// reciprocity assertion. It also applies to the original relaxing-medium
    /// field maps, without inventing an acoustic storage law for that medium.
    ///
    /// Output is in the original source order. Shape, finite inputs and the
    /// complete work budget are admitted before writing. Cancellation or an
    /// unrepresentable product may leave a prefix; discard output on error.
    pub fn source_pullback(
        &self,
        covectors: &[HarmonicFieldCovector],
        budget: HarmonicWorkBudget,
        source_gradient: &mut [C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), TransferError> {
        poll(&mut cancelled)?;
        if covectors.len() != self.receivers.len() || source_gradient.len() != self.sources.len() {
            return Err(bad(
                "harmonic source pullback shape differs from its original field",
            ));
        }
        admit_work(
            self.receivers.len() * self.sources.len() * 4,
            budget.maximum_products,
        )?;
        for covector in covectors {
            poll(&mut cancelled)?;
            if !finite(covector.pressure) || covector.particle_velocity.iter().any(|&v| !finite(v))
            {
                return Err(bad("harmonic field covector must be finite"));
            }
        }
        pullback(
            self,
            |r| {
                let v = covectors[r];
                [
                    v.pressure,
                    v.particle_velocity[0],
                    v.particle_velocity[1],
                    v.particle_velocity[2],
                ]
                .map(parts)
            },
            source_gradient,
            &mut cancelled,
        )
    }
}

impl<'a> PreparedHarmonicEnergyControl<'a> {
    /// Bind the original physical operator, targets and authored nonnegative
    /// volume quadrature. Relaxing/dispersive exterior storage is refused;
    /// equilibrium sound-speed calibration cannot authorize this energy law.
    pub fn new(
        field: &'a PreparedHarmonicObservation,
        quadrature: &'a [HarmonicVolumeQuadrature],
        targets: Option<&'a [HarmonicField]>,
        limits: HarmonicControlLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, TransferError> {
        poll(&mut cancelled)?;
        if !field.radiation_momentum {
            return Err(bad(
                "acoustic energy control requires a nondispersive resting exterior",
            ));
        }
        cap(quadrature.len(), MAX_SURFACE_QUADRATURE_POINTS)?;
        cap(quadrature.len(), limits.maximum_quadrature_points)?;
        let receivers = field.receivers.len();
        let sources = field.sources.len();
        let products = receivers * (8 * sources + 16);
        admit_work(products, limits.maximum_products)?;
        let workspace_bytes = receivers
            .checked_mul(core::mem::size_of::<Dd>() + core::mem::size_of::<Parts>())
            .and_then(|v| v.checked_add(sources.checked_mul(core::mem::size_of::<C64>())?))
            .and_then(|v| v.checked_add(core::mem::size_of::<Self>()))
            .ok_or_else(|| bad("acoustic source control workspace overflows"))?;
        cap(workspace_bytes, limits.maximum_workspace_bytes)?;
        if let Some(targets) = targets {
            if targets.len() != receivers {
                return Err(bad("acoustic target differs from original receiver order"));
            }
            for target in targets {
                poll(&mut cancelled)?;
                if !finite(target.pressure_pa)
                    || target.particle_velocity_m_s.iter().any(|&v| !finite(v))
                {
                    return Err(bad("acoustic control target must be finite"));
                }
            }
        }
        for sample in quadrature {
            poll(&mut cancelled)?;
            if sample.receiver >= receivers
                || !sample.volume_m3.is_finite()
                || sample.volume_m3 < 0.0
            {
                return Err(bad("acoustic control volume quadrature is invalid"));
            }
        }
        let mut weights = allocated(receivers, Dd::ZERO)?;
        for sample in quadrature {
            poll(&mut cancelled)?;
            weights[sample.receiver] = weights[sample.receiver] + Dd::from_f64(sample.volume_m3);
        }
        if weights.iter().any(|&w| !dd_finite(w)) {
            return Err(bad("acoustic control volume is unrepresentable"));
        }
        let covectors = allocated(receivers, [[Dd::ZERO; 2]; 4])?;
        let gradient = allocated(sources, C64::ZERO)?;
        poll(&mut cancelled)?;
        Ok(Self {
            field,
            quadrature,
            targets,
            weights,
            covectors,
            gradient,
            products,
            workspace_bytes,
        })
    }

    pub fn sources(&self) -> &[SourceId] {
        &self.field.sources
    }
    pub fn receivers(&self) -> &[ReceiverSpec] {
        &self.field.receivers
    }
    pub fn quadrature(&self) -> &[HarmonicVolumeQuadrature] {
        self.quadrature
    }
    pub fn angular_frequency_rad_s(&self) -> f64 {
        self.field.omega
    }
    pub fn producer_hash(&self) -> ContentHash {
        self.field.producer_hash
    }
    pub fn producer_uncertainty(&self) -> ProducerUncertainty {
        self.field.producer_uncertainty()
    }
    pub fn work_products(&self) -> usize {
        self.products
    }
    pub fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }

    /// Return a scalar objective and its complete original-source gradient.
    /// On failure no result is published and the objective remains reusable.
    pub fn value_gradient(
        &mut self,
        coefficients: &[C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HarmonicEnergyGradient<'_>, TransferError> {
        let energy = self.evaluate(coefficients, true, true, &mut cancelled)?;
        Ok(HarmonicEnergyGradient {
            mean_residual_energy_j: energy,
            source_gradient: &self.gradient,
        })
    }

    /// Matrix-free Hessian action on a complex Cartesian source direction.
    /// The fixed-operator quadratic has a positive-semidefinite Hessian, with
    /// no target dependence. Source constraints/regularization are not added.
    pub fn hessian_vector(
        &mut self,
        direction: &[C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<&[C64], TransferError> {
        self.evaluate(direction, false, false, &mut cancelled)?;
        Ok(&self.gradient)
    }

    fn evaluate(
        &mut self,
        coefficients: &[C64],
        subtract_target: bool,
        need_value: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, TransferError> {
        poll(cancelled)?;
        self.field.factors(coefficients, 1)?;
        let density = Dd::from_f64(self.field.medium.density);
        let speed = Dd::from_f64(self.field.medium.sound_speed);
        let half = Dd::from_f64(0.5);
        let quarter = Dd::from_f64(0.25);
        let mut energy = Dd::ZERO;
        self.covectors.fill([[Dd::ZERO; 2]; 4]);
        for (receiver, &weight) in self.weights.iter().enumerate() {
            poll(cancelled)?;
            if weight == Dd::ZERO {
                continue;
            }
            // Preserve DD source products through subtraction of the target.
            // Rounding a large forward field first would erase a weak mismatch.
            let mut residual = [[Dd::ZERO; 2]; 4];
            for (source, &coefficient) in coefficients.iter().enumerate() {
                if source % 32 == 0 {
                    poll(cancelled)?;
                }
                let transfer = self.field.transfers[source * self.field.receivers.len() + receiver];
                for component in 0..4 {
                    add_product(
                        &mut residual[component],
                        transfer[component],
                        parts(coefficient),
                        false,
                    );
                }
            }
            if subtract_target && let Some(targets) = self.targets {
                let t = targets[receiver];
                let values = [
                    t.pressure_pa,
                    t.particle_velocity_m_s[0],
                    t.particle_velocity_m_s[1],
                    t.particle_velocity_m_s[2],
                ];
                for c in 0..4 {
                    let target = parts(values[c]);
                    for part in 0..2 {
                        residual[c][part] = residual[c][part] - target[part];
                    }
                }
            }
            if residual.iter().flatten().any(|&v| !dd_finite(v)) {
                return Err(bad("acoustic control residual is unrepresentable"));
            }
            let p_over_c = residual[0].map(|v| v / speed);
            if need_value {
                let mut norm = squared(p_over_c) / density;
                for u in &residual[1..] {
                    norm = norm + density * squared(*u);
                }
                energy = energy + (weight * quarter) * norm;
            }
            for part in 0..2 {
                self.covectors[receiver][0][part] =
                    ((p_over_c[part] / density) / speed) * weight * half;
                for c in 1..4 {
                    self.covectors[receiver][c][part] =
                        (residual[c][part] * density) * weight * half;
                }
            }
            if !dd_finite(energy)
                || self.covectors[receiver]
                    .iter()
                    .flatten()
                    .any(|&v| !dd_finite(v))
            {
                return Err(bad(
                    "acoustic control energy or field derivative is unrepresentable",
                ));
            }
        }
        pullback(
            self.field,
            |r| self.covectors[r],
            &mut self.gradient,
            cancelled,
        )?;
        poll(cancelled)?;
        let value = energy.to_f64();
        if !value.is_finite() || value < 0.0 {
            return Err(bad("acoustic control energy is unrepresentable"));
        }
        Ok(value)
    }
}

fn pullback(
    field: &PreparedHarmonicObservation,
    covector: impl Fn(usize) -> Parts,
    output: &mut [C64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), TransferError> {
    let receivers = field.receivers.len();
    for (source, out) in output.iter_mut().enumerate() {
        let mut gradient = [Dd::ZERO; 2];
        for receiver in 0..receivers {
            poll(cancelled)?;
            let transfer = field.transfers[source * receivers + receiver];
            let v = covector(receiver);
            for c in 0..4 {
                add_product(&mut gradient, transfer[c], v[c], true);
            }
        }
        if gradient.iter().any(|&v| !dd_finite(v)) {
            return Err(bad("harmonic source gradient is unrepresentable"));
        }
        *out = C64::new(gradient[0].to_f64(), gradient[1].to_f64());
        if !finite(*out) {
            return Err(bad("harmonic source gradient is unrepresentable"));
        }
    }
    poll(cancelled)
}

fn add_product(sum: &mut [Dd; 2], a: C64, b: [Dd; 2], conjugate_a: bool) {
    let ar = Dd::from_f64(a.re);
    let ai = Dd::from_f64(if conjugate_a { -a.im } else { a.im });
    sum[0] = sum[0] + ar * b[0] - ai * b[1];
    sum[1] = sum[1] + ar * b[1] + ai * b[0];
}
fn parts(value: C64) -> [Dd; 2] {
    [Dd::from_f64(value.re), Dd::from_f64(value.im)]
}
fn squared(value: [Dd; 2]) -> Dd {
    value[0] * value[0] + value[1] * value[1]
}
fn dd_finite(value: Dd) -> bool {
    value.hi.is_finite() && value.lo.is_finite()
}
fn bad(message: &'static str) -> TransferError {
    TransferError::Invalid(message)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), TransferError> {
    if cancelled() {
        Err(TransferError::Cancelled)
    } else {
        Ok(())
    }
}
fn allocated<T: Clone>(count: usize, value: T) -> Result<Vec<T>, TransferError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| bad("harmonic source control allocation failed"))?;
    values.resize(count, value);
    Ok(values)
}
fn admit_work(products: usize, maximum: usize) -> Result<(), TransferError> {
    cap(maximum, MAX_PRODUCTS)?;
    if maximum == 0 {
        return Err(bad("empty source control work allowance"));
    }
    cap(products, maximum)
}

#[cfg(test)]
mod tests;
