//! Stationary index-one radiation forces on the original transport cells.
//!
//! Every actual collision transfers photon momentum, including elastic
//! scattering. The native phase sampler supplies the outgoing direction; no
//! averaged asymmetry or absorption-only pressure replaces it. This is the
//! force on stationary matter in vacuum/index-one transport. Refractive-host
//! stress, moving-medium Doppler work, thermal emission recoil and polarization
//! are not inferred by these scalar intensity histories. The optional wrench
//! query also resolves orbital torque at the native interaction positions;
//! polarization-dependent spin angular momentum remains outside this model.

use super::*;
use crate::volumes::transport::deposition::BeamPowerDeposition;
use crate::volumes::transport::scattering::field::CellScatteringGrid;

mod orbital;
pub use orbital::{
    AngularMomentumRateEstimate, CellRadiationTorque, RadiationAngularMomentumBalance,
    ScatteringWrenchDeposition,
};
use orbital::{HistoryOrbital, OrbitalTallies};

const LIGHT_M_S: f64 = 299_792_458.0;

/// Cartesian momentum per unit time, kg m/s² = N, in world coordinates.
/// Errors are fixed-sample standard errors of the signed component means;
/// they are not confidence bounds or a covariance matrix.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MomentumRateEstimate<const N: usize> {
    pub mean_n: [[f64; 3]; N],
    pub standard_error_n: [[f64; 3]; N],
}

/// Radiation force on one original X-major grid cell. Repeated visits are
/// combined within each history before estimating its sample variance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CellRadiationForce<const N: usize> {
    pub cell_index: usize,
    pub force: MomentumRateEstimate<N>,
}

/// Incident momentum rate equals escape + unresolved + material force.
/// The total force statistics retain correlations between all visited cells.
/// Unresolved radiation remains a directed photon packet, never hidden force
/// or heat. Its later interactions require increasing the collision budget.
#[derive(Clone, Debug, PartialEq)]
pub struct RadiationMomentumBalance<const N: usize> {
    pub incident_momentum_rate_n: [[f64; 3]; N],
    pub ballistic_escape: MomentumRateEstimate<N>,
    pub diffuse_escape: MomentumRateEstimate<N>,
    pub unresolved: MomentumRateEstimate<N>,
    pub force: MomentumRateEstimate<N>,
    /// Compensated incident-minus-outgoing-minus-force residual, per lane.
    pub balance_residual_n: [[f64; 3]; N],
    /// Ascending original cell indices; pure elastic scatterers are retained.
    pub cells: Vec<CellRadiationForce<N>>,
}

/// The native energy and momentum accounts from exactly the same histories.
/// Applied force may be integrated by a mechanical owner, but this stationary
/// optical solve does not debit finite-velocity work or advance a body state.
#[derive(Clone, Debug, PartialEq)]
pub struct ScatteringMomentumDeposition<const N: usize> {
    pub power: ScatteringPowerDeposition<N>,
    pub momentum: RadiationMomentumBalance<N>,
}

impl ScatteringGrid<'_> {
    /// Trace the original energy histories and their per-cell radiation force.
    ///
    /// Absorption-only lanes are exact Beer transport. Scattering uses the
    /// actual wavelength/cell phase sample. The final capped collision also
    /// samples its outgoing direction, so its recoil and unresolved momentum
    /// close together. Those final draws cannot affect another independent
    /// history; the returned power report matches the heating-only query.
    /// Force storage is bounded by `max_deposition_cells`. Cancellation or any
    /// refusal returns neither account; the medium and its sources are unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn deposit_multiple_scattering_momentum<const N: usize>(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        unit_direction: Vec3,
        wavelengths_nm: &[f64; N],
        incident_power_w: &[f64; N],
        config: ScatteringDepositionConfig,
    ) -> Result<ScatteringMomentumDeposition<N>, VolumeTransportError> {
        let mut tallies = MomentumTallies::new(config.max_deposition_cells);
        let power = self.deposit_observed(
            cx,
            origin,
            unit_direction,
            wavelengths_nm,
            incident_power_w,
            config,
            Some(&mut tallies),
        )?;
        let momentum = tallies.finish(&mut || checkpoint(cx), &power, unit_direction)?;
        Ok(ScatteringMomentumDeposition { power, momentum })
    }
}

impl<const N: usize> CellScatteringGrid<'_, '_, N> {
    /// The same joint energy/momentum query with this field's original local
    /// coefficients, phase law, wavelengths and immutable source identity.
    pub fn deposit_multiple_scattering_momentum(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        unit_direction: Vec3,
        incident_power_w: &[f64; N],
        config: ScatteringDepositionConfig,
    ) -> Result<ScatteringMomentumDeposition<N>, VolumeTransportError> {
        self.deposit_multiple_scattering_momentum_with_checkpoint(
            &mut || checkpoint(cx),
            origin,
            unit_direction,
            incident_power_w,
            config,
        )
    }

    /// Hosting-solver checkpoint variant; no partial power or force report is
    /// published when the original callback, work or spatial budget refuses.
    pub fn deposit_multiple_scattering_momentum_with_checkpoint(
        &self,
        check: &mut (impl FnMut() -> Result<(), VolumeTransportError> + ?Sized),
        origin: Point3,
        unit_direction: Vec3,
        incident_power_w: &[f64; N],
        config: ScatteringDepositionConfig,
    ) -> Result<ScatteringMomentumDeposition<N>, VolumeTransportError> {
        let mut tallies = MomentumTallies::new(config.max_deposition_cells);
        let power = deposit_prepared_observed(
            self,
            check,
            origin,
            unit_direction,
            incident_power_w,
            config,
            Some(&mut tallies),
        )?;
        let momentum = tallies.finish(check, &power, unit_direction)?;
        Ok(ScatteringMomentumDeposition { power, momentum })
    }
}

fn photon_rate(power: f64) -> Result<f64, VolumeTransportError> {
    let rate = power / LIGHT_M_S;
    if !rate.is_finite() || power > 0.0 && rate == 0.0 {
        return Err(VolumeTransportError::Unrepresentable);
    }
    Ok(rate)
}

fn components(direction: Vec3) -> [f64; 3] {
    [direction.x, direction.y, direction.z]
}

fn estimate<const N: usize>(
    moments: [[Moments; 3]; N],
) -> Result<MomentumRateEstimate<N>, VolumeTransportError> {
    let mut result = MomentumRateEstimate {
        mean_n: [[0.0; 3]; N],
        standard_error_n: [[0.0; 3]; N],
    };
    for (lane, values) in moments.iter().enumerate() {
        for (axis, value) in values.iter().enumerate() {
            (
                result.mean_n[lane][axis],
                result.standard_error_n[lane][axis],
            ) = value.estimate()?;
        }
    }
    Ok(result)
}

pub(super) struct MomentumTallies<const N: usize> {
    force: [[Moments; 3]; N],
    escape: [[[Moments; 3]; N]; 3],
    spatial: SpatialTallies<N, 3>,
    orbital: Option<OrbitalTallies<N>>,
}
impl<const N: usize> MomentumTallies<N> {
    fn new(limit: usize) -> Self {
        Self {
            force: [[Moments::default(); 3]; N],
            escape: [[[Moments::default(); 3]; N]; 3],
            spatial: SpatialTallies::new(limit),
            orbital: None,
        }
    }
    pub(super) fn exact(
        &mut self,
        lane: usize,
        origin: Point3,
        direction: Vec3,
        power: &BeamPowerDeposition<1>,
    ) -> Result<(), VolumeTransportError> {
        let mut history = HistoryMomentum::new(self, lane, 0);
        for cell in &power.cells {
            history.capture(
                cell.cell_index,
                origin,
                cell.absorbed_power_w[0],
                0.0,
                direction,
                direction,
            )?;
        }
        history.leave(0, power.transmitted_power_w[0], origin, direction)?;
        history.finish()
    }
    fn finish(
        self,
        check: &mut (impl FnMut() -> Result<(), VolumeTransportError> + ?Sized),
        power: &ScatteringPowerDeposition<N>,
        direction: Vec3,
    ) -> Result<RadiationMomentumBalance<N>, VolumeTransportError> {
        check()?;
        let force = estimate(self.force)?;
        let [ballistic_escape, diffuse_escape, unresolved] = self.escape.map(estimate);
        let (ballistic_escape, diffuse_escape, unresolved) =
            (ballistic_escape?, diffuse_escape?, unresolved?);
        let mut incident = [[0.0; 3]; N];
        let mut residual = [[0.0; 3]; N];
        for lane in 0..N {
            let rate = photon_rate(power.incident_power_w[lane])?;
            incident[lane] = components(direction).map(|d| d * rate);
            for axis in 0..3 {
                let mut balance = PowerSum::default();
                for value in [
                    incident[lane][axis],
                    -ballistic_escape.mean_n[lane][axis],
                    -diffuse_escape.mean_n[lane][axis],
                    -unresolved.mean_n[lane][axis],
                    -force.mean_n[lane][axis],
                ] {
                    balance.add(value)?;
                }
                residual[lane][axis] = balance.value();
            }
        }
        let cells =
            self.spatial
                .finish_map(check, &power.sampled_histories, |cell_index, values| {
                    CellRadiationForce {
                        cell_index,
                        force: MomentumRateEstimate {
                            mean_n: values.map(|lane| lane.map(|value| value.0)),
                            standard_error_n: values.map(|lane| lane.map(|value| value.1)),
                        },
                    }
                })?;
        check()?;
        Ok(RadiationMomentumBalance {
            incident_momentum_rate_n: incident,
            ballistic_escape,
            diffuse_escape,
            unresolved,
            force,
            balance_residual_n: residual,
            cells,
        })
    }
}

pub(super) struct HistoryMomentum<'a, const N: usize> {
    tallies: &'a mut MomentumTallies<N>,
    lane: usize,
    history: u32,
    force: [PowerSum; 3],
    escape: [[f64; 3]; 3],
    orbital: Option<HistoryOrbital>,
}
impl<'a, const N: usize> HistoryMomentum<'a, N> {
    pub(super) fn new(tallies: &'a mut MomentumTallies<N>, lane: usize, history: u32) -> Self {
        let orbital = tallies
            .orbital
            .as_ref()
            .map(|tally| HistoryOrbital::new(tally.reference));
        Self {
            tallies,
            lane,
            history,
            force: [PowerSum::default(); 3],
            escape: [[0.0; 3]; 3],
            orbital,
        }
    }
    pub(super) fn capture(
        &mut self,
        cell: usize,
        position: Point3,
        heat: f64,
        continuation: f64,
        incoming: Vec3,
        outgoing: Vec3,
    ) -> Result<(), VolumeTransportError> {
        let absorbed = photon_rate(heat)?;
        let scattered = photon_rate(continuation)?;
        let a = components(incoming);
        let b = components(outgoing);
        // P_abs*d_in + P_continue*(d_in-d_out) retains weak absorption and
        // small forward recoil without subtracting two large photon fluxes.
        let force: [f64; 3] =
            std::array::from_fn(|axis| scattered.mul_add(a[axis] - b[axis], absorbed * a[axis]));
        if force != [0.0; 3] {
            let statistics = &mut self.tallies.spatial.cell(cell)?.lanes[self.lane];
            for axis in 0..3 {
                self.force[axis].add(force[axis])?;
                statistics[axis].capture(self.history, force[axis])?;
            }
            if let Some(orbital) = &mut self.orbital {
                orbital.capture(
                    self.tallies
                        .orbital
                        .as_mut()
                        .ok_or(VolumeTransportError::Unrepresentable)?,
                    (self.lane, self.history, cell),
                    position,
                    force,
                )?;
            }
        }
        Ok(())
    }
    pub(super) fn leave(
        &mut self,
        channel: usize,
        power: f64,
        position: Point3,
        direction: Vec3,
    ) -> Result<(), VolumeTransportError> {
        let rate = photon_rate(power)?;
        self.escape[channel] = components(direction).map(|d| rate * d);
        if let Some(orbital) = &mut self.orbital {
            orbital.leave(channel, position, self.escape[channel])?;
        }
        Ok(())
    }
    pub(super) fn finish(self) -> Result<(), VolumeTransportError> {
        for axis in 0..3 {
            self.tallies.force[self.lane][axis].observe_signed(self.force[axis].value())?;
            for channel in 0..3 {
                self.tallies.escape[channel][self.lane][axis]
                    .observe_signed(self.escape[channel][axis])?;
            }
        }
        if let Some(orbital) = self.orbital {
            orbital.finish(
                self.tallies
                    .orbital
                    .as_mut()
                    .ok_or(VolumeTransportError::Unrepresentable)?,
                self.lane,
            )?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
