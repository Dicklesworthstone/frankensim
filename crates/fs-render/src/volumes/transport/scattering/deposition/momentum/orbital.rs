//! Orbital angular momentum from the original scalar photon trajectories.
//!
//! A free flight conserves (r - reference) × p because displacement is parallel
//! to p. A collision transfers (r_collision - reference) × delta_p to matter.
//! This resolves the line of action inside a cell; cell-center force moments
//! are not a substitute. Spin/polarization torque is not modeled.

use super::*;

/// Orbital angular momentum per unit time, kg m²/s² = N m.
/// Component errors are fixed-sample standard errors, not confidence bounds
/// or a covariance matrix; repeated events are combined within each history.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AngularMomentumRateEstimate<const N: usize> {
    pub mean_n_m: [[f64; 3]; N],
    pub standard_error_n_m: [[f64; 3]; N],
}

/// Torque about the report's reference point on an original X-major cell.
/// Its force is available at the same original index in the momentum report.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CellRadiationTorque<const N: usize> {
    pub cell_index: usize,
    pub torque: AngularMomentumRateEstimate<N>,
}

/// Incident orbital momentum rate equals outgoing rate plus material torque.
/// Every quantity uses the same world-coordinate reference point in meters.
/// This is extrinsic orbital angular momentum of scalar ray packets, including
/// off-axis beams, and does not infer optical spin or a phase-resolved vortex.
#[derive(Clone, Debug, PartialEq)]
pub struct RadiationAngularMomentumBalance<const N: usize> {
    pub reference_point: Point3,
    pub incident_angular_momentum_rate_n_m: [[f64; 3]; N],
    pub ballistic_escape: AngularMomentumRateEstimate<N>,
    pub diffuse_escape: AngularMomentumRateEstimate<N>,
    pub unresolved: AngularMomentumRateEstimate<N>,
    pub torque: AngularMomentumRateEstimate<N>,
    pub balance_residual_n_m: [[f64; 3]; N],
    /// Ascending original cell indices with nonzero sampled torque.
    pub cells: Vec<CellRadiationTorque<N>>,
}

/// Energy, linear momentum and orbital angular momentum of the same native
/// trajectories. Forces and torques act on stationary matter. A mechanical
/// owner may integrate them; this query does not advance bodies or debit work.
#[derive(Clone, Debug, PartialEq)]
pub struct ScatteringWrenchDeposition<const N: usize> {
    pub power: ScatteringPowerDeposition<N>,
    pub momentum: RadiationMomentumBalance<N>,
    pub orbital: RadiationAngularMomentumBalance<N>,
}

impl ScatteringGrid<'_> {
    /// Resolve native per-cell force and orbital torque about a finite world
    /// point in meters. Actual collision positions determine torque; pure
    /// absorption retains the exact source line of action through every cell.
    ///
    /// Histories, phase draws and energy estimates match the momentum query.
    /// Escaping and capped unresolved packets retain their orbital flux. The
    /// original spatial/work limits and cancellation apply to the whole result.
    /// No orbital cell storage is allocated by the existing power/force APIs.
    #[allow(clippy::too_many_arguments)]
    pub fn deposit_multiple_scattering_wrench<const N: usize>(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        unit_direction: Vec3,
        wavelengths_nm: &[f64; N],
        incident_power_w: &[f64; N],
        reference_point: Point3,
        config: ScatteringDepositionConfig,
    ) -> Result<ScatteringWrenchDeposition<N>, VolumeTransportError> {
        let mut tallies =
            MomentumTallies::with_orbital(config.max_deposition_cells, reference_point)?;
        let power = self.deposit_observed(
            cx,
            origin,
            unit_direction,
            wavelengths_nm,
            incident_power_w,
            config,
            Some(&mut tallies),
        )?;
        tallies.finish_wrench(&mut || checkpoint(cx), power, origin, unit_direction)
    }
}

impl<const N: usize> CellScatteringGrid<'_, '_, N> {
    /// Joint force/torque query using this field's original local spectra and
    /// phase laws, about the supplied finite world-coordinate point in meters.
    #[allow(clippy::too_many_arguments)]
    pub fn deposit_multiple_scattering_wrench(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        unit_direction: Vec3,
        incident_power_w: &[f64; N],
        reference_point: Point3,
        config: ScatteringDepositionConfig,
    ) -> Result<ScatteringWrenchDeposition<N>, VolumeTransportError> {
        self.deposit_multiple_scattering_wrench_with_checkpoint(
            &mut || checkpoint(cx),
            origin,
            unit_direction,
            incident_power_w,
            reference_point,
            config,
        )
    }

    /// Hosting-solver checkpoint variant. Invalid reference geometry,
    /// unrepresentable moments, budget exhaustion or cancellation publishes
    /// none of the three accounts and leaves the borrowed medium unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn deposit_multiple_scattering_wrench_with_checkpoint(
        &self,
        check: &mut (impl FnMut() -> Result<(), VolumeTransportError> + ?Sized),
        origin: Point3,
        unit_direction: Vec3,
        incident_power_w: &[f64; N],
        reference_point: Point3,
        config: ScatteringDepositionConfig,
    ) -> Result<ScatteringWrenchDeposition<N>, VolumeTransportError> {
        let mut tallies =
            MomentumTallies::with_orbital(config.max_deposition_cells, reference_point)?;
        let power = deposit_prepared_observed(
            self,
            check,
            origin,
            unit_direction,
            incident_power_w,
            config,
            Some(&mut tallies),
        )?;
        tallies.finish_wrench(check, power, origin, unit_direction)
    }
}

impl<const N: usize> MomentumTallies<N> {
    pub(super) fn with_orbital(
        limit: usize,
        reference: Point3,
    ) -> Result<Self, VolumeTransportError> {
        if ![reference.x, reference.y, reference.z]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(VolumeTransportError::Ray);
        }
        let mut tallies = Self::new(limit);
        tallies.orbital = Some(OrbitalTallies {
            reference,
            torque: [[Moments::default(); 3]; N],
            escape: [[[Moments::default(); 3]; N]; 3],
            spatial: SpatialTallies::new(limit),
        });
        Ok(tallies)
    }

    fn finish_wrench(
        mut self,
        check: &mut (impl FnMut() -> Result<(), VolumeTransportError> + ?Sized),
        power: ScatteringPowerDeposition<N>,
        origin: Point3,
        direction: Vec3,
    ) -> Result<ScatteringWrenchDeposition<N>, VolumeTransportError> {
        let orbital = self
            .orbital
            .take()
            .ok_or(VolumeTransportError::Unrepresentable)?;
        let momentum = self.finish(check, &power, direction)?;
        let orbital = orbital.finish(check, &power, origin, &momentum)?;
        Ok(ScatteringWrenchDeposition {
            power,
            momentum,
            orbital,
        })
    }
}

pub(super) struct OrbitalTallies<const N: usize> {
    pub(super) reference: Point3,
    torque: [[Moments; 3]; N],
    escape: [[[Moments; 3]; N]; 3],
    spatial: SpatialTallies<N, 3>,
}

fn angular_estimate<const N: usize>(
    moments: [[Moments; 3]; N],
) -> Result<AngularMomentumRateEstimate<N>, VolumeTransportError> {
    let value = estimate(moments)?;
    Ok(AngularMomentumRateEstimate {
        mean_n_m: value.mean_n,
        standard_error_n_m: value.standard_error_n,
    })
}

impl<const N: usize> OrbitalTallies<N> {
    fn finish(
        self,
        check: &mut (impl FnMut() -> Result<(), VolumeTransportError> + ?Sized),
        power: &ScatteringPowerDeposition<N>,
        origin: Point3,
        momentum: &RadiationMomentumBalance<N>,
    ) -> Result<RadiationAngularMomentumBalance<N>, VolumeTransportError> {
        check()?;
        let torque = angular_estimate(self.torque)?;
        let [ballistic_escape, diffuse_escape, unresolved] = self.escape.map(angular_estimate);
        let (ballistic_escape, diffuse_escape, unresolved) =
            (ballistic_escape?, diffuse_escape?, unresolved?);
        let mut incident = [[0.0; 3]; N];
        let mut residual = [[0.0; 3]; N];
        for lane in 0..N {
            incident[lane] = moment_about(
                self.reference,
                origin,
                momentum.incident_momentum_rate_n[lane],
            )?;
            for axis in 0..3 {
                let mut balance = PowerSum::default();
                for value in [
                    incident[lane][axis],
                    -ballistic_escape.mean_n_m[lane][axis],
                    -diffuse_escape.mean_n_m[lane][axis],
                    -unresolved.mean_n_m[lane][axis],
                    -torque.mean_n_m[lane][axis],
                ] {
                    balance.add(value)?;
                }
                residual[lane][axis] = balance.value();
            }
        }
        let cells =
            self.spatial
                .finish_map(check, &power.sampled_histories, |cell_index, values| {
                    CellRadiationTorque {
                        cell_index,
                        torque: AngularMomentumRateEstimate {
                            mean_n_m: values.map(|lane| lane.map(|value| value.0)),
                            standard_error_n_m: values.map(|lane| lane.map(|value| value.1)),
                        },
                    }
                })?;
        check()?;
        Ok(RadiationAngularMomentumBalance {
            reference_point: self.reference,
            incident_angular_momentum_rate_n_m: incident,
            ballistic_escape,
            diffuse_escape,
            unresolved,
            torque,
            balance_residual_n_m: residual,
            cells,
        })
    }
}

/// Compensated two-product difference for each cross-product component. The
/// FMA residual retains a small physical moment when two large products nearly
/// cancel. Nonfinite lever arms or products refuse instead of publishing NaNs.
pub(super) fn moment_about(
    reference: Point3,
    position: Point3,
    rate: [f64; 3],
) -> Result<[f64; 3], VolumeTransportError> {
    if rate == [0.0; 3] {
        return Ok([0.0; 3]);
    }
    let r = [
        position.x - reference.x,
        position.y - reference.y,
        position.z - reference.z,
    ];
    if !r.into_iter().all(f64::is_finite) {
        return Err(VolumeTransportError::Unrepresentable);
    }
    let result = std::array::from_fn(|axis| {
        let j = (axis + 1) % 3;
        let k = (axis + 2) % 3;
        let product = r[k] * rate[j];
        let correction = r[k].mul_add(rate[j], -product);
        r[j].mul_add(rate[k], -product) - correction
    });
    if !result.into_iter().all(f64::is_finite) {
        return Err(VolumeTransportError::Unrepresentable);
    }
    Ok(result)
}

pub(super) struct HistoryOrbital {
    reference: Point3,
    pub(super) torque: [PowerSum; 3],
    pub(super) escape: [[f64; 3]; 3],
}
impl HistoryOrbital {
    pub(super) fn new(reference: Point3) -> Self {
        Self {
            reference,
            torque: [PowerSum::default(); 3],
            escape: [[0.0; 3]; 3],
        }
    }
    pub(super) fn capture<const N: usize>(
        &mut self,
        tallies: &mut OrbitalTallies<N>,
        (lane, history, cell): (usize, u32, usize),
        position: Point3,
        force: [f64; 3],
    ) -> Result<(), VolumeTransportError> {
        let torque = moment_about(self.reference, position, force)?;
        if torque != [0.0; 3] {
            let statistics = &mut tallies.spatial.cell(cell)?.lanes[lane];
            for axis in 0..3 {
                self.torque[axis].add(torque[axis])?;
                statistics[axis].capture(history, torque[axis])?;
            }
        }
        Ok(())
    }
    pub(super) fn leave(
        &mut self,
        channel: usize,
        position: Point3,
        rate: [f64; 3],
    ) -> Result<(), VolumeTransportError> {
        self.escape[channel] = moment_about(self.reference, position, rate)?;
        Ok(())
    }
    pub(super) fn finish<const N: usize>(
        self,
        tallies: &mut OrbitalTallies<N>,
        lane: usize,
    ) -> Result<(), VolumeTransportError> {
        for axis in 0..3 {
            tallies.torque[lane][axis].observe_signed(self.torque[axis].value())?;
            for channel in 0..3 {
                tallies.escape[channel][lane][axis].observe_signed(self.escape[channel][axis])?;
            }
        }
        Ok(())
    }
}
