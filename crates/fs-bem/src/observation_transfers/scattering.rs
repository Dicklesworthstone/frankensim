//! Nearby prescribed sound sources observed around closed rigid obstacles.
//!
//! The original Helmholtz owner resolves the scattered wave, including
//! diffraction, and its adaptive receiver integrates that same boundary trace.
//! The exact spherical incident pressure AND gradient are then added. In
//! particular, an obstructed direct path is not discarded: its cancellation
//! by the scattered field is part of the wave solution, not a visibility test.
//!
//! This is a stationary, homogeneous, linear, one-way observation model. A
//! point source is explicit input, not an automatic replacement for a finite
//! aperture, vibrating body, or a source with acoustic loading. The common FIR
//! owner retains its frequency/tail ladder and declared reconstruction delay.

use super::*;
use crate::helmholtz::bc::{
    Admittance,
    incident::{self, IncidentField, Monopole},
};
use crate::near_field::{MAX_RECEIVERS, Options};
use std::cell::RefCell;
use std::collections::BTreeSet;

/// An independently driven point source of outward volume flow [m³/s].
#[derive(Clone, Debug, PartialEq)]
pub struct PointSourceSpec {
    /// Caller-owned source identity, distinct within the input basis.
    pub id: String,
    /// Physical point in the same metre-valued frame as obstacles and receivers.
    pub position_m: [f64; 3],
}

/// A complete source-to-receiver pressure/gradient observation of diffraction
/// around fixed, closed rigid bodies. No in-loop impedance is inferred.
pub struct ScatteringObservation<'a> {
    pub(super) surface: &'a SpherePanels,
    pub(super) receiver_geometry: Geometry<'a>,
    pub(super) medium: Medium,
    pub(super) point_sources: Vec<PointSourceSpec>,
    sources: Vec<SourceId>,
    pub(super) receivers: Vec<ReceiverSpec>,
    pub(super) highest_omega: f64,
    pub(super) minimum_clearance_m: f64,
    pub(super) quadrature: Options,
    delays_s: Vec<f64>,
}

impl<'a> ScatteringObservation<'a> {
    /// Bind the immutable physical source/obstacle/receiver geometry and a
    /// finite band. Sources and receivers must lie outside every closed body
    /// and obey actual triangle clearance. A receiver cannot coincide with a
    /// point source. Existing dense-panel and receiver-quadrature caps apply.
    ///
    /// All frequencies use the same Burton–Miller discretization; no operator
    /// switch can introduce a spurious spectral discontinuity into the FIR.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        surface: &'a SpherePanels,
        medium: Medium,
        sources: &[PointSourceSpec],
        receivers: &[ReceiverSpec],
        highest_omega: f64,
        minimum_clearance_m: f64,
        quadrature: Options,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, TransferError> {
        Self::new_admitted(
            surface,
            medium,
            sources,
            receivers,
            highest_omega,
            minimum_clearance_m,
            quadrature,
            None,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn new_admitted(
        surface: &'a SpherePanels,
        medium: Medium,
        sources: &[PointSourceSpec],
        receivers: &[ReceiverSpec],
        highest_omega: f64,
        minimum_clearance_m: f64,
        quadrature: Options,
        iterative: Option<&fs_exec::CancelGate>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, TransferError> {
        if cancelled() {
            return Err(TransferError::Cancelled);
        }
        quadrature.validate()?;
        if sources.is_empty()
            || sources.len() > MAX_RECEIVERS
            || receivers.is_empty()
            || receivers.len() > MAX_RECEIVERS
            || !highest_omega.is_finite()
            || highest_omega <= 0.0
            || !medium.density.is_finite()
            || medium.density <= 0.0
            || !medium.sound_speed.is_finite()
            || medium.sound_speed <= 0.0
            || !(medium.density * medium.sound_speed).is_finite()
            || !(highest_omega * medium.density).is_finite()
            || !(highest_omega / medium.sound_speed).is_finite()
            || !minimum_clearance_m.is_finite()
            || minimum_clearance_m <= 0.0
            || sources.iter().any(|s| s.id.is_empty() || s.id.len() > 1024)
            || receivers
                .iter()
                .any(|r| r.id.is_empty() || r.id.len() > 1024)
            || sources.iter().map(|s| &s.id).collect::<BTreeSet<_>>().len() != sources.len()
            || receivers
                .iter()
                .map(|r| &r.id)
                .collect::<BTreeSet<_>>()
                .len()
                != receivers.len()
        {
            return Err(TransferError::Invalid(
                "invalid point-source diffraction declaration",
            ));
        }
        // Refuse an unresolved upper band before geometry preparation or
        // lower-frequency solves; the same monotone gate owns every sample.
        if iterative.is_some() {
            helmholtz::admit_iterative_resolution(surface, highest_omega / medium.sound_speed)?;
        } else {
            helmholtz::admit_resolution(surface, highest_omega / medium.sound_speed)?;
        }
        let source_points: Vec<_> = sources.iter().map(|s| s.position_m).collect();
        let admit_geometry = |points: &[[f64; 3]]| match iterative {
            Some(gate) => Geometry::new_iterative(surface, points, minimum_clearance_m, gate),
            None => Geometry::new(surface, points, minimum_clearance_m),
        };
        admit_geometry(&source_points)?;
        if cancelled() {
            return Err(TransferError::Cancelled);
        }
        let points: Vec<_> = receivers.iter().map(|r| r.position_m).collect();
        let receiver_geometry = admit_geometry(&points)?;
        let mut delays_s = Vec::with_capacity(sources.len() * receivers.len());
        for source in sources {
            for receiver in receivers {
                let distance = det::sqrt(
                    source
                        .position_m
                        .iter()
                        .zip(receiver.position_m)
                        .map(|(s, r)| (s - r).powi(2))
                        .sum::<f64>(),
                );
                let delay = distance / medium.sound_speed;
                if !distance.is_finite() || distance < minimum_clearance_m || !delay.is_finite() {
                    return Err(TransferError::Invalid(
                        "point source and receiver separation is unresolved",
                    ));
                }
                // The straight-line travel time is a conservative earliest
                // arrival even when the segment intersects a rigid obstacle.
                delays_s.push(delay);
            }
        }
        if cancelled() {
            return Err(TransferError::Cancelled);
        }
        Ok(Self {
            surface,
            receiver_geometry,
            medium,
            point_sources: sources.to_vec(),
            sources: sources
                .iter()
                .map(|s| SourceId::PointVolumeFlow(s.id.clone()))
                .collect(),
            receivers: receivers.to_vec(),
            highest_omega,
            minimum_clearance_m,
            quadrature,
            delays_s,
        })
    }

    /// Complete input basis; each signal supplied to the runtime is in m³/s.
    pub fn sources(&self) -> &[SourceId] {
        &self.sources
    }

    /// Source-major, receiver-minor conservative earliest arrivals [s].
    pub fn geometric_delays_s(&self) -> Vec<f64> {
        self.delays_s.clone()
    }

    /// Bind the physical geometry, numerical band and receiver quadrature.
    /// The common transfer declaration additionally binds its FIR resolution,
    /// uncertainty, source/receiver identities and reconstruction latency.
    pub fn operator_hash(&self) -> ContentHash {
        let mut h = DomainHasher::new("fs-bem:point-source-rigid-diffraction:v1");
        for count in [
            self.surface.areas().len(),
            self.point_sources.len(),
            self.receivers.len(),
            self.quadrature.maximum_depth,
            self.quadrature.maximum_kernel_evaluations,
        ] {
            h.update(&(count as u64).to_le_bytes());
        }
        for value in [
            self.medium.density,
            self.medium.sound_speed,
            self.highest_omega,
            self.minimum_clearance_m,
            self.quadrature.relative_tolerance,
        ]
        .into_iter()
        .chain(self.surface.centroids().iter().flatten().copied())
        .chain(self.surface.normals().iter().flatten().copied())
        .chain(self.surface.areas().iter().copied())
        .chain(self.point_sources.iter().flat_map(|s| s.position_m))
        .chain(self.receivers.iter().flat_map(|r| r.position_m))
        {
            h.update(&value.to_bits().to_le_bytes());
        }
        for triangle in self.surface.triangles().into_iter().flatten() {
            for value in triangle.iter().flatten() {
                h.update(&value.to_bits().to_le_bytes());
            }
        }
        h.finalize()
    }

    /// Source-major, receiver-major `(p, dp/dx, dp/dy, dp/dz)` in the original
    /// exp(-iωt) convention, per unit prescribed source volume flow.
    ///
    /// One boundary assembly/LU serves the entire source basis, and one set of
    /// adaptive Green rows serves every receiver/source pair. Cancellation is
    /// polled between bounded stages; the inherited dense LU is not interruptible.
    pub fn frequency_sample(
        &self,
        omega: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<C64>, HelmholtzError> {
        self.frequency_sample_with_admittance(omega, Admittance::Rigid, cancelled)
    }

    pub(super) fn frequency_sample_with_admittance(
        &self,
        omega: f64,
        admittance: Admittance<'_>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<C64>, HelmholtzError> {
        if cancelled() {
            return Err(HelmholtzError::Cancelled);
        }
        if !omega.is_finite() || omega < 0.0 || omega > self.highest_omega {
            return Err(HelmholtzError::BadParameter {
                what: "diffraction observation frequency outside its band",
            });
        }
        let channels = 4 * self.point_sources.len() * self.receivers.len();
        if omega == 0.0 {
            return Ok(vec![C64::ZERO; channels]);
        }
        let k = omega / self.medium.sound_speed;
        let monopoles: Vec<_> = self
            .point_sources
            .iter()
            .map(|s| {
                [Monopole {
                    position_m: s.position_m,
                    volume_velocity_m3_s: C64::ONE,
                }]
            })
            .collect();
        let fields: Vec<_> = monopoles
            .iter()
            .map(|s| IncidentField::Monopoles(s))
            .collect();
        let solutions = incident::solve_scattering_batch(
            self.surface,
            k,
            self.medium,
            admittance,
            &fields,
            self.minimum_clearance_m,
            helmholtz::Formulation::BurtonMiller,
        )?;
        if cancelled() {
            return Err(HelmholtzError::Cancelled);
        }
        let rows = self
            .receiver_geometry
            .prepare_velocity(k, self.medium, self.quadrature)?;
        let mut output = Vec::with_capacity(channels);
        let momentum = C64::new(0.0, omega * self.medium.density);
        for (solution, field) in solutions.iter().zip(fields) {
            if cancelled() {
                return Err(HelmholtzError::Cancelled);
            }
            let scattered = rows.evaluate_velocity(&solution.solution().scattered)?;
            for ((&pressure, velocity), receiver) in scattered
                .scalar
                .pressure
                .iter()
                .zip(&scattered.particle_velocity_m_s)
                .zip(&self.receivers)
            {
                let direct = field.pressure_velocity_at(k, self.medium, receiver.position_m)?;
                output.push(pressure + direct.pressure_pa);
                output.extend(
                    (0..3).map(|axis| momentum * (velocity[axis] + direct.velocity_m_s[axis])),
                );
            }
        }
        if output.len() != channels || output.iter().any(|&v| !finite(v)) {
            return Err(HelmholtzError::BadParameter {
                what: "nonfinite or incomplete total diffraction field",
            });
        }
        if cancelled() {
            return Err(HelmholtzError::Cancelled);
        }
        Ok(output)
    }

    /// Realize actual source-flow histories through the existing converged
    /// observation FIR owner. Producer/BEM uncertainty is retained separately
    /// from temporal convergence; an exact-analytic producer claim is refused.
    pub fn bake(
        &self,
        declaration: &TransferDeclaration,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ObservationTransferBake, TransferError> {
        if declaration.sources != self.sources
            || declaration.receivers != self.receivers
            || declaration.geometric_delays_s != self.delays_s
            || declaration.producer_hash != self.operator_hash()
            || core::f64::consts::TAU * declaration.taper_stop_hz > self.highest_omega
            || matches!(
                declaration.producer_uncertainty,
                ProducerUncertainty::ExactAnalytic
            )
        {
            return Err(TransferError::Invalid(
                "cross-wired point-source diffraction declaration",
            ));
        }
        let cancel = RefCell::new(cancelled);
        bake_observation_transfers(
            |omega| self.frequency_sample(omega, || (*cancel.borrow_mut())()),
            declaration,
            || (*cancel.borrow_mut())(),
        )
        .map_err(|error| match error {
            TransferError::Helmholtz(HelmholtzError::Cancelled) => TransferError::Cancelled,
            other => other,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::TAU;

    fn medium() -> Medium {
        Medium {
            density: 1.2,
            sound_speed: 343.0,
        }
    }
    fn source() -> PointSourceSpec {
        PointSourceSpec {
            id: "prescribed-flow".into(),
            position_m: [0.02, 0.01, 0.09],
        }
    }
    fn receiver() -> ReceiverSpec {
        ReceiverSpec {
            id: "shadow-side".into(),
            position_m: [-0.03, 0.02, -0.08],
        }
    }
    fn observation<'a>(surface: &'a SpherePanels) -> ScatteringObservation<'a> {
        ScatteringObservation::new(
            surface,
            medium(),
            &[source()],
            &[receiver()],
            TAU * 2000.0,
            1e-3,
            Options::default(),
            || false,
        )
        .unwrap()
    }
    fn declaration(owner: &ScatteringObservation<'_>) -> TransferDeclaration {
        TransferDeclaration {
            producer_hash: owner.operator_hash(),
            sources: owner.sources().to_vec(),
            receivers: vec![receiver()],
            fs_hz: 8000.0,
            geometric_delays_s: owner.geometric_delays_s(),
            pass_band_upper_hz: 800.0,
            taper_stop_hz: 2000.0,
            filter_latency_samples: 128,
            fft_lengths: vec![512, 1024, 2048, 4096],
            tail_energy_tolerance: 1e-10,
            dtft_relative_tolerance: 1e-4,
            producer_uncertainty: ProducerUncertainty::Unresolved {
                source_hash: owner.operator_hash(),
            },
            maximum_frequency_entries: 4 * MAX_OBSERVATION_FFT,
        }
    }

    /// G1: pressure is the actual total diffraction field and its gradient is
    /// independently checked by physical receiver displacement, not by copying
    /// the analytic derivative formula from the implementation.
    #[test]
    fn total_diffraction_gradient_matches_displaced_pressure() {
        let surface = SpherePanels::icosphere(0.025, 0).unwrap();
        let owner = observation(&surface);
        let omega = TAU * 700.0;
        let k = omega / medium().sound_speed;
        let row = owner.frequency_sample(omega, || false).unwrap();
        let monopole = [Monopole {
            position_m: source().position_m,
            volume_velocity_m3_s: C64::ONE,
        }];
        let field = IncidentField::Monopoles(&monopole);
        let solution = incident::solve_scattering(
            &surface,
            k,
            medium(),
            Admittance::Rigid,
            field,
            1e-3,
            helmholtz::Formulation::BurtonMiller,
        )
        .unwrap();
        let point = receiver().position_m;
        let direct = field.pressure_velocity_at(k, medium(), point).unwrap();
        assert!(
            (row[0] - direct.pressure_pa).abs() > 1e-3 * direct.pressure_pa.abs(),
            "the finite obstacle must alter the physical field"
        );
        // Central parity shares vector quadrature; scalar displacement remains
        // the independent gradient check with its own admitted refinement.
        let central = solution
            .total_pressure_velocity_at(&[point], 1e-3, Options::default())
            .unwrap();
        assert!((row[0] - central.scalar.pressure[0]).abs() < 1e-12 * row[0].abs());
        let h = 2e-5;
        let mut points = Vec::with_capacity(6);
        for axis in 0..3 {
            let mut plus = point;
            let mut minus = point;
            plus[axis] += h;
            minus[axis] -= h;
            points.extend([plus, minus]);
        }
        let pressure = solution
            .total_pressure_at(&points, 1e-3, Options::default())
            .unwrap();
        for axis in 0..3 {
            let gradient =
                (pressure.pressure[2 * axis] - pressure.pressure[1 + 2 * axis]).scale(0.5 / h);
            assert!(
                (row[axis + 1] - gradient).abs() < 2e-5 * gradient.abs(),
                "axis {axis}: {:?} vs {gradient:?}",
                row[axis + 1]
            );
        }
    }

    /// G0/G4: invalid source semantics and cross-wired bakes refuse before
    /// solving; DC is zero, while a cancelled DC request is still cancelled.
    #[test]
    fn declarations_refuse_unphysical_or_cross_wired_sources() {
        let surface = SpherePanels::icosphere(0.025, 0).unwrap();
        let owner = observation(&surface);
        assert_eq!(
            owner.frequency_sample(0.0, || false).unwrap(),
            vec![C64::ZERO; 4]
        );
        assert!(matches!(
            owner.frequency_sample(0.0, || true),
            Err(HelmholtzError::Cancelled)
        ));
        assert!(owner.frequency_sample(TAU * 2001.0, || false).is_err());
        assert!(matches!(
            ScatteringObservation::new(
                &surface,
                medium(),
                &[source()],
                &[receiver()],
                TAU * 1e8,
                1e-3,
                Options::default(),
                || false,
            ),
            Err(TransferError::Helmholtz(HelmholtzError::TooCoarse { .. }))
        ));
        let mut interior = source();
        interior.position_m = [0.0; 3];
        assert!(
            ScatteringObservation::new(
                &surface,
                medium(),
                &[interior],
                &[receiver()],
                TAU * 2000.0,
                1e-3,
                Options::default(),
                || false
            )
            .is_err()
        );
        let mut singular = receiver();
        singular.position_m = source().position_m;
        assert!(
            ScatteringObservation::new(
                &surface,
                medium(),
                &[source()],
                &[singular],
                TAU * 2000.0,
                1e-3,
                Options::default(),
                || false
            )
            .is_err()
        );
        let traces = vec![vec![1.0; surface.areas().len()]];
        assert!(
            BemObservation::new(
                &surface,
                medium(),
                owner.sources(),
                &traces,
                &[receiver()],
                1e-3
            )
            .is_err(),
            "point sources are not surface traces"
        );
        let mut d = declaration(&owner);
        let original_hash = d.declaration_hash();
        d.sources[0] = SourceId::ApertureVolumeFlow(source().id);
        assert_ne!(original_hash, d.declaration_hash());
        assert!(matches!(
            owner.bake(&d, || false),
            Err(TransferError::Invalid(_))
        ));
        d = declaration(&owner);
        d.producer_uncertainty = ProducerUncertainty::ExactAnalytic;
        assert!(matches!(
            owner.bake(&d, || false),
            Err(TransferError::Invalid(_))
        ));
        assert!(matches!(
            owner.bake(&declaration(&owner), || true),
            Err(TransferError::Cancelled)
        ));
        let mut polls = 0;
        assert!(matches!(
            owner.bake(&declaration(&owner), || {
                polls += 1;
                polls == 4 // Inside the first non-DC frequency sample.
            }),
            Err(TransferError::Cancelled)
        ));
    }

    /// G1/G3: the actual BEM diffraction operator reaches a causal pressure
    /// runtime. The four-rung temporal ladder does not certify the coarse
    /// spatial mesh; that uncertainty remains explicit in the artifact.
    #[test]
    fn diffraction_bake_replays_source_flow_at_physical_receiver() {
        let surface = SpherePanels::icosphere(0.025, 0).unwrap();
        let owner = observation(&surface);
        let d = declaration(&owner);
        let bake = owner.bake(&d, || false).unwrap();
        assert_eq!(bake.ladder.len(), 4);
        assert_eq!(bake.producer_uncertainty, d.producer_uncertainty);
        let omega = TAU * 500.0;
        let direct = owner.frequency_sample(omega, || false).unwrap();
        let latency = phase(-omega * d.filter_latency_samples as f64 / d.fs_hz);
        for (kernel, sample) in bake.kernels.iter().zip(direct) {
            let reference = sample.conj() * latency;
            let measured = kernel.transfer_at(omega / d.fs_hz).unwrap();
            assert!(
                (measured - reference).abs() < 2e-4 * reference.abs().max(1.0),
                "FIR must preserve physical phase and level: {measured:?} vs {reference:?}"
            );
        }
        let prefix = (bake.geometric_delays_s[0] * bake.fs_hz).floor() as usize;
        let mut observer = bake.observer().unwrap();
        let mut peak = 0.0_f64;
        for time in 0..512 {
            let output = observer
                .step(&[if time == 0 { 1e-6 } else { 0.0 }])
                .unwrap();
            assert_eq!(output.len(), 4);
            assert!(output.iter().all(|value| value.is_finite()));
            if time < prefix {
                assert_eq!(output, vec![0.0; 4]);
            }
            peak = peak.max(output[0].abs());
        }
        assert!(
            peak > 1e-4,
            "source flow must reach receiver pressure in pascals"
        );
    }
}
