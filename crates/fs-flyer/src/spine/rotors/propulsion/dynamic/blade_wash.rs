//! Fixed-geometry samples of the other rotor's finite axial wake.
//!
//! Coefficients subtract the same receiving-disk mean used by transport.
//! Adding them to the mean incident velocity recovers each actual section's
//! external wash without counting the mean or self-induction twice.

use super::{DiskAxes, DiskFrame, DiskWash, Refusal, air_error, reject};
use crate::propcoupling::PropDisk;
use fs_airscrew::Rotor;
use fs_airscrew::dynamic_inflow::first_harmonic::HarmonicInduction;
use fs_airscrew::oblique::AzimuthQuadrature;

pub(super) struct BladeWash {
    /// Per-section mean/cosine/sine coefficients in the SOURCE disk frame.
    coefficients: [Vec<[f64; 3]>; 2],
    nonzero: [[bool; 3]; 2],
    axes: [[f64; 3]; 2],
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (a[0] * b[0] + a[1] * b[1]) + a[2] * b[2]
}

impl BladeWash {
    pub(super) fn new(
        axes: DiskAxes,
        disks: &[PropDisk; 2],
        rotor: &Rotor,
        rule: &AzimuthQuadrature,
        frames: [DiskFrame; 2],
        mean: &[[DiskWash; 2]; 2],
    ) -> Result<Self, Refusal> {
        let centers = disks.map(|disk| disk.center_m);
        let mut coefficients = [Vec::new(), Vec::new()];
        let mut nonzero = [[false; 3]; 2];
        for k in 0..2 {
            let other = 1 - k;
            let average = mean[k][other];
            let mean_modes = [
                average.mean,
                dot(average.transverse, frames[other].radial_reference()),
                dot(average.transverse, frames[other].azimuth_reference()),
            ];
            let samples = &mut coefficients[k];
            *samples = rule.section_offsets(rotor, frames[k]).map_err(air_error)?;
            // Replace offsets in place: one bounded allocation per receiver,
            // all spatial gates/projections outside the repeated force path.
            for sample in samples {
                let point = core::array::from_fn(|j| centers[k][j] + sample[j]);
                let factor = axes.wash_factors(point, centers, rotor.radius_m)?[other];
                let local = if factor == 0.0 {
                    [0.0; 3]
                } else {
                    let offset =
                        core::array::from_fn(|j| (point[j] - centers[other][j]) / rotor.radius_m);
                    [
                        factor,
                        factor * dot(offset, frames[other].radial_reference()),
                        factor * dot(offset, frames[other].azimuth_reference()),
                    ]
                };
                *sample = core::array::from_fn(|j| local[j] - mean_modes[j]);
                if !sample.iter().all(|x| x.is_finite()) {
                    return Err(reject(
                        "dynamic-propulsion-invalid",
                        "nonfinite external blade-wash coefficient",
                    ));
                }
                for j in 0..3 {
                    nonzero[k][j] |= sample[j] != 0.0;
                }
            }
        }
        Ok(Self {
            coefficients,
            nonzero,
            axes: axes.directions(),
        })
    }

    pub(super) fn active(&self, receiver: usize, source: HarmonicInduction) -> bool {
        self.nonzero[receiver]
            .into_iter()
            .zip(source.modes())
            .any(|(coefficient, value)| coefficient && value != 0.0)
    }

    pub(super) fn velocity(
        &self,
        receiver: usize,
        sample: usize,
        source: HarmonicInduction,
    ) -> [f64; 3] {
        let coefficient = self.coefficients[receiver][sample];
        let mut increment = coefficient[0] * source.mean_mps;
        if source.first_mps != [0.0; 2] {
            let first = coefficient[1] * source.first_mps[0];
            let second = coefficient[2] * source.first_mps[1];
            increment += first + second;
        }
        self.axes[1 - receiver].map(|x| increment * x)
    }
}
