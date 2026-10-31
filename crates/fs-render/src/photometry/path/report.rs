//! Interpretation of accepted photometric samples and their declared settings.

use crate::tracer::Sampler;
use fs_evidence::{NumericalCertificate, NumericalKind};
use std::io::{self, Write};

use super::{DiffuseTransportRemainder, PathIlluminanceJob, PathLuminanceJob, PointFixtureSampling};

/// Sampling uncertainty of an accepted path estimate, in the units named by
/// the measurement accessor. This does not include numerical, calibration,
/// spectral, geometry/material/medium or omitted-transport error.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SamplingUncertainty {
    /// Fewer than two IID paths cannot supply a sample variance.
    InsufficientIidSamples {
        /// Number of accepted whole paths.
        samples: u32,
    },
    /// Sample-based standard error of the mean, not a confidence interval or
    /// error bound. Zero does not rule out unsampled paths or model error.
    EstimatedIidStandardError {
        /// Standard error in the measurement accessor's stated SI units.
        standard_error: f64,
        /// Number of accepted whole paths.
        samples: u32,
    },
    /// Correlated primary samples do not justify the IID variance formula.
    CorrelatedSamples {
        /// Number of accepted whole paths; increasing it does not create an
        /// IID standard-error claim for this sampler.
        samples: u32,
    },
}

impl SamplingUncertainty {
    pub(super) fn from_estimate(
        sampler: Sampler,
        samples: u32,
        standard_error: Option<f64>,
    ) -> Self {
        match (sampler, standard_error) {
            (Sampler::Iid, Some(standard_error)) => Self::EstimatedIidStandardError {
                standard_error,
                samples,
            },
            (Sampler::Iid, None) => Self::InsufficientIidSamples { samples },
            _ => Self::CorrelatedSamples { samples },
        }
    }
}

/// Declared calculation settings copied from an admitted path job.
///
/// This is neither a complete scene identity nor an accuracy certificate.
/// Replay also requires the original scene, fixtures, medium, measurement
/// geometry and remaining options. Calibration, spectral discretization,
/// geometry/material/medium errors remain unbounded. The separate optional
/// [`super::DiffuseTransportRemainder`] addresses only a restricted finite-depth
/// model remainder; these settings alone establish no bound.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PathCalculationContext {
    /// Actual primary-sample policy of this job.
    pub sampler: Sampler,
    /// Declared seed used by its counter and scrambled primary streams.
    pub seed: u64,
    /// Declared native path vertex ceiling, not a bound on the omitted tail.
    pub max_depth: u32,
    /// CIE source-data identity and spectral sampling convention.
    pub observer_id: &'static str,
    /// Version of the shared direct photometric model.
    pub photometry_semantics_version: u32,
    /// Version of fixture-path admission, transport and accumulation.
    pub path_semantics_version: u32,
}

impl PathCalculationContext {
    pub(super) const fn new(sampler: Sampler, seed: u64, max_depth: u32) -> Self {
        Self {
            sampler,
            seed,
            max_depth,
            observer_id: crate::photometry::observer::OBSERVER_ID,
            photometry_semantics_version: crate::photometry::PHOTOMETRY_SEMANTICS_VERSION,
            path_semantics_version: super::PHOTOMETRIC_PATH_SEMANTICS_VERSION,
        }
    }
}

/// The scalar quantity to which every scalar error in a report refers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhotometricQuantity {
    /// Incident luminous flux per receiving area, in lux.
    Illuminance,
    /// Photopic radiance, in candela per square metre.
    Luminance,
}

impl PhotometricQuantity {
    /// SI unit label used in the CSV output.
    #[must_use]
    pub const fn unit(self) -> &'static str {
        match self {
            Self::Illuminance => "lx",
            Self::Luminance => "cd/m2",
        }
    }
}

/// Error accounting attached to one accepted probe or pixel, before exposure.
///
/// A model-tail bound and an estimated standard error are different objects.
/// They are never added to manufacture a total accuracy bound. In particular,
/// a zero diffuse tail or zero sample variance leaves unknown calibration,
/// discretization, geometry and physical-model errors unknown. The numerical
/// certificate is copied without promotion from the original path estimate.
///
/// The context and point-fixture policy name the actual admitted calculation.
/// The optional remainder retains its source ceiling, surface/volume
/// contraction and depth. These are provenance data, not a complete scene hash
/// or authenticated source measurements. Keep the original scene for replay.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PathErrorBudget {
    /// Row-major pixel index or original probe index.
    pub measurement_index: usize,
    /// Physical interpretation and units of `value` and scalar error fields.
    pub quantity: PhotometricQuantity,
    /// Original accepted luminance or total incident illuminance estimate.
    pub value: f64,
    /// Whole-path sampling uncertainty, never a confidence certificate.
    pub sampling: SamplingUncertainty,
    /// Original numerical claim, not strengthened by an independent model tail.
    pub numerical: NumericalCertificate,
    /// Omitted orders for an eligible diffuse surface/volume transport model.
    pub diffuse_transport_remainder: Option<DiffuseTransportRemainder>,
    /// Fixture/source calibration error. Currently unknown, not zero.
    pub calibration_error_bound: Option<f64>,
    /// Spectral discretization and observer approximation error; currently unknown.
    pub spectral_error_bound: Option<f64>,
    /// Geometry, intersections and ray-offset error; currently unknown.
    pub geometry_error_bound: Option<f64>,
    /// Physical material/medium/source-model discrepancy; currently unknown.
    pub physical_model_error_bound: Option<f64>,
    /// Original total-error bound. Unknown terms are never replaced by zero.
    pub total_error_bound: Option<f64>,
    /// Actual sampler, seed, depth, observer and transport versions.
    pub calculation: PathCalculationContext,
    /// Actual policy for sampled point-fixture connections.
    pub point_sampling: PointFixtureSampling,
}

impl PathErrorBudget {
    /// The conditional transport-tail bound in this measurement's units.
    /// `None` means unsupported or unavailable, never zero omitted light.
    #[must_use]
    pub fn transport_truncation_upper_bound(self) -> Option<f64> {
        let remainder = self.diffuse_transport_remainder?;
        Some(match self.quantity {
            PhotometricQuantity::Illuminance => remainder.illuminance_bound_lux(),
            PhotometricQuantity::Luminance => remainder.luminance_bound_cd_m2(),
        })
    }

    /// Write a header for [`Self::write_csv_row`]. No file or stdout is opened.
    pub fn write_csv_header(writer: &mut impl Write) -> io::Result<()> {
        writer.write_all(
            concat!(
                "index,unit,value,samples,sampling_uncertainty,estimated_standard_error,",
                "numerical_claim,numerical_lower,numerical_upper,transport_remainder_model,",
                "transport_truncation_upper_bound,calibration_error_bound,spectral_error_bound,",
                "geometry_error_bound,physical_model_error_bound,total_error_bound,seed,max_depth,",
                "point_sampling,connections_per_vertex,photometry_version,path_version,observer_id,",
                "max_surface_reflectance,source_luminance_ceiling_cd_m2,sampler,",
                "max_scattering_gain,contraction_upper_bound\n"
            )
            .as_bytes(),
        )
    }

    /// Stream one record with round-trip binary64 values and explicit `unknown`
    /// fields. The transport bound is the named diffuse surface/volume remainder,
    /// not a symmetric interval around `value`; standard error stays separate.
    /// Observer text is CSV-quoted. An I/O failure is returned to the caller;
    /// atomic file publication, when needed, belongs to the caller's writer.
    pub fn write_csv_row(self, writer: &mut impl Write) -> io::Result<()> {
        let (samples, sampling, standard_error) = match self.sampling {
            SamplingUncertainty::InsufficientIidSamples { samples } => {
                (samples, "insufficient-iid-samples", None)
            }
            SamplingUncertainty::EstimatedIidStandardError {
                samples,
                standard_error,
            } => (samples, "estimated-iid-standard-error", Some(standard_error)),
            SamplingUncertainty::CorrelatedSamples { samples } => {
                (samples, "correlated-samples", None)
            }
        };
        let numerical = match self.numerical.kind {
            NumericalKind::Exact => "exact",
            NumericalKind::Enclosure => "enclosure",
            NumericalKind::Estimate => "estimate",
            NumericalKind::NoClaim => "no-claim",
        };
        write!(
            writer,
            "{},{},{},{samples},{sampling},",
            self.measurement_index,
            self.quantity.unit(),
            self.value
        )?;
        write_optional(writer, standard_error)?;
        write!(writer, ",{numerical}")?;
        for endpoint in [self.numerical.lo, self.numerical.hi] {
            writer.write_all(b",")?;
            write_optional(
                writer,
                (self.numerical.kind != NumericalKind::NoClaim).then_some(endpoint),
            )?;
        }
        writer.write_all(if self.diffuse_transport_remainder.is_some() {
            b",diffuse-surface-volume-geometric-tail-v1"
        } else {
            b",unknown"
        })?;
        for bound in [
            self.transport_truncation_upper_bound(),
            self.calibration_error_bound,
            self.spectral_error_bound,
            self.geometry_error_bound,
            self.physical_model_error_bound,
            self.total_error_bound,
        ] {
            writer.write_all(b",")?;
            write_optional(writer, bound)?;
        }
        let (policy, connections) = match self.point_sampling {
            PointFixtureSampling::All => ("all", None),
            PointFixtureSampling::Importance {
                connections_per_vertex,
            } => ("importance", Some(connections_per_vertex)),
        };
        write!(
            writer,
            ",{},{},{policy},",
            self.calculation.seed,
            self.calculation.max_depth
        )?;
        match connections {
            Some(count) => write!(writer, "{count}")?,
            None => writer.write_all(b"not-applicable")?,
        }
        write!(
            writer,
            ",{},{},\"",
            self.calculation.photometry_semantics_version,
            self.calculation.path_semantics_version
        )?;
        for part in self.calculation.observer_id.split_inclusive('"') {
            writer.write_all(part.as_bytes())?;
            if part.ends_with('"') {
                writer.write_all(b"\"")?;
            }
        }
        writer.write_all(b"\",")?;
        write_optional(
            writer,
            self.diffuse_transport_remainder
                .map(|r| r.max_surface_reflectance()),
        )?;
        writer.write_all(b",")?;
        write_optional(
            writer,
            self.diffuse_transport_remainder
                .map(|r| r.source_luminance_ceiling_cd_m2()),
        )?;
        let sampler = match self.calculation.sampler {
            Sampler::Iid => "iid",
            Sampler::OwenSobol => "owen-sobol",
            Sampler::OwenSobolFullPath => "owen-sobol-full-path",
        };
        write!(writer, ",{sampler},")?;
        write_optional(
            writer,
            self.diffuse_transport_remainder.map(|r| r.max_scattering_gain()),
        )?;
        writer.write_all(b",")?;
        write_optional(
            writer,
            self.diffuse_transport_remainder
                .map(|r| r.contraction_upper_bound()),
        )?;
        writer.write_all(b"\n")
    }
}

fn write_optional(writer: &mut impl Write, value: Option<f64>) -> io::Result<()> {
    match value {
        Some(value) => write!(writer, "{value}"),
        None => writer.write_all(b"unknown"),
    }
}

impl PathLuminanceJob<'_, '_> {
    /// Copy a pixel's accepted value, uncertainty slices and calculation context.
    /// Uncomputed and out-of-range pixels return `None`. This does no tracing,
    /// sampling or scene re-admission and cannot change the accepted prefix.
    #[must_use]
    pub fn error_budget_cd_m2(&self, index: usize) -> Option<PathErrorBudget> {
        let estimate = self.pixel_estimate(index)?;
        Some(PathErrorBudget {
            measurement_index: index,
            quantity: PhotometricQuantity::Luminance,
            value: estimate.luminance_cd_m2,
            sampling: self.sampling_uncertainty_cd_m2(index)?,
            numerical: estimate.numerical,
            diffuse_transport_remainder: self.diffuse_transport_remainder(),
            calibration_error_bound: None,
            spectral_error_bound: None,
            geometry_error_bound: None,
            physical_model_error_bound: None,
            total_error_bound: estimate.total_error_bound_cd_m2,
            calculation: self.calculation_context(),
            point_sampling: self.point_sampling(),
        })
    }
}

impl PathIlluminanceJob<'_, '_> {
    /// Copy a probe's total incident lux, uncertainty slices and retained context.
    /// Uncomputed and out-of-range probes return `None`. The direct point-source
    /// term remains deterministic; its zero sampling variance is not an accuracy
    /// certificate. The original native material and medium owners are unchanged.
    #[must_use]
    pub fn error_budget_lux(&self, index: usize) -> Option<PathErrorBudget> {
        let estimate = self.probe_estimate(index)?;
        Some(PathErrorBudget {
            measurement_index: index,
            quantity: PhotometricQuantity::Illuminance,
            value: estimate.illuminance_lux,
            sampling: self.sampling_uncertainty_lux(index)?,
            numerical: estimate.numerical,
            diffuse_transport_remainder: self.diffuse_transport_remainder(),
            calibration_error_bound: None,
            spectral_error_bound: None,
            geometry_error_bound: None,
            physical_model_error_bound: None,
            total_error_bound: estimate.total_error_bound_lux,
            calculation: self.calculation_context(),
            point_sampling: self.point_sampling(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn report() -> PathErrorBudget {
        PathErrorBudget {
            measurement_index: 3,
            quantity: PhotometricQuantity::Illuminance,
            value: 12.5,
            sampling: SamplingUncertainty::EstimatedIidStandardError {
                standard_error: 0.0,
                samples: 2,
            },
            numerical: NumericalCertificate::no_claim(),
            diffuse_transport_remainder: None,
            calibration_error_bound: None,
            spectral_error_bound: None,
            geometry_error_bound: None,
            physical_model_error_bound: None,
            total_error_bound: None,
            calculation: PathCalculationContext {
                observer_id: "declared, \"observer\"",
                photometry_semantics_version: 6,
                path_semantics_version: 12,
                ..PathCalculationContext::new(Sampler::Iid, 173, 4)
            },
            point_sampling: PointFixtureSampling::All,
        }
    }

    #[test]
    fn csv_zero_sampling_error_does_not_become_a_known_physical_bound() {
        let mut bytes = Vec::new();
        report().write_csv_row(&mut bytes).unwrap();
        assert_eq!(
            String::from_utf8(bytes).unwrap(),
            concat!(
                "3,lx,12.5,2,estimated-iid-standard-error,0,no-claim,unknown,unknown,",
                "unknown,unknown,unknown,unknown,unknown,unknown,unknown,173,4,all,",
                "not-applicable,6,12,\"declared, \"\"observer\"\"\",unknown,unknown,iid,unknown,unknown\n"
            )
        );
    }

    #[test]
    fn csv_retains_zero_model_tail_but_not_a_zero_total_bound() {
        let mut value = report();
        value.quantity = PhotometricQuantity::Luminance;
        value.diffuse_transport_remainder =
            DiffuseTransportRemainder::from_envelope(4, 0.0, 0.0, fs_ivl::Interval::point(8.0));
        value.sampling = SamplingUncertainty::CorrelatedSamples { samples: 8 };
        value.calculation.sampler = Sampler::OwenSobol;
        value.point_sampling = PointFixtureSampling::Importance {
            connections_per_vertex: 2,
        };
        let mut bytes = Vec::new();
        value.write_csv_row(&mut bytes).unwrap();
        let csv = String::from_utf8(bytes).unwrap();
        assert!(csv.starts_with("3,cd/m2,12.5,8,correlated-samples,unknown,no-claim,"));
        assert!(csv.contains(
            ",diffuse-surface-volume-geometric-tail-v1,0,unknown,unknown,unknown,unknown,unknown,"
        ));
        assert!(csv.contains(",importance,2,6,12,"));
        assert!(csv.ends_with(",0,8,owen-sobol,0,0\n"));
        assert_eq!(value.transport_truncation_upper_bound(), Some(0.0));
        assert!(value.total_error_bound.is_none());
    }

    #[test]
    fn csv_keeps_volume_contraction_and_measurement_units_explicit() {
        let mut value = report();
        let remainder = DiffuseTransportRemainder::from_envelope(
            4,
            0.25,
            0.5,
            fs_ivl::Interval::point(8.0),
        )
        .unwrap();
        value.diffuse_transport_remainder = Some(remainder);
        assert_eq!(
            value.transport_truncation_upper_bound(),
            Some(remainder.illuminance_bound_lux())
        );
        value.quantity = PhotometricQuantity::Luminance;
        assert_eq!(
            value.transport_truncation_upper_bound(),
            Some(remainder.luminance_bound_cd_m2())
        );
        let mut bytes = Vec::new();
        value.write_csv_row(&mut bytes).unwrap();
        assert!(String::from_utf8(bytes).unwrap().ends_with(",0.25,8,iid,0.5,0.5\n"));
        assert!(value.total_error_bound.is_none());
    }

    #[test]
    fn csv_propagates_writer_failure() {
        struct Refuse;
        impl Write for Refuse {
            fn write(&mut self, _: &[u8]) -> io::Result<usize> {
                Err(io::Error::other("test sink refusal"))
            }
            fn flush(&mut self) -> io::Result<()> {
                Ok(())
            }
        }
        assert!(PathErrorBudget::write_csv_header(&mut Refuse).is_err());
        assert!(report().write_csv_row(&mut Refuse).is_err());
    }
}
