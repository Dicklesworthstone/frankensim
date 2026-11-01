//! Incoherent irradiance from a bounded set of actual optical branches.
//!
//! Each member supplies a launch guess and a deterministic dielectric policy.
//! The existing spectral connector locates every surface and solves every
//! wavelength, preserving physical polarized power and the original angular
//! footprint. Compensated addition then sums the independently connected powers.
//! A failure of any member refuses the complete reading; missing paths are not
//! silently assigned zero power.
//!
//! Native ordered branch signatures must be distinct. This prevents a second
//! launch guess or an equivalent policy from counting the same branch twice.
//! Distinct roots with the SAME interaction signature are conservatively refused:
//! this owner does not certify image separation or enumerate all optical paths.
//! Coatings retain their coherent internal response, but interference BETWEEN
//! the separately connected paths is outside this incoherent source model.

use super::super::{BundleIrradianceMeter, BundleIrradianceReading, BundlePhotometryError, Sum};
use super::{
    BundleBranchPolicy, ConnectedIrradianceError, ConnectedIrradianceRequest,
    ConnectedIrradianceSpectrum, CurvedPolarizedScene, PhotometricSource, PolarizedSource,
    SPECTRAL_BINS,
};
use crate::polarization::dielectric::scene::curved::bundle::connection::BundleConnectionError;
use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_geom::Vec3;

/// Fixed admission bound on explicitly selected distinct optical branches.
pub const MAX_IRRADIANCE_BRANCHES: usize = 32;

/// A local guess for one reflection/transmission history. Geometry, receiver,
/// source, polarization axes and numerical controls are shared by the request.
#[derive(Clone, Copy, Debug)]
pub struct ConnectedIrradianceBranch<'a> {
    pub initial_direction: Vec3,
    pub policy: BundleBranchPolicy<'a>,
}

#[derive(Clone, Debug, PartialEq)]
pub enum MultipathIrradianceError {
    EmptyBranches,
    TooManyBranches,
    DuplicateBranch {
        first: usize,
        second: usize,
    },
    Branch {
        branch: usize,
        error: ConnectedIrradianceError,
    },
    Photometry(BundlePhotometryError),
    NumericalRange,
}
impl core::fmt::Display for MultipathIrradianceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "multipath spectral irradiance: {self:?}")
    }
}
impl core::error::Error for MultipathIrradianceError {}

/// One complete spectral branch's contribution and measured numerical work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IrradianceBranchSummary {
    pub signature: ContentHash,
    pub reading: BundleIrradianceReading,
    pub max_endpoint_error_m: f64,
    pub max_contact_resolution_m: f64,
    pub trace_evaluations: u64,
    pub field_evaluations: u64,
    pub query_work: u64,
}

const EMPTY_SUMMARY: IrradianceBranchSummary = IrradianceBranchSummary {
    signature: ContentHash([0; 32]),
    reading: BundleIrradianceReading {
        radiant_irradiance_w_m2: 0.0,
        illuminance_lux: 0.0,
    },
    max_endpoint_error_m: 0.0,
    max_contact_resolution_m: 0.0,
    trace_evaluations: 0,
    field_evaluations: 0,
    query_work: 0,
};

/// Complete incoherent sum over the declared branch set, using fixed storage.
/// It is not a completeness claim for all images or an unresolved caustic flux.
#[derive(Clone, Debug, PartialEq)]
pub struct MultipathIrradianceSpectrum {
    bins_w_m2_nm: [f64; SPECTRAL_BINS],
    reading: BundleIrradianceReading,
    branches: [IrradianceBranchSummary; MAX_IRRADIANCE_BRANCHES],
    branch_count: usize,
    max_endpoint_error_m: f64,
    max_contact_resolution_m: f64,
    trace_evaluations: u64,
    field_evaluations: u64,
    query_work: u64,
}

impl MultipathIrradianceSpectrum {
    pub const fn spectral_irradiance_w_m2_nm(&self) -> &[f64; SPECTRAL_BINS] {
        &self.bins_w_m2_nm
    }
    pub const fn reading(&self) -> BundleIrradianceReading {
        self.reading
    }
    pub fn branches(&self) -> &[IrradianceBranchSummary] {
        &self.branches[..self.branch_count]
    }
    pub const fn max_endpoint_error_m(&self) -> f64 {
        self.max_endpoint_error_m
    }
    /// Largest contact collar, not a bound on physical irradiance error.
    pub const fn max_contact_resolution_m(&self) -> f64 {
        self.max_contact_resolution_m
    }
    pub const fn trace_evaluations(&self) -> u64 {
        self.trace_evaluations
    }
    pub const fn field_evaluations(&self) -> u64 {
        self.field_evaluations
    }
    pub const fn query_work(&self) -> u64 {
        self.query_work
    }
}

impl PhotometricSource {
    /// Sum distinct selected optical branches from this unpolarized fixture.
    /// Interface-induced polarization is preserved independently in every path.
    pub fn connect_irradiance_branches(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedPolarizedScene<'_>,
        request: ConnectedIrradianceRequest,
        branches: &[ConnectedIrradianceBranch<'_>],
    ) -> Result<MultipathIrradianceSpectrum, MultipathIrradianceError> {
        accumulate(self, cx, request, branches, |member| {
            self.connect_irradiance_with_policy(
                cx,
                scene,
                ConnectedIrradianceRequest {
                    initial_direction: member.initial_direction,
                    ..request
                },
                member.policy,
            )
        })
    }
}

impl PolarizedSource {
    /// Sum selected branches without restoring upstream lost source power or
    /// averaging away the different reflected/transmitted polarization states.
    pub fn connect_irradiance_branches(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedPolarizedScene<'_>,
        request: ConnectedIrradianceRequest,
        branches: &[ConnectedIrradianceBranch<'_>],
    ) -> Result<MultipathIrradianceSpectrum, MultipathIrradianceError> {
        accumulate(self.intensity_source(), cx, request, branches, |member| {
            self.connect_irradiance_with_policy(
                cx,
                scene,
                ConnectedIrradianceRequest {
                    initial_direction: member.initial_direction,
                    ..request
                },
                member.policy,
            )
        })
    }
}

fn accumulate(
    source: &PhotometricSource,
    cx: &Cx<'_>,
    request: ConnectedIrradianceRequest,
    branches: &[ConnectedIrradianceBranch<'_>],
    mut connect: impl FnMut(
        ConnectedIrradianceBranch<'_>,
    ) -> Result<ConnectedIrradianceSpectrum, ConnectedIrradianceError>,
) -> Result<MultipathIrradianceSpectrum, MultipathIrradianceError> {
    if branches.is_empty() {
        return Err(MultipathIrradianceError::EmptyBranches);
    }
    if branches.len() > MAX_IRRADIANCE_BRANCHES {
        return Err(MultipathIrradianceError::TooManyBranches);
    }
    let mut result = MultipathIrradianceSpectrum {
        bins_w_m2_nm: [0.0; SPECTRAL_BINS],
        reading: EMPTY_SUMMARY.reading,
        branches: [EMPTY_SUMMARY; MAX_IRRADIANCE_BRANCHES],
        branch_count: branches.len(),
        max_endpoint_error_m: 0.0,
        max_contact_resolution_m: 0.0,
        trace_evaluations: 0,
        field_evaluations: 0,
        query_work: 0,
    };
    let mut totals: [Sum; SPECTRAL_BINS] = core::array::from_fn(|_| Sum::default());
    for (branch, member) in branches.iter().copied().enumerate() {
        cx.checkpoint()
            .map_err(|_| MultipathIrradianceError::Branch {
                branch,
                error: ConnectedIrradianceError::Connection {
                    bin: 0,
                    error: BundleConnectionError::Cancelled,
                },
            })?;
        let spectrum =
            connect(member).map_err(|error| MultipathIrradianceError::Branch { branch, error })?;
        let signature = spectrum.branch_signature();
        if let Some(first) = result.branches[..branch]
            .iter()
            .position(|b| b.signature == signature)
        {
            return Err(MultipathIrradianceError::DuplicateBranch {
                first,
                second: branch,
            });
        }
        result.branches[branch] = IrradianceBranchSummary {
            signature,
            reading: spectrum.reading(),
            max_endpoint_error_m: spectrum.max_endpoint_error_m(),
            max_contact_resolution_m: spectrum.max_contact_resolution_m(),
            trace_evaluations: spectrum.trace_evaluations(),
            field_evaluations: spectrum.field_evaluations(),
            query_work: spectrum.query_work(),
        };
        for (total, value) in totals
            .iter_mut()
            .zip(spectrum.spectral_irradiance_w_m2_nm())
        {
            total.add(*value);
        }
        result.max_endpoint_error_m = result
            .max_endpoint_error_m
            .max(spectrum.max_endpoint_error_m());
        result.max_contact_resolution_m = result
            .max_contact_resolution_m
            .max(spectrum.max_contact_resolution_m());
        result.trace_evaluations = result
            .trace_evaluations
            .checked_add(spectrum.trace_evaluations())
            .ok_or(MultipathIrradianceError::NumericalRange)?;
        result.field_evaluations = result
            .field_evaluations
            .checked_add(spectrum.field_evaluations())
            .ok_or(MultipathIrradianceError::NumericalRange)?;
        result.query_work = result
            .query_work
            .checked_add(spectrum.query_work())
            .ok_or(MultipathIrradianceError::NumericalRange)?;
    }
    result.bins_w_m2_nm = totals.map(Sum::value);
    // Reuse the original calibrated quadrature and its finite-range admission.
    let mut meter =
        BundleIrradianceMeter::new(source, request.receiver, request.connection.projection);
    for (bin, value) in result.bins_w_m2_nm.iter().copied().enumerate() {
        meter
            .record_value(bin, 1_u128 << bin, value)
            .map_err(MultipathIrradianceError::Photometry)?;
    }
    result.reading = meter
        .reading()
        .map_err(MultipathIrradianceError::Photometry)?;
    Ok(result)
}

#[cfg(test)]
mod tests;
