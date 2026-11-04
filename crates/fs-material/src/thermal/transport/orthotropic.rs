//! Orthotropic heat transport from explicit principal-direction source paths.
//!
//! The caller declares that these independently selected conductivity laws are
//! principal values in the supplied fixed material basis. Scalar measurements
//! alone do not establish that symmetry or basis. Rotation preserves their
//! original units, temperature support, material identity and source receipts.
//! This is a constitutive response, not an enthalpy or mechanical evolution.

use super::{PreparedThermalTransport, ThermalTransportError};
use fs_blake3::{ContentHash, DomainHasher};

/// Source-bound orthotropic conductivity in two or three Cartesian dimensions.
///
/// Columns of `principal_to_target` are principal material directions expressed
/// in the consumer's frame. D=2 describes a declared in-plane constitutive law;
/// it does not eliminate an unknown out-of-plane flux from a 3D material.
#[derive(Clone, Debug, PartialEq)]
pub struct PreparedOrthotropicThermalTransport<const D: usize> {
    paths: [PreparedThermalTransport; D],
    principal_to_target: [[f64; D]; D],
    constant_storage: Option<[f64; 2]>,
    identity: ContentHash,
}

/// One evaluated source model, retaining its exact prepared identity.
/// The complete original receipts remain on the prepared owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OrthotropicThermalTransportSample<const D: usize> {
    identity: ContentHash,
    temperature_k: f64,
    density: f64,
    specific_heat: f64,
    conductivity: [[f64; D]; D],
    diffusivity: [[f64; D]; D],
}

impl<const D: usize> PreparedOrthotropicThermalTransport<D> {
    /// Compose already admitted principal-direction laws; no source is queried
    /// again and no directional average or isotropy assumption is introduced.
    ///
    /// Every path must have the SAME original material card, query endpoints,
    /// density source line and heat-capacity source line, including receipts.
    /// Conductivity selection may differ by principal direction. The caller
    /// owns the physical assignment of each path to the corresponding column.
    /// A finite, proper orthonormal rotation is required; it is never repaired.
    pub fn try_new(
        paths: [PreparedThermalTransport; D],
        principal_to_target: [[f64; D]; D],
    ) -> Result<Self, ThermalTransportError> {
        if D != 2 && D != 3 {
            return Err(invalid(
                "orthotropic transport requires two or three dimensions",
            ));
        }
        admit_rotation(principal_to_target)?;
        let first = &paths[0];
        for path in &paths[1..] {
            if path.material_card_identity() != first.material_card_identity()
                || path.query_endpoints() != first.query_endpoints()
                || path.source_lines()[..2] != first.source_lines()[..2]
            {
                return Err(invalid(
                    "principal paths must share the original material, context and storage sources",
                ));
            }
        }
        // A canonical line is piecewise affine. Equal original values at ALL
        // its knots prove constant storage throughout its admitted path, even
        // when conductivity has additional breakpoints. Endpoint equality alone
        // would miss an interior density/Cp excursion.
        let storage =
            core::array::from_fn(|i| first.source_lines()[i].knots[0].1.evidence.value.value);
        let constant_storage = (0..2)
            .all(|i| {
                first.source_lines()[i]
                    .knots
                    .iter()
                    .all(|(_, value)| value.evidence.value.value.to_bits() == storage[i].to_bits())
            })
            .then_some(storage);
        let mut hash = DomainHasher::new(
            "org.frankensim.fs-material.prepared-orthotropic-thermal-transport.v1",
        );
        hash.update(&(D as u64).to_le_bytes());
        for path in &paths {
            hash.update(path.identity().as_bytes());
        }
        for entry in principal_to_target.as_flattened() {
            hash.update(&entry.to_le_bytes());
        }
        Ok(Self {
            paths,
            principal_to_target,
            constant_storage,
            identity: hash.finalize(),
        })
    }

    /// Evaluate K = Q diag(k_i) Q^T and alpha = K/(rho*c), in physical SI.
    /// Reuses the original prepared scalar evaluations without allocation or
    /// database queries. Inputs are interpolated before rotation and division.
    /// No extrapolation or coefficient clipping occurs. Consumers still own
    /// admission of their represented discrete operator and condition number.
    pub fn at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<OrthotropicThermalTransportSample<D>, ThermalTransportError> {
        let first = self.paths[0].at_temperature(temperature_k)?;
        let mut principal = [0.0; D];
        principal[0] = first.conductivity_w_m_k();
        for (axis, path) in self.paths.iter().enumerate().skip(1) {
            principal[axis] = path.at_temperature(temperature_k)?.conductivity_w_m_k();
        }
        let density = first.density_kg_m3();
        let specific_heat = first.specific_heat_j_kg_k();
        let mut conductivity = [[0.0; D]; D];
        let mut diffusivity = [[0.0; D]; D];
        for row in 0..D {
            for column in 0..=row {
                let k: f64 = (0..D)
                    .map(|axis| {
                        (self.principal_to_target[row][axis]
                            * self.principal_to_target[column][axis])
                            * principal[axis]
                    })
                    .sum();
                let alpha = (k / density) / specific_heat;
                if !k.is_finite()
                    || !alpha.is_finite()
                    || (k != 0.0 && alpha == 0.0)
                    || (row == column && (k <= 0.0 || alpha <= 0.0))
                {
                    return Err(invalid("oriented thermal transport is unrepresentable"));
                }
                // Compute each symmetric entry once, without averaging unequal
                // entries or discarding the signed off-diagonal heat flux.
                conductivity[row][column] = k;
                conductivity[column][row] = k;
                diffusivity[row][column] = alpha;
                diffusivity[column][row] = alpha;
            }
        }
        Ok(OrthotropicThermalTransportSample {
            identity: self.identity,
            temperature_k,
            density,
            specific_heat,
            conductivity,
            diffusivity,
        })
    }

    /// All original directional paths, including their unmodified receipts.
    #[must_use]
    pub const fn principal_paths(&self) -> &[PreparedThermalTransport; D] {
        &self.paths
    }
    /// Original principal-to-consumer basis.
    #[must_use]
    pub const fn principal_to_target(&self) -> [[f64; D]; D] {
        self.principal_to_target
    }
    /// Identity includes ordered source paths and the physical basis.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    /// Common inclusive source-supported interval [K].
    #[must_use]
    pub fn temperature_span_k(&self) -> [f64; 2] {
        self.paths[0].temperature_span_k()
    }
    /// (rho [kg/m3], c [J/(kg K)]) only when BOTH original source models
    /// are constant over the entire prepared temperature interval.
    #[must_use]
    pub const fn constant_storage(&self) -> Option<[f64; 2]> {
        self.constant_storage
    }
}

impl<const D: usize> OrthotropicThermalTransportSample<D> {
    /// Prepared owner identity, including every original source and basis.
    #[must_use]
    pub const fn path_identity(&self) -> ContentHash {
        self.identity
    }
    /// Actual evaluation temperature [K].
    #[must_use]
    pub const fn temperature_k(&self) -> f64 {
        self.temperature_k
    }
    /// Original common density [kg/m3].
    #[must_use]
    pub const fn density_kg_m3(&self) -> f64 {
        self.density
    }
    /// Original selected specific heat [J/(kg K)].
    #[must_use]
    pub const fn specific_heat_j_kg_k(&self) -> f64 {
        self.specific_heat
    }
    /// Full conductivity tensor in the declared consumer frame [W/(m K)].
    #[must_use]
    pub const fn conductivity_w_m_k(&self) -> [[f64; D]; D] {
        self.conductivity
    }
    /// Full thermal diffusivity in that same frame [m2/s].
    #[must_use]
    pub const fn diffusivity_m2_s(&self) -> [[f64; D]; D] {
        self.diffusivity
    }
}

fn admit_rotation<const D: usize>(q: [[f64; D]; D]) -> Result<(), ThermalTransportError> {
    if q.as_flattened().iter().any(|x| !x.is_finite()) {
        return Err(invalid("finite thermal material axes required"));
    }
    let tolerance = 256.0 * f64::EPSILON;
    for row in 0..D {
        for column in 0..D {
            let dot: f64 = (0..D).map(|k| q[k][row] * q[k][column]).sum();
            if (dot - if row == column { 1.0 } else { 0.0 }).abs() > tolerance {
                return Err(invalid("orthonormal thermal material axes required"));
            }
        }
    }
    let determinant = if D == 2 {
        q[0][0] * q[1][1] - q[0][1] * q[1][0]
    } else {
        q[0][0] * (q[1][1] * q[2][2] - q[1][2] * q[2][1])
            - q[0][1] * (q[1][0] * q[2][2] - q[1][2] * q[2][0])
            + q[0][2] * (q[1][0] * q[2][1] - q[1][1] * q[2][0])
    };
    if (determinant - 1.0).abs() > tolerance {
        return Err(invalid("proper thermal material rotation required"));
    }
    Ok(())
}

fn invalid(reason: &'static str) -> ThermalTransportError {
    ThermalTransportError::Invalid(reason)
}
