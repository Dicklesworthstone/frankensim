//! Source-declared hardness apparatus and loading conditions. This is a
//! measurement selector, not a conversion between empirical hardness scales.

use fs_qty::{Dims, QtyAny};

use crate::{MatDbError, ObservationId};

/// One ordered constant-force hold in a source-declared hardness test.
/// Ramps and instrument compliance remain part of the named protocol.
#[derive(Debug, Clone, PartialEq)]
pub struct HardnessLoadStep {
    force_newtons: f64,
    dwell_seconds: f64,
}

impl HardnessLoadStep {
    /// Normalize at the source adapter first. Force must be positive and dwell
    /// nonnegative; both must be finite and carry their actual dimensions.
    pub fn new(force: QtyAny, dwell: QtyAny) -> Result<Self, MatDbError> {
        if force.dims != Dims([1, 1, -2, 0, 0, 0]) || !force.value.is_finite() || force.value <= 0.0
        {
            return Err(invalid("load must be a finite positive force in newtons"));
        }
        if dwell.dims != Dims([0, 0, 1, 0, 0, 0]) || !dwell.value.is_finite() || dwell.value < 0.0 {
            return Err(invalid(
                "dwell must be a finite nonnegative time in seconds",
            ));
        }
        Ok(Self {
            force_newtons: force.value,
            dwell_seconds: if dwell.value == 0.0 { 0.0 } else { dwell.value },
        })
    }

    /// Applied force in coherent SI.
    #[must_use]
    pub const fn force_newtons(&self) -> f64 {
        self.force_newtons
    }

    /// Duration at this force in coherent SI.
    #[must_use]
    pub const fn dwell_seconds(&self) -> f64 {
        self.dwell_seconds
    }
}

/// Exact apparatus, loading sequence, protocol and observed specimen identity.
/// The observation supplies the existing specimen/process, method, raw artifact,
/// caveats and licensed provenance. Names are source declarations: this type
/// does not certify compliance with a standard or invent missing conditions.
#[derive(Debug, Clone, PartialEq)]
pub struct HardnessTestContext {
    indenter: String,
    loading: Vec<HardnessLoadStep>,
    protocol: String,
    observation: ObservationId,
}

impl HardnessTestContext {
    /// Bound loading history, including preliminary and final holds when the
    /// source provides them. Indenter must identify its geometry/material;
    /// protocol must identify the procedure and revision being requested.
    pub fn new(
        indenter: impl Into<String>,
        loading: Vec<HardnessLoadStep>,
        protocol: impl Into<String>,
        observation: ObservationId,
    ) -> Result<Self, MatDbError> {
        let indenter = indenter.into();
        let protocol = protocol.into();
        if indenter.trim().is_empty() || protocol.trim().is_empty() {
            return Err(invalid(
                "indenter and protocol declarations must be nonblank",
            ));
        }
        if indenter.len() > 4096 || protocol.len() > 4096 {
            return Err(invalid(
                "indenter and protocol declarations exceed 4096 bytes",
            ));
        }
        if loading.is_empty() || loading.len() > 64 {
            return Err(invalid(
                "loading history must contain between 1 and 64 holds",
            ));
        }
        Ok(Self {
            indenter,
            loading,
            protocol,
            observation,
        })
    }

    /// Source-declared indenter geometry and material or exact specification.
    #[must_use]
    pub fn indenter(&self) -> &str {
        &self.indenter
    }

    /// Ordered force/dwell holds; no implicit averaging or reordering.
    #[must_use]
    pub fn loading(&self) -> &[HardnessLoadStep] {
        &self.loading
    }

    /// Exact source-declared procedure/revision, not inferred from the scale.
    #[must_use]
    pub fn protocol(&self) -> &str {
        &self.protocol
    }

    /// Existing observation identity binding specimen and test provenance.
    #[must_use]
    pub const fn observation(&self) -> ObservationId {
        self.observation
    }

    pub(crate) fn canonical_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::new();
        for text in [&self.indenter, &self.protocol] {
            bytes.extend_from_slice(&(text.len() as u64).to_le_bytes());
            bytes.extend_from_slice(text.as_bytes());
        }
        bytes.extend_from_slice(&self.observation.0.0);
        bytes.extend_from_slice(&(self.loading.len() as u64).to_le_bytes());
        for step in &self.loading {
            bytes.extend_from_slice(&step.force_newtons.to_bits().to_le_bytes());
            bytes.extend_from_slice(&step.dwell_seconds.to_bits().to_le_bytes());
        }
        bytes
    }
}

fn invalid(reason: &'static str) -> MatDbError {
    MatDbError::InvalidHardnessContext { reason }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g0_hardness_context_requires_dimensions_finite_holds_and_named_apparatus() {
        let force = |v| QtyAny::new(v, Dims([1, 1, -2, 0, 0, 0]));
        let time = |v| QtyAny::new(v, Dims([0, 0, 1, 0, 0, 0]));
        for value in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(HardnessLoadStep::new(force(value), time(10.0)).is_err());
        }
        for value in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(HardnessLoadStep::new(force(20.0), time(value)).is_err());
        }
        assert!(HardnessLoadStep::new(time(20.0), time(10.0)).is_err());
        assert!(HardnessLoadStep::new(force(20.0), force(10.0)).is_err());
        assert_eq!(
            HardnessLoadStep::new(force(20.0), time(-0.0))
                .unwrap()
                .dwell_seconds()
                .to_bits(),
            0
        );
        let step = HardnessLoadStep::new(force(20.0), time(10.0)).unwrap();
        let observation = ObservationId(fs_blake3::hash_domain("test", b"synthetic"));
        for (indenter, protocol, loading) in [
            (" ", "fixture revision 1", vec![step.clone()]),
            ("fixture pyramid", "", vec![step.clone()]),
            ("fixture pyramid", "fixture revision 1", vec![]),
            (
                "fixture pyramid",
                "fixture revision 1",
                vec![step.clone(); 65],
            ),
        ] {
            assert!(HardnessTestContext::new(indenter, loading, protocol, observation).is_err());
        }
        let context = HardnessTestContext::new(
            "fixture pyramid",
            vec![step],
            "fixture revision 1",
            observation,
        )
        .unwrap();
        assert!(
            crate::PropertyKey::new("pressure", Dims([-1, 1, -2, 0, 0, 0]))
                .with_hardness_test(context)
                .is_err()
        );
    }
}
