//! Fixed, arbitrarily oriented disks with three-dimensional actuator wash.
//!
//! This is an explicit model choice: unlike the original spanwise gate, the
//! radial distance is measured in each disk's actual plane. The retained
//! axial factor and finite 2R window are unchanged. This is still an axial,
//! quasi-steady BEMT model, not a blade-resolved crossflow or skew-wake model.

use super::{
    CoupledStep, IncidentInflow, InfluenceOperator, Panel, PropCouplingSolverSpec, PropDisk,
    Refusal, Rotor, SectionClosure, StripSpec, coupled_step_impl, hash_domain, wash_factor,
};

/// Prepared positive-thrust directions in the same frame as the panels/disks.
/// Spin-axis signs and blade handedness are separate mechanical inputs.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DiskAxes {
    directions: [[f64; 3]; 2],
}

fn invalid(message: &str) -> Refusal {
    Refusal {
        code: "spatial-disk-invalid",
        message: message.into(),
        ranked_repairs: vec!["supply finite geometry and explicit unit thrust directions".into()],
    }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    let xy = a[0] * b[0] + a[1] * b[1];
    xy + a[2] * b[2]
}

impl DiskAxes {
    /// Admit two unit directions without silently normalizing them.
    ///
    /// # Errors
    /// `spatial-disk-invalid` for nonfinite or nonunit directions.
    pub fn new(directions: [[f64; 3]; 2]) -> Result<Self, Refusal> {
        if directions.iter().any(|axis| {
            !axis.iter().all(|x| x.is_finite()) || (dot(*axis, *axis) - 1.0).abs() > 1e-12
        }) {
            return Err(invalid("disk directions must be finite unit vectors"));
        }
        Ok(Self { directions })
    }

    /// Admitted positive-thrust axes, left then right.
    #[must_use]
    pub fn directions(&self) -> [[f64; 3]; 2] {
        self.directions
    }

    /// Spatial wash factors at a physical point. Radial distance includes
    /// both transverse directions; no vertical overlap is inferred from span.
    ///
    /// # Errors
    /// `spatial-disk-invalid` for unrepresentable geometry or invalid radius.
    pub fn wash_factors(
        &self,
        point: [f64; 3],
        centers: [[f64; 3]; 2],
        radius_m: f64,
    ) -> Result<[f64; 2], Refusal> {
        if !radius_m.is_finite()
            || radius_m <= 0.0
            || !point
                .iter()
                .chain(centers.iter().flatten())
                .all(|x| x.is_finite())
        {
            return Err(invalid("disk radius and sample geometry are invalid"));
        }
        let mut factors = [0.0; 2];
        for k in 0..2 {
            let offset = core::array::from_fn(|j| point[j] - centers[k][j]);
            let axial = dot(offset, self.directions[k]);
            let radial: [f64; 3] =
                core::array::from_fn(|j| offset[j] - axial * self.directions[k][j]);
            let radial_sq = dot(radial, radial);
            if !axial.is_finite() || !radial_sq.is_finite() {
                return Err(invalid("disk-relative geometry overflow"));
            }
            // Scale first: the dimensionless law avoids overflow in dx²+R²
            // for large but representable geometry. It is the same 2R model.
            factors[k] = wash_factor(axial / radius_m, radial_sq.sqrt() / radius_m, 1.0);
        }
        Ok(factors)
    }

    pub(super) fn wash_map(
        &self,
        strips: &[StripSpec],
        panels: &[Panel],
        disks: &[PropDisk; 2],
        radius: f64,
    ) -> Result<Vec<[f64; 2]>, Refusal> {
        let centers = [disks[0].center_m, disks[1].center_m];
        strips
            .iter()
            .map(|strip| {
                let panel = &panels[strip.panel_indices[0]];
                let point = core::array::from_fn(|j| 0.5 * panel.a[j] + 0.5 * panel.b[j]);
                self.wash_factors(point, centers, radius)
            })
            .collect()
    }

    pub(super) fn velocity(
        &self,
        incident: [f64; 3],
        factors: [f64; 2],
        induction: [f64; 2],
    ) -> [f64; 3] {
        let left = factors[0] * induction[0];
        let right = factors[1] * induction[1];
        core::array::from_fn(|j| {
            let increment = left * self.directions[0][j] + right * self.directions[1][j];
            incident[j] + increment
        })
    }

    pub(super) fn axial(&self, disk: usize, incident: [f64; 3], induced: [f64; 3]) -> f64 {
        dot(
            core::array::from_fn(|j| incident[j] + induced[j]),
            self.directions[disk],
        )
    }

    pub(super) fn digest(&self, spec: &PropCouplingSolverSpec) -> String {
        let mut bytes = spec.digest().into_bytes();
        for x in self.directions.iter().flatten() {
            bytes.extend_from_slice(&x.to_bits().to_le_bytes());
        }
        hash_domain("org.frankensim.fs-flyer.spatial-disk-coupling.v1", &bytes).to_hex()
    }
}

/// Native two-way wing/BEMT solve for spatially mounted disks.
///
/// The local incident air plus native wing-induced velocity is projected
/// onto each positive-thrust axis. The scalar momentum root supplies a
/// vector wash along that axis at every actual 3D strip location. The two
/// disks may have different axes and spin rates. Only disk-average axial
/// inflow is modeled: transverse disk loading, blade harmonics, swirl, wake
/// contraction, and direct disk-to-disk interference remain outside this tier.
///
/// # Errors
/// Native wing/BEMT/coupling refusals, `spatial-disk-invalid` for geometry,
/// or `disk-inflow-invalid` for reversed/nonfinite projected throughflow.
#[allow(clippy::too_many_arguments)]
pub fn coupled_step(
    axes: &DiskAxes,
    op: &InfluenceOperator,
    panels: &[Panel],
    strips: &[StripSpec],
    closure: SectionClosure<'_>,
    rotor: &Rotor,
    disks: &[PropDisk; 2],
    freestream: [f64; 3],
    rho: f64,
    spec: &PropCouplingSolverSpec,
    warm: Option<[f64; 2]>,
    incident: IncidentInflow<'_>,
) -> Result<CoupledStep, Refusal> {
    if incident.strip_velocity_mps.len() != strips.len() {
        return Err(invalid("provide exactly one incident velocity per strip"));
    }
    coupled_step_impl(
        op,
        panels,
        strips,
        closure,
        rotor,
        disks,
        freestream,
        rho,
        spec,
        warm,
        Some(incident),
        Some(*axes),
    )
}
