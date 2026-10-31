//! Immutable incident boundary data shared by prepared and iterative solves.
use super::*;

/// Pressure and outward-normal velocity of a field regular in the bodies.
/// The binding prevents reuse at another surface, frequency or medium. It
/// does not certify that externally supplied data satisfy Helmholtz's equation.
#[derive(Clone, Debug)]
pub struct IncidentTraces {
    fingerprint: u64,
    k: f64,
    medium: Medium,
    pressure_pa: Vec<C64>,
    normal_velocity_m_s: Vec<C64>,
}

impl IncidentTraces {
    /// Sample coherent plane-wave and monopole sources on the original mesh.
    /// Clearance prevents singular placement; it is not a discretization error
    /// bound. Nearby source variations still require a spatial mesh ladder.
    pub fn from_sources(
        surface: &SpherePanels,
        k: f64,
        medium: Medium,
        sources: &[IncidentSource],
        minimum_source_clearance_m: f64,
    ) -> Result<Self, HelmholtzError> {
        Self::from_field(
            surface,
            k,
            medium,
            IncidentField::Sources(sources),
            minimum_source_clearance_m,
        )
    }

    /// Bind an analytic incident field to the original collocation points.
    /// Source admission uses component-wise winding and point-to-triangle
    /// clearance throughout the existing iterative surface range.
    pub fn from_field(
        surface: &SpherePanels,
        k: f64,
        medium: Medium,
        field: IncidentField<'_>,
        minimum_source_clearance_m: f64,
    ) -> Result<Self, HelmholtzError> {
        validate_frequency(k, medium)?;
        field.validate()?;
        if !minimum_source_clearance_m.is_finite() || minimum_source_clearance_m <= 0.0 {
            return Err(bad("incident source clearance must be positive and finite"));
        }
        validate_closed(surface)?;
        let positions: Vec<_> = match field {
            IncidentField::PlaneWave(_) => Vec::new(),
            IncidentField::Monopoles(sources) => sources.iter().map(|s| s.position_m).collect(),
            IncidentField::Sources(sources) => sources
                .iter()
                .filter_map(|s| match s {
                    IncidentSource::Monopole(s) => Some(s.position_m),
                    IncidentSource::PlaneWave(_) => None,
                })
                .collect(),
        };
        if !positions.is_empty() {
            near_field::admit_exterior_points(
                surface,
                &positions,
                minimum_source_clearance_m,
                crate::panel3d::MAX_SURFACE_PANELS,
            )?;
        }
        Self::from_admitted_field(surface, k, medium, field, None)
    }

    /// The immutable source/surface geometry must already have passed the
    /// original exterior admission. Reuse that admission over a frequency
    /// grid; retain the original field sampler and frequency/shape checks.
    pub(crate) fn from_admitted_field(
        surface: &SpherePanels,
        k: f64,
        medium: Medium,
        field: IncidentField<'_>,
        gate: Option<&fs_exec::CancelGate>,
    ) -> Result<Self, HelmholtzError> {
        validate_frequency(k, medium)?;
        field.validate()?;
        let mut pressure_pa = Vec::with_capacity(surface.areas().len());
        let mut normal_velocity_m_s = Vec::with_capacity(surface.areas().len());
        for (index, (&point, normal)) in surface
            .centroids()
            .iter()
            .zip(surface.normals())
            .enumerate()
        {
            if index % 256 == 0 && gate.is_some_and(fs_exec::CancelGate::is_requested) {
                return Err(HelmholtzError::Cancelled);
            }
            let sample = field.sample(k, medium, point)?;
            pressure_pa.push(sample.pressure_pa);
            normal_velocity_m_s.push((0..3).fold(C64::ZERO, |v, a| {
                v + sample.velocity_m_s[a].scale(normal[a])
            }));
        }
        if normal_velocity_m_s.iter().any(|&v| !finite(v)) {
            return Err(bad("nonfinite incident normal velocity"));
        }
        Ok(Self {
            fingerprint: surface_fingerprint(surface),
            k,
            medium,
            pressure_pa,
            normal_velocity_m_s,
        })
    }

    /// Bind traces transferred from a caller-owned Helmholtz field. Pressure
    /// is [Pa], velocity is [m/s] along the body's outward normal, and both use
    /// exp(-i omega t). The producer owns field regularity in every component,
    /// the Helmholtz equation and the accuracy of its boundary transfer.
    pub fn from_traces(
        surface: &SpherePanels,
        k: f64,
        medium: Medium,
        pressure_pa: &[C64],
        normal_velocity_m_s: &[C64],
    ) -> Result<Self, HelmholtzError> {
        validate_frequency(k, medium)?;
        validate_closed(surface)?;
        let n = surface.areas().len();
        if pressure_pa.len() != n || normal_velocity_m_s.len() != n {
            return Err(HelmholtzError::ShapeMismatch {
                what: "one incident pressure and normal velocity per panel",
            });
        }
        if pressure_pa
            .iter()
            .chain(normal_velocity_m_s)
            .any(|&v| !finite(v))
        {
            return Err(bad("nonfinite incident boundary traces"));
        }
        Ok(Self {
            fingerprint: surface_fingerprint(surface),
            k,
            medium,
            pressure_pa: pressure_pa.to_vec(),
            normal_velocity_m_s: normal_velocity_m_s.to_vec(),
        })
    }

    /// Incident pressure [Pa] in original panel order.
    pub fn pressure_pa(&self) -> &[C64] {
        &self.pressure_pa
    }

    /// Incident velocity [m/s] along the original outward panel normals.
    pub fn normal_velocity_m_s(&self) -> &[C64] {
        &self.normal_velocity_m_s
    }

    pub(crate) fn validate_binding(
        &self,
        surface: &SpherePanels,
        k: f64,
        medium: Medium,
    ) -> Result<(), HelmholtzError> {
        if self.fingerprint != surface_fingerprint(surface) || self.k != k || self.medium != medium
        {
            return Err(bad(
                "incident traces belong to a different surface, frequency or medium",
            ));
        }
        Ok(())
    }

    pub(crate) fn drive(&self, law: Admittance<'_>) -> Result<Vec<C64>, HelmholtzError> {
        let image = law.image(&self.pressure_pa)?;
        let inverse_rho_c = 1.0 / (self.medium.density * self.medium.sound_speed);
        let drive: Vec<_> = self
            .normal_velocity_m_s
            .iter()
            .zip(image)
            .map(|(&v, y)| -v - y.scale(inverse_rho_c))
            .collect();
        if drive.iter().any(|&v| !finite(v)) {
            return Err(bad("nonfinite incident Robin drive"));
        }
        Ok(drive)
    }

    pub(super) fn into_boundary_traces(self) -> (Vec<C64>, Vec<C64>) {
        (self.pressure_pa, self.normal_velocity_m_s)
    }
}
