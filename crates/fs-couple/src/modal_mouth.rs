//! Geometry-derived exterior loading of every retained circular-duct mode.
//!
//! A real circular mouth is a cap on the original closed exterior mesh. The
//! native BEM owner receives the signed velocity traces `psi_n / S`, where
//! `S = pi R^2`. Its transpose pressure projection is therefore already the
//! duct's complete `P = Z U` matrix in Pa s/m^3, under `exp(-i omega t)`.
//! Higher modes can have zero net volume flow and still radiate and react on
//! the bore. No division by their net flow, independent piston substitution,
//! matrix diagonalization, fit, or passivity repair occurs here.
//!
//! This is a harmonic rigid-exterior model. The mouth uses the existing m=0
//! radial basis and panel-centroid traces. Omitted azimuthal/radial modes,
//! mesh/basis convergence, the chosen BEM route's validity and interior wall
//! losses remain the respective physical owners' responsibilities.
//! [`realization`] converts explicitly supported source samples to the existing
//! complete causal MIMO load through the shared physical-source fitting owner.

/// Source-bound rational memory for all signed modal mouth coordinates.
pub mod realization;

use fs_bem::helmholtz::aperture::{
    Aperture, ApertureError, ApertureImpedance, ApertureRoute, aperture_impedance,
};
use fs_bem::helmholtz::{HelmholtzError, Medium};
use fs_bem::panel3d::SpherePanels;
use fs_bem::radiation_policy::GeometryPolicy;
use fs_duct::DuctError;
use fs_duct::modal::{AxisymmetricMouthBasis, ModalMouthImpedance, ModeCount};
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use fs_math::dd::Dd;

/// A circular port on a closed, outward-oriented exterior triangle surface.
#[derive(Clone, Debug)]
pub struct CircularMouthPatch {
    /// Centre of the physical mouth reference plane [m].
    pub centre_m: [f64; 3],
    /// Unit normal pointing from the duct into the exterior.
    pub outward_normal: [f64; 3],
    /// Radius of the original circular modal reference section [m].
    pub radius_m: f64,
    /// Unique original cap panels. Other surface panels remain rigid.
    pub panels: Vec<usize>,
}

/// Whole-binding and whole-sweep ceilings, additional to the BEM route's limits.
#[derive(Clone, Copy, Debug)]
pub struct ModalMouthBemBudget {
    /// Original exterior panel count, before geometry or trace allocation.
    pub maximum_panels: usize,
    /// Number of real coefficients in the complete modes-by-panels trace bank.
    pub maximum_trace_entries: usize,
    /// Total complex matrix entries over the requested frequency sweep. Both
    /// the unchanged BEM artifact and ready-to-use duct matrices are retained.
    pub maximum_matrix_entries: usize,
    /// Bound on `panels * modes^2 * (1 + frequencies)`, for geometry Gram and
    /// pressure projection. This excludes BEM kernel/factorization work, which
    /// remains under the explicitly selected original route's own limits.
    pub maximum_projection_terms: u64,
    /// Maximum infinity norm of `sum A_i psi_i psi_i^T / S - I`, in (0,1).
    /// This checks that the actual cap resolves the declared modal traces.
    /// It is not an acoustic solve or geometry-error certificate. The original
    /// shapes and areas are never renormalized to make this diagnostic pass.
    pub maximum_basis_gram_defect: f64,
}

/// Geometry, original radiation solve, duct load admission or cancellation.
#[derive(Debug)]
pub enum ModalMouthBemError {
    /// Incomplete physical binding, unresolved cap basis or exceeded budget.
    Invalid(&'static str),
    /// Original closed-surface validation or geometry refusal.
    Geometry(HelmholtzError),
    /// Original selected BEM route's refusal, without fallback.
    Radiation(ApertureError),
    /// Original modal basis or full Hermitian load admission refusal.
    Duct(DuctError),
    /// Cancellation before a complete binding or sweep was published.
    Cancelled,
}
impl core::fmt::Display for ModalMouthBemError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "modal mouth exterior: {self:?}")
    }
}
impl std::error::Error for ModalMouthBemError {}
impl From<DuctError> for ModalMouthBemError {
    fn from(error: DuctError) -> Self {
        Self::Duct(error)
    }
}
fn bad(what: &'static str) -> ModalMouthBemError {
    ModalMouthBemError::Invalid(what)
}
fn poll(gate: &CancelGate) -> Result<(), ModalMouthBemError> {
    if gate.is_requested() {
        Err(ModalMouthBemError::Cancelled)
    } else {
        Ok(())
    }
}

/// Immutable exterior geometry, gas and complete original modal trace bank.
pub struct BemModalMouth<'surface> {
    surface: &'surface SpherePanels,
    patch: CircularMouthPatch,
    exterior: GasState,
    basis: AxisymmetricMouthBasis,
    apertures: Vec<Aperture>,
    basis_gram_defect: f64,
    cap_area_ratio: f64,
    budget: ModalMouthBemBudget,
}

/// A complete admitted sweep and its unchanged native radiation source data.
pub struct BemModalMouthBake {
    loads: Vec<ModalMouthImpedance>,
    source: ApertureImpedance,
    basis_gram_defect: f64,
    cap_area_ratio: f64,
}
impl BemModalMouthBake {
    /// Frequency/basis/radius/gas-bound duct loads, in requested frequency order.
    pub fn loads(&self) -> &[ModalMouthImpedance] {
        &self.loads
    }
    /// Original full BEM matrices, route, mesh/shape declaration and diagnostics.
    /// `volume_gains` are modal net-flow diagnostics, not rescaling factors.
    pub fn source(&self) -> &ApertureImpedance {
        &self.source
    }
    /// Original discrete modal Gram defect, unchanged by radiation admission.
    pub fn basis_gram_defect(&self) -> f64 {
        self.basis_gram_defect
    }
    /// Original cap area divided by the declared circular section area.
    pub fn cap_area_ratio(&self) -> f64 {
        self.cap_area_ratio
    }
}

impl<'surface> BemModalMouth<'surface> {
    /// Bind a circular cap and every signed m=0 mode without changing geometry.
    /// Retained triangles must form an outward closed surface; cap triangles
    /// must lie on the declared plane, face outward and stay inside its circle.
    /// Polygonal boundary and panel quadrature error remain in the reported
    /// area ratio and Gram defect. The caller supplies the physical cap's
    /// coverage; these finite modal moments do not certify a disk tessellation.
    pub fn new(
        surface: &'surface SpherePanels,
        patch: CircularMouthPatch,
        exterior: GasState,
        modes: ModeCount,
        budget: ModalMouthBemBudget,
        gate: &CancelGate,
    ) -> Result<Self, ModalMouthBemError> {
        poll(gate)?;
        modes.validate()?;
        exterior
            .validate_coherence()
            .map_err(|_| bad("coherent exterior gas required"))?;
        let panels = surface.areas().len();
        let count = modes.0;
        let entries = panels
            .checked_mul(count)
            .ok_or_else(|| bad("modal trace size overflow"))?;
        let terms = (panels as u128) * (count as u128) * (count as u128);
        if panels == 0
            || panels > budget.maximum_panels
            || entries > budget.maximum_trace_entries
            || terms > u128::from(budget.maximum_projection_terms)
            || !budget.maximum_basis_gram_defect.is_finite()
            || budget.maximum_basis_gram_defect <= 0.0
            || budget.maximum_basis_gram_defect >= 1.0
        {
            return Err(bad(
                "complete modal mouth geometry/basis exceeds its declared limits",
            ));
        }
        let area = core::f64::consts::PI * patch.radius_m * patch.radius_m;
        let normal_square: f64 = patch.outward_normal.iter().map(|v| v * v).sum();
        if !patch.radius_m.is_finite()
            || patch.radius_m <= 0.0
            || !area.is_finite()
            || area <= 0.0
            || !area.recip().is_finite()
            || patch.centre_m.iter().any(|v| !v.is_finite())
            || !normal_square.is_finite()
            || (normal_square - 1.0).abs() > 1e-12
            || patch.panels.is_empty()
            || patch.panels.len() > panels
        {
            return Err(bad(
                "finite positive circular mouth, unit normal and cap panels required",
            ));
        }
        // Reuse the original exact-coordinate topology and orientation owner.
        // This validates the surface; it does not select or substitute a route.
        GeometryPolicy::new(surface).map_err(ModalMouthBemError::Geometry)?;
        poll(gate)?;
        let triangles = surface
            .triangles()
            .ok_or_else(|| bad("retained cap triangles required"))?;
        let basis = AxisymmetricMouthBasis::new(modes)?;
        let mut seen = vec![false; panels];
        let mut apertures: Vec<_> = (0..count)
            .map(|_| Aperture {
                shape: vec![0.0; panels],
                reference_point: patch.centre_m,
                normal: patch.outward_normal,
            })
            .collect();
        let mut gram = vec![Dd::ZERO; count * count];
        let mut cap_area = Dd::ZERO;
        let mut values = vec![0.0; count];
        for &panel in &patch.panels {
            poll(gate)?;
            if panel >= panels || seen[panel] {
                return Err(bad("unique original cap panel indices required"));
            }
            seen[panel] = true;
            let normal = surface.normals()[panel];
            if normal
                .iter()
                .zip(patch.outward_normal)
                .any(|(a, b)| (a - b).abs() > 1e-12)
            {
                return Err(bad("cap normal must point outward along the duct axis"));
            }
            for point in triangles[panel] {
                let _ = radius_fraction(point, &patch)?;
            }
            let radial = radius_fraction(surface.centroids()[panel], &patch)?;
            let weight = surface.areas()[panel] / area;
            cap_area = cap_area + Dd::from_f64(surface.areas()[panel]);
            for m in 0..count {
                values[m] = basis.shape(m, radial)?;
                let trace = values[m] / area;
                if !trace.is_finite() || (values[m] != 0.0 && trace == 0.0) {
                    return Err(bad("physical modal velocity trace is not representable"));
                }
                apertures[m].shape[panel] = trace;
            }
            for a in 0..count {
                for b in 0..count {
                    gram[a * count + b] = gram[a * count + b]
                        + Dd::from_f64(weight) * Dd::from_f64(values[a]) * Dd::from_f64(values[b]);
                }
            }
        }
        let basis_gram_defect = (0..count)
            .map(|a| {
                (0..count)
                    .map(|b| (gram[a * count + b].to_f64() - if a == b { 1.0 } else { 0.0 }).abs())
                    .sum::<f64>()
            })
            .fold(0.0, f64::max);
        let cap_area_ratio = cap_area.to_f64() / area;
        if !basis_gram_defect.is_finite()
            || !cap_area_ratio.is_finite()
            || basis_gram_defect > budget.maximum_basis_gram_defect
        {
            return Err(bad(
                "original cap quadrature does not resolve the complete modal basis; refine the mesh",
            ));
        }
        poll(gate)?;
        Ok(Self {
            surface,
            patch,
            exterior,
            basis,
            apertures,
            basis_gram_defect,
            cap_area_ratio,
            budget,
        })
    }

    /// Original physical patch, including its declared circular reference plane.
    pub fn patch(&self) -> &CircularMouthPatch {
        &self.patch
    }
    /// Fixed original duct basis; higher-mode net volume gain can vanish.
    pub fn basis(&self) -> &AxisymmetricMouthBasis {
        &self.basis
    }
    /// Original immutable exterior gas, distinct from an interior gas if needed.
    pub fn exterior(&self) -> &GasState {
        &self.exterior
    }
    /// Signed outward normal velocity per unit modal flow, in surface order.
    pub fn velocity_trace(&self, mode: usize) -> Option<&[f64]> {
        self.apertures.get(mode).map(|a| a.shape.as_slice())
    }
    /// Original discrete modal Gram defect, without basis repair.
    pub fn basis_gram_defect(&self) -> f64 {
        self.basis_gram_defect
    }
    /// Original cap area divided by the circular section area.
    pub fn cap_area_ratio(&self) -> f64 {
        self.cap_area_ratio
    }

    /// Solve the original exterior at every requested positive increasing
    /// frequency. The native batch shares its factor across all mouth modes;
    /// Galerkin also reuses geometry across the sweep. Its restricted band or
    /// any original solver refusal propagates. Full matrix load admission is
    /// pointwise Hermitian PSD, not a causal or broadband passivity certificate.
    /// No partially completed sweep is returned.
    pub fn bake(
        &self,
        omegas_rad_s: &[f64],
        route: &ApertureRoute,
        gate: &CancelGate,
    ) -> Result<BemModalMouthBake, ModalMouthBemError> {
        poll(gate)?;
        let n = self.basis.modes().0 as u128;
        let frequencies = omegas_rad_s.len() as u128;
        let entries = 2 * n * n * frequencies;
        let terms = (self.surface.areas().len() as u128) * n * n * (1 + frequencies);
        if omegas_rad_s.is_empty()
            || omegas_rad_s.iter().any(|w| !w.is_finite() || *w <= 0.0)
            || omegas_rad_s.windows(2).any(|w| w[1] <= w[0])
            || entries > self.budget.maximum_matrix_entries as u128
            || terms > u128::from(self.budget.maximum_projection_terms)
        {
            return Err(bad(
                "positive increasing complete mouth sweep within declared limits required",
            ));
        }
        let medium = Medium {
            density: self.exterior.density,
            sound_speed: self.exterior.sound_speed,
        };
        let source = aperture_impedance(
            self.surface,
            &self.apertures,
            omegas_rad_s,
            medium,
            route,
            gate,
        )
        .map_err(ModalMouthBemError::Radiation)?;
        let mut loads = Vec::with_capacity(source.rows.len());
        for row in &source.rows {
            poll(gate)?;
            // With b_n=psi_n/S, native transpose projection is exactly P_n/U_m.
            // `.acoustic()` would divide by net-flow gains and destroy this basis.
            loads.push(ModalMouthImpedance::new(
                row.omega,
                self.patch.radius_m,
                self.exterior,
                self.basis.modes(),
                row.impedance.clone(),
            )?);
        }
        poll(gate)?;
        Ok(BemModalMouthBake {
            loads,
            source,
            basis_gram_defect: self.basis_gram_defect,
            cap_area_ratio: self.cap_area_ratio,
        })
    }
}

fn radius_fraction(point: [f64; 3], patch: &CircularMouthPatch) -> Result<f64, ModalMouthBemError> {
    let delta: [f64; 3] = std::array::from_fn(|i| point[i] - patch.centre_m[i]);
    let axial: f64 = delta
        .iter()
        .zip(patch.outward_normal)
        .map(|(a, n)| a * n)
        .sum();
    let coordinate = point
        .into_iter()
        .chain(patch.centre_m)
        .map(f64::abs)
        .fold(patch.radius_m, f64::max);
    let tolerance = 128.0 * f64::EPSILON * coordinate;
    if !tolerance.is_finite()
        || tolerance > 1e-6 * patch.radius_m
        || !axial.is_finite()
        || axial.abs() > tolerance
    {
        return Err(bad(
            "mouth cap must resolve the original circular reference plane",
        ));
    }
    let radial = delta[0].hypot(delta[1]).hypot(delta[2]);
    if !radial.is_finite() || radial > patch.radius_m + tolerance {
        return Err(bad(
            "mouth cap vertex is outside the declared circular section",
        ));
    }
    // A vertex at the circle can lie one ulp outside. This only admits its
    // geometry; centroid samples must lie inside and are never clipped.
    let fraction = radial / patch.radius_m;
    if !fraction.is_finite() {
        return Err(bad("mouth radius ratio is not representable"));
    }
    Ok(fraction)
}

#[cfg(test)]
mod tests;
