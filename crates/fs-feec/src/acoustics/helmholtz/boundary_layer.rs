//! First-order thermoviscous boundary layers on the original H1 pressure space.
//!
//! With the parent's `exp(+i omega t)` convention, the Wentzell wall form is
//! `(-1+i)/2 [delta_v S + k^2 (gamma-1) delta_t C]`, where
//! `S_ij = integral(grad_T phi_i . grad_T phi_j)` and `C` is the original
//! boundary mass. This is the frequency-domain formulation of Berggren,
//! Bernland and Noreland, J. Comput. Phys. 371 (2018), 633--650,
//! <https://arxiv.org/abs/1801.04177>. It resolves tangential pressure variation;
//! replacing `S` by a scalar impedance would lose viscous dissipation.
//!
//! The gas/material owner supplies transport data. Geometry supplies the true
//! minimum curvature radius and half-gap, not the triangle diameter. A computed
//! frequency interval enforces thin layers relative to those scales and `1/k`.
//! That admission is an asymptotic condition, not a certified acoustic-error
//! bound. Faceted geometry, corner layers, moving walls, molecular relaxation,
//! finite wall heat capacity and narrow-gap LNS are not represented here.
//! Wall-patch edges have the natural zero tangential conormal derivative.

use super::{
    AcousticFeError, Coo, Csr, HelmholtzFe, HelmholtzModes, Medium, ReducedHelmholtz, cross3, dot3,
    sub3, tri_quadrature,
};
use std::collections::BTreeSet;

/// Gas transport quantities evaluated by the caller's material owner.
/// Density and sound speed come from the bound [`HelmholtzFe::medium`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoundaryLayerTransport {
    /// Dynamic shear viscosity [Pa s]. Zero disables the viscous wall term.
    pub dynamic_viscosity_pa_s: f64,
    /// Thermal conductivity [W/(m K)]. Zero disables the thermal wall term.
    pub thermal_conductivity_w_m_k: f64,
    /// Equilibrium isobaric specific heat [J/(kg K)].
    pub cp_j_kg_k: f64,
    /// Equilibrium `cp/cv`, at least one.
    pub heat_capacity_ratio: f64,
}

/// Physical geometry and a caller-selected thin-layer admission threshold.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoundaryLayerGeometry {
    /// Smallest smooth-wall curvature radius [m]; positive infinity for flat walls.
    pub minimum_curvature_radius_m: f64,
    /// Smallest half-gap to an opposing wall [m]; positive infinity if absent.
    pub minimum_half_gap_m: f64,
    /// Maximum `delta/geometry_scale` and `k delta`, in `(0, 0.05]`.
    /// This bounds expansion parameters, not the observable's numerical error.
    pub maximum_layer_ratio: f64,
}

/// Admitted layer thicknesses and expansion parameters at one frequency.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoundaryLayerScales {
    /// Viscous penetration depth `sqrt(2 mu/(rho omega))` [m].
    pub viscous_depth_m: f64,
    /// Thermal penetration depth `sqrt(2 kappa/(rho cp omega))` [m].
    pub thermal_depth_m: f64,
    /// Largest depth divided by the smaller physical geometry scale.
    pub geometry_ratio: f64,
    /// Largest depth multiplied by the acoustic wavenumber.
    pub wavelength_ratio: f64,
}

/// Mean cycle-averaged power transferred from acoustics to the stationary wall.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct BoundaryLayerPower {
    /// Viscous wall dissipation [W].
    pub viscous_w: f64,
    /// Thermal wall dissipation [W].
    pub thermal_w: f64,
}

fn invalid(what: &'static str) -> AcousticFeError {
    AcousticFeError::InvalidInput { what }
}

impl BoundaryLayerTransport {
    fn diffusivities(self, medium: Medium) -> Result<(f64, f64), AcousticFeError> {
        if !medium.density_kg_m3.is_finite()
            || medium.density_kg_m3 <= 0.0
            || !medium.sound_speed_m_s.is_finite()
            || medium.sound_speed_m_s <= 0.0
            || !self.dynamic_viscosity_pa_s.is_finite()
            || self.dynamic_viscosity_pa_s < 0.0
            || !self.thermal_conductivity_w_m_k.is_finite()
            || self.thermal_conductivity_w_m_k < 0.0
            || !self.cp_j_kg_k.is_finite()
            || self.cp_j_kg_k <= 0.0
            || !self.heat_capacity_ratio.is_finite()
            || self.heat_capacity_ratio < 1.0
        {
            return Err(invalid("boundary-layer gas transport"));
        }
        let viscous = self.dynamic_viscosity_pa_s / medium.density_kg_m3;
        let thermal = self.thermal_conductivity_w_m_k / medium.density_kg_m3 / self.cp_j_kg_k;
        if !viscous.is_finite() || !thermal.is_finite() {
            return Err(invalid("boundary-layer diffusivity overflow"));
        }
        Ok((viscous, thermal))
    }

    /// Compute the thin-layer admission interval [Hz]. Both endpoints follow
    /// from `delta <= ratio * min(radius,half_gap)` and `k delta <= ratio`.
    /// A zero-loss medium has `(0,+infinity)`. An empty band refuses.
    pub fn frequency_band_hz(
        self,
        medium: Medium,
        geometry: BoundaryLayerGeometry,
    ) -> Result<(f64, f64), AcousticFeError> {
        let (viscous, thermal) = self.diffusivities(medium)?;
        geometry.validate()?;
        let diffusion = viscous.max(thermal);
        if diffusion == 0.0 {
            return Ok((0.0, f64::INFINITY));
        }
        let scale = geometry
            .minimum_curvature_radius_m
            .min(geometry.minimum_half_gap_m);
        let ratio = geometry.maximum_layer_ratio;
        let pi = core::f64::consts::PI;
        let low = diffusion / pi / (scale * ratio).powi(2);
        let high = (ratio * medium.sound_speed_m_s).powi(2) / (4.0 * pi * diffusion);
        if low.is_nan() || !high.is_finite() || high <= 0.0 || low > high {
            return Err(invalid(
                "boundary-layer geometry has no admitted frequency band",
            ));
        }
        Ok((low, high))
    }

    /// Evaluate layer thicknesses and refuse a frequency outside the declared
    /// thin-layer geometry/wavelength band. No scalar damping is fitted.
    pub fn scales(
        self,
        medium: Medium,
        geometry: BoundaryLayerGeometry,
        frequency_hz: f64,
    ) -> Result<BoundaryLayerScales, AcousticFeError> {
        let (viscous, thermal) = self.diffusivities(medium)?;
        geometry.validate()?;
        if !frequency_hz.is_finite() || frequency_hz <= 0.0 {
            return Err(invalid("boundary-layer frequency"));
        }
        let omega = 2.0 * core::f64::consts::PI * frequency_hz;
        let viscous_depth_m = (2.0 * viscous / omega).sqrt();
        let thermal_depth_m = (2.0 * thermal / omega).sqrt();
        let deepest = viscous_depth_m.max(thermal_depth_m);
        let scale = geometry
            .minimum_curvature_radius_m
            .min(geometry.minimum_half_gap_m);
        let geometry_ratio = deepest / scale;
        let wavelength_ratio = omega / medium.sound_speed_m_s * deepest;
        if !omega.is_finite()
            || !geometry_ratio.is_finite()
            || !wavelength_ratio.is_finite()
            || geometry_ratio > geometry.maximum_layer_ratio
            || wavelength_ratio > geometry.maximum_layer_ratio
        {
            return Err(invalid(
                "frequency outside boundary-layer geometry/wavelength band; resolve thermoviscous flow",
            ));
        }
        Ok(BoundaryLayerScales {
            viscous_depth_m,
            thermal_depth_m,
            geometry_ratio,
            wavelength_ratio,
        })
    }
}

impl BoundaryLayerGeometry {
    fn validate(self) -> Result<(), AcousticFeError> {
        if self.minimum_curvature_radius_m.is_nan()
            || self.minimum_curvature_radius_m <= 0.0
            || self.minimum_half_gap_m.is_nan()
            || self.minimum_half_gap_m <= 0.0
            || !self.maximum_layer_ratio.is_finite()
            || self.maximum_layer_ratio <= 0.0
            || self.maximum_layer_ratio > 0.05
        {
            return Err(invalid(
                "boundary-layer physical geometry or thin-layer ratio",
            ));
        }
        Ok(())
    }
}

impl HelmholtzFe<'_> {
    /// Surface stiffness `integral grad_T(phi_i).grad_T(phi_j) dS` on selected
    /// original faces, using exact hierarchical derivatives and the original
    /// face quadrature. Normals/winding cannot change this bilinear form.
    /// Duplicate or nonboundary faces refuse. No volume normal derivative is
    /// retained: a basis function with a zero face trace has zero surface energy.
    pub fn boundary_tangential_stiffness(
        &self,
        faces: &[[u32; 3]],
    ) -> Result<Csr, AcousticFeError> {
        self.boundary_tangential_stiffness_controlled(faces, &|| false)
    }

    /// Surface assembly with cancellation checked before every face. A refusal
    /// publishes no partially assembled matrix.
    pub fn boundary_tangential_stiffness_controlled(
        &self,
        faces: &[[u32; 3]],
        cancelled: &dyn Fn() -> bool,
    ) -> Result<Csr, AcousticFeError> {
        if cancelled() {
            return Err(invalid("boundary tangential stiffness cancelled"));
        }
        let quad = tri_quadrature(self.order + 1);
        let mut coo = Coo::new(self.ndof(), self.ndof());
        let mut seen = BTreeSet::new();
        for &face in faces {
            if cancelled() {
                return Err(invalid("boundary tangential stiffness cancelled"));
            }
            let mut key = face;
            key.sort_unstable();
            if !seen.insert(key) {
                return Err(invalid("duplicate boundary-layer face"));
            }
            let rule = self.face_rule(face, &quad)?;
            let fi = self
                .space
                .complex
                .faces
                .binary_search(&key)
                .expect("face_rule admitted the canonical face");
            let (tet_index, omitted, _) = self.owners[fi];
            let tet = self.space.complex.tets[tet_index as usize];
            let local: Vec<usize> = (0..4).filter(|&i| i != omitted as usize).collect();
            let pa = self.positions[tet[local[0]] as usize];
            let e1 = sub3(self.positions[tet[local[1]] as usize], pa);
            let e2 = sub3(self.positions[tet[local[2]] as usize], pa);
            let twice_area = 2.0 * rule.area;
            let normal = cross3(e1, e2).map(|v| v / twice_area);
            let g1 = cross3(e2, normal).map(|v| v / twice_area);
            let g2 = cross3(normal, e1).map(|v| v / twice_area);
            let g0 = std::array::from_fn::<_, 3, _>(|i| -g1[i] - g2[i]);
            let tangent = [g0, g1, g2];
            let mut entries = Vec::new();
            for (global, function) in self.space.element_dofs(tet_index as usize) {
                let gradients: Vec<[f64; 3]> = quad
                    .iter()
                    .map(|(barycentric, _)| {
                        let mut lam = [0.0; 4];
                        for i in 0..3 {
                            lam[local[i]] = barycentric[i];
                        }
                        let derivative = function.d_lambda(lam, self.order);
                        std::array::from_fn(|axis| {
                            (0..3)
                                .map(|i| derivative[local[i]] * tangent[i][axis])
                                .sum()
                        })
                    })
                    .collect();
                if gradients.iter().any(|g| g.iter().any(|&x| x != 0.0)) {
                    entries.push((global, gradients));
                }
            }
            for (gi, a) in &entries {
                for (gj, b) in &entries {
                    let value: f64 = rule
                        .weights
                        .iter()
                        .zip(a.iter().zip(b))
                        .map(|(weight, (u, v))| weight * dot3(*u, *v))
                        .sum();
                    if !value.is_finite() {
                        return Err(invalid("boundary tangential stiffness overflow"));
                    }
                    coo.push(*gi, *gj, value);
                }
            }
        }
        if cancelled() {
            return Err(invalid("boundary tangential stiffness cancelled"));
        }
        Ok(coo.assemble())
    }
}

pub(super) struct ReducedBoundaryLayer {
    transport: BoundaryLayerTransport,
    geometry: BoundaryLayerGeometry,
    surface_stiffness: Vec<f64>,
    surface_mass: Vec<f64>,
}

impl ReducedBoundaryLayer {
    pub(super) fn new(
        fe: &HelmholtzFe<'_>,
        faces: &[[u32; 3]],
        ports: &[&[[u32; 3]]],
        transport: BoundaryLayerTransport,
        geometry: BoundaryLayerGeometry,
        project: &impl Fn(&Csr) -> Vec<f64>,
    ) -> Result<Self, AcousticFeError> {
        transport.frequency_band_hz(fe.medium, geometry)?;
        if faces.is_empty() {
            return Err(invalid("empty thermoviscous wall"));
        }
        let port_faces: BTreeSet<[u32; 3]> = ports
            .iter()
            .flat_map(|faces| faces.iter())
            .map(|face| {
                let mut key = *face;
                key.sort_unstable();
                key
            })
            .collect();
        for face in faces {
            let mut key = *face;
            key.sort_unstable();
            if port_faces.contains(&key) {
                return Err(invalid("thermoviscous wall overlaps an acoustic port"));
            }
        }
        let stiffness = fe.boundary_tangential_stiffness(faces)?;
        let mass = fe.boundary_mass(faces)?;
        let surface_stiffness = project(&stiffness);
        let surface_mass = project(&mass);
        if surface_stiffness
            .iter()
            .chain(&surface_mass)
            .any(|v| !v.is_finite())
        {
            return Err(invalid("boundary-layer projection overflow"));
        }
        Ok(Self {
            transport,
            geometry,
            surface_stiffness,
            surface_mass,
        })
    }

    fn coefficients(
        &self,
        medium: Medium,
        frequency_hz: f64,
    ) -> Result<(f64, f64), AcousticFeError> {
        let scales = self.transport.scales(medium, self.geometry, frequency_hz)?;
        let k = 2.0 * core::f64::consts::PI * frequency_hz / medium.sound_speed_m_s;
        let viscous = 0.5 * scales.viscous_depth_m;
        let thermal =
            0.5 * k * k * (self.transport.heat_capacity_ratio - 1.0) * scales.thermal_depth_m;
        if !viscous.is_finite() || !thermal.is_finite() {
            return Err(invalid("boundary-layer coefficient overflow"));
        }
        Ok((viscous, thermal))
    }

    pub(super) fn add_operator(
        &self,
        medium: Medium,
        frequency_hz: f64,
        operator: &mut [[f64; 2]],
    ) -> Result<(), AcousticFeError> {
        let (viscous, thermal) = self.coefficients(medium, frequency_hz)?;
        for ((entry, stiffness), mass) in operator
            .iter_mut()
            .zip(&self.surface_stiffness)
            .zip(&self.surface_mass)
        {
            let wall = viscous * stiffness + thermal * mass;
            entry[0] -= wall;
            entry[1] += wall;
            if entry.iter().any(|v| !v.is_finite()) {
                return Err(invalid("boundary-layer operator overflow"));
            }
        }
        Ok(())
    }
}

impl ReducedHelmholtz {
    /// Prepare the same volume/port Galerkin problem with viscous and isothermal
    /// wall boundary layers. Surface matrices are assembled and projected once;
    /// each original [`Self::response`] then adds two scalar-weighted matrices.
    ///
    /// Wall faces must be unique and disjoint from all acoustic ports. The
    /// source material owner supplies transport; geometric scale declarations
    /// must describe the real wall, not the mesh spacing. The original static
    /// enrichment and modal truncation remain controlled by the caller.
    pub fn new_with_boundary_layers(
        fe: &HelmholtzFe<'_>,
        modes: &HelmholtzModes,
        ports: &[&[[u32; 3]]],
        static_vectors: bool,
        wall_faces: &[[u32; 3]],
        transport: BoundaryLayerTransport,
        geometry: BoundaryLayerGeometry,
    ) -> Result<Self, AcousticFeError> {
        Self::build(
            fe,
            modes,
            ports,
            static_vectors,
            Some((wall_faces, transport, geometry)),
        )
    }

    /// Dissipation [W] in the wall for these reduced pressure coordinates.
    /// `Im(p* A_wall p)/(2 omega rho)` gives the physical input-work loss;
    /// viscous and thermal contributions are reported separately. The supplied
    /// coordinates must belong to this reduction, normally `response.amplitudes`.
    pub fn boundary_layer_power(
        &self,
        frequency_hz: f64,
        amplitudes: &[[f64; 2]],
    ) -> Result<BoundaryLayerPower, AcousticFeError> {
        if !frequency_hz.is_finite()
            || frequency_hz <= 0.0
            || amplitudes.len() != self.n
            || amplitudes.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "boundary-layer power frequency or coordinate shape",
            ));
        }
        let Some(wall) = &self.boundary_layer else {
            return Ok(BoundaryLayerPower::default());
        };
        let (viscous, thermal) = wall.coefficients(self.medium, frequency_hz)?;
        let quadratic = |matrix: &[f64]| -> Result<f64, AcousticFeError> {
            let mut value = 0.0;
            let mut absolute = 0.0;
            for (i, a) in amplitudes.iter().enumerate() {
                for (j, b) in amplitudes.iter().enumerate() {
                    let term = matrix[i * self.n + j] * (a[0] * b[0] + a[1] * b[1]);
                    value += term;
                    absolute += term.abs();
                }
            }
            if !value.is_finite()
                || !absolute.is_finite()
                || value < -128.0 * f64::EPSILON * absolute
            {
                return Err(invalid("boundary-layer power overflow or nonpassive form"));
            }
            Ok(value)
        };
        let scale = 1.0 / (4.0 * core::f64::consts::PI * frequency_hz * self.medium.density_kg_m3);
        let power = BoundaryLayerPower {
            viscous_w: scale * viscous * quadratic(&wall.surface_stiffness)?,
            thermal_w: scale * thermal * quadratic(&wall.surface_mass)?,
        };
        if !power.viscous_w.is_finite() || !power.thermal_w.is_finite() {
            return Err(invalid("boundary-layer power overflow"));
        }
        Ok(power)
    }
}

#[cfg(test)]
mod tests;
