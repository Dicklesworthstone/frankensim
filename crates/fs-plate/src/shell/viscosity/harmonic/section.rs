//! Original homogeneous and laminate sources for the shared harmonic shell.

use super::{
    Cx, Dd, MaterialPlane, PlaneStressHarmonic, PlateError, ReferenceDensityMoments, Resultants,
    TensorProny, TensorShellHarmonicError, bad, checkpoint, compensated_products, invalid,
    resultants,
};
use fs_material::laminate::{LaminateSection, LayerLaw};
use fs_math::c64::C64;

#[derive(Clone, Debug)]
enum Source {
    Homogeneous(TensorProny),
    Laminate {
        section: LaminateSection,
        // Only full elastic layers need a lossless tensor reduction. Original
        // Prony laws remain in the immutable source stack without being cloned.
        elastic: Vec<Option<Box<TensorProny>>>,
    },
}

/// An original full-tensor section in its facet's local `(ex, ey, ez)` frame.
/// A homogeneous section is centered; a laminate retains every layer and its
/// declared reference surface. This CST/DKT model has no shear deformation.
#[derive(Clone, Debug)]
pub struct TensorPronyShellSection {
    source: Source,
    thickness_m: f64,
    density_kg_m3: f64,
}

impl TensorPronyShellSection {
    /// Retain a centered original law, positive full thickness and density.
    /// Rotate a homogeneous full tensor with `TensorProny::rotated` beforehand.
    pub fn new(law: TensorProny, thickness_m: f64, density_kg_m3: f64) -> Result<Self, PlateError> {
        if !thickness_m.is_finite()
            || thickness_m <= 0.0
            || !density_kg_m3.is_finite()
            || density_kg_m3 <= 0.0
        {
            return Err(bad(
                "tensor shell requires positive finite thickness and reference density",
            ));
        }
        Ok(Self {
            source: Source::Homogeneous(law),
            thickness_m,
            density_kg_m3,
        })
    }

    /// Retain an oriented, possibly asymmetric or offset laminate unchanged.
    /// Each layer must supply an original full elastic or tensor-Prony law so
    /// all six material strains can be recovered and checked. Reduced-only
    /// plane-stress cards do not declare that full strain domain and are refused.
    ///
    /// The mesh is the declared reference surface, with its original normal.
    /// This section requires consistent mass to preserve signed `I1` and `I2`.
    pub fn from_laminate(section: LaminateSection) -> Result<Self, PlateError> {
        let thickness_m = section.total_thickness_m();
        let density_kg_m3 = section.i0_kg_m2() / thickness_m;
        if !density_kg_m3.is_finite() || density_kg_m3 <= 0.0 {
            return Err(bad(
                "laminate harmonic reference density is unrepresentable",
            ));
        }
        let mut elastic = Vec::with_capacity(section.layers().len());
        for layer in section.layers() {
            elastic.push(match &layer.law {
                LayerLaw::Prony(_) => None,
                LayerLaw::Elastic(law) => Some(Box::new(
                    TensorProny::from_tensor_shear(law.stiffness(), Vec::new(), law.strain_limit)
                        .map_err(|_| bad("laminate elastic layer has no admissible full strain law"))?,
                )),
                LayerLaw::PlaneStress { .. } => {
                    return Err(bad(
                        "harmonic laminate requires full tensor laws; a reduced-only layer has no full strain domain",
                    ));
                }
            });
        }
        ReferenceDensityMoments::from_laminate(&section)?;
        Ok(Self {
            source: Source::Laminate { section, elastic },
            thickness_m,
            density_kg_m3,
        })
    }

    /// Original homogeneous law, if this source is a centered single solid.
    /// A laminate has no single equivalent law; inspect `laminate()` instead.
    pub const fn law(&self) -> Option<&TensorProny> {
        match &self.source {
            Source::Homogeneous(law) => Some(law),
            Source::Laminate { .. } => None,
        }
    }

    /// Immutable original stack, including layer sources, angles and offsets.
    pub const fn laminate(&self) -> Option<&LaminateSection> {
        match &self.source {
            Source::Laminate { section, .. } => Some(section),
            Source::Homogeneous(_) => None,
        }
    }

    /// Full physical section thickness [m].
    pub const fn thickness_m(&self) -> f64 {
        self.thickness_m
    }

    /// Mean reference density [kg/m3]. Laminate inertia uses the original
    /// three moments, never a homogeneous replacement formed from this mean.
    pub const fn density_kg_m3(&self) -> f64 {
        self.density_kg_m3
    }

    pub(super) fn layer_count(&self) -> usize {
        self.laminate().map_or(1, |s| s.layers().len())
    }

    // Conservative fixed-size tensor reduction/derivative work, including
    // every retained relaxation branch. Admission is shared by reference and
    // prepared-frequency paths; a large stack cannot hide behind facet count.
    pub(super) fn material_work(&self) -> Option<usize> {
        (0..self.layer_count()).try_fold(0usize, |sum, layer| {
            self.layer_law(layer)
                .ok()?
                .branches()
                .len()
                .checked_mul(4096)
                .and_then(|v| v.checked_add(16_384))
                .and_then(|v| sum.checked_add(v))
        })
    }

    pub(super) fn density_moments(&self) -> Result<ReferenceDensityMoments, PlateError> {
        match self.laminate() {
            Some(section) => ReferenceDensityMoments::from_laminate(section),
            None => ReferenceDensityMoments::from_uniform(self.density_kg_m3),
        }
    }

    pub(super) fn layer_law(&self, layer: usize) -> Result<&TensorProny, TensorShellHarmonicError> {
        match &self.source {
            Source::Homogeneous(law) if layer == 0 => Ok(law),
            Source::Laminate { section, elastic } => match &section.layers()[layer].law {
                LayerLaw::Prony(law) => Ok(law),
                LayerLaw::Elastic(_) => elastic[layer]
                    .as_deref()
                    .ok_or_else(|| invalid("missing original elastic laminate reduction")),
                LayerLaw::PlaneStress { .. } => {
                    Err(invalid("unadmitted reduced-only laminate law"))
                }
            },
            Source::Homogeneous(_) => Err(invalid("homogeneous layer index")),
        }
    }

    fn integrate(
        &self,
        matrices: &[[f64; 9]],
        element: usize,
    ) -> Result<Resultants, TensorShellHarmonicError> {
        match self.laminate() {
            Some(section) => section
                .integrate_layer_matrices(matrices)
                .map_err(|source| TensorShellHarmonicError::Laminate { element, source }),
            None => Ok(resultants(matrices[0], self.thickness_m)?),
        }
    }

    pub(super) fn prepare(
        &self,
        element: usize,
        omega: f64,
        cx: &Cx<'_>,
    ) -> Result<(PreparedSection, Resultants, Resultants), TensorShellHarmonicError> {
        let count = self.layer_count();
        let mut layers = Vec::with_capacity(count);
        let mut storage = Vec::with_capacity(count);
        let mut loss = Vec::with_capacity(count);
        for layer in 0..count {
            checkpoint(cx)?;
            let material = self
                .layer_law(layer)?
                .prepare_plane_stress_harmonic(MaterialPlane::Axes12, omega)
                .map_err(|source| self.material_error(element, layer, source))?;
            storage.push(*material.storage_pa());
            loss.push(*material.loss_pa());
            let (z, to_principal) = if let Some(section) = self.laminate() {
                let geometry = &section.layer_integrals()[layer];
                (
                    [geometry.z_bottom_m, geometry.z_top_m],
                    section
                        .principal_strain_transform(layer)
                        .map_err(|source| TensorShellHarmonicError::Laminate { element, source })?,
                )
            } else {
                (
                    [-0.5 * self.thickness_m, 0.5 * self.thickness_m],
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                )
            };
            layers.push(PreparedLayer {
                material,
                z,
                to_principal,
            });
        }
        Ok((
            PreparedSection { layers },
            self.integrate(&storage, element)?,
            self.integrate(&loss, element)?,
        ))
    }

    pub(super) fn frequency_derivative(
        &self,
        prepared: &PreparedSection,
        element: usize,
        omega: f64,
        cx: &Cx<'_>,
    ) -> Result<(Resultants, Resultants), TensorShellHarmonicError> {
        let mut storage = Vec::with_capacity(self.layer_count());
        let mut viscosity = Vec::with_capacity(self.layer_count());
        for (layer, prepared) in prepared.layers.iter().enumerate() {
            checkpoint(cx)?;
            let (re, im) = self
                .layer_law(layer)?
                .plane_stress_complex_frequency_derivative(MaterialPlane::Axes12, omega)
                .map_err(|source| self.material_error(element, layer, source))?;
            let slope = core::array::from_fn(|i| {
                let omega = Dd::from_f64(omega);
                ((Dd::from_f64(im[i]) - Dd::from_f64(prepared.material.loss_pa()[i]) / omega)
                    / omega)
                    .to_f64()
            });
            storage.push(re);
            viscosity.push(slope);
        }
        Ok((
            self.integrate(&storage, element)?,
            self.integrate(&viscosity, element)?,
        ))
    }

    pub(super) fn material_error(
        &self,
        element: usize,
        layer: usize,
        source: fs_material::visco::ViscoError,
    ) -> TensorShellHarmonicError {
        match self.laminate() {
            Some(_) => TensorShellHarmonicError::Laminate {
                element,
                source: fs_material::laminate::LaminateError::Layer { layer, source },
            },
            None => TensorShellHarmonicError::Material { element, source },
        }
    }
}

pub(super) struct PreparedSection {
    layers: Vec<PreparedLayer>,
}

struct PreparedLayer {
    material: PlaneStressHarmonic,
    z: [f64; 2],
    to_principal: [[f64; 3]; 3],
}

impl PreparedSection {
    /// Numerical source payload plus temporary principal matrices; source
    /// laws are retained by the shell, not cloned for each prepared frequency.
    pub(super) fn payload_bytes(layers: usize) -> Option<usize> {
        layers
            .checked_mul(
                core::mem::size_of::<PreparedLayer>() + 2 * core::mem::size_of::<[f64; 9]>(),
            )
            .and_then(|v| v.checked_add(core::mem::size_of::<Self>()))
    }

    pub(super) fn layer_faces(&self, layer: usize) -> Option<[f64; 2]> {
        self.layers.get(layer).map(|prepared| prepared.z)
    }

    /// One shared original section-to-material strain recovery for both source
    /// admission and physical field observation. All transverse strains remain.
    pub(super) fn recover_strain(
        &self, layer: usize, strain: &[C64; 6], z: f64,
    ) -> Result<[C64; 6], fs_material::visco::ViscoError> {
        let prepared = self.layers.get(layer).ok_or(fs_material::visco::ViscoError::Parameters {
            what: "original shell material layer index is outside the section",
        })?;
        if !z.is_finite() || z < prepared.z[0] || z > prepared.z[1] {
            return Err(fs_material::visco::ViscoError::Parameters {
                what: "original shell observation must lie inside its material layer",
            });
        }
        let section_strain: [C64; 3] =
            core::array::from_fn(|r| strain[r] + strain[r + 3].scale(z));
        let principal = core::array::from_fn(|r| C64::new(
            compensated_products(
                (0..3).map(|i| (prepared.to_principal[r][i], section_strain[i].re)),
            ),
            compensated_products(
                (0..3).map(|i| (prepared.to_principal[r][i], section_strain[i].im)),
            ),
        ));
        prepared.material.recover_strain(&principal)
    }

    /// The original midpoint order is 01, 12, 20. Affine extrapolation gives
    /// the three vertices; their layer-face cycles bound this finite-element
    /// strain field by convexity, independently of integration-point values.
    pub(super) fn admit_facet(
        &self, section: &TensorPronyShellSection, element: usize,
        midpoint: &[[C64; 6]; 3], cx: &Cx<'_>,
    ) -> Result<(), TensorShellHarmonicError> {
        for vertex in 0..3 {
            checkpoint(cx)?;
            let strain: [C64; 6] = core::array::from_fn(|r|
                midpoint[vertex][r] + midpoint[(vertex + 2) % 3][r]
                    - midpoint[(vertex + 1) % 3][r]
            );
            self.admit(section, element, &strain, cx)?;
        }
        Ok(())
    }

    pub(super) fn admit(
        &self,
        section: &TensorPronyShellSection,
        element: usize,
        strain: &[C64; 6],
        cx: &Cx<'_>,
    ) -> Result<(), TensorShellHarmonicError> {
        for (layer, prepared) in self.layers.iter().enumerate() {
            checkpoint(cx)?;
            for z in prepared.z {
                self.recover_strain(layer, strain, z)
                    .map_err(|source| section.material_error(element, layer, source))?;
            }
        }
        Ok(())
    }
}
