//! One source-bound thin-membrane law for discrete and continuous ducts. Above thickness/radius 0.05
//! an isotropic wall takes Lame's thick-cylinder law instead (`lame`); an orthotropic wall is refused there.
//! Geometry reduction is shared with the existing graph wall. Frequency-domain
//! memory uses GeneralizedMaxwell::modulus; no fitted loss factor is introduced.
use super::{
    AxialConstraint, WallMaterialLoss, WallPatch, YOUNG_VISCOSITY_PROPERTY, invalid, property,
};
use crate::acoustic_realize::AcousticRealizeError;
use fs_material::state_point::{
    DENSITY_PROPERTY, POISSON_RATIO_PROPERTY, ResolvedMaterialStatePoint, YOUNG_MODULUS_PROPERTY,
};
use fs_math::c64::C64;
use fs_phs::{WallPin, wall_specific_impedance};
use fs_qty::{Density, Dims, DynViscosity, Pressure};
use fs_vfit::relaxation::RelaxationImpedanceSpec;

/// Largest thickness/radius for which the membrane law is used; above it an isotropic wall takes
/// Lame's thick cylinder.
const MEMBRANE_THICKNESS_OVER_RADIUS: f64 = 0.05;

mod update;
mod thermal;
mod lame;
mod orthotropic;
pub(super) mod finite_gas;

/// Immutable material, thickness, restraint and complete local loss law.
///
/// The frozen-state, small-displacement membrane may be isotropic or have
/// orthotropic principal axes aligned with the hoop and axial directions. It does not include axial shell waves, pressure prestress, thermal
/// evolution, or nonlocal wall radiation. At thickness/radius <= 0.05 the membrane law runs; above it
/// an isotropic wall takes Lame's thick-cylinder law and an orthotropic wall is refused.
#[derive(Clone, Debug, PartialEq)]
pub struct LocalMaterialWall {
    material: ResolvedMaterialStatePoint,
    thickness_m: f64,
    axial: AxialConstraint,
    loss: WallMaterialLoss,
    band_hz: [f64; 2],
    density: f64,
    young: f64,
    viscosity: f64,
    axial_factor: f64,
    orthotropic_axes: Option<[usize; 2]>,
    proportional_viscosity_time_property: Option<String>,
}

/// The same source-bound wall and its native time-domain acoustic owner.
/// Original source receipts and applicability remain attached to the operator.
#[derive(Clone, Debug)]
pub struct MaterialWallDg {
    wall: LocalMaterialWall,
    duct: fs_duct::dg::DgDuct,
}
impl MaterialWallDg {
    /// Original source, thickness, constraint and declared applicability band.
    #[must_use]
    pub const fn wall(&self) -> &LocalMaterialWall {
        &self.wall
    }
    /// Original DG owner, including reciprocal wall mechanics and all gas memory.
    #[must_use]
    pub const fn duct(&self) -> &fs_duct::dg::DgDuct {
        &self.duct
    }
    /// Move the original coupled operator into its acoustic host without
    /// rebuilding or copying its material and gas-memory coordinates.
    #[must_use]
    pub fn into_duct(self) -> fs_duct::dg::DgDuct {
        self.duct
    }
}

/// Independently sourced wall materials on the ORIGINAL DG elements.
/// The exact material receipts, temperatures, thicknesses and restraints remain
/// attached to the one gas/wall time operator; no effective material is minted.
#[derive(Clone, Debug)]
pub struct RegionalMaterialWallDg {
    walls: Vec<LocalMaterialWall>,
    band_hz: [f64; 2],
    duct: fs_duct::dg::DgDuct,
}
impl RegionalMaterialWallDg {
    /// Bind one source law per original element, in element order. Distinct
    /// material states may describe different substances or temperatures.
    /// Material interfaces must coincide with element boundaries; refinement
    /// repeats each region's source over its sub-elements explicitly.
    ///
    /// Each region preserves its mass, equilibrium stiffness, viscosity and
    /// every positive Maxwell arm in the native reciprocal DG operator.
    /// Spectra may differ in branch count and relaxation times; no common fit
    /// replaces them. Exactly homogeneous sources retain the original uniform path.
    ///
    /// The returned applicability band is the intersection of all source bands,
    /// not a filter or error certificate. Original cylindrical, frozen-state and
    /// cold-binding restrictions apply; the wall law is the membrane's up to thickness/radius 0.05
    /// and Lame's thick cylinder above it (isotropic walls only).
    pub fn new(
        walls: Vec<LocalMaterialWall>,
        duct: fs_duct::dg::DgDuct,
    ) -> Result<Self, AcousticRealizeError> {
        Self::new_controlled(walls, duct, &mut || false)
    }

    /// Cancellable source binding; failures publish no partially coupled duct.
    pub fn new_controlled(
        walls: Vec<LocalMaterialWall>,
        duct: fs_duct::dg::DgDuct,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AcousticRealizeError> {
        if cancelled() {
            return Err(AcousticRealizeError::Cancelled);
        }
        if walls.is_empty() || walls.len() != duct.mesh().elements {
            return Err(invalid(
                "regional DG wall needs one original material per element",
            ));
        }
        let mut band_hz = walls[0].band_hz;
        for wall in &walls {
            if cancelled() {
                return Err(AcousticRealizeError::Cancelled);
            }
            band_hz[0] = band_hz[0].max(wall.band_hz[0]);
            band_hz[1] = band_hz[1].min(wall.band_hz[1]);
        }
        if band_hz[0] >= band_hz[1] {
            return Err(invalid(
                "regional wall sources have no common positive frequency interval",
            ));
        }
        if walls.iter().all(|wall| *wall == walls[0]) {
            let bound = walls[0].clone().bind_dg_controlled(duct, cancelled)?;
            return Ok(Self {
                walls,
                band_hz,
                duct: bound.duct,
            });
        }
        let radius = duct.uniform_circular_radius_m().ok_or_else(|| {
            invalid("regional DG material wall requires an explicit uniform circular geometry")
        })?;
        let nodes = duct.node_positions();
        let local = duct.mesh().order + 1;
        let mut pins = Vec::with_capacity(walls.len());
        for (element, wall) in walls.iter().enumerate() {
            if cancelled() {
                return Err(AcousticRealizeError::Cancelled);
            }
            let base = element * local;
            let length = nodes[base + local - 1] - nodes[base];
            pins.push(wall.base_patch(radius, length)?.wall);
        }
        let spectra: Vec<_> = walls.iter().map(|wall| match &wall.loss {
            WallMaterialLoss::Maxwell(source) => Some(source.material()),
            _ => None,
        }).collect();
        let duct = duct
            .with_compliant_wall_relaxation_regions_controlled(&pins, &spectra, cancelled)
            .map_err(|e| match e {
                fs_duct::DuctError::Cancelled { .. } => AcousticRealizeError::Cancelled,
                error => AcousticRealizeError::Nonlinear(error.to_string()),
            })?;
        Ok(Self {
            walls,
            band_hz,
            duct,
        })
    }

    /// Original material/geometry laws, in original element order.
    #[must_use]
    pub fn walls(&self) -> &[LocalMaterialWall] {
        &self.walls
    }
    /// Common applicability interval of all original source laws [Hz].
    #[must_use]
    pub const fn band_hz(&self) -> [f64; 2] {
        self.band_hz
    }
    /// The native gas/wall operator, including retained gas fields and memory.
    #[must_use]
    pub const fn duct(&self) -> &fs_duct::dg::DgDuct {
        &self.duct
    }
    /// Move the coupled operator into the acoustic host without rebuilding its
    /// regional material, gas-memory or hereditary wall coordinates.
    #[must_use]
    pub fn into_duct(self) -> fs_duct::dg::DgDuct {
        self.duct
    }
}

impl LocalMaterialWall {
    /// Resolve the same original claims used by graph-based material walls.
    ///
    /// # Errors
    /// Missing/wrong-quantity/out-of-band material data, unsupported material
    /// dependence, nonpositive geometry or invalid isotropic coefficients.
    pub fn new(
        material: ResolvedMaterialStatePoint,
        thickness_m: f64,
        axial: AxialConstraint,
        loss: WallMaterialLoss,
        band_hz: [f64; 2],
    ) -> Result<Self, AcousticRealizeError> {
        if !thickness_m.is_finite()
            || thickness_m <= 0.0
            || !band_hz.iter().all(|f| f.is_finite() && *f > 0.0)
            || band_hz[0] >= band_hz[1]
            || !(core::f64::consts::TAU * band_hz[1]).is_finite()
        {
            return Err(invalid(
                "local material wall requires positive thickness and an ordered positive frequency band",
            ));
        }
        let density = property(&material, DENSITY_PROPERTY, Density::DIMS, band_hz)?;
        let young = property(&material, YOUNG_MODULUS_PROPERTY, Pressure::DIMS, band_hz)?;
        let nu = property(&material, POISSON_RATIO_PROPERTY, Dims([0; 6]), band_hz)?;
        if density <= 0.0 || young <= 0.0 || nu <= -1.0 || nu >= 0.5 {
            return Err(invalid("isotropic wall needs rho,E>0 and -1<nu<0.5"));
        }
        let viscosity = match &loss {
            WallMaterialLoss::Elastic | WallMaterialLoss::Maxwell(_) => 0.0,
            WallMaterialLoss::KelvinVoigt => property(
                &material,
                YOUNG_VISCOSITY_PROPERTY,
                DynViscosity::DIMS,
                band_hz,
            )?,
        };
        if viscosity < 0.0 {
            return Err(invalid("wall Young viscosity must be nonnegative"));
        }
        if let WallMaterialLoss::Maxwell(spectrum) = &loss {
            spectrum.validate_source(&material)?;
            spectrum.validate_material(young, band_hz)?;
        }
        let axial_factor = match axial {
            AxialConstraint::StressFree => 1.0,
            AxialConstraint::StrainFixed => 1.0 / (1.0 - nu * nu),
        };
        Ok(Self {
            material,
            thickness_m,
            axial,
            loss,
            band_hz,
            density,
            young,
            viscosity,
            axial_factor,
            orthotropic_axes: None,
            proportional_viscosity_time_property: None,
        })
    }

    /// Original material identity, queries and property receipts.
    #[must_use]
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }
    #[must_use]
    pub const fn thickness_m(&self) -> f64 {
        self.thickness_m
    }
    #[must_use]
    pub const fn axial_constraint(&self) -> AxialConstraint {
        self.axial
    }
    #[must_use]
    pub const fn loss(&self) -> &WallMaterialLoss {
        &self.loss
    }
    #[must_use]
    pub const fn band_hz(&self) -> [f64; 2] {
        self.band_hz
    }

    /// Bind this exact material reduction to the native DG time-domain duct.
    /// The full lateral wall responds to gas pressure; its velocity feeds back
    /// into continuity, and its mass/elastic energy is retained in DgState.
    /// Previously attached gas viscothermal memory remains a separate loss.
    ///
    /// Elastic and Kelvin-Voigt coefficients follow base_patch without a fit.
    /// Maxwell walls retain every positive source arm in the SAME DG operator;
    /// the native exact exponential advances their history together with gas
    /// and wall motion. Midpoint-network clock limits do not apply to that
    /// propagator. The original pointwise gas-field masses also couple to this
    /// wall; its inertia and spring mass use the geometric Gram matrix, never
    /// an averaged gas. The uniform-cylinder restriction applies; the wall law is the membrane's
    /// up to thickness/radius 0.05 and Lame's thick cylinder above it (isotropic walls only). This is
    /// cold compilation. `MaterialWallDg::rebind_controlled` updates an already
    /// vibrating wall and its certified propagator together; Maxwell histories
    /// require complete original source-arm correspondence on fixed geometry.
    /// The retained band is the SAME model applicability statement as the
    /// graph-wall owner, not a filter or broadband input certificate. The
    /// caller must resolve and qualify its excited frequency range separately.
    pub fn bind_dg(
        self,
        duct: fs_duct::dg::DgDuct,
    ) -> Result<MaterialWallDg, AcousticRealizeError> {
        self.bind_dg_controlled(duct, &mut || false)
    }

    /// Cancellable cold binding. Original gas/wall source refusals remain
    /// authoritative; interruption publishes no partially coupled operator.
    pub fn bind_dg_controlled(
        self,
        duct: fs_duct::dg::DgDuct,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MaterialWallDg, AcousticRealizeError> {
        if cancelled() {
            return Err(AcousticRealizeError::Cancelled);
        }
        let radius = duct.uniform_circular_radius_m().ok_or_else(|| {
            invalid("DG material wall requires an explicit uniform circular geometry")
        })?;
        let patch = self.base_patch(radius, duct.total_length())?;
        let duct = match &self.loss {
            WallMaterialLoss::Maxwell(source) => duct.with_compliant_wall_relaxation_controlled(
                patch.wall,
                source.material(),
                cancelled,
            ),
            _ => duct.with_compliant_wall_controlled(patch.wall, cancelled),
        }
        .map_err(|e| match e {
            fs_duct::DuctError::Cancelled { .. } => AcousticRealizeError::Cancelled,
            error => AcousticRealizeError::Nonlinear(error.to_string()),
        })?;
        Ok(MaterialWallDg { wall: self, duct })
    }

    /// Base mass/equilibrium/viscosity image for an actual cylindrical strip.
    /// Hereditary arms remain in `loss()`; this is not its complete impedance.
    ///
    /// At thickness/radius <= 0.05 the membrane expressions run. Above it an isotropic wall takes
    /// Lame's thick cylinder (`lame`, with the Poisson ratio of its own resolved material) and an
    /// orthotropic wall is refused. Non-positive or non-finite geometry is refused.
    pub fn base_patch(&self, radius: f64, length: f64) -> Result<WallPatch, AcousticRealizeError> {
        let h = self.thickness_m;
        if ![radius, length, h]
            .iter()
            .all(|x| x.is_finite() && *x > 0.0)
        {
            return Err(invalid(
                "local wall requires positive finite radius, length and thickness",
            ));
        }
        // (stiffness / E, surface density) from either law, in `lame::thick_cylinder`'s order.
        let (factor, surface_density) = if h / radius <= MEMBRANE_THICKNESS_OVER_RADIUS {
            ((h / radius) / radius * self.axial_factor, self.density * h)
        } else if self.orthotropic_axes.is_some() {
            return Err(invalid(
                "orthotropic wall requires thickness/radius <= 0.05; the thick-wall law is isotropic",
            ));
        } else {
            let nu = property(&self.material, POISSON_RATIO_PROPERTY, Dims([0; 6]), self.band_hz)?;
            lame::thick_cylinder(radius, h, self.density, nu, self.axial)
        };
        let patch = WallPatch {
            area_m2: core::f64::consts::TAU * radius * length,
            wall: WallPin {
                surface_density,
                stiffness_per_area: self.young * factor,
                resistance: self.viscosity * factor,
            },
        };
        patch.impedance()?;
        if self.viscosity > 0.0 && patch.wall.resistance == 0.0 {
            return Err(invalid(
                "wall material reduction lost its nonzero viscosity",
            ));
        }
        Ok(patch)
    }

    // Same source reduction as graph/DG/TMM walls, evaluated at the original
    // frustum owner's canonical radius and wetted strip length. Keep native
    // source refusals here rather than translating them through a numeric API.
    pub(crate) fn acoustic_cell_wall(&self, radius: f64, strip_length: f64)
        -> Result<fs_phs::acoustic_outlet::AcousticCellWall, AcousticRealizeError>
    {
        use fs_phs::acoustic_outlet::{AcousticCellWall, AcousticWallRelaxation};
        let patch = self.base_patch(radius, strip_length)?;
        self.load(patch)?; // Original positive-arm geometry/source admission.
        let relaxation = match &self.loss {
            WallMaterialLoss::Maxwell(source) => source.material().terms.iter()
                .filter(|(modulus, _)| *modulus > 0.0)
                .map(|&(modulus, relaxation_time_s)| AcousticWallRelaxation {
                    stiffness_per_area: (patch.wall.stiffness_per_area / self.young) * modulus,
                    relaxation_time_s,
                }).collect(),
            _ => Vec::new(),
        };
        Ok(AcousticCellWall { wall: patch.wall, relaxation })
    }

    pub(super) fn load(
        &self,
        patch: WallPatch,
    ) -> Result<RelaxationImpedanceSpec, AcousticRealizeError> {
        match &self.loss {
            WallMaterialLoss::Maxwell(spectrum) => spectrum.load(patch, self.young, self.band_hz),
            _ => RelaxationImpedanceSpec::new(patch.impedance()?, &[])
                .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string())),
        }
    }

    /// Complete surface impedance [Pa s/m] under exp(-i omega t).
    /// Storage and loss moduli come from the same source law as the graph's
    /// hereditary states, evaluated at the requested physical frequency.
    pub fn specific_impedance(&self, radius: f64, omega: f64) -> Result<C64, AcousticRealizeError> {
        if !omega.is_finite()
            || omega < core::f64::consts::TAU * self.band_hz[0]
            || omega > core::f64::consts::TAU * self.band_hz[1]
        {
            return Err(invalid(
                "material wall frequency lies outside its declared source band",
            ));
        }
        let patch = self.base_patch(radius, 1.0)?;
        let mut pin = patch.wall;
        if let WallMaterialLoss::Maxwell(spectrum) = &self.loss {
            // Retain the graph owner's admission of every positive physical arm.
            spectrum.load(patch, self.young, self.band_hz)?;
            // Reuse the complete generalized-Maxwell constitutive owner.
            let (storage, loss) = spectrum.material().modulus(omega);
            let factor = pin.stiffness_per_area / self.young;
            pin.stiffness_per_area = factor * storage;
            pin.resistance += factor * loss / omega;
        }
        wall_specific_impedance(&pin, omega)
            .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))
    }
}
