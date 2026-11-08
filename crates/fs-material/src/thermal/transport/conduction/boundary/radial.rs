//! Original caloric material faces with explicit radial resistance geometry.
//!
//! Physical radii and the selected angular footprint remain distinct from the
//! resistance-equivalent distance used by the native finite-cell heat owner.

use super::super::CaloricConductionLayer;
use super::IntegratedCaloricBoundary;
use crate::thermal::conduction::{
    LayerConductionConfig, ReservoirLayerGeometry, ReservoirLayerResponse,
};
use fs_thermal::enthalpy::conduction::CurrentConductiveBoundary;
use fs_thermal::enthalpy::conduction::implicit::boundary::{
    ImplicitBoundaryGeometry, ImplicitBoundaryLaw, ImplicitBoundaryState,
};
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};

/// An explicit radial material trace, boundary face and angular footprint.
///
/// Cylinders have insulated ends; their complete axial length is supplied.
/// Spheres use a fraction of the full solid angle. The same fraction scales
/// both the physical face area and the exact geometric conductance multiplier;
/// lateral spreading and the remaining surface channels stay outside this law.
/// No radius, extent, material volume or mass is inferred from an equivalent
/// distance, and a point source at zero radius is not admitted.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RadialCaloricBoundaryGeometry {
    coordinates: RadialCoordinates,
    area_fraction: f64,
    flux: ReservoirLayerGeometry,
    equivalent_distance_m: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum RadialCoordinates {
    Cylindrical {
        surface_radius_m: f64,
        trace_radius_m: f64,
        length_m: f64,
    },
    Spherical {
        surface_radius_m: f64,
        trace_radius_m: f64,
    },
}

impl RadialCaloricBoundaryGeometry {
    /// Either cylindrical face and its material trace, with explicit axial
    /// length and circumferential area fraction in (0, 1].
    pub fn cylindrical(
        surface_radius_m: f64,
        trace_radius_m: f64,
        length_m: f64,
        area_fraction: f64,
    ) -> Result<Self, DepositionError> {
        let flux =
            ReservoirLayerGeometry::cylindrical_face(surface_radius_m, trace_radius_m, length_m)?
                .with_area_fraction(area_fraction)?;
        Self::prepare(
            RadialCoordinates::Cylindrical {
                surface_radius_m,
                trace_radius_m,
                length_m,
            },
            area_fraction,
            flux,
        )
    }

    /// Either spherical face and its material trace, with the explicitly
    /// represented solid-angle fraction Omega/(4*pi) in (0, 1].
    pub fn spherical(
        surface_radius_m: f64,
        trace_radius_m: f64,
        area_fraction: f64,
    ) -> Result<Self, DepositionError> {
        let flux = ReservoirLayerGeometry::spherical_face(surface_radius_m, trace_radius_m)?
            .with_area_fraction(area_fraction)?;
        Self::prepare(
            RadialCoordinates::Spherical {
                surface_radius_m,
                trace_radius_m,
            },
            area_fraction,
            flux,
        )
    }

    fn prepare(
        coordinates: RadialCoordinates,
        area_fraction: f64,
        flux: ReservoirLayerGeometry,
    ) -> Result<Self, DepositionError> {
        let equivalent_distance_m = flux.surface_area_m2() / flux.conductance_multiplier_m();
        if !equivalent_distance_m.is_finite() || equivalent_distance_m <= 0.0 {
            return Err(DepositionError::Arithmetic);
        }
        Ok(Self {
            coordinates,
            area_fraction,
            flux,
            equivalent_distance_m,
        })
    }

    /// Actual selected radial face area. It never supplies material inventory.
    #[must_use]
    pub const fn surface_area_m2(self) -> f64 {
        self.flux.surface_area_m2()
    }

    /// Exact C = 1/integral(dr/A(r)), in metres, for Q = C integral(k(T) dT).
    #[must_use]
    pub const fn conductance_multiplier_m(self) -> f64 {
        self.flux.conductance_multiplier_m()
    }

    /// A_face * integral(dr/A(r)), in metres. This resistance coordinate is
    /// generally different from the physical radial trace-to-face thickness.
    #[must_use]
    pub const fn equivalent_distance_m(self) -> f64 {
        self.equivalent_distance_m
    }

    #[must_use]
    pub const fn surface_radius_m(self) -> f64 {
        match self.coordinates {
            RadialCoordinates::Cylindrical {
                surface_radius_m, ..
            }
            | RadialCoordinates::Spherical {
                surface_radius_m, ..
            } => surface_radius_m,
        }
    }

    #[must_use]
    pub const fn trace_radius_m(self) -> f64 {
        match self.coordinates {
            RadialCoordinates::Cylindrical { trace_radius_m, .. }
            | RadialCoordinates::Spherical { trace_radius_m, .. } => trace_radius_m,
        }
    }

    /// Supplied complete axial length for a cylinder; absent for a sphere.
    #[must_use]
    pub const fn cylinder_length_m(self) -> Option<f64> {
        match self.coordinates {
            RadialCoordinates::Cylindrical { length_m, .. } => Some(length_m),
            RadialCoordinates::Spherical { .. } => None,
        }
    }

    #[must_use]
    pub const fn area_fraction(self) -> f64 {
        self.area_fraction
    }

    fn same_radial_path(self, reference: Self) -> bool {
        matches!(
            (self.coordinates, reference.coordinates),
            (
                RadialCoordinates::Cylindrical { .. },
                RadialCoordinates::Cylindrical { .. }
            ) | (
                RadialCoordinates::Spherical { .. },
                RadialCoordinates::Spherical { .. }
            )
        ) && self.area_fraction == reference.area_fraction
            && (self.surface_radius_m() < self.trace_radius_m())
                == (reference.surface_radius_m() < reference.trace_radius_m())
    }
}

/// A source-integrated radial material face and its one original areal film.
///
/// The original caloric object and k(T) integral supply the material response.
/// The native boundary receives the exact radial resistance coordinate and
/// physical footprint, retaining all original heat, entropy, reference mass and
/// accepted time. A finite film may keep the actual material surface supported
/// even when its bath is outside the material temperature chart.
///
/// Geometry and the original source context are fixed during each thermal
/// solve. This adds no half-cell storage, pressure work, inferred expansion,
/// radiation or independently frozen bulk conductivity. Density-driven planar
/// thickness queries intentionally remain unsupported on this radial law.
pub struct IntegratedRadialCaloricBoundary<'l, 'a, 'm> {
    geometry: RadialCaloricBoundaryGeometry,
    prepared: IntegratedCaloricBoundary<'l, 'a, 'm>,
}

impl<'l, 'a, 'm> IntegratedRadialCaloricBoundary<'l, 'a, 'm> {
    pub fn new(
        layer: &'l CaloricConductionLayer<'a, 'm>,
        geometry: RadialCaloricBoundaryGeometry,
        surface_resistance_m2_k_w: f64,
        config: LayerConductionConfig,
    ) -> Result<Self, DepositionError> {
        let native = ImplicitBoundaryGeometry {
            area_m2: geometry.surface_area_m2(),
            half_distance_m: geometry.equivalent_distance_m(),
            surface_resistance_m2_k_w,
        };
        // Reuse original controls, native-distance admission and R''/A once.
        // No face query occurs during this preparation. Retain the exact radial
        // multiplier itself, avoiding a second A/d reconstruction in the flux.
        let mut prepared = IntegratedCaloricBoundary::new(layer, native, config)?;
        prepared.flux_geometry = geometry.flux;
        Ok(Self { geometry, prepared })
    }

    #[must_use]
    pub const fn geometry(&self) -> RadialCaloricBoundaryGeometry {
        self.geometry
    }

    /// Exact area and resistance-equivalent distance for the original observer.
    /// Original reference cell volume and mass are supplied independently.
    #[must_use]
    pub const fn native_geometry(&self) -> ImplicitBoundaryGeometry {
        self.prepared.geometry
    }

    /// Actual material-face temperature, signed outward heat and original
    /// endpoint-temperature derivatives at this fixed physical geometry.
    pub fn response_at_specific_energy(
        &self,
        specific_energy_j_kg: f64,
        bath_temperature_k: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ReservoirLayerResponse, DepositionError> {
        self.prepared.response_at_specific_energy(
            specific_energy_j_kg,
            bath_temperature_k,
            checkpoint,
        )
    }
}

impl ImplicitBoundaryLaw for IntegratedRadialCaloricBoundary<'_, '_, '_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.prepared.caloric_law()
    }

    fn state_at_specific_enthalpy(
        &self,
        specific_energy_j_kg: f64,
        bath_temperature_k: f64,
        geometry: ImplicitBoundaryGeometry,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitBoundaryState, DepositionError> {
        self.prepared.state_at_specific_enthalpy(
            specific_energy_j_kg,
            bath_temperature_k,
            geometry,
            checkpoint,
        )
    }
}

/// One explicitly selected current radial configuration on the original face.
///
/// The original radial family, inner/outer material side and angular fraction
/// remain fixed. New physical radii and cylinder length must be supplied; no
/// mechanical extent is inferred from the native resistance coordinates.
/// The original source, reference coverage, mass and areal film retain their
/// own owners. Both the current radial multiplier and R''/A_current are prepared
/// once. Source-compatible mechanics and its work require independent closure.
pub struct CurrentRadialCaloricBoundary<'r, 'l, 'a, 'm> {
    reference: &'r IntegratedRadialCaloricBoundary<'l, 'a, 'm>,
    prepared: IntegratedRadialCaloricBoundary<'l, 'a, 'm>,
}

impl<'l, 'a, 'm> IntegratedRadialCaloricBoundary<'l, 'a, 'm> {
    /// Rebind explicitly supplied radial coordinates without changing the
    /// original boundary. The native observer must receive current_geometry()
    /// on the SAME original reservoir or exposed-ambient channel.
    ///
    /// Preparation is constant work with no material-table copies or source
    /// queries. The native implicit attachment still probes original source
    /// admission at its actual cell/bath state. Every thermal trial holds this
    /// configuration fixed; a density-driven thickness tangent is not inferred.
    pub fn with_current_radial_geometry<'r>(
        &'r self,
        current: RadialCaloricBoundaryGeometry,
    ) -> Result<CurrentRadialCaloricBoundary<'r, 'l, 'a, 'm>, DepositionError> {
        if !current.same_radial_path(self.geometry) {
            return Err(DepositionError::InvalidInput(
                "current radial boundary requires its original family, material side and angular fraction",
            ));
        }
        let prepared = Self::new(
            self.prepared.layer,
            current,
            self.prepared.geometry.surface_resistance_m2_k_w,
            self.prepared.config,
        )?;
        Ok(CurrentRadialCaloricBoundary {
            reference: self,
            prepared,
        })
    }
}

impl CurrentRadialCaloricBoundary<'_, '_, '_, '_> {
    #[must_use]
    pub const fn reference_geometry(&self) -> RadialCaloricBoundaryGeometry {
        self.reference.geometry
    }

    /// Actual physical coordinates and explicitly retained angular footprint.
    #[must_use]
    pub const fn geometry(&self) -> RadialCaloricBoundaryGeometry {
        self.prepared.geometry
    }

    /// Exact current area and resistance-equivalent distance for the observer.
    /// These measures never replace the reference cell volume or mass.
    #[must_use]
    pub fn current_geometry(&self) -> CurrentConductiveBoundary {
        let geometry = self.prepared.native_geometry();
        CurrentConductiveBoundary {
            area_m2: geometry.area_m2,
            half_distance_m: geometry.half_distance_m,
        }
    }

    #[must_use]
    pub const fn native_geometry(&self) -> ImplicitBoundaryGeometry {
        self.prepared.native_geometry()
    }

    pub fn response_at_specific_energy(
        &self,
        specific_energy_j_kg: f64,
        bath_temperature_k: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ReservoirLayerResponse, DepositionError> {
        self.prepared.response_at_specific_energy(
            specific_energy_j_kg,
            bath_temperature_k,
            checkpoint,
        )
    }
}

impl ImplicitBoundaryLaw for CurrentRadialCaloricBoundary<'_, '_, '_, '_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.prepared.caloric_law()
    }

    fn supports_current_geometry(&self) -> bool {
        true
    }

    fn state_at_specific_enthalpy(
        &self,
        specific_energy_j_kg: f64,
        bath_temperature_k: f64,
        geometry: ImplicitBoundaryGeometry,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitBoundaryState, DepositionError> {
        self.prepared.state_at_specific_enthalpy(
            specific_energy_j_kg,
            bath_temperature_k,
            geometry,
            checkpoint,
        )
    }
}
