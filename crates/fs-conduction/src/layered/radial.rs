//! Generated heat in radially symmetric multilayer hollow cylinders and spheres.
//! Geometry supplies exact resistance/source moments to the existing nonlinear
//! generated-wall owner. Material integration, domain checks, cancellation and
//! the bounded shooting solve are shared with the planar implementation.

use super::geometry::{LayerShape, RadialKind};
use super::{
    LayeredWall, NonlinearGeneratedLayeredWallSolution, NonlinearGenerationOptions, checkpoint,
    finite, invalid, positive,
};
use crate::ConductionError;
use fs_exec::Cx;

mod solid;
pub use solid::{
    GeneratedSolidRadialSolution, SolidRadialGeometry, SolidRadialLayerTemperatures,
    SolidRadialLayeredBody,
};

/// Fixed radial symmetry. Cylinders have insulated ends; spheres have full
/// solid angle. The positive inner radius is a physical reservoir boundary.
/// A solid centre, axial/lateral spreading and evolving geometry are not implied.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RadialWallGeometry {
    Cylindrical { inner_radius_m: f64, length_m: f64 },
    Spherical { inner_radius_m: f64 },
}

/// Original source layers bound to increasing radii. Each layer thickness is
/// radial, each contact resistance is per ACTUAL interface area, and both films
/// are per their respective physical boundary area. No common-area assumption
/// is imposed on those physical resistances.
#[derive(Debug, Clone, PartialEq)]
pub struct RadialLayeredWall {
    wall: LayeredWall,
    geometry: RadialWallGeometry,
    radii_m: Vec<f64>,
    shapes: Vec<LayerShape>,
    reference_area_m2: f64,
    contacts_reference_m2_k_w: Vec<f64>,
    surfaces_reference_m2_k_w: [f64; 2],
}

impl RadialLayeredWall {
    /// Bind the existing source stack to a cylindrical shell of positive length.
    /// Source/process conditions retained by each table must describe this
    /// specimen and its radial conductivity; geometry creates no new claims.
    pub fn cylindrical(
        wall: LayeredWall,
        inner_radius_m: f64,
        length_m: f64,
    ) -> Result<Self, ConductionError> {
        positive("cylinder length", length_m)?;
        Self::new(
            wall,
            RadialWallGeometry::Cylindrical {
                inner_radius_m,
                length_m,
            },
        )
    }

    /// Bind the existing source stack to a complete hollow spherical shell.
    pub fn spherical(wall: LayeredWall, inner_radius_m: f64) -> Result<Self, ConductionError> {
        Self::new(wall, RadialWallGeometry::Spherical { inner_radius_m })
    }

    fn new(wall: LayeredWall, geometry: RadialWallGeometry) -> Result<Self, ConductionError> {
        let (kind, inner, area) = match geometry {
            RadialWallGeometry::Cylindrical {
                inner_radius_m,
                length_m,
            } => (
                RadialKind::Cylindrical,
                inner_radius_m,
                (2.0 * core::f64::consts::PI * inner_radius_m) * length_m,
            ),
            RadialWallGeometry::Spherical { inner_radius_m } => (
                RadialKind::Spherical,
                inner_radius_m,
                (4.0 * core::f64::consts::PI * inner_radius_m) * inner_radius_m,
            ),
        };
        positive("inner shell radius", inner)?;
        positive("inner shell area", area)?;
        let mut shapes = Vec::with_capacity(wall.layers.len());
        let mut radii_m = Vec::with_capacity(wall.layers.len() + 1);
        let mut contacts = Vec::with_capacity(wall.contact_resistances_m2_k_per_w.len());
        radii_m.push(inner);
        let mut radius = inner;
        for (index, layer) in wall.layers.iter().enumerate() {
            let shape = admit_shell_shape(kind, radius, inner, layer.thickness_m, area)?;
            radius = finite("outer layer radius", radius + layer.thickness_m)?;
            let ratio = shape.area_ratio(layer.thickness_m);
            if let Some(&resistance) = wall.contact_resistances_m2_k_per_w.get(index) {
                contacts.push(reference_resistance(resistance, ratio)?);
            }
            shapes.push(shape);
            radii_m.push(radius);
        }
        let outer_ratio = shapes.last().expect("admitted nonempty stack").area_ratio(
            wall.layers
                .last()
                .expect("admitted nonempty stack")
                .thickness_m,
        );
        let surfaces = [
            wall.surface_resistances_m2_k_per_w[0],
            reference_resistance(wall.surface_resistances_m2_k_per_w[1], outer_ratio)?,
        ];
        Ok(Self {
            wall,
            geometry,
            radii_m,
            shapes,
            reference_area_m2: area,
            contacts_reference_m2_k_w: contacts,
            surfaces_reference_m2_k_w: surfaces,
        })
    }

    /// Unchanged material tables, source receipts and physical areal resistances.
    pub fn wall(&self) -> &LayeredWall {
        &self.wall
    }
    pub const fn geometry(&self) -> RadialWallGeometry {
        self.geometry
    }
    /// Increasing physical layer boundaries, including the inner and outer face.
    pub fn radii_m(&self) -> &[f64] {
        &self.radii_m
    }
    /// Inner boundary area used ONLY to normalize the shared scalar solver.
    pub const fn reference_area_m2(&self) -> f64 {
        self.reference_area_m2
    }

    /// Solve -div(k(T) grad T)=s with uniform signed generation in each shell.
    /// Q(r) changes by the exact generated volume, while local flux is Q/A(r).
    /// The original piecewise-linear k(T) is integrated and inverted analytically;
    /// supported boundary temperatures cannot hide an unsupported interior extremum.
    ///
    /// The solver's absolute flux tolerance is referenced to the INNER area:
    /// an allowance e W/m2 means A_inner*e W in every reported heat-rate residual.
    /// All physical contacts/films use their OWN area. Temperature and
    /// source-temperature domain and bounded-work controls otherwise retain their meanings.
    /// There is no transient storage or temperature-dependent generation feedback.
    pub fn solve_temperature_dependent_with_generation(
        &self,
        cx: &Cx<'_>,
        reservoir_temperatures_k: [f64; 2],
        volumetric_generation_w_m3: &[f64],
        options: NonlinearGenerationOptions,
    ) -> Result<NonlinearGeneratedRadialWallSolution, ConductionError> {
        checkpoint(cx, 0)?;
        let normalized = self.wall.solve_generated_geometry(
            cx,
            reservoir_temperatures_k,
            self.reference_area_m2,
            volumetric_generation_w_m3,
            options,
            self.shapes.clone(),
            &self.contacts_reference_m2_k_w,
            self.surfaces_reference_m2_k_w,
        )?;
        let mut layers = Vec::with_capacity(self.shapes.len());
        for (index, (row, &shape)) in normalized
            .layer_temperatures
            .iter()
            .zip(&self.shapes)
            .enumerate()
        {
            checkpoint(cx, index)?;
            let radii = [self.radii_m[index], self.radii_m[index + 1]];
            let ratios = [shape.area_ratio(0.0), shape.area_ratio(shape.thickness())];
            let rates = [
                scaled(
                    row.face_fluxes_w_per_m2[0],
                    self.reference_area_m2,
                    "inner face heat rate",
                )?,
                scaled(
                    row.face_fluxes_w_per_m2[1],
                    self.reference_area_m2,
                    "outer face heat rate",
                )?,
            ];
            layers.push(NonlinearGeneratedRadialLayerTemperatures {
                face_radii_m: radii,
                face_areas_m2: ratios.map(|ratio| ratio * self.reference_area_m2),
                face_temperatures_k: row.face_temperatures_k,
                radial_heat_rates_w: rates,
                radial_heat_fluxes_w_per_m2: [
                    scaled(
                        row.face_fluxes_w_per_m2[0],
                        1.0 / ratios[0],
                        "inner face heat flux",
                    )?,
                    scaled(
                        row.face_fluxes_w_per_m2[1],
                        1.0 / ratios[1],
                        "outer face heat flux",
                    )?,
                ],
                generated_heat_rate_w: scaled(
                    normalized.volumetric_generation_w_m3[index] * shape.volume(),
                    self.reference_area_m2,
                    "layer generated heat",
                )?,
                minimum_temperature_k: row.minimum_temperature_k,
                minimum_radius_m: finite(
                    "minimum temperature radius",
                    radii[0] + row.minimum_position_m,
                )?,
                maximum_temperature_k: row.maximum_temperature_k,
                maximum_radius_m: finite(
                    "maximum temperature radius",
                    radii[0] + row.maximum_position_m,
                )?,
                recomputed_resistance_weighted_heat_rate_w: finite(
                    "recomputed heat rate",
                    row.recomputed_mean_flux_w_per_m2 * self.reference_area_m2,
                )?,
                heat_rate_residual_w: finite(
                    "layer heat rate residual",
                    row.mean_flux_residual_w_per_m2 * self.reference_area_m2,
                )?,
            });
        }
        Ok(NonlinearGeneratedRadialWallSolution {
            geometry: self.geometry,
            reference_area_m2: self.reference_area_m2,
            reservoir_temperatures_k,
            reservoir_heat_rates_w: normalized.reservoir_heat_rates_w,
            generated_heat_rate_w: normalized.generated_heat_rate_w,
            energy_residual_w: normalized.energy_residual_w,
            layer_temperatures: layers,
            maximum_heat_rate_residual_w: finite(
                "maximum heat rate residual",
                normalized.maximum_flux_residual_w_per_m2 * self.reference_area_m2,
            )?,
            right_boundary_residual_k: normalized.right_boundary_residual_k,
            iterations: normalized.iterations,
            inner_heat_rate_bracket_w: [
                finite(
                    "lower heat rate bracket",
                    normalized.left_flux_bracket_w_per_m2[0] * self.reference_area_m2,
                )?,
                finite(
                    "upper heat rate bracket",
                    normalized.left_flux_bracket_w_per_m2[1] * self.reference_area_m2,
                )?,
            ],
            options,
            normalized,
        })
    }
}

/// Actual local fluxes, heat rates, temperatures and radii of one radial layer.
#[derive(Debug, Clone, PartialEq)]
pub struct NonlinearGeneratedRadialLayerTemperatures {
    pub face_radii_m: [f64; 2],
    pub face_areas_m2: [f64; 2],
    pub face_temperatures_k: [f64; 2],
    /// Signed toward increasing radius; inner and outer rates differ by generation.
    pub radial_heat_rates_w: [f64; 2],
    /// Actual Q/A on each face, not the private reference-area flux.
    pub radial_heat_fluxes_w_per_m2: [f64; 2],
    pub generated_heat_rate_w: f64,
    pub minimum_temperature_k: f64,
    pub minimum_radius_m: f64,
    pub maximum_temperature_k: f64,
    pub maximum_radius_m: f64,
    /// Aref*integral(k dT)/D, using signed integral from the outer to inner face.
    /// With generation this is a resistance-weighted heat rate, not a face rate.
    pub recomputed_resistance_weighted_heat_rate_w: f64,
    pub heat_rate_residual_w: f64,
}

/// Completed steady radial solve retaining every original material/source table.
/// Numerical residuals and the final bracket are diagnostics, not interval or
/// material uncertainty certificates.
#[derive(Debug, Clone, PartialEq)]
pub struct NonlinearGeneratedRadialWallSolution {
    pub geometry: RadialWallGeometry,
    pub reference_area_m2: f64,
    pub reservoir_temperatures_k: [f64; 2],
    /// Signed heat delivered outward to the inner and outer reservoirs.
    pub reservoir_heat_rates_w: [f64; 2],
    pub generated_heat_rate_w: f64,
    pub energy_residual_w: f64,
    pub layer_temperatures: Vec<NonlinearGeneratedRadialLayerTemperatures>,
    pub maximum_heat_rate_residual_w: f64,
    pub right_boundary_residual_k: f64,
    pub iterations: usize,
    /// Bracket on signed heat toward increasing radius at the INNER boundary.
    pub inner_heat_rate_bracket_w: [f64; 2],
    pub options: NonlinearGenerationOptions,
    normalized: NonlinearGeneratedLayeredWallSolution,
}

impl NonlinearGeneratedRadialWallSolution {
    pub fn wall(&self) -> &LayeredWall {
        &self.normalized.wall
    }
    /// Original signed uniform source in each physical layer, W/m3.
    pub fn volumetric_generation_w_m3(&self) -> &[f64] {
        &self.normalized.volumetric_generation_w_m3
    }

    /// Recover the source-integrated temperature at a physical radial offset
    /// inside one layer. The offset is measured from that layer's inner face.
    pub fn temperature_at(
        &self,
        cx: &Cx<'_>,
        layer: usize,
        offset_m: f64,
    ) -> Result<f64, ConductionError> {
        self.normalized.temperature_at(cx, layer, offset_m)
    }

    /// Signed local heat flux at a physical radial offset, W/m2.
    pub fn heat_flux_at(
        &self,
        cx: &Cx<'_>,
        layer: usize,
        offset_m: f64,
    ) -> Result<f64, ConductionError> {
        checkpoint(cx, layer)?;
        let shape = self
            .normalized
            .shapes
            .get(layer)
            .ok_or_else(|| invalid("heat flux query layer index is out of range"))?;
        if !offset_m.is_finite() || offset_m < 0.0 || offset_m > shape.thickness() {
            return Err(invalid(
                "heat flux query must lie inside its physical layer",
            ));
        }
        let row = &self.normalized.layer_temperatures[layer];
        let reference = if offset_m == 0.0 {
            row.face_fluxes_w_per_m2[0]
        } else if offset_m == shape.thickness() {
            row.face_fluxes_w_per_m2[1]
        } else {
            finite(
                "local reference heat flux",
                row.face_fluxes_w_per_m2[0]
                    + self.normalized.volumetric_generation_w_m3[layer]
                        * shape.volume_length(offset_m),
            )?
        };
        scaled(
            reference,
            1.0 / shape.area_ratio(offset_m),
            "local radial heat flux",
        )
    }
}

fn reference_resistance(value: f64, area_ratio: f64) -> Result<f64, ConductionError> {
    let normalized = finite("reference-area resistance", value / area_ratio)?;
    if value != 0.0 && normalized == 0.0 {
        return Err(invalid("nonzero radial interface resistance underflowed"));
    }
    Ok(normalized)
}

fn scaled(value: f64, scale: f64, what: &'static str) -> Result<f64, ConductionError> {
    let result = finite(what, value * scale)?;
    if value != 0.0 && result == 0.0 {
        return Err(invalid("nonzero radial heat rate or flux underflowed"));
    }
    Ok(result)
}

fn admit_shell_shape(
    kind: RadialKind,
    inner: f64,
    reference: f64,
    thickness: f64,
    reference_area: f64,
) -> Result<LayerShape, ConductionError> {
    let shape = LayerShape::radial(kind, inner, reference, thickness)
        .ok_or_else(|| invalid("radial layer metrics are nonpositive or unrepresentable"))?;
    positive(
        "outer layer area",
        shape.area_ratio(thickness) * reference_area,
    )?;
    positive("physical layer volume", shape.volume() * reference_area)?;
    Ok(shape)
}
