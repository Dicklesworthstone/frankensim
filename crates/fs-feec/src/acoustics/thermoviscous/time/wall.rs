//! Actual normal heat fields in the original finite solid slabs.
//!
//! Each slab has a uniform P1 mesh through its physical thickness, tensorized
//! with its FULL original hierarchical P2 surface trace. The front coordinate
//! is the original fluid theta=T'/T0: temperature continuity is an identification,
//! not a penalty or a separate delayed boundary update. The consistent element
//! mass includes the fluid/interface-to-solid cross terms. At every depth the
//! original trace is retained, including coefficients whose front temperature
//! is constrained by an intersecting isothermal fluid face.
//!
//! C_vol*T0 integral(theta^2)/2 is stored in the SAME quadratic Hamiltonian.
//! k*T0 integral((d theta/d normal)^2) is a separate trajectory-integrated loss.
//! There is no second fluid boundary-layer loss, fitted thermal pole, lateral
//! solid conduction, multilayer homogenization or changing mean temperature.

use super::*;
use crate::acoustics::helmholtz::boundary_layer::thermal_wall::{ThermalBacking, ThermalWallDepth};
use std::collections::BTreeSet;

/// Explicit spatial resolution for the original finite wall regions.
/// The declared band includes DC. It is an applicability check, not a promise
/// that arbitrary initial fields or held drives are band limited or resolved.
#[derive(Clone, Copy, Debug)]
pub struct ThermoviscousWallDiscretization<'a> {
    /// Positive P1 element counts, in the original conductive-region order.
    /// No source region may be omitted. The complete state budget also applies.
    pub normal_elements: &'a [usize],
    /// Upper frequency of the declared normal-profile domain [Hz]. The clock
    /// must sample it; full slab thickness must admit the normal approximation
    /// down to DC. No penetration depth is substituted for physical thickness.
    pub maximum_frequency_hz: f64,
}

/// Ordering of the real stored temperature coefficients in one solid region.
#[derive(Clone, Debug)]
pub struct ConductiveWallTimeLayout {
    trace_dofs: Vec<usize>,
    depths_m: Vec<f64>,
    normal_elements: usize,
}
impl ConductiveWallTimeLayout {
    /// Original global hierarchical P2 indices, ascending and untruncated.
    pub fn trace_dofs(&self) -> &[usize] {
        &self.trace_dofs
    }
    /// Free normal-node depths [m], ascending from the fluid interface.
    /// The interface is shared with fluid state and is not duplicated here.
    /// An isothermal backing is constrained to zero and is not stored.
    pub fn depths_m(&self) -> &[f64] {
        &self.depths_m
    }
    pub const fn normal_elements(&self) -> usize {
        self.normal_elements
    }
    pub fn state_count(&self) -> usize {
        self.trace_dofs.len() * self.depths_m.len()
    }
}

pub(super) struct WallFields {
    pub layouts: Vec<ConductiveWallTimeLayout>,
    offsets: Vec<usize>,
    pub maximum_frequency_hz: Option<f64>,
}

const MAX_NORMAL_ELEMENTS: usize = 128;

impl WallFields {
    /// Admit and count before allocating any added dense matrix.
    pub fn new(
        model: &ThermoviscousFe<'_>,
        duration_s: f64,
        resolution: Option<ThermoviscousWallDiscretization<'_>>,
        mut mixed: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Self, usize), ThermoviscousTimeError> {
        let regions = model.conductive_walls();
        let mut fields = Self {
            layouts: Vec::with_capacity(regions.len()),
            offsets: Vec::with_capacity(regions.len()),
            maximum_frequency_hz: None,
        };
        let Some(resolution) = resolution else {
            if regions.is_empty() {
                return Ok((fields, mixed));
            }
            return Err(invalid(
                "finite conductive walls require explicit normal meshes and their frequency domain",
            ));
        };
        if resolution.normal_elements.len() != regions.len()
            || !resolution.maximum_frequency_hz.is_finite()
            || resolution.maximum_frequency_hz <= 0.0
            || duration_s * resolution.maximum_frequency_hz > 0.5
        {
            return Err(invalid(
                "complete wall meshes, positive frequency domain and a sampling clock required",
            ));
        }
        fields.maximum_frequency_hz = Some(resolution.maximum_frequency_hz);
        for (region, &elements) in regions.iter().zip(resolution.normal_elements) {
            poll(cancelled)?;
            let ThermalWallDepth::Slab {
                thickness_m,
                backing,
            } = region.wall.depth
            else {
                return Err(invalid(
                    "a semi-infinite wall cannot be replaced by an invented finite depth",
                ));
            };
            if elements == 0 || elements > MAX_NORMAL_ELEMENTS {
                return Err(invalid(
                    "bounded positive wall normal-element count required",
                ));
            }
            // Reuse original material/geometry validation and the worst
            // high-frequency acoustic-length check. The DC limit additionally
            // penetrates the WHOLE physical thickness, never merely one cell.
            region
                .wall
                .normal_admittance(
                    model.fe.medium(),
                    region.geometry,
                    resolution.maximum_frequency_hz,
                )
                .map_err(ThermoviscousError::from)?;
            let tangential_scale = region
                .geometry
                .minimum_curvature_radius_m
                .min(region.wall.minimum_tangential_scale_m);
            if thickness_m / tangential_scale > region.geometry.maximum_layer_ratio {
                return Err(invalid(
                    "complete finite wall violates normal-profile geometry at low frequency",
                ));
            }
            let mut trace = BTreeSet::new();
            for &face in &region.faces {
                poll(cancelled)?;
                let face = canonical(face);
                for vertex in face {
                    trace.insert(vertex as usize);
                }
                for (a, b) in [(0, 1), (0, 2), (1, 2)] {
                    let edge = model
                        .space
                        .complex
                        .edges
                        .binary_search(&[face[a], face[b]])
                        .expect("original admitted conductive boundary edge");
                    trace.insert(model.space.edge_off + edge);
                }
            }
            let free = elements - usize::from(backing == ThermalBacking::Isothermal);
            let trace_dofs: Vec<_> = trace.into_iter().collect();
            let count = free
                .checked_mul(trace_dofs.len())
                .ok_or_else(|| invalid("complete conductive wall field dimension overflow"))?;
            fields.offsets.push(mixed);
            mixed = mixed
                .checked_add(count)
                .ok_or_else(|| invalid("complete conductive wall state dimension overflow"))?;
            fields.layouts.push(ConductiveWallTimeLayout {
                trace_dofs,
                depths_m: (1..=free)
                    .map(|i| thickness_m * (i as f64 / elements as f64))
                    .collect(),
                normal_elements: elements,
            });
        }
        Ok((fields, mixed))
    }

    fn coordinate(
        &self,
        model: &ThermoviscousFe<'_>,
        region: usize,
        depth: usize,
        surface: usize,
    ) -> Option<usize> {
        let layout = &self.layouts[region];
        if depth == 0 {
            model.temperature_map[layout.trace_dofs[surface]].map(|i| model.temperature_start + i)
        } else if depth <= layout.depths_m.len() {
            Some(self.offsets[region] + (depth - 1) * layout.trace_dofs.len() + surface)
        } else {
            None // Original zero-temperature backside, eliminated exactly.
        }
    }

    pub fn assemble(
        &self,
        model: &ThermoviscousFe<'_>,
        mass: &mut [f64],
        mixed: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<f64>, ThermoviscousTimeError> {
        let mut conduction = vec![0.0; mixed * mixed];
        for (region_index, (region, layout)) in model
            .conductive_walls()
            .iter()
            .zip(&self.layouts)
            .enumerate()
        {
            poll(cancelled)?;
            let ThermalWallDepth::Slab { thickness_m, .. } = region.wall.depth else {
                unreachable!("admitted finite slab")
            };
            let (capacity, conductance) = element_coefficients(
                thickness_m,
                layout.normal_elements,
                region.wall.transport.volumetric_heat_capacity_j_m3_k,
                region.wall.transport.conductivity_w_m_k,
                model.medium.temperature_k,
            )?;
            // This is the SAME original P2 boundary quadrature as the harmonic
            // wall, before removing any isothermal front coordinates.
            let surface_mass = model
                .fe
                .boundary_mass(&region.faces)
                .map_err(ThermoviscousError::from)?;
            for depth in 0..layout.normal_elements {
                poll(cancelled)?;
                for (i, &source_i) in layout.trace_dofs.iter().enumerate() {
                    let (columns, values) = surface_mass.row(source_i);
                    for (&source_j, &surface) in columns.iter().zip(values) {
                        let Ok(j) = layout.trace_dofs.binary_search(&source_j) else {
                            continue;
                        };
                        for a in 0..2 {
                            let Some(row) = self.coordinate(model, region_index, depth + a, i)
                            else {
                                continue;
                            };
                            for b in 0..2 {
                                let Some(col) = self.coordinate(model, region_index, depth + b, j)
                                else {
                                    continue;
                                };
                                mass[row * mixed + col] +=
                                    surface * capacity * if a == b { 2.0 } else { 1.0 };
                                conduction[row * mixed + col] +=
                                    surface * conductance * if a == b { 1.0 } else { -1.0 };
                            }
                        }
                    }
                }
            }
        }
        if mass.iter().chain(&conduction).any(|v| !v.is_finite()) {
            return Err(invalid("nonfinite original solid diffusion assembly"));
        }
        Ok(conduction)
    }

    pub fn encode(
        &self,
        model: &ThermoviscousFe<'_>,
        temperatures: &[Vec<f64>],
        state: &mut [f64],
    ) -> Result<(), ThermoviscousTimeError> {
        if temperatures.len() != self.layouts.len() {
            return Err(invalid(
                "complete original conductive wall field set required",
            ));
        }
        for ((field, layout), &offset) in temperatures.iter().zip(&self.layouts).zip(&self.offsets)
        {
            if field.len() != layout.state_count() {
                return Err(invalid(
                    "complete original conductive wall temperature field required",
                ));
            }
            for (i, temperature) in field.iter().enumerate() {
                state[offset + i] = temperature / model.medium.temperature_k;
            }
        }
        Ok(())
    }

    pub fn decode(&self, model: &ThermoviscousFe<'_>, mixed: &[f64]) -> Vec<Vec<f64>> {
        self.layouts
            .iter()
            .zip(&self.offsets)
            .map(|(layout, &offset)| {
                mixed[offset..offset + layout.state_count()]
                    .iter()
                    .map(|theta| theta * model.medium.temperature_k)
                    .collect()
            })
            .collect()
    }
}

// The exact two-node P1 integrals, before tensorization with the original trace:
// M_e = capacity * [[2,1],[1,2]], K_e = conductance * [[1,-1],[-1,1]].
// Kept as a small pure kernel for direct checks against continuum heat balances.
pub(super) fn element_coefficients(
    thickness_m: f64,
    elements: usize,
    volumetric_capacity: f64,
    conductivity: f64,
    temperature_k: f64,
) -> Result<(f64, f64), ThermoviscousTimeError> {
    let dz = thickness_m / elements as f64;
    let capacity = volumetric_capacity * temperature_k * dz / 6.0;
    let conductance = conductivity * temperature_k / dz;
    if [dz, capacity, conductance]
        .into_iter()
        .any(|x| !x.is_finite() || x <= 0.0)
    {
        return Err(invalid(
            "finite positive solid normal-element storage and conductance required",
        ));
    }
    Ok((capacity, conductance))
}
