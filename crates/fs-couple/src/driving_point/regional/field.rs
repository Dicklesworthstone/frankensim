//! Midpoint discretization of the shared primitive gas-field owner.
//!
//! Only constitutive sampling is new: each cell still uses `fs-duct`'s cylinder
//! or spherical cone propagation and the same pressure/volume-flow interfaces.
//! Subdividing geometry does not introduce another loss or acoustic equation.
//! A uniform field retains the original geometry and operation order exactly.

use fs_material::gas_field::AxialGasField;

use super::{Duct, DuctError, MAX_REGIONAL_SEGMENTS, RegionalDuct, Segment, invalid};

/// Shared material/propagation refusal from the continuous duct owner.
pub use fs_duct::profile::gas_field::GasFieldDuctError;

impl RegionalDuct {
    /// Bind continuous temperature/composition and smooth radius geometry to
    /// the existing fourth-order Magnus and stable two-port owners.
    ///
    /// For a varying field, steps_per_span is the integration resolution per
    /// authored smooth interval, with the fs-duct total-work cap. No intermediate
    /// gas coefficients or radii are averaged. Uniform fields retain the exact
    /// original regional/uniform operation order. The older from_gas_field
    /// constructor remains the explicitly piecewise-constant midpoint rung.
    ///
    /// This freezes the field during rendering. Nonuniform tone holes and wall
    /// pins need separate resolved physics and refuse instead of being dropped.
    pub fn from_continuous_gas_field(
        duct: &Duct,
        field: &AxialGasField,
        steps_per_span: fs_duct::RiccatiSteps,
    ) -> Result<Self, GasFieldDuctError> {
        if duct.segments.is_empty() {
            return Err(DuctError::EmptyDuct.into());
        }
        if duct.segments.len() > MAX_REGIONAL_SEGMENTS {
            return Err(invalid("gas-field geometry exceeds MAX_REGIONAL_SEGMENTS").into());
        }
        let continuous = fs_duct::profile::gas_field::GasFieldDuct::new(
            duct, field, steps_per_span,
        )?;
        if field.uniform_state().is_some() {
            return Self::new(duct, &vec![*continuous.inlet_gas(); duct.segments.len()])
                .map_err(Into::into);
        }
        let mut offset = 0.0;
        let mut gases = Vec::with_capacity(duct.segments.len());
        for segment in &duct.segments {
            gases.push(*field.state_at(offset)?.state());
            offset += segment.axial_length();
        }
        let mut regional = Self::new(duct, &gases)?;
        regional.one_way_time_s = continuous.one_way_time_s();
        regional.inlet_gas_override = Some(*continuous.inlet_gas());
        regional.continuous = Some(continuous);
        Ok(regional)
    }

    /// The actual continuous acoustic binding, if selected at construction.
    pub fn continuous_gas_field(&self) -> Option<&fs_duct::profile::gas_field::GasFieldDuct> {
        self.continuous.as_ref()
    }

    /// Apply a supplied axial T/composition field to a bore without changing
    /// its physics code. Each authored cylinder/cone gets exactly
    /// `cells_per_segment` midpoint cells, with at most 64 total records.
    /// Increasing this axis converges the prescribed field representation;
    /// it does not solve the flow or heat transfer that created the field.
    ///
    /// Inlet source density and characteristic impedance use the exact x=0
    /// state. Terminal radiation uses the exact x=L state through the existing
    /// exterior-medium binding. The exterior can subsequently be replaced
    /// explicitly by `with_exterior_gas`; no air default is inferred.
    ///
    /// Uniform fields keep every original segment and the uniform solver's
    /// arithmetic, including existing tone holes and smooth radius profiles.
    /// Nonuniform midpoint fields admit cylinders/cones only. Their tone-hole
    /// mutual coupling and simultaneous smooth-radius integration require the
    /// continuous-field duct owner and are refused here.
    ///
    /// # Errors
    /// Missing/invalid/oversized geometry or subdivision, a length mismatch,
    /// unsupported nonuniform geometry, or an existing material/duct refusal.
    pub fn from_gas_field(
        duct: &Duct,
        field: &AxialGasField,
        cells_per_segment: usize,
    ) -> Result<Self, GasFieldDuctError> {
        if cells_per_segment == 0 || cells_per_segment > MAX_REGIONAL_SEGMENTS {
            return Err(invalid("gas-field subdivision needs 1..=64 cells per segment").into());
        }
        if duct.segments.is_empty() {
            return Err(DuctError::EmptyDuct.into());
        }
        if duct.segments.len() > MAX_REGIONAL_SEGMENTS {
            return Err(invalid("gas-field geometry exceeds MAX_REGIONAL_SEGMENTS").into());
        }
        let length: f64 = duct.segments.iter().map(Segment::axial_length).sum();
        if !length.is_finite() || length <= 0.0 || length != field.length_m() {
            return Err(invalid("gas field must cover the exact complete duct length").into());
        }
        if let Some(uniform) = field.uniform_state() {
            return Self::new(duct, &vec![*uniform.state(); duct.segments.len()])
                .map_err(Into::into);
        }
        let count = duct
            .segments
            .len()
            .checked_mul(cells_per_segment)
            .filter(|count| *count <= MAX_REGIONAL_SEGMENTS)
            .ok_or_else(|| invalid("gas-field subdivision exceeds MAX_REGIONAL_SEGMENTS"))?;
        if duct
            .segments
            .iter()
            .any(|segment| matches!(segment, Segment::Profile { .. } | Segment::ToneHole { .. }))
        {
            return Err(invalid(
                "nonuniform gas fields currently need cylinder/cone geometry without tone holes",
            )
            .into());
        }
        let inlet = *field.state_at(0.0)?.state();
        let exterior = *field.state_at(field.length_m())?.state();
        let mut cells = Vec::with_capacity(count);
        let mut positions = Vec::with_capacity(count);
        let mut offset = 0.0;
        for segment in &duct.segments {
            let length = segment.axial_length();
            if !length.is_finite() || length <= 0.0 {
                return Err(invalid("gas-field axial lengths must be positive and finite").into());
            }
            for i in 0..cells_per_segment {
                let fraction0 = i as f64 / cells_per_segment as f64;
                let fraction1 = (i + 1) as f64 / cells_per_segment as f64;
                let x0 = offset + length * fraction0;
                let x1 = offset + length * fraction1;
                let cell_length = x1 - x0;
                if !cell_length.is_finite() || cell_length <= 0.0 {
                    return Err(invalid("gas-field cell length is numerically unresolved").into());
                }
                positions.push(x0 + 0.5 * cell_length);
                cells.push(match *segment {
                    Segment::Cylinder { radius, .. } => Segment::Cylinder {
                        radius,
                        length: cell_length,
                    },
                    Segment::Cone {
                        inlet_radius,
                        outlet_radius,
                        ..
                    } => {
                        let radius = |fraction| {
                            if fraction == 0.0 {
                                inlet_radius
                            } else if fraction == 1.0 {
                                outlet_radius
                            } else {
                                inlet_radius + fraction * (outlet_radius - inlet_radius)
                            }
                        };
                        Segment::Cone {
                            inlet_radius: radius(fraction0),
                            outlet_radius: radius(fraction1),
                            length: cell_length,
                        }
                    }
                    Segment::Profile { .. } | Segment::ToneHole { .. } => {
                        unreachable!("admitted above")
                    }
                });
            }
            offset += length;
        }
        let gases = field.sample_states(&positions)?;
        let mut regional =
            Self::new(&Duct { segments: cells }, &gases)?.with_exterior_gas(exterior)?;
        regional.inlet_gas_override = Some(inlet);
        Ok(regional)
    }
}

#[cfg(test)]
pub(crate) mod tests;
