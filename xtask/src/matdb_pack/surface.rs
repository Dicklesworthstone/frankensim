//! Source-owned two-axis scalar laws for manifest v9 / material source v3.
//!
//! The tab-separated row has twelve fields:
//! `surface2d claim observations property axis0 unit0 axis1 unit1 value-unit
//! x:y:value,... i:j:k,... interpolation`. Triangle indices are zero-based
//! sample ordinals; `-` declares no triangles. Source order is retained, and
//! neither missing samples nor connectivity are inferred. Both source axes
//! need explicit manifest axis declarations, including `dimensional` when no
//! semantic kind is asserted. Existing validity, uncertainty, observation and
//! joint-statistics rows remain authoritative.

use super::*;

#[derive(Debug)]
pub(super) struct RawSurface2D {
    axes: [String; 2],
    axis_units: [String; 2],
    ordinate_unit: String,
    pub(super) points: Vec<[String; 3]>,
    triangles: Vec<[u32; 3]>,
}

pub(super) fn parse(fields: &[&str], subject: &str) -> Result<RawSurface2D, CompileError> {
    let malformed = |detail| CompileError::new("malformed_surface", subject, detail);
    if fields[9].is_empty() || fields[9].split(',').count() > fs_matdb::MAX_SURFACE_SAMPLES {
        return Err(malformed(
            "surface samples must be nonempty and fit the runtime sample limit",
        ));
    }
    let points = fields[9]
        .split(',')
        .map(|point| {
            let fields: Vec<_> = point.split(':').collect();
            let [x, y, value] = fields.as_slice() else {
                return Err(malformed("each surface sample must be x:y:value"));
            };
            Ok([
                require_number_token(x, "first surface coordinate", subject)?,
                require_number_token(y, "second surface coordinate", subject)?,
                require_number_token(value, "surface ordinate", subject)?,
            ])
        })
        .collect::<Result<Vec<_>, CompileError>>()?;
    let triangles = if fields[10] == "-" {
        Vec::new()
    } else {
        if fields[10].is_empty() || fields[10].split(',').count() > fs_matdb::MAX_SURFACE_TRIANGLES
        {
            return Err(malformed(
                "surface connectivity must fit the runtime triangle limit",
            ));
        }
        fields[10]
            .split(',')
            .map(|triangle| {
                let fields: Vec<_> = triangle.split(':').collect();
                let [a, b, c] = fields.as_slice() else {
                    return Err(malformed(
                        "each surface triangle must contain three sample indices",
                    ));
                };
                Ok([
                    parse_canonical_u32(a, "surface sample index", subject)?,
                    parse_canonical_u32(b, "surface sample index", subject)?,
                    parse_canonical_u32(c, "surface sample index", subject)?,
                ])
            })
            .collect::<Result<Vec<_>, CompileError>>()?
    };
    Ok(RawSurface2D {
        axes: [
            require_identifier(fields[4], "first surface axis", subject)?,
            require_identifier(fields[6], "second surface axis", subject)?,
        ],
        axis_units: [
            require_unit(fields[5], "first surface axis unit", subject)?,
            require_unit(fields[7], "second surface axis unit", subject)?,
        ],
        ordinate_unit: require_unit(fields[8], "surface ordinate unit", subject)?,
        points,
        triangles,
    })
}

pub(super) fn compile(
    source: &RawSurface2D,
    interpolation: InterpolationPolicy,
    kind: Option<QuantityKind>,
    manifest: &Manifest,
    subject: &str,
) -> Result<(PropertyValue, Vec<(LocalComponent, ParsedQuantity)>), CompileError> {
    if interpolation == InterpolationPolicy::ConstantWithinValidity {
        return Err(CompileError::new(
            "incompatible_interpolation",
            subject,
            "surface claims admit only linear or tabulated interpolation",
        ));
    }
    for axis in &source.axes {
        if !manifest.axis_kinds.contains_key(axis) {
            return Err(CompileError::new(
                "missing_surface_axis_declaration",
                subject,
                format!(
                    "surface axis {axis:?} needs an explicit manifest axis kind or dimensional declaration"
                ),
            ));
        }
    }
    let mut samples = Vec::with_capacity(source.points.len());
    let mut components = Vec::with_capacity(source.points.len() * 3);
    let mut dimensions = None;
    for (index, [x, y, value]) in source.points.iter().enumerate() {
        let sample = u32::try_from(index).map_err(|_| {
            CompileError::new(
                "resource_limit",
                subject,
                "surface sample index exceeds u32",
            )
        })?;
        let a = parse_declared_quantity(
            x,
            &source.axis_units[0],
            manifest.axis_kinds[&source.axes[0]],
            &format!("{subject}:surface:axis:0:{index}"),
        )?;
        let b = parse_declared_quantity(
            y,
            &source.axis_units[1],
            manifest.axis_kinds[&source.axes[1]],
            &format!("{subject}:surface:axis:1:{index}"),
        )?;
        let ordinate = parse_declared_quantity(
            value,
            &source.ordinate_unit,
            kind,
            &format!("{subject}:surface:value:{index}"),
        )?;
        let current_dims = [a.dims, b.dims, ordinate.dims];
        if dimensions.is_some_and(|expected| expected != current_dims) {
            return Err(CompileError::new(
                "surface_dims_mismatch",
                subject,
                "surface sample dimensions must remain identical",
            ));
        }
        dimensions = Some(current_dims);
        samples.push(([a.value, b.value], ordinate.value));
        components.extend([
            (LocalComponent::SurfaceCoordinate(sample, 0), a),
            (LocalComponent::SurfaceCoordinate(sample, 1), b),
            (LocalComponent::SurfaceOrdinate(sample), ordinate),
        ]);
    }
    let [x_dims, y_dims, dims] = dimensions.ok_or_else(|| {
        CompileError::new(
            "malformed_surface",
            subject,
            "surface requires at least one sample",
        )
    })?;
    Ok((
        PropertyValue::Surface2D {
            axes: source.axes.clone(),
            axis_dims: [x_dims, y_dims],
            samples,
            triangles: source.triangles.clone(),
            dims,
        },
        components,
    ))
}

/// Match original coordinate units to original validity before applying aliases.
/// The runtime owner then admits finite bounds, sample coverage and topology.
pub(super) fn bind_axes(
    axes: &mut [String; 2],
    dims: [Dims; 2],
    validities: Option<&BTreeMap<String, RawValidity>>,
    manifest: &Manifest,
    claim: &str,
) -> Result<(), CompileError> {
    let subject = format!("claim:{claim}");
    for (axis, dims) in axes.iter_mut().zip(dims) {
        let declared = validities
            .and_then(|values| values.get(axis))
            .ok_or_else(|| {
                CompileError::new(
                    "missing_surface_validity",
                    &subject,
                    "both surface axes require explicit source validity intervals",
                )
            })?;
        let lower = parse_declared_quantity(
            &declared.lower,
            &declared.unit,
            manifest.axis_kinds.get(axis).copied().flatten(),
            &format!("{subject}:validity:{axis}:lower"),
        )?;
        if lower.dims != dims {
            return Err(CompileError::new(
                "surface_validity_dims_mismatch",
                &subject,
                "surface coordinate dimensions disagree with its validity axis",
            ));
        }
        *axis = manifest.axis(axis).to_owned();
    }
    Ok(())
}
