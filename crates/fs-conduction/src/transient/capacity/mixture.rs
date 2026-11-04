//! Fixed constituent inventories in the original P1 temperature-storage owner.

use super::{MaterialHeatCapacity, VolumetricHeatCapacity, poll};
use crate::{ConductionError, ConductionMesh, material::TemperatureSpan};
use fs_exec::Cx;

/// One nonreacting constituent of a stationary, locally equilibrated composite.
///
/// Fractions refer to additive volumes of the supplied reference mesh, not
/// mass fractions. Source density and Cp remain separate per constituent:
/// `rho_mix = sum(phi_i rho_i)` and `C_vol = sum(phi_i rho_i Cp_i)`.
/// Multiplying separately averaged density and Cp is a different, incorrect
/// energy account. Fractions and original sources are sampled only at setup.
///
/// This declares one shared temperature and zero heat of mixing. It supplies
/// thermal storage only; no effective conductivity, mechanics, phase change,
/// constituent motion or experimentally validated mixture law is inferred.
#[derive(Debug, Clone, Copy)]
pub struct MaterialHeatCapacityVolumeFraction<'a> {
    /// Original constant-density/Cp source, including its temperature support.
    pub source: &'a MaterialHeatCapacity,
    /// Dimensionless P1 volume fraction at each mesh vertex, in `[0, 1]`.
    /// The constituent sum at every vertex must differ from one by at most
    /// `64 * f64::EPSILON`. Values are never clipped or renormalized.
    pub nodal_volume_fraction: &'a [f64],
}

pub(in crate::transient) struct PreparedVolumeFractions {
    pub(in crate::transient) capacities: Vec<VolumetricHeatCapacity>,
    pub(in crate::transient) spans: Vec<TemperatureSpan>,
}

pub(in crate::transient) fn prepare(
    cx: &Cx<'_>,
    mesh: &ConductionMesh,
    constituents: &[MaterialHeatCapacityVolumeFraction<'_>],
) -> Result<PreparedVolumeFractions, ConductionError> {
    poll(cx, 0)?;
    if constituents.is_empty() {
        return Err(invalid("at least one constituent is required"));
    }
    let n = mesh.vertex_count();
    for (index, constituent) in constituents.iter().enumerate() {
        if index % 512 == 0 {
            poll(cx, index)?;
        }
        if constituent.nodal_volume_fraction.len() != n {
            return Err(ConductionError::FieldLength {
                field: "constituent volume fraction",
                expected: n,
                found: constituent.nodal_volume_fraction.len(),
            });
        }
    }
    let mut capacities = Vec::with_capacity(n);
    for vertex in 0..n {
        if vertex % 512 == 0 {
            poll(cx, vertex)?;
        }
        let mut fraction_sum = 0.0;
        let mut capacity = 0.0;
        for (index, constituent) in constituents.iter().enumerate() {
            if index % 512 == 0 {
                poll(cx, vertex)?;
            }
            let fraction = constituent.nodal_volume_fraction[vertex];
            if !fraction.is_finite() || !(0.0..=1.0).contains(&fraction) {
                return Err(invalid("volume fractions must be finite and lie in [0, 1]"));
            }
            fraction_sum += fraction;
            capacity += fraction * constituent.source.value_j_per_m3_k();
        }
        if (fraction_sum - 1.0).abs() > 64.0 * f64::EPSILON {
            return Err(invalid(
                "constituent volume fractions must sum to one at every vertex",
            ));
        }
        capacities.push(VolumetricHeatCapacity::declared(capacity)?);
    }

    let mut spans = vec![TemperatureSpan::Unbounded; n];
    for (element, tet) in mesh.complex().tets.iter().enumerate() {
        if element % 512 == 0 {
            poll(cx, element)?;
        }
        for (index, constituent) in constituents.iter().enumerate() {
            if index % 512 == 0 {
                poll(cx, element)?;
            }
            if !tet
                .iter()
                .any(|&v| constituent.nodal_volume_fraction[v as usize] > 0.0)
            {
                continue;
            }
            // A nonnegative P1 constituent present at one corner is present
            // throughout the tet interior. Checking only its nonzero corners
            // could evaluate it beyond its source domain near another corner.
            for &vertex in tet {
                let span = spans[vertex as usize].intersect(constituent.source.temperature_span());
                if matches!(span, TemperatureSpan::Sampled { low, high } if low > high) {
                    return Err(invalid(
                        "active constituents have disjoint temperature support",
                    ));
                }
                spans[vertex as usize] = span;
            }
        }
    }
    poll(cx, n)?;
    Ok(PreparedVolumeFractions { capacities, spans })
}

fn invalid(what: &str) -> ConductionError {
    ConductionError::Config {
        parameter: "material volume-fraction heat capacity",
        what: what.to_owned(),
    }
}
