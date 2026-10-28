//! Source-bound axial bearings on the full fixed-end span, including afterlength.
//! The existing solid-mass axial modes and original shell eigenvectors supply
//! one signed displacement/force column. Its reaction is part of the original
//! complete linear/friction/strain solve, never a delayed body force.

use super::*;
use super::super::super::interfaces::{
    MaterialBridgeInterface, outer_material, validate_frequencies,
    validate_surface_states,
};
use crate::render::schedule::force::coupled::{
    ConstantLossConnection, ModalColumnConnection,
};
use crate::render::schedule::force::coupled::contact::friction::FrictionModalSystem;
use fs_material::state_point::ResolvedLinearInterface;
use super::super::super::shell::ShellPointGeometry;

/// Cold connection capacity, independent of instrument and physical coefficients.
pub const MAX_AXIAL_BEARINGS: usize = 16;

/// Actual contact geometry in the body mesh's global XYZ frame.
/// The complete string still has fixed endpoints; this bearing is INTERIOR.
#[derive(Clone, Debug, PartialEq)]
pub struct AxialBearingGeometry {
    pub string: usize,
    /// Fraction of the complete supported span, not the speaking length.
    pub station_fraction: f64,
    /// Original shell-assembly part; zero for a single material shell.
    pub shell_component: usize,
    pub node: usize,
    /// Positive string axis, a finite unit vector without implicit normalization.
    pub string_axis: [f64; 3],
    /// Vector [m] from the physical shell node to the bearing point.
    pub contact_arm_m: [f64; 3],
    /// Stress-free axial relative displacement [m].
    pub rest_extension_m: f64,
}
impl AxialBearingGeometry {
    /// Pure geometry admission before material I/O or modal allocation.
    pub fn validate(&self) -> Result<(), ImpactError> {
        let norm2 = self.string_axis.iter().map(|x| x*x).sum::<f64>();
        if !self.station_fraction.is_finite() || !(0.0 < self.station_fraction && self.station_fraction < 1.0)
            || !norm2.is_finite() || (norm2 - 1.0).abs() > 1e-12
            || self.contact_arm_m.iter().any(|x| !x.is_finite())
            || !self.rest_extension_m.is_finite() {
            return Err(invalid("axial bearing requires an interior full-span station, unit axis and finite arm/rest geometry"));
        }
        Ok(())
    }

    /// Physical force/moment per unit axial force: [t, r cross t]. Its dot
    /// product with [U, Theta] is t dot (U + Theta cross r), including its sign.
    #[cfg(test)]
    fn nodal_row(&self) -> Result<[f64; 6], ImpactError> {
        self.validate()?;
        self.shell_geometry().nodal_row().map_err(|e| ImpactError::Owner(e.to_string()))
    }

    fn shell_geometry(&self) -> ShellPointGeometry {
        ShellPointGeometry { shell_component: self.shell_component, node: self.node,
            direction: self.string_axis, arm_m: self.contact_arm_m }
    }
}

/// Original resolved surface law and actual geometry, retained through playback.
/// Its coefficients are obtained from the existing linear interface owner.
/// Bulk string material alone cannot supply an interface stiffness or loss.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialAxialBearing {
    pub geometry: AxialBearingGeometry,
    pub interface: ResolvedLinearInterface,
}

fn receiver<'a>(
    source: &'a MaterialStringBankSource, geometry: &AxialBearingGeometry, gate: &CancelGate,
) -> Result<(Vec<f64>, &'a ResolvedMaterialStatePoint), ImpactError> {
    geometry.validate()?;
    let point = source.shell_point(&geometry.shell_geometry(), gate)
        .map_err(|e| ImpactError::Owner(e.to_string()))?;
    Ok((point.body.shapes, point.material))
}

pub(super) fn lower(
    source: &MaterialStringBankSource, strings: &[MaterialLongitudinalString],
    components: &[ModalAcousticTimeModel], bearings: &[MaterialAxialBearing],
    gate: &CancelGate,
) -> Result<(Vec<ModalColumnConnection>, Vec<ConstantLossConnection>), ImpactError> {
    if bearings.is_empty() { return Ok((Vec::new(), Vec::new())); }
    if bearings.len() > MAX_AXIAL_BEARINGS {
        return Err(invalid("axial bearing count exceeds the declared owner capacity"));
    }
    let old_count: usize = components.iter().map(|model| model.modes().len()).sum();
    let count = strings.iter().try_fold(old_count, |n, law| {
        n.checked_add(law.model.axial_modes().len()).ok_or_else(|| invalid("axial bearing coordinate overflow"))
    })?;
    if count > source.options.network.max_modes
        || count.checked_mul(bearings.len()).is_none_or(|n| n > source.options.network.max_setup_terms) {
        return Err(invalid("axial bearing columns exceed original complete-network budgets"));
    }
    let body = components.get(source.strings.len()).ok_or_else(|| invalid("axial bearing has no original body"))?;
    let body_offset: usize = components[..source.strings.len()].iter().map(|m| m.modes().len()).sum();
    let mut columns = Vec::with_capacity(bearings.len());
    let mut anchors = Vec::new();
    for (index, bearing) in bearings.iter().enumerate() {
        if gate.is_requested() { return Err(ImpactError::Cancelled); }
        let geometry = &bearing.geometry;
        geometry.validate()?;
        let law = strings.get(geometry.string).ok_or_else(|| invalid("axial bearing names an absent material string"))?;
        let member = &source.strings[geometry.string];
        let outer = outer_material(&member.source).map_err(|e| ImpactError::Owner(e.to_string()))?;
        let (body_shapes, receiving) = receiver(source, geometry, gate)?;
        if body_shapes.len() != body.modes().len() {
            return Err(invalid("axial bearing body row differs from the original physical basis"));
        }
        validate_surface_states(&bearing.interface, outer, receiving)
            .map_err(|e| ImpactError::Owner(e.to_string()))?;
        let axial_shapes = law.model.axial_displacement_shapes(geometry.station_fraction)
            .map_err(|e| ImpactError::Owner(e.to_string()))?;
        let axial_offset = old_count + strings[..geometry.string].iter()
            .map(|s| s.model.axial_modes().len()).sum::<usize>();
        let mut column = vec![0.0; count];
        column[axial_offset..axial_offset+axial_shapes.len()].copy_from_slice(&axial_shapes);
        for (slot, shape) in column[body_offset..body_offset+body_shapes.len()].iter_mut().zip(body_shapes) {
            *slot = -shape;
        }
        columns.push(ModalColumnConnection {
            column, stiffness_n_m: bearing.interface.stiffness_n_m(),
            damping_n_s_m: bearing.interface.damping_n_s_m(),
            rest_extension_m: geometry.rest_extension_m,
        });
        if let Some(reference_omega_rad_s) = bearing.interface.constant_loss_reference_omega_rad_s() {
            anchors.push(ConstantLossConnection { connection: index, reference_omega_rad_s });
        }
    }
    Ok((columns, anchors))
}

pub(super) fn validate_bands(
    source: &MaterialStringBankSource, strings: &[MaterialLongitudinalString],
    bearings: &[MaterialAxialBearing], system: &FrictionModalSystem, gate: &CancelGate,
) -> Result<(), ImpactError> {
    let bindings: Vec<_> = bearings.iter().map(|bearing| MaterialBridgeInterface {
        string: bearing.geometry.string, interface: bearing.interface.clone(),
    }).collect();
    // The original physical prefix excludes appended gas memory. Axial
    // coordinates are physical too and are screened explicitly below.
    let physical = source.strings.len() + 1 + usize::from(source.material_cavity.is_some())
        + usize::from(source.spatial_cavity.as_ref().is_some_and(|cavity|
            cavity.coupling().total_modes() > cavity.coupling().structural_modes()));
    validate_frequencies(physical, &bindings, system, gate)
        .map_err(|e| ImpactError::Owner(e.to_string()))?;
    if !bearings.is_empty() {
        for interface in bearings.iter().map(|bearing| &bearing.interface)
            .chain(source.material_bridge_interfaces().iter().map(|binding| &binding.interface)) {
            for law in strings {
                if gate.is_requested() { return Err(ImpactError::Cancelled); }
                for &omega in law.model.axial_omegas() {
                    interface.require_angular_frequency_rad_s(omega)
                        .map_err(|e| ImpactError::Owner(e.to_string()))?;
                }
            }
        }
        // Passive spring columns add positive semidefinite stiffness. The
        // uncoupled physical minimum and complete stiffness upper bound enclose
        // the retained undamped spectrum conservatively; neither is a
        // qualification of nonlinear harmonics or omitted continuum modes.
        let lower = system.components()[..physical].iter()
            .flat_map(|component| component.modes().iter().map(|mode| mode.angular_frequency_rad_s))
            .chain(strings.iter().flat_map(|law| law.model.axial_omegas().iter().copied()))
            .fold(f64::INFINITY, f64::min);
        let upper = system.linear_frequency_upper_bound_rad_s();
        for law in strings {
            let loss = law.axial_material_loss.as_ref().map(|loss| loss.omega_band_rad_s());
            for (lo, hi) in law.elastic_frequency_band_rad_s.into_iter().chain(loss) {
                if lower < lo || upper > hi {
                    return Err(invalid("axial bearing coupled spectrum exceeds original axial material applicability"));
                }
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn physical_axis_and_arm_supply_translation_rotation_and_reciprocal_work() {
        let geometry = AxialBearingGeometry { string: 0, station_fraction: 0.8,
            shell_component: 0, node: 0, string_axis: [1.0, 0.0, 0.0],
            contact_arm_m: [0.0, 0.0, 0.02], rest_extension_m: 0.0 };
        let row = geometry.nodal_row().unwrap();
        assert_eq!(row, [1.0, 0.0, 0.0, 0.0, 0.02, 0.0]);
        let dq = [0.003, 0.007, -0.002, 0.01, 0.2, -0.03];
        let displacement = row.iter().zip(dq).map(|(b, q)| b*q).sum::<f64>();
        assert!((displacement - (0.003 + 0.02*0.2)).abs() < 1e-16);
        let force = 7.0;
        let modal_work = row.iter().zip(dq).map(|(b, q)| force*b*q).sum::<f64>();
        assert!((modal_work-force*displacement).abs() < 1e-16);
        let mut zero_arm = geometry.clone(); zero_arm.contact_arm_m = [0.0; 3];
        assert_eq!(&zero_arm.nodal_row().unwrap()[3..], &[0.0; 3]);
        for fraction in [0.0, 1.0, f64::NAN] {
            let mut invalid = geometry.clone(); invalid.station_fraction = fraction;
            assert!(invalid.validate().is_err());
        }
    }

    #[test]
    fn original_longitudinal_owner_supplies_solid_mass_shapes_and_exact_fixed_ends() {
        let build = |density| LongitudinalString::new(KcStringParams {
            length: 0.5, lin_density: density, tension: 20.0, ea: 1000.0,
        }, 2, 1, &[1, 2]).unwrap();
        let a = build(0.002); let b = build(0.008);
        let row = a.axial_displacement_shapes(0.25).unwrap();
        let scale = det::sqrt(2.0/(0.002*0.5));
        assert!((row[0] - scale*det::sqrt(0.5)).abs() < 1e-13);
        assert!((row[1] - scale).abs() < 1e-13);
        for (original, heavier) in row.iter().zip(b.axial_displacement_shapes(0.25).unwrap()) {
            assert!((heavier*2.0-original).abs() < 1e-13);
        }
        assert_eq!(a.axial_displacement_shapes(0.0).unwrap(), vec![0.0; 2]);
        assert_eq!(a.axial_displacement_shapes(1.0).unwrap(), vec![0.0; 2]);
        assert!(a.axial_displacement_shapes(1.01).is_err());
    }
}
