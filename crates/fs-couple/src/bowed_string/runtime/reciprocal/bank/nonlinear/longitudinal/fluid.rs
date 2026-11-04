//! Local axial Stokes shear in the original longitudinal string/body network.
//!
//! The cylinder owner supplies Z_parallel=2*pi*mu*z*K1(z)/K0(z). Its positive
//! viscous spectrum is half the transverse spectrum, with no instantaneous
//! added mass. Projecting a uniform local operator onto fixed-end sine modes
//! gives one reciprocal fluid branch per retained axial mode. Each branch has
//! its own real inertia and history; its reaction participates in the existing
//! exact linear/friction/finite-strain solve and work ledger.
//!
//! This shares the selected transverse gas's frozen source, radius and spectral
//! resolution. It neglects axial fluid gradients and finite-cylinder ends, as
//! required by the local slender-cylinder approximation. No axial drag is
//! inferred for banks which did not select the causal cylinder operator.
use super::*;
use crate::render::schedule::force::coupled::ModalColumnConnection;

/// Bound all new state and column allocations before material/modal assembly.
/// The native network owner subsequently admits the complete old+new link work.
pub(super) fn admit(
    source: &MaterialStringBankSource,
    resolutions: &[LongitudinalStringResolution],
    original_modes: usize,
    bearings: usize,
    gate: &CancelGate,
) -> Result<(), ImpactError> {
    let loading = source.cylinder_gas_loading();
    if loading.is_some_and(|fluid| fluid.strings().len() != resolutions.len()) {
        return Err(invalid(
            "axial gas loading lost its original material strings",
        ));
    }
    let mut count = original_modes;
    let mut links = bearings;
    for (index, resolution) in resolutions.iter().enumerate() {
        if gate.is_requested() {
            return Err(ImpactError::Cancelled);
        }
        resolution.validate()?;
        let modes = resolution.axial_modes.len();
        let memory = loading
            .map_or(Some(0), |fluid| {
                modes.checked_mul(fluid.strings()[index].flow().branches().len())
            })
            .ok_or_else(|| invalid("axial gas-memory dimension overflow"))?;
        count = count
            .checked_add(modes)
            .and_then(|n| n.checked_add(memory))
            .ok_or_else(|| invalid("complete longitudinal dimension overflow"))?;
        links = links
            .checked_add(memory)
            .ok_or_else(|| invalid("axial gas connection count overflow"))?;
    }
    let budget = source.options.network;
    let extent = links
        .checked_add(1)
        .ok_or_else(|| invalid("axial gas connection setup overflow"))?;
    let work = count
        .checked_mul(extent)
        .and_then(|n| n.checked_mul(extent));
    if count > budget.max_modes
        || links > budget.max_connections
        || work.is_none_or(|terms| terms > budget.max_setup_terms)
    {
        return Err(invalid(
            "longitudinal modes and ALL axial gas memory exceed original network budgets",
        ));
    }
    Ok(())
}

pub(super) fn tail_damping_ratio(
    source: &MaterialStringBankSource,
    string: usize,
    solid_line_mass: f64,
    omega: f64,
) -> Result<f64, ImpactError> {
    let Some(loading) = source.cylinder_gas_loading() else {
        return Ok(0.0);
    };
    let flow = loading
        .strings()
        .get(string)
        .ok_or_else(|| invalid("axial gas tail lost its original string"))?
        .flow();
    let ratio = 0.5 * flow.axial_low_tail_resistance_per_length_n_s_m2() / solid_line_mass / omega;
    if !ratio.is_finite() || ratio <= 0.0 {
        return Err(invalid(
            "axial gas tail is unrepresentable in the original solid basis",
        ));
    }
    Ok(ratio)
}

pub(super) fn append(
    source: &MaterialStringBankSource,
    strings: &[MaterialLongitudinalString],
    original: &[ModalAcousticTimeModel],
    additional: &mut Vec<ModalAcousticTimeModel>,
    connections: &mut Vec<ModalColumnConnection>,
    gate: &CancelGate,
) -> Result<(), ImpactError> {
    let Some(loading) = source.cylinder_gas_loading() else {
        return Ok(());
    };
    if loading.strings().len() != strings.len() || additional.len() != strings.len() {
        return Err(invalid(
            "axial gas must follow every original solid axial component",
        ));
    }
    let old_count: usize = original.iter().map(|c| c.modes().len()).sum();
    let axial_count: usize = additional.iter().map(|c| c.modes().len()).sum();
    let memory_count: usize = strings
        .iter()
        .zip(loading.strings())
        .map(|(law, fluid)| law.model.axial_modes().len() * fluid.flow().branches().len())
        .sum();
    // All sums/products and this complete column size were admitted before
    // construction, against the original owner (at most 4096 coordinates).
    let count = old_count + axial_count + memory_count;
    for connection in connections.iter_mut() {
        if gate.is_requested() {
            return Err(ImpactError::Cancelled);
        }
        connection.column.resize(count, 0.0);
    }
    let mut axial_offset = old_count;
    let mut fluid_offset = old_count + axial_count;
    for (index, (law, loading)) in strings.iter().zip(loading.strings()).enumerate() {
        if gate.is_requested() {
            return Err(ImpactError::Cancelled);
        }
        let component = original.len() + additional.len();
        if law.axial_fluid_component != Some(component)
            || law.axial_component != original.len() + index
        {
            return Err(invalid(
                "axial gas changed an original solid component address",
            ));
        }
        let flow = loading.flow();
        let length = loading.solid_card().length_m;
        let mode_count = law.model.axial_modes().len();
        let branch_count = flow.branches().len();
        let mut modes = Vec::with_capacity(mode_count * branch_count);
        let solid_shape = det::sqrt(2.0 / law.solid_linear_density_kg_m / length);
        for mode in 0..mode_count {
            for (branch_index, branch) in flow.axial_branches().enumerate() {
                if gate.is_requested() {
                    return Err(ImpactError::Cancelled);
                }
                let fluid_shape = det::sqrt(2.0 / branch.mass_per_length_kg_m / length);
                let resistance = branch.resistance_per_length_n_s_m2 * (0.5 * length);
                if [solid_shape, fluid_shape, resistance]
                    .iter()
                    .any(|value| !value.is_finite() || *value <= 0.0)
                {
                    return Err(invalid(
                        "axial fluid/string conjugate port is unrepresentable",
                    ));
                }
                let mut column = vec![0.0; count];
                column[axial_offset + mode] = solid_shape;
                column[fluid_offset + mode * branch_count + branch_index] = -fluid_shape;
                connections.push(ModalColumnConnection {
                    column,
                    stiffness_n_m: 0.0,
                    damping_n_s_m: resistance,
                    rest_extension_m: 0.0,
                });
                modes.push(ModalAcousticMode {
                    angular_frequency_rad_s: 0.0,
                    damping_ratio: 0.0,
                    pressure_per_modal_velocity: C64::ZERO,
                });
            }
        }
        additional.push(
            ModalAcousticTimeModel::try_new_with_free_coordinates(
                source.options.sample_rate_hz,
                modes,
                source.options.string_budget,
            )
            .map_err(|e| ImpactError::Owner(format!("axial gas modal admission: {e:?}")))?,
        );
        axial_offset += mode_count;
        fluid_offset += mode_count * branch_count;
    }
    Ok(())
}
