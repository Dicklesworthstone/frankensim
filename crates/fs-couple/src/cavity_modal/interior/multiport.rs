//! Source-coherent spatial interiors in the original duct composition algebra.
//!
//! The FE owner uses exp(+i omega t); ducts use exp(-i omega t). With the same
//! real pressure/flow coordinates and positive real power references, the
//! conversion is entrywise conjugation. Port flow stays inward on BOTH sides;
//! no extra sign flip, area factor or modal renormalization is introduced.
//! Geometry, source gas, pressure/trace resolution and feature qualification
//! remain with their original owners. This is not a fitted feature or proof
//! that a two-coordinate truncation is adequate near a real junction.

use super::{
    AcousticFeError, BoundaryLayerGeometry, BoundaryLayerTransport, ConductiveThermalWall,
    GasState, HelmholtzFe, HelmholtzModes, ThermalTransportSample, ThermalWall, ThermalWallDepth,
    boundary_layer_transport, conductive_wall_transport,
};
use fs_duct::compose::TwoPort;
use fs_feec::acoustics::helmholtz::multiport::{BoundaryMode, InteriorMultiport};
use fs_math::c64::C64;

fn invalid(what: &'static str) -> AcousticFeError {
    AcousticFeError::InvalidInput { what }
}

fn source_medium(
    fe: &HelmholtzFe<'_>,
    gas: &GasState,
) -> Result<BoundaryLayerTransport, AcousticFeError> {
    let transport = boundary_layer_transport(gas)?;
    let medium = fe.medium();
    if medium.density_kg_m3.to_bits() != gas.density.to_bits()
        || medium.sound_speed_m_s.to_bits() != gas.sound_speed.to_bits()
    {
        return Err(invalid(
            "spatial interior multiport differs from its source gas",
        ));
    }
    Ok(transport)
}

/// Bind the same coherent source gas to an original lossless spatial FE
/// multiport. Each supplied real trace retains its exact declared normalization.
pub fn interior_multiport(
    fe: &HelmholtzFe<'_>,
    modes: &HelmholtzModes,
    ports: &[BoundaryMode],
    static_vectors: bool,
    gas: &GasState,
    cancelled: &dyn Fn() -> bool,
) -> Result<InteriorMultiport, AcousticFeError> {
    source_medium(fe, gas)?;
    InteriorMultiport::new(fe, modes, ports, static_vectors, cancelled)
}

/// Spatial pressure/flow modes with the original source-gas transport and
/// tangential viscous/thermal wall operator. Wall/port overlap and physical
/// thin-layer admission remain enforced by the FE owner.
pub fn thermoviscous_multiport(
    fe: &HelmholtzFe<'_>,
    modes: &HelmholtzModes,
    ports: &[BoundaryMode],
    static_vectors: bool,
    wall_faces: &[[u32; 3]],
    gas: &GasState,
    geometry: BoundaryLayerGeometry,
    cancelled: &dyn Fn() -> bool,
) -> Result<InteriorMultiport, AcousticFeError> {
    let transport = source_medium(fe, gas)?;
    InteriorMultiport::new_with_boundary_layers(
        fe,
        modes,
        ports,
        static_vectors,
        wall_faces,
        transport,
        geometry,
        cancelled,
    )
}

/// Source-owned homogeneous solid conduction coupled to the SAME spatial
/// multiport. Original gas/solid equilibrium-temperature coherence, directional
/// conductivity and selected volumetric heat capacity are preserved.
pub fn thermally_conductive_multiport(
    fe: &HelmholtzFe<'_>,
    modes: &HelmholtzModes,
    ports: &[BoundaryMode],
    static_vectors: bool,
    wall_faces: &[[u32; 3]],
    gas: &GasState,
    geometry: BoundaryLayerGeometry,
    wall: &ThermalTransportSample,
    depth: ThermalWallDepth,
    minimum_tangential_scale_m: f64,
    cancelled: &dyn Fn() -> bool,
) -> Result<InteriorMultiport, AcousticFeError> {
    let transport = conductive_wall_transport(gas, wall)?;
    thermoviscous_multiport(
        fe,
        modes,
        ports,
        static_vectors,
        wall_faces,
        gas,
        geometry,
        cancelled,
    )?
    .with_thermal_wall(ThermalWall::Conductive(ConductiveThermalWall {
        transport,
        depth,
        minimum_tangential_scale_m,
    }))
}

/// Row-major power scattering in the duct's exp(-i omega t) convention.
/// Positive real references [Pa s/m3] are passed unchanged to the original FE
/// solve. Coordinate order and area normalization are unchanged. For existing
/// `fs_duct::compose::redheffer_star`, the caller groups left/right modes and
/// ensures matched connected references; no evanescent traveling-wave or
/// Euclidean-passivity assumption is added by this conversion.
pub fn duct_scattering(
    interior: &InteriorMultiport,
    frequency_hz: f64,
    reference_resistance: &[f64],
    cancelled: &dyn Fn() -> bool,
) -> Result<Vec<C64>, AcousticFeError> {
    let mut scattering = interior.scattering(frequency_hz, reference_resistance, cancelled)?;
    for row in scattering.chunks_mut(interior.port_mode_count()) {
        if cancelled() {
            return Err(invalid("spatial interior duct conversion cancelled"));
        }
        for value in row {
            *value = value.conj();
        }
    }
    Ok(scattering)
}

/// Export exactly two explicitly retained coordinates to the original duct
/// `TwoPort`, ready for `compose_chain` or `TwoPort::star`. Both use the same
/// positive real reference. No extra modes are dropped to obtain this shape.
/// References are coordinate choices, not inferred characteristic impedances.
pub fn duct_two_port(
    interior: &InteriorMultiport,
    frequency_hz: f64,
    reference_impedance: f64,
    cancelled: &dyn Fn() -> bool,
) -> Result<TwoPort, AcousticFeError> {
    if interior.port_mode_count() != 2 {
        return Err(invalid(
            "duct TwoPort requires exactly two retained FE port modes",
        ));
    }
    let s = duct_scattering(interior, frequency_hz, &[reference_impedance; 2], cancelled)?;
    TwoPort::new([s[0], s[1], s[2], s[3]], reference_impedance)
        .map_err(|_| invalid("invalid spatial interior duct scattering block"))
}

#[cfg(test)]
mod tests;
