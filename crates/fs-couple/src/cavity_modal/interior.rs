//! Source-gas binding for the three-dimensional thermoviscous FE interior.
//!
//! Gas calorics and transport stay in `fs-material`; boundary-layer geometry,
//! original pressure traces and the frequency solve stay in `fs-feec`. This
//! adapter only binds their numerical views to the same coherent gas state.

use fs_feec::acoustics::helmholtz::boundary_layer::thermal_wall::{
    ConductiveThermalWall, ThermalWall, ThermalWallDepth, WallThermalTransport,
};
use fs_feec::acoustics::helmholtz::boundary_layer::{
    BoundaryLayerGeometry, BoundaryLayerTransport,
};
use fs_feec::acoustics::helmholtz::{
    AcousticFeError, HelmholtzFe, HelmholtzModes, ReducedHelmholtz,
};
use fs_material::gas::GasState;
use fs_material::thermal::transport::ThermalTransportSample;

/// Read wall transport from one coherent source gas; no air constants or
/// independent heat-capacity/Prandtl approximations are introduced.
pub fn boundary_layer_transport(gas: &GasState) -> Result<BoundaryLayerTransport, AcousticFeError> {
    gas.validate_coherence()
        .map_err(|_| AcousticFeError::InvalidInput {
            what: "thermoviscous interior requires a coherent source gas",
        })?;
    Ok(BoundaryLayerTransport {
        dynamic_viscosity_pa_s: gas.dynamic_viscosity,
        thermal_conductivity_w_m_k: gas.thermal_conductivity,
        cp_j_kg_k: gas.specific_heat_cp,
        heat_capacity_ratio: gas.gamma,
    })
}

/// Attach thermoviscous walls to the original FE problem using the same gas
/// that supplied its density and sound speed. A different volume medium refuses
/// before reduction, even if each individual state would otherwise be valid.
/// The caller retains the FE owner's geometric admission and approximation
/// limits; this adapter introduces no additional physical model or solver.
pub fn thermoviscous_interior(
    fe: &HelmholtzFe<'_>,
    modes: &HelmholtzModes,
    ports: &[&[[u32; 3]]],
    static_vectors: bool,
    wall_faces: &[[u32; 3]],
    gas: &GasState,
    geometry: BoundaryLayerGeometry,
) -> Result<ReducedHelmholtz, AcousticFeError> {
    let transport = boundary_layer_transport(gas)?;
    let medium = fe.medium();
    if medium.density_kg_m3.to_bits() != gas.density.to_bits()
        || medium.sound_speed_m_s.to_bits() != gas.sound_speed.to_bits()
    {
        return Err(AcousticFeError::InvalidInput {
            what: "thermoviscous wall gas differs from the original FE volume medium",
        });
    }
    ReducedHelmholtz::new_with_boundary_layers(
        fe,
        modes,
        ports,
        static_vectors,
        wall_faces,
        transport,
        geometry,
    )
}

/// Bind a source-evaluated solid to the same gas equilibrium temperature.
/// The supplied sample retains its original prepared-material identity; the
/// caller owns the source path, caloric constraint and wall-normal direction.
/// No heat capacity is converted and no material constants are inferred.
pub fn conductive_wall_transport(
    gas: &GasState,
    wall: &ThermalTransportSample,
) -> Result<WallThermalTransport, AcousticFeError> {
    boundary_layer_transport(gas)?;
    if gas.temperature.to_bits() != wall.temperature_k().to_bits() {
        return Err(AcousticFeError::InvalidInput {
            what: "conductive acoustic wall and gas require one equilibrium temperature",
        });
    }
    Ok(WallThermalTransport {
        conductivity_w_m_k: wall.conductivity_w_m_k(),
        volumetric_heat_capacity_j_m3_k: wall.volumetric_heat_capacity_j_m3_k(),
    })
}

/// Original three-dimensional pressure FE with conjugate gas/solid thermal
/// boundary layers and the original viscous tangential surface form. Material
/// preparation occurs outside this builder; every response only evaluates the
/// exact normal diffusion profile of a homogeneous solid at fixed equilibrium.
/// Geometry, backing condition and local-layer admission remain explicit.
pub fn thermally_conductive_interior(
    fe: &HelmholtzFe<'_>,
    modes: &HelmholtzModes,
    ports: &[&[[u32; 3]]],
    static_vectors: bool,
    wall_faces: &[[u32; 3]],
    gas: &GasState,
    geometry: BoundaryLayerGeometry,
    wall: &ThermalTransportSample,
    depth: ThermalWallDepth,
    minimum_tangential_scale_m: f64,
) -> Result<ReducedHelmholtz, AcousticFeError> {
    let transport = conductive_wall_transport(gas, wall)?;
    thermoviscous_interior(fe, modes, ports, static_vectors, wall_faces, gas, geometry)?
        .with_thermal_wall(ThermalWall::Conductive(ConductiveThermalWall {
            transport,
            depth,
            minimum_tangential_scale_m,
        }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_feec::acoustics::helmholtz::Medium;
    use fs_feec::fixtures::single_tet;
    use fs_material::gas::GasSpec;

    mod conductive;

    const GEOMETRY: BoundaryLayerGeometry = BoundaryLayerGeometry {
        minimum_curvature_radius_m: f64::INFINITY,
        minimum_half_gap_m: 0.05,
        maximum_layer_ratio: 0.05,
    };

    fn basis(fe: &HelmholtzFe<'_>) -> HelmholtzModes {
        let n = fe.ndof();
        HelmholtzModes {
            lambda: vec![0.0; n],
            frequency_hz: vec![0.0; n],
            residual: vec![0.0; n],
            shapes: (0..n)
                .map(|i| {
                    let mut p = vec![0.0; n];
                    p[i] = 1.0;
                    p
                })
                .collect(),
            below_low: 0,
            below_high: n,
        }
    }

    #[test]
    fn original_dry_and_humid_source_gases_reach_full_fe_wall_work() {
        let (complex, points) = single_tet();
        let points: Vec<_> = points.iter().map(|p| p.map(|x| 0.2 * x)).collect();
        let gases = [
            GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap(),
            GasState::try_new_moist_air(310.15, 101325.0, 0.9).unwrap(),
        ];
        let mut observed = Vec::new();
        for gas in gases {
            let transport = boundary_layer_transport(&gas).unwrap();
            for (actual, expected) in [
                (transport.dynamic_viscosity_pa_s, gas.dynamic_viscosity),
                (
                    transport.thermal_conductivity_w_m_k,
                    gas.thermal_conductivity,
                ),
                (transport.cp_j_kg_k, gas.specific_heat_cp),
                (transport.heat_capacity_ratio, gas.gamma),
            ] {
                assert_eq!(actual.to_bits(), expected.to_bits());
            }
            let medium = Medium {
                density_kg_m3: gas.density,
                sound_speed_m_s: gas.sound_speed,
            };
            let fe = HelmholtzFe::new(&complex, &points, 1, medium).unwrap();
            let faces = fe.boundary_faces();
            let source = &faces[..1];
            let model = thermoviscous_interior(
                &fe,
                &basis(&fe),
                &[source],
                false,
                &faces[1..],
                &gas,
                GEOMETRY,
            )
            .unwrap();
            let response = model.response(300.0, 0, 0.01, &[[0.0, 0.0]]).unwrap();
            let power = model
                .boundary_layer_power(300.0, &response.amplitudes)
                .unwrap();
            let input = 0.5 * response.inflow * response.mean_pressure[0][0];
            assert!(power.viscous_w > 0.0 && power.thermal_w > 0.0);
            assert!((input - power.viscous_w - power.thermal_w).abs() < 1e-9 * input.abs());
            observed.push(response.mean_pressure[0]);
        }
        assert_ne!(observed[0], observed[1]);
    }

    #[test]
    fn mixed_or_incoherent_volume_and_wall_gases_refuse() {
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
        let (complex, positions) = single_tet();
        let medium = Medium {
            density_kg_m3: gas.density,
            sound_speed_m_s: f64::from_bits(gas.sound_speed.to_bits() + 1),
        };
        let fe = HelmholtzFe::new(&complex, &positions, 1, medium).unwrap();
        let faces = fe.boundary_faces();
        assert!(
            thermoviscous_interior(
                &fe,
                &basis(&fe),
                &[&faces[..1]],
                false,
                &faces[1..],
                &gas,
                GEOMETRY
            )
            .is_err()
        );
        let corrupt = GasState {
            specific_heat_cp: gas.specific_heat_cp * 1.01,
            ..gas
        };
        assert!(boundary_layer_transport(&corrupt).is_err());
    }
}
