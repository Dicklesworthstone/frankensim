use super::*;
use fs_evidence::ValidityDomain;
use fs_feec::acoustics::helmholtz::boundary_layer::thermal_wall::ThermalBacking;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::thermal::transport::{PreparedThermalTransport, ThermalTransportSourceQuery};
use fs_qty::{Density, Dims, Temperature};

fn material(conductivity: f64, specific_heat: f64) -> PreparedThermalTransport {
    let roles = [
        ("density", Density::DIMS, 1000.0),
        ("specific_heat", Dims([2, 0, -2, -1, 0, 0]), specific_heat),
        (
            "normal_conductivity",
            Dims([1, 1, -3, -1, 0, 0]),
            conductivity,
        ),
    ];
    let mut claims = ClaimSet::new();
    for (name, dims, value) in roles {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(250.0, value), (350.0, value)],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    fs_qty::QuantitySpec::dimensional(Temperature::DIMS),
                    250.0,
                    350.0,
                ),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic conductive acoustic wall".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "manufactured verification fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    PreparedThermalTransport::try_new(
        &card,
        ThermalTransportSourceQuery {
            density: &PropertyKey::new(roles[0].0, roles[0].1),
            specific_heat: &PropertyKey::new(roles[1].0, roles[1].1),
            conductivity: &PropertyKey::new(roles[2].0, roles[2].1),
            temperature_axis: "T",
            lower: &QueryPoint::new()
                .with_quantity(
                    "T",
                    fs_qty::QuantitySpec::dimensional(Temperature::DIMS),
                    250.0,
                )
                .unwrap(),
            upper: &QueryPoint::new()
                .with_quantity(
                    "T",
                    fs_qty::QuantitySpec::dimensional(Temperature::DIMS),
                    350.0,
                )
                .unwrap(),
            selections: [ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly); 3],
        },
    )
    .unwrap()
}

#[test]
fn original_material_source_changes_conductive_wall_and_full_fe_work() {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let (complex, points) = single_tet();
    let points: Vec<_> = points.iter().map(|p| p.map(|x| 0.2 * x)).collect();
    let fe = HelmholtzFe::new(
        &complex,
        &points,
        1,
        Medium {
            density_kg_m3: gas.density,
            sound_speed_m_s: gas.sound_speed,
        },
    )
    .unwrap();
    let faces = fe.boundary_faces();
    let mut pressures = Vec::new();
    for (conductivity, specific_heat) in [(0.01, 1.0), (0.03, 1.0), (0.01, 3.0)] {
        let source = material(conductivity, specific_heat);
        let sample = source.at_temperature(gas.temperature).unwrap();
        assert_eq!(sample.path_identity(), source.identity());
        let transport = conductive_wall_transport(&gas, &sample).unwrap();
        assert_eq!(
            transport.conductivity_w_m_k.to_bits(),
            conductivity.to_bits()
        );
        assert_eq!(
            transport.volumetric_heat_capacity_j_m3_k.to_bits(),
            (1000.0 * specific_heat).to_bits()
        );
        let model = thermally_conductive_interior(
            &fe,
            &basis(&fe),
            &[&faces[..1]],
            false,
            &faces[1..],
            &gas,
            GEOMETRY,
            &sample,
            ThermalWallDepth::Slab {
                thickness_m: 1e-4,
                backing: ThermalBacking::Adiabatic,
            },
            0.1,
        )
        .unwrap();
        let response = model.response(300.0, 0, 0.01, &[[0.0, 0.0]]).unwrap();
        let loss = model
            .boundary_layer_power(300.0, &response.amplitudes)
            .unwrap();
        let input = 0.5 * response.inflow * response.mean_pressure[0][0];
        assert!(loss.thermal_w > 0.0 && loss.viscous_w > 0.0);
        assert!((input - loss.thermal_w - loss.viscous_w).abs() < 1e-9 * input.abs());
        pressures.push(response.mean_pressure[0]);
    }
    assert_ne!(pressures[0], pressures[1]);
    assert_ne!(pressures[0], pressures[2]);
}

#[test]
fn conductive_source_requires_one_gas_solid_equilibrium_temperature() {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let source = material(0.01, 1.0);
    assert!(conductive_wall_transport(&gas, &source.at_temperature(294.15).unwrap()).is_err());
    let sample = source.at_temperature(gas.temperature).unwrap();
    let corrupt = GasState {
        specific_heat_cp: gas.specific_heat_cp * 1.01,
        ..gas
    };
    assert!(conductive_wall_transport(&corrupt, &sample).is_err());
}
