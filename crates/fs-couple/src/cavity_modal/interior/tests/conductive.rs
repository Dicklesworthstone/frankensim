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

#[test]
fn molecular_pressure_heat_field_retains_original_solid_source_and_work() {
    use crate::cavity_modal::interior::molecular::{MolecularInterior, MolecularThermalWall};
    use fs_feec::acoustics::helmholtz::multiport::BoundaryMode;
    use fs_material::relaxation::{
        InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw,
        RelaxingMedium, SpeciesRelaxation,
    };
    use fs_math::c64::C64;
    use fs_thermochem::SpeciesId;
    let spec = GasSpec {
        gamma: 1.25,
        ..GasSpec::dry_air_ussa1976()
    };
    let gas = GasState::try_new(&spec, 293.15, 101325.0).unwrap();
    let species = SpeciesId::new("verification-gas").unwrap();
    let state = RelaxationAbsorption::new(vec![SpeciesRelaxation {
        species: species.clone(),
        processes: vec![RelaxationProcess {
            reservoir: InternalReservoir::Rotational {
                degrees_of_freedom: 2,
            },
            rate: RelaxationRateLaw::PartnerIndependent {
                p_tau_ref_pa_s: 0.001 * gas.pressure,
                t_ref_k: gas.temperature,
                exponent: 0.0,
                t_min_k: 200.0,
                t_max_k: 400.0,
            },
        }],
    }])
    .unwrap()
    .evaluate(&RelaxingMedium::new(gas, vec![(species, 1.0)]).unwrap())
    .unwrap();
    let (complex, points) = single_tet();
    let points: Vec<_> = points.iter().map(|p| p.map(|v| 0.2 * v)).collect();
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
    let ports = [BoundaryMode::uniform(&faces[..1])];
    let mut pressures = Vec::new();
    for (conductivity, capacity) in [(0.01, 1.0), (0.03, 1.0), (0.01, 3.0)] {
        let source = material(conductivity, capacity);
        let sample = source.at_temperature(gas.temperature).unwrap();
        let model = MolecularInterior::new(
            &fe,
            &basis(&fe),
            &ports,
            &faces[1..],
            state.clone(),
            MolecularThermalWall::Conductive {
                sample: &sample,
                depth: ThermalWallDepth::Slab {
                    thickness_m: 1e-4,
                    backing: ThermalBacking::Adiabatic,
                },
                minimum_tangential_scale_m: 0.1,
            },
            GEOMETRY,
            &|| false,
        )
        .unwrap();
        assert_eq!(
            model.wall_source().unwrap().path_identity(),
            source.identity()
        );
        let response = model
            .prepare(300.0, &|| false)
            .unwrap()
            .response(&[C64::from_re(1e-6)], &|| false)
            .unwrap();
        let input = response.acoustic.input_power_w();
        assert!(
            response.dissipation.molecular_w > 0.0 && response.dissipation.solid_conduction_w > 0.0
        );
        assert!((input - response.dissipation.total_w()).abs() < 2e-9 * input);
        pressures.push(response.acoustic.pressure_pa[0]);
        let wrong = source.at_temperature(294.15).unwrap();
        assert!(
            MolecularInterior::new(
                &fe,
                &basis(&fe),
                &ports,
                &faces[1..],
                state.clone(),
                MolecularThermalWall::Conductive {
                    sample: &wrong,
                    depth: ThermalWallDepth::HalfSpace,
                    minimum_tangential_scale_m: 0.1
                },
                GEOMETRY,
                &|| false
            )
            .is_err()
        );
    }
    assert_ne!(pressures[0], pressures[1]);
    assert_ne!(pressures[0], pressures[2]);
}

fn resolved_wall_fixture<'c>(
    mesh: &'c fs_rep_mesh::TetComplex,
    points: &[[f64; 3]],
    gas: &GasState,
) -> (
    fs_feec::acoustics::thermoviscous::ThermoviscousFe<'c>,
    Vec<[u32; 3]>,
) {
    use fs_feec::acoustics::helmholtz::multiport::BoundaryMode;
    use fs_feec::acoustics::thermoviscous::{
        ThermoviscousBoundary, ThermoviscousFe, ThermoviscousLimits,
    };
    let fe = HelmholtzFe::new(
        mesh,
        points,
        2,
        Medium {
            density_kg_m3: gas.density,
            sound_speed_m_s: gas.sound_speed,
        },
    )
    .unwrap();
    let faces = fe.boundary_faces();
    let ports: Vec<_> = faces
        .iter()
        .copied()
        .filter(|face| face.iter().all(|&i| points[i as usize][0] == 0.0))
        .collect();
    let walls: Vec<_> = faces.into_iter().filter(|f| !ports.contains(f)).collect();
    let model = ThermoviscousFe::new(
        mesh,
        points,
        super::super::thermoviscous::source_medium(gas, 0.0).unwrap(),
        &ThermoviscousBoundary {
            no_slip_faces: walls.clone(),
            isothermal_faces: vec![],
            ports: vec![BoundaryMode::uniform(&ports)],
        },
        ThermoviscousLimits::default(),
        &mut || false,
    )
    .unwrap();
    (model, walls)
}

#[test]
fn resolved_finite_walls_keep_original_material_and_joint_molecular_work() {
    use super::super::thermoviscous::{ThermoviscousInterior, wall::MaterialWallRegion};
    use fs_feec::fixtures::kuhn_cube;
    use fs_material::relaxation::{RelaxingMedium, input::read_relaxation};
    use fs_math::c64::C64;

    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let law = read_relaxation(
        &mut "declared 1 synthetic 1 vibration 1 900 1 power 101.325 293.15 0 200 500"
            .split_ascii_whitespace(),
    )
    .unwrap();
    let source = law
        .evaluate(&RelaxingMedium::new(gas, vec![(law.species()[0].species.clone(), 1.0)]).unwrap())
        .unwrap();
    let (mesh, unit) = kuhn_cube(2);
    let points: Vec<_> = unit.iter().map(|p| p.map(|x| 0.02 * x)).collect();
    let mut flows = Vec::new();
    for (conductivity, capacity) in [(0.01, 1.0), (0.03, 1.0), (0.01, 3.0)] {
        let prepared = material(conductivity, capacity);
        let sample = prepared.at_temperature(gas.temperature).unwrap();
        let (fe, walls) = resolved_wall_fixture(&mesh, &points, &gas);
        let model = ThermoviscousInterior::from_relaxation(fe, source.clone(), [100.0, 1000.0])
            .unwrap()
            .with_conductive_walls(
                &[MaterialWallRegion {
                    faces: &walls,
                    sample: &sample,
                    depth: ThermalWallDepth::Slab {
                        thickness_m: 1e-4,
                        backing: ThermalBacking::Adiabatic,
                    },
                    minimum_tangential_scale_m: 0.02,
                    geometry: GEOMETRY,
                }],
                &mut || false,
            )
            .unwrap();
        assert_eq!(model.wall_sources()[0].path_identity(), prepared.identity());
        assert_eq!(model.inner().conductive_walls().len(), 1);
        let response = model
            .response(300.0, &[C64::new(1.0, 0.2)], &mut || false)
            .unwrap();
        let power = response.power;
        println!("resolved source solid k={conductivity} cp={capacity} power={power:?}");
        assert!(power.viscous_w > 0.0 && power.conductive_w > 0.0);
        assert!(power.caloric_w > 0.0 && power.solid_conductive_w > 0.0);
        assert!(power.balance_residual_w().abs() < 2e-8 * power.input_w);
        flows.push(response.inward_flow_m3_s[0]);
    }
    for changed in &flows[1..] {
        assert!((*changed - flows[0]).abs() > 1e-7 * flows[0].abs());
    }
}

#[test]
fn resolved_solid_binding_refuses_wrong_temperature_and_unowned_numeric_wall() {
    use super::super::thermoviscous::{ThermoviscousInterior, wall::MaterialWallRegion};
    use fs_feec::acoustics::thermoviscous::wall::ConductiveWallRegion;
    use fs_feec::fixtures::kuhn_cube;

    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let (mesh, unit) = kuhn_cube(1);
    let points: Vec<_> = unit.iter().map(|p| p.map(|x| 0.02 * x)).collect();
    let prepared = material(0.01, 1.0);
    let wrong = prepared.at_temperature(294.15).unwrap();
    let (fe, walls) = resolved_wall_fixture(&mesh, &points, &gas);
    let model = ThermoviscousInterior::from_gas(fe, gas, [100.0, 1000.0]).unwrap();
    assert!(
        model
            .with_conductive_walls(
                &[MaterialWallRegion {
                    faces: &walls,
                    sample: &wrong,
                    depth: ThermalWallDepth::HalfSpace,
                    minimum_tangential_scale_m: 0.02,
                    geometry: GEOMETRY,
                }],
                &mut || false,
            )
            .is_err()
    );

    let sample = prepared.at_temperature(gas.temperature).unwrap();
    let (fe, walls) = resolved_wall_fixture(&mesh, &points, &gas);
    let numerical = fe
        .with_conductive_walls(
            &[ConductiveWallRegion {
                faces: walls,
                wall: ConductiveThermalWall {
                    transport: conductive_wall_transport(&gas, &sample).unwrap(),
                    depth: ThermalWallDepth::HalfSpace,
                    minimum_tangential_scale_m: 0.02,
                },
                geometry: GEOMETRY,
                temperature_k: gas.temperature,
            }],
            &mut || false,
        )
        .unwrap();
    assert!(ThermoviscousInterior::from_gas(numerical, gas, [100.0, 1000.0]).is_err());
}
