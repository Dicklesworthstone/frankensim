//! Synthetic source cards exercise the original shell/material/FE owners.
//! These are discrete physical regressions, not measured-shell validation.
use super::super::{source_fluid_medium, source_medium};
use super::*;
use crate::shear_shell::{
    ResolvedShearShell, ShellProjectionBudget, ShellReductionOptions, ShellRegionMaterial,
    TRANSVERSE_SHEAR_VISCOSITY_PROPERTIES, compile_material_shear_shell,
    connections::{
        MaterialShellConnection, ShellInterfacePoint, ShellInterfaceSide, stiffness_per_area_dims,
        viscosity_per_area_dims,
    },
    pressure::ShellPressureSurfaces,
};
use crate::thin_plate::{PlateMaterialModel, PlateRegionMaterial, PlateThicknessConstraint};
use fs_evidence::ValidityDomain;
use fs_feec::acoustics::helmholtz::Medium;
use fs_feec::acoustics::thermoviscous::{ThermoviscousField, time::ThermoviscousTimeLimits};
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
    UncertaintyModel,
};
use fs_material::gas::GasSpec;
use fs_material::plane_stress_viscosity::PLANE_STRESS_VISCOSITY_PROPERTIES;
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_interface_state_point, resolve_material_state_point,
};
use fs_plate::{
    PlateRegion,
    shell::{
        ShellMesh,
        mitc::{
            ConstraintFrame, MitcBudget, NodalDirectors, ShellConstraint, ShellMass,
            connections::ShellConnectionBudget,
        },
    },
};
use fs_qty::{Density, Dims, DynViscosity, Pressure};
use fs_rep_mesh::TetComplex;

mod relaxation;

const BAND: [f64; 2] = [20.0, 10_000.0];
const POINTS: [[f64; 3]; 6] = [
    [0.0, 0.0, 0.0],
    [0.02, 0.0, 0.0],
    [0.02, 0.02, 0.0],
    [0.0, 0.02, 0.0],
    [0.01, 0.01, 0.0],
    [0.01, 0.01, 0.01],
];
const TRIANGLES: [[usize; 3]; 4] = [[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]];

fn claims(rows: &[(&str, Dims, f64)]) -> (ClaimSet, Vec<ScalarPropertyRequirement>) {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for &(name, dims, value) in rows {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained()
                    .with("T", 280.0, 310.0)
                    .with("omega", 1.0, 1e6),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic conforming shell/fluid specimen".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap(),
        );
    }
    (claims, requirements)
}

fn material(anchored: bool) -> ResolvedMaterialStatePoint {
    let mut rows = vec![
        ("density", Density::DIMS, 800.0),
        ("young_modulus_1", Pressure::DIMS, 1e7),
        ("young_modulus_2", Pressure::DIMS, 4e6),
        ("poisson_ratio_12", Dims::NONE, 0.25),
        ("shear_modulus_12", Pressure::DIMS, 1.8e6),
        ("shear_modulus_13", Pressure::DIMS, 8e5),
        ("shear_modulus_23", Pressure::DIMS, 3e5),
    ];
    rows.extend(
        PLANE_STRESS_VISCOSITY_PROPERTIES
            .into_iter()
            .zip([100.0, 70.0, 30.0, 20.0, 4.0, -3.0])
            .map(|(name, value)| (name, DynViscosity::DIMS, value)),
    );
    rows.extend(
        TRANSVERSE_SHEAR_VISCOSITY_PROPERTIES
            .into_iter()
            .zip([40.0, 80.0, 5.0])
            .map(|(name, value)| (name, DynViscosity::DIMS, value)),
    );
    if anchored {
        rows.push((
            crate::shear_shell::CONSTANT_LOSS_REFERENCE_OMEGA,
            crate::shear_shell::CONSTANT_LOSS_REFERENCE_OMEGA_DIMS,
            1000.0,
        ));
    }
    let (claims, requirements) = claims(&rows);
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic orthotropic shell".into(),
            phase: "solid".into(),
            process: "physical regression".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_material_state_point(
        &card,
        &QueryPoint::new()
            .with("T", 293.15)
            .unwrap()
            .with("omega", 1000.0)
            .unwrap(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn source(
    permutation: [usize; 5],
    side: ShellFace,
    anchored: bool,
    gate: &CancelGate,
) -> (
    ResolvedShearShell,
    NodalDirectors,
    ShellReductionOptions,
    Vec<Vec<f64>>,
) {
    let mut points = vec![[0.0; 3]; 5];
    for (old, &new) in permutation.iter().enumerate() {
        points[new] = POINTS[old];
    }
    let mesh = ShellMesh::new(
        points,
        TRIANGLES
            .iter()
            .map(|tri| tri.map(|n| permutation[n]))
            .collect(),
    )
    .unwrap();
    let budget = MitcBudget {
        max_nodes: 5,
        max_triangles: 4,
    };
    let material = material(anchored);
    let shell = compile_material_shear_shell(
        mesh,
        vec![ShellRegionMaterial {
            region: PlateRegionMaterial {
                region: PlateRegion {
                    name: "panel".into(),
                    triangle_indices: (0..4).collect(),
                },
                material: material.clone(),
                model: PlateMaterialModel::Orthotropic12,
                thickness_constraint: PlateThicknessConstraint::FixedThickness(0.0002),
                material_angle_rad: 0.31,
            },
            shear_correction: 5.0 / 6.0,
            axis1_global: [1.0, 0.0, 0.0],
        }],
        budget,
        gate,
    )
    .unwrap();
    let mut fixed = material.material().clone();
    fixed.process = "immobile interface counterpart".into();
    let (interface_claims, requirements) = claims(&[
        ("spring", stiffness_per_area_dims(), 1e8),
        ("dashpot", viscosity_per_area_dims(), 5000.0),
    ]);
    let card = InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: material.material().clone(),
            texture_frame: "center".into(),
        },
        SurfaceSpec {
            material: fixed.clone(),
            texture_frame: "fixed".into(),
        },
        SystemContext {
            medium: "synthetic".into(),
            third_body: None,
            environment: "fixed".into(),
            history: "linear".into(),
        },
        interface_claims,
        vec![],
    )
    .unwrap();
    let interface = resolve_interface_state_point(
        &card,
        &material.query(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let shell = shell
        .with_interfaces(
            vec![MaterialShellConnection {
                name: "center dashpot".into(),
                interface,
                stiffness_property: "spring".into(),
                viscosity_property: "dashpot".into(),
                area_m2: 1e-5,
                direction: [1.0, 0.0, 0.0],
                left: ShellInterfaceSide::Region {
                    name: "panel".into(),
                    points: vec![ShellInterfacePoint {
                        node: permutation[4],
                        weight: 1.0,
                        offset_m: [0.0; 3],
                    }],
                },
                right: ShellInterfaceSide::Fixed { material: fixed },
            }],
            ShellConnectionBudget {
                max_connections: 1,
                max_terms: 1,
                max_work: 100_000,
            },
            gate,
        )
        .unwrap();
    let selected: Vec<_> = (0..4)
        .map(|triangle| OrientedShellFacet { triangle, side })
        .collect();
    let shell = shell
        .with_pressure_surfaces(
            ShellPressureSurfaces {
                exterior: selected.clone(),
                cavity: selected,
            },
            gate,
        )
        .unwrap();
    let directors = NodalDirectors::new(shell.mesh(), vec![[0.0, 0.0, 1.0]; 5], 0.9).unwrap();
    let mut constraints: Vec<_> = (0..4)
        .map(|old| ShellConstraint {
            node: permutation[old],
            mask: [true; 6],
            frame: ConstraintFrame::Global,
        })
        .collect();
    constraints.push(ShellConstraint {
        node: permutation[4],
        mask: [false, false, true, true, true, true],
        frame: ConstraintFrame::Global,
    });
    let options = ShellReductionOptions {
        constraints,
        mass: ShellMass::Consistent,
        drilling_alpha: 1e-3,
        eigenvalue_window: (1e-6, 1e10),
        n_modes: 2,
        assembly_budget: budget,
        projection_budget: ShellProjectionBudget {
            max_modes: 2,
            max_projection_terms: 10_000_000,
            reference_resistance_n_s_m: 7.0,
        },
    };
    let ports = [0, 1]
        .map(|component| {
            let mut row = vec![0.0; 30];
            row[6 * permutation[4] + component] = 1.0;
            row
        })
        .to_vec();
    (shell, directors, options, ports)
}

fn image(permutation: [usize; 5], side: ShellFace, anchored: bool) -> MaterialShellImage {
    let gate = CancelGate::new_clock_free();
    let (source, directors, options, ports) = source(permutation, side, anchored, &gate);
    source.reduce(&directors, &options, &ports, &gate).unwrap()
}

fn fluid(permutation: [usize; 6]) -> (TetComplex, Vec<[f64; 3]>) {
    let mut positions = vec![[0.0; 3]; 6];
    for (old, &new) in permutation.iter().enumerate() {
        positions[new] = POINTS[old];
    }
    let mesh = TetComplex::from_tets(
        6,
        TRIANGLES
            .iter()
            .map(|tri| {
                [
                    permutation[tri[0]] as u32,
                    permutation[tri[1]] as u32,
                    permutation[tri[2]] as u32,
                    permutation[5] as u32,
                ]
            })
            .collect(),
    );
    (mesh, positions)
}
fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}
fn controls() -> ShellWallControls {
    ShellWallControls {
        maximum_source_nodes: 5,
        maximum_source_triangles: 4,
        maximum_source_interfaces: 1,
        maximum_force_ports: 2,
        maximum_modes: 2,
        maximum_damping_channels: 3,
        relaxation: MovingWallRelaxationLimits {
            maximum_channels: 32,
            maximum_workspace_bytes: 1_000_000,
            maximum_projection_terms: 10_000,
        },
        maximum_setup_bytes: 1_000_000,
        maximum_setup_work: 10_000_000,
        sample_rate_hz: 48_000,
        modal: ModalAcousticTimeBudget::audible_reference(),
        fluid: ThermoviscousLimits::default(),
    }
}
fn map(image: &MaterialShellImage, positions: &[[f64; 3]]) -> Vec<Option<u32>> {
    image
        .source()
        .mesh()
        .nodes
        .iter()
        .map(|point| Some(positions.iter().position(|p| p == point).unwrap() as u32))
        .collect()
}
fn geometry<'c>(
    mesh: &'c TetComplex,
    positions: &[[f64; 3]],
    medium: ThermoviscousMedium,
    thermo: ThermoviscousThermodynamics,
) -> HelmholtzFe<'c> {
    HelmholtzFe::new(
        mesh,
        positions,
        2,
        Medium {
            density_kg_m3: medium.density_kg_m3,
            sound_speed_m_s: thermo.sound_speed_m_s(medium).unwrap(),
        },
    )
    .unwrap()
}
fn boundary(mesh: &TetComplex) -> ThermoviscousBoundary {
    // Only exterior faces: each base facet plus the four triangular sidewalls.
    let mut faces = Vec::new();
    for face in &mesh.faces {
        if mesh
            .tets
            .iter()
            .filter(|tet| face.iter().all(|v| tet.contains(v)))
            .count()
            == 1
        {
            faces.push(*face);
        }
    }
    ThermoviscousBoundary {
        no_slip_faces: faces,
        isothermal_faces: vec![],
        ports: vec![],
    }
}
fn assembly<'s, 'c>(
    image: &'s MaterialShellImage,
    mesh: &'c TetComplex,
    positions: &[[f64; 3]],
) -> ShellWallAssembly<'s, 'c> {
    let medium = source_medium(&gas(), 0.0).unwrap();
    let thermo = ThermoviscousThermodynamics::ideal_gas(medium).unwrap();
    ShellWallAssembly::new(
        image,
        geometry(mesh, positions, medium, thermo),
        &map(image, positions),
        medium,
        thermo,
        &boundary(mesh),
        controls(),
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}

fn manual<'c>(
    image: &MaterialShellImage,
    mesh: &'c TetComplex,
    positions: &[[f64; 3]],
) -> ThermoviscousFe<'c> {
    let medium = source_medium(&gas(), 0.0).unwrap();
    let thermo = ThermoviscousThermodynamics::ideal_gas(medium).unwrap();
    let mapped = map(image, positions);
    let shapes = image
        .physical_modes()
        .iter()
        .map(|mode| {
            let mut shape = vec![[0.0; 3]; mesh.vertex_count + mesh.edges.len()];
            for (node, &vertex) in mapped.iter().enumerate() {
                shape[vertex.unwrap() as usize] =
                    [mode[6 * node], mode[6 * node + 1], mode[6 * node + 2]];
            }
            shape
        })
        .collect();
    let faces = image
        .source()
        .mesh()
        .tris
        .iter()
        .map(|tri| canonical(tri.map(|node| mapped[node].unwrap())))
        .collect();
    let modes = image
        .radiators()
        .iter()
        .map(|mode| ModalAcousticMode {
            angular_frequency_rad_s: mode.omega,
            damping_ratio: 0.0,
            pressure_per_modal_velocity: C64::ZERO,
        })
        .collect();
    let mut channels: Vec<Vec<f64>> = image
        .damping()
        .connections(0)
        .iter()
        .map(|connection| {
            connection
                .left
                .shapes
                .iter()
                .map(|&value| connection.damping_n_s_m.sqrt() * value)
                .collect()
        })
        .collect();
    channels.extend(
        image
            .interface_damping_connections(0)
            .iter()
            .map(|connection| {
                connection
                    .left
                    .shapes
                    .iter()
                    .map(|&value| connection.damping_n_s_m.sqrt() * value)
                    .collect::<Vec<_>>()
            }),
    );
    ThermoviscousFe::new_with_moving_wall(
        mesh,
        positions,
        medium,
        thermo,
        &boundary(mesh),
        MovingWall {
            faces,
            velocity_shapes: shapes,
            structure: ModalAcousticTimeModel::try_new(
                48_000,
                modes,
                ModalAcousticTimeBudget::audible_reference(),
            )
            .unwrap(),
        },
        controls().fluid,
        &mut || false,
    )
    .unwrap()
    .with_moving_wall_damping_channels(channels, &mut || false)
    .unwrap()
}
fn close(a: &[C64], b: &[C64], tolerance: f64) {
    assert_eq!(a.len(), b.len());
    let scale = b.iter().map(|x| x.abs()).fold(1e-30, f64::max);
    assert!(
        a.iter()
            .zip(b)
            .all(|(x, y)| (*x - *y).abs() < tolerance * scale)
    );
}
fn response(
    source: &ShellFluidInterior<'_, '_>,
    gate: &CancelGate,
) -> (ThermoviscousField, Vec<C64>) {
    let mut force = vec![C64::ZERO; source.image().physical_modes().len()];
    source
        .project_forces(
            &[C64::new(1e-4, 2e-5), C64::new(-3e-5, 4e-5)],
            &mut force,
            gate,
        )
        .unwrap();
    let field = source
        .interior()
        .response_with_wall_forces(1000.0, &[], &force, &mut || false)
        .unwrap();
    (field, force)
}

#[test]
fn original_shell_matches_native_loading_and_retains_zero_volume_shear_motion() {
    let gate = CancelGate::new_clock_free();
    let image = image([0, 1, 2, 3, 4], ShellFace::Negative, false);
    let (mesh, positions) = fluid([0, 1, 2, 3, 4, 5]);
    let source = assembly(&image, &mesh, &positions)
        .bind_gas(gas(), BAND, &gate)
        .unwrap();
    assert!(std::ptr::eq(source.image(), &image));
    assert!(image.physical_modes().iter().flatten().any(|v| *v < 0.0));
    assert!(
        image
            .cavity_volume_shapes_m2_per_sqrt_kg()
            .iter()
            .all(|v| *v == 0.0)
    );
    let wall = source.interior().inner().moving_wall().unwrap();
    for (r, shape) in wall.velocity_shapes.iter().enumerate() {
        for node in 0..5 {
            assert_eq!(
                shape[node],
                image.physical_modes()[r][6 * node..6 * node + 3]
            );
        }
        assert!(shape[mesh.vertex_count..].iter().all(|v| *v == [0.0; 3]));
        assert_eq!(
            wall.structure.modes()[r].angular_frequency_rad_s,
            image.radiators()[r].omega
        );
        assert_eq!(wall.structure.damping_rate_per_s(r), 0.0);
    }
    let actual_c = source
        .interior()
        .inner()
        .moving_wall_damping_matrix_per_s()
        .unwrap();
    let interface = &image.source().interfaces()[0];
    let column = &image.interface_columns()[0];
    for i in 0..2 {
        for j in 0..2 {
            let expected = image.damping().matrix_per_s()[2 * i + j]
                + interface.law().damping_n_s_m * column[i] * column[j];
            assert!((actual_c[2 * i + j] - expected).abs() < 2e-11 * expected.abs().max(1e-30));
        }
    }
    let (field, force) = response(&source, &gate);
    let original = manual(&image, &mesh, &positions);
    let expected = original
        .prepare(1000.0, &mut || false)
        .unwrap()
        .response_with_wall_forces(&[], &force, &mut || false)
        .unwrap();
    close(&field.pressure_pa, &expected.pressure_pa, 2e-10);
    close(
        &field.moving_wall.modal_velocity_m_sqrt_kg_per_s,
        &expected.moving_wall.modal_velocity_m_sqrt_kg_per_s,
        2e-10,
    );
    close(
        &field.moving_wall.fluid_force_n_per_sqrt_kg,
        &expected.moving_wall.fluid_force_n_per_sqrt_kg,
        2e-10,
    );
    assert!(
        field
            .moving_wall
            .fluid_force_n_per_sqrt_kg
            .iter()
            .any(|f| f.abs() > 1e-12)
    );
    assert!(field.power.viscous_w > 0.0);
    assert!(field.power.structural_damping_w > 0.0);
    assert!(field.power.balance_residual_w().abs() < 2e-8 * field.power.input_w);
    let physical_work: f64 = image
        .force_ports()
        .iter()
        .zip([C64::new(1e-4, 2e-5), C64::new(-3e-5, 4e-5)])
        .map(|(row, amplitude)| {
            let velocity = row
                .iter()
                .zip(&field.moving_wall.modal_velocity_m_sqrt_kg_per_s)
                .fold(C64::ZERO, |v, (&shape, &u)| v + u.scale(shape));
            0.5 * (amplitude.conj() * velocity).re
        })
        .sum();
    assert!((physical_work - field.moving_wall.external_input_w).abs() < 2e-12 * physical_work);

    let time = source
        .interior()
        .prepare_time(1e-5, ThermoviscousTimeLimits::default(), &mut || false)
        .unwrap();
    let reference_time = original
        .prepare_time(1e-5, ThermoviscousTimeLimits::default(), &mut || false)
        .unwrap();
    let real_force: Vec<_> = force.iter().map(|f| f.re).collect();
    let step = time
        .step(&time.zero_state(), &[], &real_force, &mut || false)
        .unwrap();
    let reference = reference_time
        .step(&reference_time.zero_state(), &[], &real_force, &mut || {
            false
        })
        .unwrap();
    for (a, b) in step
        .state
        .wall_velocity_m_sqrt_kg_per_s
        .iter()
        .zip(&reference.state.wall_velocity_m_sqrt_kg_per_s)
    {
        assert!((a - b).abs() < 2e-10 * b.abs().max(1e-30));
    }
    assert!(step.work.structural_damping_heat_j > 0.0);
    assert!(
        (step.work.structural_damping_heat_j / reference.work.structural_damping_heat_j - 1.0)
            .abs()
            < 2e-9
    );
    assert!(step.work.balance_residual_j().abs() < 2e-8 * step.work.supplied_j());
}

#[test]
fn source_and_fluid_relabeling_preserve_physical_work_and_fields() {
    let gate = CancelGate::new_clock_free();
    let first = image([0, 1, 2, 3, 4], ShellFace::Negative, false);
    let second = image([2, 0, 4, 1, 3], ShellFace::Negative, false);
    let (mesh_a, points_a) = fluid([0, 1, 2, 3, 4, 5]);
    let permutation = [3, 5, 1, 0, 4, 2];
    let (mesh_b, points_b) = fluid(permutation);
    let a = assembly(&first, &mesh_a, &points_a)
        .bind_gas(gas(), BAND, &gate)
        .unwrap();
    let b = assembly(&second, &mesh_b, &points_b)
        .bind_gas(gas(), BAND, &gate)
        .unwrap();
    let (a, _) = response(&a, &gate);
    let (b, _) = response(&b, &gate);
    let reordered: Vec<_> = permutation.iter().map(|&new| b.pressure_pa[new]).collect();
    close(&a.pressure_pa, &reordered, 2e-7);
    for (old, &new) in permutation.iter().enumerate() {
        close(&a.velocity_m_s[old], &b.velocity_m_s[new], 2e-7);
    }
    assert!((a.power.structural_damping_w / b.power.structural_damping_w - 1.0).abs() < 2e-7);
    assert!((a.power.input_w / b.power.input_w - 1.0).abs() < 2e-7);
}

#[test]
fn generic_material_binding_and_force_projection_preserve_the_original_owners() {
    let gate = CancelGate::new_clock_free();
    let image = image([0, 1, 2, 3, 4], ShellFace::Negative, false);
    let (mesh, positions) = fluid([0, 1, 2, 3, 4, 5]);
    let material = gas();
    let (medium, thermo) = source_fluid_medium(&material, material.pressure, 0.0).unwrap();
    let source = ShellWallAssembly::new(
        &image,
        geometry(&mesh, &positions, medium, thermo),
        &map(&image, &positions),
        medium,
        thermo,
        &boundary(&mesh),
        controls(),
        &gate,
    )
    .unwrap()
    .bind_material(material, BAND, &gate)
    .unwrap();
    assert_eq!(source.interior().inner().medium(), medium);
    assert_eq!(source.interior().inner().thermodynamics(), thermo);
    assert_eq!(
        source.interior().material_source().density,
        material.density
    );
    let mut forces = [C64::new(7.0, 8.0); 2];
    let before = forces;
    assert!(
        source
            .project_forces(&[C64::ONE], &mut forces, &gate)
            .is_err()
    );
    assert_eq!(forces, before);
    assert!(
        source
            .project_forces(&[C64::new(f64::NAN, 0.0), C64::ONE], &mut forces, &gate)
            .is_err()
    );
    assert_eq!(forces, before);
    gate.request();
    assert!(matches!(
        source.project_forces(&[C64::ONE; 2], &mut forces, &gate),
        Err(ShellWallError::Cancelled)
    ));
    assert_eq!(forces, before);
    let continued = source.into_interior();
    assert!(continued.inner().moving_wall().is_some());
    assert!(
        continued
            .inner()
            .moving_wall_damping_matrix_per_s()
            .is_some()
    );
}

#[test]
fn stale_or_incomplete_geometry_unsupported_laws_and_resources_refuse() {
    let gate = CancelGate::new_clock_free();
    let image = image([0, 1, 2, 3, 4], ShellFace::Negative, false);
    let (mesh, positions) = fluid([0, 1, 2, 3, 4, 5]);
    let medium = source_medium(&gas(), 0.0).unwrap();
    let thermo = ThermoviscousThermodynamics::ideal_gas(medium).unwrap();
    let build = |source,
                 xyz: &[[f64; 3]],
                 mapping: &[Option<u32>],
                 boundary: &ThermoviscousBoundary,
                 limits| {
        ShellWallAssembly::new(
            source,
            geometry(&mesh, xyz, medium, thermo),
            mapping,
            medium,
            thermo,
            boundary,
            limits,
            &gate,
        )
        .map(|_| ())
    };
    let mapped = map(&image, &positions);
    let mut missing = mapped.clone();
    missing[4] = None;
    assert!(build(&image, &positions, &missing, &boundary(&mesh), controls()).is_err());
    let mut welded = mapped.clone();
    welded[4] = welded[0];
    assert!(build(&image, &positions, &welded, &boundary(&mesh), controls()).is_err());
    let mut stale = positions.clone();
    stale[4][0] += 1e-8;
    assert!(build(&image, &stale, &mapped, &boundary(&mesh), controls()).is_err());
    let mut absent_face = boundary(&mesh);
    absent_face.no_slip_faces.retain(|face| *face != [0, 1, 4]);
    assert!(build(&image, &positions, &mapped, &absent_face, controls()).is_err());
    for changed in 0..6 {
        let mut limits = controls();
        match changed {
            0 => limits.maximum_modes = 1,
            1 => limits.maximum_damping_channels = 2,
            2 => limits.maximum_setup_bytes = 1,
            3 => limits.maximum_setup_work = 1,
            4 => limits.maximum_source_interfaces = 0,
            _ => limits.maximum_force_ports = 1,
        }
        assert!(matches!(
            build(&image, &positions, &mapped, &boundary(&mesh), limits),
            Err(ShellWallError::Budget { .. })
        ));
    }
    let wrong_side = self::image([0, 1, 2, 3, 4], ShellFace::Positive, false);
    assert!(
        build(
            &wrong_side,
            &positions,
            &mapped,
            &boundary(&mesh),
            controls()
        )
        .is_err()
    );
    let anchored = self::image([0, 1, 2, 3, 4], ShellFace::Negative, true);
    assert!(build(&anchored, &positions, &mapped, &boundary(&mesh), controls()).is_err());
    assert!(
        assembly(&image, &mesh, &positions)
            .bind_gas(gas(), [0.01, 1000.0], &gate)
            .is_err()
    );
    gate.request();
    assert!(matches!(
        build(&image, &positions, &mapped, &boundary(&mesh), controls()),
        Err(ShellWallError::Cancelled)
    ));
}
