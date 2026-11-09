//! Original material memory reaches distributed reflection and physical work.
//! Synthetic constitutive data exercise the coupling, not measured accuracy.

use fs_couple::bernoulli_aperture::wall::material::{
    AxialConstraint, LocalMaterialWall, WallMaterialLoss, relaxation::WallRelaxation,
};
use fs_couple::cavity_modal::interior::heterogeneous::AxialGasEmbedding;
use fs_couple::cavity_modal::interior::heterogeneous::mechanical::{
    CylindricalMaterialBoundary, MaterialMechanicalError, MaterialMechanicalHarmonic,
    MaterialMechanicalHarmonicError, MaterialMechanicalTime,
};
use fs_couple::cavity_modal::interior::heterogeneous::porous::{
    AcousticVolumeRegion, AcousticVolumeRegionState, PorousRegionsInterior,
};
use fs_couple::cavity_modal::interior::heterogeneous::regions::FluidRegionsInterior;
use fs_evidence::ValidityDomain;
use fs_feec::acoustics::helmholtz::heterogeneous::time::mechanical::{
    MechanicalBoundary, MechanicalTimeLimits,
};
use fs_feec::acoustics::helmholtz::heterogeneous::{
    BoundaryAdmittance, EquivalentFluidControls, HeterogeneousError, VelocityTrace,
};
use fs_feec::acoustics::helmholtz::{HelmholtzFe, Medium};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::fluid_acoustics::FluidAcoustics;
use fs_material::gas::{GasSpec, GasState};
use fs_material::porous::{PorousCard, PorousFluidLaw};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_math::c64::C64;
use fs_phs::acoustic_outlet::AcousticWallRelaxation;
use fs_phs::{WallPin, wall_specific_impedance};
use fs_qty::{Density, Dims, Pressure};

const RADIUS: f64 = 0.1;
const THICKNESS: f64 = 0.0005;
const DT: f64 = 1e-5;
const BAND: [f64; 2] = [1.0, 10_000.0];
const PLACEHOLDER: Medium = Medium {
    density_kg_m3: 7.0,
    sound_speed_m_s: 123.0,
};

fn source(thickness: f64, time_scale: f64, band: [f64; 2]) -> LocalMaterialWall {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, dims, value) in [
        ("density", Density::DIMS, 1000.0),
        ("young_modulus", Pressure::DIMS, 2e7),
        ("poisson_ratio", Dims([0; 6]), 0.3),
        ("e1", Pressure::DIMS, 1e7),
        ("t1", Dims([0, 0, 1, 0, 0, 0]), 0.002 * time_scale),
        ("e2", Pressure::DIMS, 4e7),
        ("t2", Dims([0, 0, 1, 0, 0, 0]), 0.02 * time_scale),
        ("e0", Pressure::DIMS, 0.0),
        ("t0", Dims([0, 0, 1, 0, 0, 0]), 0.03),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained()
                    .with("T", 290.0, 310.0)
                    .with("omega", 0.0, core::f64::consts::TAU * 100_000.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic cylindrical membrane fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic membrane".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let snapshot = resolve_material_state_point(
        &card,
        &QueryPoint::new()
            .with("T", 300.0)
            .unwrap()
            .with("omega", 100.0)
            .unwrap(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let memory = WallRelaxation::from_state_point(
        &snapshot,
        &[
            ("e1".into(), "t1".into()),
            ("e0".into(), "t0".into()),
            ("e2".into(), "t2".into()),
        ],
        band,
        0.25,
        0.4,
    )
    .unwrap();
    LocalMaterialWall::new(
        snapshot,
        thickness,
        AxialConstraint::StressFree,
        WallMaterialLoss::Maxwell(memory),
        band,
    )
    .unwrap()
}

// A short faceted cylindrical side, with the remaining cuts rigid.
// At its vertices the radial direction differs from the facet normal by 0.05 rad.
fn mesh() -> (fs_rep_mesh::TetComplex, Vec<[f64; 3]>) {
    let (complex, mut points) = fs_feec::fixtures::kuhn_cube(1);
    for point in &mut points {
        point[0] = RADIUS * 0.05_f64.cos() - 0.02 * (1.0 - point[0]);
        point[1] = RADIUS * 0.05_f64.sin() * (2.0 * point[1] - 1.0);
        point[2] *= 0.2;
    }
    (complex, points)
}

fn faces(geometry: &HelmholtzFe<'_>) -> Vec<[u32; 3]> {
    geometry
        .boundary_faces()
        .into_iter()
        .filter(|face| {
            face.iter()
                .all(|&v| geometry.positions()[v as usize][0] == RADIUS * 0.05_f64.cos())
        })
        .collect()
}

fn boundary<'s>(
    source: &'s LocalMaterialWall,
    faces: Vec<[u32; 3]>,
) -> CylindricalMaterialBoundary<'s> {
    CylindricalMaterialBoundary {
        faces,
        source,
        radius_m: RADIUS,
        relative_geometry_tolerance: 0.06,
        axis: AxialGasEmbedding::new([0.0; 3], [0.0, 0.0, 1.0]).unwrap(),
    }
}

fn close(actual: &[f64], expected: &[f64]) {
    let scale = expected.iter().map(|v| v.abs()).fold(1e-30_f64, f64::max);
    assert_eq!(actual.len(), expected.len());
    for (&a, &b) in actual.iter().zip(expected) {
        assert!((a - b).abs() < 2e-10 * scale, "{a} != {b}");
    }
}

/// G1: independently evaluated original constitutive impedance agrees with
/// every native pressure/wall coordinate away from mechanical poles.
#[test]
fn porous_volume_and_original_material_wall_share_pressure_motion_and_heat() {
    let wall = source(THICKNESS, 1.0, BAND);
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100_000.0).unwrap();
    let porous = PorousFluidLaw::try_new(PorousCard::melamine_uf_nguyen2024(), 1.81e-5).unwrap();
    let (complex, positions) = mesh();
    let geometry = HelmholtzFe::new(&complex, &positions, 2, PLACEHOLDER).unwrap();
    let selected = faces(&geometry);
    let regions = vec![0; complex.tets.len()];
    let sources = [AcousticVolumeRegion::Porous {
        law: &porous,
        fluid: &gas,
    }];
    let frequency_hz = 800.0;
    let model = MaterialMechanicalHarmonic::new(
        &geometry,
        frequency_hz,
        &sources,
        &regions,
        &[],
        vec![boundary(&wall, selected.clone())],
        &[],
        EquivalentFluidControls::default(),
        &mut || false,
    )
    .unwrap();
    assert!(core::ptr::eq(model.boundaries()[0].source, &wall));
    assert_eq!(model.wall_band_hz(), BAND);
    assert_eq!(model.interior().element_regions(), regions);
    let state = match model.interior().region_state(0).unwrap() {
        AcousticVolumeRegionState::Porous(state) => state,
        _ => panic!("original porous material was replaced"),
    };
    assert_eq!(
        state,
        porous
            .effective_fluid(&gas, core::f64::consts::TAU * frequency_hz)
            .unwrap()
    );
    let original = model.prepared().wall(0).unwrap();
    assert_eq!(original.relaxation.len(), 2);
    let WallMaterialLoss::Maxwell(memory) = model.boundaries()[0].source.loss() else {
        panic!()
    };
    assert_eq!(memory.material().terms.len(), 3);
    assert_eq!(memory.source_identity(), Some(wall.material().identity()));

    let minimum_x = positions
        .iter()
        .map(|point| point[0])
        .fold(f64::INFINITY, f64::min);
    let drive: Vec<VelocityTrace> = geometry
        .boundary_faces()
        .into_iter()
        .filter(|face| face.iter().all(|&v| positions[v as usize][0] == minimum_x))
        .map(|face| (face, [C64::new(0.001, -0.0003); 3]))
        .collect();
    let actual = model.prepared().solve(&drive, &mut || false).unwrap();
    // The source's original complex law uses exp(-i omega t); conjugate once
    // for this independently assembled exp(+i omega t) Robin reference.
    let impedance = wall
        .specific_impedance(RADIUS, core::f64::consts::TAU * frequency_hz)
        .unwrap()
        .conj();
    let reference = PorousRegionsInterior::new(
        &geometry,
        frequency_hz,
        &sources,
        &regions,
        &[BoundaryAdmittance {
            faces: selected,
            admittance_m_pa_s: impedance.recip(),
        }],
        &[],
        EquivalentFluidControls::default(),
        &mut || false,
    )
    .unwrap();
    let expected = reference.operator().solve(&drive, &mut || false).unwrap();
    let scale = expected
        .pressure_pa
        .iter()
        .map(|p| p.abs())
        .fold(0.0_f64, f64::max);
    for (&a, &b) in actual.pressure_pa.iter().zip(&expected.pressure_pa) {
        assert!((a - b).abs() < 2e-10 * scale, "{a:?} versus {b:?}");
    }
    for (&dof, &velocity) in model
        .prepared()
        .wall_dofs(0)
        .unwrap()
        .iter()
        .zip(&actual.walls[0].normal_velocity_m_s)
    {
        let expected_velocity = actual.pressure_pa[dof] / impedance;
        assert!((velocity - expected_velocity).abs() < 2e-10 * scale / impedance.abs());
    }
    assert_eq!(actual.walls[0].relaxation_memory.len(), 2);
    assert!(actual.wall_accounting.relaxation_dissipation_w > 0.0);
    assert!(actual.density_dissipation_w > 0.0 && actual.compressional_dissipation_w > 0.0);
    assert_eq!(actual.boundary_dissipation_w, 0.0);
    let power_scale = actual.input_power_w.abs();
    assert!(
        (actual.wall_accounting.dissipation_w() - expected.boundary_dissipation_w).abs()
            < 2e-10 * power_scale
    );
    assert!(
        (actual.density_dissipation_w - expected.density_dissipation_w).abs() < 2e-10 * power_scale
    );
    assert!(
        (actual.compressional_dissipation_w - expected.compressional_dissipation_w).abs()
            < 2e-10 * power_scale
    );
    assert!(actual.power_balance_residual_w().abs() < 2e-10 * power_scale);
    assert!(actual.fluid_power_balance_residual_w().abs() < 2e-10 * power_scale);
    assert!(actual.wall_accounting.power_balance_residual_w().abs() < 2e-10 * power_scale);
}

/// G0/G1: original source bandwidth and actual membrane geometry remain binding.
#[test]
fn porous_mechanical_material_binding_preserves_frequency_geometry_and_cancellation_refusals() {
    let wall = source(THICKNESS, 1.0, [100.0, 2000.0]);
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100_000.0).unwrap();
    let (complex, positions) = mesh();
    let geometry = HelmholtzFe::new(&complex, &positions, 1, PLACEHOLDER).unwrap();
    let regions = vec![0; complex.tets.len()];
    let sources = [AcousticVolumeRegion::Fluid(&gas)];
    let good = boundary(&wall, faces(&geometry));
    for frequency_hz in [0.0, 99.0, 2001.0, f64::NAN] {
        assert!(matches!(
            MaterialMechanicalHarmonic::new(
                &geometry,
                frequency_hz,
                &sources,
                &regions,
                &[],
                vec![good.clone()],
                &[],
                EquivalentFluidControls::default(),
                &mut || false,
            ),
            Err(MaterialMechanicalHarmonicError::Frequency { boundary: 0, .. })
        ));
    }
    let mut wrong_geometry = good.clone();
    wrong_geometry.radius_m *= 1.2;
    assert!(matches!(
        MaterialMechanicalHarmonic::new(
            &geometry,
            800.0,
            &sources,
            &regions,
            &[],
            vec![wrong_geometry],
            &[],
            EquivalentFluidControls::default(),
            &mut || false,
        ),
        Err(MaterialMechanicalHarmonicError::Wall(
            MaterialMechanicalError::Geometry { .. }
        ))
    ));
    assert!(matches!(
        MaterialMechanicalHarmonic::new(
            &geometry,
            800.0,
            &sources,
            &regions,
            &[],
            vec![good.clone()],
            &[],
            EquivalentFluidControls::default(),
            &mut || true,
        ),
        Err(MaterialMechanicalHarmonicError::Wall(
            MaterialMechanicalError::Interior(HeterogeneousError::Cancelled,)
        ))
    ));
    let thick = source(0.006, 1.0, BAND);
    assert!(matches!(
        MaterialMechanicalHarmonic::new(
            &geometry,
            800.0,
            &sources,
            &regions,
            &[],
            vec![boundary(&thick, good.faces)],
            &[],
            EquivalentFluidControls::default(),
            &mut || false,
        ),
        Err(MaterialMechanicalHarmonicError::Wall(
            MaterialMechanicalError::Source { .. }
        ))
    ));
}

/// G1/G3: original source identity, complex reflection and full transient parity.
#[test]
fn complete_material_spectrum_reaches_reflection_pressure_and_wall_work() {
    let source = source(THICKNESS, 1.0, BAND);
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100_000.0).unwrap();
    let (complex, positions) = mesh();
    let geometry = HelmholtzFe::new(&complex, &positions, 1, PLACEHOLDER).unwrap();
    let fluid = FluidRegionsInterior::new(
        &geometry,
        &[&gas],
        &vec![0; complex.tets.len()],
        3,
        &mut || false,
    )
    .unwrap();
    let selected = faces(&geometry);
    let bound = MaterialMechanicalTime::new(
        fluid.operator(),
        DT,
        &[],
        vec![boundary(&source, selected.clone())],
        MechanicalTimeLimits::default(),
        &mut || false,
    )
    .unwrap();
    assert!(std::ptr::eq(bound.boundaries()[0].source, &source));
    assert_eq!(bound.band_hz(), BAND);
    let WallMaterialLoss::Maxwell(memory) = bound.boundaries()[0].source.loss() else {
        panic!()
    };
    assert_eq!(memory.material().terms.len(), 3); // The exact-zero source arm is still recorded.
    assert_eq!(memory.source_identity(), Some(source.material().identity()));
    let wall = bound.prepared().wall(0).unwrap();
    assert_eq!(wall.relaxation.len(), 2);
    let fluid_state = gas.acoustic_state().unwrap();
    let z0 = C64::from_re(fluid_state.density_kg_m3() * fluid_state.sound_speed_m_s());
    for frequency in [10.0, 100.0, 400.0, 2000.0] {
        let omega = core::f64::consts::TAU * frequency;
        let base = wall_specific_impedance(&wall.wall, omega).unwrap();
        let actual = wall.relaxation.iter().fold(base, |z, arm| {
            z + C64::from_re(arm.stiffness_per_area)
                / C64 {
                    re: 1.0 / arm.relaxation_time_s,
                    im: -omega,
                }
        });
        let expected = source.specific_impedance(RADIUS, omega).unwrap();
        assert!((actual - expected).abs() < 5e-13 * expected.abs());
        assert!(
            (actual - base).abs() > 1.0,
            "memory must alter complex wall impedance"
        );
        let reflection = (actual - z0) / (actual + z0);
        let reference = (expected - z0) / (expected + z0);
        assert!((reflection - reference).abs() < 1e-12);
        assert!(reflection.abs() < 1.0);
    }

    // Independent synthetic coefficient image, not the adapter's output.
    let factor = THICKNESS / (RADIUS * RADIUS);
    let numerical = MechanicalBoundary {
        faces: selected,
        wall: WallPin {
            surface_density: 1000.0 * THICKNESS,
            stiffness_per_area: 2e7 * factor,
            resistance: 0.0,
        },
        relaxation: vec![
            AcousticWallRelaxation {
                stiffness_per_area: 1e7 * factor,
                relaxation_time_s: 0.002,
            },
            AcousticWallRelaxation {
                stiffness_per_area: 4e7 * factor,
                relaxation_time_s: 0.02,
            },
        ],
    };
    let direct = fluid
        .operator()
        .prepare_time_with_walls(
            DT,
            &[],
            &[numerical],
            MechanicalTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let n = fluid.operator().mass().nrows();
    let mut actual = bound
        .prepared()
        .state(&vec![0.0; n], &vec![1.0; n])
        .unwrap();
    let mut expected = direct.state(&vec![0.0; n], &vec![1.0; n]).unwrap();
    let initial = bound.prepared().energy(&actual).unwrap().total_j();
    let mut heat = 0.0;
    for _ in 0..80 {
        let a = bound.prepared().step(&actual, &[], &mut || false).unwrap();
        let b = direct.step(&expected, &[], &mut || false).unwrap();
        close(a.state.pressure_pa(), b.state.pressure_pa());
        close(
            a.state.wall_displacement_m(0).unwrap(),
            b.state.wall_displacement_m(0).unwrap(),
        );
        for arm in 0..2 {
            close(
                a.state.wall_relaxation_memory(0, arm).unwrap(),
                b.state.wall_relaxation_memory(0, arm).unwrap(),
            );
        }
        assert!(a.work.balance_residual_j().abs() < 1e-9 * initial);
        heat += a.work.wall_relaxation_dissipation_j;
        actual = a.state;
        expected = b.state;
    }
    let final_energy = bound.prepared().energy(&actual).unwrap();
    assert!(final_energy.wall_relaxation_j > 0.0 && heat > 0.0);
    assert!((final_energy.total_j() + heat - initial).abs() < 2e-9 * initial);
}

/// Original source bands and independent histories remain attached per patch.
#[test]
fn separate_original_wall_sources_keep_their_bands_and_all_arms() {
    let a = source(THICKNESS, 1.0, [1.0, 3000.0]);
    let b = source(THICKNESS, 2.0, [20.0, 5000.0]);
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100_000.0).unwrap();
    let (complex, positions) = mesh();
    let geometry = HelmholtzFe::new(&complex, &positions, 1, PLACEHOLDER).unwrap();
    let fluid = FluidRegionsInterior::new(
        &geometry,
        &[&gas],
        &vec![0; complex.tets.len()],
        3,
        &mut || false,
    )
    .unwrap();
    let faces = faces(&geometry);
    assert_eq!(faces.len(), 2);
    let bound = MaterialMechanicalTime::new(
        fluid.operator(),
        DT,
        &[],
        vec![boundary(&a, vec![faces[0]]), boundary(&b, vec![faces[1]])],
        MechanicalTimeLimits::default(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(bound.band_hz(), [20.0, 3000.0]);
    assert!(std::ptr::eq(bound.boundaries()[1].source, &b));
    assert_eq!(bound.prepared().wall_count(), 2);
    assert_eq!(
        bound.prepared().wall(1).unwrap().relaxation[0].relaxation_time_s,
        0.004
    );
    let other = source(THICKNESS, 1.0, [4000.0, 5000.0]);
    assert!(matches!(
        MaterialMechanicalTime::new(
            fluid.operator(),
            DT,
            &[],
            vec![
                boundary(&a, vec![faces[0]]),
                boundary(&other, vec![faces[1]])
            ],
            MechanicalTimeLimits::default(),
            &mut || false
        ),
        Err(MaterialMechanicalError::Geometry { .. })
    ));
}

/// Geometry, constitutive accuracy and native resource refusals stay authoritative.
#[test]
fn incompatible_geometry_underresolved_memory_and_native_limits_refuse() {
    let original = source(THICKNESS, 1.0, BAND);
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100_000.0).unwrap();
    let (complex, positions) = mesh();
    let geometry = HelmholtzFe::new(&complex, &positions, 1, PLACEHOLDER).unwrap();
    let fluid = FluidRegionsInterior::new(
        &geometry,
        &[&gas],
        &vec![0; complex.tets.len()],
        3,
        &mut || false,
    )
    .unwrap();
    let good = boundary(&original, faces(&geometry));
    let limits = MechanicalTimeLimits::default();
    let mut wrong_radius = good.clone();
    wrong_radius.radius_m *= 2.0;
    let mut unresolved_facets = good.clone();
    unresolved_facets.relative_geometry_tolerance = 0.001;
    let mut cap = good.clone();
    cap.faces = geometry
        .boundary_faces()
        .into_iter()
        .filter(|face| face.iter().all(|&v| positions[v as usize][2] == 0.0))
        .collect();
    for invalid in [wrong_radius, unresolved_facets, cap] {
        assert!(matches!(
            MaterialMechanicalTime::new(
                fluid.operator(),
                DT,
                &[],
                vec![invalid],
                limits,
                &mut || false
            ),
            Err(MaterialMechanicalError::Geometry { .. })
        ));
    }
    assert!(matches!(
        MaterialMechanicalTime::new(
            fluid.operator(),
            0.0006,
            &[],
            vec![good.clone()],
            limits,
            &mut || false
        ),
        Err(MaterialMechanicalError::Source { .. })
    ));
    let no_memory = MechanicalTimeLimits {
        max_wall_relaxation_dofs: 0,
        ..limits
    };
    assert!(matches!(
        MaterialMechanicalTime::new(
            fluid.operator(),
            DT,
            &[],
            vec![good.clone()],
            no_memory,
            &mut || false
        ),
        Err(MaterialMechanicalError::Interior(
            HeterogeneousError::Limit { .. }
        ))
    ));
    assert!(matches!(
        MaterialMechanicalTime::new(
            fluid.operator(),
            DT,
            &[],
            vec![good.clone()],
            limits,
            &mut || true
        ),
        Err(MaterialMechanicalError::Interior(
            HeterogeneousError::Cancelled
        ))
    ));
    let thick = source(0.006, 1.0, BAND);
    assert!(matches!(
        MaterialMechanicalTime::new(
            fluid.operator(),
            DT,
            &[],
            vec![boundary(&thick, good.faces)],
            limits,
            &mut || false
        ),
        Err(MaterialMechanicalError::Source { .. })
    ));
}
