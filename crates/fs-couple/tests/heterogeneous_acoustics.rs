//! Original source packs reach spatial acoustic operators and driven pressure.
//! Synthetic constitutive data exercise the binding, not measured-gas accuracy.

use fs_couple::cavity_modal::interior::heterogeneous::{
    AxialGasEmbedding, AxialGasInterior, AxialGasInteriorError,
};
use fs_feec::acoustics::helmholtz::heterogeneous::{BoundaryAdmittance, HarmonicLimits};
use fs_feec::acoustics::helmholtz::{HelmholtzFe, Medium};
use fs_feec::fixtures::{kuhn_cube, single_tet};
use fs_matdb::QueryPoint;
use fs_material::gas_field::{AxialGasField, PrimitiveGasState};
use fs_material::gas_species::input::ThermochemicalPhase;
use fs_material::thermochemical_gas::packs::from_packs;
use fs_math::c64::C64;
use fs_sparse::Csr;
use fs_thermochem::{Composition, CompositionBasis};

#[path = "../src/render/schedule/reed/field_input/tests/support.rs"]
mod gas_support;

fn primitive(temperature_k: f64, a: f64) -> PrimitiveGasState {
    PrimitiveGasState {
        temperature_k,
        composition: Composition::new(CompositionBasis::MoleFraction, vec![a, 1.0 - a]).unwrap(),
    }
}

fn field(stations: Vec<f64>, states: Vec<PrimitiveGasState>, cellwise: bool) -> AxialGasField {
    let a = gas_support::packs("A", 0.028, 3.5, 200.0);
    let b = gas_support::packs("B", 0.044, 4.5, 200.0);
    let model = from_packs(
        &[a.binding(), b.binding()],
        primitive(300.0, 0.5).composition,
        None,
    )
    .unwrap();
    let make = if cellwise {
        AxialGasField::new_piecewise_constant
    } else {
        AxialGasField::new
    };
    make(
        stations,
        100_000.0,
        states,
        model,
        QueryPoint::new(),
        ThermochemicalPhase::FrozenGas,
    )
    .unwrap()
}

fn axis() -> AxialGasEmbedding {
    AxialGasEmbedding::new([0.0; 3], [1.0, 0.0, 0.0]).unwrap()
}

// Deliberately unrelated: only geometry is borrowed by the new source owner.
const PLACEHOLDER: Medium = Medium {
    density_kg_m3: 1000.0,
    sound_speed_m_s: 1500.0,
};

fn dense(matrix: &Csr) -> Vec<f64> {
    let mut values = vec![0.0; matrix.nrows() * matrix.ncols()];
    for row in 0..matrix.nrows() {
        let (columns, entries) = matrix.row(row);
        for (&column, &entry) in columns.iter().zip(entries) {
            values[row * matrix.ncols() + column] = entry;
        }
    }
    values
}

#[test]
fn uniform_original_source_recovers_physical_matrices_without_reference_fluid_scaling() {
    let (complex, positions) = single_tet();
    let geometry = HelmholtzFe::new(&complex, &positions, 2, PLACEHOLDER).unwrap();
    let source = field(vec![0.0, 1.0], vec![primitive(300.0, 0.75); 2], false);
    let bound = AxialGasInterior::new(&geometry, &source, axis(), 4, &mut || false).unwrap();
    assert!(std::ptr::eq(bound.source(), &source));
    assert!(std::ptr::eq(bound.operator().geometry(), &geometry));
    let gas = source.uniform_state().unwrap().state();
    let probe = bound.probe(0, [0.25; 4]).unwrap();
    assert_eq!(probe.medium(), Medium {
        density_kg_m3: gas.density,
        sound_speed_m_s: gas.sound_speed,
    });
    assert!(matches!(bound.probe(1, [0.25; 4]), Err(AxialGasInteriorError::Interior(_))));
    assert!(matches!(bound.probe(0, [0.2; 4]), Err(AxialGasInteriorError::Interior(_))));
    for (physical, geometric, scale) in [
        (
            bound.operator().stiffness(),
            geometry.stiffness(),
            1.0 / gas.density,
        ),
        (
            bound.operator().mass(),
            geometry.mass(),
            1.0 / (gas.density * gas.sound_speed * gas.sound_speed),
        ),
    ] {
        let expected = dense(geometric)
            .into_iter()
            .map(|value| value * scale)
            .collect::<Vec<_>>();
        let norm = expected
            .iter()
            .copied()
            .map(f64::abs)
            .fold(0.0_f64, f64::max);
        for (actual, expected) in dense(physical).into_iter().zip(expected) {
            assert!((actual - expected).abs() < 2e-12 * norm);
        }
    }
}

#[test]
fn primitive_gradient_changes_driven_pressure_and_retains_boundary_work() {
    let (complex, positions) = kuhn_cube(2);
    let geometry = HelmholtzFe::new(&complex, &positions, 1, PLACEHOLDER).unwrap();
    let cold = primitive(300.0, 0.875);
    let sources = [
        field(vec![0.0, 1.0], vec![cold.clone(); 2], false),
        field(vec![0.0, 1.0], vec![cold, primitive(700.0, 0.125)], false),
    ];
    let faces = geometry.boundary_faces();
    let source_faces: Vec<_> = faces
        .iter()
        .copied()
        .filter(|face| face.iter().all(|&v| positions[v as usize][0] == 0.0))
        .collect();
    let loads = source_faces
        .iter()
        .map(|&face| (face, [C64::from_re(0.001); 3]))
        .collect::<Vec<_>>();
    let mut responses = Vec::new();
    for source in &sources {
        let bound = AxialGasInterior::new(&geometry, source, axis(), 5, &mut || false).unwrap();
        let prepared = bound
            .operator()
            .prepare(
                90.0,
                &[BoundaryAdmittance {
                    faces: faces.clone(),
                    admittance_m_pa_s: C64::from_re(0.002),
                }],
                HarmonicLimits::default(),
                &mut || false,
            )
            .unwrap();
        let result = prepared.solve(&loads, &mut || false).unwrap();
        assert!(result.input_power_w > 0.0);
        assert!(result.kinetic_energy_j > 0.0 && result.compressional_energy_j > 0.0);
        assert!(result.power_balance_residual_w().abs() < 1e-9 * result.input_power_w);
        for element in [0, geometry.tetrahedra().len() / 2] {
            let probe = bound.probe(element, [0.25; 4]).unwrap();
            let original = source.state_at(axis().coordinate_m(probe.position_m()).unwrap()).unwrap();
            assert_eq!(probe.medium(), Medium {
                density_kg_m3: original.state().density,
                sound_speed_m_s: original.state().sound_speed,
            });
            let sample = probe.harmonic(&result.pressure_pa, prepared.frequency_hz()).unwrap();
            let p = bound.operator().pressure_at(&result.pressure_pa, element, [0.25; 4]).unwrap();
            assert_eq!(sample.pressure_pa, p);
        }
        responses.push(result.pressure_pa);
    }
    let norm = responses[0].iter().map(|p| p.norm_sq()).sum::<f64>();
    let difference = responses[0]
        .iter()
        .zip(&responses[1])
        .map(|(&a, &b)| (a - b).norm_sq())
        .sum::<f64>();
    assert!(
        difference > 1e-4 * norm,
        "spatial source must affect the actual pressure solve"
    );
    let middle = sources[1].state_at(0.5).unwrap();
    let ends = [
        sources[1].state_at(0.0).unwrap(),
        sources[1].state_at(1.0).unwrap(),
    ];
    let averaged_density = 0.5 * (ends[0].state().density + ends[1].state().density);
    assert!((middle.state().density / averaged_density - 1.0).abs() > 0.01);
}

#[test]
fn material_jumps_require_conforming_tetrahedra_but_equal_cells_do_not() {
    let source = field(
        vec![0.0, 0.5, 1.0],
        vec![
            // The inlet-only state also differs from the first physical cell.
            primitive(450.0, 0.5),
            primitive(300.0, 1.0),
            primitive(600.0, 0.0),
        ],
        true,
    );
    let (single, single_positions) = single_tet();
    let coarse = HelmholtzFe::new(&single, &single_positions, 1, PLACEHOLDER).unwrap();
    assert!(matches!(
        AxialGasInterior::new(&coarse, &source, axis(), 3, &mut || false),
        Err(AxialGasInteriorError::UnresolvedJump { station_m: 0.5, .. })
    ));
    let (complex, positions) = kuhn_cube(2);
    let conforming = HelmholtzFe::new(&complex, &positions, 1, PLACEHOLDER).unwrap();
    let bound = AxialGasInterior::new(&conforming, &source, axis(), 3, &mut || false).unwrap();
    let coefficients: Vec<_> = positions.iter().map(|point| C64::new(1.0, -point[0])).collect();
    let mut sampled = [0_usize; 3];
    for (element, tet) in conforming.tetrahedra().iter().enumerate() {
        let center_x = tet.iter().map(|&v| 0.25 * positions[v as usize][0]).sum::<f64>();
        let original = source.state_at(center_x).unwrap();
        for (face_kind, station) in [0.0, 0.5, 1.0].into_iter().enumerate() {
            if tet.iter().filter(|&&v| positions[v as usize][0] == station).count() != 3 {
                continue;
            }
            let barycentric = tet.map(|v| if positions[v as usize][0] == station {
                1.0 / 3.0
            } else {
                0.0
            });
            let probe = bound.probe(element, barycentric).unwrap();
            assert_eq!(probe.medium(), Medium {
                density_kg_m3: original.state().density,
                sound_speed_m_s: original.state().sound_speed,
            });
            let sample = probe.harmonic(&coefficients, 1.0).unwrap();
            let expected_v = 1.0 / (core::f64::consts::TAU * original.state().density);
            assert!((sample.particle_velocity_m_s[0].re / expected_v - 1.0).abs() < 1e-12);
            assert!((sample.active_intensity_w_m2[0] / (0.5 * expected_v) - 1.0).abs() < 1e-12);
            if station == 0.0 || (station == 0.5 && center_x > 0.5) {
                let endpoint = source.state_at(station).unwrap();
                assert_ne!(probe.medium().density_kg_m3, endpoint.state().density,
                    "one-sided FE trace must not adopt the source endpoint's other side");
            }
            sampled[face_kind] += 1;
        }
    }
    assert!(sampled.into_iter().all(|count| count > 0));
    let equal = field(vec![0.0, 0.5, 1.0], vec![primitive(300.0, 1.0); 3], true);
    assert!(AxialGasInterior::new(&coarse, &equal, axis(), 3, &mut || false).is_ok());
}

#[test]
fn original_material_failure_inside_the_mesh_is_preserved() {
    let a = gas_support::packs("A", 0.028, 3.5, 200.0);
    let b = gas_support::packs("B", 0.044, 4.5, 500.0);
    let model = from_packs(
        &[a.binding(), b.binding()],
        primitive(300.0, 0.5).composition,
        None,
    )
    .unwrap();
    let source = AxialGasField::new(
        vec![0.0, 1.0],
        100_000.0,
        vec![primitive(300.0, 1.0), primitive(700.0, 0.0)],
        model,
        QueryPoint::new(),
        ThermochemicalPhase::FrozenGas,
    )
    .unwrap();
    let (complex, positions) = single_tet();
    let geometry = HelmholtzFe::new(&complex, &positions, 1, PLACEHOLDER).unwrap();
    assert!(matches!(
        AxialGasInterior::new(&geometry, &source, axis(), 3, &mut || false),
        Err(AxialGasInteriorError::Material(_))
    ));
}

#[test]
fn physical_embedding_source_domain_and_cancellation_are_explicit() {
    let embedding = AxialGasEmbedding::new([2.0, 3.0, 4.0], [0.0, 0.0, -1e300]).unwrap();
    assert_eq!(embedding.coordinate_m([7.0, 9.0, 3.5]).unwrap(), 0.5);
    assert!(AxialGasEmbedding::new([0.0; 3], [0.0; 3]).is_err());
    assert!(axis().coordinate_m([f64::NAN, 0.0, 0.0]).is_err());
    let (complex, positions) = single_tet();
    let geometry = HelmholtzFe::new(&complex, &positions, 1, PLACEHOLDER).unwrap();
    let short = field(vec![0.0, 0.75], vec![primitive(300.0, 0.5); 2], false);
    assert!(matches!(
        AxialGasInterior::new(&geometry, &short, axis(), 3, &mut || false),
        Err(AxialGasInteriorError::OutsideSource { .. })
    ));
    assert!(matches!(
        AxialGasInterior::new(&geometry, &short, axis(), 3, &mut || true),
        Err(AxialGasInteriorError::Cancelled)
    ));
}
