//! Original fluid sources on conforming interfaces, with independent P1 integrals.

use std::cell::Cell;
use std::error::Error;

use fs_couple::cavity_modal::interior::heterogeneous::regions::{
    FluidRegionsInterior, FluidRegionsInteriorError,
};
use fs_feec::acoustics::helmholtz::{HelmholtzFe, Medium};
use fs_feec::fixtures::two_tets;
use fs_material::MaterialError;
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;
use fs_qty::{Pressure, Temperature};
use fs_sparse::Csr;
use fs_thermochem::water::if97::{If97Region, If97State};

// The borrowed geometry's medium is deliberately unrelated to either source.
const PLACEHOLDER: Medium = Medium {
    density_kg_m3: 7.0,
    sound_speed_m_s: 123.0,
};

fn entry(matrix: &Csr, row: usize, column: usize) -> f64 {
    let (columns, values) = matrix.row(row);
    columns
        .iter()
        .position(|&value| value == column)
        .map_or(0.0, |i| values[i])
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * expected.abs().max(f64::MIN_POSITIVE),
        "{actual} != {expected}"
    );
}

#[test]
fn gas_water_interface_assembles_original_physical_states_and_shared_pressure_dofs() {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100_000.0).unwrap();
    let water = If97State::at_pressure_temperature(
        Pressure::new(100_000.0),
        Temperature::new(300.0),
        If97Region::Region1,
    )
    .unwrap();
    let sources: [&dyn FluidAcoustics; 2] = [&gas, &water];
    let (complex, positions) = two_tets();
    let geometry = HelmholtzFe::new(&complex, &positions, 1, PLACEHOLDER).unwrap();
    let bound = FluidRegionsInterior::new(&geometry, &sources, &[0, 1], 3, &mut || false).unwrap();
    assert!(std::ptr::eq(bound.sources()[0], sources[0]));
    assert!(std::ptr::eq(bound.sources()[1], sources[1]));
    assert!(std::ptr::eq(bound.operator().geometry(), &geometry));
    assert_eq!(bound.element_regions(), &[0, 1]);
    let a = gas.acoustic_state().unwrap();
    let b = water.acoustic_state().unwrap();
    assert_eq!(bound.region_state(0), Some(a));
    assert_eq!(bound.region_state(1), Some(b));
    assert_eq!(bound.region_state(2), None);

    // The same geometrical point on the shared face has two explicit material
    // traces. Velocity uses its original density, never the geometry placeholder
    // or the other side's gas/liquid state.
    let coefficients: Vec<_> = positions.iter().map(|point| {
        C64::new(1.0, 1.0 - point.iter().sum::<f64>())
    }).collect();
    for (element, original) in [a, b].into_iter().enumerate() {
        let barycentric = geometry.tetrahedra()[element].map(|v| {
            if (1..=3).contains(&v) { 1.0 / 3.0 } else { 0.0 }
        });
        let probe = bound.probe(element, barycentric).unwrap();
        assert_eq!(probe.medium(), Medium {
            density_kg_m3: original.density_kg_m3(),
            sound_speed_m_s: original.sound_speed_m_s(),
        });
        let sample = probe.harmonic(&coefficients, 1.0).unwrap();
        let expected_v = 1.0 / (core::f64::consts::TAU * original.density_kg_m3());
        near(sample.pressure_pa.re, 1.0);
        for axis in 0..3 {
            near(sample.particle_velocity_m_s[axis].re, expected_v);
            near(sample.active_intensity_w_m2[axis], 0.5 * expected_v);
        }
    }

    // Volumes 1/6 and 1/3; the shared face is x+y+z=1. Vertex 0
    // belongs only to gas, vertex 4 only to water, vertices 1..3 to both.
    // Exact P1 integrals: int(phi_i²)=V/10, int(phi_i phi_j)=V/20.
    let k = bound.operator().stiffness();
    let m = bound.operator().mass();
    let inv_a = 1.0 / a.density_kg_m3();
    let inv_b = 1.0 / b.density_kg_m3();
    let beta_a = a.adiabatic_compressibility_per_pa();
    let beta_b = b.adiabatic_compressibility_per_pa();
    near(entry(k, 0, 0), 0.5 * inv_a);
    near(entry(k, 4, 4), 0.25 * inv_b);
    near(entry(k, 1, 1), inv_a / 6.0 + inv_b / 4.0);
    near(entry(k, 1, 4), -inv_b / 12.0);
    near(entry(m, 0, 0), beta_a / 60.0);
    near(entry(m, 4, 4), beta_b / 30.0);
    near(entry(m, 1, 1), beta_a / 60.0 + beta_b / 30.0);
    near(entry(m, 1, 4), beta_b / 60.0);
    assert_eq!(entry(k, 0, 4), 0.0);
    assert_eq!(entry(m, 0, 4), 0.0);
    for row in 0..5 {
        let sum: f64 = (0..5).map(|col| entry(k, row, col)).sum();
        assert!(sum.abs() < 1e-12 * entry(k, row, row));
    }
}

struct CountedSource {
    state: Cell<FluidAcousticState>,
    queries: Cell<usize>,
    refuse: bool,
}

impl FluidAcoustics for CountedSource {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        self.queries.set(self.queries.get() + 1);
        if self.refuse {
            return Err(MaterialError::Parameters {
                what: "explicit test source refusal".into(),
            });
        }
        Ok(self.state.get())
    }
}

fn counted(refuse: bool) -> CountedSource {
    CountedSource {
        state: Cell::new(FluidAcousticState::try_new(1.25, 340.0).unwrap()),
        queries: Cell::new(0),
        refuse,
    }
}

#[test]
fn one_query_per_used_source_freezes_coherent_coefficients_across_high_order_elements() {
    let used = counted(false);
    let unused = counted(true);
    let sources: [&dyn FluidAcoustics; 2] = [&used, &unused];
    let (complex, positions) = two_tets();
    let geometry = HelmholtzFe::new(&complex, &positions, 3, PLACEHOLDER).unwrap();
    let bound = FluidRegionsInterior::new(&geometry, &sources, &[0, 0], 5, &mut || false).unwrap();
    let frozen = used.state.get();
    used.state
        .set(FluidAcousticState::try_new(900.0, 1200.0).unwrap());
    assert_eq!(bound.region_state(0), Some(frozen));
    assert_eq!(bound.region_state(1), None);
    for (physical, geometric, scale) in [
        (
            bound.operator().stiffness(),
            geometry.stiffness(),
            1.0 / frozen.density_kg_m3(),
        ),
        (
            bound.operator().mass(),
            geometry.mass(),
            frozen.adiabatic_compressibility_per_pa(),
        ),
    ] {
        let norm = (0..geometry.ndof())
            .flat_map(|row| geometric.row(row).1)
            .map(|value| value.abs() * scale)
            .fold(0.0_f64, f64::max);
        for row in 0..geometry.ndof() {
            for col in 0..geometry.ndof() {
                assert!(
                    (entry(physical, row, col) - scale * entry(geometric, row, col)).abs()
                        < 1e-11 * norm
                );
            }
        }
    }
    let mut coefficients = vec![C64::ZERO; geometry.ndof()];
    for (vertex, point) in positions.iter().enumerate() {
        coefficients[vertex] = C64::new(1.0, -point[0]);
    }
    for element in 0..geometry.tetrahedra().len() {
        let probe = bound.probe(element, [0.25; 4]).unwrap();
        assert_eq!(probe.medium(), Medium {
            density_kg_m3: frozen.density_kg_m3(),
            sound_speed_m_s: frozen.sound_speed_m_s(),
        });
        let sample = probe.harmonic(&coefficients, 1.0).unwrap();
        near(sample.particle_velocity_m_s[0].re,
            1.0 / (core::f64::consts::TAU * frozen.density_kg_m3()));
    }
    assert!(matches!(bound.probe(2, [0.25; 4]), Err(FluidRegionsInteriorError::Interior(_))));
    assert!(matches!(bound.probe(0, [f64::NAN; 4]), Err(FluidRegionsInteriorError::Interior(_))));
    assert_eq!(used.queries.get(), 1);
    assert_eq!(unused.queries.get(), 0);
}

#[test]
fn invalid_regions_and_quadrature_refuse_before_query_and_original_material_error_survives() {
    let source = counted(false);
    let (complex, positions) = two_tets();
    let geometry = HelmholtzFe::new(&complex, &positions, 1, PLACEHOLDER).unwrap();
    assert!(matches!(
        FluidRegionsInterior::new(&geometry, &[&source], &[0], 3, &mut || false),
        Err(FluidRegionsInteriorError::ElementCount {
            expected: 2,
            supplied: 1
        })
    ));
    assert!(matches!(
        FluidRegionsInterior::new(&geometry, &[&source], &[0, 7], 3, &mut || false),
        Err(FluidRegionsInteriorError::UnknownRegion {
            tetrahedron: 1,
            region: 7,
            source_count: 1
        })
    ));
    assert!(matches!(
        FluidRegionsInterior::new(&geometry, &[&source], &[0, 0], 2, &mut || false),
        Err(FluidRegionsInteriorError::Interior(_))
    ));
    assert_eq!(source.queries.get(), 0);
    let mut gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100_000.0).unwrap();
    gas.sound_speed *= 1.01;
    let error = FluidRegionsInterior::new(&geometry, &[&source, &gas], &[1, 0], 3, &mut || false)
        .err()
        .expect("incoherent original gas must refuse");
    assert!(matches!(
        error,
        FluidRegionsInteriorError::Material { region: 1, .. }
    ));
    assert!(error.source().is_some());
    assert_eq!(source.queries.get(), 0);
}

#[test]
fn cancellation_before_and_during_assembly_never_returns_a_partial_binding() {
    let source = counted(false);
    let (complex, positions) = two_tets();
    let geometry = HelmholtzFe::new(&complex, &positions, 1, PLACEHOLDER).unwrap();
    assert!(matches!(
        FluidRegionsInterior::new(&geometry, &[&source], &[0, 0], 3, &mut || true),
        Err(FluidRegionsInteriorError::Cancelled)
    ));
    assert_eq!(source.queries.get(), 0);
    assert!(matches!(
        FluidRegionsInterior::new(&geometry, &[&source], &[0, 0], 3, &mut || source
            .queries
            .get()
            > 0),
        Err(FluidRegionsInteriorError::Cancelled)
    ));
    assert_eq!(source.queries.get(), 1);
}
