//! G1/G3: original JCA(L) material responses in full 3-D pressure volumes.

use fs_couple::cavity_modal::interior::heterogeneous::porous::{
    AcousticVolumeRegion, AcousticVolumeRegionState, PorousRegionsInterior,
    PorousRegionsInteriorError,
};
use fs_feec::acoustics::helmholtz::heterogeneous::{
    EquivalentFluidControls, HeterogeneousError, VelocityTrace,
};
use fs_feec::acoustics::helmholtz::{HelmholtzFe, Medium};
use fs_feec::fixtures::{masked_cube_grid, single_tet};
use fs_material::MaterialError;
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_material::fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState};
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::gas::{GasSpec, GasState};
use fs_material::porous::{PorousCard, PorousFluidLaw};
use fs_math::c64::C64;
use fs_rep_mesh::TetComplex;
use std::cell::Cell;

#[path = "porous_interior/molecular.rs"]
mod molecular;

fn gas(temperature_k: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature_k, 101325.0).unwrap()
}

fn law() -> PorousFluidLaw {
    PorousFluidLaw::try_new(PorousCard::melamine_uf_nguyen2024(), 1.81e-5).unwrap()
}

fn medium(gas: &GasState) -> Medium {
    Medium { density_kg_m3: gas.density, sound_speed_m_s: gas.sound_speed }
}

fn brick(nx: usize, length: f64) -> (TetComplex, Vec<[f64; 3]>) {
    let (complex, mut positions) = masked_cube_grid(nx, 1, 1, &|_, _, _| true);
    for point in &mut positions {
        point[0] *= length / nx as f64;
        point[1] *= 0.012;
        point[2] *= 0.015;
    }
    (complex, positions)
}

fn inlet(geometry: &HelmholtzFe<'_>) -> Vec<[u32; 3]> {
    geometry.boundary_faces().into_iter().filter(|face| face.iter()
        .all(|&v| geometry.positions()[v as usize][0] == 0.0)).collect()
}

fn drive(faces: &[[u32; 3]]) -> Vec<VelocityTrace> {
    faces.iter().map(|&face| (face, [C64::ONE; 3])).collect()
}

fn mean_pressure(geometry: &HelmholtzFe<'_>, faces: &[[u32; 3]], field: &[C64]) -> C64 {
    let (load, area) = geometry.boundary_load(faces).unwrap();
    load.iter().zip(field).fold(C64::ZERO, |sum, (&b, &p)| sum + p.scale(b / area))
}

#[test]
fn rigid_backed_porous_volume_converges_to_original_jcal_layer_impedance() {
    let gas = gas(293.15);
    let law = law();
    let frequency = 1000.0;
    let omega = core::f64::consts::TAU * frequency;
    let length = 0.05;
    let expected = law.surface_impedance_rigid_fluid(&gas, omega, length).unwrap();
    let source = [AcousticVolumeRegion::Porous { law: &law, fluid: &gas }];
    let mut errors = Vec::new();
    for nx in [2, 4, 8] {
        let (complex, positions) = brick(nx, length);
        let geometry = HelmholtzFe::new(&complex, &positions, 2, medium(&gas)).unwrap();
        let model = PorousRegionsInterior::new(
            &geometry,
            frequency,
            &source,
            &vec![0; complex.tets.len()],
            &[],
            &[],
            EquivalentFluidControls::default(),
            &mut || false,
        ).unwrap();
        let faces = inlet(&geometry);
        let response = model.operator().solve(&drive(&faces), &mut || false).unwrap();
        let actual = mean_pressure(&geometry, &faces, &response.pressure_pa);
        errors.push((actual - expected).abs() / expected.abs());
        let snapshot = match model.region_state(0).unwrap() {
            AcousticVolumeRegionState::Porous(response) => response,
            _ => panic!("porous material was replaced"),
        };
        assert_eq!(snapshot, law.effective_fluid(&gas, omega).unwrap());
        assert!(response.density_dissipation_w > 0.0 && response.compressional_dissipation_w > 0.0);
        assert_eq!(response.boundary_dissipation_w, 0.0);
        assert!(response.power_balance_residual_w().abs() < 2e-10 * response.input_power_w);
    }
    assert!(errors[1] < errors[0] / 4.0 && errors[2] < errors[1] / 4.0, "{errors:?}");
    assert!(errors[2] < 5e-4, "{errors:?}");
}

#[test]
fn conforming_free_fluid_and_porous_layer_use_the_same_pressure_and_area_convention() {
    let gas = gas(310.15);
    let law = law();
    let frequency = 800.0;
    let omega = core::f64::consts::TAU * frequency;
    let (complex, positions) = brick(12, 0.1);
    let geometry = HelmholtzFe::new(&complex, &positions, 2, medium(&gas)).unwrap();
    let regions: Vec<_> = complex.tets.iter().map(|tet| {
        let x: f64 = tet.iter().map(|&v| positions[v as usize][0]).sum::<f64>() / 4.0;
        usize::from(x > 0.05)
    }).collect();
    let sources = [
        AcousticVolumeRegion::Fluid(&gas),
        AcousticVolumeRegion::Porous { law: &law, fluid: &gas },
    ];
    let model = PorousRegionsInterior::new(
        &geometry,
        frequency,
        &sources,
        &regions,
        &[],
        &[],
        EquivalentFluidControls::default(),
        &mut || false,
    ).unwrap();
    let faces = inlet(&geometry);
    let response = model.operator().solve(&drive(&faces), &mut || false).unwrap();
    let z_layer = law.surface_impedance_rigid_fluid(&gas, omega, 0.05).unwrap();
    let z_gas = gas.density * gas.sound_speed;
    let (s, c) = (omega * 0.05 / gas.sound_speed).sin_cos();
    let expected = (z_layer.scale(c) + C64::new(0.0, z_gas * s))
        / (C64::from_re(c) + C64::new(0.0, s / z_gas) * z_layer);
    let actual = mean_pressure(&geometry, &faces, &response.pressure_pa);
    assert!((actual - expected).abs() < 4e-4 * expected.abs(), "{actual:?} versus {expected:?}");
    for (&region, power) in regions.iter().zip(&response.element_power) {
        if region == 0 {
            assert_eq!(power.total_dissipation_w(), 0.0);
        } else {
            assert!(power.total_dissipation_w() > 0.0);
        }
    }
    assert!(response.power_balance_residual_w().abs() < 3e-10 * response.input_power_w);
}

struct CountedGas<'a> {
    gas: &'a GasState,
    queries: Cell<usize>,
}

impl FluidAcoustics for CountedGas<'_> {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        self.gas.acoustic_state()
    }
}

impl FluidTransport for CountedGas<'_> {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        self.gas.transport_state()
    }
}

impl ThermoviscousFluid for CountedGas<'_> {
    fn thermoviscous_state(&self) -> Result<ThermoviscousFluidState, MaterialError> {
        self.queries.set(self.queries.get() + 1);
        self.gas.thermoviscous_state()
    }
}

#[test]
fn sources_are_sampled_once_after_admission_and_keep_native_refusals() {
    let gas = gas(293.15);
    let counted = CountedGas { gas: &gas, queries: Cell::new(0) };
    let law = law();
    let (complex, mut positions) = single_tet();
    for p in &mut positions { *p = p.map(|x| x * 0.01); }
    let geometry = HelmholtzFe::new(&complex, &positions, 2, medium(&gas)).unwrap();
    let sources = [AcousticVolumeRegion::Porous { law: &law, fluid: &counted }];
    assert!(matches!(
        PorousRegionsInterior::new(
            &geometry,
            1000.0,
            &sources,
            &[1],
            &[],
            &[],
            EquivalentFluidControls::default(),
            &mut || false,
        ),
        Err(PorousRegionsInteriorError::UnknownRegion { .. }),
    ));
    assert_eq!(counted.queries.get(), 0);
    assert!(matches!(
        PorousRegionsInterior::new(
            &geometry,
            1000.0,
            &sources,
            &[0],
            &[],
            &[],
            EquivalentFluidControls { max_volume_points: 0, ..EquivalentFluidControls::default() },
            &mut || false,
        ),
        Err(PorousRegionsInteriorError::Interior(HeterogeneousError::Limit { .. })),
    ));
    assert_eq!(counted.queries.get(), 0);
    let model = PorousRegionsInterior::new(
        &geometry,
        1000.0,
        &sources,
        &[0],
        &[],
        &[],
        EquivalentFluidControls::default(),
        &mut || false,
    ).unwrap();
    assert_eq!(counted.queries.get(), 1);
    let faces = geometry.boundary_faces();
    let first = model.operator().solve(&drive(&faces[..1]), &mut || false).unwrap();
    assert_eq!(
        model.operator().solve(&drive(&faces[..1]), &mut || false).unwrap().pressure_pa,
        first.pressure_pa,
    );
    assert_eq!(counted.queries.get(), 1);
    let unsupported = law.reference_card().valid_hz.0 / 2.0;
    assert!(matches!(
        PorousRegionsInterior::new(
            &geometry,
            unsupported,
            &sources,
            &[0],
            &[],
            &[],
            EquivalentFluidControls::default(),
            &mut || false,
        ),
        Err(PorousRegionsInteriorError::Material { region: 0, .. }),
    ));
    assert_eq!(
        counted.queries.get(),
        1,
        "unsupported frequency must not query a replacement fluid",
    );
}
