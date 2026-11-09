use super::*;
use crate::curle3d::quadrupole::flux::tests::molecular::{
    GRADIENT, affine_flow, affine_velocity, law, stress,
};
use fs_flux::compressible::spatial::mixture::molecular as multi;

const CONTROL_BOX: PermeableBox3 = PermeableBox3 { lower: [2; 3], upper: [3; 3] };

fn check_stress_loading(
    inviscid: &PreparedPermeableSurface3,
    full: &PreparedPermeableSurface3,
    tensor: [[f64; 3]; 3],
) {
    assert_eq!(full.faces(), inviscid.faces());
    assert_eq!(full.source_positions_m(), inviscid.source_positions_m());
    let mut total = [0.0; 3];
    for ((face, original), current) in inviscid.faces().iter()
        .zip(inviscid.source_components()).zip(full.source_components())
    {
        // Viscous force has no mass source. Native conserved momentum remains
        // the exact input to both original surface Q channels.
        assert_eq!(current[0], original[0]);
        for (component, row) in full.config().grid_to_acoustic.iter().enumerate() {
            let tau: f64 = (0..3).map(|a| row[a] * tensor[a][face.side.axis()]).sum();
            let expected = -tau * face.side.sign() * face.area_m2;
            close(current[component + 1] - original[component + 1], expected, 2e-11);
            total[component] += current[component + 1] - original[component + 1];
        }
    }
    // Constant stress has zero divergence: equal opposite faces cancel,
    // while retaining all individually nonzero signed shear/normal sources.
    for resultant in total {
        close(resultant, 0.0, 2e-11);
    }
}

#[test]
fn g1_native_viscous_surface_loading_has_full_stress_sign_area_and_rotation() {
    let (mu, bulk) = (0.12, 0.07);
    let mut flow = affine_flow(mu, bulk, GRADIENT);
    let settings = config(flow.flow().gas());
    let count = CONTROL_BOX.face_count(flow.flow().grid()).unwrap();
    let before = flow.flow().states().to_vec();
    let mut inviscid = PreparedPermeableSurface3::from_frozen(
        flow.flow(), CONTROL_BOX, settings, controls(count, 0), || false,
    ).unwrap();
    let mut full = PreparedPermeableSurface3::from_molecular_frozen(
        &mut flow, CONTROL_BOX, settings, controls(count, 0), || false,
    ).unwrap();
    assert!(full.includes_viscous_stress());
    assert!(!inviscid.includes_viscous_stress());
    check_stress_loading(&inviscid, &full, stress(mu, bulk, GRADIENT));
    assert_eq!(flow.flow().states(), before);
    assert_eq!(full.sample_time_s(), flow.flow().time().value());
    assert_eq!(
        full.sample_frozen(flow.flow(), || false),
        Err(CompressibleVolumeError3::SourceMismatch),
    );
    assert_eq!(
        inviscid.sample_molecular_frozen(&mut flow, || false),
        Err(CompressibleVolumeError3::SourceMismatch),
    );
}

#[test]
fn g0_zero_viscosity_preserves_surface_sources_and_dd_keeps_residual_force() {
    let mut flow = affine_flow(0.0, 0.0, GRADIENT);
    let settings = config(flow.flow().gas());
    let count = CONTROL_BOX.face_count(flow.flow().grid()).unwrap();
    let inviscid = PreparedPermeableSurface3::from_frozen(
        flow.flow(), CONTROL_BOX, settings, controls(count, 0), || false,
    ).unwrap();
    let full = PreparedPermeableSurface3::from_molecular_frozen(
        &mut flow, CONTROL_BOX, settings, controls(count, 0), || false,
    ).unwrap();
    assert_eq!(full.source_components(), inviscid.source_components());

    // An O(1) pressure residual survives cancellation of O(1e16) advective
    // and viscous loading. Subtracting tau AFTER a rounded inviscid force
    // would incorrectly return zero.
    let mut settings = settings;
    settings.grid_to_acoustic = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let face = PermeableFace3 {
        adjacent_cells: [0, 1],
        side: BoundarySide::XMax,
        area_m2: 1.0,
    };
    let mut primitive = Primitive::default();
    primitive.density = 1.0;
    primitive.velocity = [1e8, 0.0, 0.0];
    primitive.pressure = settings.reference_pressure_pa + 1.0;
    let states = [(primitive, [1e8, 0.0, 0.0]); 2];
    let original = surface_flux(&face, settings, states, None).unwrap();
    assert_eq!(
        surface_flux(&face, settings, states, Some([0.0; 3])).unwrap(),
        original,
    );
    let resolved = surface_flux(&face, settings, states, Some([1e16, 0.0, 0.0])).unwrap();
    assert_eq!(resolved[0], original[0]);
    assert_eq!(&resolved[1..], &[1.0, 0.0, 0.0]);
}

#[test]
fn g1_mixture_surface_retains_original_composition_and_viscous_traction() {
    let mechanism = cantera_h2o2(HydrogenCalorics::Original).unwrap();
    let gas = MixtureGas::new(
        mechanism.thermodynamics().clone(),
        CaloricDomain::new(Temperature::new(300.0), Temperature::new(900.0)),
        || false,
    ).unwrap();
    let settings = config(&frozen(6));
    let grid = CartesianGrid::new(5, 5, 5, 0.2, 0.3, 0.4).unwrap();
    let species = gas.species_count();
    let mut ambient = vec![0.0; species];
    ambient[6] = settings.acoustic.density_kg_m3();
    let (mu, bulk) = (0.12, 0.07);
    let transport = law(gas.source().clone(), mu, bulk);
    let mut states = vec![0.0; 125 * (species + 4)];
    for (cell, row) in states.chunks_exact_mut(species + 4).enumerate() {
        let x = grid.cell_center(cell).unwrap();
        let mut partial = vec![0.0; species];
        partial[0] = 0.15 + 0.02 * x[0];
        partial[6] = 1.7 - partial[0];
        state_into(&gas, &partial, affine_velocity(x, GRADIENT), settings.reference_pressure_pa, row)
            .unwrap();
    }
    let bytes = MixtureEuler3d::required_workspace_bytes(&gas, &grid).unwrap();
    let native = MixtureEuler3d::new(gas, grid, states, PERIODIC, bytes, || false).unwrap();
    let source = multi::SpatialTransport::new(native.gas().clone(), transport, || false).unwrap();
    let bytes = multi::MolecularTransport3d::required_workspace_bytes(&native, &source).unwrap();
    let operator = multi::MolecularTransport3d::new(&native, source, bytes, || false).unwrap();
    let bytes = multi::MolecularEuler3d::required_workspace_bytes(&native, &operator).unwrap();
    let mut flow = multi::MolecularEuler3d::new(native, operator, bytes, || false).unwrap();
    let before = flow.flow().states().to_vec();
    let count = CONTROL_BOX.face_count(&grid).unwrap();
    let inviscid = PreparedPermeableSurface3::from_mixture(
        flow.flow(), &ambient, CONTROL_BOX, settings, controls(count, species), || false,
    ).unwrap();
    let mut full = PreparedPermeableSurface3::from_molecular_mixture(
        &mut flow, &ambient, CONTROL_BOX, settings, controls(count, species),
        multi::DissipationControl::default(), || false,
    ).unwrap();
    assert_eq!(full.reference_partial_densities_kg_m3(), ambient);
    check_stress_loading(&inviscid, &full, stress(mu, bulk, GRADIENT));
    assert_eq!(flow.flow().states(), before);
    assert_eq!(
        full.sample_mixture(flow.flow(), || false),
        Err(CompressibleVolumeError3::SourceMismatch),
    );
    flow.step(Time::new(1e-7), multi::MolecularStepControl::default(), || false).unwrap();
    let accepted = flow.flow().states().to_vec();
    full.sample_molecular_mixture(&mut flow, multi::DissipationControl::default(), || false).unwrap();
    assert_eq!(full.sample_time_s(), flow.flow().time().value());
    assert_eq!(flow.flow().states(), accepted);
}
