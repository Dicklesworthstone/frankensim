use super::*;
use fs_flux::compressible::spatial::{mixture::molecular as multi, molecular as single};
use fs_gas::multicomponent::dissipation::{
    ConstantTransport, TransportCoefficients, TransportDomain,
};
use fs_qty::DynViscosity;
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::transport::kinetic::ThermalConductivityQuantityV1;
use fs_thermochem::{ChargeVector, ElementId, ElementalMatrix};

pub(in crate::curle3d::quadrupole::flux) const GRADIENT: [[f64; 3]; 3] = [[2.0, 3.0, -4.0], [-1.0, -0.5, 5.0], [6.0, -2.0, 1.5]];
const ROTATION: [[f64; 3]; 3] =
    [[0.48, -0.6, 0.64], [0.36, 0.8, 0.48], [-0.8, 0.0, 0.6]];
const CENTER: usize = 2 + 5 * (2 + 5 * 2);

pub(in crate::curle3d::quadrupole::flux) fn law(source: IdealGasEquilibriumV1, mu: f64, bulk: f64) -> ConstantTransport {
    let n = source.models().len();
    let mut binary = vec![1e-4; n * n];
    for i in 0..n {
        binary[i * n + i] = 0.0;
    }
    ConstantTransport::new(
        source,
        "manufactured-acoustic-stress",
        1,
        TransportDomain {
            minimum_temperature_k: 300.0,
            maximum_temperature_k: 900.0,
            minimum_pressure_pa: 1e3,
            maximum_pressure_pa: 1e7,
        },
        TransportCoefficients {
            dynamic_viscosity: DynViscosity::new(mu),
            bulk_viscosity: DynViscosity::new(bulk),
            thermal_conductivity: ThermalConductivityQuantityV1::new(0.0),
        },
        binary,
    )
    .unwrap()
}

pub(in crate::curle3d::quadrupole::flux) fn frozen_law(gas: &FrozenGas, mu: f64, bulk: f64) -> ConstantTransport {
    let models = gas.source().components().to_vec();
    let ids: Vec<_> = models.iter().map(|m| m.species().clone()).collect();
    // The test's sole species is the original pinned N2 caloric source.
    let source = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("N").unwrap()],
            ids.clone(),
            vec![vec![2]],
        )
        .unwrap(),
        ChargeVector::new(ids, vec![0]).unwrap(),
        models,
    )
    .unwrap();
    law(source, mu, bulk)
}

pub(in crate::curle3d::quadrupole::flux) fn affine_velocity(x: [f64; 3], gradient: [[f64; 3]; 3]) -> [f64; 3] {
    let center = [0.5, 0.75, 1.0];
    let base = [0.1, -0.2, 0.3];
    std::array::from_fn(|i| {
        base[i] + (0..3).map(|a| gradient[i][a] * (x[a] - center[a])).sum::<f64>()
    })
}

pub(in crate::curle3d::quadrupole::flux) fn flow_with_velocity(
    grid: CartesianGrid,
    mu: f64,
    bulk: f64,
    velocity: impl Fn([f64; 3]) -> [f64; 3],
) -> single::MolecularEuler3d<ConstantTransport> {
    let gas = frozen(6);
    let transport = frozen_law(&gas, mu, bulk);
    let count = grid.shape().iter().product();
    let states = (0..count)
        .map(|cell| {
            State::from_primitive(
                &gas,
                1.0,
                velocity(grid.cell_center(cell).unwrap()),
                gas.gas_constant() * 500.0,
            )
            .unwrap()
        })
        .collect();
    let flow = Euler3d::new(gas, grid, states, PERIODIC, || false).unwrap();
    single::MolecularEuler3d::new(flow, transport, || false).unwrap()
}

pub(in crate::curle3d::quadrupole::flux) fn affine_flow(
    mu: f64,
    bulk: f64,
    gradient: [[f64; 3]; 3],
) -> single::MolecularEuler3d<ConstantTransport> {
    flow_with_velocity(
        CartesianGrid::new(5, 5, 5, 0.2, 0.3, 0.4).unwrap(),
        mu,
        bulk,
        |x| affine_velocity(x, gradient),
    )
}

// Independent continuum Newtonian law, including all transverse derivatives
// and authored bulk viscosity. Manufactured affine fields make it exact.
pub(in crate::curle3d::quadrupole::flux) fn stress(mu: f64, bulk: f64, gradient: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let dilation: f64 = (0..3).map(|i| gradient[i][i]).sum();
    std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            mu * (gradient[i][j] + gradient[j][i])
                + if i == j { (bulk - 2.0 * mu / 3.0) * dilation } else { 0.0 }
        })
    })
}

fn rotated(stress: [[f64; 3]; 3]) -> [f64; 6] {
    [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)].map(|(i, j)| {
        (0..3)
            .flat_map(|a| (0..3).map(move |b| ROTATION[i][a] * stress[a][b] * ROTATION[j][b]))
            .sum()
    })
}

#[test]
fn g1_original_full_stress_has_physical_face_sign_bulk_and_acoustic_rotation() {
    let (mu, bulk) = (0.12, 0.07);
    let mut flow = affine_flow(mu, bulk, GRADIENT);
    let original = flow.flow().states().to_vec();
    let mut settings = config(flow.flow().gas());
    settings.grid_to_acoustic = ROTATION;
    let volume = flow.flow().grid().cell_volume();
    let expected = stress(mu, bulk, GRADIENT);
    {
        let field = flow.stress_field(|| false).unwrap();
        assert_eq!(field.time_s(), 0.0);
        assert_eq!(field.momentum(CENTER), Some(original[CENTER].momentum));
        assert_eq!(field.momentum(original.len()), None);
        assert_eq!(field.face_traction_pa(3, [2; 3]), None);
        assert_eq!(field.face_traction_pa(0, [6, 2, 2]), None);
        for axis in 0..3 {
            let traction = field.face_traction_pa(axis, [2; 3]).unwrap();
            for component in 0..3 {
                close(traction[component], expected[component][axis]);
            }
        }
        let observed = field.cell_stress_pa(CENTER).unwrap();
        for (k, (i, j)) in [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
            .into_iter()
            .enumerate()
        {
            close(observed[k], expected[i][j]);
        }
    }
    let bytes = PreparedCompressibleVolume3::required_workspace_bytes(1).unwrap();
    let inviscid =
        PreparedCompressibleVolume3::from_frozen(flow.flow(), &[CENTER], settings, bytes, || false)
            .unwrap();
    let full = PreparedCompressibleVolume3::from_molecular_frozen(
        &mut flow, &[CENTER], settings, bytes, || false,
    )
    .unwrap();
    assert!(full.includes_viscous_stress());
    assert!(!inviscid.includes_viscous_stress());
    for (k, tau) in rotated(expected).into_iter().enumerate() {
        close(
            full.stress_moments_n_m()[0][k] - inviscid.stress_moments_n_m()[0][k],
            -tau * volume,
        );
    }
    assert_eq!(flow.flow().states(), original);
    assert_eq!(flow.flow().time().value(), 0.0);
    assert!(flow.last_step().is_none());
}

#[test]
fn g0_rigid_motion_dilatation_and_zero_transport_have_their_physical_limits() {
    let rotation = [[0.0, -2.0, 3.0], [2.0, 0.0, -4.0], [-3.0, 4.0, 0.0]];
    let dilation = [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]];
    for (mu, bulk, gradient, expected) in [
        (0.12, 0.07, [[0.0; 3]; 3], [0.0; 6]),
        (0.12, 0.07, rotation, [0.0; 6]),
        (0.12, 0.07, dilation, [0.42, 0.42, 0.42, 0.0, 0.0, 0.0]),
        (0.0, 0.07, GRADIENT, [0.21, 0.21, 0.21, 0.0, 0.0, 0.0]),
        (0.0, 0.0, GRADIENT, [0.0; 6]),
    ] {
        let mut flow = affine_flow(mu, bulk, gradient);
        let field = flow.stress_field(|| false).unwrap();
        for (actual, target) in field.cell_stress_pa(CENTER).unwrap().into_iter().zip(expected) {
            close(actual, target);
        }
    }
    let mut flow = affine_flow(0.0, 0.0, GRADIENT);
    let settings = config(flow.flow().gas());
    let bytes = PreparedCompressibleVolume3::required_workspace_bytes(1).unwrap();
    let inviscid =
        PreparedCompressibleVolume3::from_frozen(flow.flow(), &[CENTER], settings, bytes, || false)
            .unwrap();
    let full = PreparedCompressibleVolume3::from_molecular_frozen(
        &mut flow, &[CENTER], settings, bytes, || false,
    )
    .unwrap();
    assert_eq!(full.stress_moments_n_m(), inviscid.stress_moments_n_m());
}

#[test]
fn g1_smooth_periodic_shear_stress_refines_at_second_order() {
    let mu = 0.12;
    let wave = core::f64::consts::TAU;
    let mut errors = Vec::new();
    for n in [8, 16, 32] {
        let grid = CartesianGrid::new(n, 2, 2, 1.0 / n as f64, 0.25, 0.25).unwrap();
        let mut flow =
            flow_with_velocity(grid, mu, 0.07, |x| [0.0, (wave * x[0]).sin(), 0.0]);
        let field = flow.stress_field(|| false).unwrap();
        let mut error = 0.0;
        for cell in 0..n {
            let x = grid.cell_center(cell).unwrap()[0];
            let tensor = field.cell_stress_pa(cell).unwrap();
            let expected = mu * wave * (wave * x).cos();
            error += (tensor[3] - expected).powi(2);
            for k in [0, 1, 2, 4, 5] {
                close(tensor[k], 0.0);
            }
        }
        errors.push((error / n as f64).sqrt());
    }
    for pair in errors.windows(2) {
        let order = (pair[0] / pair[1]).log2();
        assert!((1.9..2.1).contains(&order), "{errors:?}, order={order}");
    }
}

#[test]
fn g1_multicomponent_native_stress_keeps_species_owner_and_global_face_axes() {
    let mechanism = cantera_h2o2(HydrogenCalorics::Original).unwrap();
    let gas = MixtureGas::new(
        mechanism.thermodynamics().clone(),
        CaloricDomain::new(Temperature::new(300.0), Temperature::new(900.0)),
        || false,
    )
    .unwrap();
    let (mu, bulk) = (0.12, 0.07);
    let transport = law(gas.source().clone(), mu, bulk);
    let grid = CartesianGrid::new(5, 5, 5, 0.2, 0.3, 0.4).unwrap();
    let species = gas.species_count();
    let mut states = vec![0.0; 125 * (species + 4)];
    for (cell, row) in states.chunks_exact_mut(species + 4).enumerate() {
        let x = grid.cell_center(cell).unwrap();
        let mut densities = vec![0.0; species];
        densities[0] = 0.1 + 0.1 * x[0];
        densities[6] = 0.7 - densities[0];
        state_into(&gas, &densities, affine_velocity(x, GRADIENT), 100000.0, row).unwrap();
    }
    let bytes = MixtureEuler3d::required_workspace_bytes(&gas, &grid).unwrap();
    let native = MixtureEuler3d::new(gas, grid, states, PERIODIC, bytes, || false).unwrap();
    let source = multi::SpatialTransport::new(native.gas().clone(), transport, || false).unwrap();
    let bytes = multi::MolecularTransport3d::required_workspace_bytes(&native, &source).unwrap();
    let operator = multi::MolecularTransport3d::new(&native, source, bytes, || false).unwrap();
    let bytes = multi::MolecularEuler3d::required_workspace_bytes(&native, &operator).unwrap();
    let mut flow = multi::MolecularEuler3d::new(native, operator, bytes, || false).unwrap();
    let original = flow.flow().states().to_vec();
    let mut settings = config(&frozen(6));
    settings.grid_to_acoustic = ROTATION;
    let expected = stress(mu, bulk, GRADIENT);
    {
        let field = flow.stress_field(multi::DissipationControl::default(), || false).unwrap();
        let start = CENTER * (species + 4) + species;
        assert_eq!(
            field.momentum(CENTER).unwrap(),
            [original[start], original[start + 1], original[start + 2]],
        );
        for axis in 0..3 {
            let traction = field.face_traction_pa(axis, [2; 3]).unwrap();
            for component in 0..3 {
                close(traction[component], expected[component][axis]);
            }
        }
    }
    let bytes = PreparedCompressibleVolume3::required_workspace_bytes(1).unwrap();
    let inviscid =
        PreparedCompressibleVolume3::from_mixture(flow.flow(), &[CENTER], settings, bytes, || false)
            .unwrap();
    let mut full = PreparedCompressibleVolume3::from_molecular_mixture(
        &mut flow, &[CENTER], settings, bytes, multi::DissipationControl::default(), || false,
    )
    .unwrap();
    assert_eq!(full.mixture_source(), Some(flow.flow().gas()));
    for (k, tau) in rotated(expected).into_iter().enumerate() {
        close(
            full.stress_moments_n_m()[0][k] - inviscid.stress_moments_n_m()[0][k],
            -tau * grid.cell_volume(),
        );
    }
    assert_eq!(flow.flow().states(), original);
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

#[test]
fn g4_full_stress_mode_and_late_cancellation_preserve_source_and_published_sample() {
    let mut flow = affine_flow(0.12, 0.07, GRADIENT);
    let settings = config(flow.flow().gas());
    let cells = [CENTER];
    let bytes = PreparedCompressibleVolume3::required_workspace_bytes(cells.len()).unwrap();
    assert!(matches!(
        PreparedCompressibleVolume3::from_molecular_frozen(
            &mut flow, &cells, settings, bytes - 1, || false,
        ),
        Err(CompressibleVolumeError3::Workspace { .. }),
    ));
    let mut sample = PreparedCompressibleVolume3::from_molecular_frozen(
        &mut flow, &cells, settings, bytes, || false,
    )
    .unwrap();
    let mut reference = PreparedCompressibleVolume3::from_molecular_frozen(
        &mut flow, &cells, settings, bytes, || false,
    )
    .unwrap();
    let mut inviscid =
        PreparedCompressibleVolume3::from_frozen(flow.flow(), &cells, settings, bytes, || false)
            .unwrap();
    assert_eq!(sample.sample_frozen(flow.flow(), || false), Err(CompressibleVolumeError3::SourceMismatch));
    assert_eq!(
        inviscid.sample_molecular_frozen(&mut flow, || false),
        Err(CompressibleVolumeError3::SourceMismatch),
    );
    let previous = sample.stress_moments_n_m().to_vec();
    flow.step(Time::new(1e-7), single::MolecularStepControl::default(), || false).unwrap();
    let states = flow.flow().states().to_vec();
    let report = flow.last_step().copied();
    let time = flow.flow().time();
    let mut calls = 0;
    reference.sample_molecular_frozen(&mut flow, || {
        calls += 1;
        false
    }).unwrap();
    assert!(calls > 3);
    let mut polls = 0;
    assert_eq!(
        sample.sample_molecular_frozen(&mut flow, || {
            polls += 1;
            polls == calls
        }),
        Err(CompressibleVolumeError3::Cancelled),
    );
    assert_eq!(sample.sample_time_s(), 0.0);
    assert_eq!(sample.stress_moments_n_m(), previous);
    assert_eq!(flow.flow().states(), states);
    assert_eq!(flow.flow().time(), time);
    assert_eq!(flow.last_step().copied(), report);
    sample.sample_molecular_frozen(&mut flow, || false).unwrap();
    assert_eq!(sample.stress_moments_n_m(), reference.stress_moments_n_m());
    assert_eq!(sample.sample_time_s(), time.value());
}
