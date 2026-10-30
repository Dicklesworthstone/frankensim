use super::*;
use fs_feec::acoustics::helmholtz::Medium;
use fs_feec::fixtures::single_tet;
use fs_material::gas::{GasSpec, GasState};
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw, RelaxingMedium,
    SpeciesRelaxation,
};
use fs_thermochem::SpeciesId;

const GEOMETRY: BoundaryLayerGeometry = BoundaryLayerGeometry {
    minimum_curvature_radius_m: f64::INFINITY,
    minimum_half_gap_m: 0.05,
    maximum_layer_ratio: 0.05,
};

pub(super) fn source(gas: GasState, tau: Option<f64>) -> RelaxationState {
    let species = SpeciesId::new("verification-gas").unwrap();
    let processes = tau
        .map(|tau| {
            vec![RelaxationProcess {
                reservoir: InternalReservoir::Rotational {
                    degrees_of_freedom: 2,
                },
                rate: RelaxationRateLaw::PartnerIndependent {
                    p_tau_ref_pa_s: tau * gas.pressure,
                    t_ref_k: gas.temperature,
                    exponent: 0.0,
                    t_min_k: 200.0,
                    t_max_k: 400.0,
                },
            }]
        })
        .unwrap_or_default();
    RelaxationAbsorption::new(vec![SpeciesRelaxation {
        species: species.clone(),
        processes,
    }])
    .unwrap()
    .evaluate(&RelaxingMedium::new(gas, vec![(species, 1.0)]).unwrap())
    .unwrap()
}

#[test]
fn original_molecular_source_drives_spatial_compression_and_duct_coordinates() {
    let spec = GasSpec {
        gamma: 1.25,
        ..GasSpec::dry_air_ussa1976()
    };
    let gas = GasState::try_new(&spec, 293.15, 101325.0).unwrap();
    let (mesh, points) = single_tet();
    let points: Vec<_> = points.iter().map(|p| p.map(|v| 0.2 * v)).collect();
    let fe = HelmholtzFe::new(
        &mesh,
        &points,
        1,
        Medium {
            density_kg_m3: gas.density,
            sound_speed_m_s: gas.sound_speed,
        },
    )
    .unwrap();
    let faces = fe.boundary_faces();
    let basis = HelmholtzModes {
        lambda: vec![0.0],
        frequency_hz: vec![0.0],
        residual: vec![0.0],
        shapes: vec![vec![1.0; fe.ndof()]],
        below_low: 0,
        below_high: 1,
    };
    let ports = [
        BoundaryMode::uniform(&faces[..1]),
        BoundaryMode::uniform(&faces[1..2]),
    ];
    let mut pressures = Vec::new();
    for tau in [None, Some(0.0002), Some(0.002)] {
        let state = source(gas, tau);
        let expected = state
            .caloric_response(core::f64::consts::TAU * 300.0)
            .unwrap()
            .bulk_modulus_pa
            .conj()
            * C64::new(
                0.0,
                -1.0 / (core::f64::consts::TAU * 300.0 * (0.2_f64.powi(3) / 6.0)),
            );
        let model = MolecularInterior::new(
            &fe,
            &basis,
            &ports,
            &faces[2..],
            state.clone(),
            MolecularThermalWall::Adiabatic,
            GEOMETRY,
            &|| false,
        )
        .unwrap();
        assert_eq!(model.source(), &state);
        assert!(model.wall_source().is_none());
        let response = model
            .prepare(300.0, &|| false)
            .unwrap()
            .response(&[C64::from_re(1e-7), C64::ZERO], &|| false)
            .unwrap();
        let z = response.acoustic.pressure_pa[0].scale(1e7);
        assert!((z - expected).abs() / expected.abs() < 2e-12);
        if tau.is_some() {
            let input = response.acoustic.input_power_w();
            assert!(response.dissipation.molecular_w > 0.0);
            assert!((input - response.dissipation.total_w()).abs() < 1e-10 * input);
        }
        let r = 1e6;
        let s = model.duct_scattering(300.0, &[r, r], &|| false).unwrap();
        let z_duct = expected.conj();
        let reflected = z_duct / (z_duct.scale(2.0) + C64::from_re(r));
        // Equal-pressure two-port rank-one Z: S11=-R/(2Z+R), S12=2Z/(2Z+R).
        assert!((s[1] - reflected.scale(2.0)).abs() < 2e-12);
        assert!((s[0] - (reflected.scale(2.0) - C64::ONE)).abs() < 2e-12);
        assert!(model.duct_two_port(300.0, r, &|| false).is_ok());
        pressures.push(z);
    }
    assert!(pressures.windows(2).all(|p| (p[1] - p[0]).abs() > 1e-3));
    let wrong = HelmholtzFe::new(
        &mesh,
        &points,
        1,
        Medium {
            density_kg_m3: gas.density * 1.001,
            sound_speed_m_s: gas.sound_speed,
        },
    )
    .unwrap();
    assert!(
        MolecularInterior::new(
            &wrong,
            &basis,
            &ports,
            &faces[2..],
            source(gas, Some(0.001)),
            MolecularThermalWall::Adiabatic,
            GEOMETRY,
            &|| false
        )
        .is_err()
    );
}
