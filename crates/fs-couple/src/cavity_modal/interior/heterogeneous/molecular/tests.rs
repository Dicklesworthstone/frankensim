//! Original public molecular sources, frozen FE mass and causal full fields.
use super::*;
use fs_material::gas::{GasSpec, GasState};
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw, RelaxingMedium,
    SpeciesRelaxation, VibrationalMode,
};
use fs_rep_mesh::TetComplex;
use fs_thermochem::SpeciesId;

const X: [[f64; 3]; 5] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [0.0, 0.0, -1.0],
];
const GEOMETRY_MEDIUM: Medium = Medium {
    density_kg_m3: 123.0,
    sound_speed_m_s: 456.0,
};
fn source(temperature: f64, pressure: f64, theta: f64) -> RelaxationState {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature, pressure).unwrap();
    let id = SpeciesId::new("explicit-fixture").unwrap();
    let medium = RelaxingMedium::new(gas, vec![(id.clone(), 1.0)]).unwrap();
    RelaxationAbsorption::new(vec![SpeciesRelaxation {
        species: id,
        processes: vec![RelaxationProcess {
            reservoir: InternalReservoir::Vibrational(vec![VibrationalMode {
                theta_k: theta,
                degeneracy: 1,
            }]),
            rate: RelaxationRateLaw::PartnerIndependent {
                p_tau_ref_pa_s: 101325.0 * 0.002,
                t_ref_k: 293.15,
                exponent: 0.0,
                t_min_k: 250.0,
                t_max_k: 400.0,
            },
        }],
    }])
    .unwrap()
    .evaluate(&medium)
    .unwrap()
}
fn entry(m: &fs_sparse::Csr, i: usize, j: usize) -> f64 {
    let (c, v) = m.row(i);
    c.binary_search(&j).map_or(0.0, |at| v[at])
}

#[test]
fn source_partition_builds_original_frozen_compressibility_and_density_once() {
    let mesh = TetComplex::from_tets(5, vec![[0, 1, 2, 3], [0, 2, 1, 4]]);
    let g = HelmholtzFe::new(&mesh, &X, 1, GEOMETRY_MEDIUM).unwrap();
    let a = source(293.15, 101325.0, 600.0);
    let b = source(353.15, 101325.0, 1600.0);
    let regions = [
        MolecularBulkRegion {
            tetrahedra: vec![0],
            source: &a,
        },
        MolecularBulkRegion {
            tetrahedra: vec![1],
            source: &b,
        },
    ];
    let interior = MolecularBulkInterior::new(&g, &regions, 3, &mut || false).unwrap();
    assert!(core::ptr::eq(interior.sources()[0].source, &a));
    assert_eq!(interior.compressibility_laws()[1].source(), &b);
    assert_eq!(interior.region_for_tetrahedron(1), Some(1));
    for (source, apex) in [(&a, 3), (&b, 4)] {
        let cv = source.frozen_cv_j_kg_k();
        let cp = cv + source.gas().specific_gas_constant;
        let beta = cv / (cp * source.gas().pressure);
        assert!((entry(interior.model().mass(), apex, apex) / (beta / 60.0) - 1.0).abs() < 2e-14);
        // Each apex shape has |grad phi|=1, hence K_aa=V/rho.
        assert!(
            (entry(interior.model().stiffness(), apex, apex)
                / (1.0 / (6.0 * source.gas().density))
                - 1.0)
                .abs()
                < 2e-14
        );
        assert!(
            (beta
                / source
                    .caloric_response(0.0)
                    .unwrap()
                    .bulk_modulus_pa
                    .re
                    .recip()
                - 1.0)
                .abs()
                > 1e-3
        );
    }
    let owner = interior
        .prepare_time(1e-5, &[], &[], &[], BulkTimeLimits::default(), &mut || {
            false
        })
        .unwrap();
    assert_eq!(owner.arm_count(), 2);
    assert_eq!(owner.arm_dofs(0).unwrap().len(), 4);
    assert!(MolecularBulkInterior::new(&g, &regions[..1], 3, &mut || false).is_err());
    let overlap = [regions[0].clone(), regions[0].clone()];
    assert!(MolecularBulkInterior::new(&g, &overlap, 3, &mut || false).is_err());
    let high_pressure = source(353.15, 100000.0, 1600.0);
    let incompatible = [
        regions[0].clone(),
        MolecularBulkRegion {
            tetrahedra: vec![1],
            source: &high_pressure,
        },
    ];
    assert!(MolecularBulkInterior::new(&g, &incompatible, 3, &mut || false).is_err());
}

#[test]
fn full_pressure_fields_converge_to_independent_closed_gas_relaxation() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let g = HelmholtzFe::new(&mesh, &X[..4], 1, GEOMETRY_MEDIUM).unwrap();
    let gas = source(293.15, 101325.0, 600.0);
    let interior = MolecularBulkInterior::new(
        &g,
        &[MolecularBulkRegion {
            tetrahedra: vec![0],
            source: &gas,
        }],
        3,
        &mut || false,
    )
    .unwrap();
    let cv = gas.frozen_cv_j_kg_k();
    let r = gas.gas().specific_gas_constant;
    let c = gas.processes()[0].heat_capacity_j_kg_k;
    let beta_frozen = cv / ((cv + r) * gas.gas().pressure);
    let beta_equilibrium = (cv + c) / ((cv + c + r) * gas.gas().pressure);
    // Direct one-reservoir Cv inversion at constant compression, independent
    // of the new pressure realization and its midpoint algebra.
    let relaxation_time = gas.processes()[0].relaxation_time_s * cv / (cv + c);
    let initial_pressure = 10.0;
    let equilibrium = initial_pressure * beta_frozen / beta_equilibrium;
    let duration = 0.003;
    let expected = equilibrium
        + (initial_pressure - equilibrium) * fs_math::det::exp(-duration / relaxation_time);
    let mut errors = Vec::new();
    for steps in [20, 40, 80] {
        let owner = interior
            .prepare_time(
                duration / steps as f64,
                &[],
                &[],
                &[],
                BulkTimeLimits::default(),
                &mut || false,
            )
            .unwrap();
        let mut state = owner
            .state(
                &[0.0; 4],
                &[initial_pressure; 4],
                &[vec![0.0; 4]],
                &mut || false,
            )
            .unwrap();
        let e0 = owner.energy(&state, &mut || false).unwrap().total_j();
        let mut heat = 0.0;
        for _ in 0..steps {
            let next = owner
                .step_projected(&state, &[0.0; 4], &mut || false)
                .unwrap();
            assert!(next.work.molecular_dissipation_j >= 0.0);
            heat += next.work.molecular_dissipation_j;
            state = next.state;
        }
        let e1 = owner.energy(&state, &mut || false).unwrap().total_j();
        assert!((e1 + heat - e0).abs() < 2e-11 * e0);
        let pressure = state.pressure_pa().iter().sum::<f64>() / 4.0;
        assert!(
            state
                .pressure_pa()
                .iter()
                .all(|&p| (p - pressure).abs() < 2e-10)
        );
        errors.push((pressure - expected).abs());
    }
    assert!(
        errors[0] / errors[1] > 3.8 && errors[1] / errors[2] > 3.8,
        "{errors:?}"
    );
    assert!(errors[2] < 5e-5);
    eprintln!("original-source closed-volume pressure errors {errors:?}");
}
