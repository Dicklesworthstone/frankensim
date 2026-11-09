//! G1/G3: public FE adapter with independent exact P1 tetrahedron integrals.
use super::*;
use crate::acoustics::helmholtz::{AcousticFeError, HelmholtzFe, Medium};
use crate::acoustics::helmholtz::heterogeneous::{HeterogeneousHelmholtzFe, MechanicalHarmonicLimits};
use crate::acoustics::helmholtz::heterogeneous::time::ResistiveBoundary;
use fs_rep_mesh::TetComplex;

#[test]
fn original_heterogeneous_mass_and_robin_heat_set_the_complex_cavity_mode() {
    let points = [[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]];
    let mesh = TetComplex::from_tets(4, vec![[0,1,2,3]]);
    for (rho, c) in [(1.0,1.0),(2.0,3.0)] {
        let medium = Medium { density_kg_m3: rho, sound_speed_m_s: c };
        let geometry = HelmholtzFe::new(&mesh, &points, 1, medium).unwrap();
        let model = HeterogeneousHelmholtzFe::new(&geometry, 3,
            &mut |_,_| Ok::<_,AcousticFeError>(medium), &mut || false).unwrap();
        let admittance = 0.01;
        let robin = ResistiveBoundary { faces: vec![[0,1,2]], admittance_m_pa_s: admittance };
        let owner = model.prepare_with_coupled_mechanical_walls(0.1, &[robin], &[], &[[1,2,3]],
            MechanicalHarmonicLimits::default(), &mut || false).unwrap();
        assert_eq!(owner.resonance_coordinate_count(),1);
        assert_eq!(owner.resonance_free_pressure_dofs(), &[0]);
        // Integral |grad(phi0)|²/rho=1/(2rho), integral phi0²/(rho*c²)=1/(60rho*c²),
        // and integral_boundary Y*phi0²=Y/12. No production matrix is the oracle.
        let decay = 2.5*rho*c*c*admittance;
        let expected = C64::new(-decay, (30.0*c*c-decay*decay).sqrt());
        let options = MechanicalPoleOptions { spectral_scale: expected.abs(), step_tolerance: 1e-12,
            residual_tolerance: 1e-12, max_iterations: 24, max_dimension: 16,
            max_dense_bytes: 2<<20, max_work: 1_000_000 };
        let limits = MechanicalPoleLimits { max_evaluations: 64, max_payload_bytes: 2<<20,
            physical_residual_tolerance: 1e-11 };
        let pole = owner.resonance_pole(expected.scale(0.99), &[C64::ONE], &options, limits, &mut || false).unwrap();
        assert!((pole.rate_per_s-expected).abs() < 2e-10);
        assert_eq!(&pole.pressure_pa[1..], &[C64::ZERO;3]);
        assert!(pole.walls.is_empty());
        assert!(pole.physical_residual < 1e-11);
    }
}
