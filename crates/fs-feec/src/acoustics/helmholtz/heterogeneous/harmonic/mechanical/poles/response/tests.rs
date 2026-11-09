//! Public source projection against independent P1 volume/face integrals.
use super::*;
use crate::acoustics::helmholtz::{AcousticFeError, HelmholtzFe, Medium};
use crate::acoustics::helmholtz::heterogeneous::{
    HeterogeneousHelmholtzFe, MechanicalHarmonicLimits, MechanicalPoleOptions,
};
use crate::acoustics::helmholtz::heterogeneous::time::ResistiveBoundary;
use fs_rep_mesh::TetComplex;

#[test]
fn public_pressure_residue_uses_original_face_motion_and_release_constraints() {
    let points = [[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]];
    let mesh = TetComplex::from_tets(4, vec![[0,1,2,3]]);
    let medium = Medium { density_kg_m3: 2.0, sound_speed_m_s: 3.0 };
    let geometry = HelmholtzFe::new(&mesh, &points, 1, medium).unwrap();
    let model = HeterogeneousHelmholtzFe::new(&geometry, 3,
        &mut |_,_| Ok::<_,AcousticFeError>(medium), &mut || false).unwrap();
    let admittance = 0.01;
    let robin = ResistiveBoundary {
        faces: vec![[0,1,2]], admittance_m_pa_s: admittance,
    };
    let owner = model.prepare_with_coupled_mechanical_walls(0.1, &[robin],
        &[], &[[1,2,3]], MechanicalHarmonicLimits::default(), &mut || false).unwrap();
    let rho = medium.density_kg_m3;
    let c = medium.sound_speed_m_s;
    let decay = 2.5*rho*c*c*admittance;
    let expected = C64::new(-decay, (30.0*c*c-decay*decay).sqrt());
    let pole_options = MechanicalPoleOptions {
        spectral_scale: expected.abs(), step_tolerance: 1e-12,
        residual_tolerance: 1e-12, max_iterations: 24, max_dimension: 16,
        max_dense_bytes: 2<<20, max_work: 1_000_000,
    };
    let limits = MechanicalPoleLimits {
        max_evaluations: 64, max_payload_bytes: 2<<20,
        physical_residual_tolerance: 1e-11,
    };
    let pole = owner.resonance_pole(expected, &[C64::ONE], &pole_options,
        limits, &mut || false).unwrap();
    let options = MechanicalPoleResponseOptions {
        spectral_scale: expected.abs(), pole_residual_tolerance: 1e-11,
        pole_step_tolerance: 1e-11, linear_residual_tolerance: 1e-11,
        max_dimension: 16, max_vectors: 1, max_dense_bytes: 2<<20,
        max_work: 1_000_000,
    };
    let prepared = owner.prepare_resonance_response(&pole, options,
        limits, &mut || false).unwrap();
    let speed = C64::new(0.7, -0.3);
    let residue = prepared.forcing_residue(&[([0,1,2], [speed;3])],
        &[], &[], &[], &mut || false).unwrap();
    // Integral phi0 = area/3=1/6 on the driven face; M00=1/(60*rho*c^2),
    // C00=Y/12. The physical force is s*speed/6.
    let slope = expected.scale(2.0/(60.0*rho*c*c))+C64::from_re(admittance/12.0);
    let expected_residue = expected*speed.scale(1.0/6.0)/slope;
    assert!((residue.pressure_pa_per_s[0]-expected_residue).abs()
        < 2e-9*expected_residue.abs());
    assert_eq!(&residue.pressure_pa_per_s[1..], &[C64::ZERO;3]);
    assert!(residue.walls.is_empty());
    assert!(prepared.forcing_residue(&[([1,2,3], [speed;3])],
        &[], &[], &[], &mut || false).is_err());
    assert!(prepared.forcing_residue(&[], &[], &[], &[], &mut || true).is_err());
    let zero = prepared.forcing_residue(&[], &[], &[], &[], &mut || false).unwrap();
    assert!(zero.coordinates.iter().all(|&z| z == C64::ZERO));
}
