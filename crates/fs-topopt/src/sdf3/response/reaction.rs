//! Numerical actuator/support-force and moment targets for response design.
use fs_cutfem::elastic3::adaptive::dirichlet::PrescribedMotion3;

/// One scalar Nitsche-consistent embedded reaction observation. The mode acts
/// only on the operator's selected embedded supports. Constant unit vectors
/// select forces; rotation modes select moments. Targets/scales must use the
/// resulting physical units. Modes are pure, reference-configuration laws and
/// must not change with density. Values are not actuator-energy certificates.
#[derive(Clone, Copy)]
pub struct ReactionTarget3<'a> {
    pub mode: &'a PrescribedMotion3<'a>,
    pub target: f64,
    pub scale: f64,
    pub weight: f64,
}
