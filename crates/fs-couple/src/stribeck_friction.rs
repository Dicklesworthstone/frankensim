//! Compatibility surface for the shared constitutive friction owner.
//!
//! Scalar modal and bowed-string consumers retain the same types and calls.
//! Stable-root continuation, its numerical budgets and retained tangential
//! history now belong to `fs-tribo`, so contact/body consumers can use the same
//! implementation through `TangentialLaw` without depending on this coupler.

pub use fs_tribo::stribeck_friction::{
    ContactBranch, HystereticFriction, HystereticFrictionState, ImplicitFrictionError,
    ImplicitFrictionStep, SetValuedFrictionStep, StribeckFriction, StribeckLoadLine,
    solve_implicit_friction,
};
