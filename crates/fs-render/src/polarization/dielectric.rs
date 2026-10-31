//! Polarized dielectric optics and frame-aware ray-power transport.
//!
//! Interface matrices use +Q along s, p chosen so (s,p,propagation) is
//! right-handed, U=2 Re(E_s conj(E_p)), V=2 Im(E_s conj(E_p)), and
//! exp(-i omega t). Reflection reverses p at normal incidence. Transmission
//! matrices transfer power; the legacy sampling adapter includes its radiance
//! eta-squared factor. Only lossless isotropic interfaces are admitted, with
//! native Beer-Lambert absorption as a separate bulk effect.
//!
//! [`transport`] carries the frame and current material across successive
//! physical ray interactions. It is forward power transport, not a replacement
//! for the scalar camera-path integrator.

mod interface;
pub use interface::{
    PolarizedDielectricSample, PolarizedFresnel, attenuate_stokes,
    fresnel_dielectric_polarized, sample_smooth_dielectric_polarized,
    spectral_fresnel_dielectric_polarized,
};
pub mod transport;
