//! Reusable drivers for the executable physical scan-design examples.
//! Numerical models and optimization kernels remain in the original fs crates.
pub mod population;
pub mod scenarios;
/// Explicit caller-supplied geometry and thermal data for measurement fitting.
pub mod specimen;
/// Explicit multi-heater schedules and selected authored power controls.
pub mod heaters;
/// File-declared physical requirements for measured scan programs.
pub mod requirements;
