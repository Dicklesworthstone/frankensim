//! The Smith & Woodhouse (2000) plastic thermal law of rosin is a render law in
//! [`crate::rosin_thermal`] since 2026-11-06. This module keeps its benchmark
//! names for the X-Bw2 bakeoff and the g1 tests.

pub use crate::rosin_thermal::{
    MAX_HISTORY_STEPS, RosinContact, SmithWoodhouse, SmithWoodhouseParams, SmithWoodhouseState,
    SteadySlidingFit, ThermalSolid, half_space_unit_weight,
};
