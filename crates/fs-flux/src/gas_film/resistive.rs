//! Compatibility surface for reciprocal quasistatic Reynolds-film resistance.
//!
//! The single constitutive implementation lives in `fs-tribo`, alongside other
//! solver-independent interface laws. Mechanical coupling can consume that
//! lower owner directly without importing the continuum flow solver stack.
//! Public types and numerical behavior are shared exactly with this surface.

pub use fs_tribo::resistive_film::{
    FilmCell, FilmChannel, FilmError, FilmLimits, FilmReport, GapPort, MAX_CELLS, MAX_PORTS,
    ResistiveFilm,
};

#[cfg(test)]
mod tests;
