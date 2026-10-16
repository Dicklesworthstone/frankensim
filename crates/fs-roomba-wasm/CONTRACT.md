# CONTRACT: fs-roomba-wasm

## Purpose and layer

Layer: **L6 WASM boundary**. This crate exposes the US 6,594,844 Roomba
museum composition owned by `fs-mbd::roomba` to browser consumers. It performs
packet admission and stable JSON-envelope serialization only. Generic
constant-twist differential-drive law remains in `fs-mbd::planar_drive`; the
source-bounded optical redirect and display non-penetration composition remain
in `fs-mbd::roomba`.

## Public boundary

`roomba_step(&[f64]) -> String` accepts packet version 1:

| Index | Field |
| --- | --- |
| 0 | packet version (`1`) |
| 1 | fixed `dt_s` |
| 2 | wheel speed, m/s |
| 3 | turn rate, rad/s |
| 4–5 | room width/height, m |
| 6 | optical sensor height, inches |
| 7 | explicit wall distance, inches, or `-1` to derive it from geometry |
| 8 | optical subsystem enabled (`0` or `1`) |
| 9–11 | previous x/y pose, m, and heading, rad |
| 12 | mode (`0` spiral, `1` straight, `2` turn, `3` backup) |
| 13 | time in mode, s |
| 14 | exact unsigned 32-bit random seed |
| 15–17 | left wheel, right wheel, and side-brush angles, rad |
| 18… | zero or more collider quadruples `(x_m, y_m, width_m, height_m)` |

Success returns `{"ok":{...}}` with the complete browser tape state, stable
mode/reason strings, and a contact index (`-1` clear, `-2` room boundary, or
the zero-based collider index). Refusal returns
`{"refusal":{"code", "message", "repairs"}}` and never a partial `ok`
payload.

## Invariants and refusal

- Packet length, version, enum/boolean encodings, integer seed, collider stride,
  finite values, room geometry, and fixed-step bounds are admitted before a
  result is serialized.
- At most 64 colliders are forwarded to the bounded owner.
- A successful output contains only finite numeric JSON values.
- The boundary is stateless. Deterministic replay is owned by the caller's
  explicit previous-state packet and the deterministic `fs-mbd` owner.

## Cancellation and resources

The call is synchronous, bounded by 64 colliders, performs no I/O, starts no
threads or tasks, and has no partial commit. Cancellation is therefore the
call boundary: callers may decline to start the next fixed step.

## Unsafe boundary and dependencies

No unsafe code. The crate depends only on `fs-mbd` plus target-specific
`wasm-bindgen`; it owns no duplicate physics law.

## No-claim boundaries

This boundary does not claim collision forces, impulse/restitution, friction,
traction, wheel slip, motor or battery performance, dust pickup, coverage,
localization accuracy, or hardware validation. The room projection and
contextual cleaning path are pedagogical display mechanics, while the optical
field intersection and change-of-direction logic are the source-bounded patent
lane.
