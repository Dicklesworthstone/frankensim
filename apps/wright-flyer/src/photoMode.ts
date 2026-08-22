// Photo mode (presentation-plane only): one body class toggled by the
// P key. The class hides every HUD surface and applies a period plate
// grade + frame border via CSS filters on the renderer canvas — no
// render-pipeline coupling, so it composes with any QoS tier and any
// post-processing stack. State is trivially serializable: active|not.

export interface PhotoModeToggle {
  /** True after the toggle if photo mode is now ON. */
  readonly active: boolean;
}

const PHOTO_CLASS = "wf-photo";

/** Toggle photo mode on a document root. Returns the new state. Pure
 * DOM-class mutation; safe to call repeatedly. */
export function togglePhotoMode(root: Document): PhotoModeToggle {
  const active = !root.body.classList.contains(PHOTO_CLASS);
  root.body.classList.toggle(PHOTO_CLASS, active);
  return { active };
}

/** Current state without toggling. */
export function photoModeActive(root: Document): boolean {
  return root.body.classList.contains(PHOTO_CLASS);
}
