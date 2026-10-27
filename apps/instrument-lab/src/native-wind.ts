import "./native-wind.css";
// One canonical native fixture, bundled verbatim by Vite; no copied physics defaults.
import exampleSource from "../../../crates/fs-couple/examples/reed-two-medium-tone-hole.performance?raw";
import { nativeWindStatus, renderNativeWind, type NativeWindStatus } from "./native-wind-client.mjs";
import { el } from "./util";

const EXAMPLE_NAME = "reed-two-medium-tone-hole.performance";
const EXAMPLE_URL = "https://github.com/Dicklesworthstone/frankensim/blob/main/crates/fs-couple/examples/reed-two-medium-tone-hole.performance";
const LOCAL_FILE_LIMIT = 1_048_576; // Bound file reading even on a static host without a renderer.
const message = (error: unknown): string => error instanceof Error ? error.message : String(error);

/** A source editor and native render client. Synthesis stays in the Rust renderer. */
export function createNativeWindPanel(): HTMLElement {
  const source = el("textarea", { id: "native-wind-source", class: "nw-source", rows: "15",
    spellcheck: "false", autocapitalize: "off", autocomplete: "off", wrap: "off",
    "aria-describedby": "native-wind-source-help" });
  source.value = exampleSource;
  const sourceInfo = el("span", { class: "small muted" });
  const connection = el("span", { class: "nw-connection", role: "status", "aria-live": "polite" }, "Checking local renderer…");
  const result = el("p", { class: "nw-result small", role: "status", "aria-live": "polite", "aria-atomic": "true" });
  const render = el("button", { type: "button", class: "btn primary" }, "Render WAV");
  const cancel = el("button", { type: "button", class: "btn" }, "Cancel render");
  const check = el("button", { type: "button", class: "btn" }, "Check renderer");
  const open = el("button", { type: "button", class: "btn" }, "Open performance…");
  const reset = el("button", { type: "button", class: "btn" }, "Reset example");
  const file = el("input", { type: "file", accept: ".performance,.txt,text/plain",
    "aria-label": "Open a native wind performance" });
  file.hidden = true;
  const audio = el("audio", { controls: "", preload: "metadata", "aria-label": "Your native wind render" });
  const download = el("a", { class: "btn", download: "wind-render.wav" }, "Download WAV");
  const renderedFrom = el("p", { class: "small muted" });
  const playback = el("div", { class: "nw-playback" }, audio, download, renderedFrom);
  playback.hidden = true;
  cancel.hidden = true;

  let status: NativeWindStatus | null = null;
  let checking: AbortController | null = null;
  let rendering: AbortController | null = null;
  let loadingFile = false;
  let previewUrl: string | null = null;
  let renderedSource: string | null = null;
  let sourceName = EXAMPLE_NAME;

  const controls = (): void => {
    const active = rendering !== null;
    render.disabled = active || loadingFile || checking !== null || !status?.available || status.busy || !source.value.trim();
    cancel.hidden = !active;
    check.disabled = active || checking !== null;
    open.disabled = active || loadingFile;
    reset.disabled = active || loadingFile;
    source.readOnly = active || loadingFile;
    source.setAttribute("aria-busy", String(active));
    render.textContent = active ? "Rendering…" : "Render WAV";
  };
  const edited = (): void => {
    sourceInfo.textContent = `${sourceName} · ${new TextEncoder().encode(source.value).byteLength.toLocaleString()} bytes`;
    if (renderedSource !== null) {
      renderedFrom.textContent = source.value === renderedSource
        ? "This WAV was rendered from the source currently shown."
        : "This WAV is from your previous source. Render again to hear your edits.";
    }
    controls();
  };
  const clearPreview = (): void => {
    audio.pause();
    audio.removeAttribute("src");
    audio.load();
    download.removeAttribute("href");
    playback.hidden = true;
    renderedSource = null;
    if (previewUrl !== null) URL.revokeObjectURL(previewUrl);
    previewUrl = null;
  };
  const showResult = (text: string, refused = false): void => {
    result.textContent = text;
    result.classList.toggle("nw-refused", refused);
  };

  const checkConnection = async (): Promise<void> => {
    checking?.abort();
    const controller = new AbortController();
    checking = controller;
    status = null;
    connection.textContent = "Checking local renderer…";
    controls();
    const timeout = window.setTimeout(() => controller.abort(), 5000);
    try {
      const next = await nativeWindStatus(controller.signal);
      if (checking !== controller) return;
      status = next;
      connection.textContent = !next.available ? "Native renderer unavailable"
        : next.busy ? "Local renderer is busy" : "Local renderer ready";
      if (!next.available) {
        showResult("Open the local Instrument Lab server to render and listen. You can edit the source here.");
      } else if (next.busy) {
        showResult("Another render is running. Check the renderer again when it finishes.");
      } else {
        showResult(`Ready to render. The local render limit is ${Math.ceil(next.timeoutMs / 1000)} seconds.`);
      }
    } catch {
      if (checking !== controller) return;
      connection.textContent = "Native renderer unavailable";
      showResult("This host has no native renderer connection. Open the local Instrument Lab server to render and listen.");
    } finally {
      window.clearTimeout(timeout);
      if (checking === controller) {
        checking = null;
        controls();
      }
    }
  };

  render.addEventListener("click", () => {
    if (!status?.available || status.busy || rendering || checking || loadingFile) return;
    const connected = status;
    const submitted = source.value;
    const controller = new AbortController();
    rendering = controller;
    audio.pause();
    showResult("Rendering your performance from its material, gas and geometry data…");
    controls();
    let timedOut = false;
    const timeout = window.setTimeout(() => {
      timedOut = true;
      controller.abort();
    }, connected.timeoutMs + 5000);
    void (async () => {
      try {
        const wave = await renderNativeWind(submitted, connected, controller.signal);
        if (rendering !== controller || controller.signal.aborted) return;
        const nextUrl = URL.createObjectURL(wave);
        clearPreview();
        previewUrl = nextUrl;
        renderedSource = submitted;
        audio.src = nextUrl;
        download.href = nextUrl;
        download.download = `${sourceName.replace(/\.[^.]*$/, "") || "wind"}-render.wav`;
        playback.hidden = false;
        edited();
        showResult("Render complete. Press play to listen, or download the WAV.");
      } catch (error) {
        if (rendering !== controller) return;
        showResult(controller.signal.aborted
          ? timedOut ? "The render timed out. Check the local renderer before trying again." : "Render cancelled."
          : message(error), !controller.signal.aborted || timedOut);
      } finally {
        window.clearTimeout(timeout);
        if (rendering === controller) {
          rendering = null;
          controls();
        }
      }
    })();
  });
  cancel.addEventListener("click", () => {
    rendering?.abort();
    showResult("Cancelling render…");
  });
  check.addEventListener("click", () => void checkConnection());
  source.addEventListener("input", edited);
  open.addEventListener("click", () => file.click());
  reset.addEventListener("click", () => {
    source.value = exampleSource;
    sourceName = EXAMPLE_NAME;
    edited();
  });
  file.addEventListener("change", () => {
    const selected = file.files?.[0];
    file.value = "";
    if (!selected || rendering || loadingFile) return;
    const limit = Math.min(LOCAL_FILE_LIMIT, status?.maximumSourceBytes ?? LOCAL_FILE_LIMIT);
    if (selected.size > limit) {
      showResult(`Choose a performance file no larger than ${limit.toLocaleString()} bytes.`, true);
      return;
    }
    loadingFile = true;
    controls();
    void selected.arrayBuffer().then((bytes) => {
      const text = new TextDecoder("utf-8", { fatal: true }).decode(bytes);
      source.value = text;
      sourceName = selected.name;
      showResult("Performance loaded. Edit the source or render it.");
      edited();
    }).catch(() => {
      showResult("This performance could not be read as UTF-8 text. Your previous source is still in the editor.", true);
    }).finally(() => {
      loadingFile = false;
      controls();
    });
  });
  window.addEventListener("pagehide", () => {
    checking?.abort();
    rendering?.abort();
    clearPreview();
  });
  window.addEventListener("pageshow", (event) => {
    if (event.persisted) void checkConnection();
  });

  const panel = el("section", { class: "panel nw-panel", "aria-labelledby": "native-wind-title" },
    el("div", { class: "panel-head" },
      el("div", {}, el("h3", { id: "native-wind-title" }, "Render a wind performance"),
        el("p", { class: "small muted" }, "Edit the source, run the native physics, and listen to your result.")),
      connection),
    el("p", { class: "nw-intro" }, "Starter: two-medium tone hole — a short boundary diagnostic with a reed and three tube sections in different gas regions. Change the gas mixtures, temperature, geometry or gestures in the description, or open a self-contained reed or flue performance."),
    el("details", { class: "nw-help small" },
      el("summary", {}, "About the source and local renderer"),
      el("p", { id: "native-wind-source-help" }, "ambient sets the inlet and chimney gas; a segment’s gas entry overrides its local medium; exterior sets the surrounding gas. The interior and exterior must share a static pressure because this boundary model omits mean flow. This short numerical fixture is not a measured clarinet or oboe."),
      el("p", {}, "The local renderer accepts self-contained frankensim-reed-performance-v1 and frankensim-flue-dg-performance-v1 descriptions. File-backed material and valve packs need the command-line workflow."),
      el("p", {}, "Rendering uses the local music_render executable. The hosted review page can show and edit the source; audio is available when you open it through the local renderer server."),
      el("p", {}, el("a", { href: EXAMPLE_URL, target: "_blank", rel: "noopener" }, "Canonical example"),
        " · ", el("a", { href: "https://github.com/Dicklesworthstone/frankensim/tree/main/apps/instrument-lab#native-wind-renderer", target: "_blank", rel: "noopener" }, "Local renderer setup"))),
    el("div", { class: "nw-source-head" }, el("label", { for: "native-wind-source" }, "Performance source"),
      sourceInfo),
    source,
    el("div", { class: "btn-row nw-actions" }, render, cancel, open, reset, check, file),
    result, playback);
  edited();
  void checkConnection();
  return panel;
}
