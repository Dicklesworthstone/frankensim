// Entry of the instruments-only review build (instruments.html): the site's fonts and stylesheet, then the section itself.
// The full page gets the same pieces from main.ts and a second script tag; this entry pulls in nothing else.
import "@fontsource-variable/fraunces/full.css";
import "@fontsource-variable/inter";
import "@fontsource/jetbrains-mono/400.css";
import "@fontsource/jetbrains-mono/500.css";
import "./styles.css";
import "./blind-ab";
import "./candidates";
import "./instruments";
