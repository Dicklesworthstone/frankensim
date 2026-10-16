// Development entry for the live piano section (live-piano-dev.html only).
import '../styles.css';
import { mountLivePiano } from './pianoSection';

const root = document.getElementById('live-piano');
if (root) mountLivePiano(root);
