// Development entry for the live drum section (live-drum-dev.html only).
import '../styles.css';
import { mountLiveDrum } from './drumSection';

const root = document.getElementById('live-drum');
if (root) mountLiveDrum(root);
