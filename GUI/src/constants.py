APP_TITLE = "TEOAE Cosmic Analyzer"
APP_SUBTITLE = (
    "Explore patients, templates, SNR, reproducibility, and waveform matching "
    "in a fantasy-space interface."
)

DEFAULT_SNR_THRESHOLD = 3.0
DEFAULT_REPEAT_CORR_THRESHOLD = 0.35
DEFAULT_BANDPASS_LOW = 1000.0
DEFAULT_BANDPASS_HIGH = 3800.0
DEFAULT_FILTER_ORDER = 4
DEFAULT_RESPONSE_START = 0.004
DEFAULT_RESPONSE_END = 0.016
DEFAULT_NOISE_START = 0.017
DEFAULT_NOISE_END = 0.020
DEFAULT_MAX_LAG = 20
DEFAULT_HEAD_LEN = 100

COSMIC_CSS = """
<style>
:root {
  --bg1: #050816;
  --bg2: #0b1029;
  --bg3: #1a0f35;
  --card: rgba(10, 14, 38, 0.68);
  --line: rgba(255,255,255,0.12);
  --text: #e8ecff;
  --muted: #b8c1ff;
}
html, body, [data-testid="stAppViewContainer"], .stApp {
  background:
    radial-gradient(circle at 20% 20%, rgba(120,80,255,0.20), transparent 18%),
    radial-gradient(circle at 80% 15%, rgba(0,200,255,0.16), transparent 20%),
    radial-gradient(circle at 70% 80%, rgba(255,120,180,0.12), transparent 18%),
    linear-gradient(180deg, var(--bg1), var(--bg2) 48%, var(--bg3));
  color: var(--text);
}
[data-testid="stHeader"] { background: rgba(0,0,0,0); }
[data-testid="stSidebar"] {
  background: linear-gradient(180deg, rgba(6,10,25,0.96), rgba(18,10,40,0.94));
  border-right: 1px solid var(--line);
}
.cosmic-title {
  font-size: 2.5rem;
  font-weight: 800;
  letter-spacing: 0.04em;
  background: linear-gradient(90deg, #eaf1ff, #88d7ff, #c18cff, #eaf1ff);
  -webkit-background-clip: text;
  -webkit-text-fill-color: transparent;
  margin-bottom: 0.25rem;
}
.cosmic-subtitle { color: #b8c1ff; margin-bottom: 1.2rem; }
[data-testid="stMetric"] {
  background: var(--card);
  border: 1px solid var(--line);
  border-radius: 18px;
  padding: 0.8rem;
}
.stPlotlyChart, [data-testid="stDataFrame"] {
  background: var(--card);
  border: 1px solid var(--line);
  border-radius: 18px;
  padding: 0.5rem;
}
.cosmic-card {
  background: var(--card);
  border: 1px solid var(--line);
  border-radius: 20px;
  padding: 1rem 1.1rem;
}
.stars, .stars2, .stars3 {
  position: fixed;
  inset: 0;
  pointer-events: none;
  z-index: 0;
}
.stars {
  background-image:
    radial-gradient(2px 2px at 20px 30px, rgba(255,255,255,0.7), transparent),
    radial-gradient(1.5px 1.5px at 140px 80px, rgba(136,215,255,0.8), transparent),
    radial-gradient(1.5px 1.5px at 250px 200px, rgba(255,255,255,0.7), transparent),
    radial-gradient(2px 2px at 360px 120px, rgba(193,140,255,0.7), transparent);
  background-size: 520px 320px;
  animation: drift 60s linear infinite;
  opacity: 0.8;
}
.stars2 {
  background-image:
    radial-gradient(2px 2px at 50px 150px, rgba(255,255,255,0.5), transparent),
    radial-gradient(1px 1px at 200px 60px, rgba(136,215,255,0.6), transparent),
    radial-gradient(1px 1px at 300px 220px, rgba(255,255,255,0.5), transparent);
  background-size: 620px 420px;
  animation: drift 90s linear infinite reverse;
  opacity: 0.55;
}
.stars3 {
  background-image:
    radial-gradient(4px 4px at 100px 120px, rgba(255,255,255,0.16), transparent),
    radial-gradient(3px 3px at 420px 90px, rgba(136,215,255,0.18), transparent),
    radial-gradient(5px 5px at 580px 260px, rgba(193,140,255,0.15), transparent);
  background-size: 760px 560px;
  animation: pulse 8s ease-in-out infinite;
}
.planet {
  position: fixed;
  border-radius: 50%;
  pointer-events: none;
  z-index: 0;
}
.planet-a {
  width: 220px; height: 220px; right: -40px; top: 90px;
  background:
    radial-gradient(circle at 30% 30%, rgba(255,255,255,0.35), transparent 22%),
    radial-gradient(circle at 35% 40%, #7d54ff, #3a185f 65%, rgba(0,0,0,0.18) 100%);
}
.planet-b {
  width: 140px; height: 140px; left: -30px; bottom: 100px;
  background:
    radial-gradient(circle at 30% 30%, rgba(255,255,255,0.28), transparent 22%),
    radial-gradient(circle at 35% 40%, #5fe2ff, #0d5078 65%, rgba(0,0,0,0.18) 100%);
}
.ship {
  position: fixed;
  top: 18%; left: -10%; width: 120px; height: 40px;
  z-index: 0; opacity: 0.55; animation: fly 28s linear infinite;
}
.ship:before {
  content: ''; position: absolute; left: 0; top: 11px; width: 90px; height: 18px;
  background: linear-gradient(90deg, rgba(255,255,255,0.3), rgba(136,215,255,0.5), rgba(193,140,255,0.3));
  clip-path: polygon(0 60%, 20% 20%, 75% 20%, 100% 50%, 75% 80%, 20% 80%);
  border-radius: 10px;
}
.ship:after {
  content: ''; position: absolute; left: 82px; top: 15px; width: 55px; height: 8px;
  background: linear-gradient(90deg, rgba(136,215,255,0.6), rgba(255,255,255,0));
  filter: blur(2px);
}
@keyframes drift { from { transform: translateY(0px) translateX(0px);} to { transform: translateY(-120px) translateX(-60px);} }
@keyframes pulse { 0%,100% { opacity: 0.25; transform: scale(1);} 50% { opacity: 0.45; transform: scale(1.08);} }
@keyframes fly { 0% { transform: translateX(0vw) translateY(0px) scale(0.95); opacity: 0;} 6% { opacity: 0.55;} 50% { transform: translateX(60vw) translateY(-40px) scale(1.05);} 100% { transform: translateX(125vw) translateY(-70px) scale(0.9); opacity: 0;} }
</style>
<div class="stars"></div>
<div class="stars2"></div>
<div class="stars3"></div>
<div class="planet planet-a"></div>
<div class="planet planet-b"></div>
<div class="ship"></div>
"""